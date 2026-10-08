#!/usr/bin/env python3
"""Plan the ONT whole-genome pass from the HPRC release 2 ONT file index.

Every R10.4.1 sample with Lai et al. or 1000 Genomes panel truth gets one input
BAM, named <sample>__R1041. Its runs (one PromethION flow cell each) are added
in run order with the same depth rule as the PacBio pass (select_inputs.py rule
8): stop at TARGET_DEPTH, leave out a run that would pass MAX_DEPTH once the
total is above MIN_DEPTH, then shed the smallest runs while the total stays
above MIN_DEPTH. Depth is the HPRC index value per run.

Writes, under data/ont/:
  HPRC_ONT_release2.index.csv    the HPRC index as downloaded
  merge_plan.tsv                 input_bam, sample, prep, filename, size_gb, cov_x, note
  manifest.tsv                   the same columns as data/manifest.tsv, for the scorer
  download_manifest.tsv          sample, filename, s3_uri, bytes (for 1_download.sbatch)
  stages/stageN.*                the same three files for each stage below

  build_ont_plan.py [--refresh]   --refresh downloads the index again
"""
import csv
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = lambda *p: os.path.join(REPO, "data", *p)
OUT = D("ont")
INDEX_URL = ("https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/"
             "main/data_tables/sequencing_data/data_ont_release2_v1.0.index.csv")
INDEX = os.path.join(OUT, "HPRC_ONT_release2.index.csv")
S3_HTTPS = "https://s3-us-west-2.amazonaws.com/"

MIN_DEPTH = 20.0
TARGET_DEPTH = 30.0
MAX_DEPTH = 40.0

# Stage 0 is the pilot. Stage 1 is the rest of the Lai set (four-field truth).
# Stage 2 is ten 1000 Genomes samples across the range of total depth in the
# index. Stage 3 is every other sample.
STAGES = {
    0: ["NA19240"],
    1: ["HG02818", "HG03486", "NA18906", "NA20129", "NA21309"],
    2: ["NA18944", "HG00290", "NA19036", "NA20346", "HG01786",
        "HG00146", "NA20809", "HG01167", "NA18565", "HG03583"],
}


def truth_samples():
    s = set()
    for key in ("lai_answer_key.tsv", "kgp_answer_key.tsv"):
        for line in open(D(key)):
            if line.strip():
                s.add(line.split("\t", 1)[0])
    return s


def head_bytes(uri):
    req = urllib.request.Request(S3_HTTPS + uri[len("s3://"):], method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers["Content-Length"])


def plan_sample(runs):
    runs = sorted(runs, key=lambda r: r["filename"])
    taken, total = [], 0.0
    for r in runs:
        c = float(r["coverage"])
        if taken and total >= TARGET_DEPTH:
            continue
        if taken and total >= MIN_DEPTH and total + c > MAX_DEPTH:
            continue
        taken.append(r)
        total += c
    while total > MAX_DEPTH and len(taken) > 1:
        smallest = min(taken, key=lambda r: float(r["coverage"]))
        if total - float(smallest["coverage"]) < MIN_DEPTH:
            break
        taken.remove(smallest)
        total -= float(smallest["coverage"])
    note = "%d of %d runs (~%.0fx); %s" % (len(taken), len(runs), total, taken[0]["basecaller_model"])
    return taken, note


def write(prefix, entries):
    with open(prefix + "merge_plan.tsv", "w") as fh:
        fh.write("input_bam\tsample\tprep\tfilename\tsize_gb\tcov_x\tnote\n")
        for e in entries:
            for r in e["runs"]:
                fh.write("\t".join([e["id"], e["sample"], "R1041", r["filename"],
                                    "%.1f" % (r["bytes"] / 1e9), r["coverage"], e["note"]]) + "\n")
    with open(prefix + "manifest.tsv", "w") as fh:
        fh.write("input_bam\tsample\tcohort\tlibrary\tinstrument\tfilename\ts3_uri\tbytes\tcov_x\n")
        for e in entries:
            for r in e["runs"]:
                cohort = r["path"].split("/")[4]
                fh.write("\t".join([e["id"], e["sample"], cohort, r["library_id"], r["instrument_model"],
                                    r["filename"], r["path"], str(r["bytes"]), r["coverage"]]) + "\n")
    with open(prefix + "download_manifest.tsv", "w") as fh:
        for e in entries:
            for r in e["runs"]:
                fh.write("\t".join([e["sample"], r["filename"], r["path"], str(r["bytes"])]) + "\n")


def main(argv):
    os.makedirs(os.path.join(OUT, "stages"), exist_ok=True)
    if "--refresh" in argv or not os.path.exists(INDEX):
        urllib.request.urlretrieve(INDEX_URL, INDEX)

    truth = truth_samples()
    by = {}
    for r in csv.DictReader(open(INDEX)):
        if r["sequencing_chemistry"] == "R1041" and r["sample_id"] in truth:
            by.setdefault(r["sample_id"], []).append(r)

    entries = []
    for s in sorted(by):
        if len({r["basecaller_model"] for r in by[s]}) > 1:
            sys.exit("%s mixes basecaller models; plan it by hand" % s)
        runs, note = plan_sample(by[s])
        entries.append({"id": s + "__R1041", "sample": s, "runs": runs, "note": note})

    runs = [r for e in entries for r in e["runs"]]
    with ThreadPoolExecutor(16) as pool:
        for r, b in zip(runs, pool.map(head_bytes, [r["path"] for r in runs])):
            r["bytes"] = b

    write(os.path.join(OUT, ""), entries)
    staged = {s for v in STAGES.values() for s in v}
    stages = dict(STAGES)
    stages[3] = [e["sample"] for e in entries if e["sample"] not in staged]
    for n, samples in stages.items():
        missing = set(samples) - set(by)
        if missing:
            sys.exit("stage %d samples not in the plan: %s" % (n, " ".join(sorted(missing))))
        write(os.path.join(OUT, "stages", "stage%d." % n), [e for e in entries if e["sample"] in samples])

    for n, samples in stages.items():
        sel = [e for e in entries if e["sample"] in samples]
        tb = sum(r["bytes"] for e in sel for r in e["runs"]) / 1e12
        print("stage %d: %d input BAMs, %d runs, %.2f TB" % (n, len(sel), sum(len(e["runs"]) for e in sel), tb))
    print("all: %d input BAMs, %d runs, %.2f TB" % (len(entries), len(runs), sum(r["bytes"] for r in runs) / 1e12))


if __name__ == "__main__":
    main(sys.argv[1:])
