#!/usr/bin/env python3
"""Plan ONT stage 4: the Lai samples whose R10.4.1 runs are not in the release 2 index.

The release 2 index (build_ont_plan.py) stops at January 2024. HPRC has since
sequenced most of the Lai samples again on R10.4.1 (ultra-long runs, dorado
sup 5.0.0 or 5.2.0). Those runs are listed only in the per-submission data
tables of github.com/human-pangenomics/HPRC_metadata. This script reads those
tables, takes every Lai sample that stages 0 to 3 did not type, and plans one
input BAM per sample, <sample>__R1041, with the same depth rule as
build_ont_plan.py.

Some table paths point to the wrong cohort folder (HPRC for HPRC_PLUS), and
one run is only under submissions/. Each run is therefore resolved with a HEAD
request against the table path, the other cohort folder, and the submission
path from the release table. A run that resolves nowhere stops the script.

The samples are split into batches so that only a few batches are on scratch
at a time (run_batches.sh). Batch 1 is the smallest sample, as a pilot. The
other batches are filled in sample order up to --batch_tb of downloads.

Writes, under data/ont/:
  hprc_metadata/                 the HPRC_METADATA tables as downloaded
  stages/stage4.*                merge_plan, manifest and download_manifest for all of stage 4
  stages/stage4/bNN.*            the same three files for each batch

  build_ont_stage4.py [--refresh] [--batch_tb 1.0]
"""
import argparse
import csv
import os
import sys
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from build_ont_plan import MAX_DEPTH, MIN_DEPTH, TARGET_DEPTH, plan_sample

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = lambda *p: os.path.join(REPO, "data", *p)
OUT = D("ont")
META = D("ont", "hprc_metadata")
RAW = "https://raw.githubusercontent.com/human-pangenomics/HPRC_metadata/main/submissions/"
S3_HTTPS = "https://human-pangenomics.s3.us-west-2.amazonaws.com/"

# Every R10.4.1 submission folder in HPRC_metadata. Year2/Year3 R10 overlap the
# release 2 index; their samples were typed in stages 0 to 3 and drop out below.
SUBMISSIONS = {
    "UCSC_HPRC_ONT_Year1_R10": ["data_table", "release"],
    "UCSC_HPRC_PLUS_nanopore_R10": ["data_table", "release"],
    "UCSC_HPRC_nanopore_Year2_R10": ["data_table"],
    "UCSC_HPRC_nanopore_Year3_R10": ["data_table"],
    "UCSC_HPRC_nanopore_Year5": ["data_table"],
}


def fetch(refresh):
    os.makedirs(META, exist_ok=True)
    paths = []
    for sub, kinds in SUBMISSIONS.items():
        for kind in kinds:
            name = "%s_%s.csv" % (sub, kind)
            p = os.path.join(META, name)
            if refresh or not os.path.exists(p):
                urllib.request.urlretrieve(RAW + sub + "/" + name, p)
            paths.append((kind, p))
    return paths


def read_csv(p):
    with open(p, newline="") as fh:
        return [{k.strip(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(fh)]


def lai_samples():
    return {l.split("\t", 1)[0] for l in open(D("lai_answer_key.tsv")) if l.strip()}


def typed_samples():
    """Samples already planned in stages 0 to 3 (data/ont/merge_plan.tsv)."""
    with open(D("ont", "merge_plan.tsv"), newline="") as fh:
        return {r["sample"] for r in csv.DictReader(fh, delimiter="\t")}


def head(uri):
    req = urllib.request.Request(S3_HTTPS + uri[len("s3://human-pangenomics/"):], method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return int(r.headers["Content-Length"])
    except urllib.error.HTTPError as e:
        if e.code in (403, 404):
            return None
        raise


def candidates(r, submitted, folders):
    p = r["path"]
    out = [p]
    for a, b in (("/working/HPRC/", "/working/HPRC_PLUS/"), ("/working/HPRC_PLUS/", "/working/HPRC/")):
        if a in p:
            out.append(p.replace(a, b))
    out += submitted.get(r["filename"], [])
    out += ["%s%s/basecalled/%s" % (f, r["sample_id"], r["filename"]) for f in folders]
    return out


def resolve(r, submitted, folders):
    for uri in candidates(r, submitted, folders):
        b = head(uri)
        if b:
            return uri, b
    return None, None


def cohort(uri):
    parts = uri.split("/")
    if parts[3] == "working":
        return parts[4]
    return "HPRC_PLUS" if "HPRC_PLUS" in parts[4] else "HPRC"


def write(prefix, entries):
    os.makedirs(os.path.dirname(prefix), exist_ok=True)
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
                fh.write("\t".join([e["id"], e["sample"], cohort(r["uri"]), r["library_id"], r["instrument_model"],
                                    r["filename"], r["uri"], str(r["bytes"]), r["coverage"]]) + "\n")
    with open(prefix + "download_manifest.tsv", "w") as fh:
        for e in entries:
            for r in e["runs"]:
                fh.write("\t".join([e["sample"], r["filename"], r["uri"], str(r["bytes"])]) + "\n")


def batches(entries, batch_tb):
    size = lambda e: sum(r["bytes"] for r in e["runs"])
    pilot = min(entries, key=size)
    out, cur, tot = [[pilot]], [], 0
    for e in entries:
        if e is pilot:
            continue
        if cur and tot + size(e) > batch_tb * 1e12:
            out.append(cur)
            cur, tot = [], 0
        cur.append(e)
        tot += size(e)
    if cur:
        out.append(cur)
    return out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--refresh", action="store_true", help="download the HPRC_metadata tables again")
    ap.add_argument("--batch_tb", type=float, default=1.0, help="downloads per batch, TB (default 1.0)")
    a = ap.parse_args(argv)

    tables = fetch(a.refresh)
    # The release tables have no header: submission path, working path.
    rows, submitted, folders = [], {}, set()
    for kind, p in tables:
        if kind == "data_table":
            rows += read_csv(p)
            continue
        for row in csv.reader(open(p, newline="")):
            for v in (x.strip() for x in row):
                if v.startswith("s3://") and "/submissions/" in v:
                    submitted.setdefault(os.path.basename(v), []).append(v)
                    folders.add("/".join(v.split("/")[:5]) + "/")

    want = lai_samples() - typed_samples()
    by = {}
    for r in rows:
        if "r10.4.1" in r["basecaller_model"] and r["sample_id"] in want and r["filename"].endswith(".bam"):
            by.setdefault(r["sample_id"], {})[r["filename"]] = r
    by = {s: list(v.values()) for s, v in by.items()}

    entries = []
    for s in sorted(by):
        if len({r["basecaller_model"] for r in by[s]}) > 1:
            sys.exit("%s mixes basecaller models; plan it by hand" % s)
        runs, note = plan_sample(by[s])
        entries.append({"id": s + "__R1041", "sample": s, "runs": runs, "note": note})

    runs = [r for e in entries for r in e["runs"]]
    with ThreadPoolExecutor(16) as pool:
        for r, (uri, b) in zip(runs, pool.map(lambda r: resolve(r, submitted, folders), runs)):
            if not uri:
                sys.exit("not on S3 at any known path: %s" % r["filename"])
            r["uri"], r["bytes"] = uri, b

    write(os.path.join(OUT, "stages", "stage4."), entries)
    bs = batches(entries, a.batch_tb)
    for i, b in enumerate(bs, 1):
        write(os.path.join(OUT, "stages", "stage4", "b%02d." % i), b)

    print("depth rule %g/%g/%gx; Lai samples without R10.4.1 runs: %s" % (
        MIN_DEPTH, TARGET_DEPTH, MAX_DEPTH, " ".join(sorted(want - set(by))) or "none"))
    for i, b in enumerate(bs, 1):
        tb = sum(r["bytes"] for e in b for r in e["runs"]) / 1e12
        print("b%02d: %2d input BAMs, %2d runs, %.2f TB  %s" % (
            i, len(b), sum(len(e["runs"]) for e in b), tb, " ".join(e["sample"] for e in b)))
    print("stage 4: %d input BAMs, %d runs, %.2f TB" % (len(entries), len(runs), sum(r["bytes"] for r in runs) / 1e12))


if __name__ == "__main__":
    main(sys.argv[1:])
