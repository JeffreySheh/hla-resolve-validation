#!/usr/bin/env python3
"""Plan ONT stage 4: the Lai samples whose R10.4.1 runs are not in the release 2 index.

The release 2 index (build_ont_plan.py) stops at January 2024. HPRC has since
sequenced most of the Lai samples again on R10.4.1 (ultra-long runs, dorado
sup 5.0.0 or 5.2.0). Those runs are listed only in the per-submission data
tables of github.com/human-pangenomics/HPRC_metadata. This script reads those
tables, takes every Lai sample that stages 0 to 3 did not type, and plans one
input BAM per sample, <sample>__R1041, with the same depth rule as
build_ont_plan.py.

Only runs that HPRC has released to working/ are used. Some table paths point
to the wrong cohort folder (HPRC for HPRC_PLUS), so each run is resolved with a
HEAD request against the table path and the other cohort folder. A run found
in neither is not released and is left out before the depth rule is applied.

With --write_batches, the samples are also split into batches so that only a
few batches are on scratch at a time (run_batches.sh). Batch 1 is the smallest
sample, as a pilot. The other batches are filled in sample order up to
--batch_tb of downloads. The batch files in the repo are the ones that ran.

Writes, under data/ont/:
  hprc_metadata/                 the HPRC_METADATA data tables as downloaded
  stages/stage4.*                merge_plan, manifest and download_manifest for all of stage 4
  stages/stage4/bNN.*            the same three files for each batch (--write_batches)

  build_ont_stage4.py [--refresh] [--write_batches] [--batch_tb 1.0]
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
SUBMISSIONS = [
    "UCSC_HPRC_ONT_Year1_R10",
    "UCSC_HPRC_PLUS_nanopore_R10",
    "UCSC_HPRC_nanopore_Year2_R10",
    "UCSC_HPRC_nanopore_Year3_R10",
    "UCSC_HPRC_nanopore_Year5",
]


def fetch(refresh):
    os.makedirs(META, exist_ok=True)
    paths = []
    for sub in SUBMISSIONS:
        name = sub + "_data_table.csv"
        p = os.path.join(META, name)
        if refresh or not os.path.exists(p):
            urllib.request.urlretrieve(RAW + sub + "/" + name, p)
        paths.append(p)
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


def candidates(r):
    p = r["path"]
    out = [p]
    for a, b in (("/working/HPRC/", "/working/HPRC_PLUS/"), ("/working/HPRC_PLUS/", "/working/HPRC/")):
        if a in p:
            out.append(p.replace(a, b))
    return [u for u in out if "/working/" in u]


def resolve(r):
    for uri in candidates(r):
        b = head(uri)
        if b:
            return uri, b
    return None, None


def cohort(uri):
    return uri.split("/")[4]


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
    ap.add_argument("--write_batches", action="store_true", help="also write stages/stage4/bNN.*")
    ap.add_argument("--batch_tb", type=float, default=1.0, help="downloads per batch, TB (default 1.0)")
    a = ap.parse_args(argv)

    rows = [r for p in fetch(a.refresh) for r in read_csv(p)]
    want = lai_samples() - typed_samples()
    found = {}
    for r in rows:
        if "r10.4.1" in r["basecaller_model"] and r["sample_id"] in want and r["filename"].endswith(".bam"):
            found[r["filename"]] = r

    # Released runs only, resolved before the depth rule so that a missing run is not planned.
    by = {}
    with ThreadPoolExecutor(16) as pool:
        for r, (uri, b) in zip(found.values(), pool.map(resolve, found.values())):
            if uri:
                r["uri"], r["bytes"] = uri, b
                by.setdefault(r["sample_id"], []).append(r)

    entries = []
    for s in sorted(by):
        if len({r["basecaller_model"] for r in by[s]}) > 1:
            sys.exit("%s mixes basecaller models; plan it by hand" % s)
        runs, note = plan_sample(by[s])
        entries.append({"id": s + "__R1041", "sample": s, "runs": runs, "note": note})
    runs = [r for e in entries for r in e["runs"]]

    write(os.path.join(OUT, "stages", "stage4."), entries)
    bs = batches(entries, a.batch_tb) if a.write_batches else []
    for i, b in enumerate(bs, 1):
        write(os.path.join(OUT, "stages", "stage4", "b%02d." % i), b)

    print("depth rule %g/%g/%gx; Lai samples without released R10.4.1 runs: %s" % (
        MIN_DEPTH, TARGET_DEPTH, MAX_DEPTH, " ".join(sorted(want - set(by))) or "none"))
    for i, b in enumerate(bs, 1):
        tb = sum(r["bytes"] for e in b for r in e["runs"]) / 1e12
        print("b%02d: %2d input BAMs, %2d runs, %.2f TB  %s" % (
            i, len(b), sum(len(e["runs"]) for e in b), tb, " ".join(e["sample"] for e in b)))
    print("stage 4: %d input BAMs, %d runs, %.2f TB" % (len(entries), len(runs), sum(r["bytes"] for r in runs) / 1e12))


if __name__ == "__main__":
    main(sys.argv[1:])
