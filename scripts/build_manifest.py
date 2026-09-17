#!/usr/bin/env python3
"""Write data/manifest.tsv: every input file behind every typed BAM.

Joins the plan (data/merge_plan.tsv, written once by select_inputs.py) to the
S3 object each file came from, so the benchmark can be reproduced from a fixed
list rather than by re-running the selection rules against a bucket that keeps
changing. The S3 listing is fetched per sample over HTTPS (no credentials).

    python3 scripts/build_manifest.py

Columns: input_bam, sample, cohort, library, instrument, filename, s3_uri, bytes,
cov_x. `cohort` is the S3 collection the file sits under (HPRC or HPRC_PLUS),
i.e. who sequenced it. `instrument` is read from the movie name: m84 = Revio,
m64 and m54 = Sequel II / IIe. An input BAM that mixes both is "Revio + Sequel II".
"""

import csv
import os
import re
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLAN = os.path.join(REPO_ROOT, "data", "merge_plan.tsv")
OUT = os.path.join(REPO_ROOT, "data", "manifest.tsv")

S3 = "https://human-pangenomics.s3.amazonaws.com/?list-type=2&prefix=working/%s/%s/raw_data/PacBio_HiFi/"
COHORTS = ("HPRC", "HPRC_PLUS")
KEY = re.compile(r"<Key>([^<]+\.bam)</Key>(?:(?!</Contents>).)*?<Size>(\d+)</Size>")


def instrument(filename):
    m = re.match(r"(?:.*_)?(m\d{2})", filename)
    return {"m84": "Revio", "m64": "Sequel II", "m54": "Sequel II"}.get(m.group(1) if m else "", "unknown")


def listing(samples):
    """{(sample, filename): (cohort, uri, bytes)} from the bucket."""
    seen = {}
    for s in sorted(samples):
        for cohort in COHORTS:
            xml = subprocess.run(["curl", "-s", S3 % (cohort, s)], capture_output=True, text=True).stdout
            for key, size in KEY.findall(xml):
                fname = key.split("/")[-1]
                if (s, fname) in seen:
                    continue  # first cohort listed wins, as the download did
                seen[(s, fname)] = (cohort, "s3://human-pangenomics/" + key, size)
        if not any(k[0] == s for k in seen):
            print("WARNING: no S3 objects listed for %s" % s, file=sys.stderr)
    return seen


def main():
    plan = list(csv.DictReader(open(PLAN), delimiter="\t"))
    samples = {r["sample"] for r in plan}
    s3 = listing(samples)

    by_entry = {}
    for r in plan:
        by_entry.setdefault(r["input_bam"], []).append(r)

    rows, missing = [], []
    for entry in sorted(by_entry):
        files = by_entry[entry]
        insts = sorted({instrument(r["filename"]) for r in files})
        inst = " + ".join(insts) if len(insts) > 1 else insts[0]
        for r in files:
            hit = s3.get((r["sample"], r["filename"]))
            if hit is None:
                missing.append((r["sample"], r["filename"]))
                cohort, uri, size = "", "", ""
            else:
                cohort, uri, size = hit
            rows.append([entry, r["sample"], cohort, r["prep"], inst, r["filename"], uri, size, r["cov_x"]])

    # entries outside the plan (the HG002 GIAB and HPRC_PLUS single cells) are
    # kept from the existing manifest
    extra = []
    if os.path.exists(OUT):
        extra = [r for r in csv.reader(open(OUT), delimiter="\t") if r and r[0] not in by_entry and r[0] != "input_bam"]
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["input_bam", "sample", "cohort", "library", "instrument", "filename", "s3_uri", "bytes", "cov_x"])
        w.writerows(rows)
        w.writerows(extra)
    print("wrote %s: %d files, %d entries, %d samples" % (OUT, len(rows), len(by_entry), len(samples)))
    if missing:
        print("no S3 object found for %d files:" % len(missing))
        for s, f in missing:
            print("  %s\t%s" % (s, f))


if __name__ == "__main__":
    main()
