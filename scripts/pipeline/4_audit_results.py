#!/usr/bin/env python3
"""Collect HLA-Resolve results from a typing run into one file.

    python3 scripts/pipeline/4_audit_results.py <results root> > results/audit.txt

Each directory under the root is one input BAM (<sample>__<library>), typed by
3_type.sbatch, with HLA-Resolve's own <sample>/ nested inside:

    <root>/<input_bam>/<sample>/hla_typing_results/allele_output.csv
    <root>/<input_bam>/<sample>.hla_resolve.log

Writes three tab-separated blocks to stdout: the allele calls, a coverage
summary per input BAM, and the per-gene ARS depth, all read from the run's own
output files. An input BAM with no result yet is listed with empty fields.

Stdlib only, Python 3.6+.
"""

import csv
import os
import re
import sys

GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
COLS = ["HLA-%s_%s" % (g, k) for g in GENES for k in ("1", "2")]

# coverage table in the run log: "HLA-A   gene   45.2   100.0%  99.0%  98.0%"
COV_RE = re.compile(r"^(HLA-[A-Z0-9]+)\s+(gene|CDS|ARS)\s+([0-9.]+)\s")


def norm(a):
    a = (a or "").strip()
    return a[4:] if a.startswith("HLA-") else a


def read_alleles(path):
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    return rows[0] if rows else None


def parse_log(path):
    """{gene: {"gene": x, "CDS": x, "ARS": x}} from the log's coverage table."""
    cov = {}
    if not os.path.exists(path):
        return cov
    with open(path, errors="replace") as fh:
        for line in fh:
            m = COV_RE.match(line)
            if m:
                cov.setdefault(norm(m.group(1)), {})[m.group(2)] = float(m.group(3))
    return cov


def mean(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def fmt(x, nd=2):
    return "" if x is None else ("%." + str(nd) + "f") % x


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    root = sys.argv[1]
    entries = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)))

    allele_rows, cov_rows, ars_rows = [], [], []
    for e in entries:
        sample = e.split("__", 1)[0]
        row = read_alleles(os.path.join(root, e, sample, "hla_typing_results", "allele_output.csv"))
        allele_rows.append([e] + ([norm(row.get(c)) for c in COLS] if row else [""] * 16))

        cov = parse_log(os.path.join(root, e, sample + ".hla_resolve.log"))
        ars = [cov.get(g, {}).get("ARS") for g in GENES]
        cov_rows.append([e,
                         fmt(mean([cov.get(g, {}).get("gene") for g in GENES])),
                         fmt(mean([cov.get(g, {}).get("CDS") for g in GENES])),
                         fmt(mean(ars)),
                         fmt(min(v for v in ars if v is not None) if any(v is not None for v in ars) else None)])
        ars_rows.append([e] + [fmt(v, 1) for v in ars])

    def block(title, header, rows):
        print("=== %s ===" % title)
        print("\t".join(header))
        for r in rows:
            print("\t".join(r))
        print("=== END %s ===" % title)
        print()

    print("# HLA-Resolve results audit")
    print("# root: %s" % root)
    print("# entries: %d" % len(entries))
    print()
    block("BLOCK1 ALLELE CALLS", ["input_bam"] + COLS, allele_rows)
    block("BLOCK2 COVERAGE SUMMARY", ["input_bam", "mean_gene_depth", "mean_CDS_depth", "mean_ARS_depth", "min_ARS_depth"], cov_rows)
    block("BLOCK3 PER-GENE ARS DEPTH", ["input_bam"] + GENES, ars_rows)


if __name__ == "__main__":
    main()
