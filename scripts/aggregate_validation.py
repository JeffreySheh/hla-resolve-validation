#!/usr/bin/env python3
"""Aggregate the two WGS validations into one per-sample table and one summary.

Datasets
  hprc  data/hprc_calls.csv   43 samples, truth = Lai et al. (data/lai_answer_key.tsv),
                              eight genes, one- to four-field resolution.
                              40 rows are Supplementary Table 17 of the manuscript;
                              HG002, HG01258 and HG03579 were verified separately.
  1kgp  results/kgp_calls.csv 89 samples, truth = 1000 Genomes HLA panel
                              (data/kgp_answer_key.tsv), five genes, one- and
                              two-field resolution. Written by test_against_1kgp.py.
                              Instrument comes from data/kgp_input_manifest.csv
                              (build_kgp_manifest.py).

Each dataset is graded against its own truth with test_against_truth.grade_gene.
The combined row covers only what both truths resolve: one- and two-field
concordance on A, B, C, DQB1 and DRB1. Two samples occur in both datasets from
the same sequencing library (HG00733, NA19240); the combined row keeps the hprc
copy and drops the 1kgp copy so no library is counted twice.

Outputs
  results/validation_by_sample.csv   one row per sample x dataset
  results/validation_summary.txt     per-dataset and combined tables
"""

import csv
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_against_truth as T  # noqa: E402

REPO_ROOT = os.path.dirname(HERE)
D = lambda *p: os.path.join(REPO_ROOT, *p)  # noqa: E731

DATASETS = [
    # name, calls file, key file, depths graded, genes with truth
    ("hprc", D("data", "hprc_calls.csv"), D("data", "lai_answer_key.tsv"),
     [1, 2, 3, 4], T.GENES),
    ("1kgp", D("results", "kgp_calls.csv"), D("data", "kgp_answer_key.tsv"),
     [1, 2], ["A", "B", "C", "DQB1", "DRB1"]),
]
SHARED_GENES = ["A", "B", "C", "DQB1", "DRB1"]
# same library in both datasets -> keep the hprc copy in the combined row
SAME_LIBRARY = {"HG00733", "NA19240"}

KGP_MANIFEST = D("data", "kgp_input_manifest.csv")   # instrument for 1kgp rows
BY_SAMPLE = D("results", "validation_by_sample.csv")
SUMMARY = D("results", "validation_summary.txt")


def load_calls(path):
    calls, cov, inst = {}, {}, {}
    for row in csv.DictReader(open(path, newline="")):
        s = row["sample"].strip()
        calls[s] = {g: [row["HLA-%s_1" % g].strip(), row["HLA-%s_2" % g].strip()]
                    for g in T.GENES}
        cov[s] = row.get("HLA_coverage_x", "").strip()
        inst[s] = row.get("instrument", "").strip()
    return calls, cov, inst


def grade(key, calls, genes, n):
    h = d = 0
    for g in genes:
        a, b = T.grade_gene(key[g], calls[g], n)
        h += a
        d += b
    return h, d


def fmt(h, d):
    return "%d/%d = %.2f%%" % (h, d, 100.0 * h / d) if d else "n/a"


def main():
    T.selftest()
    rows = []
    tables = {}
    combined = {n: [0, 0] for n in (1, 2)}
    combined_n = 0

    manifest_inst = {r["sample"]: r["instrument"]
                     for r in csv.DictReader(open(KGP_MANIFEST, newline=""))}

    for name, calls_path, key_path, depths, genes in DATASETS:
        key = T.load_key(key_path)
        if name == "hprc":
            T.apply_exclusions(key)
        calls, cov, inst = load_calls(calls_path)
        if name == "1kgp":
            inst = {s: manifest_inst.get(s, "") for s in calls}
        per_gene = {n: defaultdict(lambda: [0, 0]) for n in depths}
        called = expected = 0
        for s in sorted(calls):
            typed = sum(1 for g in T.GENES if any(T.is_called(a) for a in calls[s][g]))
            n_called = sum(1 for g in T.GENES for a in calls[s][g] if T.is_called(a))
            called += n_called
            expected += 16
            row = {"sample": s, "dataset": name, "instrument": inst[s],
                   "HLA_coverage_x": cov[s], "genes_typed": "%d/8" % typed,
                   "alleles_called": "%d/16" % n_called}
            for n in (1, 2, 3, 4):
                if n in depths:
                    h, d = grade(key[T.key_sample(s)], calls[s], genes, n)
                    row["%df" % n] = "%d/%d" % (h, d)
                    for g in genes:
                        a, b = T.grade_gene(key[T.key_sample(s)][g], calls[s][g], n)
                        per_gene[n][g][0] += a
                        per_gene[n][g][1] += b
                else:
                    row["%df" % n] = ""
            rows.append(row)
            if not (name == "1kgp" and s in SAME_LIBRARY):
                combined_n += 1
                for n in (1, 2):
                    h, d = grade(key[T.key_sample(s)], calls[s], SHARED_GENES, n)
                    combined[n][0] += h
                    combined[n][1] += d
        tables[name] = (len(calls), called, expected, depths, genes, per_gene)

    cols = ["sample", "dataset", "instrument", "HLA_coverage_x", "genes_typed",
            "alleles_called", "1f", "2f", "3f", "4f"]
    with open(BY_SAMPLE, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    out = []
    for name, (n_samples, called, expected, depths, genes, per_gene) in tables.items():
        out.append("%s  (%d samples, %d genes with truth)" % (name.upper(), n_samples, len(genes)))
        out.append("  call rate  %s  (all eight genes, key not consulted)" % fmt(called, expected))
        out.append("  " + "".join(g.rjust(9) for g in genes) + "   " + "TOTAL".rjust(18))
        for n in depths:
            cells = "".join(("%d/%d" % tuple(per_gene[n][g])).rjust(9) for g in genes)
            h = sum(v[0] for v in per_gene[n].values())
            d = sum(v[1] for v in per_gene[n].values())
            out.append("  %df%s   %s" % (n, cells, fmt(h, d).rjust(18)))
        out.append("")
    out.append("COMBINED  (%d samples; A, B, C, DQB1, DRB1; %s counted once, hprc copy kept)"
               % (combined_n, ", ".join(sorted(SAME_LIBRARY))))
    for n in (1, 2):
        out.append("  %df  %s" % (n, fmt(*combined[n])))
    open(SUMMARY, "w").write("\n".join(out) + "\n")
    print("\n".join(out))
    print("\nPer sample: %s" % BY_SAMPLE)
    print("Summary:    %s" % SUMMARY)


if __name__ == "__main__":
    main()
