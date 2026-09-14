#!/usr/bin/env python3
"""Grade hla_resolve WGS calls against the 1000 Genomes Project HLA panel.

Truth:  data/kgp_answer_key.tsv
        Converted from 20181129_HLA_types_full_1000_Genomes_Project_panel.txt
        (Abi-Rached et al. 2018, PolyPheMe on exome data, IMGT 3.28). Five
        genes (A, B, C, DQB1, DRB1) at two-field resolution, so only 1-field
        and 2-field concordance are evaluated. DPA1, DPB1 and DQA1 have no
        truth and are excluded from concordance; they still count toward the
        call rate, which never consults the key.

        Notation in the source and how the key records it:
          blank / "None"   no truth            -> empty cell (excluded)
          02:02*           exome-limited call  -> DQB1*02:02, asterisk dropped;
                                                  listed in data/kgp_star_alleles.tsv
          03:01/04         ambiguous 2nd field -> A*03:01/A*03:04 (any option matches)
          49:01 50:01      ambiguous 1st field -> B*49:01/B*50:01

Calls:  data/audit_1kgp.txt
        Raw hla_resolve output collected on the cluster by audit_jsheh.py
        (allele_output.csv per sample, coverage from the run log). This script
        writes the calls for every sample that has 1kGP truth to
        results/kgp_calls.csv and grades that file. A sample whose run produced
        no allele_output.csv is kept with zero calls so it counts against the
        call rate; it cannot affect concordance.

Grading rules, outputs and the self-test are those of test_against_truth.py;
this script only swaps the key, the calls file and the output paths.

After the standard report it prints concordance again restricted to genes
whose ARS (peptide-binding domain) depth reached MIN_GENE_DEPTH, read from
BLOCK3 of the audit. Call rate is unaffected. Gene x sample combos dropped by
the gate are listed so the two tables can be reconciled.
"""

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import test_against_truth as T  # noqa: E402

REPO_ROOT = os.path.dirname(HERE)
AUDIT_PATH = os.path.join(REPO_ROOT, "data", "audit_1kgp.txt")
KEY_PATH = os.path.join(REPO_ROOT, "data", "kgp_answer_key.tsv")
CALLS_PATH = os.path.join(REPO_ROOT, "results", "kgp_calls.csv")
MISSES_PATH = os.path.join(REPO_ROOT, "results", "kgp_misses.txt")
PER_SAMPLE_PATH = os.path.join(REPO_ROOT, "results", "kgp_concordance_by_sample.txt")

COLS = ["HLA-%s_%d" % (g, k) for g in T.GENES for k in (1, 2)]

# Per-gene ARS depth a gene must reach for its calls to enter the depth-gated
# concordance table. hla_resolve's own gate is 8x; 10x is the reporting cut.
MIN_GENE_DEPTH = 10.0


def read_block(path, name):
    """Rows of one '=== NAME ===' block of the audit file, keyed by sample."""
    lines = open(path).read().split("\n")
    start = lines.index("=== %s ===" % name)
    header = lines[start + 1].split("\t")
    rows = {}
    for line in lines[start + 2:]:
        if line.startswith("=== END"):
            break
        fields = line.split("\t")
        if len(fields) == len(header):
            rows[fields[0]] = dict(zip(header, fields))
    return rows


def write_calls(audit_path, key_samples, out_path):
    """Build the results CSV that test_against_truth.load_calls reads.

    Only samples with 1kGP truth are written. A run with no allele_output.csv
    appears as a row of no-calls with coverage NA.
    """
    alleles = read_block(audit_path, "BLOCK1 ALLELE CALLS")
    coverage = read_block(audit_path, "BLOCK2 COVERAGE SUMMARY")

    samples = sorted(s for s in alleles if s in key_samples)
    failed = []
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "HLA_coverage_x"] + COLS)
        for s in samples:
            depth = coverage.get(s, {}).get("mean_gene_depth", "").strip()
            calls = [alleles[s][c] for c in COLS]
            if all(c == "NA" for c in calls):
                failed.append(s)
                depth = "NA"
            w.writerow([s, depth] + calls)
    return samples, failed


def gene_depths(audit_path):
    """{sample: {gene: ARS depth}} from BLOCK3; genes without a value are absent."""
    out = {}
    for s, row in read_block(audit_path, "BLOCK3 PER-GENE ARS DEPTH").items():
        out[s] = {g: float(v) for g, v in row.items()
                  if g != "sample" and v.strip()}
    return out


def depth_gated_report(depths):
    """Concordance over genes at or above MIN_GENE_DEPTH only."""
    key = T.load_key(T.KEY_PATH)
    calls, _ = T.load_calls(T.CALLS_PATH)
    per_gene = {n: {g: [0, 0] for g in T.GENES} for n in T.DEPTHS}
    dropped = []
    for s in sorted(calls):
        for g in T.GENES:
            if not any(T.is_called(a) for a in calls[s][g]):
                continue                       # no-call: call-rate, not here
            d = depths.get(s, {}).get(g)
            if d is None or d < MIN_GENE_DEPTH:
                dropped.append((s, g, d))
                continue
            for n in T.DEPTHS:
                h, dd = T.grade_gene(key[s][g], calls[s][g], n)
                per_gene[n][g][0] += h
                per_gene[n][g][1] += dd

    print("\nCONCORDANCE, genes with ARS depth >= %gx only" % MIN_GENE_DEPTH)
    print("  %d called gene x sample combos below the gate excluded" % len(dropped))
    print("\n      " + "".join(g.rjust(9) for g in T.GENES) + "   " + "TOTAL".rjust(18))
    for n in T.DEPTHS:
        cells = "".join(("%d/%d" % tuple(per_gene[n][g])).rjust(9) for g in T.GENES)
        hits = sum(v[0] for v in per_gene[n].values())
        denom = sum(v[1] for v in per_gene[n].values())
        total = "%d/%d = %.2f%%" % (hits, denom, 100.0 * hits / denom) if denom else "n/a"
        print("  %df  %s   %s" % (n, cells, total.rjust(18)))
    if dropped:
        print("\n  excluded by the gate:")
        for s, g, d in dropped:
            print("    %-8s HLA-%-5s %s" % (s, g, "no depth" if d is None else "%.1fx" % d))


def main():
    key_samples = {row.split("\t")[0] for row in open(KEY_PATH) if row.strip()}
    samples, failed = write_calls(AUDIT_PATH, key_samples, CALLS_PATH)

    print("1kGP validation")
    print("  audit: %s" % AUDIT_PATH)
    print("  %d samples in the audit have 1kGP truth -> %s" % (len(samples), CALLS_PATH))
    if failed:
        print("  %d run(s) produced no allele_output.csv (kept as no-calls): %s"
              % (len(failed), " ".join(failed)))

    T.EXCLUDED_ALLELES = set()      # the Lai exclusion does not apply here
    T.KEY_PATH = KEY_PATH
    T.CALLS_PATH = CALLS_PATH
    T.MISSES_PATH = MISSES_PATH
    T.PER_SAMPLE_PATH = PER_SAMPLE_PATH
    T.DEPTHS = [1, 2]
    T.main()

    depth_gated_report(gene_depths(AUDIT_PATH))


if __name__ == "__main__":
    main()
