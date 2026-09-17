#!/usr/bin/env python3
"""
Grade 1000G-panel samples at 2-field on the 5 panel genes (A,B,C,DQB1,DRB1).
Panel ambiguity shorthand "18:01/02" means {18:01,18:02}; "24:02" means {24:02}.

Usage: python grade_1kg.py <results_1kg_dir> <panel.txt> <output.csv>
"""
import os, sys, csv

GENES = ["A", "B", "C", "DQB1", "DRB1"]
COLS  = {"A": (3, 4), "B": (5, 6), "C": (7, 8), "DQB1": (9, 10), "DRB1": (11, 12)}
ALLELE_COLS = ["HLA-%s_%d" % (g, i) for g in GENES for i in (1, 2)]

def two(a):                                 # -> first 2 fields, or None
    a = a.split("*")[-1]
    f = a.split(":")
    return ":".join(f[:2]) if len(f) >= 2 else None

def expand(cell):                           # truth cell -> set of 2-field options
    cell = cell.strip()
    if not cell:
        return set()
    p = cell.split("/")
    f1 = p[0].split(":")[0]
    return {two(x) for x in [p[0]] + [f1 + ":" + x for x in p[1:]]} - {None}

def load_panel(path):
    truth = {}
    for i, line in enumerate(open(path)):
        c = line.rstrip("\n").split("\t")
        if i and len(c) >= 13:
            truth[c[2].strip()] = {g: [expand(c[a]), expand(c[b])] for g, (a, b) in COLS.items()}
    return truth

def coverage(log):                          # avg of the 8 'HLA-X gene' depths
    if not os.path.isfile(log):
        return ""
    d = [float(p[2]) for p in (l.split() for l in open(log))
         if len(p) >= 3 and p[0].startswith("HLA-") and p[1] == "gene"]
    return str(round(sum(d) / len(d), 1)) if d else ""

def grade_gene(calls, truth_slots):         # -> (hits, gradable) at 2-field
    calls = [c for c in (two(x) for x in calls) if c]
    slots = [s for s in truth_slots if s]   # only slots with a truth allele count
    hits = 0
    left = list(calls)
    for s in slots:
        for c in list(left):
            if c in s:
                hits += 1
                left.remove(c)
                break
    return hits, len(slots)

def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python grade_1kg.py <results_1kg_dir> <panel.txt> <output.csv>")
    results_dir, panel, out = sys.argv[1:4]
    truth = load_panel(panel)
    rows, no_truth = [], []

    for s in sorted(os.listdir(results_dir)):
        apath = os.path.join(results_dir, s, s, "hla_typing_results", "allele_output.csv")
        if not os.path.isfile(apath):
            continue
        if s not in truth:
            no_truth.append(s); continue
        calls = next(csv.DictReader(open(apath)), None)
        if calls is None:
            continue
        typed = hits = den = 0
        for g in GENES:
            pair = [calls.get("HLA-%s_1" % g, ""), calls.get("HLA-%s_2" % g, "")]
            if any(x.strip() for x in pair):
                typed += 1
            h, d = grade_gene(pair, truth[s][g])
            hits += h; den += d
        row = {"sample": s, "dataset": "1KG", "HLA_coverage_x":
               coverage(os.path.join(results_dir, s, s + ".hla_resolve.log")),
               "genes_typed": "%d/5" % typed, "concordance_2field": "%d/%d" % (hits, den)}
        row.update({c: calls.get(c, "") for c in ALLELE_COLS})
        rows.append(row)

    header = ["sample", "dataset", "HLA_coverage_x", "genes_typed", "concordance_2field"] + ALLELE_COLS
    w = csv.DictWriter(open(out, "w", newline=""), fieldnames=header)
    w.writeheader(); w.writerows(rows)
    print("Wrote %d samples to %s" % (len(rows), out))
    if no_truth:
        print("Not in panel (%d): %s" % (len(no_truth), ", ".join(no_truth)))

if __name__ == "__main__":
    main()
