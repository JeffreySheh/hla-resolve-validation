#!/usr/bin/env python3
"""Grade every typed entry against its reference typings, both truths at once.

    python3 scripts/grade.py ~/Desktop/audit_all.txt [more audit files...]

An entry is one sample x one library (<sample>__<library>); its calls come from
the audit written by audit_results.py. The truth for an entry is chosen by
sample: the Lai et al. key (eight genes, up to four fields) when the sample is
in it, otherwise the 1000 Genomes panel (HLA-A, -B, -C, -DQB1, -DRB1, two
fields). Grading semantics are those of test_against_truth.py (unordered
pairing, denominator set by the truth's resolution, no-calls count against
call rate and not against concordance); the functions are imported from there.

Reads   data/manifest.tsv           entry -> sample, cohort, library, instrument
        data/lai_answer_key.tsv     Lai et al. reference typings
        data/kgp_answer_key.tsv     1000 Genomes panel, converted
        data/kgp_star_alleles.tsv   panel alleles flagged exome-only
Writes  results/calls.tsv           one row per entry: metadata, depth, call
                                    rate, concordance at each field, alleles
        results/misses.tsv          every discordant gene x entry with depth
        results/summary.txt         call rate and concordance tables per truth,
                                    a combined table, and a depth-gated table
"""

import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_against_truth as T  # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(REPO_ROOT, "data", "manifest.tsv")
LAI_KEY = os.path.join(REPO_ROOT, "data", "lai_answer_key.tsv")
KGP_KEY = os.path.join(REPO_ROOT, "data", "kgp_answer_key.tsv")
KGP_STAR = os.path.join(REPO_ROOT, "data", "kgp_star_alleles.tsv")
RESULTS = os.path.join(REPO_ROOT, "results")

GENES = T.GENES
TRUTHS = {
    "Lai": {"genes": GENES, "depths": [1, 2, 3, 4]},
    "1kGP": {"genes": ["A", "B", "C", "DQB1", "DRB1"], "depths": [1, 2]},
}
COMBINED_GENES = ["A", "B", "C", "DQB1", "DRB1"]
GATE_DEPTH = 10.0      # per-gene ARS depth for the gated table
DROPOUT_DEPTH = 20.0   # a miss at a gene below this is flagged as likely allelic dropout


# --------------------------------------------------------------------------
# input
# --------------------------------------------------------------------------

def read_block(path, name):
    rows, inside = [], False
    for line in open(path):
        line = line.rstrip("\n")
        if line == "=== %s ===" % name:
            inside = True
            continue
        if line == "=== END %s ===" % name:
            break
        if inside and line and not line.startswith("#"):
            rows.append(line.split("\t"))
    header, body = rows[0], rows[1:]
    return [dict(zip(header, r)) for r in body]


def load_audits(paths):
    """calls[entry][gene] = [a1, a2]; depth[entry] = coverage row; ars[entry][gene] = float."""
    calls, depth, ars = {}, {}, {}
    for p in paths:
        for r in read_block(p, "BLOCK1 ALLELE CALLS"):
            e = r["sample"]
            row = {g: [r["HLA-%s_1" % g], r["HLA-%s_2" % g]] for g in GENES}
            if not any(T.is_called(a) for g in GENES for a in row[g]):
                continue  # entry still running: audit found the directory but no calls
            calls[e] = row
        for r in read_block(p, "BLOCK2 COVERAGE SUMMARY"):
            depth[r["sample"]] = r
        for r in read_block(p, "BLOCK3 PER-GENE ARS DEPTH"):
            ars[r["sample"]] = {g: float(r[g]) if r.get(g, "") not in ("", "NA") else None for g in GENES}
    return calls, depth, ars


def load_manifest():
    meta = {}
    for r in csv.DictReader(open(MANIFEST), delimiter="\t"):
        m = meta.setdefault(r["entry"], {"sample": r["sample"], "cohort": r["cohort"],
                                          "library": r["library"], "instrument": r["instrument"],
                                          "files": []})
        m["files"].append(r["filename"])
    return meta


def load_star():
    star = set()
    if os.path.exists(KGP_STAR):
        for r in csv.DictReader(open(KGP_STAR), delimiter="\t"):
            star.add((r["sample"], r["gene"]))
    return star


# --------------------------------------------------------------------------
# grading
# --------------------------------------------------------------------------

def grade_entry(truth_row, call_row, genes, depths):
    """{(gene, n): (hits, denom)}"""
    out = {}
    for g in genes:
        for n in depths:
            out[(g, n)] = T.grade_gene(truth_row[g], call_row[g], n)
    return out


def fmt(h, d):
    return "%d/%d" % (h, d) if d else ""


def pct(h, d):
    return "%6.2f%%" % (100.0 * h / d) if d else "     --"


def table(title, genes, depths, cell):
    """cell(gene, n) -> (hits, denom). Returns the printed block as a string."""
    lines = [title, "  %6s" % "" + "".join("%9s" % g for g in genes) + "%20s" % "TOTAL"]
    for n in depths:
        h_t = d_t = 0
        row = "  %4df  " % n
        for g in genes:
            h, d = cell(g, n)
            h_t += h
            d_t += d
            row += "%9s" % fmt(h, d)
        row += "%12s = %s" % (fmt(h_t, d_t), pct(h_t, d_t).strip())
        lines.append(row)
    return "\n".join(lines)


def main(paths):
    calls, depth, ars = load_audits(paths)
    meta = load_manifest()
    lai = T.load_key(LAI_KEY)
    kgp = T.load_key(KGP_KEY)
    star = load_star()
    T.apply_exclusions(lai)

    os.makedirs(RESULTS, exist_ok=True)
    # only entries the manifest defines; a typed directory the plan later
    # dropped (a thin library superseded by a deeper one) is left out
    not_planned = sorted(e for e in calls if e not in meta)
    entries = sorted(e for e in calls if e in meta)
    truth_of, graded, skipped = {}, {}, []
    for e in entries:
        s = T.key_sample(e)
        if s in lai:
            truth_of[e] = ("Lai", lai[s])
        elif s in kgp:
            truth_of[e] = ("1kGP", kgp[s])
        else:
            skipped.append(e)
            continue
        name, key = truth_of[e]
        graded[e] = grade_entry(key, calls[e], TRUTHS[name]["genes"], TRUTHS[name]["depths"])

    # ---- calls.tsv --------------------------------------------------------
    cols = (["entry", "sample", "cohort", "library", "instrument", "truth", "n_files",
             "mean_gene_depth", "mean_ARS_depth", "min_ARS_depth"]
            + ["ARS_%s" % g for g in GENES]
            + ["genes_typed", "alleles_called", "1f", "2f", "3f", "4f", "ABCDQB1DRB1_1f", "ABCDQB1DRB1_2f"]
            + ["HLA-%s_%d" % (g, i) for g in GENES for i in (1, 2)])
    with open(os.path.join(RESULTS, "calls.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for e in entries:
            m = meta.get(e, {})
            name = truth_of.get(e, ("", None))[0]
            d = depth.get(e, {})
            a = ars.get(e, {})
            typed = sum(1 for g in GENES if any(T.is_called(x) for x in calls[e][g]))
            called = sum(1 for g in GENES for x in calls[e][g] if T.is_called(x))
            gr = graded.get(e, {})
            per_field = []
            for n in (1, 2, 3, 4):
                h = sum(v[0] for (g, k), v in gr.items() if k == n)
                dd = sum(v[1] for (g, k), v in gr.items() if k == n)
                per_field.append(fmt(h, dd))
            # the five genes both references cover, for the combined table
            for n in (1, 2):
                h = sum(gr.get((g, n), (0, 0))[0] for g in COMBINED_GENES)
                dd = sum(gr.get((g, n), (0, 0))[1] for g in COMBINED_GENES)
                per_field.append(fmt(h, dd))
            w.writerow([e, T.key_sample(e), m.get("cohort", ""), m.get("library", ""), m.get("instrument", ""),
                        name, len(m.get("files", [])),
                        d.get("mean_gene_depth", ""), d.get("mean_ARS_depth", ""), d.get("min_ARS_depth", ""),
                        *["%.1f" % a[g] if a.get(g) is not None else "" for g in GENES],
                        "%d/8" % typed, "%d/16" % called, *per_field,
                        *[calls[e][g][i] for g in GENES for i in (0, 1)]])

    # ---- misses.tsv -------------------------------------------------------
    n_miss = 0
    with open(os.path.join(RESULTS, "misses.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["entry", "sample", "truth", "gene", "field", "hits", "denom", "ARS_depth",
                    "call_1", "call_2", "truth_1", "truth_2", "note"])
        for e in entries:
            if e not in graded:
                continue
            name, key = truth_of[e]
            for g in TRUTHS[name]["genes"]:
                # report the shallowest field at which the gene first misses
                for n in TRUTHS[name]["depths"]:
                    h, d = graded[e][(g, n)]
                    if d and h < d:
                        dep = ars.get(e, {}).get(g)
                        notes = []
                        if dep is not None and dep < DROPOUT_DEPTH:
                            notes.append("ARS depth under %gx" % DROPOUT_DEPTH)
                        if name == "1kGP" and (T.key_sample(e), g) in star:
                            notes.append("panel allele flagged exome-only (*)")
                        w.writerow([e, T.key_sample(e), name, g, n, h, d,
                                    "%.1f" % dep if dep is not None else "",
                                    *calls[e][g], *key[g], "; ".join(notes)])
                        n_miss += 1
                        break

    # ---- summary.txt ------------------------------------------------------
    out = []
    out.append("hla_resolve WGS validation: %d entries graded from %d audit file(s)" % (len(graded), len(paths)))
    if skipped:
        out.append("no reference typing for %d entries (not graded): %s" % (len(skipped), " ".join(skipped)))
    if not_planned:
        out.append("typed but not in the manifest (ignored): %s" % " ".join(not_planned))
    for name in ("Lai", "1kGP"):
        es = [e for e in entries if truth_of.get(e, ("",))[0] == name]
        if not es:
            continue
        genes, depths = TRUTHS[name]["genes"], TRUTHS[name]["depths"]
        samples = {T.key_sample(e) for e in es}
        out.append("")
        out.append("== %s truth: %d entries, %d samples ==" % (name, len(es), len(samples)))
        called = sum(1 for e in es for g in genes for x in calls[e][g] if T.is_called(x))
        expected = 2 * len(genes) * len(es)
        untyped = [(e, g) for e in es for g in genes if not any(T.is_called(x) for x in calls[e][g])]
        out.append("call rate: %d/%d = %.2f%% (%d gene x entry untyped)" % (called, expected, 100.0 * called / expected, len(untyped)))
        for e, g in untyped:
            dep = ars.get(e, {}).get(g)
            out.append("  %-45s HLA-%-5s ARS %s" % (e, g, "%.1fx" % dep if dep is not None else "NA"))
        cell = lambda g, n, es=es: (sum(graded[e][(g, n)][0] for e in es), sum(graded[e][(g, n)][1] for e in es))
        out.append(table("concordance among called alleles:", genes, depths, cell))
        if name == "1kGP":
            gated = lambda g, n, es=es: (
                sum(graded[e][(g, n)][0] for e in es if (ars.get(e, {}).get(g) or 0) >= GATE_DEPTH),
                sum(graded[e][(g, n)][1] for e in es if (ars.get(e, {}).get(g) or 0) >= GATE_DEPTH))
            out.append(table("concordance, genes with ARS depth >= %gx only:" % GATE_DEPTH, genes, depths, gated))
        deps = sorted(float(depth[e]["mean_gene_depth"]) for e in es if depth.get(e, {}).get("mean_gene_depth"))
        if deps:
            out.append("mean HLA gene depth: mean %.1fx, median %.1fx, range %.1f-%.1fx" % (
                sum(deps) / len(deps), deps[len(deps) // 2], deps[0], deps[-1]))
    es = [e for e in entries if e in graded]
    out.append("")
    out.append("== combined, %d entries, HLA-A/-B/-C/-DQB1/-DRB1 at 1-2 fields ==" % len(es))
    cell = lambda g, n, es=es: (sum(graded[e][(g, n)][0] for e in es), sum(graded[e][(g, n)][1] for e in es))
    out.append(table("concordance among called alleles:", COMBINED_GENES, [1, 2], cell))
    out.append("")
    out.append("misses: %d gene x entry -> results/misses.tsv" % n_miss)
    text = "\n".join(out) + "\n"
    open(os.path.join(RESULTS, "summary.txt"), "w").write(text)
    print(text)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
