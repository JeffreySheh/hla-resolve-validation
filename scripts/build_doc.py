#!/usr/bin/env python3
"""Write docs/wgs_validation.md for the hla_resolve repository from the graded run.

    python3 scripts/build_doc.py [output.md] [--version 0.9.8] [--imgt 3.64.0]

Reads   results/calls.tsv     one row per entry (grade.py)
        results/misses.tsv    discordant gene x entry (grade.py)
        data/manifest.tsv     input files per entry with S3 URI (build_manifest.py)
Writes  the markdown document (default: ../hla_resolve/docs/wgs_validation.md)

Every number in the document is recomputed here from calls.tsv, so the text
cannot drift from the tables.
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import test_against_truth as T  # noqa: E402

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CALLS = os.path.join(REPO_ROOT, "results", "calls.tsv")
MISSES = os.path.join(REPO_ROOT, "results", "misses.tsv")
MANIFEST = os.path.join(REPO_ROOT, "data", "manifest.tsv")
DEFAULT_OUT = os.path.join(os.path.dirname(REPO_ROOT), "hla_resolve", "docs", "wgs_validation.md")

LAI_URL = "https://doi.org/10.1016/j.csbj.2024.03.030"
KGP_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt"
GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
KGP_GENES = ["A", "B", "C", "DQB1", "DRB1"]
FIELDS = ["1f", "2f", "3f", "4f"]
# One clause on how the panel was read; set to "" to drop it.
PANEL_NOTE = "Exome-based (IPD-IMGT/HLA 3.28)"
# Discordant-allele table: fields to list (empty tuple = no table).
MISS_FIELDS = (1, 2, 3)


def frac(s):
    """'x/y' -> (x, y); '' -> (0, 0)."""
    if not s:
        return 0, 0
    a, b = s.split("/")
    return int(a), int(b)


def pct(h, d):
    return "%.1f%%" % (100.0 * h / d) if d else "–"


def sum_field(rows, f):
    h = d = 0
    for r in rows:
        a, b = frac(r[f])
        h += a
        d += b
    return h, d


def stats(xs):
    xs = sorted(xs)
    n = len(xs)
    med = xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2
    return sum(xs) / n, med, xs[0], xs[-1]


def truth_link(t):
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default=DEFAULT_OUT)
    ap.add_argument("--version", default="0.9.8")
    ap.add_argument("--imgt", default="3.64.0")
    args = ap.parse_args()

    calls = list(csv.DictReader(open(CALLS), delimiter="\t"))
    keys = {"Lai": T.load_key(os.path.join(REPO_ROOT, "data", "lai_answer_key.tsv")),
            "1kGP": T.load_key(os.path.join(REPO_ROOT, "data", "kgp_answer_key.tsv"))}
    excluded = {(s_, g, a) for s_, g, a in T.EXCLUDED_ALLELES}
    misses = list(csv.DictReader(open(MISSES), delimiter="\t"))
    manifest = defaultdict(list)
    for r in csv.DictReader(open(MANIFEST), delimiter="\t"):
        manifest[r["entry"]].append(r)

    # group entries by sample, samples in order of first appearance (sorted)
    by_sample = OrderedDict()
    for r in sorted(calls, key=lambda r: (r["sample"], r["entry"])):
        by_sample.setdefault(r["sample"], []).append(r)
    lai = [r for r in calls if r["truth"] == "Lai"]
    kgp = [r for r in calls if r["truth"] == "1kGP"]
    n_lai_s = len({r["sample"] for r in lai})
    n_kgp_s = len({r["sample"] for r in kgp})
    cohorts = sorted({r["cohort"] for r in calls})
    depths = [float(r["mean_gene_depth"]) for r in calls if r["mean_gene_depth"]]
    d_mean, d_med, d_min, d_max = stats(depths)

    L = []
    w = L.append
    w("# WGS Validation")
    w("")
    w("HLA-Resolve v%s (IPD-IMGT/HLA release %s) was run on public whole-genome PacBio reads for **%d sequencing libraries from %d samples**. "
      "All libraries were sequenced by the Human Pangenome Reference Consortium (HPRC; cohorts %s) except the two Genome in a Bottle (GIAB) HG002 libraries. "
      "Every library was typed independently. Libraries sequenced across multiple SMRT cells were combined to a target genome-wide coverage depth of 30×. "
      "Samples with two libraries appear twice as independent measurements."
      % (args.version, args.imgt, len(calls), len(by_sample), ", ".join(c for c in cohorts if c != "GIAB")))
    w("")
    w("Allele calls were compared to two reference sets:")
    w("")
    w("- **Lai et al. 2024** ([%s](%s); [Supplementary File 6](Lai_Supplementary-6.xlsx))" % (LAI_URL.split("//")[1], LAI_URL))
    w("  - %d HPRC release 1 samples, %d libraries" % (n_lai_s, len(lai)))
    w("  - HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1")
    w("  - Resolution: up through four-field (full) resolution")
    w("- **1000 Genomes Project HLA panel** ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](%s))" % KGP_URL)
    w("  - %d samples, %d libraries" % (n_kgp_s, len(kgp)))
    w("  - HLA-A, -B, -C, -DQB1, -DRB1")
    w("  - Resolution: two-field resolution")
    if PANEL_NOTE:
        w("  - %s" % PANEL_NOTE)
    w("")
    w("The following command was used for each sample:")
    w("")
    w("```bash")
    w("hla_resolve --input_file <LIBRARY_uBAM> --sample_name <SAMPLE> \\")
    w("            --platform pacbio --scheme WGS \\")
    w("            --output_dir <OUTPUT_DIR> --threads <THREADS>")
    w("```")
    w("")
    w("## Scoring")
    w("")
    w("Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution. "
      "A call resolved to fewer fields than the reference is considered discordant at that level of resolution. "
      "By default, genes do not type if their antigen recognition site coverage depth sits below 8x, so call rate was measured as the proportion of total possible allele calls emitted by HLA-Resolve.")
    w("")

    # ---- summary tables ---------------------------------------------------
    for name, rows, genes, fields in (("Lai et al. reference (%d libraries, %d samples)" % (len(lai), n_lai_s), lai, GENES, FIELDS),
                                      ("1000 Genomes panel (%d libraries, %d samples)" % (len(kgp), n_kgp_s), kgp, KGP_GENES, FIELDS[:2])):
        w("## %s" % name)
        w("")
        called = sum(1 for r in rows for g in genes for i in (1, 2) if r["HLA-%s_%d" % (g, i)].strip())
        expected = 2 * len(genes) * len(rows)
        w("| Metric | Value |")
        w("|---|---:|")
        w("| Total possible allele calls | %d |" % expected)
        w("| Alleles called | %d |" % called)
        w("| Call rate | %s |" % pct(called, expected))
        w("")
        # reference alleles that cannot be evaluated: absent from the key, or excluded
        tname = "Lai" if rows is lai else "1kGP"
        absent = defaultdict(int)
        dropped = []
        for r in rows:
            k = keys[tname][r["sample"]]
            for g in genes:
                for a in k[g]:
                    if not T.is_called(a):
                        absent[g] += 1
                for s_, g_, a_ in excluded:
                    if s_ == r["sample"] and g_ == g:
                        dropped.append("%s %s" % (s_, a_))
        w("| Resolution | Concordance | Concordant alleles | Alleles evaluated |")
        w("|---|---:|---:|---:|")
        for f in fields:
            h, d = sum_field(rows, f)
            w("| %s-field | %s | %d | %d |" % (f[0], pct(h, d), h, d))
        w("")
        notes = []
        if dropped:
            notes.append("%d uncertain reference allele%s excluded (%s)" % (len(set(dropped)), "s were" if len(set(dropped)) > 1 else " was", ", ".join(sorted(set(dropped)))))
        if absent:
            notes.append("the reference has missing " + " and ".join(
                "HLA-%s calls for %d samples (%d alleles)" % (g, n // 2, n) for g, n in sorted(absent.items())))
        if tname == "Lai":
            notes.append("not all reference alleles have third and fourth field values")
        if notes:
            w("Alleles evaluated is lower than alleles called because " + ", and ".join(notes) + ".")
            w("")
    both = [r for r in calls if r["truth"] in ("Lai", "1kGP")]
    w("## All libraries")
    w("")
    w("Both reference sets pooled, %d libraries. The number of alleles evaluated falls with resolution because the 1000 Genomes panel stops at two fields and covers five genes, and because a reference allele is evaluated only at the fields it specifies." % len(both))
    w("")
    called = sum(1 for r in both for g in (GENES if r["truth"] == "Lai" else KGP_GENES) for i in (1, 2) if r["HLA-%s_%d" % (g, i)].strip())
    expected = sum(2 * len(GENES if r["truth"] == "Lai" else KGP_GENES) for r in both)
    w("| Metric | Value |")
    w("|---|---:|")
    w("| Total possible allele calls | %d |" % expected)
    w("| Alleles called | %d |" % called)
    w("| Call rate | %s |" % pct(called, expected))
    w("")
    w("| Resolution | Concordance | Concordant alleles | Alleles evaluated |")
    w("|---|---:|---:|---:|")
    for f in FIELDS:
        h, d = sum_field(both, f)
        w("| %s-field | %s | %d | %d |" % (f[0], pct(h, d), h, d))
    w("")

    # ---- per-library table ------------------------------------------------
    w("## Results by library")
    w("")
    w("One row per sequencing library, grouped by sample.")
    w("")
    w("- The 3- and 4-field columns are blank where the reference is two-field.")
    w("- HLA coverage is the mean depth across the eight classical HLA genes from the HLA-Resolve run log (all %d libraries: mean %.1f×, median %.1f×, range %.1f–%.1f×)." % (len(calls), d_mean, d_med, d_min, d_max))
    w("")
    w("| Sample | Library | Cohort | Instrument | Reference | HLA&nbsp;coverage | Genes&nbsp;typed | Alleles&nbsp;called | 1-field | 2-field | 3-field | 4-field |")
    w("|---|---|---|---|---|---:|---|---|---|---|---|---|")
    for s, rows in by_sample.items():
        for r in rows:
            f3 = r["3f"] if r["truth"] == "Lai" else ""
            f4 = r["4f"] if r["truth"] == "Lai" else ""
            w("| %s | %s | %s | %s | %s | %s× | %s | %s | %s | %s | %s | %s |" % (
                s, r["library"], r["cohort"], r["instrument"], truth_link(r["truth"]),
                "%.1f" % float(r["mean_gene_depth"]), r["genes_typed"], r["alleles_called"], r["1f"], r["2f"], f3, f4))
    w("")

    # ---- misses -------------------------------------------------------------
    shown = [m for m in misses if int(m["field"]) in MISS_FIELDS]
    if MISS_FIELDS and shown:
        w("## Discordant alleles")
        w("")
        w("Gene × library combinations discordant at %s, with the shallowest discordant field, the gene's mean ARS depth, the calls and the reference alleles."
          % (" or ".join("%d" % f for f in MISS_FIELDS) + " fields" if len(MISS_FIELDS) > 1 else "%d field" % MISS_FIELDS[0]))
        w("")
        w("| Sample | Library | Gene | Field | ARS&nbsp;depth | Called | Reference | Note |")
        w("|---|---|---|---|---:|---|---|---|")
        lib_of = {r["entry"]: r["library"] for r in calls}
        for m in sorted(shown, key=lambda m: (int(m["field"]), m["sample"], m["gene"])):
            w("| %s | %s | %s | %s | %s× | %s / %s | %s / %s | %s |" % (
                m["sample"], lib_of.get(m["entry"], ""), m["gene"], m["field"], m["ARS_depth"],
                m["call_1"], m["call_2"], m["truth_1"], m["truth_2"], m["note"]))
        w("")

    # ---- input files --------------------------------------------------------
    w("## Input files")
    w("")
    w("Each library was typed from one unaligned BAM made by concatenating the files listed. HPRC files are public on the AWS Open Data registry and need no credential. Build a path by substituting the Cohort, Sample and Files values from the table into:")
    w("")
    w("```")
    w("s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/PacBio_HiFi/<FILE>")
    w("```")
    w("")
    ex_sample = next(iter(by_sample))
    ex_entry = by_sample[ex_sample][0]["entry"]
    ex_files = manifest.get(ex_entry, [])
    ex_cohort = by_sample[ex_sample][0]["cohort"]
    ex_key = ex_files[0]["s3_uri"] if ex_files and ex_files[0]["s3_uri"].startswith("s3://") else "s3://human-pangenomics/working/%s/%s/raw_data/PacBio_HiFi/%s" % (ex_cohort, ex_sample, ex_files[0]["filename"] if ex_files else "<FILE>")
    w("For example, the first file listed for %s (cohort %s) is at:" % (ex_sample, ex_cohort))
    w("")
    w("```")
    w(ex_key)
    w("```")
    w("")
    w("Files are public and need no credentials:")
    w("")
    w("```bash")
    w("aws s3 cp --no-sign-request \\")
    w("  %s ." % ex_key)
    w("```")
    w("")
    w("Some files sit one folder deeper (`primrose/`, `wMods/`). The manifest below gives the exact key.")
    w("")
    w("| Sample | Library | Cohort | Files |")
    w("|---|---|---|---|")
    for s, rows in by_sample.items():
        for r in rows:
            files = manifest.get(r["entry"], [])
            names = "<br>".join(f["s3_uri"].split("/PacBio_HiFi/")[-1] if f["s3_uri"] else f["filename"] for f in files)
            w("| %s | %s | %s | %s |" % (s, r["library"], r["cohort"], names))
    w("")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    open(args.out, "w").write("\n".join(L) + "\n")
    print("wrote %s: %d libraries, %d samples, %d discordant gene x library" % (args.out, len(calls), len(by_sample), len(misses)))


if __name__ == "__main__":
    main()
