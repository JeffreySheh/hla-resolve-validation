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
from collections import OrderedDict, defaultdict

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
    return {"Lai": "[Lai][lai]", "1kGP": "[1kGP][kgp]"}.get(t, t)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", nargs="?", default=DEFAULT_OUT)
    ap.add_argument("--version", default="0.9.8")
    ap.add_argument("--imgt", default="3.64.0")
    args = ap.parse_args()

    calls = list(csv.DictReader(open(CALLS), delimiter="\t"))
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
      "Every library was typed on its own: a sample with two libraries appears twice, as two independent measurements."
      % (args.version, args.imgt, len(calls), len(by_sample), ", ".join(c for c in cohorts if c != "GIAB")))
    w("")
    w("Allele calls were compared to two reference sets:")
    w("")
    w("- **Lai et al. 2024** ([%s](%s)) for %d HPRC year-1 samples (%d libraries): eight classical genes, up to four fields." % (LAI_URL.split("//")[1], LAI_URL, n_lai_s, len(lai)))
    w("- **1000 Genomes Project HLA panel** ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](%s)) for %d samples (%d libraries): HLA-A, -B, -C, -DQB1 and -DRB1 at two fields. "
      "The panel is exome-based (PolyPheMe, IPD-IMGT/HLA 3.28); alleles it lists with an alternative (`A*03:01/A*03:04`) were accepted as equivalent, and its exome-only flag (`*`) was kept." % (KGP_URL, n_kgp_s, len(kgp)))
    w("")
    w("## Input selection")
    w("")
    w("Input files were chosen by a fixed set of rules (`select_inputs.py` in the validation repository) and are listed under [Input files](#input-files):")
    w("")
    w("1. One library per entry, identified from the HPRC file index (`library_ID`) or the BAM read-group `LB` tag. A sample sequenced from two libraries gives two entries.")
    w("2. Excluded: files under 1 GB (barcode cross-talk), files whose read group names another sample, non-WGS preparations (Fiber-seq, AMPure, dilution tests), and second encodings of the same SMRT cell.")
    w("3. SMRT cells were added in run order to an estimated 30× genome-wide depth (HPRC index coverage, or file size × a per-chemistry rate). A cell that would take the library past 40× was left out once it was above 20×. A library under 20× on its own was typed only when it was the sample's only library.")
    w("4. Cells were concatenated with `samtools cat`; no duplicate removal.")
    w("")
    w("```bash")
    w("hla_resolve --input_file <LIBRARY_uBAM> --sample_name <SAMPLE> \\")
    w("            --platform pacbio --scheme WGS \\")
    w("            --output_dir <OUTPUT_DIR> --threads 24")
    w("```")
    w("")
    w("## Grading")
    w("")
    w("Call rate counts every allele HLA-Resolve emitted against two per gene. Concordance is measured among alleles called: a gene the tool did not type (antigen recognition site under 8× mean depth) counts against call rate, not concordance. "
      "An allele was evaluated at a field resolution only if the reference specified it to that resolution; a call resolved to fewer fields than the reference is discordant. The two alleles of a gene were paired with the reference by the assignment that maximizes agreement.")
    w("")
    w("Coverage is the mean depth across the eight classical HLA genes from the HLA-Resolve run log. Across all %d libraries: mean %.1f×, median %.1f×, range %.1f–%.1f×." % (len(calls), d_mean, d_med, d_min, d_max))
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
        w("| Resolution | Concordance | Concordant alleles | Alleles evaluated |")
        w("|---|---:|---:|---:|")
        for f in fields:
            h, d = sum_field(rows, f)
            w("| %s-field | %s | %d | %d |" % (f[0], pct(h, d), h, d))
        w("")
    both = [r for r in calls if r["truth"] in ("Lai", "1kGP")]
    w("## All libraries, HLA-A, -B, -C, -DQB1, -DRB1")
    w("")
    w("The five genes both references cover, over all %d libraries, at one and two fields." % len(both))
    w("")
    w("| Resolution | Concordance | Concordant alleles | Alleles evaluated |")
    w("|---|---:|---:|---:|")
    for f in (1, 2):
        h, d = sum_field(both, "ABCDQB1DRB1_%df" % f)
        w("| %d-field | %s | %d | %d |" % (f, pct(h, d), h, d))
    w("")

    # ---- per-library table ------------------------------------------------
    w("## Results by library")
    w("")
    w("One row per sequencing library, grouped by sample. 3- and 4-field columns are blank where the reference is two-field.")
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
    w("## Discordant alleles")
    w("")
    w("Every gene × library discordant at any field, with the shallowest discordant field, the gene's mean ARS depth, the calls and the reference alleles.")
    w("")
    w("| Sample | Library | Gene | Field | ARS&nbsp;depth | Called | Reference | Note |")
    w("|---|---|---|---|---:|---|---|---|")
    lib_of = {r["entry"]: r["library"] for r in calls}
    for m in sorted(misses, key=lambda m: (int(m["field"]), m["sample"], m["gene"])):
        w("| %s | %s | %s | %s | %s× | %s / %s | %s / %s | %s |" % (
            m["sample"], lib_of.get(m["entry"], ""), m["gene"], m["field"], m["ARS_depth"],
            m["call_1"], m["call_2"], m["truth_1"], m["truth_2"], m["note"]))
    w("")

    # ---- input files --------------------------------------------------------
    w("## Input files")
    w("")
    w("Each library was typed from one unaligned BAM made by concatenating the files listed. HPRC files are public on the AWS Open Data registry and need no credentials; the S3 URI for each file is in `data/manifest.tsv` of the validation repository. Build a path from the table with:")
    w("")
    w("```")
    w("s3://human-pangenomics/working/<COHORT>/<SAMPLE>/raw_data/PacBio_HiFi/<FILE>")
    w("```")
    w("")
    w("Some files sit one folder deeper (`primrose/`, `wMods/`); `aws s3 ls --recursive` or the manifest gives the exact key.")
    w("")
    w("| Sample | Library | Cohort | Files |")
    w("|---|---|---|---|")
    for s, rows in by_sample.items():
        for r in rows:
            files = manifest.get(r["entry"], [])
            names = "<br>".join(f["s3_uri"].split("/PacBio_HiFi/")[-1] if f["s3_uri"] else f["filename"] for f in files)
            w("| %s | %s | %s | %s |" % (s, r["library"], r["cohort"], names))
    w("")

    w("[lai]: %s" % LAI_URL)
    w("[kgp]: %s" % KGP_URL)
    w("")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    open(args.out, "w").write("\n".join(L) + "\n")
    print("wrote %s: %d libraries, %d samples, %d discordant gene x library" % (args.out, len(calls), len(by_sample), len(misses)))


if __name__ == "__main__":
    main()
