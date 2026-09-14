#!/usr/bin/env python3
"""
Audit hla_resolve raw output on the cluster.

Reads each sample's hla_typing_results/allele_output.csv plus the coverage
data, so the result is independent of however the results spreadsheet was
assembled.

Usage (on the cluster):
    python3 audit_jsheh.py                                    # default root
    python3 audit_jsheh.py <root>                             # other tree
    python3 audit_jsheh.py <root> v0.9.0                      # nested layout

The optional second argument is a per-sample subdirectory to prefer -- Jeffrey's
v0.9.0 re-run lives at <root>/<sample>/v0.9.0/hla_typing_results/ while the
older run sits at <root>/<sample>/hla_typing_results/. Samples lacking the
subdirectory fall back to the sample directory, and BLOCK4 records which tree
each sample was read from, so a half-updated tree cannot silently mislead.

Coverage comes from mosdepth's regions.bed.gz when present (authoritative, and
not truncated by re-runs) and from the log's coverage table otherwise.

Stdlib only, Python 3.6+. Writes nothing, reads only.
"""

import csv
import gzip
import os
import re
import sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/hb/groups/cornejo_training/jsheh"
SUBDIR = sys.argv[2] if len(sys.argv) > 2 else ""

GENES = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
COLS = ["HLA-%s_%s" % (g, k) for g in GENES for k in ("1", "2")]

# "HLA-A gene 45.2 100.0% 99.0% 98.0%"
COV_RE = re.compile(r"^(HLA-[A-Z0-9]+)\s+(gene|CDS|ARS)\s+([0-9.]+)\s")
INSUF_RE = re.compile(r"^Gene\s+(\S+)\s+has insufficient (.+?) coverage")
SUMMARY_NAME = ".mosdepth.summary.txt"


def norm(a):
    a = (a or "").strip()
    return a[4:] if a.startswith("HLA-") else a


def sample_dirs(root):
    out = []
    try:
        for name in sorted(os.listdir(root)):
            if os.path.isdir(os.path.join(root, name)):
                out.append(name)
    except OSError as e:
        sys.exit("ERROR: cannot list %s (%s)" % (root, e))
    return out


def resolve_dir(root, sample):
    """Pick the directory that actually holds this sample's results.

    Three layouts occur in practice and they are not distinguishable from the
    root alone, so each candidate is tested for hla_typing_results/allele_output.csv
    rather than assumed:

      <root>/<sample>/<SUBDIR>/   explicit, e.g. Jeffrey's v0.9.0 re-run
      <root>/<sample>/<sample>/  hla_resolve's own nesting -- sample_manager.py
                                 joins sample_ID onto --output_dir, so passing
                                 --output_dir <dir>/<sample> yields <sample>/<sample>
      <root>/<sample>/           flattened, e.g. a hand copy

    Returns (path, which) where `which` is recorded in BLOCK4.
    """
    base = os.path.join(root, sample)
    marker = os.path.join("hla_typing_results", "allele_output.csv")

    cands = []
    if SUBDIR:
        cands.append((os.path.join(base, SUBDIR), SUBDIR))
    cands.append((os.path.join(base, sample), "<sample>/<sample> (hla_resolve nesting)"))
    cands.append((base, "flat"))

    for path, which in cands:
        if os.path.exists(os.path.join(path, marker)):
            return path, which

    # Nothing had results; fall back so the caller still reports a useful note.
    if SUBDIR and os.path.isdir(os.path.join(base, SUBDIR)):
        return os.path.join(base, SUBDIR), "%s (no allele_output.csv)" % SUBDIR
    return base, "flat (no allele_output.csv)"


def read_alleles(sample_dir):
    path = os.path.join(sample_dir, "hla_typing_results", "allele_output.csv")
    if not os.path.exists(path):
        return None, "MISSING allele_output.csv"
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except Exception as e:
        return None, "UNREADABLE (%s)" % e
    if not rows:
        return None, "EMPTY"
    return rows[0], ("" if len(rows) == 1 else "%d rows, used first" % len(rows))


def find_log(sample_dir, sample):
    parent = os.path.dirname(sample_dir)
    for cand in (os.path.join(sample_dir, "%s.hla_resolve.log" % sample),
                 os.path.join(sample_dir, "hla_resolve.log"),
                 os.path.join(parent, "%s.hla_resolve.log" % sample)):
        if os.path.exists(cand):
            return cand
    return None


def parse_log(path):
    cov, insuf = {}, []
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                m = COV_RE.match(line)
                if m:
                    gene = norm(m.group(1))
                    cov.setdefault(gene, {})[m.group(2)] = float(m.group(3))
                    continue
                m = INSUF_RE.match(line)
                if m:
                    insuf.append("%s(%s)" % (m.group(1),
                                             m.group(2).split("[")[0].strip()))
    except Exception as e:
        return {}, ["LOG UNREADABLE: %s" % e]
    return cov, insuf


def read_regions_bed(sample_dir, sample):
    """Per-gene depths from mosdepth's regions.bed.gz. Returns (cov, lengths)."""
    path = None
    for cand in (os.path.join(sample_dir, "mosdepth", sample + ".regions.bed.gz"),
                 os.path.join(sample_dir, sample + ".regions.bed.gz")):
        if os.path.exists(cand):
            path = cand
            break
    if path is None:
        return {}, {}
    sums, lens = {}, {}
    try:
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                start, end, name, depth = int(f[1]), int(f[2]), f[3], float(f[4])
                length = end - start
                if name.endswith("_ARS"):
                    gene, kind = name[:-4], "ARS"
                elif "_CDS_" in name:
                    gene, kind = name.split("_CDS_")[0], "CDS"
                else:
                    gene, kind = name.split("_")[0], "gene"
                gene = norm(gene)
                sums.setdefault(gene, {}).setdefault(kind, 0.0)
                lens.setdefault(gene, {}).setdefault(kind, 0)
                sums[gene][kind] += depth * length     # length-weighted, so the
                lens[gene][kind] += length             # pooled CDS exons are right
    except Exception:
        return {}, {}
    cov = {}
    for gene in sums:
        cov[gene] = {}
        for kind in sums[gene]:
            L = lens[gene][kind]
            cov[gene][kind] = sums[gene][kind] / L if L else 0.0
    return cov, lens


def mosdepth_total(sample_dir, sample):
    for base in (os.path.join(sample_dir, "mosdepth", sample + SUMMARY_NAME),
                 os.path.join(sample_dir, sample + SUMMARY_NAME)):
        if not os.path.exists(base):
            continue
        try:
            with open(base) as fh:
                for line in fh:
                    f = line.split("\t")
                    if f and f[0].strip() in ("total_region", "total"):
                        return f[3].strip()
        except Exception:
            return "ERR"
    return ""


def mean(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def fmt(x, nd=2):
    return "" if x is None else ("%." + str(nd) + "f") % x


def main():
    samples = sample_dirs(ROOT)
    print("# hla_resolve raw-output audit")
    print("# root: %s" % ROOT)
    print("# per-sample subdirectory preferred: %s" % (SUBDIR or "(none)"))
    print("# sample directories found: %d" % len(samples))
    print()

    allele_rows, cov_rows, ars_rows, notes = [], [], [], []

    for s in samples:
        sdir, which = resolve_dir(ROOT, s)
        if which not in ("", SUBDIR):
            notes.append("%s\ttree\tread from: %s" % (s, which))

        row, anote = read_alleles(sdir)
        if anote:
            notes.append("%s\tallele_output\t%s" % (s, anote))
        if row is None:
            allele_rows.append([s] + ["NA"] * 16)
        else:
            vals = []
            for c in COLS:
                v = row.get(c)
                if v is None:
                    v = row.get(c.replace("HLA-", ""))
                vals.append(norm(v) or "NA")
            allele_rows.append([s] + vals)

        cov, lens = read_regions_bed(sdir, s)
        source = "regions.bed.gz" if cov else ""
        insuf = []

        lpath = find_log(sdir, s)
        if lpath:
            logcov, insuf = parse_log(lpath)
            if not cov:
                cov = logcov
                source = "log" if cov else ""
        elif not cov:
            notes.append("%s\tlog\tMISSING <sample>.hla_resolve.log" % s)

        if not cov:
            notes.append("%s\tcoverage\tno per-gene depths available" % s)
            cov_rows.append([s, "0", "", "", "", "", "", "", ""])
            ars_rows.append([s] + [""] * 8)
            continue
        if insuf:
            notes.append("%s\tcoverage\t%s" % (s, "; ".join(insuf)))

        num = sum(cov[x]["gene"] * lens[x]["gene"]
                  for x in GENES if x in lens and "gene" in lens[x])
        den = sum(lens[x]["gene"] for x in GENES if x in lens and "gene" in lens[x])
        weighted = (num / den) if den else None

        ars_vals = [cov.get(x, {}).get("ARS") for x in GENES]
        present = [v for v in ars_vals if v is not None]
        cov_rows.append([s, str(len(cov)),
                         fmt(mean([cov.get(x, {}).get("gene") for x in GENES])),
                         fmt(weighted),
                         fmt(mean([cov.get(x, {}).get("CDS") for x in GENES])),
                         fmt(mean([cov.get(x, {}).get("ARS") for x in GENES])),
                         fmt(min(present) if present else None),
                         mosdepth_total(sdir, s), source])
        ars_rows.append([s] + [fmt(v, 1) for v in ars_vals])

    def block(title, header, rows):
        print("=== %s ===" % title)
        print("\t".join(header))
        for r in rows:
            print("\t".join(str(x) for x in r))
        print("=== END %s ===" % title)
        print()

    block("BLOCK1 ALLELE CALLS", ["sample"] + COLS, allele_rows)
    block("BLOCK2 COVERAGE SUMMARY",
          ["sample", "n_genes", "mean_gene_depth", "len_weighted_gene_depth",
           "mean_CDS_depth", "mean_ARS_depth", "min_ARS_depth",
           "mosdepth_total_region", "source"], cov_rows)
    block("BLOCK3 PER-GENE ARS DEPTH", ["sample"] + GENES, ars_rows)

    print("=== BLOCK4 NOTES ===")
    for n in (notes or ["(none)"]):
        print(n)
    print("=== END BLOCK4 NOTES ===")
    print()
    print("# done: %d samples" % len(samples))


if __name__ == "__main__":
    main()
