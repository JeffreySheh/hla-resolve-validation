# HLA-Resolve WGS validation

Benchmark of [HLA-Resolve](https://github.com/matthewglasenapp/hla_resolve) v0.9.8 (IPD-IMGT/HLA release 3.64.0) on public whole-genome PacBio reads for 143 sequencing libraries from 129 samples. All libraries were sequenced by the Human Pangenome Reference Consortium (HPRC; cohorts HPRC, HPRC_PLUS) except the two Genome in a Bottle (GIAB) HG002 libraries. Every library was typed independently. Libraries sequenced across multiple SMRT cells were combined to a target genome-wide coverage depth of 30×. Samples with two libraries appear twice as independent measurements. The results are in the HLA-Resolve repository under [`docs/wgs_validation.md`](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/wgs_validation.md). This repository holds the scripts and data that produced them.

## Sample and library selection

Inputs were chosen by `scripts/select_inputs.py` from the HPRC PacBio HiFi data on S3 (`working/HPRC` and `working/HPRC_PLUS`), and the chosen files are listed in `data/manifest.tsv`. The goals were to type one sequencing library at a time, at a depth typical of a single whole-genome run (30x), using only standard WGS data.

**Included**

- One input BAM per sequencing library, identified from the HPRC file index or the BAM read group. A sample with two libraries at 20× or more is typed twice.
- SMRT cells from the same sequencing library were combined to reach a target depth of 30x but not exceed 40x.

**Excluded**

- Non-WGS libraries: Fiber-seq, AMPure, and dilution tests.
- Barcode cross-talk: files under 1 GB, or whose read group names another sample.
- Duplicate encodings of the same cell (`.5mc`, `.fixed_kinetics`, `.pm_v1.4.0` copies).
- Libraries below 20× when the same sample has a deeper library. A sample whose only library is below 20× was still typed.

A few HPRC metadata entries are wrong or missing (e.g. NA19909's four preps are all labelled as one library). `data/merge_overrides.tsv` lists these hand fixes, with a note on each.

## Scoring

Allele calls were compared to two reference sets:

- **Lai et al. 2024** ([doi.org/10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030); Supplementary File 6)
  - 42 HPRC release 1 samples, 49 libraries
  - HLA-A, -B, -C, -DPA1, -DPB1, -DQA1, -DQB1, -DRB1
  - Resolution: up through four-field (full) resolution
- **1000 Genomes Project HLA panel** ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt))
  - 87 samples, 94 libraries
  - HLA-A, -B, -C, -DQB1, -DRB1
  - Resolution: two-field resolution
  - Exome-based (IPD-IMGT/HLA 3.28)

Concordance is reported among alleles called. An allele was evaluated at a given field resolution only if the reference specified it to that resolution. A call resolved to fewer fields than the reference is considered discordant at that level of resolution. By default, genes do not type if their antigen recognition site coverage depth sits below 8x, so call rate was measured as the proportion of total possible allele calls emitted by HLA-Resolve. The two alleles of a gene are paired with the reference by the assignment that maximizes agreement; a reference ambiguity (`A*03:01/A*03:04`) accepts either allele.

## Layout

```
data/
  manifest.tsv                     sequencing read file manifest
  merge_plan.tsv                   assigns each read file to a final input BAM for HLA typing
  merge_overrides.tsv              metadata corrections used by select_inputs.py
  library_by_file.tsv              library and sample tags parsed from every downloaded BAM header
  samples_all.txt                  List of the 129 HPRC samples typed
  HPRC_PacBio_HiFi.file.index.csv  Downloaded HPRC PacBio HiFi file index
  lai_answer_key.tsv               Lai et al. reference typings (from their Supplementary File 6)
  kgp_answer_key.tsv               1000 Genomes panel reference typings
  kgp_star_alleles.tsv             panel alleles the 1000 Genomes file marks as exome-only (`*`)
scripts/
  select_inputs.py                 applies the selection rules and writes data/merge_plan.tsv
  build_manifest.py                creates data/manifest.tsv from merge_plan.tsv and the S3 listing
  pipeline/
    1_download.sbatch              downloads the read files (fetch_worker.sh handles one file)
    fetch_worker.sh
    2_merge.sbatch                 concatenates SMRT cells into the final input BAMs for HLA typing
    3_type.sbatch                  runs HLA-Resolve on each input BAM
    4_audit_results.py             collects allele calls and coverage into results/audit.txt
    5_test_against_truth.py        compares allele calls to the reference typings
results/
  audit.txt                        allele calls and coverage from every HLA-Resolve run
  calls.tsv                        HLA typing results per input BAM: calls, coverage, concordance
  misses.tsv                       discordant HLA allele calls
  summary.txt                      call rate and concordance tables
archive/                           the original pipeline and its results
```

## Workflow

Start from `data/manifest.tsv`. It was made by downloading a larger set of files per sample and applying the selection rules with `scripts/select_inputs.py`; the selection is already done, so reproducing the results starts with the download.

1. **Download read files.** `scripts/pipeline/1_download.sbatch` downloads the files in `data/manifest.tsv`.
2. **Merge SMRT cells.** `scripts/pipeline/2_merge.sbatch` concatenates read files from the same library sequenced across different SMRT cells when needed to reach ~30× genome-wide coverage depth. The files to be concatenated are specified by `data/merge_plan.tsv`, which was created during input selection (see above).
3. **Perform HLA typing.** `scripts/pipeline/3_type.sbatch` runs HLA-Resolve on each finalized BAM.
4. **Collect results.** `scripts/pipeline/4_audit_results.py` gathers the allele calls and coverage from every run into one file, `results/audit.txt`.
5. **Evaluate HLA allele concordance.** `scripts/pipeline/5_test_against_truth.py` compares the HLA allele calls from HLA-Resolve to the reference typings and writes `results/calls.tsv`, `results/misses.tsv`, and `results/summary.txt`.
