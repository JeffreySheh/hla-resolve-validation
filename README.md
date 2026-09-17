# HLA-Resolve WGS validation

Benchmark of [HLA-Resolve](https://github.com/matthewglasenapp/hla_resolve) v0.9.8 (IPD-IMGT/HLA release 3.64.0) on public whole-genome PacBio reads for 143 sequencing libraries from 129 samples. All libraries were sequenced by the Human Pangenome Reference Consortium (HPRC; cohorts HPRC, HPRC_PLUS) except the two Genome in a Bottle (GIAB) HG002 libraries. Every library was typed independently. Libraries sequenced across multiple SMRT cells were combined to a target genome-wide coverage depth of 30×. Samples with two libraries appear twice as independent measurements. The results are in the HLA-Resolve repository under [`docs/wgs_validation.md`](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/wgs_validation.md). This repository holds the scripts and data that produced them.

## Sample and library selection

Inputs were chosen by `scripts/select_inputs.py` from the HPRC PacBio HiFi data on S3 (`working/HPRC` and `working/HPRC_PLUS`), and the chosen files are listed in `data/manifest.tsv`. The goals were to type one sequencing library at a time, at a depth typical of a single whole-genome run (30x), using only standard WGS data.

**Included**

- One entry per sequencing library, identified from the HPRC file index or the BAM read group. A sample with two libraries at 20× or more is typed twice.
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
  manifest.tsv                     every input file per entry: sample, cohort, library, instrument, filename, S3/HTTPS URI, bytes, index coverage
  merge_plan.tsv                   select_inputs.py output the manifest was built from
  merge_overrides.tsv              metadata corrections used by select_inputs.py
  library_by_file.tsv              BAM header scan (LB/PU/SM) of every downloaded file, input to select_inputs.py
  samples_all.txt                  the 129 HPRC-sequenced samples
  HPRC_PacBio_HiFi.file.index.csv  HPRC metadata snapshot used for library and coverage
  s3_listing.tsv                   cached bucket listing behind manifest.tsv
  lai_answer_key.tsv               Lai et al. reference typings (from their Supplementary File 6)
  kgp_answer_key.tsv               1000 Genomes panel, converted (ambiguities expanded, `*` flags stripped into kgp_star_alleles.tsv, corrupted numeric cells decoded)
scripts/
  select_inputs.py                 rules -> data/merge_plan.tsv
  build_manifest.py                merge_plan.tsv + S3 listing -> data/manifest.tsv
  pipeline/1_download.sbatch, fetch_worker.sh   download from a manifest (verified, resumable)
  pipeline/2_merge.sbatch          samtools cat per entry, from a plan
  pipeline/3_type.sbatch           hla_resolve per entry
  audit_results.py                 results tree -> one audit file (calls, coverage, per-gene ARS depth)
  grade.py                         audit -> results/ (both references); imports test_against_truth.py
  test_against_truth.py            grading functions
  build_doc.py                     results -> docs/wgs_validation.md in the hla_resolve repository
results/
  audit.txt                        raw audit of the typing run
  calls.tsv                        one row per entry: metadata, depth, per-gene ARS depth, call rate, concordance by field, alleles
  misses.tsv                       every discordant gene x entry
  summary.txt                      call rate and concordance tables
archive/                           the original pipeline and its results (see History)
```

## Workflow

Steps 1 and 2 were run once; their outputs are committed. A re-run starts at step 3 from the committed manifest.

1. **Select inputs.** `select_inputs.py` applies the rules above to the header scan of every downloaded BAM (`data/library_by_file.tsv`) and the HPRC file index, and writes `data/merge_plan.tsv`: one row per input file, grouped into entries.
2. **Build the manifest.** `build_manifest.py` looks up the S3 object behind each planned file and writes `data/manifest.tsv`, the fixed list of inputs with cohort, library, instrument, URI, and size.
3. **Download.** `1_download.sbatch` fetches the manifest's files, one array task per sample, and leaves a `.done` marker next to each verified file.
4. **Merge.** `2_merge.sbatch` concatenates each entry's files with `samtools cat` into one unaligned BAM (`merged/<entry>.hifi.bam`). An entry with a single file is symlinked.
5. **Type.** `3_type.sbatch` runs HLA-Resolve on each merged BAM. Results land in `results/<entry>/<sample>/hla_typing_results/`, with the run log beside them.
6. **Audit.** `audit_results.py` reads every result directory and writes one file, `results/audit.txt`, with the allele calls, coverage summary, and per-gene ARS depth.
7. **Grade.** `grade.py` scores the audit against both reference sets and writes `results/calls.tsv` (one row per entry), `results/misses.tsv`, and `results/summary.txt`.
8. **Document.** `build_doc.py` turns the graded results into `docs/wgs_validation.md` in the HLA-Resolve repository.
