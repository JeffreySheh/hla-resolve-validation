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

Steps in order. Steps 1–2 were run once and their outputs are committed; a re-run starts at step 3 from the committed manifest.

| # | Step | Script | Reads | Writes |
|---|------|--------|-------|--------|
| 1 | Select inputs by the rules above | `select_inputs.py data/samples_all.txt` | `data/library_by_file.tsv` (BAM header scan), `data/HPRC_PacBio_HiFi.file.index.csv`, `data/merge_overrides.tsv` | `data/merge_plan.tsv` (entry, sample, library, filename, size, depth, note) |
| 2 | Record the S3 object behind every planned file | `build_manifest.py` | `data/merge_plan.tsv`, HPRC bucket listing (cached in `data/s3_listing.tsv`) | `data/manifest.tsv` |
| 3 | Download | `pipeline/1_download.sbatch` + `fetch_worker.sh`, one array task per sample | a 4-column cut of `data/manifest.tsv` (sample, filename, URI, bytes) | `<scratch>/pacbio_bam/<sample>/<file>` with a `.done` marker per verified file |
| 4 | Merge | `pipeline/2_merge.sbatch`, one array task per entry | a frozen copy of `data/merge_plan.tsv`, the downloaded files | `<scratch>/merged/<entry>.hifi.bam` (+ `.done`); a single-file entry is a symlink |
| 5 | Type | `pipeline/3_type.sbatch`, one array task per entry | the merged BAMs | `<scratch>/results/<entry>/<sample>/hla_typing_results/allele_output.csv`, run log `<entry>/<sample>.hla_resolve.log` |
| 6 | Audit | `audit_results.py <results root>` | the results tree | `results/audit.txt` (allele calls, coverage summary, per-gene ARS depth, notes) |
| 7 | Grade | `grade.py results/audit.txt` | `results/audit.txt`, `data/manifest.tsv`, the two answer keys | `results/calls.tsv`, `results/misses.tsv`, `results/summary.txt` |
| 8 | Document | `build_doc.py` | `results/calls.tsv`, `results/misses.tsv`, `data/manifest.tsv` | `docs/wgs_validation.md` in the HLA-Resolve repository |

Commands as run on the UCSC cluster (`$P` scratch, `$R` this repository):

```bash
P=/mnt/beegfs/cluster/scratch/<user>; R=$P/hla-resolve-validation
X=--exclude=node-29,node-31,node-33,node-34,node-35,node-36,node-37,node-38,node-39,node-40,node-113,node-115   # windfall nodes that are slow or often preempted

# 3. download (the two GIAB HG002 files are HTTPS: wget -c <url>, then touch .<file>.done beside it)
awk -F'\t' 'NR>1 && $7 ~ /^s3:/ {print $2"\t"$6"\t"$7"\t"$8}' $R/data/manifest.tsv > $P/pacbio_bam/pacbio_manifest.tsv
N=$(cut -f1 $P/pacbio_bam/pacbio_manifest.tsv | sort -u | wc -l)
OUTDIR=$P/pacbio_bam WORKER=$R/scripts/pipeline/fetch_worker.sh sbatch --array=1-$N $X $R/scripts/pipeline/1_download.sbatch

# 4. merge; PLAN is a frozen copy because a task maps its array index to the Nth entry of PLAN when it starts
cp $R/data/merge_plan.tsv $P/plan.tsv; PLAN=$P/plan.tsv
N=$(tail -n +2 $PLAN | cut -f1 | sort -u | wc -l)
MERGEDIR=$P/merged SEARCH_ROOTS="$P/pacbio_bam" sbatch --array=1-$N%40 $X $R/scripts/pipeline/2_merge.sbatch

# 5. type (conda env with hla_resolve active; 24 CPU, ~3 h and ~27 GB peak per entry)
MERGEDIR=$P/merged RESULTS=$P/results sbatch --array=1-$N%10 $X $R/scripts/pipeline/3_type.sbatch

# 6-8. audit, grade, document
python3 $R/scripts/audit_results.py $P/results > $R/results/audit.txt
python3 $R/scripts/grade.py $R/results/audit.txt
python3 $R/scripts/build_doc.py
```

Notes from the run: the shared filesystem, not CPU, sets the pace (about 0.5 GB/s in total; merging 9 TB took ~8 h at 40 concurrent streams, typing ~3 h per entry at 10–20 concurrent). Windfall jobs get preempted and the array scripts use `--no-requeue`, so a failed task is resubmitted by index. A typing task that dies at import (NFS hiccup on the conda env) retries itself.

## History

The first version of this pipeline (Jeffrey Sheh, 2026; kept under `archive/`) downloaded every BAM per sample, merged the Revio cells if there were at least three and the Sequel II cells otherwise, typed each sample once, and graded against Lai (40 samples) and the 1000 Genomes panel (89 samples). An audit of the BAM headers found the reasons for the re-run:

- the merge rule keyed on the movie-name prefix and skipped `m54*` Sequel IIe cells (16 samples affected; one sample typed from 28,786 reads);
- it had no notion of library: five samples merged two or three library preparations, one of them Fiber-seq;
- barcode cross-talk stubs (sub-GB files carrying another sample's reads) entered merges;
- for 13 samples the Sequel IIe cells were in an S3 subfolder a non-recursive listing never showed, so they ran from one thin Revio cell (6–22×);
- ten HPRC samples ran at over 40× from four merged cells;
- the 1000 Genomes grader counted no-calls as misses, dropped `*`-flagged alleles, and paired alleles greedily.

The download step and worker are unchanged from that version and were used for the re-run. The merge and typing scripts were replaced by plan-driven versions, the grader by `grade.py`, and the sample-level results by library-level ones.
