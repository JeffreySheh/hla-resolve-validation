# HLA-Resolve WGS validation

Benchmark of [HLA-Resolve](https://github.com/matthewglasenapp/hla_resolve) on public PacBio whole-genome sequencing: 146 sequencing libraries from 130 samples, typed on one tool version and graded against two reference typings. The published summary is [`docs/wgs_validation.md`](https://github.com/matthewglasenapp/hla_resolve/blob/main/docs/wgs_validation.md) in the HLA-Resolve repository; this repository holds the input manifest, the pipeline that produced it, and the graded results.

## Design

- **Unit of the benchmark is a sequencing library** (one sample × one library preparation), named `<sample>__<library>`. A sample sequenced from two libraries is typed twice and reported twice; the two runs are independent measurements.
- **Two reference sets.** Lai et al. 2024 ([doi:10.1016/j.csbj.2024.03.030](https://doi.org/10.1016/j.csbj.2024.03.030)) for the HPRC release 1 samples: eight classical genes, up to four fields. The 1000 Genomes Project HLA panel ([20181129_HLA_types_full_1000_Genomes_Project_panel.txt](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/20181129_HLA_types_full_1000_Genomes_Project_panel.txt)) for the 1000 Genomes samples HPRC sequenced in later years: HLA-A, -B, -C, -DQB1, -DRB1 at two fields. A sample in both sets is graded against Lai.
- **Input selection is rule-based and recorded.** The rules ran once (`scripts/select_inputs.py`) and their output is the fixed file list in `data/manifest.tsv`; reproduction uses the manifest, not the rules.
- **One tool version** for every library (v0.9.8, IPD-IMGT/HLA 3.64.0), recorded in each run log.

### Selection rules

Applied per sample to every BAM under `working/HPRC/<sample>/raw_data/PacBio_HiFi/` and `working/HPRC_PLUS/…` on the HPRC bucket (listed recursively — some cells sit in `primrose/` or `wMods/` subfolders):

1. Library identity from the HPRC file index (`library_ID`, movie suffix stripped), else the BAM read-group `LB` tag, else "untagged" for old `.ccs.bam` files without tags.
2. Excluded: downloads that did not verify; files under 1 GB (barcode cross-talk stubs); files whose `LB`/`SM` names another sample; non-WGS preparations (Fiber-seq, AMPure, dilution tests); a second encoding of the same SMRT cell (`.5mc`, `.fixed_kinetics`, `.pm_v1.4.0`).
3. The deepest library is the sample's primary entry. Every other library at ≥20× on its own is an additional entry. A library under 20× is typed only when it is the sample's only one (flagged).
4. Depth cap: cells of a library are added in run order to an estimated 30× genome-wide (HPRC index coverage per file, or file size × a per-chemistry rate; pinned by override for six HPRC_PLUS pools the index does not cover). Once above 20×, a cell that would pass 40× is left out.
5. Cells are concatenated with `samtools cat`; no duplicate removal, no kinetics stripping.

Hand corrections to the metadata are in `data/merge_overrides.tsv` (each with a note).

### Grading

Call rate counts every allele emitted against two per gene per library; HLA-Resolve does not type a gene whose antigen recognition site is under 8× mean depth, and such genes count against call rate, not concordance. Concordance is measured among alleles called. An allele is evaluated at a field resolution only if the reference specifies it to that resolution; a call resolved to fewer fields than the reference is discordant. The two alleles of a gene are paired with the reference by the assignment that maximizes agreement; a reference ambiguity (`A*03:01/A*03:04`) accepts either allele. One Lai allele (`HG01358 DRB1*04:92#`, an uncertain call) is excluded.

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

## Reproducing the run

Paths below are the ones used on the UCSC cluster; `$P` is a scratch directory, `$R` this repository.

```bash
P=/mnt/beegfs/cluster/scratch/<user>; R=$P/hla-resolve-validation
X=--exclude=node-29,node-31,node-33,node-34,node-35,node-36,node-37,node-38,node-39,node-40,node-113,node-115   # windfall nodes that are slow or often preempted

# 1. download: the S3 files of the manifest, one array task per sample
awk -F'\t' 'NR>1 && $7 ~ /^s3:/ {print $2"\t"$6"\t"$7"\t"$8}' $R/data/manifest.tsv > $P/pacbio_bam/pacbio_manifest.tsv
N=$(cut -f1 $P/pacbio_bam/pacbio_manifest.tsv | sort -u | wc -l)
OUTDIR=$P/pacbio_bam WORKER=$R/scripts/pipeline/fetch_worker.sh sbatch --array=1-$N $X $R/scripts/pipeline/1_download.sbatch
# the two GIAB HG002 files are HTTPS (wget -c <url>; touch .<file>.done beside it)

# 2. merge, one array task per entry; PLAN is a frozen copy (a task maps its index to the Nth entry of PLAN when it starts)
cp $R/data/merge_plan.tsv $P/plan.tsv; PLAN=$P/plan.tsv
N=$(tail -n +2 $PLAN | cut -f1 | sort -u | wc -l)
MERGEDIR=$P/merged SEARCH_ROOTS="$P/pacbio_bam" sbatch --array=1-$N%40 $X $R/scripts/pipeline/2_merge.sbatch

# 3. type (conda env with hla_resolve active; 24 CPU, ~3 h and ~27 GB peak per entry)
MERGEDIR=$P/merged RESULTS=$P/results sbatch --array=1-$N%10 $X $R/scripts/pipeline/3_type.sbatch

# 4. audit, grade, document
python3 $R/scripts/audit_results.py $P/results > $R/results/audit.txt
python3 $R/scripts/grade.py $R/results/audit.txt
python3 $R/scripts/build_doc.py            # writes ../hla_resolve/docs/wgs_validation.md
```

The selection itself is re-run with `python3 scripts/select_inputs.py data/samples_all.txt` (needs a fresh `data/library_by_file.tsv` from the BAM headers) followed by `python3 scripts/build_manifest.py`.

Practical notes from the run: the shared filesystem, not CPU, sets the pace (≈0.5 GB/s total; merging 9 TB took ~8 h at 40 concurrent streams, typing ~3 h per entry at 10–20 concurrent); windfall jobs are preempted and the array scripts use `--no-requeue`, so failed tasks are resubmitted by index; a task that dies at import (NFS hiccup on the conda env) is retried by `3_type.sbatch`.

## History

The first version of this pipeline (Jeffrey Sheh, 2026; kept under `archive/`) downloaded every BAM per sample, merged the Revio cells if there were at least three and the Sequel II cells otherwise, typed each sample once, and graded against Lai (40 samples) and the 1000 Genomes panel (89 samples). An audit of the BAM headers found the reasons for the re-run:

- the merge rule keyed on the movie-name prefix and skipped `m54*` Sequel IIe cells (16 samples affected; one sample typed from 28,786 reads);
- it had no notion of library: five samples merged two or three library preparations, one of them Fiber-seq;
- barcode cross-talk stubs (sub-GB files carrying another sample's reads) entered merges;
- for 13 samples the Sequel IIe cells were in an S3 subfolder a non-recursive listing never showed, so they ran from one thin Revio cell (6–22×);
- ten HPRC samples ran at over 40× from four merged cells;
- the 1000 Genomes grader counted no-calls as misses, dropped `*`-flagged alleles, and paired alleles greedily.

The download step and worker are unchanged from that version and were used for the re-run. The merge and typing scripts were replaced by plan-driven versions, the grader by `grade.py`, and the sample-level results by library-level ones.
