# WGS validation re-run: plan

Branch `matt/1kgp-validation`. Cluster clone: `/mnt/beegfs/cluster/scratch/mglasena/hla-resolve-validation`.
Read-only source data: `/mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam`. All new writes go under `/mnt/beegfs/cluster/scratch/mglasena/`.

## Goal

One benchmark of HLA-Resolve on PacBio WGS, built from a recorded manifest of input files, graded against two references (Lai et al. for the HPRC year-1 samples, the 1000 Genomes panel for the rest), and published as `hla_resolve/docs/wgs_validation.md`.

## Rules (reported in the manuscript methods; applied once)

Implemented in `scripts/select_inputs.py`, corrections in `data/merge_overrides.tsv`.

1. Unit of the benchmark is an **entry**: one sample x one library prep, id `<sample>__<prep>`. A sample sequenced from two preps is two entries — independent measurements.
2. Library identity from the HPRC file index (`library_ID`, movie suffix stripped); BAM `LB` where the index has no record; older untagged `.ccs.bam` treated as one prep per sample.
3. Excluded before selection: unverified downloads; files under 1 GB (barcode cross-talk); files whose `LB`/`SM` names another sample; non-WGS preps (Fiber-seq, Ampure, dilution tests); second encodings of a movie (`fixed_kinetics`, `pm_v1.4.0`, `.5mc` vs plain).
4. Depth: cells added in chronological order to an estimated 30x genome-wide; a cell that would exceed 40x is left out once above 20x; a prep under 20x on its own is not an entry, except that the sample's deepest prep may merge in the next one when it is below 20x (flagged).
5. Depth estimate: HPRC index coverage per file, or bytes x a per-file-type rate; pinned by override for the six HPRC_PLUS pools the index does not cover.
6. Merge with `samtools cat`; no duplicate removal, no kinetics stripping.

## Steps

| # | Step | Script | Status |
|---|------|--------|--------|
| 1 | Download the 40 Sequel IIe cells for 13 thin samples | `1_download.sbatch` + `fetch_worker.sh` (Jeffrey's) with `data/download_manifest.tsv` | running, job 5575328 / 5575341 |
| 2 | Re-scan BAM headers over both roots -> `data/library_by_file.tsv` | one-liner, see below | after 1 |
| 3 | Regenerate the plan | `select_inputs.py data/samples_all.txt` -> `data/merge_plan.tsv` | after 2 |
| 4 | Merge every entry | `2_merge_plan.sbatch`, array over entries | after 3 |
| 5 | Type every entry, one hla_resolve version | `3_type_plan.sbatch`, array over entries | after 4 |
| 6 | Audit outputs -> one tsv | `audit_results.py <results root>` | after 5 |
| 7 | Grade both truths, aggregate | `test_against_truth.py`, `test_against_1kgp.py`, `aggregate_validation.py` | after 6 |
| 8 | Write the doc | `build_doc.py` (to write) -> `hla_resolve/docs/wgs_validation.md` | after 7 |
| 9 | Consolidate the repo | see below | after 8 validates |
| 10 | PR to Jeffrey's `main`; message with the rule changes | | after 9 |

Header re-scan (step 2), run with the `hla_resolve` env active:

```bash
for J in /mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam /mnt/beegfs/cluster/scratch/mglasena/pacbio_bam; do
for d in $J/*/; do s=$(basename $d); [[ $s == logs || $s == merged ]] && continue
  for b in $d*.bam; do [[ -f $b ]] || continue; f=$(basename $b)
    samtools view -H "$b" | awk -F'\t' -v s="$s" -v f="$f" -v sz="$(stat -c %s "$b")" -v dn="$([[ -f $d.$f.done ]] && echo done || echo nodone)" \
      '$1=="@RG"{lb="";pu="";sm="";for(i=2;i<=NF;i++){if($i~/^LB:/)lb=substr($i,4);if($i~/^PU:/)pu=substr($i,4);if($i~/^SM:/)sm=substr($i,4)}; print s"\t"f"\t"sz"\t"dn"\t"lb"\t"pu"\t"sm; exit}
       END{if(!NR)print s"\t"f"\t"sz"\t"dn"\t\t\t"}'
  done
done; done > ~/library_by_file.tsv
```

Merge and type (steps 4-5), after `git pull` on the cluster:

```bash
P=/mnt/beegfs/cluster/scratch/mglasena; R=$P/hla-resolve-validation
cp $R/data/merge_plan.tsv $P/plan_$(date +%Y%m%d).tsv; PLAN=$P/plan_$(date +%Y%m%d).tsv   # frozen copy: a task maps its array index to the Nth entry of PLAN when it starts, so PLAN must not change while the array runs; never point it at the repo file
X=--exclude=node-29,node-31,node-33,node-34,node-35,node-36,node-37,node-38,node-39,node-40,node-113,node-115   # bad, slow, or often-preempted nodes, see below
N=$(tail -n +2 $PLAN | cut -f1 | sort -u | wc -l)
MERGEDIR=$P/merged SEARCH_ROOTS="$P/pacbio_bam:/mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam" \
  sbatch --array=1-$N $X --output=$P/logs/merge_%A_%a.out $R/scripts/pipeline/2_merge_plan.sbatch
# when all merges report ok:
MERGEDIR=$P/merged RESULTS=$P/results \
  sbatch --array=1-$N%15 $X --output=$P/logs/type_%A_%a.out $R/scripts/pipeline/3_type_plan.sbatch
```

Use `/mnt/beegfs/...` paths in anything that runs on a node; `/hb/scratch` is a login-node symlink.

Nodes to exclude on windfall: node-29, node-33, node-34, node-35, node-36 (Jeffrey's list); node-31, node-37 (slow S3, 2026-09-14); node-38, node-39, node-40 (lab-colibri owners preempt windfall jobs often); node-113, node-115 (crashed mid-run 2026-09-14). Add to the list when a new one shows up.

## Doc layout (step 8)

One table, one row per entry, grouped by sample:

| Column | Values |
|--------|--------|
| `sample` | |
| `library` | prep name |
| `cohort` | `HPRC`, `HPRC_PLUS`, `GIAB` — who sequenced, from the S3 prefix |
| `instrument` | Revio / Sequel II |
| `truth` | `Lai` (doi.org/10.1016/j.csbj.2024.03.030) or `1kGP` (ftp.1000genomes.ebi.ac.uk/.../20181129_HLA_types_full_1000_Genomes_Project_panel.txt) |
| depth, genes typed, alleles called, 1f-4f | 3f/4f blank for 1kGP truth |

Summary tables per truth; a combined 1-2 field row on A, B, C, DQB1, DRB1 over all entries. Input files per entry from the manifest. Truth caveats stated: 1kGP is exome-based PolyPheMe at IMGT 3.28, five genes, two fields; `*`-flagged alleles; `C*17:71` postdates the panel.

## Consolidation (step 9)

Selection runs once; what ships is its output plus deterministic reproduction.

```
data/
  manifest.tsv                     entry, sample, cohort, library, instrument, filename, s3_uri
  merge_overrides.tsv
  lai_answer_key.tsv  kgp_answer_key.tsv  kgp_star_alleles.tsv
  HPRC_PacBio_HiFi.file.index.csv  metadata snapshot the selection used
scripts/
  select_inputs.py                 the rules -> data/manifest.tsv
  download.sbatch + fetch_worker.sh   reads manifest.tsv            (Jeffrey's, kept)
  merge.sbatch                     reads manifest.tsv, samtools cat  (from 2_merge_plan)
  type.sbatch                      per entry                         (from 3_type_plan)
  audit_results.py
  grade.py                         both truths, entry-aware          (test_against_truth + test_against_1kgp)
  build_doc.py
results/
  calls.tsv  concordance_by_entry.tsv  misses.txt  summary.txt
README.md                          what changed from the original pipeline and why
```

Delete after step 8 validates: `build_kgp_manifest.py`, `library_by_file.tsv`, `movies_1kgp.tsv`, `merged_readgroups.tsv`, `pacbio_manifest_jsheh.tsv`, `audit_1kgp.txt`, `samples_*.txt`, `merge_plan_excluded.tsv`, `download_manifest.tsv`, `hprc_calls.csv`, `2_merge.sbatch`, `3_type.sbatch`, `grade_1kg.py`, `hla_validation_final.csv`, `hla_validation_1kg.csv`, `methods.md` (old rule), `matt_download_pilot/`.

## Divergence from the original

Kept: download step and worker. Replaced: chemistry-count merge rule -> library-aware plan; per-sample typing -> per-entry; greedy grader that counts no-calls as misses -> unordered best-assignment grader with call rate separate. The original stays on `main` history.

## Open decisions

- hla_resolve version for the re-run (whatever is in the cluster env at submit time; record it).
- Whether a per-gene depth gate (10x / 15x) is reported alongside the ungated concordance.
- HG005 (Lai truth, only an 8x single cell): exclude.
