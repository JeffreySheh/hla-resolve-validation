#!/bin/bash
#SBATCH --job-name=hla_type
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --no-requeue
# STEP 3 (plan-driven): run hla_resolve on one merged entry.
#
#   PLAN=.../data/merge_plan.tsv MERGEDIR=/mnt/beegfs/.../merged RESULTS=/mnt/beegfs/.../results \
#   sbatch --array=1-N --output=.../logs/type_%A_%a.out 3_type_plan.sbatch
#
# Array index N maps to the Nth entry in the plan (sorted, unique), the same
# order 2_merge_plan.sbatch uses. Output goes to $RESULTS/<entry>/; hla_resolve
# nests its own <sample>/ inside that, so results land at
#   $RESULTS/<entry>/<sample>/hla_typing_results/allele_output.csv
# and the log at $RESULTS/<entry>/<sample>.hla_resolve.log.
#
# --no-requeue and a fresh output directory: a run restarted over partial
# output has produced silently wrong calls before (stale samtools sort spills).
# 80 GB covers the largest merged inputs seen (~40 GB peak on a 41 GB BAM).

set -uo pipefail
PLAN=${PLAN:?set PLAN to merge_plan.tsv}
MERGEDIR=${MERGEDIR:?set MERGEDIR}
RESULTS=${RESULTS:?set RESULTS}
THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}

ENTRY=$(tail -n +2 "$PLAN" | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$ENTRY" ] || { echo "array index past end"; exit 0; }
SAMPLE=${ENTRY%%__*}
BAM="$MERGEDIR/$ENTRY.hifi.bam"
OUT="$RESULTS/$ENTRY"
log() { echo "[$(date +%H:%M:%S)] $ENTRY: $*"; }

[ -f "$MERGEDIR/.$ENTRY.done" ] && [ -s "$BAM" ] || { log "merged BAM missing or unverified: $BAM"; exit 1; }
if [ -s "$OUT/$SAMPLE/hla_typing_results/allele_output.csv" ]; then log "already typed"; exit 0; fi
rm -rf "$OUT"; mkdir -p "$OUT"

command -v hla_resolve >/dev/null 2>&1 || { log "hla_resolve not on PATH (activate the env before sbatch)"; exit 1; }
log "hla_resolve $(hla_resolve --version 2>/dev/null | head -1) on $(readlink -f "$BAM")"

# The conda env lives on NFS home; when many tasks import it in the same
# second a task can die at import with FileNotFoundError under site-packages.
# Such a death comes within seconds, so an attempt that ends inside 2 minutes
# is retried after a pause, up to three attempts.
for attempt in 1 2 3; do
  t0=$(date +%s)
  hla_resolve \
    --input_file "$BAM" \
    --sample_name "$SAMPLE" \
    --platform pacbio \
    --scheme WGS \
    --output_dir "$OUT" \
    --threads "$THREADS" \
    --clean_up
  rc=$?
  [ $rc -eq 0 ] && break
  [ $(( $(date +%s) - t0 )) -lt 120 ] || break
  log "attempt $attempt died within 2 min (rc=$rc); retrying in 30 s"
  sleep 30; rm -rf "$OUT"; mkdir -p "$OUT"
done
[ $rc -eq 0 ] && [ -s "$OUT/$SAMPLE/hla_typing_results/allele_output.csv" ] \
  && log "ok" || { log "FAILED rc=$rc"; exit 1; }
