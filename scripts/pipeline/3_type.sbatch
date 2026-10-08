#!/bin/bash
#SBATCH --job-name=hla_type
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --no-requeue
# STEP 3: run hla_resolve on one input BAM.
#
#   PLAN=.../data/merge_plan.tsv MERGEDIR=/mnt/beegfs/.../merged RESULTS=/mnt/beegfs/.../results \
#   sbatch --array=1-N --output=.../logs/type_%A_%a.out 3_type.sbatch
#
# Array index N maps to the Nth input BAM in the plan (sorted, unique), the same
# order 2_merge.sbatch uses. Output goes to $RESULTS/<input_bam>/; hla_resolve
# nests its own <sample>/ inside that, so results land at
#   $RESULTS/<input_bam>/<sample>/hla_typing_results/allele_output.csv
# and the log at $RESULTS/<input_bam>/<sample>.hla_resolve.log.
#
# --no-requeue and a fresh output directory: a run restarted over partial
# output has produced silently wrong calls before (stale samtools sort spills).
# 80 GB covers the largest merged inputs seen (~40 GB peak on a 41 GB BAM).

set -uo pipefail
PLAN=${PLAN:?set PLAN to merge_plan.tsv}
MERGEDIR=${MERGEDIR:?set MERGEDIR}
RESULTS=${RESULTS:?set RESULTS}
THREADS=${THREADS:-$SLURM_CPUS_PER_TASK}
PLATFORM=${PLATFORM:-pacbio}
SUFFIX=${SUFFIX:-$([ "$PLATFORM" = ont ] && echo ont || echo hifi)}

ENTRY=$(tail -n +2 "$PLAN" | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$ENTRY" ] || { echo "array index past end"; exit 0; }
SAMPLE=${ENTRY%%__*}
BAM="$MERGEDIR/$ENTRY.$SUFFIX.bam"
OUT="$RESULTS/$ENTRY"
log() { echo "[$(date +%H:%M:%S)] $ENTRY: $*"; }

[ -f "$MERGEDIR/.$ENTRY.done" ] && [ -s "$BAM" ] || { log "merged BAM missing or unverified: $BAM"; exit 1; }
if [ -s "$OUT/$SAMPLE/hla_typing_results/allele_output.csv" ]; then log "already typed"; exit 0; fi
rm -rf "$OUT"; mkdir -p "$OUT"

command -v hla_resolve >/dev/null 2>&1 || { log "hla_resolve not on PATH (activate the env before sbatch)"; exit 1; }
log "hla_resolve $(hla_resolve --version 2>/dev/null | head -1) on $(readlink -f "$BAM")"

# ONT: hla_resolve blocks --platform ont with WGS until it is validated, and
# HLA_RESOLVE_ALLOW_ONT_WGS lifts the block for this pass. The Clair3 model
# follows the basecaller in the @RG DS field, for example
# dna_r10.4.1_e8.2_400bps_sup@v4.3.0 -> r1041_e82_400bps_sup_v430.
EXTRA=()
if [ "$PLATFORM" = ont ]; then
  export HLA_RESOLVE_ALLOW_ONT_WGS=1
  BC=$(samtools view -H "$BAM" | grep -m1 '^@RG' | grep -o 'basecall_model=[^ 	]*' | cut -d= -f2)
  case "$BC" in
    dna_r10.4.1_e8.2_400bps_sup@v4.1.0) MODEL=r1041_e82_400bps_sup_v410 ;;
    dna_r10.4.1_e8.2_400bps_sup@v4.2.0) MODEL=r1041_e82_400bps_sup_v420 ;;
    dna_r10.4.1_e8.2_400bps_sup@v4.3.0) MODEL=r1041_e82_400bps_sup_v430 ;;
    dna_r10.4.1_e8.2_400bps_sup@v5.0.0) MODEL=r1041_e82_400bps_sup_v500 ;;
    *) log "unknown basecall model '$BC' in the @RG header"; exit 1 ;;
  esac
  log "basecall model $BC, Clair3 model $MODEL"
  EXTRA=(--clair3_model "$MODEL")
fi

# The conda env lives on NFS home; when many tasks import it in the same
# second a task can die at import with FileNotFoundError under site-packages.
# Such a death comes within seconds, so an attempt that ends inside 2 minutes
# is retried after a pause, up to three attempts.
for attempt in 1 2 3; do
  t0=$(date +%s)
  hla_resolve \
    --input_file "$BAM" \
    --sample_name "$SAMPLE" \
    --platform "$PLATFORM" \
    --scheme WGS \
    --output_dir "$OUT" \
    --threads "$THREADS" \
    ${EXTRA[@]+"${EXTRA[@]}"}
  rc=$?
  [ $rc -eq 0 ] && break
  [ $(( $(date +%s) - t0 )) -lt 120 ] || break
  log "attempt $attempt died within 2 min (rc=$rc); retrying in 30 s"
  sleep 30; rm -rf "$OUT"; mkdir -p "$OUT"
done
[ $rc -eq 0 ] && [ -s "$OUT/$SAMPLE/hla_typing_results/allele_output.csv" ] \
  && log "ok" || { log "FAILED rc=$rc"; exit 1; }
