#!/bin/bash
# Submit ONT stage 4 in batches, so that only about two batches are on scratch
# at a time. Run on the login node from the cluster clone, with the hla_resolve
# env active:
#
#   BATCHES="b01" bash scripts/pipeline/run_batches.sh             # the pilot
#   BATCHES="b02 b03 b04" bash scripts/pipeline/run_batches.sh     # the rest
#
# Each batch is four jobs. Tasks follow their own sample (aftercorr), and a task
# whose sample failed upstream is cancelled, not left pending:
#   download (array)  ->  merge (array, deletes the raw runs it merged)
#                     ->  type (array)  ->  cleanup (one job, deletes the
#                         merged BAM of every sample that typed with no error)
# Batch N+1 downloads only after batch N has merged and batch N-1 is cleaned up.
# A sample that fails keeps its inputs, and its batch report says why:
#   $W/<batch>.cleanup.tsv
#
# Settings (environment):
#   W              work root (default /mnt/beegfs/cluster/scratch/mglasena/ont_wgs_stage4)
#   BATCHES        batch names from data/ont/stages/stage4/ (default: all)
#   AFTER          a dependency for the first batch, for example afterany:12345
#   CLAIR3_MODEL   default r1041_e82_400bps_sup_v500 (the SIF has no sup v5.2.0 model)
#   SLURM_ARGS     partition flags for every job (default: the lab partitions)
set -euo pipefail

V=$(cd "$(dirname "$0")/../.." && pwd)
P=$V/scripts/pipeline
S=$V/data/ont/stages/stage4
W=${W:-/mnt/beegfs/cluster/scratch/mglasena/ont_wgs_stage4}
RAW=$W/raw; MERGED=$W/merged; RESULTS=$W/results; LOGS=$W/logs
CLAIR3_MODEL=${CLAIR3_MODEL:-r1041_e82_400bps_sup_v500}
SLURM_ARGS=${SLURM_ARGS:---partition=lab-colibri,lab-colibri-hmem --account=pi-jkoc --qos=pi-jkoc}
BATCHES=${BATCHES:-$(cd "$S" && ls b*.merge_plan.tsv | sed 's/\.merge_plan\.tsv$//')}

command -v hla_resolve >/dev/null || { echo "activate the hla_resolve env first"; exit 1; }
mkdir -p "$RAW" "$MERGED" "$RESULTS" "$LOGS"

join() { local IFS=,; echo "$*"; }
sub() { sbatch --parsable --kill-on-invalid-dep=yes $SLURM_ARGS "$@"; }

M_PREV=""; C_PREV=""; C_PREV2=""
for B in $BATCHES; do
    PLAN=$S/$B.merge_plan.tsv; DL=$S/$B.download_manifest.tsv
    [ -s "$PLAN" ] && [ -s "$DL" ] || { echo "no plan for $B"; exit 1; }
    N=$(tail -n +2 "$PLAN" | cut -f1 | sort -u | wc -l)

    DEPS=()
    [ -z "$M_PREV" ] && [ -n "${AFTER:-}" ] && DEPS+=("$AFTER")
    [ -n "$M_PREV" ] && DEPS+=("afterany:$M_PREV")
    [ -n "$C_PREV2" ] && DEPS+=("afterany:$C_PREV2")
    DEP=$(join ${DEPS[@]+"${DEPS[@]}"})

    # The download script asks for --exclusive, which waits a long time on
    # busy nodes, so it is removed here as in the earlier stages.
    D=$(sed '/^#SBATCH --exclusive/d' "$P/1_download.sbatch" | sub --array=1-"$N" \
        --job-name="dl_$B" --output="$LOGS/$B.download_%A_%a.out" ${DEP:+--dependency=$DEP} \
        --export=ALL,OUTDIR="$RAW",BUILD_MANIFEST=0,MANIFEST="$DL",WORKER="$P/fetch_worker.sh")
    M=$(sub --array=1-"$N" --job-name="merge_$B" --output="$LOGS/$B.merge_%A_%a.out" \
        --dependency=aftercorr:"$D" \
        --export=ALL,PLAN="$PLAN",MERGEDIR="$MERGED",SEARCH_ROOTS="$RAW",SUFFIX=ont,DELETE_RAW="$RAW" \
        "$P/2_merge.sbatch")
    T=$(sub --array=1-"$N" --job-name="type_$B" --output="$LOGS/$B.type_%A_%a.out" \
        --dependency=aftercorr:"$M" \
        --export=ALL,PLAN="$PLAN",MERGEDIR="$MERGED",RESULTS="$RESULTS",PLATFORM=ont,CLAIR3_MODEL="$CLAIR3_MODEL" \
        "$P/3_type.sbatch")
    C=$(sub --job-name="clean_$B" --output="$LOGS/$B.cleanup_%j.out" --dependency=afterany:"$T" \
        --export=ALL,PLAN="$PLAN",MERGEDIR="$MERGED",RAWDIR="$RAW",RESULTS="$RESULTS",REPORT="$W/$B.cleanup.tsv",SUFFIX=ont \
        "$P/6_cleanup.sbatch")
    echo "$B ($N samples): download $D, merge $M, type $T, cleanup $C${DEP:+, waits on $DEP}"
    M_PREV=$M; C_PREV2=$C_PREV; C_PREV=$C
done
