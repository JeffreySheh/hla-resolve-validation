#!/bin/bash
#SBATCH --job-name=hla_download
#SBATCH --output=/mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam/logs/download_%A_%a.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --requeue
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
# Run with: sbatch --array=1-44 --exclude=node-29,node-31,node-33,node-34,node-35,node-36,node-37 ~/1_download.sbatch

[ -f /software/scripts/utility/requeue.feature ] && source /software/scripts/utility/requeue.feature
set -uo pipefail

export OUTDIR=${OUTDIR:-/mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam}
mkdir -p "$OUTDIR"

find_beside() {
    local name=$1 d
    for d in "$OUTDIR" "${SLURM_SUBMIT_DIR:-$PWD}" "$PWD" "$HOME"; do
        [ -f "$d/$name" ] && { echo "$d/$name"; return 0; }
    done
    return 1
}

SAMPLE_LIST=${SAMPLE_LIST:-$(find_beside samples.tsv)} || { echo "cannot find samples.tsv"; exit 1; }
WORKER=${WORKER:-$(find_beside fetch_worker.sh)} || { echo "cannot find fetch_worker.sh"; exit 1; }
MANIFEST=${MANIFEST:-$OUTDIR/pacbio_manifest.tsv}

BUCKET=${BUCKET:-human-pangenomics}
PREFIX=${PREFIX:-working}
SUBDIR=${SUBDIR:-raw_data/PacBio_HiFi}

NJOBS=${NJOBS:-4}
S3_CONC=${S3_CONC:-8}

export FORCE=${FORCE:-0}
export VERIFY_MD5=${VERIFY_MD5:-0}

if ! command -v aws >/dev/null 2>&1; then
    if ! command -v module >/dev/null 2>&1; then
        for init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh /etc/profile.d/modules.sh; do
            [ -f "$init" ] && { source "$init"; break; }
        done
    fi
    module load "${AWS_MODULE:-aws}" 2>/dev/null || module load aws/2.13.15 2>/dev/null
fi
command -v aws >/dev/null 2>&1 || { echo "FATAL: aws not on PATH"; exit 1; }
echo "aws: $(command -v aws)"

export QUICKCHECK=1
command -v samtools >/dev/null 2>&1 || module load "${SAMTOOLS_MODULE:-samtools}" 2>/dev/null
if command -v samtools >/dev/null 2>&1; then
    echo "samtools: $(command -v samtools)"
else
    echo "samtools: not found, using BGZF EOF check"
    export QUICKCHECK=0
fi

build_manifest() {
    echo "building manifest from $SAMPLE_LIST ..."
    local TMP; TMP=$(mktemp -d)
    while IFS=$'\t' read -r sample cohort want; do
        [ -n "${sample:-}" ] || continue
        case "$sample" in \#*) continue ;; esac
        (
          try_cohorts="${cohort:-HPRC HPRC_PLUS}"
          [ -n "$cohort" ] && try_cohorts="$cohort"
          for c in $try_cohorts; do
            out=$(aws s3 ls --no-sign-request "s3://$BUCKET/$PREFIX/$c/$sample/$SUBDIR/" 2>/dev/null \
              | awk -v s="$sample" -v c="$c" \
                  '$4 ~ /\.bam$/ && $4 !~ /fail_reads/ {print s"\t"c"\t"$4"\t"$3}')
            [ -n "$out" ] && { printf '%s\n' "$out"; break; }
          done
        ) > "$TMP/$sample.ls" &
    done < "$SAMPLE_LIST"
    wait

    local MISSING=0
    : > "$TMP/rows"
    while IFS=$'\t' read -r sample cohort want; do
        [ -n "${sample:-}" ] || continue
        case "$sample" in \#*) continue ;; esac
        if [ ! -s "$TMP/$sample.ls" ]; then
            echo "ERROR: no bams found for $sample" >&2; MISSING=1; continue
        fi
        cat "$TMP/$sample.ls" >> "$TMP/rows"
    done < "$SAMPLE_LIST"

    [ "$MISSING" = 0 ] || { echo "aborting, some samples not found" >&2; rm -rf "$TMP"; return 1; }

    sort -t"$(printf '\t')" -k4,4nr "$TMP/rows" \
      | awk -F'\t' -v b="$BUCKET" -v p="$PREFIX" -v d="$SUBDIR" \
          '{printf "%s\t%s\ts3://%s/%s/%s/%s/%s/%s\t%s\n", $1, $3, b, p, $2, $1, d, $3, $4}' \
      > "$MANIFEST.tmp"
    mv -f "$MANIFEST.tmp" "$MANIFEST"

    local n tot
    n=$(wc -l < "$TMP/rows"); tot=$(awk -F'\t' '{t+=$4} END{print t}' "$TMP/rows")
    awk -v n="$n" -v t="$tot" 'BEGIN{printf "built manifest: %d files, %.2f TB\n", n, t/1e12}' >&2
    rm -rf "$TMP"
}

if [ "${SLURM_ARRAY_TASK_ID:-1}" = "${SLURM_ARRAY_TASK_MIN:-1}" ]; then
    if [ ! -s "$MANIFEST" ] || [ "$SAMPLE_LIST" -nt "$MANIFEST" ]; then
        build_manifest || exit 1
    fi
else
    for _ in $(seq 1 60); do [ -s "$MANIFEST" ] && break; sleep 5; done
fi
[ -s "$MANIFEST" ] || { echo "manifest never appeared"; exit 1; }
echo "manifest: $MANIFEST"

THIS_SAMPLE=$(cut -f1 "$MANIFEST" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID:-1}p")
[ -n "$THIS_SAMPLE" ] || { echo "array index past end"; exit 0; }
echo "array task ${SLURM_ARRAY_TASK_ID:-1} -> $THIS_SAMPLE"

WORK_MANIFEST=$OUTDIR/.manifest_row.${SLURM_ARRAY_TASK_ID:-1}.$$
awk -F'\t' -v s="$THIS_SAMPLE" '$1==s' "$MANIFEST" > "$WORK_MANIFEST"
trap 'rm -f "$WORK_MANIFEST"' EXIT

TOTAL=$(wc -l < "$WORK_MANIFEST")
[ "$TOTAL" -gt 0 ] || { echo "no files for $THIS_SAMPLE"; exit 1; }

while IFS=$'\t' read -r s f u sz; do
    rm -f "$OUTDIR/$s/$f".part*
done < "$WORK_MANIFEST"

echo "concurrency: $NJOBS files x $S3_CONC requests"

export AWS_CONFIG_FILE=$OUTDIR/.aws_config_job.${SLURM_ARRAY_TASK_ID:-1}
{
  printf '[default]\n'
  printf 'retry_mode = adaptive\n'
  printf 'max_attempts = 10\n'
  printf 's3 =\n'
  printf '    max_concurrent_requests = %s\n' "$S3_CONC"
  printf '    multipart_threshold = 64MB\n'
  printf '    multipart_chunksize = 64MB\n'
} > "$AWS_CONFIG_FILE"
export AWS_EC2_METADATA_DISABLED=true

awk -F'\t' '{print $1, $2, $3, $4}' "$WORK_MANIFEST" \
  | xargs -P "$NJOBS" -n 4 bash "$WORKER"

NDONE=0
while IFS=$'\t' read -r s f u sz; do
    [ -f "$OUTDIR/$s/.$f.done" ] && [ -s "$OUTDIR/$s/$f" ] && NDONE=$(( NDONE + 1 ))
done < "$WORK_MANIFEST"

echo "$THIS_SAMPLE: $NDONE / $TOTAL files verified"
if [ "$NDONE" -lt "$TOTAL" ]; then
    echo "some files missing, rerun to retry"
    exit 1
fi
echo "all files present and verified"
exit 0
