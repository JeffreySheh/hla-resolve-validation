#!/bin/bash
#SBATCH --job-name=download_hprc
#SBATCH --output=/hb/scratch/mglasena/pacbio_bam/logs/pb_dl_%A_%a.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --requeue
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --array=1-17

# Download 69 HPRC PacBio HiFi ccs.bam files (1.98 TB) onto scratch, in parallel.
#
# Safe to rerun any number of times. Verified files are skipped, so a rerun after
# a timeout or failure only fetches what is still missing.
#
# The logs directory must exist before you submit, since Slurm opens the output
# file itself:  mkdir -p /hb/scratch/mglasena/pacbio_bam/logs
#
#   sbatch download_pacbio.sbatch                      fetch everything outstanding
#   STATUS_ONLY=1 bash download_pacbio.sbatch          report progress, download nothing
#   SAMPLES="HG00438 NA19240" sbatch download_pacbio.sbatch    restrict to some samples
#   FORCE=1 SAMPLES=HG00438 sbatch download_pacbio.sbatch      redownload even if done

[ -f /software/scripts/utility/requeue.feature ] && source /software/scripts/utility/requeue.feature
set -uo pipefail

export OUTDIR=${OUTDIR:-/hb/scratch/mglasena/pacbio_bam}

# Find the manifest and worker next to wherever you keep them, no env vars needed.
find_beside() {
    local name=$1 d
    for d in "$OUTDIR" "${SLURM_SUBMIT_DIR:-$PWD}" "$PWD" "$HOME"; do
        [ -f "$d/$name" ] && { echo "$d/$name"; return 0; }
    done
    return 1
}
MANIFEST=${MANIFEST:-$(find_beside pacbio_manifest.tsv)} || { echo "cannot find pacbio_manifest.tsv"; exit 1; }
WORKER=${WORKER:-$(find_beside fetch_bam.sh)} || { echo "cannot find fetch_bam.sh"; exit 1; }

# Compute nodes have a 1 GbE external NIC (~125 MB/s), and one stream already pulls
# ~45 MB/s, so 3 to 4 concurrent files saturate a node. More than that just adds
# contention. Throughput scales by adding nodes, not by adding streams per node.
NJOBS=${NJOBS:-4}         # files in flight
S3_CONC=${S3_CONC:-8}     # parallel range requests per file
SAMPLES=${SAMPLES:-}      # optional space separated whitelist
STATUS_ONLY=${STATUS_ONLY:-0}

export FORCE=${FORCE:-0}
export VERIFY_MD5=${VERIFY_MD5:-0}   # size + BGZF check is enough; set 1 for full md5

# aws lives behind a module, and batch shells do not get the module function for free.
if ! command -v aws >/dev/null 2>&1; then
    if ! command -v module >/dev/null 2>&1; then
        for init in /usr/share/lmod/lmod/init/bash /etc/profile.d/lmod.sh /etc/profile.d/modules.sh; do
            [ -f "$init" ] && { source "$init"; break; }
        done
    fi
    module load "${AWS_MODULE:-aws}" 2>/dev/null || module load aws/2.13.15 2>/dev/null
fi
command -v aws >/dev/null 2>&1 || { echo "FATAL: aws not on PATH and module load failed"; exit 1; }
echo "aws:      $(command -v aws)"

# samtools is optional. Without it the worker falls back to checking the BGZF EOF
# marker directly, which catches truncation just as well.
export QUICKCHECK=1
command -v samtools >/dev/null 2>&1 || module load "${SAMTOOLS_MODULE:-samtools}" 2>/dev/null
if command -v samtools >/dev/null 2>&1; then
    echo "samtools: $(command -v samtools)"
else
    echo "samtools: not found, using BGZF EOF check instead"
    export QUICKCHECK=0
fi

mkdir -p "$OUTDIR"
echo "manifest: $MANIFEST"
echo "worker:   $WORKER"

# In array mode each task takes one sample, so tasks spread across nodes and each gets
# its own NIC. Set SAMPLES explicitly to override.
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && [ -z "$SAMPLES" ]; then
    SAMPLES=$(cut -f1 "$MANIFEST" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
    [ -n "$SAMPLES" ] || { echo "array index ${SLURM_ARRAY_TASK_ID} is past the end of the sample list"; exit 0; }
    echo "array task ${SLURM_ARRAY_TASK_ID} -> $SAMPLES"
fi

# Apply the optional sample filter.
WORK_MANIFEST=$MANIFEST
if [ -n "$SAMPLES" ]; then
    WORK_MANIFEST=$OUTDIR/.manifest_subset.$$
    awk -F'\t' -v want="$SAMPLES" 'BEGIN{n=split(want,a," ");for(i=1;i<=n;i++)k[a[i]]=1} k[$1]' \
        "$MANIFEST" > "$WORK_MANIFEST"
    trap 'rm -f "$WORK_MANIFEST"' EXIT
    echo "filtered to samples: $SAMPLES"
fi

TOTAL=$(wc -l < "$WORK_MANIFEST")
[ "$TOTAL" -gt 0 ] || { echo "no files selected"; exit 1; }

is_done() { [ -f "$OUTDIR/$1/.$2.done" ] && [ -s "$OUTDIR/$1/$2" ]; }

report() {
    local ndone=0 nmiss=0 bytes=0
    while IFS=$'\t' read -r s f u sz; do
        if is_done "$s" "$f"; then
            ndone=$(( ndone + 1 )); bytes=$(( bytes + sz ))
        else
            nmiss=$(( nmiss + 1 )); [ "${1:-}" = verbose ] && echo "  outstanding: $s/$f"
        fi
    done < "$WORK_MANIFEST"
    echo "complete: $ndone / $TOTAL   outstanding: $nmiss   ($(( bytes / 1000000000 )) GB on disk)"
    return "$nmiss"
}

echo "=== before ==="
report
df -h "$OUTDIR" | tail -1

if [ "$STATUS_ONLY" = 1 ]; then
    report verbose
    exit 0
fi

# Clear partial files from a previous run so nothing stale is adopted. Scoped to this
# job's own manifest, so sharding across nodes with SAMPLES does not let one job delete
# another's in-flight downloads. The shards must not overlap.
# The glob matters: `aws s3 cp` downloads into its own temp name <dest>.<random>, so a
# cancelled job leaves $f.part.XXXXXXXX behind. Matching only $f.part orphans those
# forever and they are tens of GB each.
while IFS=$'\t' read -r s f u sz; do
    rm -f "$OUTDIR/$s/$f".part*
done < "$WORK_MANIFEST"

echo "concurrency: $NJOBS files x $S3_CONC requests"

# Job-local aws config so we do not touch ~/.aws
export AWS_CONFIG_FILE=$OUTDIR/.aws_config_job
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

START=$SECONDS
awk -F'\t' '{print $1, $2, $3, $4}' "$WORK_MANIFEST" \
  | xargs -P "$NJOBS" -n 4 bash "$WORKER"
ELAPSED=$(( SECONDS - START ))

echo
echo "=== after (${ELAPSED}s) ==="
report verbose
NMISS=$?
du -sh "$OUTDIR" 2>/dev/null

if [ "$NMISS" -gt 0 ]; then
    echo "rerun this script to retry the outstanding files"
    exit 1
fi
echo "all files present and verified"
exit 0
