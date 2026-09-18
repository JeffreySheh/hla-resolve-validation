#!/bin/bash
#SBATCH --job-name=hla_type
#SBATCH --output=/mnt/beegfs/cluster/scratch/jsheh/hla/pacbio_bam/logs/type_%A_%a.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --requeue
#SBATCH --cpus-per-task=12
#SBATCH --mem=30G
#SBATCH --time=2-00:00:00
# STEP 3 of 3: run hla_resolve on each merged bam. Reads any SAMPLE_LIST.
#   SAMPLE_LIST=~/samples_17.tsv sbatch --array=1-17 ~/3_type.sbatch

[ -f /software/scripts/utility/requeue.feature ] && source /software/scripts/utility/requeue.feature
set -eo pipefail
mkdir -p /mnt/beegfs/cluster/scratch/jsheh/hla/singularity_cache

WORK=/mnt/beegfs/cluster/scratch/jsheh/hla
MERGEDIR=${MERGEDIR:-$WORK/pacbio_bam/merged}
RESULTS=${RESULTS:-$WORK/results_1kg}
SAMPLE_LIST=${SAMPLE_LIST:-$HOME/samples_17.tsv}
[ -f "$SAMPLE_LIST" ] || { echo "no sample list: $SAMPLE_LIST"; exit 1; }

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '[:space:]')
[ -n "$SAMPLE" ] || exit 0

INPUT="$MERGEDIR/$SAMPLE.hifi.bam"
[ -s "$INPUT" ] || { echo "ERROR: no merged bam for $SAMPLE (run step 2 first): $INPUT"; exit 1; }

if [ -d "$RESULTS/$SAMPLE" ] && [ -n "$(ls -A "$RESULTS/$SAMPLE" 2>/dev/null)" ]; then
    echo "$SAMPLE already typed, skipping"; exit 0
fi

samtools view -c "$INPUT" >/dev/null 2>&1 || { echo "ERROR: $INPUT unreadable"; exit 1; }

mkdir -p "$RESULTS/$SAMPLE"
echo "typing $SAMPLE from $INPUT"
hla_resolve \
  --input_file "$INPUT" \
  --sample_name "$SAMPLE" \
  --platform pacbio \
  --scheme WGS \
  --output_dir "$RESULTS/$SAMPLE" \
  --threads 12 \
  --clean_up

echo "$SAMPLE: typing complete"
