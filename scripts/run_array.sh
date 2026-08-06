#!/bin/bash
#SBATCH --job-name=hla_array
#SBATCH --output=hla_array_%A_%a.out
#SBATCH --partition=lab-colibri
#SBATCH --account=pi-jkoc
#SBATCH --qos=pi-jkoc
#SBATCH --mem=30G
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH --array=1-17%2

set -eo pipefail
source ~/hla_env.sh

WORK=/hb/scratch/jsheh/hla
MANIFEST=data/samples.tsv

# pull THIS task's line from the manifest by array index
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
SAMPLE=$(echo "$LINE" | cut -f1)
BUCKET=$(echo "$LINE" | cut -f2)
FILES=$(echo "$LINE" | cut -f3)

DATA=$WORK/data/$SAMPLE
mkdir -p "$DATA" "$WORK/results/$SAMPLE"
cd "$DATA"

# download this sample's files from S3
B=s3://human-pangenomics/working/$BUCKET/$SAMPLE/raw_data/PacBio_HiFi
for f in $FILES ; do
    [ -s "$f" ] || aws s3 cp --no-sign-request "$B/$f" .
done

# merge same-library runs, then type
samtools merge -@ 12 -f $SAMPLE.merged.hifi.bam $FILES
samtools view -c -@ 12 $SAMPLE.merged.hifi.bam

hla_resolve \
  --input_file "$DATA/$SAMPLE.merged.hifi.bam" \
  --sample_name $SAMPLE \
  --platform pacbio \
  --scheme WGS \
  --output_dir "$WORK/results/$SAMPLE" \
  --threads 12
