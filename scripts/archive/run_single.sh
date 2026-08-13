#!/bin/bash
# Single-sample hla_resolve run. Fill in the swap block, then sbatch.
# Two modes: MERGE (multiple same-library runs) or SINGLE (one deep run, no merge).

#SBATCH --job-name=hla_SAMPLE
#SBATCH --output=SAMPLE_%j.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --mem=60G
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00

set -eo pipefail
source ~/hla_env.sh

# ===== SWAP PER SAMPLE =====
SAMPLE=HG00438
BUCKET=HPRC                 # HPRC or HPRC_PLUS
FILES="m64043_200710_174426.ccs.bam m64043_200711_235708.ccs.bam"   # one or many
# ===========================

WORK=/hb/scratch/jsheh/hla
DATA=$WORK/data/$SAMPLE
mkdir -p "$DATA" "$WORK/results/$SAMPLE"
cd "$DATA"

B=s3://human-pangenomics/working/$BUCKET/$SAMPLE/raw_data/PacBio_HiFi
for f in $FILES ; do
    [ -s "$f" ] || aws s3 cp --no-sign-request "$B/$f" .
done

# --- MERGE mode (multiple same-library files): ---
samtools merge -@ 24 -f $SAMPLE.merged.hifi.bam $FILES
INPUT="$DATA/$SAMPLE.merged.hifi.bam"

# --- SINGLE mode (one deep run, no merge): comment out the two lines above
#     and set INPUT directly, e.g.:
# INPUT="$DATA/$(echo $FILES | awk '{print $1}')"

hla_resolve \
  --input_file "$INPUT" \
  --sample_name $SAMPLE \
  --platform pacbio \
  --scheme WGS \
  --output_dir "$WORK/results/$SAMPLE" \
  --threads 24
