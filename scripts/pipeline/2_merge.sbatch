#!/bin/bash
#SBATCH --job-name=hla_merge
#SBATCH --output=/hb/scratch/jsheh/hla/pacbio_bam/logs/merge_%A_%a.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --requeue
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
# STEP 2: pick ONE chemistry per sample, combine -> one bam. Reads any SAMPLE_LIST.
#   SAMPLE_LIST=~/samples_17.tsv sbatch --array=1-17 ~/2_merge.sbatch

set -uo pipefail
OUTDIR=${OUTDIR:-/hb/scratch/jsheh/hla/pacbio_bam}
MERGEDIR=$OUTDIR/merged
REVIO_MIN=${REVIO_MIN:-3}
FORCE_SQ2=${FORCE_SQ2:-" HG01952 HG01978 "}
SAMPLE_LIST=${SAMPLE_LIST:-$HOME/samples_17.tsv}
[ -f "$SAMPLE_LIST" ] || { echo "no sample list: $SAMPLE_LIST"; exit 1; }
mkdir -p "$MERGEDIR"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '[:space:]')
[ -n "$SAMPLE" ] || exit 0
command -v samtools >/dev/null || { echo "FATAL: no samtools (source ~/hla_env.sh)"; exit 1; }

OUT="$MERGEDIR/$SAMPLE.hifi.bam"; DONE="$MERGEDIR/.$SAMPLE.done"; TMP="$OUT.tmp"; HDR="$MERGEDIR/.$SAMPLE.hdr"
log() { echo "[$(date +%H:%M:%S)] $SAMPLE: $*"; }
[ -f "$DONE" ] && [ -s "$OUT" ] && { log "already done"; exit 0; }
rm -f "$DONE" "$TMP" "$HDR"

REVIO=(); SQ2=()
for f in "$OUTDIR/$SAMPLE"/*.bam; do
    [ -s "$f" ] && [ -f "$OUTDIR/$SAMPLE/.$(basename "$f").done" ] || continue
    case "$(basename "$f")" in m84*) REVIO+=("$f");; m64*) SQ2+=("$f");; esac
done

if [[ "$FORCE_SQ2" == *" $SAMPLE "* ]] && [ ${#SQ2[@]} -gt 0 ]; then CHEM="SequelII(forced)"; INPUTS=("${SQ2[@]}")
elif [ ${#REVIO[@]} -ge "$REVIO_MIN" ]; then CHEM=Revio; INPUTS=("${REVIO[@]}")
elif [ ${#SQ2[@]} -gt 0 ]; then CHEM=SequelII; INPUTS=("${SQ2[@]}")
elif [ ${#REVIO[@]} -gt 0 ]; then CHEM="Revio(thin)"; INPUTS=("${REVIO[@]}")
else log "no usable files"; exit 1; fi
NIN=${#INPUTS[@]}
log "$CHEM: ${#REVIO[@]} Revio / ${#SQ2[@]} SqII -> using $NIN"

KIN=$(samtools view "${INPUTS[0]}" 2>/dev/null | head -500 | grep -oE '(fi|fp|ri|rp|ip|pw):' | head -1)

if [ "$NIN" -eq 1 ] && [ -z "$KIN" ]; then
    cp "${INPUTS[0]}" "$TMP"
elif [ -n "$KIN" ]; then
    log "stripping kinetics during merge"
    samtools merge -@8 -f - "${INPUTS[@]}" | samtools view -b -@8 -x fi -x fp -x ri -x rp -x ip -x pw -o "$TMP" -
else
    { samtools view -H "${INPUTS[0]}" | grep '^@HD'
      for f in "${INPUTS[@]}"; do samtools view -H "$f" | grep '^@RG'; done | awk '!s[$0]++'
      samtools view -H "${INPUTS[0]}" | grep -vE '^@(HD|RG)'; } > "$HDR"
    if [ "$(grep -c '^@RG' "$HDR")" = "$NIN" ]; then
        samtools cat -h "$HDR" -o "$TMP" "${INPUTS[@]}"
    else
        samtools merge -@8 -f "$TMP" "${INPUTS[@]}"
    fi
fi
[ -s "$TMP" ] || { log "merge produced nothing"; exit 1; }

samtools quickcheck -u "$TMP" || { log "quickcheck FAILED"; rm -f "$TMP"; exit 1; }
[ "$(samtools view -H "$TMP" | grep -c '^@RG')" = "$NIN" ] || { log "read-group count wrong"; rm -f "$TMP"; exit 1; }

mv -f "$TMP" "$OUT"; rm -f "$HDR"; touch "$DONE"
log "ok: $CHEM, $NIN files -> $(du -h --apparent-size "$OUT" | cut -f1)"
