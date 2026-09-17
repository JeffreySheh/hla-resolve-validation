#!/bin/bash
#SBATCH --job-name=hla_merge
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --no-requeue
# STEP 2: merge exactly the files data/merge_plan.tsv lists for
# one input BAM (a sample x library prep, named <sample>__<prep>) into one unaligned
# BAM with samtools cat. No chemistry rule here;
# the decision was made by scripts/select_inputs.py and is reviewable in the
# plan. Input files are looked up in every root in SEARCH_ROOTS (read-only);
# output goes under MERGEDIR only.
#
#   PLAN=.../data/merge_plan.tsv MERGEDIR=/mnt/beegfs/.../merged \
#   SEARCH_ROOTS="/mnt/beegfs/.../mglasena/pacbio_bam:/mnt/beegfs/.../jsheh/hla/pacbio_bam" \
#   sbatch --array=1-N --output=.../logs/merge_%A_%a.out 2_merge.sbatch
#
# Array index N maps to the Nth input BAM in the plan (sorted, unique). Output is
# <input_bam>.hifi.bam; type it with --sample_name <sample> and an output dir
# named after it. Use
# absolute /mnt/beegfs paths: /hb/scratch symlinks are not visible on nodes.

set -uo pipefail
PLAN=${PLAN:?set PLAN to merge_plan.tsv}
MERGEDIR=${MERGEDIR:?set MERGEDIR}
SEARCH_ROOTS=${SEARCH_ROOTS:?set SEARCH_ROOTS, colon-separated}
mkdir -p "$MERGEDIR"
command -v samtools >/dev/null 2>&1 || module load "${SAMTOOLS_MODULE:-samtools}" 2>/dev/null
command -v samtools >/dev/null 2>&1 || { echo "FATAL: no samtools"; exit 1; }

ENTRY=$(tail -n +2 "$PLAN" | cut -f1 | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$ENTRY" ] || { echo "array index past end"; exit 0; }
SAMPLE=${ENTRY%%__*}
log() { echo "[$(date +%H:%M:%S)] $ENTRY: $*"; }

OUT="$MERGEDIR/$ENTRY.hifi.bam"; DONE="$MERGEDIR/.$ENTRY.done"
TMP="$OUT.tmp"; HDR="$MERGEDIR/.$ENTRY.hdr"; LIST="$MERGEDIR/.$ENTRY.inputs"
[ -f "$DONE" ] && [ -s "$OUT" ] && { log "already done"; exit 0; }
rm -f "$DONE" "$TMP" "$HDR"

# resolve every planned file against the search roots; a verified download
# has a .done marker beside it
INPUTS=(); : > "$LIST"
while IFS=$'\t' read -r entry s prep fname rest; do
    [ "$entry" = "$ENTRY" ] || continue
    found=""
    IFS=: read -ra ROOTS <<< "$SEARCH_ROOTS"
    for root in "${ROOTS[@]}"; do
        f="$root/$SAMPLE/$fname"
        if [ -s "$f" ] && [ -f "$root/$SAMPLE/.$fname.done" ]; then found="$f"; break; fi
    done
    [ -n "$found" ] || { log "MISSING or unverified: $fname"; exit 1; }
    INPUTS+=("$found"); echo "$found" >> "$LIST"
done < "$PLAN"
NIN=${#INPUTS[@]}
[ "$NIN" -gt 0 ] || { log "no files in plan"; exit 1; }
log "$NIN input(s) from plan"

if [ "$NIN" -eq 1 ]; then
    ln -sf "${INPUTS[0]}" "$OUT"
    log "single file: symlinked"
else
    # one @HD, every @RG, the rest of the header from the first file
    { samtools view -H "${INPUTS[0]}" | grep '^@HD'
      for f in "${INPUTS[@]}"; do samtools view -H "$f" | grep '^@RG'; done | awk '!s[$0]++'
      samtools view -H "${INPUTS[0]}" | grep -vE '^@(HD|RG)'; } > "$HDR"
    NRG=$(grep -c '^@RG' "$HDR")
    [ "$NRG" -ge "$NIN" ] || { log "header has $NRG @RG for $NIN inputs"; rm -f "$HDR"; exit 1; }
    samtools cat -h "$HDR" -o "$TMP" "${INPUTS[@]}" || { log "samtools cat failed"; rm -f "$TMP"; exit 1; }
    samtools quickcheck -u "$TMP" || { log "quickcheck FAILED"; rm -f "$TMP"; exit 1; }
    mv -f "$TMP" "$OUT"; rm -f "$HDR"
fi
touch "$DONE"
log "ok: $NIN files -> $(du -h --apparent-size -L "$OUT" | cut -f1)"
