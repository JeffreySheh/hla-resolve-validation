#!/bin/bash
#SBATCH --job-name=hla_cleanup
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
#SBATCH --no-requeue
# STEP 6: after a batch is typed, delete the inputs of every input BAM that
# typed with no error, and keep the inputs of every one that did not.
#
#   PLAN=.../stage4/b02.merge_plan.tsv MERGEDIR=.../merged RAWDIR=.../raw \
#   RESULTS=.../results REPORT=.../b02.cleanup.tsv sbatch 6_cleanup.sbatch
#
# An input BAM typed with no error when all three are true:
#   $RESULTS/<input_bam>/.typed_ok exists (3_type.sbatch writes it on exit 0)
#   allele_output.csv is not empty
#   the hla_resolve log ends the run with "(status: ok)"
# Then the merged BAM (or its symlink) and the planned runs under RAWDIR are
# deleted, and .<input_bam>.cleaned is left in MERGEDIR. Results are kept.
# REPORT gets one line per input BAM: input_bam, action, reason.

set -uo pipefail
PLAN=${PLAN:?set PLAN}
MERGEDIR=${MERGEDIR:?set MERGEDIR}
RAWDIR=${RAWDIR:?set RAWDIR}
RESULTS=${RESULTS:?set RESULTS}
REPORT=${REPORT:?set REPORT}
SUFFIX=${SUFFIX:-ont}
for d in "$MERGEDIR" "$RAWDIR" "$RESULTS"; do
    [ -d "$d" ] && [ "$(readlink -f "$d")" != / ] || { echo "FATAL: bad directory '$d'"; exit 1; }
done
log() { echo "[$(date +%H:%M:%S)] $*"; }

: > "$REPORT"
NKEPT=0
for ENTRY in $(tail -n +2 "$PLAN" | cut -f1 | sort -u); do
    SAMPLE=${ENTRY%%__*}
    OUT="$RESULTS/$ENTRY"
    CSV="$OUT/$SAMPLE/hla_typing_results/allele_output.csv"
    LOGF="$OUT/$SAMPLE.hla_resolve.log"
    why=""
    [ -f "$OUT/.typed_ok" ] || why="no .typed_ok (typing failed, cancelled or never ran)"
    [ -z "$why" ] && [ ! -s "$CSV" ] && why="allele_output.csv missing or empty"
    [ -z "$why" ] && ! grep -q "(status: ok)" "$LOGF" 2>/dev/null && why="log has no '(status: ok)'"
    if [ -n "$why" ]; then
        printf '%s\tkept\t%s\n' "$ENTRY" "$why" >> "$REPORT"
        log "$ENTRY: kept, $why"
        NKEPT=$(( NKEPT + 1 ))
        continue
    fi

    rm -f "$MERGEDIR/$ENTRY.$SUFFIX.bam" "$MERGEDIR/.$ENTRY.done" "$MERGEDIR/.$ENTRY.inputs" "$MERGEDIR/.$ENTRY.hdr"
    n=0
    while IFS=$'\t' read -r entry s prep fname rest; do
        [ "$entry" = "$ENTRY" ] || continue
        f="$RAWDIR/$SAMPLE/$fname"
        [ -e "$f" ] && n=$(( n + 1 ))
        rm -f "$f" "$RAWDIR/$SAMPLE/.$fname.done"
    done < "$PLAN"
    rmdir "$RAWDIR/$SAMPLE" 2>/dev/null
    touch "$MERGEDIR/.$ENTRY.cleaned"
    printf '%s\tdeleted\tmerged BAM and %d raw run(s)\n' "$ENTRY" "$n" >> "$REPORT"
    log "$ENTRY: deleted merged BAM and $n raw run(s)"
done
log "done: $NKEPT kept, report in $REPORT"
