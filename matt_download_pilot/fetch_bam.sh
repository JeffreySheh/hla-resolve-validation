#!/bin/bash
# Worker: download one HiFi ccs.bam from S3, verify, mark done.
# Called by download_pacbio.sbatch via xargs. Args: sample filename s3_uri expected_bytes
set -uo pipefail

SAMPLE=$1
FNAME=$2
URI=$3
EXPECT=$4

DEST="$OUTDIR/$SAMPLE"
OUT="$DEST/$FNAME"
DONE="$DEST/.$FNAME.done"
TMP="$OUT.part"

log() { printf '[%s] %s %s: %s\n' "$(date +%H:%M:%S)" "$SAMPLE" "$FNAME" "$*"; }

# Every complete BGZF file ends in this 28 byte empty-block marker, so its absence
# means truncation. Cheap second guard alongside quickcheck.
BGZF_EOF=1f8b08040000000000ff0600424302001b0003000000000000000000

verify() {
    local f=$1 sz tail_hex want got
    sz=$(stat -c %s "$f" 2>/dev/null) || return 1
    [ "$sz" = "$EXPECT" ] || { log "size $sz != expected $EXPECT"; return 1; }
    # -u is required: these are unaligned ccs.bam with no @SQ lines, and without it
    # quickcheck returns 8 ("no targets in header") on every single one of them.
    if [ "$QUICKCHECK" = 1 ]; then
        samtools quickcheck -u "$f" || { log "quickcheck failed"; return 1; }
    fi
    tail_hex=$(tail -c 28 "$f" | od -An -tx1 -v | tr -d ' \n')
    [ "$tail_hex" = "$BGZF_EOF" ] || { log "BGZF EOF marker missing, truncated"; return 1; }
    if [ "$VERIFY_MD5" = 1 ]; then
        want=$(aws s3 cp --no-sign-request --quiet "$URI.md5" - 2>/dev/null | awk '{print $1}')
        got=$(md5sum "$f" | awk '{print $1}')
        [ -n "$want" ] && [ "$want" != "$got" ] && { log "md5 $got != $want"; return 1; }
    fi
    return 0
}

[ "${FORCE:-0}" = 1 ] && rm -f "$DONE" "$OUT"

# A done marker only counts if the file it vouches for is still there.
if [ -f "$DONE" ] && [ -s "$OUT" ]; then
    exit 0
fi
rm -f "$DONE"

mkdir -p "$DEST"

# Adopt a complete file left behind by an earlier run.
if [ -s "$OUT" ] && verify "$OUT"; then
    touch "$DONE"
    log "already present, verified"
    exit 0
fi

for attempt in 1 2 3; do
    rm -f "$TMP"
    log "download attempt $attempt"
    if aws s3 cp --no-sign-request --only-show-errors "$URI" "$TMP" && verify "$TMP"; then
        mv -f "$TMP" "$OUT"
        touch "$DONE"
        log "ok ($(( EXPECT / 1000000 )) MB)"
        exit 0
    fi
    sleep $(( attempt * 20 ))
done

rm -f "$TMP"
log "FAILED after 3 attempts"
exit 1
