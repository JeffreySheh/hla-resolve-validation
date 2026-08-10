#!/bin/bash
#SBATCH --job-name=merge_hprc
#SBATCH --output=/hb/scratch/mglasena/pacbio_bam/logs/merge_%A_%a.out
#SBATCH --partition=windfall
#SBATCH --account=windfall
#SBATCH --qos=windfall
#SBATCH --requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-17
#
# Combine each sample's per-movie ccs.bam into one unaligned bam for hla_resolve.
#
# Uses `samtools cat` with a constructed header rather than `samtools merge`.
#
#   `samtools merge` decompresses and recompresses every record. Measured on this
#   cluster it ran at ~15 MB/s per task, which is a 4 to 5 hour job for 2 TB.
#
#   `samtools cat` copies compressed BGZF blocks untouched, so it runs at I/O speed.
#   Its catch is that it takes the header from the first file only, which would drop
#   the @RG lines of every later movie and leave those reads pointing at read groups
#   that do not exist. Passing -h with a header containing all of them fixes that.
#
# Verified equivalent on real HG00438 data: identical read names, sequences,
# qualities, and read-to-read-group mapping. Only tag order differs, which is
# not meaningful in BAM.
#
# Safe to rerun. Verified outputs are skipped.
#
#   sbatch merge_bams.sbatch                       combine everything outstanding
#   STATUS_ONLY=1 bash merge_bams.sbatch           report, do nothing
#   SAMPLES="HG00438" sbatch merge_bams.sbatch     one sample
#   DELETE_INPUTS=1 sbatch merge_bams.sbatch       remove per-movie files after verify
#   VERIFY_COUNTS=1 sbatch merge_bams.sbatch       also reconcile read counts (slow)
#   USE_MERGE=1 sbatch merge_bams.sbatch           fall back to samtools merge

[ -f /software/scripts/utility/requeue.feature ] && source /software/scripts/utility/requeue.feature
set -uo pipefail

export OUTDIR=${OUTDIR:-/hb/scratch/mglasena/pacbio_bam}
MERGEDIR=${MERGEDIR:-$OUTDIR/merged}

MANIFEST=""
for d in "$OUTDIR" "${SLURM_SUBMIT_DIR:-$PWD}" "$PWD" "$HOME"; do
    [ -f "$d/pacbio_manifest.tsv" ] && { MANIFEST="$d/pacbio_manifest.tsv"; break; }
done
[ -n "$MANIFEST" ] || { echo "cannot find pacbio_manifest.tsv"; exit 1; }

THREADS=${THREADS:-${SLURM_CPUS_PER_TASK:-8}}
SAMPLES=${SAMPLES:-}
STATUS_ONLY=${STATUS_ONLY:-0}
DELETE_INPUTS=${DELETE_INPUTS:-0}
VERIFY_COUNTS=${VERIFY_COUNTS:-0}
USE_MERGE=${USE_MERGE:-0}

mkdir -p "$MERGEDIR"

BGZF_EOF=1f8b08040000000000ff0600424302001b0003000000000000000000

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ] && [ -z "$SAMPLES" ]; then
    SAMPLES=$(cut -f1 "$MANIFEST" | sort -u | sed -n "${SLURM_ARRAY_TASK_ID}p")
    [ -n "$SAMPLES" ] || { echo "array index past end of sample list"; exit 0; }
fi
[ -n "$SAMPLES" ] || SAMPLES=$(cut -f1 "$MANIFEST" | sort -u)

if [ "$STATUS_ONLY" = 1 ]; then
    n=0; ok=0
    for s in $SAMPLES; do
        n=$((n+1))
        if [ -f "$MERGEDIR/.$s.done" ] && [ -s "$MERGEDIR/$s.hifi.bam" ]; then
            ok=$((ok+1))
            printf '  %-9s DONE   %s\n' "$s" "$(du -h --apparent-size "$MERGEDIR/$s.hifi.bam" 2>/dev/null | cut -f1)"
        elif [ -s "$MERGEDIR/$s.hifi.bam.tmp" ]; then
            printf '  %-9s running %s so far\n' "$s" "$(du -h --apparent-size "$MERGEDIR/$s.hifi.bam.tmp" 2>/dev/null | cut -f1)"
        else
            printf '  %-9s pending\n' "$s"
        fi
    done
    echo "combined: $ok / $n"
    [ "$ok" = "$n" ]
    exit $?
fi

# Only the work path needs samtools. STATUS_ONLY above deliberately runs without it,
# so progress can be checked from a plain shell with no conda env active.
command -v samtools >/dev/null 2>&1 || { echo "FATAL: samtools not on PATH"; exit 1; }

RC=0
for SAMPLE in $SAMPLES; do
    OUT="$MERGEDIR/$SAMPLE.hifi.bam"
    DONE="$MERGEDIR/.$SAMPLE.done"
    TMP="$OUT.tmp"
    HDR="$MERGEDIR/.$SAMPLE.header.sam"

    log() { printf '[%s] %s: %s\n' "$(date +%H:%M:%S)" "$SAMPLE" "$*"; }

    if [ -f "$DONE" ] && [ -s "$OUT" ]; then
        log "already combined, skipping"
        continue
    fi
    rm -f "$DONE" "$TMP" "$HDR"

    INPUTS=()
    IN_BYTES=0
    MISSING=0
    while IFS=$'\t' read -r s f u sz; do
        [ "$s" = "$SAMPLE" ] || continue
        if [ -f "$OUTDIR/$s/.$f.done" ] && [ -s "$OUTDIR/$s/$f" ]; then
            INPUTS+=("$OUTDIR/$s/$f")
            IN_BYTES=$(( IN_BYTES + sz ))
        else
            log "input not downloaded or unverified: $f"; MISSING=1
        fi
    done < "$MANIFEST"
    [ "$MISSING" = 0 ] || { log "SKIP, inputs incomplete"; RC=1; continue; }
    NIN=${#INPUTS[@]}

    # Kinetics roughly double file size and nothing downstream uses them. The HPRC raw
    # release has them stripped already, so this normally finds nothing. Note that
    # stripping requires decoding every record, so it forces the slow merge path.
    KIN=""
    for f in "${INPUTS[@]}"; do
        found=$(samtools view "$f" 2>/dev/null | head -1000 | awk '
            { for(i=12;i<=NF;i++){ split($i,a,":"); t=a[1]
                if(t=="fi"||t=="fp"||t=="ri"||t=="rp"||t=="ipd"||t=="pw") k[t]=1 } }
            END { for(t in k) printf "%s ", t }')
        [ -n "$found" ] && { KIN="$found"; break; }
    done
    [ -n "$KIN" ] && log "kinetics tags present ($KIN), forcing decode path to strip them" \
                  || log "no kinetics tags found, nothing to strip"

    if [ "$VERIFY_COUNTS" = 1 ]; then
        IN_READS=0
        for f in "${INPUTS[@]}"; do
            c=$(samtools view -c -@ "$THREADS" "$f") || c=-1
            [ "$c" -lt 0 ] && { IN_READS=-1; break; }
            IN_READS=$(( IN_READS + c ))
        done
        log "input reads: $IN_READS"
    fi

    if [ "$NIN" -eq 1 ] && [ -z "$KIN" ]; then
        log "single movie, copying"
        cp "${INPUTS[0]}" "$TMP" || { log "copy FAILED"; rm -f "$TMP"; RC=1; continue; }

    elif [ -n "$KIN" ] || [ "$USE_MERGE" = 1 ]; then
        log "combining $NIN movies with samtools merge (decode path)"
        if [ -n "$KIN" ]; then
            samtools merge -@ "$THREADS" -f - "${INPUTS[@]}" \
              | samtools view -b -@ "$THREADS" -x fi -x fp -x ri -x rp -x ipd -x pw -o "$TMP" - \
              || { log "merge+strip FAILED"; rm -f "$TMP"; RC=1; continue; }
        else
            samtools merge -@ "$THREADS" -f "$TMP" "${INPUTS[@]}" \
              || { log "merge FAILED"; rm -f "$TMP"; RC=1; continue; }
        fi

    else
        # Header: @HD from the first file, then every distinct @RG across all movies,
        # then the remaining header lines from the first file.
        {
            samtools view -H "${INPUTS[0]}" | grep '^@HD' || true
            for f in "${INPUTS[@]}"; do samtools view -H "$f" | grep '^@RG' || true; done | awk '!seen[$0]++'
            samtools view -H "${INPUTS[0]}" | grep -v '^@HD' | grep -v '^@RG' || true
        } > "$HDR"
        nrg_hdr=$(grep -c '^@RG' "$HDR")
        if [ "$nrg_hdr" != "$NIN" ]; then
            log "header build produced $nrg_hdr @RG for $NIN movies, falling back to samtools merge"
            samtools merge -@ "$THREADS" -f "$TMP" "${INPUTS[@]}" \
              || { log "merge FAILED"; rm -f "$TMP"; RC=1; continue; }
        else
            log "combining $NIN movies with samtools cat, $nrg_hdr read groups in header"
            samtools cat -h "$HDR" -o "$TMP" "${INPUTS[@]}" \
              || { log "cat FAILED"; rm -f "$TMP"; RC=1; continue; }
        fi
    fi

    # verify
    fail=""
    samtools quickcheck -u "$TMP" || fail="quickcheck"
    if [ -z "$fail" ]; then
        hex=$(tail -c 28 "$TMP" | od -An -tx1 -v | tr -d ' \n')
        [ "$hex" = "$BGZF_EOF" ] || fail="bgzf-eof"
    fi
    if [ -z "$fail" ]; then
        nrg=$(samtools view -H "$TMP" | grep -c '^@RG')
        [ "$nrg" = "$NIN" ] || fail="read groups: header has $nrg, expected $NIN"
    fi
    # cat copies blocks verbatim, so the output should be within a hair of the inputs.
    # This catches a truncated copy without reading the data back.
    if [ -z "$fail" ] && [ -z "$KIN" ] && [ "$USE_MERGE" = 0 ] && [ "$NIN" -gt 1 ]; then
        outb=$(stat -c %s "$TMP")
        lo=$(( IN_BYTES / 100 * 97 ))
        [ "$outb" -ge "$lo" ] || fail="size: got $outb, expected near $IN_BYTES"
    fi
    if [ -z "$fail" ] && [ "$VERIFY_COUNTS" = 1 ] && [ "$IN_READS" -ge 0 ]; then
        out_reads=$(samtools view -c -@ "$THREADS" "$TMP")
        [ "$out_reads" = "$IN_READS" ] || fail="read count: got $out_reads, expected $IN_READS"
    fi

    if [ -n "$fail" ]; then
        log "VERIFY FAILED ($fail)"; rm -f "$TMP"; RC=1; continue
    fi

    mv -f "$TMP" "$OUT"
    rm -f "$HDR"
    touch "$DONE"
    log "ok, $NIN movies -> $(du -h --apparent-size "$OUT" | cut -f1), $NIN read groups preserved"

    if [ "$DELETE_INPUTS" = 1 ]; then
        for f in "${INPUTS[@]}"; do rm -f "$f"; done
        rm -f "$OUTDIR/$SAMPLE"/.*.done
        log "removed $NIN per-movie input files"
    fi
done

exit $RC
