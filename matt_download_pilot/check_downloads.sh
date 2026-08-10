#!/bin/bash
# Progress and health check for the HPRC PacBio download.
#
#   check_downloads.sh          one-shot report
#   check_downloads.sh -r       also measure per-node rate over 60s (spots slow NICs)
#   check_downloads.sh -w       refresh every 30s until everything is complete
#
# Reads only the manifest and the download directory, so it is safe to run at any
# time from anywhere, including while jobs are running.
set -uo pipefail

OUTDIR=${OUTDIR:-/hb/scratch/mglasena/pacbio_bam}
RATE=0
WATCH=0
CLEAN=0
while getopts "rwch" o; do
    case $o in
        r) RATE=1 ;;
        w) WATCH=1 ;;
        c) CLEAN=1 ;;
        h) sed -n '2,10p' "$0"; exit 0 ;;
        *) exit 1 ;;
    esac
done

MANIFEST=""
for d in "$OUTDIR" "$(dirname "$0")" "$PWD" "$HOME"; do
    [ -f "$d/pacbio_manifest.tsv" ] && { MANIFEST="$d/pacbio_manifest.tsv"; break; }
done
[ -n "$MANIFEST" ] || { echo "cannot find pacbio_manifest.tsv"; exit 1; }

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# Array task index -> sample, the same 1-based sorted order the sbatch script uses.
cut -f1 "$MANIFEST" | sort -u | awk '{print NR"\t"$0}' > "$TMP/idx"

# Report apparent size, not allocated blocks. Plain `du` under-reports on beegfs by a
# few percent, which makes verified samples look incomplete against the manifest.
DUFLAGS=""
du --apparent-size -sk "$TMP" >/dev/null 2>&1 && DUFLAGS="--apparent-size"
dir_kb() { du $DUFLAGS -sk "$1" 2>/dev/null | cut -f1; }

# sample -> node, resolved through the array task index
node_map() {
    squeue -u "$USER" -h -o '%i %N' 2>/dev/null | grep '_' \
      | awk '{n=$1; sub(/^.*_/,"",n); print n"\t"$2}' \
      | join -t"$(printf '\t')" -1 1 -2 1 -o 1.2,2.2 <(sort -k1,1 "$TMP/idx") - 2>/dev/null
}

report() {
    node_map | sort > "$TMP/nodes" 2>/dev/null || : > "$TMP/nodes"

    awk -F'\t' '{b[$1]+=$4; n[$1]++} END{for(s in b) printf "%s\t%s\t%s\n", s, b[s], n[s]}' \
        "$MANIFEST" | sort > "$TMP/targets"

    : > "$TMP/rows"
    while IFS=$'\t' read -r s tgt_b nfiles; do
        ndone=$(ls "$OUTDIR/$s"/.*.done 2>/dev/null | wc -l | tr -d ' ')
        kb=$(dir_kb "$OUTDIR/$s"); [ -n "$kb" ] || kb=0
        node=$(awk -F'\t' -v s="$s" '$1==s{print $2; exit}' "$TMP/nodes")
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "${node:--}" "$ndone" "$nfiles" "$kb" "$tgt_b" >> "$TMP/rows"
    done < "$TMP/targets"

    awk -F'\t' '
    BEGIN { printf "%-9s %-9s %8s %10s %10s %6s  %s\n", "SAMPLE","NODE","FILES","HAVE","TARGET","PCT","" }
    {
        have=$5*1024; tgt=$6
        pct=(tgt>0)? int(have*100/tgt) : 0; if (pct>100) pct=100
        bar=""; for(i=0;i<int(pct/5);i++) bar=bar"#"
        printf "%-9s %-9s %4s/%-3s %9.1fG %9.1fG %5d%%  %s\n", $1,$2,$3,$4, have/1073741824, tgt/1073741824, pct, bar
        d+=$3; n+=$4; H+=have; T+=tgt
    }
    END {
        printf "\nTOTAL     verified %d/%d files   %.1f / %.1f GiB   %d%%\n", d,n,H/1073741824,T/1073741824,(T>0)?H*100/T:0
        exit (d==n)?0:1
    }' "$TMP/rows"
    rc=$?

    nparts=$(find "$OUTDIR" -name '*.ccs.bam.part*' 2>/dev/null | wc -l | tr -d ' ')
    [ "$nparts" -gt 0 ] && echo "  .part files present (in flight or orphaned): $nparts"
    bad=$(grep -l 'FAILED\|attempt [23]\|!=\|No space\|truncated' "$OUTDIR"/logs/*.out 2>/dev/null)
    if [ -n "$bad" ]; then
        echo "  WARNING: retries or errors found in:"
        echo "$bad" | sed 's/^/    /'
    fi
    return $rc
}

rate_check() {
    echo
    echo "measuring per-node rate over 60s ..."
    : > "$TMP/before"
    while IFS=$'\t' read -r i s; do
        printf '%s\t%s\n' "$s" "$(dir_kb "$OUTDIR/$s")" >> "$TMP/before"
    done < "$TMP/idx"
    sleep 60
    node_map | sort > "$TMP/nodes2" 2>/dev/null || : > "$TMP/nodes2"
    printf '  %-9s %-9s %8s\n' SAMPLE NODE 'MB/s'
    total=0
    while IFS=$'\t' read -r s before; do
        now=$(dir_kb "$OUTDIR/$s"); [ -n "$now" ] || now=0
        [ -n "$before" ] || before=0
        mbs=$(( (now - before) * 1024 / 60 / 1000000 ))
        [ "$mbs" -le 0 ] && continue
        node=$(awk -F'\t' -v s="$s" '$1==s{print $2; exit}' "$TMP/nodes2")
        printf '  %-9s %-9s %7s\n' "$s" "${node:--}" "$mbs"
        total=$(( total + mbs ))
    done < "$TMP/before"
    echo "  aggregate: $total MB/s"
    echo "  ~100 MB/s per node is a saturated 1 GbE link and is what you want."
    echo "  ~11 MB/s means a 100 Mbit node, cancel that task and add the node to --exclude."
}

# Delete leftover aws temp files, but only for samples where every file is verified,
# so a running download is never touched.
cleanup() {
    echo
    echo "removing orphaned .part files from fully verified samples ..."
    awk -F'\t' '{n[$1]++} END{for(s in n) printf "%s\t%s\n", s, n[s]}' "$MANIFEST" | sort > "$TMP/cnt"
    freed_kb=0
    while IFS=$'\t' read -r s nfiles; do
        ndone=$(ls "$OUTDIR/$s"/.*.done 2>/dev/null | wc -l | tr -d ' ')
        npart=$(find "$OUTDIR/$s" -maxdepth 1 -name '*.ccs.bam.part*' 2>/dev/null | wc -l | tr -d ' ')
        [ "$npart" -gt 0 ] || continue
        if [ "$ndone" != "$nfiles" ]; then
            printf '  %-9s SKIP, still downloading (%s/%s verified, %s part files)\n' "$s" "$ndone" "$nfiles" "$npart"
            continue
        fi
        kb=$(find "$OUTDIR/$s" -maxdepth 1 -name '*.ccs.bam.part*' -exec du -sk {} + 2>/dev/null | awk '{t+=$1} END{print t+0}')
        find "$OUTDIR/$s" -maxdepth 1 -name '*.ccs.bam.part*' -delete 2>/dev/null
        freed_kb=$(( freed_kb + kb ))
        printf '  %-9s removed %s orphan file(s), %.1f GiB\n' "$s" "$npart" "$(awk -v k="$kb" 'BEGIN{print k/1048576}')"
    done < "$TMP/cnt"
    awk -v k="$freed_kb" 'BEGIN{printf "  reclaimed %.1f GiB total\n", k/1048576}'
}

if [ "$WATCH" = 1 ]; then
    while true; do
        clear; date
        if report; then echo; echo "ALL COMPLETE"; exit 0; fi
        sleep 30
    done
fi

date
report; RC=$?
[ "$CLEAN" = 1 ] && cleanup
[ "$RATE" = 1 ] && rate_check
# 0 only when every file is verified, so this is usable in a conditional
exit $RC
