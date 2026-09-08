#!/bin/bash
# Build pacbio_manifest.tsv from a sample list, validating every object against S3.
#
#   ./make_manifest.sh samples.tsv > pacbio_manifest.tsv
#
# Input is tab separated, one sample per line:
#   column 1  sample name          HG00438
#   column 2  cohort               HPRC or HPRC_PLUS  (the S3 subdirectory)
#   column 3  optional, space separated movie filenames to restrict to.
#             Leave it off to take every .bam in the sample's PacBio_HiFi directory.
#
# Output is tab separated: sample, filename, full s3 uri, exact size in bytes,
# sorted largest first so the longest downloads start earliest.
#
# Every path is confirmed to exist and its size comes from S3 itself, so the
# downloader can do a byte-exact check without a second network round trip.

# Ignore this line

set -uo pipefail

SAMPLES=${1:-}
[ -n "$SAMPLES" ] && [ -f "$SAMPLES" ] || { echo "usage: $0 samples.tsv > pacbio_manifest.tsv" >&2; exit 1; }
command -v aws >/dev/null 2>&1 || { echo "aws not on PATH, try: module load aws" >&2; exit 1; }

BUCKET=${BUCKET:-human-pangenomics}
PREFIX=${PREFIX:-working}
SUBDIR=${SUBDIR:-raw_data/PacBio_HiFi}

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

# List every sample's directory in parallel, one aws call each.
while IFS=$'\t' read -r sample cohort want; do
    [ -n "${sample:-}" ] || continue
    case "$sample" in \#*) continue ;; esac
    (
      aws s3 ls --no-sign-request "s3://$BUCKET/$PREFIX/$cohort/$sample/$SUBDIR/" 2>/dev/null \
        | awk -v s="$sample" -v c="$cohort" '$4 ~ /\.bam$/ {print s"\t"c"\t"$4"\t"$3}'
    ) > "$TMP/$sample.ls" &
done < "$SAMPLES"
wait

MISSING=0
: > "$TMP/rows"
while IFS=$'\t' read -r sample cohort want; do
    [ -n "${sample:-}" ] || continue
    case "$sample" in \#*) continue ;; esac
    if [ ! -s "$TMP/$sample.ls" ]; then
        echo "ERROR: no bams found for $sample under $cohort" >&2; MISSING=1; continue
    fi
    if [ -n "${want:-}" ]; then
        for f in $want; do
            line=$(awk -F'\t' -v f="$f" '$3==f' "$TMP/$sample.ls")
            if [ -z "$line" ]; then
                echo "ERROR: $sample requested $f but it is not in S3" >&2; MISSING=1; continue
            fi
            echo "$line" >> "$TMP/rows"
        done
    else
        cat "$TMP/$sample.ls" >> "$TMP/rows"
    fi
done < "$SAMPLES"

[ "$MISSING" = 0 ] || { echo "aborting, some requested files were not found" >&2; exit 1; }

sort -t"$(printf '\t')" -k4,4nr "$TMP/rows" \
  | awk -F'\t' -v b="$BUCKET" -v p="$PREFIX" -v d="$SUBDIR" \
      '{printf "%s\t%s\ts3://%s/%s/%s/%s/%s/%s\t%s\n", $1, $3, b, p, $2, $1, d, $3, $4}'

n=$(wc -l < "$TMP/rows"); tot=$(awk -F'\t' '{t+=$4} END{print t}' "$TMP/rows")
awk -v n="$n" -v t="$tot" 'BEGIN{printf "built manifest: %d files, %.2f TB\n", n, t/1e12}' >&2
