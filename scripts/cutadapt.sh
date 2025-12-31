#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# PORTFOLIO-SAFE TEMPLATE
# - This script is anonymized and uses generic paths.
# - Set RAW_DIR to the directory containing your *_R1_001.fastq.gz / *_R2_001.fastq.gz pairs.
# - Set FWD_PRIMER / REV_PRIMER (or edit defaults) before running.
# ------------------------------------------------------------------------------

set -euo pipefail
RAW_DIR="${RAW_DIR:-raw_fastq}"
cd "$RAW_DIR"
FWD="${FWD_PRIMER:-CTTGGTCATTTAGAGGAAGTAA}"  # example default
REV="${REV_PRIMER:-GCTGCGTTCTTCATCGATGC}"   # example default

mkdir -p trimmed logs

summary="logs/cutadapt_summary.tsv"
echo -e "sample\ttotal_pairs\twith_fwd\twith_rev\ttoo_short\twritten" > "$summary"

shopt -s nullglob
for r1 in *_R1_001.fastq.gz; do
  r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  base="${r1%_R1_001.fastq.gz}"

  if [[ ! -f "$r2" ]]; then
    echo "WARN: missing R2 for $r1 -> expected $r2; skipping" >&2
    continue
  fi

  log="logs/${base}_cutadapt.log"

  # Trim primers when detected (not required to be at position 1)
  # Keep pairs even if only one read had a primer (common in ITS)
  cutadapt \
    -g "${FWD}" \
    -G "${REV}" \
    -e 0.15 \
    --minimum-length 50 \
    --pair-filter=any \
    -j 0 \
    -o "trimmed/${base}_R1_trim.fastq.gz" \
    -p "trimmed/${base}_R2_trim.fastq.gz" \
    "$r1" "$r2" | tee "$log" >/dev/null

  # Parse key numbers from the log (cutadapt output is stable for these lines)
  total_pairs=$(awk -F': +' '/Total read pairs processed:/ {gsub(/,/, "", $2); print $2}' "$log")
  with_fwd=$(awk -F': +' '/Read 1 with adapter:/ {gsub(/[(),%]/, "", $2); split($2,a," "); gsub(/,/, "", a[1]); print a[1]}' "$log")
  with_rev=$(awk -F': +' '/Read 2 with adapter:/ {gsub(/[(),%]/, "", $2); split($2,a," "); gsub(/,/, "", a[1]); print a[1]}' "$log")
  too_short=$(awk -F': +' '/Pairs that were too short:/ {gsub(/,/, "", $2); split($2,a," "); print a[1]}' "$log")
  written=$(awk -F': +' '/Pairs written \(passing filters\):/ {gsub(/,/, "", $2); split($2,a," "); print a[1]}' "$log")

  echo -e "${base}\t${total_pairs:-NA}\t${with_fwd:-NA}\t${with_rev:-NA}\t${too_short:-NA}\t${written:-NA}" >> "$summary"
done

echo "Done. Trimmed reads in: trimmed/"
echo "Logs + summary in: logs/ (see $summary)"