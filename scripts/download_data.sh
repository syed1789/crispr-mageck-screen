#!/bin/bash
# Fetches the HL60/KBM7 pooled CRISPR knockout screen count table (Wang et al.,
# Science 2014, "Genetic Screens in Human Cells Using the CRISPR-Cas9 System";
# read counts re-hosted as the canonical MAGeCK RRA tutorial dataset).
# Source: https://github.com/davidliwei/mageck2-demo (BSD-3-Clause demo data)
set -euo pipefail
cd "$(dirname "$0")/../data"

URL="https://raw.githubusercontent.com/davidliwei/mageck2-demo/main/demo1_rra_from_count_table/sample.txt"
OUT="sgrna_counts_hl60_kbm7.txt"

if [ -f "$OUT" ]; then
  echo "Already downloaded: $OUT"
  exit 0
fi

echo "Downloading sgRNA count table from $URL"
curl -L --fail -o "$OUT" "$URL"

lines=$(wc -l < "$OUT")
echo "Downloaded $OUT ($lines lines)"
