#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
MANIFEST="$ROOT/data/metadata/dataset_manifest.tsv"
OUT_DIR="$ROOT/data/raw/geo_soft"
LOG_FILE="$ROOT/Log/download_log.txt"

mkdir -p "$OUT_DIR"

download_one() {
  local gse="$1"
  local prefix="${gse%???}nnn"
  local url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/soft/${gse}_family.soft.gz"
  local out="$OUT_DIR/${gse}_family.soft.gz"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $gse $url" | tee -a "$LOG_FILE"
  curl -fL "$url" -o "$out"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE  $gse -> $out" | tee -a "$LOG_FILE"
}

tail -n +2 "$MANIFEST" | cut -f1 | while read -r dataset_id; do
  if [[ -f "$OUT_DIR/${dataset_id}_family.soft.gz" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP  $dataset_id already downloaded" | tee -a "$LOG_FILE"
    continue
  fi
  download_one "$dataset_id"
done

echo "All GEO SOFT downloads completed."
