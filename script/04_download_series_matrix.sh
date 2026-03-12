#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
MANIFEST="$ROOT/data/metadata/dataset_manifest.tsv"
OUT_DIR="$ROOT/data/raw/series_matrix"
LOG_FILE="$ROOT/Log/series_matrix_download.log"

mkdir -p "$OUT_DIR"

download_matrix_files() {
  local gse="$1"
  local prefix="${gse%???}nnn"
  local base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}/${gse}/matrix/"
  local listing
  local files_text

  listing="$(curl -s "$base_url")"
  files_text="$(printf "%s" "$listing" | rg -o '[A-Za-z0-9._-]+_series_matrix\.txt\.gz' | sort -u || true)"

  if [[ -z "$files_text" ]]; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP  $gse no series matrix found" | tee -a "$LOG_FILE"
    return
  fi

  while IFS= read -r file; do
    [[ -z "$file" ]] && continue
    local out="$OUT_DIR/$file"
    local url="${base_url}${file}"

    if [[ -f "$out" ]]; then
      echo "[$(date '+%Y-%m-%d %H:%M:%S')] SKIP  $file already downloaded" | tee -a "$LOG_FILE"
      continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] START $file $url" | tee -a "$LOG_FILE"
    curl -fL "$url" -o "$out"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] DONE  $file -> $out" | tee -a "$LOG_FILE"
  done <<EOF
$files_text
EOF
}

tail -n +2 "$MANIFEST" | cut -f1 | while read -r dataset_id; do
  download_matrix_files "$dataset_id"
done

echo "All series matrix downloads completed."
