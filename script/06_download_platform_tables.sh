#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
OUT_DIR="${ROOT_DIR}/data/raw/platform_tables"
LOG_FILE="${ROOT_DIR}/Log/06_download_platform_tables.log"

mkdir -p "${OUT_DIR}"

platforms=(
  "GPL571"
  "GPL17586"
  "GPL24120"
  "GPL22945"
  "GPL23126"
)

{
  echo "[platform_table_download] $(date '+%Y-%m-%d %H:%M:%S')"

  for gpl in "${platforms[@]}"; do
    url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gpl}&targ=self&form=text&view=data"
    out_file="${OUT_DIR}/${gpl}_platform_table.txt"

    echo "Downloading ${gpl} -> ${out_file}"
    curl -L --fail --silent --show-error "${url}" -o "${out_file}"
    wc -l "${out_file}"
  done
} | tee "${LOG_FILE}"
