#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
ARCHIVE_DIR="${ROOT_DIR}/data/processed/single_cell_archives/GSE266146"
OUT_DIR="${ROOT_DIR}/data/processed/single_cell_matrices/GSE266146"
LOG_FILE="${ROOT_DIR}/Log/gse266146_expand_log.txt"

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

echo "# GSE266146 expansion log" > "${LOG_FILE}"
echo "# started: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"

find "${ARCHIVE_DIR}" -type f -name "*.zip" | sort | while read -r zip_path; do
  sample_id="$(echo "${zip_path}" | awk -F'/' '{print $(NF-2)}')"
  base_name="$(basename "${zip_path}" .zip)"
  target_dir="${OUT_DIR}/${sample_id}/${base_name}"
  marker_file="${target_dir}/.expanded_ok"

  mkdir -p "${target_dir}"

  if [[ -f "${marker_file}" ]]; then
    echo "SKIP\t${sample_id}\t${base_name}" | tee -a "${LOG_FILE}"
    continue
  fi

  echo "EXPAND\t${sample_id}\t${base_name}" | tee -a "${LOG_FILE}"
  unzip -qo "${zip_path}" -d "${target_dir}"

  inner_zip="$(find "${target_dir}" -maxdepth 1 -type f -name "*_MEX.zip" | head -n 1 || true)"
  if [[ -n "${inner_zip}" ]]; then
    mkdir -p "${target_dir}/mex"
    unzip -qo "${inner_zip}" -d "${target_dir}/mex"
  fi

  touch "${marker_file}"
done

echo "# finished: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"
echo "Expansion log written to ${LOG_FILE}"
