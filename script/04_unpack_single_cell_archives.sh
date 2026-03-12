#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
RAW_DIR="${ROOT_DIR}/data/raw/single_cell"
OUT_DIR="${ROOT_DIR}/data/processed/single_cell_archives"
LOG_FILE="${ROOT_DIR}/Log/single_cell_unpack_log.txt"

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

echo "# single-cell archive unpack log" > "${LOG_FILE}"
echo "# started: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"

find "${RAW_DIR}" -type f -name "*.tar.gz" | sort | while read -r archive_path; do
  rel_path="${archive_path#${RAW_DIR}/}"
  dataset_id="$(echo "${rel_path}" | cut -d'/' -f1)"
  sample_id="$(echo "${rel_path}" | cut -d'/' -f2)"
  base_name="$(basename "${archive_path}" .tar.gz)"
  target_dir="${OUT_DIR}/${dataset_id}/${sample_id}/${base_name}"
  marker_file="${target_dir}/.unpacked_ok"

  mkdir -p "${target_dir}"

  if [[ -f "${marker_file}" ]]; then
    echo "SKIP\t${dataset_id}\t${sample_id}\t${base_name}" | tee -a "${LOG_FILE}"
    continue
  fi

  if ! tar -tzf "${archive_path}" >/dev/null 2>&1; then
    echo "WARN\t${dataset_id}\t${sample_id}\t${base_name}\tarchive not ready or corrupted" | tee -a "${LOG_FILE}"
    continue
  fi

  echo "UNPACK\t${dataset_id}\t${sample_id}\t${base_name}" | tee -a "${LOG_FILE}"
  tar -xzf "${archive_path}" -C "${target_dir}"
  touch "${marker_file}"
done

echo "# finished: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"
echo "Unpack log written to ${LOG_FILE}"
