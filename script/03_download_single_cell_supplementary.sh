#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
MANIFEST="${ROOT_DIR}/data/metadata/single_cell_supplementary_manifest.tsv"
OUT_DIR="${ROOT_DIR}/data/raw/single_cell"
LOG_FILE="${ROOT_DIR}/Log/single_cell_download_log.txt"

DATASET_FILTER="${1:-ALL}"
MODE="${2:-analysis_ready}"

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${LOG_FILE}")"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Manifest not found: ${MANIFEST}" >&2
  exit 1
fi

echo "# single-cell supplementary download log" > "${LOG_FILE}"
echo "# started: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"
echo "# dataset_filter=${DATASET_FILTER}" >> "${LOG_FILE}"
echo "# mode=${MODE}" >> "${LOG_FILE}"

tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r dataset_id sample_id sample_title sample_type source_name platform_id characteristic_summary relation_summary supplementary_file file_name file_type; do
  if [[ "${DATASET_FILTER}" != "ALL" && "${dataset_id}" != "${DATASET_FILTER}" ]]; then
    continue
  fi

  keep_file=0
  case "${MODE}" in
    analysis_ready)
      case "${dataset_id}" in
        GSE279086)
          if [[ "${file_name}" == *processed* ]]; then
            keep_file=1
          fi
          ;;
        GSE209781|GSE266146)
          keep_file=1
          ;;
      esac
      ;;
    processed_only)
      if [[ "${file_name}" == *processed* ]]; then
        keep_file=1
      fi
      ;;
    all)
      keep_file=1
      ;;
    *)
      echo "Unsupported mode: ${MODE}" >&2
      exit 1
      ;;
  esac

  if [[ "${keep_file}" -ne 1 ]]; then
    continue
  fi

  target_dir="${OUT_DIR}/${dataset_id}/${sample_id}"
  target_path="${target_dir}/${file_name}"
  tmp_path="${target_path}.part"
  mkdir -p "${target_dir}"

  if [[ -s "${target_path}" ]]; then
    echo "SKIP\t${dataset_id}\t${sample_id}\t${file_name}" | tee -a "${LOG_FILE}"
    continue
  fi

  echo "GET\t${dataset_id}\t${sample_id}\t${file_name}" | tee -a "${LOG_FILE}"
  rm -f "${tmp_path}"
  if curl -L --fail --retry 5 --retry-all-errors --retry-delay 2 --connect-timeout 30 -o "${tmp_path}" "${supplementary_file}" && [[ -s "${tmp_path}" ]]; then
    mv "${tmp_path}" "${target_path}"
  else
    rm -f "${tmp_path}"
    echo "ERR\t${dataset_id}\t${sample_id}\t${file_name}" | tee -a "${LOG_FILE}"
  fi
done

echo "# finished: $(date '+%Y-%m-%d %H:%M:%S')" >> "${LOG_FILE}"
echo "Download log written to ${LOG_FILE}"
