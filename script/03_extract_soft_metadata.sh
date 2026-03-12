#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
ROOT="${ROOT_DIR}"
LOG_FILE="$ROOT/Log/soft_metadata_extract.log"

python3 "$ROOT/script/python/01_extract_soft_metadata.py" | tee "$LOG_FILE"
