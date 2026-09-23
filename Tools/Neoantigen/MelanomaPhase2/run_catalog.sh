#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 OUTPUT_DIR [INPUT_ROOT]" >&2
  echo "INPUT_ROOT defaults to the current directory and must contain inputs/ from SOURCE_LOCK.json." >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
OUTPUT_DIR="$1"
INPUT_ROOT="${2:-$PWD}"
exec "${PYTHON_BIN:-python3}" -B "$SCRIPT_DIR/build_phase2_catalog.py" \
  --input-root "$INPUT_ROOT" --output-dir "$OUTPUT_DIR"
