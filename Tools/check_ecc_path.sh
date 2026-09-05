#!/usr/bin/env bash
# The real scoped Apple module/CLI, numerical path checks, pinned external
# comparisons and cache/CLI tests. No fake stores or generated interfaces.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$(mktemp -d "${TMPDIR:-/tmp}/numivivo-ecc-path.XXXXXX")}" 
if [[ $# -gt 1 ]]; then echo 'usage: check_ecc_path.sh [output-directory]' >&2; exit 2; fi
if [[ "$(uname -s)" != Darwin ]]; then echo 'Run this production CLI/store check on macOS.' >&2; exit 2; fi
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
if [[ -n "${VIVO_PATH_REFERENCE_DIR:-}" ]]; then
  REFERENCE="$(cd "$VIVO_PATH_REFERENCE_DIR" && pwd)"
else
  REFERENCE="$OUT/oracle"
  if [[ -n "${VIVO_ORACLE_PYTHON:-}" ]]; then PYTHON="$VIVO_ORACLE_PYTHON"; else
    python3 -m venv "$OUT/oracle-venv"
    PYTHON="$OUT/oracle-venv/bin/python"
    "$PYTHON" -m pip install 'numpy==1.26.4' 'scipy==1.13.1' 'h5py==3.11.0' 'pyscf==2.8.0'
  fi
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$PYTHON" "$ROOT/ReferenceAdapters/PySCF/compare_ecc_path.py" export "$REFERENCE"
  "$PYTHON" -m pip freeze > "$OUT/oracle-python-lock.txt"
fi
if [[ -n "${VIVO_SCOPED_CLI_DIR:-}" ]]; then
  CLI="$(cd "$VIVO_SCOPED_CLI_DIR" && pwd)"
else
  CLI="$OUT/scoped-cli"
  bash "$ROOT/Tools/NativeChemistry/build_scoped_cli.sh" "$CLI" 2>&1 | tee "$OUT/cli-build.log"
fi
shasum -a 256 "$ROOT/Tools/NativeChemistry/ECCPathChecks.swift" \
  "$ROOT/ReferenceAdapters/PySCF/compare_ecc_path.py" \
  "$ROOT/Tools/NativeChemistry/check_ecc_path_cli.py" > "$OUT/check-source-sha256.txt"
swiftc -swift-version 6 -O -parse-as-library -I "$CLI" -L "$CLI" -lNumiVivoKit -framework Accelerate \
  "$ROOT/Tools/NativeChemistry/ECCPathChecks.swift" -o "$OUT/ecc-path-checks" 2>&1 | tee "$OUT/test-build.log"
"$OUT/ecc-path-checks" "$OUT/native" 2>&1 | tee "$OUT/native-checks.log"
python3 "$ROOT/ReferenceAdapters/PySCF/compare_ecc_path.py" compare "$OUT/native" "$REFERENCE" "$OUT/pyscf-comparison.json"
python3 "$ROOT/Tools/NativeChemistry/check_ecc_path_cli.py" "$CLI/numivivo-chemistry" "$REFERENCE" "$OUT/cli-checks" \
  2>&1 | tee "$OUT/cli-checks.log"
echo "Shared ECC path observations: $OUT"
