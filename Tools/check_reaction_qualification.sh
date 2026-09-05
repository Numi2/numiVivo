#!/usr/bin/env bash
# Native stationary/Hessian/RRHO work, independent references, and real CLI/store.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$(mktemp -d "${TMPDIR:-/tmp}/numivivo-reaction.XXXXXX")}" 
if [[ $# -gt 1 ]]; then echo 'usage: check_reaction_qualification.sh [output-directory]' >&2; exit 2; fi
if [[ "$(uname -s)" != Darwin ]]; then echo 'Run the production CLI/store qualification on macOS.' >&2; exit 2; fi
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
if [[ -n "${VIVO_REACTION_REFERENCE_DIR:-}" ]]; then
  REFERENCE="$(cd "$VIVO_REACTION_REFERENCE_DIR" && pwd)"
else
  REFERENCE="$OUT/oracle"
  if [[ -n "${VIVO_ORACLE_PYTHON:-}" ]]; then PYTHON="$VIVO_ORACLE_PYTHON"; else
    python3 -m venv "$OUT/oracle-venv";PYTHON="$OUT/oracle-venv/bin/python"
    "$PYTHON" -m pip install 'numpy==1.26.4' 'scipy==1.13.1' 'h5py==3.11.0' 'pyscf==2.8.0'
  fi
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    "$PYTHON" "$ROOT/ReferenceAdapters/PySCF/reaction_qualification_oracle.py" export "$REFERENCE"
  "$PYTHON" -m pip freeze > "$OUT/oracle-python-lock.txt"
fi
CLI="$OUT/scoped-cli"
bash "$ROOT/Tools/NativeChemistry/build_scoped_cli.sh" "$CLI" 2>&1 | tee "$OUT/cli-build.log"
shasum -a 256 "$ROOT/Tools/NativeChemistry/ReactionQualificationChecks.swift" \
 "$ROOT/Tools/NativeChemistry/check_reaction_cli.py" "$ROOT/ReferenceAdapters/PySCF/reaction_qualification_oracle.py" > "$OUT/test-sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -I "$CLI" -L "$CLI" -lNumiVivoKit -framework Accelerate \
 "$ROOT/Tools/NativeChemistry/ReactionQualificationChecks.swift" -o "$OUT/reaction-checks" 2>&1 | tee "$OUT/test-build.log"
"$OUT/reaction-checks" "$REFERENCE/reference.json" "$OUT/native" 2>&1 | tee "$OUT/native-checks.log"
python3 "$ROOT/ReferenceAdapters/PySCF/reaction_qualification_oracle.py" compare "$OUT/native" "$REFERENCE" "$OUT/pyscf-comparison.json"
python3 "$ROOT/Tools/NativeChemistry/check_reaction_cli.py" "$CLI/numivivo-chemistry" "$OUT/cli-checks" 2>&1 | tee "$OUT/cli-checks.log"
echo "Reaction qualification observations: $OUT"
