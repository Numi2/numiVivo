#!/usr/bin/env bash
# Reproducible executable conformance, not an oracle-only export.
# Output contains exact source/fixture identities and failing check reports.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:-$(mktemp -d "${TMPDIR:-/tmp}/numivivo-conformance.XXXXXX")}" 
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
if [[ "$OUT" == "$ROOT" || "$OUT" == "$ROOT/Sources"* || "$OUT" == "$ROOT/Tools"* ]]; then
  echo 'Use an output directory outside source and tool directories.' >&2; exit 2
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ -n "${VIVO_ORACLE_PYTHON:-}" ]] && ! command -v "$VIVO_ORACLE_PYTHON" >/dev/null 2>&1; then
  echo "Configured oracle interpreter is not executable: $VIVO_ORACLE_PYTHON" >&2; exit 2
fi
# An existing pinned oracle environment may be selected explicitly. Otherwise
# create an isolated conformance-only environment, never a production dependency.
if [[ -n "${VIVO_REFERENCE_DIR:-}" ]]; then
  REFERENCE="$(cd "$VIVO_REFERENCE_DIR" && pwd)"
else
  if [[ -z "${VIVO_ORACLE_PYTHON:-}" ]]; then
    python3 -m venv "$OUT/oracle-environment"
    VIVO_ORACLE_PYTHON="$OUT/oracle-environment/bin/python"
    "$VIVO_ORACLE_PYTHON" -m pip install numpy==1.26.4 scipy==1.13.1 h5py==3.11.0 pyscf==2.8.0 > "$OUT/oracle-install.log" 2>&1
  fi
  "$VIVO_ORACLE_PYTHON" -m pip freeze > "$OUT/oracle-lock.txt"
  REFERENCE="$OUT/reference"
  "$VIVO_ORACLE_PYTHON" "$ROOT/ReferenceAdapters/PySCF/export_native_reference.py" "$REFERENCE" > "$OUT/oracle.log" 2>&1
fi
# Never silently test fewer systems/methods because a reference disappeared.
for file in manifest.json h2.json lih.json water.json; do test -f "$REFERENCE/$file"; done
python3 "$ROOT/ReferenceAdapters/PySCF/verify_reference_manifest.py" "$REFERENCE"
mkdir -p "$OUT/source" "$OUT/native"
FILES=()
for directory in QM ManyBody Embedding; do
  for file in "$ROOT/Sources/NumiVivoKit/$directory/"*.swift; do FILES+=("$file"); done
 done
for file in "$ROOT/Sources/NumiVivoKit/QMEnv/"VivoCPCM*.swift; do FILES+=("$file"); done
FILES+=("$ROOT/Sources/NumiVivoKit/QMEnv/VivoSmoothCPCM.swift")
SNAPSHOT=()
for file in "${FILES[@]}"; do
  relative="${file#"$ROOT/"}"; mkdir -p "$OUT/source/$(dirname "$relative")"
  cp "$file" "$OUT/source/$relative"; SNAPSHOT+=("$OUT/source/$relative")
done
for test in ExternalConformance ECCIntegrationChecks SmoothSolventChecks CorrelatedChemistryChecks; do
  if [[ "$test" == CorrelatedChemistryChecks ]]; then
    # The original standalone regression remains unchanged. Add only the import
    # needed to compile that driver against this reusable production module.
    printf '@testable import NumiVivoNumerics\n' > "$OUT/source/$test.swift"
    cat "$ROOT/Tools/NativeChemistry/$test.swift" >> "$OUT/source/$test.swift"
  else
    cp "$ROOT/Tools/NativeChemistry/$test.swift" "$OUT/source/$test.swift"
  fi
done
if command -v shasum >/dev/null 2>&1; then
  (cd "$OUT/source" && find . -name '*.swift' -type f -exec shasum -a 256 {} \; | sort) > "$OUT/source-sha256.txt"
else
  (cd "$OUT/source" && find . -name '*.swift' -type f -exec sha256sum {} \; | sort) > "$OUT/source-sha256.txt"
fi
cp "$REFERENCE/manifest.json" "$OUT/reference-manifest.json"
SWIFTC="${SWIFTC:-swiftc}"
"$SWIFTC" --version > "$OUT/compiler.txt"; uname -srm > "$OUT/platform.txt"
FLAGS=(-swift-version 6 -O -whole-module-optimization -parse-as-library)
EXT=so
if [[ "$(uname -s)" == Darwin ]]; then FLAGS+=(-framework Accelerate); EXT=dylib; fi
"$SWIFTC" "${FLAGS[@]}" -enable-testing -emit-module -emit-library -module-name NumiVivoNumerics \
  "${SNAPSHOT[@]}" -emit-module-path "$OUT/native/NumiVivoNumerics.swiftmodule" \
  -o "$OUT/native/libNumiVivoNumerics.$EXT" > "$OUT/build-numerics.log" 2>&1
for test in ExternalConformance ECCIntegrationChecks SmoothSolventChecks CorrelatedChemistryChecks; do
  "$SWIFTC" -swift-version 6 -O -parse-as-library -D NUMIVIVO_SCOPED_NUMERICS \
    -I "$OUT/native" -L "$OUT/native" -lNumiVivoNumerics -Xlinker -rpath -Xlinker "$OUT/native" \
    "$OUT/source/$test.swift" -o "$OUT/native/$test" > "$OUT/build-$test.log" 2>&1
done
"$OUT/native/ExternalConformance" "$REFERENCE" "$OUT/external-results.json" | tee "$OUT/external.log"
"$OUT/native/ECCIntegrationChecks" "$OUT/ecc" | tee "$OUT/ecc.log"
"$OUT/native/SmoothSolventChecks" "$OUT/solvent" | tee "$OUT/solvent.log"
"$OUT/native/CorrelatedChemistryChecks" "$OUT/correlated" | tee "$OUT/correlated.log"
# The independent ECC reconstruction uses only source Hamiltonians and solver
# declarations; native energies/matrices are observations to compare, not inputs
# to an oracle fit. It also checks the exported bath-projector reconstruction.
if [[ -n "${VIVO_ORACLE_PYTHON:-}" ]]; then
  "$VIVO_ORACLE_PYTHON" "$ROOT/ReferenceAdapters/PySCF/compare_ecc_result.py" "$OUT/ecc" | tee "$OUT/ecc-oracle.log"
fi
echo "Conformance artifacts: $OUT"
