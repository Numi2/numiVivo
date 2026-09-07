#!/usr/bin/env bash
# Focused native execution with exact source/compiler provenance. Not a full
# Apple package, GPU test, chemical benchmark or measured speedup claim.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${1:?usage: check_property_refinement.sh output-directory}"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
case "$OUT/" in "$ROOT/"|"$ROOT/Sources/"*|"$ROOT/Tools/"*) echo 'Use an output directory outside source.' >&2; exit 2;; esac
FILES=()
for directory in QM ManyBody Embedding Refinement; do
  for file in "$ROOT/Sources/NumiVivoKit/$directory/"*.swift; do FILES+=("$file"); done
done
for file in "$ROOT/Sources/NumiVivoKit/QMEnv/"VivoCPCM*.swift; do FILES+=("$file"); done
FILES+=("$ROOT/Sources/NumiVivoKit/QMEnv/VivoSmoothCPCM.swift"
        "$ROOT/Sources/NumiVivoKit/Structure/VivoMolecularStructure.swift"
        "$ROOT/Sources/NumiVivoKit/Geometry/VivoCartesianGeometry.swift"
        "$ROOT/Tools/NativeChemistry/PropertyRefinementChecks.swift")
FLAGS=(-swift-version 6 -Onone -whole-module-optimization -parse-as-library)
if [[ "$(uname -s)" == Darwin ]]; then FLAGS+=(-framework Accelerate); fi
"${SWIFTC:-swiftc}" --version > "$OUT/compiler.txt"; uname -srm > "$OUT/platform.txt"
if command -v shasum >/dev/null; then shasum -a 256 "${FILES[@]}" > "$OUT/source-sha256.txt"
else sha256sum "${FILES[@]}" > "$OUT/source-sha256.txt"; fi
"${SWIFTC:-swiftc}" "${FLAGS[@]}" "${FILES[@]}" -o "$OUT/property-refinement-checks" > "$OUT/build.log" 2>&1
"$OUT/property-refinement-checks" "$OUT/results" | tee "$OUT/run.log"
