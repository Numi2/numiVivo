#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
FILES=()
for name in VivoSparseCounts VivoSingleCellAnalysis VivoOmicsSourceDecoder VivoOmicsLinearStatistics VivoOmicsNegativeBinomial VivoOmicsNBDevianceMoments VivoOmicsNBAdaptiveMoments VivoOmicsNBSupport VivoOmicsNBCohort VivoPseudobulkDifferentialExpression; do
 FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
done
MAIN="$ROOT/Tools/Omics/NegativeBinomial/QuasiLikelihood/Moments/Adaptive/Main.swift"
shasum -a 256 "${FILES[@]}" "$MAIN" > "$OUT/sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -I "$ROOT/Sources/CNumiVivoZlib" "${FILES[@]}" "$MAIN" -o "$OUT/nb-adaptive-moments"
