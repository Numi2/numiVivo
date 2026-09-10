#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
FILES=()
for name in VivoSparseCounts VivoSingleCellAnalysis VivoOmicsSourceDecoder VivoOmicsLinearStatistics VivoOmicsNegativeBinomial VivoOmicsNBSupport VivoOmicsNBCohort VivoPseudobulkDifferentialExpression; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
done
shasum -a 256 "${FILES[@]}" "$ROOT/Tools/Omics/NegativeBinomial/HeldOutRisk/Main.swift" > "$OUT/sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -I "$ROOT/Sources/CNumiVivoZlib" \
  "${FILES[@]}" "$ROOT/Tools/Omics/NegativeBinomial/HeldOutRisk/Main.swift" -o "$OUT/nb-heldout"
