#!/bin/bash
# Compile the actual Omics and artifact owners, without the unrelated Metal product.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
FILES=()
for name in VivoSparseCounts VivoMatrixMarketCounts VivoOmicsLinearStatistics VivoOmicsNegativeBinomial VivoOmicsNBCohort VivoOmicsSourceDecoder VivoPseudobulkDifferentialExpression VivoSingleCellAnalysis VivoSingleCellArtifacts VivoSingleCellCampaign VivoSingleCellCampaignIO VivoSingleCellReduction VivoSingleCellNeighbors VivoSingleCellClustering VivoSingleCellEmbedding VivoSingleCellIntegration VivoSingleCellCohortAnalysis VivoSingleCellExamples VivoSingleCellExchange VivoSingleCellProcessing VivoHDF5 VivoSingleCellH5AD VivoH5ADPseudobulk VivoH5ADElements VivoH5ADAnnotations; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
done
for name in VivoArtifactPrimitives CanonicalArtifact VivoArtifactStore VivoRootedFileStore; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Artifacts/$name.swift")
done
shasum -a 256 "${FILES[@]}" "$ROOT/Tools/Omics/H5AD/Main.swift" > "$OUT/sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -enable-testing -module-name NumiVivoKit \
  -I "$ROOT/Sources/CNumiVivoZlib" -emit-module -emit-library -static "${FILES[@]}" \
  -emit-module-path "$OUT/NumiVivoKit.swiftmodule" -o "$OUT/libNumiVivoKit.a"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tools/Omics/H5AD/Main.swift" -o "$OUT/h5ad-check"
