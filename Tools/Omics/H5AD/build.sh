#!/bin/bash
# Compile the actual Omics and artifact owners, without the unrelated Metal product.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
FILES=()
for name in VivoSparseCounts VivoBufferedCountRecords VivoOmicsFileSnapshot VivoH5ADCountStore VivoMultiAssay VivoH5ADCountAccess VivoMultiAssayH5MUImport VivoMultiAssayTenX VivoMultiAssayH5MU VivoMultiAssayIO VivoMatrixMarketCounts VivoOmicsLinearStatistics VivoOmicsNegativeBinomial VivoOmicsNBCohort VivoOmicsNBSupport VivoOmicsSourceDecoder VivoPseudobulkDifferentialExpression VivoSingleCellAnalysis VivoSingleCellArtifacts VivoSingleCellCampaign VivoSingleCellCampaignIO VivoSingleCellReduction VivoSingleCellPrograms VivoSingleCellReference VivoSingleCellReferenceIO VivoPerturbation VivoPerturbationIO VivoComposition VivoCompositionIO VivoTargetKernel VivoTargetKernelIO VivoH5ADReduction VivoH5ADPCA VivoFrozenPCAProjection VivoH5ADPCAQuery VivoWindowedPCANeighbors VivoHNSWNeighbors VivoPCAGraphStore VivoPCANeighborBundle VivoSingleCellNeighbors VivoClusteringGraph VivoPCAGraphClustering VivoSingleCellClustering VivoEmbeddingSchedule VivoPCAGraphEmbedding VivoSingleCellEmbedding VivoSingleCellIntegration VivoIntegrationMatrix VivoMNNIntegration VivoPCAMNNIntegration VivoPCAIntegration VivoSingleCellCohortAnalysis VivoSingleCellExamples VivoSingleCellExchange VivoSingleCellProcessing VivoHDF5 VivoSingleCellH5AD VivoH5ADPseudobulk VivoH5ADElements VivoH5ADAnnotations VivoH5ADProjectionIO VivoH5ADProjection VivoOmicsNBDevianceMoments VivoOmicsNBAdaptiveMoments VivoOmicsRobustLowess VivoOmicsNBQLGlobalScale VivoOmicsNBAbundance VivoOmicsQLSpecialFunctions VivoOmicsPrecisionLowess VivoOmicsQLModeration VivoOmicsNBQLInference VivoOmicsNBQLCohort; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
done
for name in VivoArtifactPrimitives CanonicalArtifact VivoArtifactStore VivoRootedFileStore; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Artifacts/$name.swift")
done
xcrun clang++ -std=c++23 -O3 -I "$ROOT/Sources/NumiVivoCore/include" -c "$ROOT/Sources/NumiVivoCore/OmicsHNSW.cpp" -o "$OUT/OmicsHNSW.o"
xcrun clang++ -std=c++23 -O3 -I "$ROOT/Sources/NumiVivoCore/include" -c "$ROOT/Sources/NumiVivoCore/OmicsGaussian.cpp" -o "$OUT/OmicsGaussian.o"
shasum -a 256 "$ROOT/Sources/NumiVivoCore/OmicsGaussian.cpp" "$ROOT/Sources/NumiVivoCore/include/NumiVivoCore/NumiVivoOmicsGaussian.h" "$ROOT/Sources/NumiVivoCore/OmicsHNSW.cpp" "$ROOT/Sources/NumiVivoCore/include/NumiVivoCore/NumiVivoOmicsHNSW.h" "$ROOT/Sources/NumiVivoCore/ThirdParty/hnswlib/"*.h "${FILES[@]}" "$ROOT/Tools/Omics/H5AD/Main.swift" > "$OUT/sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -enable-testing -module-name NumiVivoKit \
  -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -emit-module -emit-library -static "${FILES[@]}" \
  -emit-module-path "$OUT/NumiVivoKit.swiftmodule" -o "$OUT/libNumiVivoKit.a"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tools/Omics/H5AD/Main.swift" "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/h5ad-check"
