#!/bin/bash
# Compile the actual Omics and artifact owners, including optional count Metal kernels.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
FILES=()
for name in VivoAccessibilityTFIDF VivoAccessibilityTFIDFIO VivoSparseCounts VivoBufferedCountRecords VivoOmicsFileSnapshot VivoH5ADCountStore VivoMetalCountNormalization VivoMultiAssay VivoH5ADCountAccess VivoMultiAssayH5MUImport VivoMultiAssayTenX VivoMultiAssayVisium VivoMultiAssayH5MU VivoMultiAssayIO VivoMatrixMarketCounts VivoOmicsLinearStatistics VivoOmicsNegativeBinomial VivoCountObservation VivoCountObservationCalibration VivoPairedCountMoments VivoCountRateLikelihood VivoJointCountResponse VivoAdaptiveJointCountResponse VivoOmicsNBCohort VivoOmicsNBSupport VivoOmicsSourceDecoder VivoOmicsDesignObservation VivoExpressionReportJSON VivoFileExpression VivoPseudobulkDifferentialExpression VivoSingleCellAnalysis VivoSingleCellArtifacts VivoSingleCellCampaign VivoSingleCellCampaignIO VivoSingleCellReduction VivoSingleCellPrograms VivoSingleCellReference VivoReferenceLogistic VivoCellTypistReference VivoCellTypistReferenceIO VivoSingleCellReferenceIO VivoPerturbation VivoPerturbationIO VivoDurationPerturbation VivoDurationPerturbationIO VivoComposition VivoCompositionIO VivoTargetKernel VivoTargetKernelIO VivoH5ADReduction VivoH5ADPCA VivoFrozenPCAProjection VivoH5ADPCAQuery VivoWindowedPCANeighbors VivoHNSWNeighbors VivoPCAGraphStore VivoPCANeighborBundle VivoSingleCellNeighbors VivoClusteringGraph VivoPCAGraphClustering VivoSingleCellClustering VivoEmbeddingSchedule VivoPCAGraphEmbedding VivoSingleCellEmbedding VivoSingleCellIntegration VivoIntegrationMatrix VivoMNNIntegration VivoPCAMNNIntegration VivoPCAIntegration VivoSingleCellCohortAnalysis VivoSingleCellExamples VivoSingleCellExchange VivoSingleCellProcessing VivoHDF5 VivoSingleCellH5AD VivoH5ADPseudobulk VivoCountStreamPseudobulk VivoFileCellAxis VivoFileCountStream VivoH5ADElements VivoH5ADAnnotations VivoH5ADProjectionIO VivoH5ADLegacyIO VivoH5ADProjection VivoOmicsNBDevianceMoments VivoOmicsNBAdaptiveMoments VivoOmicsRobustLowess VivoOmicsNBQLGlobalScale VivoOmicsNBAbundance VivoOmicsQLSpecialFunctions VivoOmicsPrecisionLowess VivoOmicsQLModeration VivoOmicsNBQLInference VivoOmicsNBQLCohort; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
done
for name in VivoArtifactPrimitives CanonicalArtifact VivoArtifactStore VivoRootedFileStore; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Artifacts/$name.swift")
done
CLI=()
if [[ "${2:-}" == "--with-cli" ]]; then
  # Compile the actual product router and identity boundary, with their owning
  # document/error definitions. No generated interfaces or command stand-ins.
  for name in VivoPerturbationAggregateBatch VivoH5ADPrograms; do
    FILES+=("$ROOT/Sources/NumiVivoKit/Omics/$name.swift")
  done
  FILES+=("$ROOT/Sources/NumiVivoKit/QM/VivoElectronicTypes.swift"
    "$ROOT/Sources/NumiVivoKit/Kinetics/VivoCovalentKinetics.swift"
    "$ROOT/Sources/NumiVivoKit/Kinetics/VivoKineticsDocumentIO.swift")
  CLI=("$ROOT/Sources/NumiVivoCLI/VivoSingleCellCLICommands.swift"
    "$ROOT/Sources/NumiVivoCLI/VivoWorkflowCLIImplementation.swift"
    "$ROOT/Tools/Omics/H5AD/CLIMain.swift")
fi
xcrun clang++ -std=c++23 -O3 -I "$ROOT/Sources/NumiVivoCore/include" -c "$ROOT/Sources/NumiVivoCore/OmicsHNSW.cpp" -o "$OUT/OmicsHNSW.o"
xcrun clang++ -std=c++23 -O3 -I "$ROOT/Sources/NumiVivoCore/include" -c "$ROOT/Sources/NumiVivoCore/OmicsGaussian.cpp" -o "$OUT/OmicsGaussian.o"
shasum -a 256 "$ROOT/Sources/NumiVivoCore/OmicsGaussian.cpp" "$ROOT/Sources/NumiVivoCore/include/NumiVivoCore/NumiVivoOmicsGaussian.h" "$ROOT/Sources/NumiVivoCore/OmicsHNSW.cpp" "$ROOT/Sources/NumiVivoCore/include/NumiVivoCore/NumiVivoOmicsHNSW.h" "$ROOT/Sources/NumiVivoCore/ThirdParty/hnswlib/"*.h "${FILES[@]}" "$ROOT/Tools/Omics/H5AD/Main.swift" > "$OUT/sources.sha256"
if [[ "${2:-}" == "--with-cli" ]]; then
  shasum -a 256 "${CLI[@]}" >> "$OUT/sources.sha256"
fi
swiftc -swift-version 6 -O -parse-as-library -enable-testing -module-name NumiVivoKit \
  -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -emit-module -emit-library -static "${FILES[@]}" \
  -emit-module-path "$OUT/NumiVivoKit.swiftmodule" -o "$OUT/libNumiVivoKit.a"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tools/Omics/H5AD/Main.swift" "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/h5ad-check"

if [[ "${2:-}" == "--with-cli" ]]; then
  swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
    "${CLI[@]}" "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/numivivo-omics"
fi
