#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
swift --version
swiftc -swift-version 6 -parse-as-library -I "$ROOT/Sources/CNumiVivoZlib" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSparseCounts.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoMatrixMarketCounts.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSingleCellAnalysis.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSingleCellCampaign.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSingleCellCampaignIO.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsSourceDecoder.swift" \
  "$ROOT/Sources/NumiVivoKit/Artifacts/VivoRootedFileStore.swift" \
  "$ROOT/Tools/Omics/main.swift" -o "$WORK/omics-checks"
"$WORK/omics-checks"
swiftc -swift-version 6 -parse-as-library \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsLinearStatistics.swift" \
  "$ROOT/Tools/Omics/statistics-smoke.swift" -o "$WORK/statistics-checks"
"$WORK/statistics-checks"
