#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
swift --version
swiftc -swift-version 6 -parse-as-library \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSparseCounts.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoMatrixMarketCounts.swift" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoSingleCellAnalysis.swift" \
  "$ROOT/Tools/Omics/main.swift" -o "$WORK/omics-checks"
"$WORK/omics-checks"
