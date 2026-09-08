#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
swift --version
swiftc -swift-version 6 -parse-as-library \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsLinearStatistics.swift" \
  "$ROOT/Tools/Omics/statistics-smoke.swift" -o "$WORK/statistics-checks"
"$WORK/statistics-checks"
swiftc -swift-version 6 -parse-as-library -I "$ROOT/Sources/CNumiVivoZlib" \
  "$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsSourceDecoder.swift" \
  "$ROOT/Tools/Omics/gzip-smoke.swift" -o "$WORK/gzip-checks"
"$WORK/gzip-checks"
