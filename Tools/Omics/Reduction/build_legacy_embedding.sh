#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
LIB="${1:?scoped H5AD build directory}"
OUT="${2:?output executable}"
xcrun swiftc -swift-version 6 -O -parse-as-library -I "$LIB" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$LIB" -lNumiVivoKit \
  "$ROOT/Tools/Omics/Reduction/LegacyEmbeddingReference.swift" "$ROOT/Tools/Omics/Reduction/LegacyEmbeddingMain.swift" "$LIB/OmicsHNSW.o" -lc++ -o "$OUT"
