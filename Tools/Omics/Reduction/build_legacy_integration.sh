#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
LIB="${1:?scoped H5AD build directory}"
OUT="${2:?output executable}"
xcrun swiftc -swift-version 6 -O -parse-as-library -I "$LIB" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$LIB" -lNumiVivoKit \
  "$ROOT/Tools/Omics/Reduction/LegacyIntegrationReference.swift" "$ROOT/Tools/Omics/Reduction/LegacyIntegrationMain.swift" "$LIB/OmicsHNSW.o" "$LIB/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT"
