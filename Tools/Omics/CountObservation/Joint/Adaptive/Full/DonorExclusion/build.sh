#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../../../../.." && pwd)"
OUT="${1:?fresh scoped build directory}"
TARGET="${2:-$OUT/donor-excluded-calibration}"
test ! -e "$TARGET"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" \
 -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" \
 -L "$OUT" -lNumiVivoKit \
 "$ROOT/Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Calibrate.swift" \
 "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$TARGET"
