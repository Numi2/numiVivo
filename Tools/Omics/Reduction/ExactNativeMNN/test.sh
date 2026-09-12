#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?new output directory}"
mkdir "$OUT"
xcrun clang++ -std=c++23 -O1 -g -fsanitize=address,undefined \
  -I "$ROOT/Sources/NumiVivoCore/include" \
  "$ROOT/Tools/Omics/Reduction/ExactNativeMNN/Test.cpp" \
  "$ROOT/Sources/NumiVivoCore/OmicsMNN.cpp" -o "$OUT/exact-mnn-tests"
"$OUT/exact-mnn-tests"
