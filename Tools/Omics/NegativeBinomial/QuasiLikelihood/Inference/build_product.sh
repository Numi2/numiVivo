#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../../.." && pwd)"
OUT="${1:?output directory}"
LIBRARY="${2:?qualified H5AD library build directory}"
mkdir -p "$OUT"
shasum -a 256 "$LIBRARY/libNumiVivoKit.a" "$LIBRARY/OmicsHNSW.o" "$ROOT/Tools/Omics/NegativeBinomial/QuasiLikelihood/Inference/Product.swift" > "$OUT/product-sources.sha256"
swiftc -swift-version 6 -O -parse-as-library -I "$LIBRARY" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$LIBRARY" -lNumiVivoKit \
  "$ROOT/Tools/Omics/NegativeBinomial/QuasiLikelihood/Inference/Product.swift" "$LIBRARY/OmicsHNSW.o" -lc++ -o "$OUT/ql-product"
