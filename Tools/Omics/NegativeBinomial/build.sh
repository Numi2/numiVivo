#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:?output directory}"
bash "$ROOT/Tools/Omics/H5AD/build.sh" "$OUT"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" \
  -L "$OUT" -lNumiVivoKit "$ROOT/Tools/Omics/NegativeBinomial/Main.swift" -o "$OUT/nb-check"
shasum -a 256 "$ROOT/Tools/Omics/NegativeBinomial/Main.swift" >> "$OUT/sources.sha256"
