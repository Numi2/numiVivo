#!/bin/bash
# Compile the exact FP64 QR/NB owners for bounded profile diagnostics.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?output directory}"
mkdir -p "$OUT"
swiftc --version > "$OUT/toolchain.txt"
git -C "$ROOT" rev-parse HEAD > "$OUT/owner-commit.txt"
FILES=("$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsLinearStatistics.swift"
       "$ROOT/Sources/NumiVivoKit/Omics/VivoOmicsNegativeBinomial.swift"
       "$ROOT/Tools/Omics/NegativeBinomial/ProfileAudit/Main.swift")
shasum -a 256 "${FILES[@]}" > "$OUT/sources.sha256"
swiftc -swift-version 6 -O -parse-as-library "${FILES[@]}" -o "$OUT/profile-audit"
shasum -a 256 "$OUT/profile-audit" > "$OUT/binary.sha256"
