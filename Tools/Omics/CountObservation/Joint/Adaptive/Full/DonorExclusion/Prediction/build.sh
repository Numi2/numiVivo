#!/bin/bash
set -euo pipefail
BUILD="${1:?qualified scoped build directory}"
SOURCE="${2:?matching owner source snapshot}"
TARGET="${3:?new executable path}"
test ! -e "$TARGET"
swiftc -swift-version 6 -O -parse-as-library -I "$BUILD" \
 -I "$SOURCE/Sources/CNumiVivoZlib" -I "$SOURCE/Sources/NumiVivoCore/include" \
 -L "$BUILD" -lNumiVivoKit "$(dirname "$0")/Run.swift" \
 "$BUILD/OmicsHNSW.o" "$BUILD/OmicsGaussian.o" \
 -framework Accelerate -lc++ -o "$TARGET"
