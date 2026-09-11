#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?qualified runtime directory}"
FRAMEWORKS="$(xcode-select -p)/Platforms/MacOSX.platform/Developer/Library/Frameworks"
PLUGINS="$(xcode-select -p)/Toolchains/XcodeDefault.xctoolchain/usr/lib/swift/host/plugins/testing"
swiftc -plugin-path "$PLUGINS" -F "$FRAMEWORKS" -Xlinker -rpath -Xlinker "$FRAMEWORKS" -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tests/NumiVivoIntegrationTests/MultiAssayTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/MultiAssayH5MUImportTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/MultiAssayVisiumTests.swift" \
  "$ROOT/Tools/Omics/Reduction/GaussianKernel/TestMain.swift" \
  "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/visium-tests"
: "${NUMIVIVO_HDF5_LIBRARY:?qualified HDF5 library}"
NUMIVIVO_TEST_HDF5=1 "$OUT/visium-tests"
