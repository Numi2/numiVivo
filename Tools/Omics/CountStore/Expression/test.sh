#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?path to built scoped runtime}"
DEVELOPER="$(xcode-select -p)"
FRAMEWORKS="$DEVELOPER/Platforms/MacOSX.platform/Developer/Library/Frameworks"
PLUGINS="$DEVELOPER/Toolchains/XcodeDefault.xctoolchain/usr/lib/swift/host/plugins/testing"
xcrun swiftc -swift-version 6 -parse-as-library -F "$FRAMEWORKS" \
  -Xfrontend -plugin-path -Xfrontend "$PLUGINS" \
  -Xlinker -rpath -Xlinker "$FRAMEWORKS" \
  -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellCountStreamTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellFileAxisTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellCohortTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellNBCohortTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellNBSupportTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellFileExpressionTests.swift" \
  "$ROOT/Tools/Omics/Reduction/GaussianKernel/TestMain.swift" \
  "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" "$OUT/OmicsMNN.o" -framework Accelerate -lc++ -o "$OUT/file-expression-tests"
"$OUT/file-expression-tests"
