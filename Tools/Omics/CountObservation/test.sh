#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
OUT="${1:?built scoped runtime directory}"
DEVELOPER="$(xcode-select -p)"
FRAMEWORKS="$DEVELOPER/Platforms/MacOSX.platform/Developer/Library/Frameworks"
PLUGINS="$DEVELOPER/Toolchains/XcodeDefault.xctoolchain/usr/lib/swift/host/plugins/testing"
xcrun swiftc -swift-version 6 -parse-as-library -F "$FRAMEWORKS" \
  -Xfrontend -plugin-path -Xfrontend "$PLUGINS" -Xlinker -rpath -Xlinker "$FRAMEWORKS" \
  -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellNegativeBinomialTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellPerturbationTests.swift" \
  "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellCountObservationTests.swift" \
  "$ROOT/Tools/Omics/Reduction/GaussianKernel/TestMain.swift" \
  "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/count-observation-tests"
"$OUT/count-observation-tests"
xcrun swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
  "$ROOT/Tools/Omics/CountObservation/Check.swift" "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/count-observation-check"
