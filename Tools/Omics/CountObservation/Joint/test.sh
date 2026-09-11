#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../../../.." && pwd)"
OUT="${1:?built scoped runtime}"
D="$(xcode-select -p)"
F="$D/Platforms/MacOSX.platform/Developer/Library/Frameworks"
P="$D/Toolchains/XcodeDefault.xctoolchain/usr/lib/swift/host/plugins/testing"
swiftc -swift-version 6 -parse-as-library -F "$F" -Xfrontend -plugin-path -Xfrontend "$P" -Xlinker -rpath -Xlinker "$F" \
 -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellNegativeBinomialTests.swift" \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellPerturbationTests.swift" \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellCountObservationTests.swift" \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellCountCalibrationTests.swift" \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellPairedCountTests.swift" \
 "$ROOT/Tests/NumiVivoIntegrationTests/SingleCellJointCountTests.swift" \
 "$ROOT/Tools/Omics/Reduction/GaussianKernel/TestMain.swift" \
 "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/joint-tests"
"$OUT/joint-tests"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -I "$ROOT/Sources/CNumiVivoZlib" -I "$ROOT/Sources/NumiVivoCore/include" -L "$OUT" -lNumiVivoKit \
 "$ROOT/Tools/Omics/CountObservation/Joint/Run.swift" "$OUT/OmicsHNSW.o" "$OUT/OmicsGaussian.o" -framework Accelerate -lc++ -o "$OUT/joint-counts"
