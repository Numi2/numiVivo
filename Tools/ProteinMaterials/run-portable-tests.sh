#!/bin/sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT HUP INT TERM
mkdir -p "$TMP/Sources/NumiVivoKit" "$TMP/Tests/NumiVivoIntegrationTests"
cat > "$TMP/Package.swift" <<'PACKAGE'
// swift-tools-version: 6.0
import PackageDescription
let package = Package(name: "ProteinStressPortable", targets: [
    .target(name: "NumiVivoKit"),
    .testTarget(name: "NumiVivoIntegrationTests", dependencies: ["NumiVivoKit"])
])
PACKAGE
# Compile the actual numerical source and exact existing vector declaration.
# This is a portable observable test, NOT the full Apple package or GPU runtime.
{
    printf 'import Foundation\n'
    sed -n '/^public struct VivoVector3D:/,/^}/p' "$ROOT/Sources/NumiVivoKit/Structure/VivoMolecularStructure.swift"
} > "$TMP/Sources/NumiVivoKit/VivoVector3D.swift"
cp "$ROOT/Sources/NumiVivoKit/ProteinMaterials/VivoProteinStressObservables.swift" "$TMP/Sources/NumiVivoKit/"
cp "$ROOT/Tests/NumiVivoIntegrationTests/ProteinStressObservableTests.swift" "$TMP/Tests/NumiVivoIntegrationTests/"
swift test --package-path "$TMP" -c release
