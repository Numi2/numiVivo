#!/usr/bin/env python3
"""Run the exact native binder owner/tests without Apple-only package dependencies."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]

def main() -> None:
    with tempfile.TemporaryDirectory(prefix="numivivo-binder-") as tmp:
        root = Path(tmp)
        (root / "Sources/NumiVivoKit").mkdir(parents=True)
        (root / "Tests/BinderTests").mkdir(parents=True)
        for source in sorted((ROOT / "Sources/NumiVivoKit/Binder").glob("*.swift")):
            shutil.copy2(source, root / "Sources/NumiVivoKit" / source.name)
        for test in sorted((ROOT / "Tests/NumiVivoIntegrationTests").glob("VivoBinder*Tests.swift")):
            shutil.copy2(test, root / "Tests/BinderTests" / test.name)
        (root / "Package.swift").write_text('''// swift-tools-version: 6.0
import PackageDescription
let package = Package(name: "BinderCheck", targets: [
    .target(name: "NumiVivoKit"),
    .testTarget(name: "BinderTests", dependencies: ["NumiVivoKit"])
], swiftLanguageModes: [.v6])
''')
        subprocess.run(["swift", "test", "--package-path", str(root), "--jobs", "2"], check=True)

if __name__ == "__main__":
    main()
