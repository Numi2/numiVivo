#!/usr/bin/env python3
"""Build the exact native binder owners in a retained, portable Swift package.

The OpenSSL/canonical-JSON facade is harness-only; artifact primitives are native. This does not qualify the full
Apple package, production artifact store, GPU backend or biological prediction.
"""
import argparse
import json
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def prepare(root: Path) -> None:
    if root.exists():
        raise FileExistsError(root)
    library = root / "Sources/NumiVivoKit"
    cli = root / "Sources/BinderCLI"
    tests = root / "Tests/BinderTests"
    for folder in (library, cli, tests):
        folder.mkdir(parents=True)
    sources = [ROOT / "Tools/BinderBenchmark/PortableJSON.swift",
               ROOT / "Sources/NumiVivoKit/Artifacts/VivoArtifactPrimitives.swift"]
    sources += [ROOT / "Sources/NumiVivoKit/Structure" / name for name in (
        "VivoMolecularStructure.swift", "VivoStructureValidator.swift", "VivoMolecularInterface.swift")]
    sources += sorted((ROOT / "Sources/NumiVivoKit/Binder").glob("*.swift"))
    for source in sources:
        shutil.copy2(source, library / source.name)
    shutil.copy2(ROOT / "Sources/NumiVivoCLI/VivoBinderCLICommands.swift", cli)
    (cli / "main.swift").write_text(
        'import Foundation\nexit(VivoBinderCLICommands().run(arguments: Array(CommandLine.arguments.dropFirst())))\n'
    )
    for source in sorted((ROOT / "Tests/NumiVivoIntegrationTests").glob("VivoBinder*Tests.swift")):
        shutil.copy2(source, tests / source.name)
    (root / "Package.swift").write_text('''// swift-tools-version: 6.0
import PackageDescription
let package = Package(name: "BinderCheck", targets: [
    .target(name: "NumiVivoKit", linkerSettings: [.linkedLibrary("crypto", .when(platforms: [.linux]))]),
    .executableTarget(name: "BinderCLI", dependencies: ["NumiVivoKit"]),
    .testTarget(name: "BinderTests", dependencies: ["NumiVivoKit"])
], swiftLanguageModes: [.v6])
''')


def build(root: Path, test: bool = False) -> Path:
    prepare(root)
    subprocess.run(["swift", "build", "--package-path", str(root), "--jobs", "2"], check=True)
    if test:
        subprocess.run(["swift", "test", "--package-path", str(root), "--jobs", "2"], check=True)
    bin_path = subprocess.check_output(
        ["swift", "build", "--package-path", str(root), "--show-bin-path"], text=True
    ).strip()
    binary = Path(bin_path) / "BinderCLI"
    if not binary.is_file():
        raise FileNotFoundError(binary)
    return binary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path, help="new retained package directory")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    binary = build(args.package.absolute(), args.test)
    print(json.dumps({"binary": str(binary), "scope": "portable exact-source binder subset"}))


if __name__ == "__main__":
    main()
