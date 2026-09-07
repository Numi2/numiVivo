#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
out="${1:-/tmp/numivivo-adaptive-chemistry}"
mkdir -p "$out"
cd "$root"
# Retain the actual production SplitMix64 implementation without importing the
# optimizer's unrelated model/UI types. Its complete source is hashed below.
python3 - "$out" <<'PY'
import pathlib,sys
source=pathlib.Path('Sources/NumiVivoKit/Calibration/VivoAdaptiveEnsembleOptimizer.swift').read_text()
start=source.index('struct VivoSplitMix64:')
pathlib.Path(sys.argv[1],'ProductionSplitMix64.swift').write_text('import Foundation\n'+source[start:])
PY
sources=(Sources/NumiVivoKit/QM/VivoElectronicTypes.swift Sources/NumiVivoKit/QM/VivoQMDenseAlgebra.swift
 Sources/NumiVivoKit/Chemistry/VivoMultistateMBAR.swift Sources/NumiVivoKit/AdaptiveChemistry/*.swift)
{ for source in "${sources[@]}" Sources/NumiVivoKit/Calibration/VivoAdaptiveEnsembleOptimizer.swift Tools/NativeChemistry/AdaptiveChemistryChecks.swift; do shasum -a 256 "$source"; done; } > "$out/source-sha256.txt"
swiftc --version > "$out/compiler.txt"
uname -a > "$out/platform.txt"
swiftc -swift-version 6 -Onone -whole-module-optimization -parse-as-library \
 "${sources[@]}" "$out/ProductionSplitMix64.swift" Tools/NativeChemistry/AdaptiveChemistryChecks.swift \
 -o "$out/adaptive-checks" 2>&1 | tee "$out/build.log"
"$out/adaptive-checks" "$out" 2>&1 | tee "$out/run.log"
