#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [[ $# != 1 && $# != 3 ]]; then
  echo "usage: bash Tools/Posterior/run_portable_checks.sh <new-output-directory> [--rng-fragment <extracted-production-RNG.swift>]" >&2
  exit 64
fi
mkdir -p "$1"
OUT="$(cd "$1" && pwd)"
if [[ -e "$OUT/manifest.json" ]]; then echo "choose a new output directory; retained verification is immutable" >&2; exit 64; fi
WORK="$(mktemp -d "${TMPDIR:-/tmp}/numivivo-posterior.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
if [[ $# == 3 ]]; then
  [[ "$2" == --rng-fragment ]] || { echo "unknown option" >&2; exit 64; }
  RNG_SOURCE="$(cd "$(dirname "$3")" && pwd)/$(basename "$3")"
  cp "$RNG_SOURCE" "$WORK/Random.swift"
  RNG_MODE="provided verbatim production RNG excerpt; no claim to compiling its parent optimizer"
else
  RNG_SOURCE="$ROOT/Sources/NumiVivoKit/Calibration/VivoAdaptiveEnsembleOptimizer.swift"
  python3 - "$RNG_SOURCE" "$WORK/Random.swift" <<'PY'
import sys
from pathlib import Path
text=Path(sys.argv[1]).read_text()
marker='public struct VivoSplitMix64:'
if text.count(marker)!=1: raise SystemExit('ambiguous RNG extraction')
fragment=text[text.index(marker):]
# The public RNG is the final declaration in the current optimizer source.
# Reject a moved/extended extraction instead of silently including other code.
depth=0; end=None
for i,c in enumerate(fragment):
    if c=='{': depth+=1
    elif c=='}':
        depth-=1
        if depth==0: end=i+1; break
if end is None or fragment[end:].strip(): raise SystemExit('RNG boundary changed; update extraction explicitly')
Path(sys.argv[2]).write_text('import Foundation\n'+fragment)
PY
  RNG_MODE="verbatim public RNG extracted from production optimizer; other optimizer code excluded"
fi
cd "$ROOT"
SWIFTC="${SWIFTC:-swiftc}"
"$SWIFTC" --version | tee "$OUT/toolchain.log"
SOURCES=(
  Sources/NumiVivoKit/Kinetics/VivoCovalentKinetics.swift
  Sources/NumiVivoKit/Pharmacology/VivoTargetEngagementReference.swift
  Sources/NumiVivoKit/Pharmacology/VivoFiniteDrugReactions.swift
  Sources/NumiVivoKit/Pharmacology/VivoFiniteDrugWorkflow.swift
  Sources/NumiVivoKit/Calibration/VivoPosteriorModel.swift
  Sources/NumiVivoKit/Calibration/VivoPosteriorNumerics.swift
  Sources/NumiVivoKit/Calibration/VivoTemperedPosteriorSampler.swift
  Sources/NumiVivoKit/Calibration/VivoTargetEngagementStudy.swift
  Sources/NumiVivoKit/Calibration/VivoGaussianAssay.swift
  Sources/NumiVivoKit/Calibration/VivoTargetAssayModel.swift
  Sources/NumiVivoKit/Calibration/VivoTargetPosteriorProblem.swift
  Sources/NumiVivoKit/Calibration/VivoTargetPosteriorPrediction.swift
  Sources/NumiVivoKit/Calibration/VivoTargetPosteriorSensitivity.swift
  Sources/NumiVivoKit/Calibration/VivoTargetExperimentDesign.swift
  Sources/NumiVivoKit/Calibration/VivoAggregateOccupancyBenchmark.swift
)
case "$(uname -s)" in
 Linux) LINK=(-lcrypto); LIB=libNumiVivoPortable.so ;;
 Darwin) LINK=(); LIB=libNumiVivoPortable.dylib ;;
 *) echo "unsupported harness platform" >&2; exit 64 ;;
esac
"$SWIFTC" -O -whole-module-optimization -swift-version 6 -warnings-as-errors -enable-testing \
  -emit-library -emit-module -module-name NumiVivoPortable \
  Tools/Posterior/PortableSupport.swift "$WORK/Random.swift" "${SOURCES[@]}" "${LINK[@]}" \
  -emit-module-path "$WORK/NumiVivoPortable.swiftmodule" -o "$WORK/$LIB" 2>&1 | tee "$OUT/compile-library.log"
export NUMIVIVO_CHECK_OUTPUT="$OUT"
for CHECK in PortableChecks DomainChecks; do
  "$SWIFTC" -Onone -parse-as-library -swift-version 6 -warnings-as-errors -D NUMIVIVO_PORTABLE_SEPARATE \
    -I "$WORK" -L "$WORK" -lNumiVivoPortable -Xlinker -rpath -Xlinker "$WORK" \
    "Tools/Posterior/$CHECK.swift" -o "$WORK/$CHECK" 2>&1 | tee "$OUT/compile-$CHECK.log"
  "$WORK/$CHECK" "$OUT/$CHECK.json" | tee "$OUT/$CHECK.log"
done
python3 Tools/Posterior/write_manifest.py "$OUT" "$RNG_SOURCE" "$WORK/Random.swift" "$RNG_MODE" "${SOURCES[@]}"
echo "Retained selected-source results: $OUT"
