#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK="$(mktemp -d "${TMPDIR:-/tmp}/numivivo-kinetics.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT
SWIFTC="${SWIFTC:-swiftc}"
"$SWIFTC" --version
"$SWIFTC" -O -swift-version 6 -warnings-as-errors \
  "$ROOT/Sources/NumiVivoKit/Kinetics/VivoCovalentKinetics.swift" \
  "$ROOT/Sources/NumiVivoKit/Kinetics/VivoTargetEngagementProgramSource.swift" \
  "$ROOT/Sources/NumiVivoKit/Pharmacology/VivoTargetEngagementReference.swift" \
  "$ROOT/Sources/NumiVivoKit/Chemistry/VivoTransitionStateRate.swift" \
  "$ROOT/Sources/NumiVivoKit/Calibration/VivoTargetEngagementStudy.swift" \
  "$ROOT/Tools/Kinetics/ReferenceChecks.swift" -o "$WORK/reference-checks"
if [[ $# -eq 0 ]]; then
  "$WORK/reference-checks"
elif [[ $# -eq 2 && "$1" == "--generate-fixtures" ]]; then
  "$WORK/reference-checks" "$2"
else
  echo "usage: bash Tools/Kinetics/run_reference_checks.sh [--generate-fixtures <directory>]" >&2
  exit 64
fi
