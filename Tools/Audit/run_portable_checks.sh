#!/usr/bin/env bash
# Narrow repository-source checks. No package build, Metal, network or mocks.
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cd "$root"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/numivivo-audit.XXXXXX")
trap 'rm -rf "$scratch"' EXIT
cxx=${CXX:-clang++}
swiftc=${SWIFTC:-swiftc}
command -v "$cxx" >/dev/null
command -v "$swiftc" >/dev/null
"$cxx" --version | head -n 1
"$swiftc" --version
flags=(-std=c++23 -O1 -g -Wall -Wextra -fno-omit-frame-pointer)
case "${NVIVO_AUDIT_SANITIZERS:-1}" in
  1) flags+=(-fsanitize=address,undefined) ;;
  0) ;;
  *) echo 'NVIVO_AUDIT_SANITIZERS must be 0 or 1' >&2; exit 64 ;;
esac
core=Sources/NumiVivoCore
"$cxx" "${flags[@]}" -I "$core/include" -I "$core" \
  Tools/Audit/PackMutationChecks.cpp \
  "$core/Pack.cpp" "$core/PackValidation.cpp" "$core/SourceSemantics.cpp" \
  "$core/Units.cpp" "$core/Diagnostics.cpp" "$core/Json.cpp" "$core/Hash.cpp" \
  -o "$scratch/native-checks"
"$scratch/native-checks"
"$cxx" -std=c++23 -fsyntax-only -Wall -Wextra -Wconversion -I "$core/include" \
  "$core/Pack.cpp" "$core/PackValidation.cpp" "$core/SourceSemantics.cpp" "$core/Safety.cpp"
"$swiftc" -swift-version 6 Sources/NumiVivoKit/Artifacts/VivoRootedFileStore.swift \
  Tools/Audit/RootedFileStoreChecks.swift -o "$scratch/rooted-checks"
"$scratch/rooted-checks"
"$swiftc" -swift-version 6 Sources/NumiVivoKit/Coupling/VivoCouplingIterationController.swift \
  Tools/Audit/CouplingIterationChecks.swift -o "$scratch/coupling-checks"
"$scratch/coupling-checks"
# Syntax parsing is deliberately reported separately from typechecking.
for source in \
  Sources/NumiVivoKit/Artifacts/VivoProgramPack.swift \
  Sources/NumiVivoKit/Artifacts/VivoArtifactStore.swift \
  Sources/NumiVivoKit/Physiology/VivoPhysiologyModelValidator.swift \
  Sources/NumiVivoKit/Physiology/VivoPhysiologyRuntime.swift \
  Sources/NumiVivoKit/Physiology/VivoMolecularPhysiologyBridge.swift \
  Sources/NumiVivoKit/Physiology/VivoMolecularPhysiologyCoordinator.swift \
  Sources/NumiVivoKit/Physiology/VivoPhysiologyResumeCheckpoint.swift; do
  "$swiftc" -frontend -parse "$source"
done
printf '%s\n' 'Selected Swift syntax checks passed: 7 (not Apple SDK typechecking)' \
  'Portable audit checks complete; Apple/GPU validation remains separate.'
