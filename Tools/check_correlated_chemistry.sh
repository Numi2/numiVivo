#!/usr/bin/env bash
# Scoped portable check. Does not build or qualify the complete NumiVivo package.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ $# -gt 1 ]]; then echo "usage: $0 [output-directory]" >&2; exit 2; fi
OUT="${1:-$(mktemp -d "${TMPDIR:-/tmp}/numivivo-correlated.XXXXXX")}" 
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
SWIFTC="${SWIFTC:-swiftc}"
SOURCES=(
  Sources/NumiVivoKit/QM/VivoElectronicTypes.swift
  Sources/NumiVivoKit/QM/VivoQMDenseAlgebra.swift
  Sources/NumiVivoKit/QM/VivoGaussianIntegrals.swift
  Sources/NumiVivoKit/QM/VivoHartreeFock.swift
  Sources/NumiVivoKit/QM/VivoLDAFunctional.swift
  Sources/NumiVivoKit/QM/VivoRestrictedLDA.swift
  Sources/NumiVivoKit/QM/VivoGaussianElectrostaticPotential.swift
  Sources/NumiVivoKit/QMEnv/VivoCPCM.swift
  Sources/NumiVivoKit/QMEnv/VivoCPCMOperator.swift
  Sources/NumiVivoKit/QMEnv/VivoCPCMRHF.swift
  Sources/NumiVivoKit/ManyBody/VivoEmbeddedHamiltonian.swift
  Sources/NumiVivoKit/ManyBody/VivoConfigurationInteraction.swift
  Sources/NumiVivoKit/ManyBody/VivoCoupledCluster.swift
  Sources/NumiVivoKit/ManyBody/VivoCASSCF.swift
  Sources/NumiVivoKit/ManyBody/VivoManyBodySolver.swift
  Sources/NumiVivoKit/Embedding/VivoFragmentNumberMatching.swift
  Tools/NativeChemistry/CorrelatedChemistryChecks.swift
)
FILES=()
for source in "${SOURCES[@]}"; do FILES+=("$ROOT/$source"); done
"$SWIFTC" --version > "$OUT/compiler.txt"
uname -srm > "$OUT/platform.txt"
# Hash the actual compilation inputs, independent of git availability.
if command -v shasum >/dev/null 2>&1; then
  (cd "$ROOT" && shasum -a 256 "${SOURCES[@]}") > "$OUT/source-sha256.txt"
elif command -v sha256sum >/dev/null 2>&1; then
  (cd "$ROOT" && sha256sum "${SOURCES[@]}") > "$OUT/source-sha256.txt"
else
  echo "shasum or sha256sum is required to record source identities" >&2; exit 2
fi
FLAGS=(-swift-version 6 -O -parse-as-library)
if [[ "$(uname -s)" == Darwin ]]; then FLAGS+=(-framework Accelerate); fi
"$SWIFTC" "${FLAGS[@]}" "${FILES[@]}" -o "$OUT/correlated-checks"
"$OUT/correlated-checks" "$OUT" | tee "$OUT/results.log"
echo "Results and source identities: $OUT"
