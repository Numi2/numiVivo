#!/usr/bin/env bash
# Build selected production sources and the real electronic/path command routers.
# This is not the complete NumiVivo product. No stand-ins or generated interfaces.
# Static linkage keeps the executing binary hash bound to the numerical code.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${1:?usage: build_scoped_cli.sh output-directory}"
if [[ "$(uname -s)" != Darwin ]]; then
  echo 'The scoped integration check uses Apple CryptoKit and Accelerate; run on macOS.' >&2; exit 2
fi
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
FILES=()
for directory in QM ManyBody Embedding; do
  for file in "$ROOT/Sources/NumiVivoKit/$directory/"*.swift; do FILES+=("$file"); done
 done
for file in "$ROOT/Sources/NumiVivoKit/QMEnv/"VivoCPCM*.swift; do FILES+=("$file"); done
FILES+=("$ROOT/Sources/NumiVivoKit/QMEnv/VivoSmoothCPCM.swift")
for name in VivoArtifactPrimitives CanonicalArtifact VivoArtifactStore VivoRootedFileStore; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Artifacts/$name.swift")
done
FILES+=("$ROOT/Sources/NumiVivoKit/JSONValue.swift")
for name in VivoChemistryWorkflow VivoElectronicWorkflow VivoElectronicWorkflowOperations VivoAdvancedElectronicWorkflow VivoAdvancedChemistryOperations VivoElectronicRequestDispatch VivoECCDMETOperation VivoMolecularECCPathWorkflow; do
  FILES+=("$ROOT/Sources/NumiVivoKit/Workflow/$name.swift")
done
CLI=("$ROOT/Sources/NumiVivoCLI/VivoElectronicCLICommands.swift" "$ROOT/Sources/NumiVivoCLI/VivoECCPathCLICommands.swift" "$ROOT/Tools/NativeChemistry/ScopedCLIMain.swift")
shasum -a 256 "${FILES[@]}" "${CLI[@]}" > "$OUT/sources.sha256"
swiftc --version > "$OUT/compiler.txt"; uname -srm > "$OUT/platform.txt"
swiftc -swift-version 6 -O -whole-module-optimization -parse-as-library \
  -module-name NumiVivoKit -emit-module -emit-library -static -framework Accelerate \
  "${FILES[@]}" -emit-module-path "$OUT/NumiVivoKit.swiftmodule" -o "$OUT/libNumiVivoKit.a"
swiftc -swift-version 6 -O -parse-as-library -I "$OUT" -L "$OUT" -lNumiVivoKit \
  "${CLI[@]}" -o "$OUT/numivivo-chemistry"
"$OUT/numivivo-chemistry" chemistry-help
"$OUT/numivivo-chemistry" chemistry-path-help
