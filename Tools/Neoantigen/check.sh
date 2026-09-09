#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
BUILD=$(mktemp -d "${TMPDIR:-/tmp}/numivivo-neoantigen.XXXXXX")
trap 'rm -rf "$BUILD"' EXIT
LINK=()
if [[ "$(uname -s)" == "Linux" ]]; then LINK=(-lcrypto); fi
swiftc -swift-version 6 -O \
  Tools/Posterior/PortableSupport.swift \
  Sources/NumiVivoKit/Neoantigen/VivoNeoantigenWorkbench.swift \
  Sources/NumiVivoKit/Neoantigen/VivoNeoantigenExample.swift \
  Sources/NumiVivoKit/Neoantigen/VivoNeoantigenHTML.swift \
  Tools/Neoantigen/Check.swift "${LINK[@]}" -o "$BUILD/neoantigen-check"
"$BUILD/neoantigen-check" "$@"
