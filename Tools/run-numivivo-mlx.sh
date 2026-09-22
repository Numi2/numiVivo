#!/usr/bin/env bash
# Build and run the MLX-backed NumiVivo CLI with the Metal resource bundle
# produced by Xcode. SwiftPM does not produce that bundle for mlx-swift 0.31.6.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
derived_data="${NUMIVIVO_XCODE_DERIVED_DATA:-$repo_root/.build/xcode-mlx}"
configuration="${NUMIVIVO_XCODE_CONFIGURATION:-Release}"

if (( $# == 0 )); then
    printf 'usage: %s <numivivo arguments...>\n' "$0" >&2
    exit 64
fi

if [[ "$(uname -s)" != Darwin || "$(uname -m)" != arm64 ]]; then
    printf 'MLX cell-response commands require an Apple-silicon macOS host.\n' >&2
    exit 69
fi

if ! command -v xcodebuild >/dev/null 2>&1; then
    printf 'xcodebuild is required to build the MLX Metal runtime bundle.\n' >&2
    exit 69
fi

if [[ "$configuration" != "Release" && "$configuration" != "Debug" ]]; then
    printf 'NUMIVIVO_XCODE_CONFIGURATION must be Release or Debug.\n' >&2
    exit 64
fi

xcodebuild build -quiet \
    -scheme NumiVivo-Package \
    -configuration "$configuration" \
    -destination 'platform=macOS,arch=arm64' \
    -derivedDataPath "$derived_data" \
    -skipPackagePluginValidation \
    CODE_SIGNING_ALLOWED=NO

runtime_dir="$derived_data/Build/Products/$configuration"
cli="$runtime_dir/numivivo"
mlx_metallib="$runtime_dir/mlx-swift_Cmlx.bundle/Contents/Resources/default.metallib"

if [[ ! -x "$cli" || ! -f "$mlx_metallib" ]]; then
    printf 'Xcode did not produce the MLX runtime bundle.\n' >&2
    exit 70
fi

exec "$cli" "$@"
