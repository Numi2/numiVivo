#!/usr/bin/env python3
"""Fail before package resolution when the selected Apple compiler cannot build MLX.

This checks build prerequisites, not GPU execution or scientific correctness.
The MLX 0.31.6 pin requires Swift 6.3. No fallback compiler or package downgrade
is selected automatically. Both Swift launch paths used by the workflows are
checked so an unrelated PATH toolchain cannot silently override Xcode.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
from typing import Any

MINIMUM_SWIFT = (6, 3, 0)


def swift_version(text: str) -> tuple[int, int, int]:
    match = re.search(r"\bSwift version (\d+)\.(\d+)(?:\.(\d+))?\b", text)
    if not match:
        raise ValueError("Could not identify the Swift compiler version")
    return tuple(int(v or 0) for v in match.groups())


def validate(record: dict[str, Any]) -> list[str]:
    errors = []
    if record["system"] != "Darwin":
        errors.append("The complete NumiVivo product requires an Apple build host")
    if record["machine"] != "arm64":
        errors.append("This CI qualification requires an Apple Silicon host")
    if not record["developer_directory_exists"]:
        errors.append("The explicitly selected DEVELOPER_DIR does not exist")
    versions = []
    for label in ("swift", "xcrun_swift"):
        result = record["commands"][label]
        if result["returncode"] != 0:
            errors.append(label + " failed: " + result["output"].strip())
            continue
        try:
            version = swift_version(result["output"])
        except ValueError as exc:
            errors.append(label + ": " + str(exc))
            continue
        versions.append(version)
        if version < MINIMUM_SWIFT:
            errors.append(label + " is older than Swift 6.3 required by MLX 0.31.6")
    if len(versions) == 2 and versions[0] != versions[1]:
        errors.append("PATH Swift and Xcode Swift differ; repair toolchain selection")
    for label in ("xcode", "metal"):
        result = record["commands"][label]
        if result["returncode"] != 0 or not result["output"].strip():
            errors.append(label + " is unavailable: " + result["output"].strip())
    return errors


def command(arguments: list[str]) -> dict[str, Any]:
    try:
        result = subprocess.run(arguments, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, timeout=60,
                                check=False)
        return {"argv": arguments, "returncode": result.returncode, "output": result.stdout}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"argv": arguments, "returncode": -1, "output": str(exc)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    developer = os.environ.get("DEVELOPER_DIR", "")
    record: dict[str, Any] = {
        "system": platform.system(), "machine": platform.machine(),
        "developer_directory": developer,
        "developer_directory_exists": bool(developer) and Path(developer).is_dir(),
        "minimum_swift": list(MINIMUM_SWIFT),
        "commands": {
            "swift": command(["swift", "--version"]),
            "xcrun_swift": command(["xcrun", "swift", "--version"]),
            "xcode": command(["xcodebuild", "-version"]),
            "metal": command(["xcrun", "--sdk", "macosx", "--find", "metal"]),
        },
    }
    errors = validate(record)
    record["errors"] = errors
    record["status"] = "failed" if errors else "build-prerequisites-passed"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    for error in errors:
        print("Toolchain preflight: " + error, file=sys.stderr)
    if not errors:
        print("Apple build prerequisites passed; no GPU execution is implied")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
