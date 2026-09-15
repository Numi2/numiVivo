#!/usr/bin/env python3
"""Render a deterministic status view of the existing benchmark registry.

This program does not run, validate, or promote a benchmark.  It only exposes
the claims, scope, limitations, and exact retained evidence paths declared in
``benchmark-manifest.json`` so operational compatibility is not conflated with
biological or predictive qualification.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def load_manifest(path: Path) -> tuple[dict[str, Any], str]:
    raw = path.read_bytes()
    try:
        manifest = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError(f"invalid JSON in {path}: {error}") from error
    if not isinstance(manifest, dict):
        raise ValueError("benchmark manifest must be a JSON object")
    if manifest.get("schema") != "numivivo.org/benchmark-manifest/v1":
        raise ValueError("unsupported benchmark manifest schema")
    if not isinstance(manifest.get("benchmarks"), list):
        raise ValueError("benchmark manifest must contain a benchmarks array")
    return manifest, hashlib.sha256(raw).hexdigest()


def declared_entry(entry: dict[str, Any]) -> dict[str, Any]:
    required = ("id", "status", "evidenceClass", "source", "scope", "result",
                "limitations", "evidencePaths")
    missing = [field for field in required if field not in entry]
    if missing:
        raise ValueError(f"benchmark {entry.get('id', '<unknown>')} lacks {', '.join(missing)}")
    if not isinstance(entry["id"], str) or not entry["id"]:
        raise ValueError("benchmark id must be a nonempty string")
    if not isinstance(entry["evidencePaths"], list) or not all(
        isinstance(path, str) and path for path in entry["evidencePaths"]
    ):
        raise ValueError(f"benchmark {entry['id']} has invalid evidencePaths")
    if not isinstance(entry["limitations"], list) or not all(
        isinstance(limit, str) and limit for limit in entry["limitations"]
    ):
        raise ValueError(f"benchmark {entry['id']} has invalid limitations")
    return {
        "id": entry["id"],
        "status": entry["status"],
        "evidenceClass": entry["evidenceClass"],
        "source": entry["source"],
        "scope": entry["scope"],
        "claim": entry["result"],
        "evidencePaths": entry["evidencePaths"],
        "limitations": entry["limitations"],
    }


def make_summary(manifest: dict[str, Any], manifest_sha256: str,
                 selected_ids: set[str] | None = None) -> dict[str, Any]:
    entries = [declared_entry(entry) for entry in manifest["benchmarks"]]
    ids = [entry["id"] for entry in entries]
    if len(ids) != len(set(ids)):
        raise ValueError("benchmark manifest has duplicate ids")
    known = set(ids)
    unknown = (selected_ids or set()) - known
    if unknown:
        raise ValueError("unknown benchmark id(s): " + ", ".join(sorted(unknown)))
    if selected_ids:
        entries = [entry for entry in entries if entry["id"] in selected_ids]
    entries.sort(key=lambda entry: entry["id"])
    return {
        "mode": "declared-evidence-status-only",
        "manifest": {
            "schema": manifest["schema"],
            "updated": manifest.get("updated"),
            "status": manifest.get("status"),
            "sha256": manifest_sha256,
        },
        "purpose": manifest.get("purpose"),
        "biologicalOutcomePrediction": manifest.get("biologicalOutcomePrediction"),
        "benchmarks": entries,
    }


def markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Declared benchmark evidence status",
        "",
        "This view is derived from the existing benchmark manifest. It does not rerun or "
        "independently validate any benchmark.",
        "",
        f"- Manifest SHA-256: `{summary['manifest']['sha256']}`",
        f"- Manifest updated: `{summary['manifest']['updated']}`",
        f"- Suite status: `{summary['manifest']['status']}`",
        "",
        "| Benchmark | Status | Evidence class | Declared claim | Exact evidence paths |",
        "| --- | --- | --- | --- | --- |",
    ]
    for entry in summary["benchmarks"]:
        claim = str(entry["claim"]).replace("|", "\\|")
        paths = "<br>".join(entry["evidencePaths"])
        lines.append(
            f"| {entry['id']} | {entry['status']} | {entry['evidenceClass']} | "
            f"{claim} | {paths} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=Path(__file__).with_name("benchmark-manifest.json"),
        help="existing benchmark-manifest.json to render",
    )
    parser.add_argument("--benchmark", action="append", default=[],
                        help="render one declared benchmark id; repeatable")
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    arguments = parser.parse_args(argv)
    try:
        manifest, digest = load_manifest(arguments.manifest)
        summary = make_summary(manifest, digest, set(arguments.benchmark))
    except (OSError, ValueError) as error:
        print(f"benchmark-status: {error}", file=sys.stderr)
        return 2
    if arguments.format == "markdown":
        sys.stdout.write(markdown(summary))
    else:
        sys.stdout.write(canonical_json(summary) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
