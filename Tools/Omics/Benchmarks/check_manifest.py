#!/usr/bin/env python3
"""Validate the checked-in public benchmark registry without downloading data.

The registry is a scope and provenance contract, not a substitute for running a
benchmark. This checker verifies that declared source identities, replicate
semantics, references, limitations and executable evidence paths are present.
Large source files and archived receipts are intentionally not required locally.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse


SCHEMA = "numivivo.org/benchmark-manifest/v1"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_ENTRY_KEYS = {
    "id", "status", "evidenceClass", "modality", "organism", "source",
    "scope", "design", "expectedBiology", "references", "evidencePaths", "result", "limitations",
}


def fail(message: str) -> None:
    raise ValueError(message)


def require(condition: bool, message: str) -> None:
    if not condition:
        fail(message)


def url(value: object, label: str) -> None:
    require(isinstance(value, str), f"{label} must be a string")
    parsed = urlparse(value)
    require(parsed.scheme in {"http", "https"} and bool(parsed.netloc), f"{label} must be an HTTP(S) URL")


def sha(value: object, label: str) -> None:
    require(isinstance(value, str) and SHA256.fullmatch(value) is not None, f"{label} must be a lowercase SHA-256")


def positive_integer(value: object, label: str, allow_zero: bool = False) -> None:
    require(isinstance(value, int) and not isinstance(value, bool), f"{label} must be an integer")
    require(value >= (0 if allow_zero else 1), f"{label} must be {'nonnegative' if allow_zero else 'positive'}")


def git_path_exists(root: Path, relative: str) -> bool:
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        return False
    if (root / path).exists():
        return True
    result = subprocess.run(
        ["git", "cat-file", "-e", f"HEAD:{relative}"],
        cwd=root,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.returncode == 0


def validate(manifest: dict, root: Path) -> dict:
    require(manifest.get("schema") == SCHEMA, "unsupported benchmark manifest schema")
    require(re.fullmatch(r"20[0-9]{2}-[0-9]{2}-[0-9]{2}", str(manifest.get("updated", ""))) is not None,
            "updated must be an ISO date")
    outcome = manifest.get("biologicalOutcomePrediction")
    require(isinstance(outcome, dict), "biologicalOutcomePrediction must be an object")
    require(outcome.get("status") == "not-established", "outcome status must remain not-established")
    for key in ("availableInputs", "missingEvidence"):
        values = outcome.get(key)
        require(isinstance(values, list) and values and all(isinstance(v, str) and v.strip() for v in values),
                f"outcome {key} must contain non-empty strings")

    entries = manifest.get("benchmarks")
    require(isinstance(entries, list) and entries, "benchmarks must be a non-empty list")
    ids: set[str] = set()
    classes = {"measured", "technical-measurement"}
    for index, entry in enumerate(entries):
        prefix = f"benchmarks[{index}]"
        require(isinstance(entry, dict), f"{prefix} must be an object")
        require(set(entry) == REQUIRED_ENTRY_KEYS, f"{prefix} has unexpected or missing keys")
        identifier = entry["id"]
        require(isinstance(identifier, str) and re.fullmatch(r"[a-z0-9][a-z0-9-]+", identifier) is not None,
                f"{prefix}.id is invalid")
        require(identifier not in ids, f"duplicate benchmark id: {identifier}")
        ids.add(identifier)
        require(entry["status"] in {"completed", "ineligible", "partial"}, f"{prefix}.status is invalid")
        require(entry["evidenceClass"] in classes, f"{prefix}.evidenceClass is invalid")
        for key in ("modality", "organism", "result"):
            require(isinstance(entry[key], str) and entry[key].strip(), f"{prefix}.{key} must be non-empty")

        source = entry["source"]
        require(isinstance(source, dict), f"{prefix}.source must be an object")
        require(isinstance(source.get("study"), str) and source["study"].strip(), f"{prefix}.source.study is required")
        for key in ("publication", "download"):
            if key in source:
                url(source[key], f"{prefix}.source.{key}")
        if "sha256" in source:
            sha(source["sha256"], f"{prefix}.source.sha256")
        if "artifactSha256" in source:
            require(isinstance(source["artifactSha256"], list) and source["artifactSha256"],
                    f"{prefix}.source.artifactSha256 must be non-empty")
            for artifact in source["artifactSha256"]:
                require(isinstance(artifact, dict) and isinstance(artifact.get("kind"), str),
                        f"{prefix}.source.artifactSha256 entry is invalid")
                sha(artifact.get("sha256"), f"{prefix}.source.artifactSha256.sha256")

        scope = entry["scope"]
        require(isinstance(scope, dict), f"{prefix}.scope must be an object")
        for key in ("cells", "features", "nonzeros", "donors", "libraries", "sites"):
            if key in scope:
                positive_integer(scope[key], f"{prefix}.scope.{key}", allow_zero=(key == "donors"))
        require(isinstance(scope.get("replicateSemantics"), str) and scope["replicateSemantics"].strip(),
                f"{prefix}.scope.replicateSemantics is required")

        design = entry["design"]
        require(isinstance(design, dict), f"{prefix}.design must be an object")
        for key in ("comparison", "formula", "holdout"):
            require(isinstance(design.get(key), str) and design[key].strip(), f"{prefix}.design.{key} is required")
        for key in ("expectedBiology", "references", "evidencePaths", "limitations"):
            values = entry[key]
            require(isinstance(values, list) and values and all(isinstance(v, str) and v.strip() for v in values),
                    f"{prefix}.{key} must contain non-empty strings")
        for evidence_path in entry["evidencePaths"]:
            require(git_path_exists(root, evidence_path), f"{prefix} evidence path is not in HEAD: {evidence_path}")
        require(any(word in " ".join(entry["limitations"]).lower() for word in ("no ", "not ", "remain", "unqualified", "insufficient")),
                f"{prefix}.limitations must retain at least one boundary statement")

    return {"schema": SCHEMA, "benchmarks": len(entries), "ids": sorted(ids),
            "outcomePrediction": outcome["status"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=Path(__file__).with_name("benchmark-manifest.json"))
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    try:
        manifest = json.loads(args.manifest.read_text())
        result = validate(manifest, args.repo_root.resolve())
    except (OSError, json.JSONDecodeError, ValueError) as error:
        print(f"benchmark manifest: FAIL: {error}", file=sys.stderr)
        return 1
    print(json.dumps({"status": "passed", **result}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
