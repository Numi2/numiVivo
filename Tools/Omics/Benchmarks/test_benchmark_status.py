#!/usr/bin/env python3
"""Focused regression checks for benchmark_status.py."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE = Path(__file__).with_name("benchmark_status.py")
SPEC = importlib.util.spec_from_file_location("benchmark_status", MODULE)
assert SPEC is not None and SPEC.loader is not None
STATUS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STATUS)


def benchmark(identifier: str) -> dict[str, object]:
    return {
        "id": identifier,
        "status": "completed",
        "evidenceClass": "measured",
        "source": {"study": identifier},
        "scope": {"cells": 1},
        "result": f"declared {identifier} claim",
        "limitations": ["not a promotion claim"],
        "evidencePaths": [f"evidence/{identifier}.json"],
    }


class BenchmarkStatusTests(unittest.TestCase):
    def test_summary_is_sorted_and_preserves_claim_boundaries(self) -> None:
        manifest = {
            "schema": "numivivo.org/benchmark-manifest/v1",
            "updated": "2026-09-15",
            "status": "experimental-real-data-suite",
            "purpose": "test only",
            "biologicalOutcomePrediction": {"status": "not-established"},
            "benchmarks": [benchmark("zeta"), benchmark("alpha")],
        }
        raw = json.dumps(manifest).encode()
        result = STATUS.make_summary(
            manifest, hashlib.sha256(raw).hexdigest(), {"alpha"}
        )
        self.assertEqual(result["mode"], "declared-evidence-status-only")
        self.assertEqual([item["id"] for item in result["benchmarks"]], ["alpha"])
        self.assertEqual(result["benchmarks"][0]["claim"], "declared alpha claim")
        self.assertEqual(result["benchmarks"][0]["evidencePaths"], ["evidence/alpha.json"])

    def test_invalid_or_unknown_entries_fail_closed(self) -> None:
        manifest = {
            "schema": "numivivo.org/benchmark-manifest/v1",
            "benchmarks": [benchmark("only")],
        }
        with self.assertRaisesRegex(ValueError, "unknown benchmark"):
            STATUS.make_summary(manifest, "0" * 64, {"missing"})
        incomplete = benchmark("bad")
        incomplete.pop("limitations")
        with self.assertRaisesRegex(ValueError, "lacks limitations"):
            STATUS.make_summary(
                {"schema": "numivivo.org/benchmark-manifest/v1", "benchmarks": [incomplete]},
                "0" * 64,
            )

    def test_load_manifest_hashes_the_actual_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "manifest.json"
            raw = b'{"schema":"numivivo.org/benchmark-manifest/v1","benchmarks":[]}'
            path.write_bytes(raw)
            manifest, digest = STATUS.load_manifest(path)
        self.assertEqual(manifest["benchmarks"], [])
        self.assertEqual(digest, hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    unittest.main()
