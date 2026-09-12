#!/usr/bin/env python3
"""Portable contract tests. Fabricated scalar tables are NOT biological evidence."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
import evaluate as e


class AdapterTests(unittest.TestCase):
    def scores(self, value):
        return {**{m: value for m in e.METRICS}, "avg_score": value}

    def table(self, directory, value=-0.2):
        path = Path(directory) / "scores.csv"
        fields = ["metric", "from_baseline", "from_replicate", "anchor_source", "real_bundle_id", "real_bundle_digest", "anchor_digest"]
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
            for metric in (*e.METRICS, "avg_score"):
                w.writerow(dict(metric=metric, from_baseline="" if metric == "avg_score" else "0.1",
                                from_replicate=value, anchor_source="real_bundle", real_bundle_id="synthetic-test-only",
                                real_bundle_digest="a" * 64, anchor_digest="b" * 64))
        return path

    def test_empty_run_cannot_score(self):
        self.assertIsNone(e.aggregate(["A", "B", "C"], {})["score"])

    def test_partial_run_cannot_score(self):
        r = e.aggregate(["A", "B", "C"], {"A": self.scores(0.9), "B": self.scores(0.8)})
        self.assertIsNone(r["score"]); self.assertEqual(r["missingContexts"], ["C"])

    def test_equal_context_mean_and_no_global_clipping(self):
        r = e.aggregate(["A", "B", "C"], {"A": self.scores(-2), "B": self.scores(2), "C": self.scores(3)})
        self.assertEqual(r["score"], 1); self.assertEqual(r["contextScores"]["A"]["avg_score"], -2)

    def test_foreign_context_and_duplicate_context_rejected(self):
        with self.assertRaises(ValueError): e.aggregate(["A", "B", "C"], {"D": self.scores(1)})
        with self.assertRaises(ValueError): e.aggregate(["A", "A", "C"], {})

    def test_missing_metric_rejected(self):
        s = self.scores(1); del s[e.METRICS[0]]
        with self.assertRaises(ValueError): e.aggregate(["A", "B", "C"], {"A": s})

    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError): e.aggregate(["A", "B", "C"], {"A": self.scores(float("nan"))})

    def test_wrong_metric_mean_rejected(self):
        s = self.scores(0.5); s["avg_score"] = 1
        with self.assertRaises(ValueError): e.aggregate(["A", "B", "C"], {"A": s})

    def test_official_output_column_and_negative_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertAlmostEqual(e.read_scores(self.table(tmp))["avg_score"], -0.2)

    def test_diagnostic_average_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.table(tmp); p.write_text(p.read_text().replace("avg_score,,", "avg_score,0.8,"))
            with self.assertRaises(ValueError): e.read_scores(p)

    def test_missing_anchor_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.table(tmp); p.write_text(p.read_text().replace("real_bundle,", "diagnostic,"))
            with self.assertRaises(ValueError): e.read_scores(p)

    def test_duplicate_or_extra_metric_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self.table(tmp); p.write_text(p.read_text() + p.read_text().splitlines()[1] + "\n")
            with self.assertRaises(ValueError): e.read_scores(p)

    def test_duplicate_json_fields_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "x.json"; p.write_text('{"phase":"validation","phase":"final-test"}')
            with self.assertRaises(ValueError): e.read_json(p)

    def test_nonfinite_json_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "x.json"; p.write_text('{"score":NaN}')
            with self.assertRaises(ValueError): e.read_json(p)

    def test_snapshot_symlink_rejected_and_bytes_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); (p / "a").write_bytes(b"original")
            h = e.digest(p / "a", p / "b")
            self.assertEqual(h, e.digest(p / "b")); self.assertEqual((p / "b").read_bytes(), b"original")
            (p / "link").symlink_to(p / "a")
            with self.assertRaises(OSError): e.digest(p / "link")
            with self.assertRaises(FileExistsError): e.digest(p / "a", p / "b")

    def test_pins_match_native_contract_and_lock(self):
        here = Path(__file__).resolve().parent
        lock = e.read_json(here / "evaluator-lock.json")
        self.assertEqual(lock["commit"], e.COMMIT); self.assertEqual(lock["version"], e.VERSION)
        self.assertEqual(lock["ruleVersion"], e.RULE_VERSION); self.assertEqual(tuple(lock["scoredMetrics"]), e.METRICS)
        swift = (here.parents[2] / "Sources/NumiVivoKit/Omics/VivoArc2026Contract.swift").read_text()
        self.assertIn(e.COMMIT, swift)
        for name in e.METRICS: self.assertIn('"' + name + '"', swift)


if __name__ == "__main__":
    unittest.main()
