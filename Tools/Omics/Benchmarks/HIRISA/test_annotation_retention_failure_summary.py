#!/usr/bin/env python3
"""Controls for the failure-targeting annotation-retention ledger."""
import unittest

from summarize_annotation_retention_failures import build_summary, freeze_extension


def comparison(
    stratum,
    label,
    name,
    *,
    sufficient=True,
    sensitive=True,
    rare=False,
    failure=False,
    baseline=0.8,
):
    return {
        "stratum": stratum,
        "label": label,
        "labelName": name,
        "rare": rare,
        "sufficientSupport": sufficient,
        "controlSensitive": sensitive,
        "baselineMeanRecall": baseline,
        "erasureMeanRecall": 0.1,
        "candidateMeanRecall": baseline - (0.2 if failure else 0.01),
        "meanRecallLoss": 0.2 if failure else 0.01,
        "maximumFoldRecallLoss": 0.2 if failure else 0.01,
        "marginFailure": failure,
    }


class FailureSummaryChecks(unittest.TestCase):
    def test_classifies_shared_new_resolved_and_unsupported_records(self):
        ledger = {
            "labels": ["A", "B", "C", "D"],
            "strata": [
                {"preparation": "P0", "treatment": "T0"},
                {"preparation": "P1", "treatment": "T1"},
            ],
        }
        native = {
            "matrix": "native",
            "matrixSHA256": "baseline-hash",
            "comparisons": [
                comparison(0, 0, "A", failure=True),
                comparison(0, 1, "B"),
                comparison(1, 2, "C", failure=True),
                comparison(1, 3, "D", sufficient=False, failure=True),
            ],
        }
        candidate = {
            "matrix": "condition-f1",
            "comparisons": [
                comparison(0, 0, "A", failure=True),
                comparison(0, 1, "B", failure=True),
                comparison(1, 2, "C"),
                comparison(1, 3, "D", sufficient=False, failure=False),
            ],
        }
        summary = build_summary(
            ledger, native, candidate,
            {"ledger": "l", "native": "n", "candidate": "c"},
        )
        self.assertEqual(summary["counts"]["allComparisons"], 4)
        self.assertEqual(summary["counts"]["supportedControlSensitive"], 3)
        self.assertEqual(summary["counts"]["unsupported"], 1)
        self.assertEqual(summary["counts"]["sharedMarginFailures"], 1)
        self.assertEqual(summary["counts"]["newCandidateMarginFailures"], 1)
        self.assertEqual(summary["counts"]["resolvedByCandidate"], 1)
        self.assertEqual(summary["counts"]["unsupportedNoPromotionDecision"], 1)
        self.assertEqual(
            [record["outcome"] for record in summary["records"]],
            [
                "shared_margin_failure",
                "new_candidate_margin_failure",
                "resolved_by_candidate",
                "unsupported_no_promotion_decision",
            ],
        )

    def test_rejects_misaligned_baseline_metrics(self):
        ledger = {
            "labels": ["A"],
            "strata": [{"preparation": "P", "treatment": "T"}],
        }
        native = {"comparisons": [comparison(0, 0, "A")]}
        candidate = {"comparisons": [comparison(0, 0, "A", baseline=0.7)]}
        with self.assertRaisesRegex(ValueError, "baselineMeanRecall"):
            build_summary(ledger, native, candidate, {})

    def test_allows_only_the_pinned_candidate_matrix_freeze_extension(self):
        base = {
            "cells": 10,
            "matrices": {"native": {"SHA256": "native-score", "path": "native"}},
        }
        candidate = {
            "cells": 10,
            "matrices": {
                "native": {"SHA256": "native-score", "path": "native"},
                "condition-f1": {"SHA256": "candidate-score", "path": "candidate"},
            },
        }
        provenance = freeze_extension(
            base, candidate, "condition-f1", "candidate-score"
        )
        self.assertEqual(
            provenance["relationship"], "candidate_freeze_extends_native_freeze"
        )
        bad = dict(candidate)
        bad["cells"] = 11
        with self.assertRaisesRegex(ValueError, "outside matrices"):
            freeze_extension(base, bad, "condition-f1", "candidate-score")


if __name__ == "__main__":
    unittest.main()
