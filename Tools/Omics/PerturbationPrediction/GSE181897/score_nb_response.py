"""Score native negative-binomial response bundles against the frozen source aggregate.

This scorer intentionally uses only the native aggregate report and prediction
reports. It does not infer labels or inspect outcomes while fitting. The
negative-binomial estimate is scored on the full shared panel and, separately,
on the feature indices identified by the fitted NB model.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def load(path: Path) -> Any:
    return json.loads(path.read_text())


def mean(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot average an empty sequence")
    return float(sum(values) / len(values))


def rms(values: list[float]) -> float:
    return math.sqrt(mean([value * value for value in values]))


def group_key(group: dict[str, Any]) -> tuple[str, str]:
    donor = str(group["donorID"]).removeprefix("GSE181897:exp_id:")
    condition = str(group["condition"]).removeprefix("source-code:")
    return donor, condition


parser = argparse.ArgumentParser()
parser.add_argument("--source-aggregate", type=Path, required=True)
parser.add_argument("--evidence", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)

source_report_path = args.source_aggregate / "report.json"
source = load(source_report_path)
pseudobulk = source["pseudobulk"]
source_features = pseudobulk["featureIDs"]
matrix = pseudobulk["matrix"]
groups = pseudobulk["groups"]
if matrix["cellCount"] != len(groups):
    raise ValueError("aggregate row count does not match group count")
if matrix["featureCount"] != len(source_features):
    raise ValueError("aggregate feature count does not match feature IDs")

source_index: dict[str, int] = {}
for index, feature in enumerate(source_features):
    if feature in source_index:
        raise ValueError(f"duplicate source feature {feature}")
    source_index[feature] = index

lookup = {group_key(group): index for index, group in enumerate(groups)}
if len(lookup) != len(groups):
    raise ValueError("duplicate source donor/condition group")

row_offsets = matrix["rowOffsets"]
feature_indices = matrix["featureIndices"]
counts = matrix["counts"]
if len(row_offsets) != len(groups) + 1:
    raise ValueError("invalid aggregate row offsets")
if row_offsets[-1] != len(counts) or len(counts) != len(feature_indices):
    raise ValueError("invalid aggregate sparse vectors")


source_log_cache: dict[int, tuple[list[float], int]] = {}


def source_log_row(row: int, panel: list[str]) -> tuple[list[float], int]:
    cached = source_log_cache.get(row)
    if cached is not None:
        return cached
    start, end = row_offsets[row], row_offsets[row + 1]
    total = sum(counts[start:end])
    if total <= 0:
        raise ValueError(f"empty source library at row {row}")
    row_counts = dict(zip(feature_indices[start:end], counts[start:end]))
    values = [0.0] * len(panel)
    for offset, panel_feature in enumerate(panel):
        source_feature = panel_feature.split("|", 1)[-1]
        source_column = source_index.get(source_feature)
        if source_column is None:
            continue
        value = row_counts.get(source_column, 0)
        values[offset] = math.log1p(value * 1_000_000.0 / total)
    cached = (values, total)
    source_log_cache[row] = cached
    return cached


origins = ("Kang", "HIRISA")
baselines = ("noChange", "meanResponse", "medianResponse", "contextRidge", "negativeBinomialEffect")
scores: list[dict[str, Any]] = []
checks: list[dict[str, Any]] = []

for origin in origins:
    model_path = args.evidence / f"model-{origin}-nb" / "model.json"
    model = load(model_path)
    panel = model["featureIDs"]
    if len(panel) != 11_800:
        raise ValueError(f"{origin} has unexpected panel length")
    if not all(feature.startswith("symbol|") for feature in panel):
        raise ValueError(f"{origin} panel is not symbol-qualified")
    nb = model.get("negativeBinomial")
    if not isinstance(nb, dict):
        raise ValueError(f"{origin} model has no negative-binomial component")
    tested = list(nb["testedFeatureIndices"])
    tested_set = set(tested)
    if sorted(tested) != tested or len(tested_set) != len(tested):
        raise ValueError(f"{origin} NB feature indices are not sorted and unique")
    if any(index < 0 or index >= len(panel) for index in tested):
        raise ValueError(f"{origin} NB feature index is out of range")

    seen: list[str] = []
    origin_scores: list[dict[str, Any]] = []
    for report_path in sorted(args.evidence.glob(f"prediction-{origin}-*-nb/report.json")):
        report = load(report_path)
        if report["featureIDs"] != panel:
            raise ValueError(f"{report_path} panel differs from model")
        for prediction in report["predictions"]:
            donor = str(prediction["group"]["donorID"]).removeprefix("GSE181897:exp_id:")
            seen.append(donor)
            control_row = lookup[(donor, "C")]
            treated_row = lookup[(donor, "B")]
            control, control_total = source_log_row(control_row, panel)
            truth, _ = source_log_row(treated_row, panel)
            native_control = prediction["control"]
            if len(native_control) != len(panel):
                raise ValueError(f"{report_path} control vector length differs")
            control_difference = max(abs(a - b) for a, b in zip(native_control, control))
            if control_difference > 1e-10:
                raise ValueError(f"{report_path} control differs by {control_difference}")
            if int(prediction["libraryCounts"]) != control_total:
                raise ValueError(f"{report_path} control library total differs")
            estimates = {estimate["baseline"]: estimate for estimate in prediction["estimates"]}
            if set(estimates) != set(baselines):
                raise ValueError(f"{report_path} does not contain the five expected baselines")
            for baseline in baselines:
                estimate = estimates[baseline]
                predicted = [float(value) for value in estimate["predictedTreated"]]
                if len(predicted) != len(panel):
                    raise ValueError(f"{report_path} predicted vector length differs")
                errors = [value - actual for value, actual in zip(predicted, truth)]
                identified_errors = [errors[index] for index in tested_set]
                raw = estimate["unclippedResponse"]
                record: dict[str, Any] = {
                    "origin": origin,
                    "donor": donor,
                    "baseline": baseline,
                    "controlCells": len(groups[control_row]["sourceCellIndices"]),
                    "treatedCells": len(groups[treated_row]["sourceCellIndices"]),
                    "availableFeatures": len(tested) if baseline == "negativeBinomialEffect" else len(panel),
                    "clippedFeatures": sum(1 for value in raw if float(value) < 0),
                    "responseMAE": mean([abs(value) for value in errors]),
                    "responseRMSE": rms(errors),
                }
                if baseline == "negativeBinomialEffect":
                    record["identifiedResponseMAE"] = mean([abs(value) for value in identified_errors])
                    record["identifiedResponseRMSE"] = rms(identified_errors)
                    available = estimate.get("availableFeatureIndices")
                    if available != tested:
                        raise ValueError(f"{report_path} NB available feature list differs from model")
                origin_scores.append(record)

    expected_donors = sorted(
        (donor for donor, condition in lookup if condition == "C" and (donor, "B") in lookup),
        key=int,
    )
    if sorted(seen, key=int) != expected_donors or len(seen) != len(set(seen)):
        raise ValueError(f"{origin} donor coverage differs from source cohort")
    scores.extend(origin_scores)
    checks.append(
        {
            "origin": origin,
            "donors": len(seen),
            "predictions": len(seen),
            "features": len(panel),
            "identifiedFeatures": len(tested),
            "modelSHA256": sha256(model_path),
        }
    )

summaries: list[dict[str, Any]] = []
for origin in origins:
    rows = [row for row in scores if row["origin"] == origin]
    by_baseline = {
        baseline: {
            "meanMAE": mean([row["responseMAE"] for row in rows if row["baseline"] == baseline]),
            "meanRMSE": mean([row["responseRMSE"] for row in rows if row["baseline"] == baseline]),
        }
        for baseline in baselines
    }
    nb_rows = [row for row in rows if row["baseline"] == "negativeBinomialEffect"]
    by_donor = {
        donor: {row["baseline"]: row["responseRMSE"] for row in rows if row["donor"] == donor}
        for donor in sorted({row["donor"] for row in rows}, key=int)
    }
    summaries.append(
        {
            "origin": origin,
            "donors": len(by_donor),
            "mean": by_baseline,
            "negativeBinomialIdentifiedMeanMAE": mean([row["identifiedResponseMAE"] for row in nb_rows]),
            "negativeBinomialIdentifiedMeanRMSE": mean([row["identifiedResponseRMSE"] for row in nb_rows]),
            "negativeBinomialBetterThanNoChange": sum(
                row["responseRMSE"] < by_donor[row["donor"]]["noChange"] for row in nb_rows
            ),
            "negativeBinomialBetterThanMean": sum(
                row["responseRMSE"] < by_donor[row["donor"]]["meanResponse"] for row in nb_rows
            ),
            "negativeBinomialGainPercent": 100.0
            * (1.0 - by_baseline["negativeBinomialEffect"]["meanRMSE"] / by_baseline["noChange"]["meanRMSE"]),
        }
    )

manifest = {
    "status": "passed-independent-source-aggregate-scoring",
    "sourceAggregateSHA256": sha256(source_report_path),
    "sourceRows": len(groups),
    "sourceFeatures": len(source_features),
    "panelFeatures": 11_800,
    "origins": list(origins),
    "predictions": len(scores) // len(baselines),
    "estimateVectors": len(scores),
    "checks": checks,
}
(args.out / "scores.json").write_text(json.dumps(scores, indent=2, sort_keys=True) + "\n")
(args.out / "summary.json").write_text(json.dumps(summaries, indent=2, sort_keys=True) + "\n")
(args.out / "checks.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
print(json.dumps(summaries, indent=2, sort_keys=True))
