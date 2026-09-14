#!/usr/bin/env python3
"""Independently score the frozen cross-study native NB2 predictions.

The scorer reads only the frozen source count matrices, cohort metadata, native
model/report JSON, and the fold manifest.  It reconstructs source log1p-CPM
rows without NumPy, verifies native training pseudobulk rows, checks the
legacy baselines and NB response transformation, and scores every held-out
treated row.  It never supplies held-out treated counts to model fitting.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import struct
import zipfile
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def reshape(values: list[Any], shape: tuple[int, ...]) -> Any:
    if len(shape) == 1:
        if len(values) != shape[0]:
            raise ValueError("invalid one-dimensional array length")
        return values
    if len(shape) != 2:
        raise ValueError(f"unsupported NPY shape {shape}")
    rows, columns = shape
    if len(values) != rows * columns:
        raise ValueError("invalid two-dimensional array length")
    return [values[row * columns : (row + 1) * columns] for row in range(rows)]


def load_npy(payload: bytes) -> Any:
    if payload[:6] != b"\x93NUMPY":
        raise ValueError("not a NumPy array")
    major, minor = payload[6], payload[7]
    if (major, minor) == (1, 0):
        header_length = struct.unpack_from("<H", payload, 8)[0]
        data_offset = 10 + header_length
    elif (major, minor) in ((2, 0), (3, 0)):
        header_length = struct.unpack_from("<I", payload, 8)[0]
        data_offset = 12 + header_length
    else:
        raise ValueError(f"unsupported NPY version {(major, minor)}")
    header = ast.literal_eval(payload[10 if major == 1 else 12 : data_offset].decode("latin1").strip())
    dtype = header["descr"]
    shape = tuple(header["shape"])
    if header["fortran_order"]:
        raise ValueError("Fortran-order arrays are not supported")
    count = math.prod(shape)
    raw = payload[data_offset:]
    if dtype in ("<u8", "|u8"):
        if len(raw) < count * 8:
            raise ValueError("truncated uint64 array")
        values = [struct.unpack_from("<Q", raw, offset * 8)[0] for offset in range(count)]
    elif dtype.startswith("<U"):
        width = int(dtype[2:])
        item_size = width * 4
        if len(raw) < count * item_size:
            raise ValueError("truncated Unicode array")
        values = [
            raw[offset * item_size : (offset + 1) * item_size].decode("utf-32le").rstrip("\x00")
            for offset in range(count)
        ]
    elif dtype.startswith("|S"):
        width = int(dtype[2:])
        if len(raw) < count * width:
            raise ValueError("truncated byte-string array")
        values = [
            raw[offset * width : (offset + 1) * width].split(b"\x00", 1)[0].decode("utf-8")
            for offset in range(count)
        ]
    else:
        raise ValueError(f"unsupported NPY dtype {dtype}")
    return reshape(values, shape)


def load_npz(path: Path) -> dict[str, Any]:
    with zipfile.ZipFile(path) as archive:
        return {
            name.removesuffix(".npy"): load_npy(archive.read(name))
            for name in archive.namelist()
            if name.endswith(".npy")
        }


def max_difference(actual: list[float], expected: list[float]) -> float:
    if len(actual) != len(expected):
        raise ValueError(f"vector lengths differ: {len(actual)} != {len(expected)}")
    return max((abs(float(a) - float(b)) for a, b in zip(actual, expected)), default=0.0)


def require_close(actual: list[float], expected: list[float], tolerance: float, label: str) -> float:
    difference = max_difference(actual, expected)
    if difference > tolerance:
        raise ValueError(f"{label} differs by {difference} (limit {tolerance})")
    return difference


def arithmetic_mean(values: list[float]) -> float:
    if not values:
        raise ValueError("cannot average an empty sequence")
    return float(sum(values) / len(values))


def root_mean_square(values: list[float]) -> float:
    return math.sqrt(arithmetic_mean([value * value for value in values]))


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return float(ordered[middle])
    return float((ordered[middle - 1] + ordered[middle]) / 2.0)


def pearson(left: list[float], right: list[float]) -> float | None:
    left_mean = arithmetic_mean(left)
    right_mean = arithmetic_mean(right)
    centered_left = [value - left_mean for value in left]
    centered_right = [value - right_mean for value in right]
    denominator = math.sqrt(
        sum(value * value for value in centered_left) * sum(value * value for value in centered_right)
    )
    if denominator == 0:
        return None
    return float(sum(a * b for a, b in zip(centered_left, centered_right)) / denominator)


def solve_many(matrix: list[list[float]], right_hand_side: list[list[float]]) -> list[list[float]]:
    """Solve A X = B for the small donor-by-donor ridge system."""
    size = len(matrix)
    if size == 0 or len(right_hand_side) != size:
        raise ValueError("invalid linear-system dimensions")
    columns = len(right_hand_side[0])
    augmented = [list(matrix[row]) + list(right_hand_side[row]) for row in range(size)]
    for pivot in range(size):
        pivot_row = max(range(pivot, size), key=lambda row: abs(augmented[row][pivot]))
        if abs(augmented[pivot_row][pivot]) < 1e-12:
            raise ValueError("singular ridge system")
        augmented[pivot], augmented[pivot_row] = augmented[pivot_row], augmented[pivot]
        scale = augmented[pivot][pivot]
        augmented[pivot] = [value / scale for value in augmented[pivot]]
        for row in range(size):
            if row == pivot:
                continue
            factor = augmented[row][pivot]
            if factor == 0:
                continue
            augmented[row] = [
                value - factor * pivot_value for value, pivot_value in zip(augmented[row], augmented[pivot])
            ]
    return [row[size : size + columns] for row in augmented]


def source_log_rows(
    counts: list[list[int]], panel_columns: list[int]
) -> tuple[list[list[float]], list[int]]:
    logs: list[list[float]] = []
    totals: list[int] = []
    for row in counts:
        total = int(sum(row))
        if total <= 0:
            raise ValueError("empty source library")
        totals.append(total)
        logs.append([math.log1p(row[column] * 1_000_000.0 / total) for column in panel_columns])
    return logs, totals


def reconstruct_model(
    model: dict[str, Any],
    counts: list[list[int]],
    logs: list[list[float]],
    sample_records: list[dict[str, Any]],
    training_indices: list[int],
    panel_columns: list[int],
    label: str,
) -> dict[str, Any]:
    donors = sorted({sample_records[index]["donorID"] for index in training_indices})
    if model["trainingDonors"] != donors:
        raise ValueError(f"{label} training donor list differs")
    pairs: list[tuple[int, int]] = []
    for donor in donors:
        controls = [
            index
            for index in training_indices
            if sample_records[index]["donorID"] == donor and sample_records[index]["condition"] == "control"
        ]
        treated = [
            index
            for index in training_indices
            if sample_records[index]["donorID"] == donor and sample_records[index]["condition"] == "IFNB"
        ]
        if len(controls) != 1 or len(treated) != 1:
            raise ValueError(f"{label} donor {donor} is not a single paired control/IFNB")
        pairs.append((controls[0], treated[0]))

    paired_counts = [
        [counts[control][column] + counts[treated][column] for column in panel_columns]
        for control, treated in pairs
    ]
    selected = [
        feature
        for feature in range(len(panel_columns))
        if sum(row[feature] for row in paired_counts) >= 10
        and sum(row[feature] > 0 for row in paired_counts) >= 2
    ]
    if model["selectedFeatureIndices"] != selected:
        raise ValueError(f"{label} selected feature indices differ")

    # Logs are already in the explicit panel order; panel_columns is used
    # only for indexing the source count matrix above.
    controls = [list(logs[control]) for control, _ in pairs]
    responses = [
        [treated_value - control_value for treated_value, control_value in zip(logs[treated], logs[control])]
        for control, treated in pairs
    ]
    centers = [arithmetic_mean([row[feature] for row in controls]) for feature in selected]
    scales: list[float] = []
    for offset, feature in enumerate(selected):
        variance = arithmetic_mean([(row[feature] - centers[offset]) ** 2 for row in controls])
        scale = math.sqrt(variance)
        if scale == 0 or all(row[feature] == controls[0][feature] for row in controls):
            scale = 1.0
        scales.append(scale)
    contexts = [
        [
            (row[feature] - centers[offset]) / scales[offset] / math.sqrt(len(selected))
            for offset, feature in enumerate(selected)
        ]
        for row in controls
    ]
    mean_response = [
        arithmetic_mean([row[feature] for row in responses]) for feature in range(len(panel_columns))
    ]
    median_response = [
        median([row[feature] for row in responses]) for feature in range(len(panel_columns))
    ]
    kernel = [
        [sum(left * right for left, right in zip(contexts[row], contexts[column])) for column in range(len(donors))]
        for row in range(len(donors))
    ]
    right_hand_side = [
        [responses[row][feature] - mean_response[feature] for feature in range(len(panel_columns))]
        for row in range(len(donors))
    ]
    dual = solve_many(
        [[kernel[row][column] + (1.0 if row == column else 0.0) for column in range(len(donors))]
         for row in range(len(donors))],
        right_hand_side,
    )

    differences = {
        "selectedFeatureIndices": selected == model["selectedFeatureIndices"],
        "contextCenters": require_close(model["contextCenters"], centers, 1e-8, f"{label}/contextCenters"),
        "contextScales": require_close(model["contextScales"], scales, 1e-8, f"{label}/contextScales"),
        "meanResponse": require_close(model["meanResponse"], mean_response, 1e-8, f"{label}/meanResponse"),
        "medianResponse": require_close(model["medianResponse"], median_response, 1e-8, f"{label}/medianResponse"),
    }
    model_contexts = model["contexts"]
    if len(model_contexts) != len(contexts):
        raise ValueError(f"{label} context row count differs")
    differences["contexts"] = max(
        (max_difference(actual, expected) for actual, expected in zip(model_contexts, contexts)),
        default=0.0,
    )
    if differences["contexts"] > 1e-8:
        raise ValueError(f"{label}/contexts differs by {differences['contexts']}")
    model_dual = model["dualCoefficients"]
    if len(model_dual) != len(dual):
        raise ValueError(f"{label} dual row count differs")
    differences["dualCoefficients"] = max(
        (max_difference(actual, expected) for actual, expected in zip(model_dual, dual)),
        default=0.0,
    )
    if differences["dualCoefficients"] > 1e-8:
        raise ValueError(f"{label}/dualCoefficients differs by {differences['dualCoefficients']}")
    return {
        "selected": selected,
        "centers": centers,
        "scales": scales,
        "contexts": contexts,
        "mean": mean_response,
        "median": median_response,
        "dual": dual,
        "differences": differences,
    }


def sparse_training_check(
    report_path: Path,
    study: str,
    source_features: list[str],
    counts: list[list[int]],
    sample_records: list[dict[str, Any]],
    training_indices: list[int],
    label: str,
) -> None:
    pseudobulk = load_json(report_path)["pseudobulk"]
    if pseudobulk["featureIDs"] != source_features:
        raise ValueError(f"{label} training feature IDs differ")
    matrix = pseudobulk["matrix"]
    groups = pseudobulk["groups"]
    if len(groups) != len(training_indices):
        raise ValueError(f"{label} training group count differs")
    offsets = matrix["rowOffsets"]
    indices = matrix["featureIndices"]
    values = matrix["counts"]
    if len(offsets) != len(groups) + 1 or offsets[-1] != len(indices) or len(indices) != len(values):
        raise ValueError(f"{label} invalid sparse matrix")
    sample_rows = {record["id"]: row for row, record in enumerate(sample_records)}
    reported_ids: list[str] = []
    for row, group in enumerate(groups):
        sample_ids = group["sampleIDs"]
        if len(sample_ids) != 1 or sample_ids[0] not in sample_rows:
            raise ValueError(f"{label} unknown training sample")
        sample_id = sample_ids[0]
        reported_ids.append(sample_id)
        source_row = sample_rows[sample_id]
        begin, end = offsets[row], offsets[row + 1]
        actual_indices = indices[begin:end]
        actual_values = values[begin:end]
        expected = [(column, value) for column, value in enumerate(counts[source_row]) if value]
        if list(zip(actual_indices, actual_values)) != expected:
            raise ValueError(f"{label} sparse training row differs for {sample_id}")
    expected_ids = [sample_records[index]["id"] for index in training_indices]
    if sorted(reported_ids) != sorted(expected_ids):
        raise ValueError(f"{label} training sample coverage differs")


parser = argparse.ArgumentParser()
parser.add_argument("--inputs", type=Path, required=True)
parser.add_argument("--native", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)

input_freeze = load_json(args.inputs / "input-freeze.json")
for relative, digest in input_freeze["files"].items():
    if sha256(args.inputs / relative) != digest:
        raise ValueError(f"input freeze mismatch: {relative}")
prediction_freeze = load_json(args.native / "prediction-freeze.json")
native_manifest = load_json(args.native / "manifest.json")
if prediction_freeze["completedFolds"] != 26 or prediction_freeze["failures"]:
    raise ValueError("native prediction freeze is incomplete")
if native_manifest["status"] != "passed-complete-native-nb2":
    raise ValueError("native NB2 manifest is not a complete pass")

panel = load_json(args.inputs / "panel.json")
folds = load_json(args.inputs / "folds.json")
cohorts = {
    study: load_json(args.inputs / f"{study}-cohort.json")
    for study in ("Kang", "HIRISA")
}
source_data: dict[str, dict[str, Any]] = {}
for study, cohort in cohorts.items():
    arrays = load_npz(args.inputs / f"{study}-source-counts.npz")
    counts = arrays["counts"]
    feature_ids = arrays["featureIDs"]
    sample_ids = arrays["sampleIDs"]
    records = cohort["samples"]
    if sample_ids != [record["id"] for record in records]:
        raise ValueError(f"{study} source sample IDs differ from cohort")
    if len(counts) != len(records) or len(counts[0]) != len(feature_ids):
        raise ValueError(f"{study} source matrix shape differs from metadata")
    if [int(sum(row)) for row in counts] != [int(value) for value in cohort["libraryCounts"]]:
        raise ValueError(f"{study} source library totals differ from metadata")
    suffix_index: dict[str, int] = {}
    for index, feature in enumerate(feature_ids):
        suffix = str(feature).split("|", 1)[-1]
        if suffix in suffix_index:
            raise ValueError(f"{study} duplicate source feature suffix {suffix}")
        suffix_index[suffix] = index
    panel_columns = []
    for feature in panel:
        suffix = feature.split("|", 1)[-1]
        if suffix not in suffix_index:
            raise ValueError(f"{study} missing shared feature {feature}")
        panel_columns.append(suffix_index[suffix])
    logs, totals = source_log_rows(counts, panel_columns)
    source_data[study] = {
        "counts": counts,
        "featureIDs": feature_ids,
        "records": records,
        "sampleRows": {record["id"]: row for row, record in enumerate(records)},
        "panelColumns": panel_columns,
        "logs": logs,
        "totals": totals,
        "sha256": sha256(args.inputs / f"{study}-source-counts.npz"),
    }

baselines = ("noChange", "meanResponse", "medianResponse", "contextRidge", "negativeBinomialEffect")
scores: list[dict[str, Any]] = []
checks: list[dict[str, Any]] = []
fold_seen: set[str] = set()
for fold in folds:
    fold_id = fold["id"]
    fold_root = args.native / fold_id
    if fold_id in fold_seen:
        raise ValueError(f"duplicate fold {fold_id}")
    fold_seen.add(fold_id)
    train = fold["trainingStudy"]
    query = fold["queryStudy"]
    train_data = source_data[train]
    query_data = source_data[query]
    model_path = fold_root / "model/model.json"
    report_path = fold_root / "prediction/report.json"
    training_report_path = fold_root / "model/training/report.json"
    for path in (model_path, report_path, training_report_path):
        relative = str(path.relative_to(args.native))
        expected_digest = prediction_freeze["files"].get(relative)
        if expected_digest is None or sha256(path) != expected_digest:
            raise ValueError(f"prediction freeze mismatch: {relative}")
    model = load_json(model_path)
    report = load_json(report_path)
    if model["featureIDs"] != panel or report["featureIDs"] != panel:
        raise ValueError(f"{fold_id} panel differs")
    if len(report["predictions"]) != 1:
        raise ValueError(f"{fold_id} expected one held-out prediction")
    nb = model.get("negativeBinomial")
    if not isinstance(nb, dict):
        raise ValueError(f"{fold_id} has no NB2 component")
    tested = list(nb["testedFeatureIndices"])
    if tested != sorted(set(tested)) or any(index < 0 or index >= len(panel) for index in tested):
        raise ValueError(f"{fold_id} invalid NB2 feature indices")
    if len(nb["logEffects"]) != len(panel) or len(nb["statuses"]) != len(panel):
        raise ValueError(f"{fold_id} NB2 arrays differ from panel")

    training_indices = list(fold["trainingIndices"])
    sparse_training_check(
        training_report_path,
        train,
        train_data["featureIDs"],
        train_data["counts"],
        train_data["records"],
        training_indices,
        fold_id,
    )
    reconstructed = reconstruct_model(
        model,
        train_data["counts"],
        train_data["logs"],
        train_data["records"],
        training_indices,
        train_data["panelColumns"],
        fold_id,
    )

    target_records = query_data["records"]
    query_record = target_records[fold["queryIndex"]]
    treated_record = target_records[fold["scoringIndex"]]
    if (
        query_record["condition"] != "control"
        or treated_record["condition"] != "IFNB"
        or not (
            query_record["donorID"]
            == treated_record["donorID"]
            == fold["heldOutDonor"]
        )
    ):
        raise ValueError(f"{fold_id} query/scoring metadata is inconsistent")
    prediction = report["predictions"][0]
    group = prediction["group"]
    if group["donorID"] != fold["heldOutDonor"] or group["sampleIDs"] != [query_record["id"]]:
        raise ValueError(f"{fold_id} prediction group differs from query metadata")
    control = query_data["logs"][fold["queryIndex"]]
    truth = query_data["logs"][fold["scoringIndex"]]
    require_close(prediction["control"], control, 1e-10, f"{fold_id}/control")
    if int(prediction["libraryCounts"]) != query_data["totals"][fold["queryIndex"]]:
        raise ValueError(f"{fold_id} control library total differs")

    selected = reconstructed["selected"]
    context = reconstructed["contexts"]
    query_context = [
        (control[feature] - reconstructed["centers"][offset])
        / reconstructed["scales"][offset]
        / math.sqrt(len(selected))
        for offset, feature in enumerate(selected)
    ]
    ridge_weights = [sum(query_context[offset] * row[offset] for offset in range(len(selected))) for row in context]
    ridge = [
        reconstructed["mean"][feature]
        + sum(ridge_weights[donor] * reconstructed["dual"][donor][feature] for donor in range(len(context)))
        for feature in range(len(panel))
    ]
    # The training and query cohorts share the explicit panel order.  The
    # expression above is deliberately written in panel coordinates.
    estimates = {estimate["baseline"]: estimate for estimate in prediction["estimates"]}
    if set(estimates) != set(baselines):
        raise ValueError(f"{fold_id} expected five response baselines")
    references: dict[str, list[float]] = {
        "noChange": [0.0] * len(panel),
        "meanResponse": reconstructed["mean"],
        "medianResponse": reconstructed["median"],
        "contextRidge": ridge,
    }
    nb_raw: list[float] = []
    for feature, status in enumerate(nb["statuses"]):
        effect = nb["logEffects"][feature]
        if status == "tested" and effect is not None:
            baseline_cpm = math.expm1(control[feature])
            nb_raw.append(math.log1p(baseline_cpm * math.exp(float(effect))) - control[feature])
        else:
            nb_raw.append(0.0)
    references["negativeBinomialEffect"] = nb_raw
    observed_change = [actual - baseline for actual, baseline in zip(truth, control)]
    for baseline in baselines:
        estimate = estimates[baseline]
        raw = [float(value) for value in estimate["unclippedResponse"]]
        expected_raw = references[baseline]
        require_close(raw, expected_raw, 2e-8, f"{fold_id}/{baseline}/unclippedResponse")
        expected_treated = [max(0.0, control[feature] + raw[feature]) for feature in range(len(panel))]
        actual_treated = [float(value) for value in estimate["predictedTreated"]]
        require_close(actual_treated, expected_treated, 2e-8, f"{fold_id}/{baseline}/predictedTreated")
        expected_response = [actual - baseline for actual, baseline in zip(expected_treated, control)]
        require_close(
            [float(value) for value in estimate["predictedResponse"]],
            expected_response,
            2e-8,
            f"{fold_id}/{baseline}/predictedResponse",
        )
        expected_sum = sum(math.expm1(value) for value in expected_treated)
        if abs(float(estimate["impliedCPMSum"]) - expected_sum) > 1e-4:
            raise ValueError(f"{fold_id}/{baseline}/impliedCPMSum differs")
        errors = [actual - observed for actual, observed in zip(actual_treated, truth)]
        identified_errors = [errors[index] for index in tested]
        if baseline == "negativeBinomialEffect" and estimate.get("availableFeatureIndices") != tested:
            raise ValueError(f"{fold_id} NB available feature list differs")
        scores.append(
            {
                "fold": fold_id,
                "mode": fold["mode"],
                "trainingStudy": train,
                "queryStudy": query,
                "donor": fold["heldOutDonor"],
                "baseline": baseline,
                "features": len(panel),
                "identifiedFeatures": len(tested) if baseline == "negativeBinomialEffect" else len(panel),
                "responseRMSE": root_mean_square(errors),
                "responseMAE": arithmetic_mean([abs(value) for value in errors]),
                "responsePearson": pearson(
                    [actual - baseline for actual, baseline in zip(actual_treated, control)],
                    observed_change,
                ),
                "identifiedResponseRMSE": root_mean_square(identified_errors),
                "identifiedResponseMAE": arithmetic_mean([abs(value) for value in identified_errors]),
                "clippedFeatures": sum(1 for value in raw if value < 0),
            }
        )
    checks.append(
        {
            "fold": fold_id,
            "trainingStudy": train,
            "queryStudy": query,
            "mode": fold["mode"],
            "donor": fold["heldOutDonor"],
            "modelSHA256": sha256(model_path),
            "testedFeatures": len(tested),
            "modelMaximumDifference": max(
                value
                for key, value in reconstructed["differences"].items()
                if key != "selectedFeatureIndices"
            ),
        }
    )

if len(fold_seen) != len(folds) or len(folds) != 26:
    raise ValueError("unexpected fold coverage")

summaries: list[dict[str, Any]] = []
for target in ("Kang", "HIRISA"):
    for mode in ("cross", "within"):
        subset = [row for row in scores if row["queryStudy"] == target and row["mode"] == mode]
        by_baseline = {
            baseline: {
                "meanMAE": arithmetic_mean([row["responseMAE"] for row in subset if row["baseline"] == baseline]),
                "meanRMSE": arithmetic_mean([row["responseRMSE"] for row in subset if row["baseline"] == baseline]),
                "meanIdentifiedMAE": arithmetic_mean(
                    [row["identifiedResponseMAE"] for row in subset if row["baseline"] == baseline]
                ),
                "meanIdentifiedRMSE": arithmetic_mean(
                    [row["identifiedResponseRMSE"] for row in subset if row["baseline"] == baseline]
                ),
            }
            for baseline in baselines
        }
        nb_rows = [row for row in subset if row["baseline"] == "negativeBinomialEffect"]
        donor_scores = {
            donor: {row["baseline"]: row["responseRMSE"] for row in subset if row["donor"] == donor}
            for donor in sorted({row["donor"] for row in subset})
        }
        summaries.append(
            {
                "queryStudy": target,
                "mode": mode,
                "donors": len(donor_scores),
                "mean": by_baseline,
                "negativeBinomialIdentifiedMeanMAE": arithmetic_mean(
                    [row["identifiedResponseMAE"] for row in nb_rows]
                ),
                "negativeBinomialIdentifiedMeanRMSE": arithmetic_mean(
                    [row["identifiedResponseRMSE"] for row in nb_rows]
                ),
                "negativeBinomialBetterThanNoChange": sum(
                    row["responseRMSE"] < donor_scores[row["donor"]]["noChange"] for row in nb_rows
                ),
                "negativeBinomialBetterThanMean": sum(
                    row["responseRMSE"] < donor_scores[row["donor"]]["meanResponse"] for row in nb_rows
                ),
                "negativeBinomialGainPercent": 100.0
                * (
                    1.0
                    - by_baseline["negativeBinomialEffect"]["meanRMSE"]
                    / by_baseline["noChange"]["meanRMSE"]
                ),
            }
        )

checks_manifest = {
    "status": "passed-independent-source-npz-scoring",
    "sourceInputFreezeSHA256": sha256(args.inputs / "input-freeze.json"),
    "predictionFreezeSHA256": sha256(args.native / "prediction-freeze.json"),
    "sourceCountsSHA256": {study: source_data[study]["sha256"] for study in source_data},
    "folds": len(folds),
    "predictions": len(scores) // len(baselines),
    "estimateVectors": len(scores),
    "panelFeatures": len(panel),
    "checks": checks,
}
write_json(args.out / "scores.json", scores)
write_json(args.out / "summary.json", summaries)
write_json(args.out / "checks.json", checks_manifest)
print(json.dumps(summaries, indent=2, sort_keys=True))
