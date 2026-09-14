#!/usr/bin/env python3
"""Run a provisional Adamson held-target sensitivity experiment.

The deposited control labels are still unresolved against the primary roster.
This driver therefore declares its control pool in the output provenance and
never changes the authoritative role audit.  It uses the selected native
pseudobulk report, aggregates all deposited guide rows for each candidate
target, runs the existing native target-kernel CLI, freezes every prediction,
then scores the frozen reports against the held aggregate.
"""
import argparse
import gzip
import hashlib
import json
import math
import shutil
import subprocess
import time
from pathlib import Path


NAMESPACE = "GO-direct-BP-MF-CC-MyGene-20260906"
CONTROL_LABELS = ("63(mod)_pBA580", "Gal4-4(mod)_pBA582")
CONTEXT = {
    "id": "Adamson2016-K562-upr-provisional-role-sensitivity",
    "organism": "NCBITaxon:9606",
    "featureNamespace": "Ensembl-gene",
    "perturbationNamespace": "Adamson2016-guide-target",
    "countUnit": "umiCount",
}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    return {"bytes": list(bytes.fromhex(sha(path)))}


def write(path, value):
    path.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    )


def compress(source, destination):
    with source.open("rb") as raw, destination.open("wb") as output:
        with gzip.GzipFile(filename="", fileobj=output, mode="wb", mtime=0) as compressed:
            shutil.copyfileobj(raw, compressed)


def copy_fingerprint(value):
    return {"bytes": list(value["bytes"])}


def merge_rows(matrix, rows):
    """Merge sorted CSR rows into one sorted positive-count row."""
    offsets = matrix["rowOffsets"]
    indices = matrix["featureIndices"]
    counts = matrix["counts"]
    merged = {}
    for row in rows:
        for cursor in range(offsets[row], offsets[row + 1]):
            column = indices[cursor]
            merged[column] = merged.get(column, 0) + counts[cursor]
    columns = sorted(merged)
    return columns, [merged[column] for column in columns]


def subset(training, keep):
    original = training["matrix"]
    offsets, columns, values = [0], [], []
    for row in [0, *[index + 1 for index in keep]]:
        start, end = original["rowOffsets"][row : row + 2]
        columns.extend(original["featureIndices"][start:end])
        values.extend(original["counts"][start:end])
        offsets.append(len(values))
    selection = dict(
        training["selection"],
        targets=[training["selection"]["targets"][index] for index in keep],
    )
    return dict(
        training,
        selection=selection,
        matrix=dict(
            original,
            cellCount=len(keep) + 1,
            rowOffsets=offsets,
            featureIndices=columns,
            counts=values,
        ),
    )


def row_cpm(matrix, row, feature_count):
    start, end = matrix["rowOffsets"][row : row + 2]
    total = sum(matrix["counts"][start:end])
    if total <= 0:
        raise AssertionError(f"empty aggregate row {row}")
    result = [0.0] * feature_count
    for cursor in range(start, end):
        result[matrix["featureIndices"][cursor]] = (
            matrix["counts"][cursor] / total * 1_000_000.0
        )
    return result


def descriptor(entry, annotations):
    target = entry["target"]
    terms = annotations.get(target, {}).get("terms", [])
    if entry["status"] == "supported":
        status = "available"
    elif entry["status"] in ("no-usable-GO-annotations", "annotation-not-found"):
        status = "noData"
        terms = []
    else:
        status = "unresolvedIdentity"
        terms = []
    value = {"targetID": target, "status": status, "terms": sorted(terms)}
    if "sourceFeatureID" in entry:
        value["featureID"] = entry["sourceFeatureID"]
    return value


def run_command(binary, command, log_path, commands, expected=0):
    before = time.monotonic()
    with log_path.open("w") as log:
        result = subprocess.run(
            ["/usr/bin/time", "-l", str(binary.resolve()), *map(str, command)],
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    commands.append(
        {
            "name": log_path.stem,
            "returnCode": result.returncode,
            "expected": expected,
            "elapsedSeconds": time.monotonic() - before,
        }
    )
    if result.returncode != expected:
        raise RuntimeError(f"native command failed: {log_path}")


def metric(predicted, actual):
    if len(predicted) != len(actual):
        raise AssertionError("metric vector dimensions")
    errors = [a - b for a, b in zip(predicted, actual)]
    rmse = math.sqrt(sum(error * error for error in errors) / len(errors))
    mae = sum(abs(error) for error in errors) / len(errors)
    mean_a = sum(actual) / len(actual)
    mean_b = sum(predicted) / len(predicted)
    centered_a = [value - mean_a for value in actual]
    centered_b = [value - mean_b for value in predicted]
    denominator = math.sqrt(
        sum(value * value for value in centered_a)
        * sum(value * value for value in centered_b)
    )
    correlation = (
        sum(a * b for a, b in zip(centered_a, centered_b)) / denominator
        if denominator > 0
        else None
    )
    nonzero = [index for index, value in enumerate(actual) if value != 0]
    sign_accuracy = (
        sum((predicted[index] == 0) == (actual[index] == 0) or predicted[index] * actual[index] > 0 for index in nonzero)
        / len(nonzero)
        if nonzero
        else None
    )
    return {
        "rmse": rmse,
        "mae": mae,
        "correlation": correlation,
        "signAccuracyOnNonzeroActual": sign_accuracy,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("binary", "cohort-report", "cohort-receipt", "coverage", "annotations", "out"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--limit", type=int, default=0, help="run only the first N held targets")
    args = parser.parse_args()
    if args.limit < 0:
        raise ValueError("--limit must be nonnegative")
    args.out.mkdir(parents=True, exist_ok=False)

    report = json.loads(args.cohort_report.read_text())
    receipt = json.loads(args.cohort_receipt.read_text())
    coverage = json.loads(args.coverage.read_text())
    with gzip.open(args.annotations, "rt") as stream:
        annotations = json.load(stream)
    pseudobulk = report["pseudobulk"]
    groups = pseudobulk["groups"]
    matrix = pseudobulk["matrix"]
    by_condition = {group["condition"]: index for index, group in enumerate(groups)}
    if len(by_condition) != len(groups):
        raise AssertionError("duplicate selected pseudobulk conditions")
    for label in CONTROL_LABELS:
        if label not in by_condition:
            raise AssertionError(f"missing provisional control candidate: {label}")
    feature_ids = pseudobulk["featureIDs"]
    feature_index = set(feature_ids)

    selected = []
    for entry in sorted(coverage, key=lambda item: item["target"]):
        guides = [guide for guide in entry["sourceGuides"] if guide in by_condition]
        if not guides:
            continue
        if "sourceFeatureID" in entry and entry["sourceFeatureID"] not in feature_index:
            raise AssertionError(f"descriptor feature is absent from cohort: {entry['target']}")
        selected.append((entry, guides))
    if len(selected) < 2:
        raise AssertionError("too few selected target prefixes")

    control_columns, control_values = merge_rows(
        matrix, [by_condition[label] for label in CONTROL_LABELS]
    )
    target_rows = []
    for entry, guides in selected:
        columns, values = merge_rows(matrix, [by_condition[guide] for guide in guides])
        target_rows.append((entry, guides, columns, values))

    selection = {
        "schemaVersion": 1,
        "context": CONTEXT,
        "controlCondition": "provisional-control-pool-63-Gal4",
        "targets": [{"id": entry["target"], "condition": entry["target"]} for entry, _, _, _ in target_rows],
        "provenance": (
            "Provisional Adamson role-sensitivity run. The selected source labels "
            "63(mod)_pBA580 and Gal4-4(mod)_pBA582 are pooled as a declared "
            "control candidate set; the primary label-to-control mapping remains "
            "unresolved. All deposited guides for each candidate prefix are "
            "aggregated. This output is not an authoritative role assignment or "
            "biological qualification."
        ),
    }
    rows = [(control_columns, control_values)] + [
        (columns, values) for _, _, columns, values in target_rows
    ]
    row_offsets, feature_indices, counts = [0], [], []
    for columns, values in rows:
        feature_indices.extend(columns)
        counts.extend(values)
        row_offsets.append(len(counts))
    training = {
        "schemaVersion": 1,
        "selection": selection,
        "source": copy_fingerprint(receipt["source"]),
        "sourceReport": copy_fingerprint(receipt["report"]),
        "evidence": report["metadata"]["evidence"],
        "featureIDs": feature_ids,
        "matrix": {
            "cellCount": len(rows),
            "featureCount": len(feature_ids),
            "rowOffsets": row_offsets,
            "featureIndices": feature_indices,
            "counts": counts,
        },
    }
    write(args.out / "selection.json", selection)
    write(args.out / "training.json", training)
    compress(args.out / "training.json", args.out / "training.json.gz")
    write(
        args.out / "aggregate-manifest.json",
        {
            "status": "prepared-provisional-role-sensitivity",
            "sourceReportSHA256": sha(args.cohort_report),
            "sourceReceiptSHA256": sha(args.cohort_receipt),
            "coverageSHA256": sha(args.coverage),
            "annotationsSHA256": sha(args.annotations),
            "controlCandidates": list(CONTROL_LABELS),
            "controlMappingAuthoritative": False,
            "selectedCells": len(report["metadata"]["cells"]),
            "selectedGuideGroups": len(groups),
            "selectedTargetPrefixes": len(target_rows),
            "supportedTargetPrefixes": sum(entry["status"] == "supported" for entry, *_ in target_rows),
            "featureCount": len(feature_ids),
            "trainingMatrixNonzeros": len(counts),
            "targetGuides": [
                {"target": entry["target"], "selectedGuides": guides}
                for entry, guides, _, _ in target_rows
            ],
        },
    )

    descriptors = {
        entry["target"]: descriptor(entry, annotations)
        for entry, _, _, _ in target_rows
    }
    all_targets = [entry["target"] for entry, _, _, _ in target_rows]
    commands = []

    fold_targets = all_targets[: args.limit] if args.limit else all_targets
    prediction_fingerprints = []
    started = time.monotonic()
    for held_index, held in enumerate(fold_targets):
        fold = args.out / "folds" / held
        fold.mkdir(parents=True)
        keep = [index for index, target in enumerate(all_targets) if target != held]
        selected_training = subset(training, keep)
        if held in [target["id"] for target in selected_training["selection"]["targets"]]:
            raise AssertionError("held target leaked into training selection")

        mutated = dict(training, matrix=dict(training["matrix"], counts=list(training["matrix"]["counts"])))
        start, end = mutated["matrix"]["rowOffsets"][held_index + 1 : held_index + 3]
        mutated["matrix"]["counts"][start:end] = [1] * (end - start)
        if subset(mutated, keep) != selected_training:
            raise AssertionError("held-out outcome mutation changed training")

        scratch = args.out / "scratch"
        scratch.mkdir()
        write(scratch / "training.json", selected_training)
        plan = {
            "schemaVersion": 1,
            "context": CONTEXT,
            "descriptorNamespace": NAMESPACE,
            "source": fingerprint(args.annotations),
            "provenance": selection["provenance"],
            "targets": [descriptors[target] for target in all_targets if target != held],
            "regularization": 1,
            "maximumWork": 200_000_000,
        }
        query = {
            "schemaVersion": 1,
            "context": CONTEXT,
            "descriptorNamespace": NAMESPACE,
            "source": fingerprint(args.annotations),
            "queries": [{"id": held, "descriptor": descriptors[held]}],
            "maximumWork": 200_000_000,
        }
        write(fold / "plan.json", plan)
        write(fold / "query-input.json", query)
        run_command(
            args.binary,
            ["singlecell-target-kernel-fit", scratch / "training.json", "--plan", fold / "plan.json", "--output", scratch / "model"],
            fold / "fit.log",
            commands,
        )
        run_command(
            args.binary,
            ["singlecell-target-kernel-verify", scratch / "model"],
            fold / "verify.log",
            commands,
        )
        run_command(
            args.binary,
            ["singlecell-target-kernel-predict", scratch / "model", "--plan", fold / "query-input.json", "--output", scratch / "prediction"],
            fold / "predict.log",
            commands,
        )
        run_command(
            args.binary,
            ["singlecell-target-kernel-prediction-verify", scratch / "prediction"],
            fold / "prediction-verify.log",
            commands,
        )
        report_path = scratch / "prediction" / "report.json"
        query_path = scratch / "prediction" / "query.json"
        receipt_path = scratch / "prediction" / "receipt.json"
        compress(report_path, fold / "report.json.gz")
        shutil.copy2(query_path, fold / "query.json")
        shutil.copy2(receipt_path, fold / "receipt.json")
        result = json.loads(report_path.read_text())
        prediction_fingerprints.append(
            {
                "target": held,
                "status": result["queries"][0]["status"],
                "reportSHA256": sha(fold / "report.json.gz"),
                "querySHA256": sha(query_path),
                "receiptSHA256": sha(receipt_path),
                "heldOutcomeMutationSelectedCountsExact": True,
            }
        )
        if held_index == 0:
            retained = args.out / "first-model"
            retained.mkdir()
            for path in (scratch / "model").iterdir():
                compress(path, retained / (path.name + ".gz"))
        shutil.rmtree(scratch)
        write(args.out / "commands.json", commands)
        write(args.out / "predictions-frozen.json", prediction_fingerprints)
        print(f"{held_index + 1}/{len(fold_targets)} {held}: prediction frozen", flush=True)

    write(args.out / "commands.json", commands)
    write(args.out / "predictions-frozen.json", prediction_fingerprints)

    # Scoring starts only after every requested fold has a retained, verified
    # prediction report and its immutable query/receipt hashes are recorded.
    control_cpm = row_cpm(training["matrix"], 0, len(feature_ids))
    control_log = [math.log1p(value) for value in control_cpm]
    held_by_target = {entry["target"]: (index + 1, entry) for index, (entry, _, _, _) in enumerate(target_rows)}
    scores = []
    for frozen in prediction_fingerprints:
        target = frozen["target"]
        entry = held_by_target[target][1]
        if frozen["status"] != "predicted":
            scores.append({"target": target, "status": frozen["status"], "descriptorStatus": descriptor(entry, annotations)["status"]})
            continue
        with gzip.open(args.out / "folds" / target / "report.json.gz", "rt") as stream:
            prediction = json.load(stream)
        actual_cpm = row_cpm(training["matrix"], held_by_target[target][0], len(feature_ids))
        actual_response = [math.log1p(value) - baseline for value, baseline in zip(actual_cpm, control_log)]
        query = prediction["queries"][0]
        if query["status"] != "predicted":
            raise AssertionError("frozen report status changed")
        methods = {baseline["method"]: baseline["expression"] for baseline in prediction["baselines"]}
        methods["GO-term-kernel"] = query["expression"]
        methods["matched-supported-shuffle"] = query["shuffledExpression"]
        metrics = {
            name: metric([value - baseline for value, baseline in zip(expression, control_log)], actual_response)
            for name, expression in methods.items()
        }
        scores.append({"target": target, "status": "scored", "descriptorStatus": descriptor(entry, annotations)["status"], "metrics": metrics})

    scored = [record for record in scores if record["status"] == "scored" and record["descriptorStatus"] == "available"]
    aggregate = {}
    for method in ("GO-term-kernel", "noChange", "meanSingleResponse", "meanSupportedResponse", "matched-supported-shuffle"):
        values = [record["metrics"][method]["rmse"] for record in scored]
        aggregate[method] = {"targets": len(values), "meanRMSE": sum(values) / len(values) if values else None}
    go_rmse = aggregate["GO-term-kernel"]["meanRMSE"]
    comparators = [aggregate[name]["meanRMSE"] for name in ("noChange", "meanSingleResponse", "meanSupportedResponse")]
    outcome = {
        "status": "provisional-role-sensitivity-scored",
        "controlCandidates": list(CONTROL_LABELS),
        "controlMappingAuthoritative": False,
        "predictionQualification": False,
        "elapsedSeconds": time.monotonic() - started,
        "foldsRequested": len(fold_targets),
        "foldsScored": len(scored),
        "scores": scores,
        "aggregate": aggregate,
        "goLowerThanAllThreeComparators": go_rmse is not None and all(value is not None and go_rmse < value for value in comparators),
        "note": "Scores are conditional on the declared provisional pooled control and do not resolve experimental roles or establish biological transfer.",
    }
    write(args.out / "scores.json", outcome)
    write(args.out / "receipt.json", {
        "status": "passed-provisional-role-sensitivity",
        "binarySHA256": sha(args.binary),
        "driverSHA256": sha(Path(__file__)),
        "coverageSHA256": sha(args.coverage),
        "annotationsSHA256": sha(args.annotations),
        "cohortReportSHA256": sha(args.cohort_report),
        "cohortReceiptSHA256": sha(args.cohort_receipt),
        "trainingSHA256": sha(args.out / "training.json"),
        "predictionsFrozen": True,
        "foldsRequested": len(fold_targets),
        "foldsScored": len(scored),
        "controlMappingAuthoritative": False,
        "predictionQualification": False,
        "elapsedSeconds": time.monotonic() - started,
    })
    if list(args.out.rglob(".numivivo-target-kernel-*")):
        raise AssertionError("native staging directory remained")
    print(json.dumps(outcome["aggregate"], sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
