#!/usr/bin/env python3
"""Execute a fixed six-fold score-only comparison on a pinned published table.

No tuning, structural/MD features or biological qualification. All protocols are
written before importing/scoring. Laboratory endpoints remain separate. Outputs
are immutable and include exact inputs, executable and independently checked metrics.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from check_native import close, metric_reference
from fetch_source import BLOB, FEATURES, REVISION, TARGETS, URL


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def dump(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n")


def verify_import(raw: bytes, imported: dict, assay: str) -> None:
    rows = list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
    selected = {r["uuid"]: r for r in rows if r["target"] in TARGETS}
    records = imported["dataset"]["records"]
    if len(records) != len(selected) or len(selected) != sum(r["target"] in TARGETS for r in rows):
        raise AssertionError("source/import candidate identity mismatch")
    expected_outcomes = {"binder": "binder", "non_binder": "nonBinder",
                         "not_tested": "notTested", "expression_failure": "expressionFailure"}
    if imported["dataset"]["sourceSHA256"] != sha(raw):
        raise AssertionError("source hash mismatch")
    for record in records:
        source = selected[record["id"]]
        if record["sourceFields"] != source or record["target"] != source["target"]:
            raise AssertionError("source fields changed")
        outcome = source[f"{assay}_binding"]
        if record["outcome"] != expected_outcomes.get(outcome, "inconclusive"):
            raise AssertionError("assay outcome changed")
        if record["rawOutcome"] != (outcome or "<missing>"):
            raise AssertionError("raw assay outcome changed")
        if record["leakageGroup"] != "exact-sequence:" + source["sequence"]:
            raise AssertionError("sequence group changed")
        expected = {name: float(source[name]) for name in FEATURES if source[name] != ""}
        if record["features"] != expected:
            raise AssertionError("feature values changed")


def verify_report(imported: dict, report: dict) -> None:
    records = {r["id"]: r for r in imported["dataset"]["records"]}
    plan, model = report["plan"], report["model"]
    test_targets, train_targets = set(plan["testTargets"]), set(plan["trainingTargets"])
    test_groups = {r["leakageGroup"] for r in records.values() if r["target"] in test_targets}
    names = plan["modelFeatures"]
    required = set(names + [plan["baselineFeature"]])
    def eligible(row: dict) -> bool:
        return row["outcome"] in ("binder", "nonBinder") and required <= row["features"].keys()
    expected_train = sorted(r["id"] for r in records.values()
                            if r["target"] in train_targets and r["leakageGroup"] not in test_groups and eligible(r))
    if report["trainingIDs"] != expected_train or model["featureNames"] != names:
        raise AssertionError("training population or feature identity mismatch")
    train = [records[i] for i in expected_train]
    for j, name in enumerate(names):
        mean = sum(r["features"][name] for r in train) / len(train)
        variance = sum((r["features"][name] - mean) ** 2 for r in train) / len(train)
        close(model["means"][j], mean)
        close(model["scales"][j], math.sqrt(variance) if variance > 1e-20 else 1)
    def vector(row: dict) -> list[float]:
        return [(row["features"][name] - model["means"][j]) / model["scales"][j] for j, name in enumerate(names)]
    def predict(row: dict) -> float:
        z = model["intercept"] + sum(a * b for a, b in zip(vector(row), model["coefficients"]))
        return 1 / (1 + math.exp(-z)) if z >= 0 else math.exp(z) / (1 + math.exp(z))
    errors = [predict(r) - float(r["outcome"] == "binder") for r in train]
    gradient = [sum(errors) / len(train)]
    for j, weight in enumerate(model["coefficients"]):
        gradient.append(sum(e * vector(r)[j] for e, r in zip(errors, train)) / len(train) + plan["ridgePenalty"] * weight)
    if max(map(abs, gradient)) > 1.01e-7:
        raise AssertionError("independent fit-gradient check failed")
    prevalence = sum(r["outcome"] == "binder" for r in train) / len(train)
    if {t["target"] for t in report["targets"]} != test_targets:
        raise AssertionError("test-target identity mismatch")
    for target in report["targets"]:
        expected_ids = sorted(r["id"] for r in records.values() if r["target"] == target["target"] and eligible(r))
        if expected_ids != target["candidateIDs"]:
            raise AssertionError("test population changed")
        labels = [float(records[i]["outcome"] == "binder") for i in expected_ids]
        if labels != target["labels"]:
            raise AssertionError("test outcomes changed")
        if len(target["modelProbabilities"]) != len(labels) or len(target["baselineScores"]) != len(labels):
            raise AssertionError("prediction dimensions changed")
        for identity, p, score in zip(expected_ids, target["modelProbabilities"], target["baselineScores"]):
            close(p, predict(records[identity]))
            close(score, records[identity]["features"][plan["baselineFeature"]])
        metric_reference(labels, target["baselineScores"], target["baseline"], plan["topK"])
        metric_reference(labels, target["modelProbabilities"], target["learned"], plan["topK"])
        metric_reference(labels, [prevalence] * len(labels), target["trainingPrevalence"], plan["topK"])


def execute(source: Path, binary: Path, output: Path) -> dict:
    source, binary, output = source.absolute(), binary.absolute(), output.absolute()
    if output.exists() or not output.parent.is_dir():
        raise ValueError("output must be new and its parent must exist")
    if binary.is_symlink() or not binary.is_file():
        raise ValueError("binary must be a regular file")
    source_file = source / "source.csv"
    if source_file.is_symlink() or not source_file.is_file() or source_file.stat().st_size > 64 * 1024 * 1024:
        raise ValueError("source must be a bounded regular file")
    raw = source_file.read_bytes()
    if len(raw) > 64 * 1024 * 1024:
        raise ValueError("source exceeds 64 MiB")
    blob = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
    if blob != BLOB:
        raise ValueError("source does not match the frozen published Git blob")
    protocol = {"schemaVersion": 1, "sourceRevision": REVISION, "sourceGitBlobSHA1": BLOB,
                "sourceSHA256": sha(raw), "sourceURL": URL, "targets": TARGETS,
                "assays": ["adaptyv", "twist"], "modelFeatures": FEATURES,
                "baselineFeature": "ipsae_min_boltz2", "topK": 10, "ridgePenalty": 0.1,
                "holdout": "one complete target at a time; exact-sequence overlap purged",
                "evidenceStatus": "fixed-retrospective-development; not prospective validation"}
    with tempfile.TemporaryDirectory(prefix=".binder-campaign-", dir=output.parent) as tmp:
        staging = Path(tmp) / "campaign"
        staging.mkdir()
        (staging / "source.csv").write_bytes(raw)
        if (source / "SOURCE.json").is_file():
            shutil.copy2(source / "SOURCE.json", staging / "SOURCE.json")
        shutil.copy2(binary, staging / "native-cli")
        implementation = sha((staging / "native-cli").read_bytes())
        plans = staging / "plans"; plans.mkdir()
        # All six plans are constructed before any import or outcome analysis.
        for assay in protocol["assays"]:
            dump(staging / f"import-{assay}.json", {"schemaVersion": 1, "assay": assay,
                                                   "targets": TARGETS, "features": FEATURES})
            for target in TARGETS:
                plan = {"schemaVersion": 1, "sourceSHA256": sha(raw),
                        "trainingTargets": [t for t in TARGETS if t != target], "testTargets": [target],
                        "baselineFeature": protocol["baselineFeature"], "modelFeatures": FEATURES,
                        "topK": 10, "ridgePenalty": 0.1}
                dump(plans / f"{assay}-{target}.json", plan)
        protocol["planSHA256"] = {p.name: sha(p.read_bytes()) for p in sorted(plans.glob("*.json"))}
        dump(staging / "protocol.json", protocol)
        frozen = sha((staging / "protocol.json").read_bytes())
        logs = staging / "logs"; logs.mkdir()
        def run(label: str, *args: object) -> None:
            completed = subprocess.run([str(staging / "native-cli"), *map(str, args)], capture_output=True, text=True, timeout=1800)
            (logs / f"{label}.stdout").write_text(completed.stdout)
            (logs / f"{label}.stderr").write_text(completed.stderr)
            if completed.returncode:
                raise RuntimeError(f"{label} failed ({completed.returncode}): {completed.stderr}")
        summaries = []
        for assay in protocol["assays"]:
            imported_dir = staging / f"input-{assay}"
            run(f"import-{assay}", "binder-import", staging / "source.csv", staging / f"import-{assay}.json", imported_dir)
            imported = json.loads((imported_dir / "imported.json").read_text())
            verify_import(raw, imported, assay)
            for target in TARGETS:
                label = f"{assay}-{target}"
                result = staging / f"result-{label}"
                run(label, "binder-evaluate", imported_dir, plans / f"{label}.json", result)
                run(f"verify-{label}", "binder-verify", result)
                report = json.loads((result / "report.json").read_text())
                verify_report(imported, report)
                item = report["targets"][0]
                summaries.append({"assay": assay, "target": target, "trainingCount": len(report["trainingIDs"]),
                                  "excludedCount": len(report["excluded"]), "baseline": item["baseline"],
                                  "scoreEnsemble": item["learned"], "trainingPrevalence": item["trainingPrevalence"]})
        if sha((staging / "protocol.json").read_bytes()) != frozen:
            raise AssertionError("protocol changed")
        summary = {"schemaVersion": 1, "status": "completed", "completedAt": datetime.now(timezone.utc).isoformat(),
                   "implementationSHA256": implementation, "protocolSHA256": frozen,
                   "sourceSHA256": sha(raw), "sourceRows": len(list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))),
                   "folds": summaries, "independentNumericalVerification": "passed",
                   "limitations": ["Score-only comparison; no Numi physical features or demonstrated physical uplift.",
                                   "Retrospective development; related designs/targets are not independent replications.",
                                   "Exact-sequence groups do not establish homology independence.",
                                   "Assay endpoints remain separate; missing outcomes are not failures.",
                                   "The supplied executable may be the portable subset, not the full Apple product."]}
        dump(staging / "summary.json", summary)
        files = {str(p.relative_to(staging)): sha(p.read_bytes()) for p in sorted(staging.rglob("*")) if p.is_file()}
        dump(staging / "manifest.json", {"schemaVersion": 1, "files": files})
        if output.exists():
            raise FileExistsError(output)
        os.rename(staging, output)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="directory from fetch_source.py")
    parser.add_argument("binary", type=Path, help="numivivo or portable BinderCLI executable")
    parser.add_argument("output", type=Path, help="new campaign directory")
    args = parser.parse_args()
    print(json.dumps(execute(args.source, args.binary, args.output), sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
