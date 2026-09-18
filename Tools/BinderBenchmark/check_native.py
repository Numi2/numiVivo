#!/usr/bin/env python3
"""Test exact binder Swift/CLI source with the existing portable SHA/JSON facade.
No full Apple package, production artifact-store or biological qualification is implied.
"""
from pathlib import Path
import csv
import hashlib
import json
import math
import random
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]


def close(a: float, b: float, tolerance: float = 1e-10) -> None:
    if not math.isclose(a, b, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError((a, b))


def metric_reference(labels: list[float], scores: list[float], result: dict, k: int) -> None:
    positives = [s for y, s in zip(labels, scores) if y == 1]
    negatives = [s for y, s in zip(labels, scores) if y == 0]
    if positives and negatives:
        auc = sum(float(a > b) + 0.5 * (a == b) for a in positives for b in negatives)
        close(result["auroc"], auc / (len(positives) * len(negatives)))
    else:
        assert result.get("auroc") is None
    k = min(k, len(labels))
    cut = sorted(scores, reverse=True)[k - 1]
    above = [y for y, s in zip(labels, scores) if s > cut]
    tied = [y for y, s in zip(labels, scores) if s == cut]
    expected = sum(above) + (k - len(above)) * sum(tied) / len(tied)
    close(result["expectedHitsAtK"], expected)
    close(result["precisionAtK"], expected / k)
    if positives:
        ap = 0.0
        for threshold in set(scores):
            newly_positive = sum(y for y, s in zip(labels, scores) if s == threshold)
            selected = [y for y, s in zip(labels, scores) if s >= threshold]
            ap += newly_positive / len(positives) * sum(selected) / len(selected)
        close(result["averagePrecision"], ap)
    if "brier" in result:
        close(result["brier"], sum((y - s) ** 2 for y, s in zip(labels, scores)) / len(labels))


def cli_check(binary: Path, root: Path) -> None:
    rng = random.Random(183)
    names = ["ipsae_min_boltz2", "ipsae_min_ptxv2", "ipsae_min_ef2full"]
    source, config = root / "source.csv", root / "config.json"
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    with source.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["uuid", "target", "sequence", "adaptyv_binding"] + names)
        for i in range(90):
            x = [rng.random() for _ in names]
            prob = 1 / (1 + math.exp(-(-1 + 3 * x[0] - x[1] + x[2])))
            outcome = "binder" if rng.random() < prob else "non_binder"
            writer.writerow([f"synthetic-{i:03}", "TRAIN" if i < 60 else "TEST",
                             "ACD" + alphabet[i // 20] + alphabet[i % 20], outcome] + x)
    config.write_text(json.dumps({"schemaVersion": 1, "assay": "adaptyv", "targets": ["TRAIN", "TEST"], "features": names}))
    def run(*args: str, expected: int = 0) -> subprocess.CompletedProcess:
        result = subprocess.run([str(binary), *map(str, args)], capture_output=True, text=True)
        if result.returncode != expected:
            raise AssertionError((args, result.returncode, result.stderr))
        return result
    imported, result, plan = root / "input", root / "result", root / "plan.json"
    run("binder-import", source, config, imported)
    subprocess.run(["python3", str(ROOT / "Tools/BinderBenchmark/make_plan.py"), str(imported),
                    "--test-target", "TEST", "--output", str(plan)], check=True)
    run("binder-evaluate", imported, plan, result)
    run("binder-verify", result)
    run("binder-import", expected=64)
    run("binder-import", source, config, imported, expected=65)
    receipt = json.loads((result / "receipt.json").read_text())
    assert receipt["implementationSHA256"] == hashlib.sha256(binary.read_bytes()).hexdigest()
    for name, digest in receipt["files"].items():
        assert hashlib.sha256((result / name).read_bytes()).hexdigest() == digest
    assert (result / "source.csv").read_bytes() == source.read_bytes()
    report = json.loads((result / "report.json").read_text())
    rows = {r["id"]: r for r in json.loads((result / "imported.json").read_text())["dataset"]["records"]}
    train = [rows[i] for i in report["trainingIDs"]]
    model = report["model"]
    for j, name in enumerate(names):
        mean = sum(r["features"][name] for r in train) / len(train)
        scale = math.sqrt(sum((r["features"][name] - mean) ** 2 for r in train) / len(train))
        close(model["means"][j], mean); close(model["scales"][j], scale)
    def vector(row: dict) -> list[float]:
        return [(row["features"][name] - model["means"][j]) / model["scales"][j] for j, name in enumerate(names)]
    def predict(row: dict) -> float:
        z = model["intercept"] + sum(a * b for a, b in zip(vector(row), model["coefficients"]))
        return 1 / (1 + math.exp(-z))
    errors = [predict(r) - float(r["outcome"] == "binder") for r in train]
    gradient = [sum(errors) / len(train)]
    for j, w in enumerate(model["coefficients"]):
        gradient.append(sum(e * vector(r)[j] for e, r in zip(errors, train)) / len(train) + 0.1 * w)
    assert max(map(abs, gradient)) <= 1.01e-7
    prevalence = sum(r["outcome"] == "binder" for r in train) / len(train)
    for target in report["targets"]:
        for identity, p in zip(target["candidateIDs"], target["modelProbabilities"]):
            assert rows[identity]["target"] == "TEST"
            close(p, predict(rows[identity]))
        metric_reference(target["labels"], target["baselineScores"], target["baseline"], 10)
        metric_reference(target["labels"], target["modelProbabilities"], target["learned"], 10)
        metric_reference(target["labels"], [prevalence] * len(target["labels"]), target["trainingPrevalence"], 10)
    (result / "report.json").write_text("{}")
    run("binder-verify", result, expected=65)
    print("PASS: synthetic CLI import/evaluate/replay, SHA-256, independent Python fit/metric checks, tamper/overwrite rejection")


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="numivivo-binder-") as tmp:
        root = Path(tmp)
        (root / "Sources/NumiVivoKit").mkdir(parents=True)
        (root / "Tests/BinderTests").mkdir(parents=True)
        (root / "Sources/BinderCLI").mkdir(parents=True)
        shutil.copy2(ROOT / "Tools/Posterior/PortableSupport.swift", root / "Sources/NumiVivoKit")
        shutil.copy2(ROOT / "Sources/NumiVivoCLI/VivoBinderCLICommands.swift", root / "Sources/BinderCLI")
        (root / "Sources/BinderCLI/main.swift").write_text("import Foundation\nexit(VivoBinderCLICommands().run(arguments: Array(CommandLine.arguments.dropFirst())))\n")
        for source in sorted((ROOT / "Sources/NumiVivoKit/Binder").glob("*.swift")):
            shutil.copy2(source, root / "Sources/NumiVivoKit" / source.name)
        for test in sorted((ROOT / "Tests/NumiVivoIntegrationTests").glob("VivoBinder*Tests.swift")):
            shutil.copy2(test, root / "Tests/BinderTests" / test.name)
        (root / "Package.swift").write_text('''// swift-tools-version: 6.0
import PackageDescription
let package = Package(name: "BinderCheck", targets: [
    .target(name: "NumiVivoKit", linkerSettings: [.linkedLibrary("crypto", .when(platforms: [.linux]))]),
    .executableTarget(name: "BinderCLI", dependencies: ["NumiVivoKit"]),
    .testTarget(name: "BinderTests", dependencies: ["NumiVivoKit"])
], swiftLanguageModes: [.v6])
''')
        subprocess.run(["swift", "test", "--package-path", str(root), "--jobs", "2"], check=True)
        cli_check(root / ".build/debug/BinderCLI", root)

if __name__ == "__main__":
    main()
