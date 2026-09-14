#!/usr/bin/env python3
"""Run the opt-in native negative-binomial route on frozen IFNB folds.

The source fold plans and H5AD files are immutable inputs.  This driver writes
an explicit copy of each training plan with ``responseModel`` set to
``negativeBinomial`` and invokes only the native fit/predict/verifier commands.
Held-out treated rows are never passed to the native process.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


parser = argparse.ArgumentParser()
parser.add_argument("--inputs", type=Path, required=True)
parser.add_argument("--binary", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)

freeze = json.loads((args.inputs / "input-freeze.json").read_text())
for relative, digest in freeze["files"].items():
    if sha256(args.inputs / relative) != digest:
        raise ValueError(f"input freeze mismatch: {relative}")

execution = {
    "binary": str(args.binary.resolve()),
    "binarySHA256": sha256(args.binary),
    "inputFreezeSHA256": sha256(args.inputs / "input-freeze.json"),
    "runnerSHA256": sha256(Path(__file__)),
    "platform": platform.platform(),
    "hdf5Library": os.environ.get("NUMIVIVO_HDF5_LIBRARY"),
    "startedUnix": time.time(),
}
write(args.out / "execution.json", execution)

commands: list[dict[str, Any]] = []


def invoke(arguments: list[Path | str], log: Path) -> bool:
    started = time.time()
    completed = subprocess.run(
        [str(args.binary.resolve()), *(str(value) for value in arguments)],
        capture_output=True,
        text=True,
        env=os.environ.copy(),
    )
    log.write_text(completed.stdout + completed.stderr)
    commands.append(
        {
            "arguments": [str(value) for value in arguments],
            "exitCode": completed.returncode,
            "elapsedSeconds": time.time() - started,
            "log": str(log.relative_to(args.out)),
        }
    )
    write(args.out / "commands.json", commands)
    return completed.returncode == 0


folds = json.loads((args.inputs / "folds.json").read_text())
failures: list[str] = []
completed_folds = 0
for fold_index, fold in enumerate(folds):
    source = args.inputs / fold["id"]
    destination = args.out / fold["id"]
    destination.mkdir()
    original_plan = json.loads((source / "training.json").read_text())
    negative_binomial_plan = dict(original_plan)
    negative_binomial_plan["responseModel"] = "negativeBinomial"
    write(destination / "training-nb.json", negative_binomial_plan)

    stages = [
        (
            [
                "singlecell-perturbation-fit",
                source / "training.h5ad",
                "--plan",
                destination / "training-nb.json",
                "--output",
                destination / "model",
            ],
            destination / "fit.log",
            "fit",
        ),
        (
            ["singlecell-perturbation-verify", destination / "model"],
            destination / "verify-model.log",
            "verify-model",
        ),
        (
            [
                "singlecell-perturbation-predict",
                source / "query.h5ad",
                "--plan",
                source / "query.json",
                "--reference",
                destination / "model",
                "--output",
                destination / "prediction",
            ],
            destination / "predict.log",
            "predict",
        ),
        (
            ["singlecell-perturbation-prediction-verify", destination / "prediction"],
            destination / "verify-prediction.log",
            "verify-prediction",
        ),
    ]
    fold_ok = True
    for arguments, log, stage in stages:
        if not invoke(arguments, log):
            failures.append(f"{fold['id']}:{stage}")
            fold_ok = False
            break
    if fold_ok and fold_index == 0:
        repeat = destination / "repeat"
        repeat_ok = invoke(
            [
                "singlecell-perturbation-predict",
                source / "query.h5ad",
                "--plan",
                source / "query.json",
                "--reference",
                destination / "model",
                "--output",
                repeat,
            ],
            destination / "repeat.log",
        )
        if not repeat_ok:
            failures.append(f"{fold['id']}:repeat-predict")
            fold_ok = False
        elif (repeat / "report.json").read_bytes() != (destination / "prediction/report.json").read_bytes():
            failures.append(f"{fold['id']}:repeat-bytes")
            fold_ok = False
    if fold_ok:
        completed_folds += 1
    print(json.dumps({"completed": completed_folds, "fold": fold["id"], "ok": fold_ok}), flush=True)

files = {
    str(path.relative_to(args.out)): sha256(path)
    for path in sorted(args.out.rglob("*"))
    if path.is_file()
}
write(
    args.out / "prediction-freeze.json",
    {
        "schemaVersion": 1,
        "files": files,
        "commands": len(commands),
        "completedFolds": completed_folds,
        "totalFolds": len(folds),
        "scoringStarted": False,
        "failures": failures,
        "completedUnix": time.time(),
    },
)
status = "passed-complete-native-nb2" if not failures and completed_folds == len(folds) else "failed-native-nb2"
manifest = {
    "status": status,
    "inputFreezeSHA256": sha256(args.inputs / "input-freeze.json"),
    "runnerSHA256": sha256(Path(__file__)),
    "binarySHA256": sha256(args.binary),
    "folds": len(folds),
    "completedFolds": completed_folds,
    "commands": len(commands),
    "failures": failures,
    "predictionFreezeSHA256": sha256(args.out / "prediction-freeze.json"),
}
write(args.out / "manifest.json", manifest)
print(json.dumps(manifest, indent=2, sort_keys=True), flush=True)
if failures:
    raise SystemExit(1)
