#!/usr/bin/env python3
"""Run a fresh, matched native and Scanpy raw-count-to-graph comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


ORDER = (
    ("cpu", "scanpy", "metal"),
    ("scanpy", "metal", "cpu"),
    ("metal", "cpu", "scanpy"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def output_digest(root: Path, excluded_names: set[str] | None = None) -> str:
    excluded_names = excluded_names or set()
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name in excluded_names:
            continue
        relative = path.relative_to(root).as_posix().encode()
        digest.update(relative)
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256(path)))
    return digest.hexdigest()


def timed_command(
    command: list[str], log_path: Path, environment: dict[str, str]
) -> tuple[float, int]:
    started = time.monotonic()
    with log_path.open("w") as log:
        completed = subprocess.run(
            ["/usr/bin/time", "-l", *command],
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
            text=True,
            check=False,
        )
    elapsed = time.monotonic() - started
    text = log_path.read_text()
    matches = re.findall(r"(\d+)\s+maximum resident set size", text)
    if not matches:
        raise RuntimeError(f"time output has no maximum resident set size: {log_path}")
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}; see {log_path}"
        )
    return elapsed, max(int(value) for value in matches)


def checked_command(
    command: list[str], log_path: Path, environment: dict[str, str]
) -> None:
    with log_path.open("a") as log:
        completed = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=environment,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"verification failed ({completed.returncode}): {' '.join(command)}; see {log_path}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--pca-plan", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--scanpy-script", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cohort", default=None, help="Human-readable cohort identifier for the protocol")
    parser.add_argument("--cells-expected", type=int, default=None, help="Expected cell count recorded in the protocol")
    parser.add_argument("--repetitions", type=int, default=3)
    args = parser.parse_args()
    if args.cells_expected is not None and args.cells_expected <= 0:
        raise SystemExit("--cells-expected must be positive")
    for path in (args.source, args.pca_plan, args.binary, args.scanpy_script):
        if not path.is_file():
            raise SystemExit(f"missing input: {path}")
    if args.repetitions != 3:
        raise SystemExit("the predeclared protocol requires exactly three repetitions")
    if args.output.exists():
        raise SystemExit(f"output already exists; use a new path to retain prior evidence: {args.output}")
    args.output.mkdir(parents=True)

    graph_base = {
        "schemaVersion": 1,
        "inputKind": "fitted",
        "storage": "json",
        "neighbors": {
            "neighbors": 20,
            "maximumDistancePairs": 500_000_000,
            "representation": "pca",
        },
        "execution": {
            "workers": 1,
            "queryBlockRows": 128,
            "candidateBlockRows": 8_192,
        },
    }
    graph_cpu_plan = json.loads(json.dumps(graph_base))
    graph_metal_plan = json.loads(json.dumps(graph_base))
    graph_metal_plan["execution"]["backend"] = "metalFP32"
    shutil.copy2(args.pca_plan, args.output / "pca-plan.json")
    write_json(args.output / "graph-cpu-plan.json", graph_cpu_plan)
    write_json(args.output / "graph-metal-plan.json", graph_metal_plan)

    native_environment = os.environ.copy()
    native_environment["NUMIVIVO_HDF5_LIBRARY"] = (
        "/Users/n/numivivo-integration-reference-py/lib/python3.13/"
        "site-packages/h5py/.dylibs/libhdf5.320.0.0.dylib"
    )
    scanpy_environment = os.environ.copy()
    scanpy_environment.update(
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        NUMBA_NUM_THREADS="1",
        PYTHONHASHSEED="0",
    )
    protocol = {
        "cohort": args.cohort or args.source.stem,
        "cellsExpected": args.cells_expected,
        "repetitions": 3,
        "order": [list(item) for item in ORDER],
        "nativePipeline": "singlecell-h5ad-pca then singlecell-pca-neighbors",
        "nativeVerification": "pca and graph verify commands run after each timed pipeline",
        "scanpyPipeline": "raw sparse H5AD -> normalize_total(10000) -> log1p -> Seurat HVG(2000,20 bins) -> ARPACK PCA(20,seed 7) -> exact sklearn neighbors(20)",
        "timingBoundary": "sum of fresh command wall times, including startup, source reads and output writes; verification is outside the timing interval",
        "memoryBoundary": "maximum resident set size reported by /usr/bin/time -l across timed commands",
        "nativeGraphPlan": graph_base,
        "sourceSHA256": sha256(args.source),
        "pcaPlanSHA256": sha256(args.pca_plan),
        "binarySHA256": sha256(args.binary),
        "scanpyDriverSHA256": sha256(args.scanpy_script),
        "biologicalQualification": False,
    }
    write_json(args.output / "protocol.json", protocol)

    results: list[dict[str, object]] = []
    for repetition in range(args.repetitions):
        for backend in ORDER[repetition]:
            name = f"rep-{repetition + 1}-{backend}"
            root = args.output / name
            root.mkdir()
            if backend == "scanpy":
                elapsed, rss = timed_command(
                    [sys.executable, str(args.scanpy_script), str(args.source), str(root)],
                    args.output / f"{name}.log",
                    scanpy_environment,
                )
                result = {
                    "repetition": repetition + 1,
                    "backend": backend,
                    "pipelineWallSeconds": elapsed,
                    "peakRSSBytes": rss,
                    "verified": True,
                    "outputSHA256": output_digest(root),
                    "deterministicOutputSHA256": output_digest(root, {"report.json"}),
                }
            else:
                pca_root = root / "pca"
                graph_root = root / "graph"
                pca_plan_path = args.output / "pca-plan.json"
                graph_plan_path = args.output / f"graph-{backend}-plan.json"
                pca_seconds, pca_rss = timed_command(
                    [
                        str(args.binary),
                        "singlecell-h5ad-pca",
                        str(args.source),
                        "--plan",
                        str(pca_plan_path),
                        "--output",
                        str(pca_root),
                    ],
                    args.output / f"{name}-pca.log",
                    native_environment,
                )
                checked_command(
                    [str(args.binary), "singlecell-h5ad-pca-verify", str(pca_root)],
                    args.output / f"{name}-pca.log",
                    native_environment,
                )
                graph_seconds, graph_rss = timed_command(
                    [
                        str(args.binary),
                        "singlecell-pca-neighbors",
                        str(pca_root),
                        "--plan",
                        str(graph_plan_path),
                        "--output",
                        str(graph_root),
                    ],
                    args.output / f"{name}-graph.log",
                    native_environment,
                )
                checked_command(
                    [
                        str(args.binary),
                        "singlecell-pca-neighbors-verify",
                        str(graph_root),
                    ],
                    args.output / f"{name}-graph.log",
                    native_environment,
                )
                result = {
                    "repetition": repetition + 1,
                    "backend": backend,
                    "pcaWallSeconds": pca_seconds,
                    "graphWallSeconds": graph_seconds,
                    "pipelineWallSeconds": pca_seconds + graph_seconds,
                    "peakRSSBytes": max(pca_rss, graph_rss),
                    "verified": True,
                    "outputSHA256": output_digest(root),
                    "deterministicOutputSHA256": output_digest(root),
                }
            results.append(result)
            write_json(args.output / "timings.json", {"results": results})

    consistency = {}
    for backend in ("cpu", "metal", "scanpy"):
        hashes = [str(item["outputSHA256"]) for item in results if item["backend"] == backend]
        deterministic_hashes = [
            str(item["deterministicOutputSHA256"])
            for item in results
            if item["backend"] == backend
        ]
        consistency[backend] = {
            "repetitions": len(hashes),
            "allOutputHashesMatch": len(set(hashes)) == 1,
            "allDeterministicArtifactHashesMatch": len(set(deterministic_hashes)) == 1,
        }
    write_json(args.output / "timings.json", {"results": results, "repetitionConsistency": consistency})

    comparison_path = args.output / "comparison.json"
    comparison_log = args.output / "comparison.log"
    compare_command = [
        sys.executable,
        str(Path(__file__).with_name("compare.py")),
        "--source",
        str(args.source),
        "--scanpy",
        str(args.output / "rep-3-scanpy"),
        "--native-cpu",
        str(args.output / "rep-3-cpu"),
        "--native-metal",
        str(args.output / "rep-3-metal"),
        "--output",
        str(comparison_path),
    ]
    with comparison_log.open("w") as log:
        completed = subprocess.run(
            compare_command,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=scanpy_environment,
            text=True,
            check=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"numerical comparison failed; see {comparison_log}")
    comparison = json.loads(comparison_path.read_text())
    if args.cells_expected is not None and comparison.get("cells") != args.cells_expected:
        raise RuntimeError(
            f"source cell count {comparison.get('cells')} differs from --cells-expected {args.cells_expected}"
        )
    summary = {
        "status": "PASS",
        "protocolSHA256": sha256(args.output / "protocol.json"),
        "timing": results,
        "repetitionConsistency": consistency,
        "comparison": str(comparison_path),
        "biologicalQualification": False,
    }
    write_json(args.output / "summary.json", summary)
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
