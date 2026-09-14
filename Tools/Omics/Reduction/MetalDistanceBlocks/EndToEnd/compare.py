#!/usr/bin/env python3
"""Compare the final native graphs with the independent Scanpy graph."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse


def read_json(path: Path):
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def native_graph(root: Path) -> tuple[dict, np.ndarray, np.ndarray, sparse.csr_matrix]:
    graph = read_json(root / "graph.json")
    cells = len(graph["cells"])
    indices = np.asarray(graph["neighborIndices"], dtype=np.int64)
    distances = np.asarray(graph["neighborDistances"], dtype=np.float64)
    if indices.size != distances.size or indices.size % cells != 0:
        raise SystemExit(f"invalid native neighbor shape in {root}")
    width = indices.size // cells
    indices = indices.reshape(cells, width)
    distances = distances.reshape(cells, width)
    matrix = sparse.csr_matrix(
        (
            np.asarray(graph["weights"], dtype=np.float64),
            np.asarray(graph["columnIndices"], dtype=np.int64),
            np.asarray(graph["rowOffsets"], dtype=np.int64),
        ),
        shape=(cells, cells),
    )
    matrix.sort_indices()
    if matrix.nnz != len(graph["weights"]):
        raise SystemExit(f"invalid native graph offsets in {root}")
    return graph, indices, distances, matrix


def cell_identity(graph: dict) -> tuple[list[str], list[str]]:
    return (
        [str(item["barcode"]) for item in graph["cells"]],
        [str(item["sampleID"]) for item in graph["cells"]],
    )


def selected_features(pca_root: Path) -> list[str]:
    model = read_json(pca_root / "model.json")
    features = model["features"]
    return [str(features[index]["featureID"]) for index in model["selectedFeatureIndices"]]


def compare_neighbors(
    native_indices: np.ndarray,
    native_distances: np.ndarray,
    reference: sparse.csr_matrix,
) -> dict:
    if native_indices.shape[0] != reference.shape[0]:
        raise SystemExit("native and Scanpy cell counts differ")
    missing = 0
    extra = 0
    reference_entries = 0
    max_distance_error = 0.0
    rows_with_duplicate_native = 0
    for row in range(reference.shape[0]):
        start, end = reference.indptr[row : row + 2]
        reference_by_cell = {
            int(index): float(value)
            for index, value in zip(reference.indices[start:end], reference.data[start:end])
            if int(index) != row
        }
        candidate_by_cell = {}
        for index, value in zip(native_indices[row], native_distances[row]):
            index = int(index)
            if index == row:
                continue
            if index in candidate_by_cell:
                rows_with_duplicate_native += 1
            candidate_by_cell[index] = float(value)
        reference_cells = set(reference_by_cell)
        candidate_cells = set(candidate_by_cell)
        missing += len(reference_cells - candidate_cells)
        extra += len(candidate_cells - reference_cells)
        reference_entries += len(reference_cells)
        for index in reference_cells & candidate_cells:
            max_distance_error = max(
                max_distance_error,
                abs(reference_by_cell[index] - candidate_by_cell[index]),
            )
    return {
        "referenceNonSelfNeighborEntries": reference_entries,
        "missingNeighbors": missing,
        "extraNeighbors": extra,
        "neighborSetMatch": missing == 0 and extra == 0 and rows_with_duplicate_native == 0,
        "rowsWithDuplicateNativeNeighbors": rows_with_duplicate_native,
        "maximumCommonNeighborDistanceDifference": max_distance_error,
    }


def compare_edges(native: sparse.csr_matrix, reference: sparse.csr_matrix) -> dict:
    if native.shape != reference.shape:
        raise SystemExit("native and Scanpy graph shapes differ")
    native.sort_indices()
    reference.sort_indices()
    n = native.shape[0]
    native_keys = np.repeat(np.arange(n, dtype=np.int64), np.diff(native.indptr)) * n + native.indices
    reference_keys = (
        np.repeat(np.arange(n, dtype=np.int64), np.diff(reference.indptr)) * n + reference.indices
    )
    common, native_at, reference_at = np.intersect1d(
        native_keys, reference_keys, assume_unique=True, return_indices=True
    )
    if len(common):
        weight_error = float(
            np.max(np.abs(native.data[native_at] - reference.data[reference_at]))
        )
    else:
        weight_error = 0.0
    return {
        "nativeEdges": int(native.nnz),
        "scanpyEdges": int(reference.nnz),
        "commonEdges": int(len(common)),
        "missingEdges": int(len(reference_keys) - len(common)),
        "extraEdges": int(len(native_keys) - len(common)),
        "exactEdgeCoordinates": len(common) == len(native_keys) == len(reference_keys),
        "maximumCommonFuzzyWeightDifference": weight_error,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--scanpy", type=Path, required=True)
    parser.add_argument("--native-cpu", type=Path, required=True)
    parser.add_argument("--native-metal", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"comparison output already exists: {args.output}")

    scanpy_identity = read_json(args.scanpy / "identity.json")
    reference_distances = sparse.load_npz(args.scanpy / "distances.npz").tocsr()
    reference_connectivities = sparse.load_npz(args.scanpy / "connectivities.npz").tocsr()
    reference_distances.sort_indices()
    reference_connectivities.sort_indices()

    comparisons = {}
    selected = {}
    identities = {}
    for label, root in (("cpu", args.native_cpu), ("metal", args.native_metal)):
        graph, indices, distances, native_connectivities = native_graph(root / "graph")
        barcodes, sample_ids = cell_identity(graph)
        native_selected = selected_features(root / "pca")
        identities[label] = {
            "cellBarcodesMatch": barcodes == scanpy_identity["cellBarcodes"],
            "sampleIDsMatch": sample_ids == scanpy_identity["sampleIDs"],
            "cells": len(barcodes),
        }
        selected[label] = {
            "nativeCount": len(native_selected),
            "scanpyCount": len(scanpy_identity["selectedFeatureIDs"]),
            "orderedMatch": native_selected == scanpy_identity["selectedFeatureIDs"],
            "setMatch": set(native_selected) == set(scanpy_identity["selectedFeatureIDs"]),
        }
        comparisons[label] = {
            **compare_neighbors(indices, distances, reference_distances),
            **compare_edges(native_connectivities, reference_connectivities),
        }

    gate = all(
        item["cellBarcodesMatch"] and item["sampleIDsMatch"] for item in identities.values()
    ) and all(item["setMatch"] for item in selected.values()) and all(
        item["neighborSetMatch"] and item["exactEdgeCoordinates"] for item in comparisons.values()
    )
    result = {
        "status": "PASS" if gate else "FAIL",
        "scope": "final repetition raw-count-to-graph numerical comparison; timing and memory are separate",
        "sourceSHA256": sha256(args.source),
        "cells": int(reference_distances.shape[0]),
        "identities": identities,
        "selectedFeatures": selected,
        "comparisons": comparisons,
        "biologicalQualification": False,
        "matchedTimingClaim": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
