#!/usr/bin/env python3
"""Run the independently specified raw-count-to-neighbor Scanpy pipeline."""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.update(
    OPENBLAS_NUM_THREADS="1",
    OMP_NUM_THREADS="1",
    VECLIB_MAXIMUM_THREADS="1",
    NUMBA_NUM_THREADS="1",
)

# The driver is named scanpy.py for a discoverable command path. Remove its
# directory before importing the installed Scanpy package so the driver cannot
# shadow that package.
if sys.path and Path(sys.path[0] or ".").resolve() == Path(__file__).resolve().parent:
    sys.path.pop(0)

import anndata as ad
import numpy as np
import scanpy as sc
from scipy import sparse
from threadpoolctl import threadpool_limits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise SystemExit(f"output already contains files: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)

    sc.settings.n_jobs = 1
    with threadpool_limits(limits=1):
        started = time.perf_counter()
        data = ad.read_h5ad(args.source)
        if not sparse.issparse(data.X):
            raise SystemExit("source count matrix is dense")
        data.X = data.X.astype(np.float64)
        sc.pp.normalize_total(data, target_sum=10_000)
        sc.pp.log1p(data)
        sc.pp.highly_variable_genes(
            data, flavor="seurat", n_top_genes=2_000, n_bins=20
        )
        sc.pp.pca(
            data,
            n_comps=20,
            svd_solver="arpack",
            dtype="float64",
            random_state=7,
        )
        sc.pp.neighbors(
            data,
            n_neighbors=20,
            n_pcs=20,
            use_rep="X_pca",
            knn=True,
            method="umap",
            metric="euclidean",
            transformer="sklearn",
            random_state=7,
        )
        elapsed = time.perf_counter() - started

        distances = data.obsp["distances"].tocsr()
        connectivities = data.obsp["connectivities"].tocsr()
        distances.sort_indices()
        connectivities.sort_indices()
        sparse.save_npz(args.output / "distances.npz", distances)
        sparse.save_npz(args.output / "connectivities.npz", connectivities)
        np.save(args.output / "scores.npy", np.asarray(data.obsm["X_pca"], dtype=np.float64))

        if "native_sample" not in data.obs:
            raise SystemExit("source is missing native_sample identity column")
        selected = data.var_names[data.var["highly_variable"]].astype(str).tolist()
        identity = {
            "cellBarcodes": [str(value) for value in data.obs_names.tolist()],
            "sampleIDs": [str(value) for value in data.obs["native_sample"].tolist()],
            "selectedFeatureIDs": selected,
        }
        (args.output / "identity.json").write_text(json.dumps(identity, indent=2) + "\n")
        report = {
            "pipeline": "raw sparse H5AD -> normalize_total -> log1p -> Seurat HVG -> ARPACK PCA -> exact sklearn neighbors",
            "cells": int(data.n_obs),
            "sourceFeatures": int(data.n_vars),
            "selectedFeatures": len(selected),
            "components": 20,
            "normalizationTarget": 10_000,
            "hvg": {"flavor": "seurat", "nTopGenes": 2_000, "nBins": 20},
            "pca": {"svdSolver": "arpack", "dtype": "float64", "randomState": 7},
            "neighbors": {
                "nNeighbors": 20,
                "nPcs": 20,
                "method": "umap",
                "metric": "euclidean",
                "transformer": "sklearn",
                "randomState": 7,
            },
            "distanceEntries": int(distances.nnz),
            "connectivityEntries": int(connectivities.nnz),
            "pipelineSecondsInsideScript": elapsed,
            "sparseInput": True,
            "biologicalQualification": False,
        }
        (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
