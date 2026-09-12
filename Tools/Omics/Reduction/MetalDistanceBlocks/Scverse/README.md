# Complete Kang graph construction versus exact Scanpy

On the same frozen **24,673-cell, 20-PC Kang representation**, native Metal FP32
constructs the neighbor/fuzzy graph in median **2.353 s**, versus **4.622 s** for
Scanpy's exact sklearn-backed path: **49.1% lower elapsed time (1.96x)** in this
measured stage. Native CPU FP64 takes 5.372 s. All 468,787 non-self neighbors and
706,016 connectivity coordinates match Scanpy in both native backends.

## Matching scope and numerical evidence

The [official Scanpy neighbor API](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.neighbors.html)
supports the explicit sklearn transformer. The installed source confirms this
selects brute-force KNeighborsTransformer. The experiment fixes `n_neighbors=20`,
`use_rep="X_pca"`, `n_pcs=20`, Euclidean metric, UMAP fuzzy connectivity, seed 7,
and one CPU worker/thread. It does not use Scanpy's approximate default route.
Native execution also uses one CPU worker, query tiles 128 and candidate tiles
8192. Metal uses the published native owner; no substitute GPU algorithm is used.

All three repetitions per backend reproduce their graph outputs exactly.
Maximum common neighbor-distance difference versus Scanpy is 1.74e-13 for CPU
and 1.48e-6 for Metal. Maximum fuzzy-weight difference is 4.52e-6 CPU and 1.53e-5
Metal. Matching topology is not bitwise weight equality. This adds an ecosystem
numerical comparison; source-type preservation remains the separate
[Kang qualification](../Kang/README.md), not new biological generalization.

## Timing boundaries and every run

| Backend | Graph-stage seconds | Median | Total process seconds |
| --- | --- | ---: | --- |
| Native CPU FP64 | 5.372, 5.451, 5.365 | 5.372 | 6.328, 5.967, 5.881 |
| Native Metal FP32 | 2.375, 2.352, 2.353 | 2.353 | 2.892, 2.875, 2.872 |
| Scanpy exact | 18.632, 4.594, 4.622 | 4.622 | 75.810, 6.175, 6.210 |

The predeclared order is CPU/Scanpy/Metal, Scanpy/Metal/CPU, Metal/CPU/Scanpy.
Every first run is retained. The first Scanpy process built a Matplotlib font
cache before the graph timer; its first graph call was also slower. No separate
profile establishes how much graph time each initialization step caused.

The graph timer includes reading/validating the original score records and
constructing neighbors plus fuzzy connectivity. Native metadata/identity decoding
and Python imports precede the timer; serialization follows it. Scanpy constructs
its AnnData container inside the timer and uses a sparse empty count axis plus
resident PCA scores. Native runs the resident-result graph route while reading
scores through its bounded tiles. Neither creates a dense cells-by-genes matrix.
These are concrete implementation scopes, not identical internal memory layouts.

Process peak RSS ranges 323.0–326.2 MB CPU, 331.6–337.8 MB Metal and 511.4–575.8 MB
Scanpy, including runtime/import/setup costs. Results are from exercised caches
on one shared physical M4 Pro desktop, with three fixed-order repetitions. There
is no confidence interval, cold-cache claim, approximate-Scanpy comparison or
universal GPU speedup. **Do not compare this stage timing with full NumiVivo CLI
publication:** the latter additionally reconstructs upstream PCA and publishes
verified bundles. The earlier 30.4% full-CLI CPU/Metal result has its own scope.

## Reproduction and evidence

Versions: Scanpy 1.12.4, umap-learn 0.5.12, scikit-learn 1.9.1, NumPy 2.5.3,
SciPy 1.18.1, Numba 0.67.0 and AnnData 0.13.3.post0. The isolated environment's
resolved package list and install log are retained. `Graph.swift` links the
qualified native library; `scanpy_graph.py` invokes the actual Scanpy API.
`run.py` binds source/binary/script identities before execution. `check.py`
checks every sparse neighbor and edge, graph repetitions and timing records.

`verify.py` verifies retained evidence, not fresh execution. The archive retains
all source drivers, protocols, package versions, logs and comparisons. Complete
native graphs, Scanpy sparse matrices, binary and linked native dependencies
remain externally hash-bound at `/Users/n/numivivo-metal-scverse-20260912`.
This closes the stated exact-scverse comparison for this stage/cohort; broader
workflows, other data sizes and the remaining GPU algorithms stay unqualified.
