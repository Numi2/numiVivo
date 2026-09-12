# Independent raw-count-to-graph comparison

The full Kang dataset (24,673 cells) was processed independently through Scanpy
1.12.4: original sparse H5AD counts → total-count normalization to 10,000 → log1p
→ Seurat-flavor 2,000 highly variable features / 20 mean bins → centered float64
ARPACK PCA / 20 components / seed 7 → exact sklearn neighbors / k=20 including self
→ UMAP fuzzy connectivities. Native PCA scores and feature selections were **not**
inputs to this reference computation.

All source cell barcodes and sample IDs are verified in order against each native
graph. Independent feature membership is checked against the native model. Sparse
counts stay sparse, with resident sparse storage in the reference; PCA scores and
loadings are dense low-dimensional arrays.

| Complete graph comparison | Native CPU FP64 | Native Metal FP32 |
| --- | ---: | ---: |
| Reference non-self neighbors | 468,787 | 468,787 |
| Missing / extra neighbors | 0 / 0 | 0 / 0 |
| Common edge coordinates | 706,016 / 706,016 | 706,016 / 706,016 |
| Maximum common neighbor-distance difference | 7.66e-11 | 1.48e-6 |
| Maximum fuzzy-weight difference | 4.51e-6 | 1.53e-5 |

Both pass the pre-run exact neighbor-set and edge-coordinate gate. Weight differences
are reported rather than rounded away. This extends the earlier fixed-PC graph
comparison to independently selected features and independently calculated PCA.
Native graphs come from the qualified complete Kang graph outputs; these were not
recomputed in this turn. Native PCA output equivalence across the later CPU
optimization was verified in the [borrowed traversal evidence](../BorrowedTraversal/README.md).

This is a single real-cohort numerical comparison. It does not qualify biological
labels, batch integration, perturbation prediction, UMAP coordinates or cross-library
clustering. It is not a full-workflow speed benchmark: the recorded Scanpy interval
includes reading and graph computation but excludes imports, provenance checks and
output writes, while native publication has different boundaries.

`protocol.json` and `run.py` retain every parameter and command path; `identity.py`
checks row identities and selected features. `results.json` binds native graph and
source hashes. `manifest.json` binds scripts, outputs, independent score array and
sparse NPZ graph files. Large arrays remain at
`/Users/n/numivivo-scverse-chain-20260912`. The reference uses the existing
`/Users/n/numivivo-metal-scverse-20260912/env/bin/python`, with OPENBLAS_NUM_THREADS,
OMP_NUM_THREADS, VECLIB_MAXIMUM_THREADS and NUMBA_NUM_THREADS all 1; package versions
are recorded in the [PCA comparison](../ScversePCA/versions.json). Reproduction
requires the original source/model and qualified native graph paths and a fresh
output directory. Full source tables and biological metadata are not invented or
replaced with synthetic fixtures.
