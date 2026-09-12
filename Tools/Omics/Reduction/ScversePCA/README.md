# Matched-feature sparse PCA against Scanpy

Executed 2026-09-12 on physical M4 Pro using Scanpy 1.12.4 and its sparse ARPACK
PCA implementation. All 24,673 Kang and 13,863 Hagai cells are included. The native
outputs are from the borrowed-buffer implementation published in `105022c1`.

The independent reference reads each original H5AD, retains the sparse count matrix,
normalizes all source features to the native target 10,000, applies log1p, then selects
the exact 2,000 native-selected features in the recorded order. Scanpy performs
centered float64 PCA with 20 components, ARPACK and seed 7. This checks PCA given
the selected features; it does **not** independently qualify highly-variable-feature
selection. No dense cells × genes matrix is constructed; the reference does keep
the sparse matrix resident and produces dense cells × 20 and features × 20 arrays.

Cell barcodes and sample IDs are checked against the explicit native import mapping.
The first driver incorrectly assumed every barcode was the H5AD row name and stopped
before Hagai PCA. Hagai instead uses its `barcode` column, while row names include
sample prefixes. The failed attempt and its successful Kang result are retained.
The corrected driver verifies the mapped columns without changing data or native outputs.

| Full cohort | Minimum loading-subspace singular value | Maximum relative variance difference | Maximum sign-aligned score difference |
| --- | ---: | ---: | ---: |
| Kang | 0.9999999999999987 | 5.11e-15 | 2.09e-10 |
| Hagai | 0.9999999999999991 | 3.33e-15 | 1.14e-9 |

Both pass the pre-run subspace threshold 0.99999 and relative variance threshold 1e-5.
Loading differences after sign alignment are at most 1.84e-11 and 1.14e-10;
normalization-center differences are at most 9.33e-14 and 3.11e-14, respectively.
Signs are aligned because PCA eigenvectors have arbitrary orientation. Close subspaces
and variances establish numerical agreement, not biological prediction accuracy.

The corrected single-run Scanpy PCA-only times are 0.252 s Kang and 0.304 s Hagai.
Reading, sparse normalization and feature slicing took a separate 1.163 s and 0.896 s.
These exclude interpreter/import startup, native HVG selection, source hashing,
publication and output writing. They are **not comparable to native full-CLI medians**,
and no native/scverse speed ratio is claimed. A matched stage benchmark and memory
comparison remain necessary before a PCA performance claim against scverse.

`run.py` records the exact API call and runtime paths. `protocol.json`, `versions.json`,
`results.json` and logs retain execution evidence. `manifest.json` SHA-256 binds these
files, the reference NPZ arrays and the installed Scanpy PCA source snapshot. Large
reference arrays and that source snapshot remain external at
`/Users/n/numivivo-pca-scverse-20260912`; native model, score and loading hashes are in
`results.json`. Run with the existing
`/Users/n/numivivo-metal-scverse-20260912/env/bin/python`, setting OPENBLAS_NUM_THREADS,
OMP_NUM_THREADS, VECLIB_MAXIMUM_THREADS and NUMBA_NUM_THREADS to 1, in a fresh output
workspace after updating the driver's root. The original full datasets remain required.
