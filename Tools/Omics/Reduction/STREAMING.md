# Streamed H5AD normalization, feature selection and PCA

Add `"reduction": {}` to the existing `singlecell-h5ad-pseudobulk` plan alongside
its explicit `mapping` and `contrasts`. Optional settings are:

```json
"reduction": {
  "normalizationTarget": 10000,
  "pca": {"highlyVariableFeatures": 2000, "components": 20, "maximumBasis": 128},
  "maximumCacheBytes": 2000000000,
  "maximumEntryVisits": 2000000000
}
```

Publish with `singlecell-h5ad-pseudobulk SOURCE --plan PLAN --output BUNDLE`;
reconstruct with `singlecell-h5ad-pseudobulk-verify BUNDLE`. Plans without reduction
retain their encoding and behavior. The bundle archives the original H5AD,
plan, report and receipt. Replay reconstructs the scratch representation.

## Ownership and storage

`VivoH5ADPseudobulk` owns the snapshot, QC, publication and reconstruction.
`VivoH5ADReduction` adds two source scans after QC: all-feature moments/HVG
selection, then selected-entry cache writing. Library-size normalization uses
all genes and retains zero-count cells. `VivoSingleCellReduction` shares its
original Seurat-style HVG selection and centered Krylov/Jacobi PCA between the
resident sparse representation and streamed operator. No cells-by-genes or
full feature covariance matrix is materialized.

The private cache uses 16-byte little-endian records: UInt32 row, UInt32 selected
column and Float64 log-normalized value. Writes use a buffer capped at 1 MiB;
PCA reads an explicit POSIX read-only private mmap. File size, bounds and values
are checked before arithmetic traversal. Scratch creation is exclusive, rejects
symlinks and uses mode 0600. Normal completion and thrown errors unlink scratch;
the existing publication owner cleans rejected staging. Abrupt process death
is not covered by the cleanup checks.

The report records cache SHA256, bytes, selected entry count, three source passes
and actual PCA entry visits. Before cache writing, byte and work bounds reject
oversized requests. Visits equal selected entries times
`2 * min(maximumBasis, selectedFeatures) + 3 * components`. PCA rank, residual
and orthogonality gates remain unchanged. The source scanner retains its
1 GiB / 100 million nonzero bounds; the resident count projection retains its
five-million-entry default. Cache bytes are capped at two billion and the
entry-visit default is two billion (configurable up to twenty billion).

## Full experimental matrix qualification

Both prepared sources retain the complete deposited count scope established by
the existing [benchmark suite](../Benchmarks/README.md). Reference code uses
AnnData/Scanpy, all source genes and exact cell identities, with sparse ARPACK
PCA. It materializes a sparse reference matrix; it is not the native storage path.

| Dataset | Cells × genes | Source nonzeros | Selected cache bytes | Aligned PCA score error |
| --- | ---: | ---: | ---: | ---: |
| Baron human, CSR | 8,569 × 20,125 | 16,171,764 | 18,358,208 | 2.780e-12 |
| Hagai mouse, CSC | 13,863 × 22,048 | 32,848,185 | 53,794,864 | 1.515e-11 |

Both select exactly the same 2,000 genes as Scanpy and fit 20 components. Maximum
native relative residuals are 9.785e-12 and 2.981e-11. All-feature means,
variances, bins and dispersions also pass explicit tolerances, including
Scanpy's float32 normalized-dispersion representation.

Each production workflow passes eight assertions: publish, reconstruct, exact
repeat, absence of retained scratch, cache/work/basis rejection and staging
cleanup. A separate numerical fixture checks exact CSR/CSC PCA equality and
zero-cell/zero-feature handling. The full Kang resident integration, graph,
clustering and UMAP report is exactly unchanged. All 37 Swift tests in ten suites,
11 existing Kang workflow assertions and 18 general CLI assertions pass.
The initial Swift initializer compilation failure is retained alongside the
successful scoped build, tests and full release build.

Observed production publication took 7.64 seconds for Baron and 14.72 seconds
for Hagai, with `/usr/bin/time -l` peak footprints 282,149,728 and 303,694,664 bytes
on this laptop. These are single-run observations, not CPU/scverse performance
comparisons. Commands, runtime/source hashes, references and rejection stderr
are retained in [evidence](evidence/2026-09-09-streaming/source-state.json).

```sh
python check_streamed_cli.py --binary /path/to/numivivo --prepared /data/prepared --out /new/native
python check_reference.py --h5ad /data/prepared/prepared.h5ad --report /new/native/bundle/report.json --out /new/reference.json
python check_streamed_formats.py --binary /path/to/numivivo --out /new/formats
```

## Remaining scope

Metadata, QC, aggregates, moments, basis, scores, loadings and JSON reports remain
resident. This route currently stops at PCA; streamed graph/integration wiring,
million-cell execution, parallel kernels and GPU performance remain unqualified.
Numerical agreement on these matrices does not establish cell annotation,
biological preservation, disease effects or held-out perturbation prediction.

Optional [fixed expression programs](../Programs/README.md) share the archived
source and add one separate scan. When programs and PCA are enabled together,
there are four scans in total; the reduction component still accounts for its
three scans, including QC.
