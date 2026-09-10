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
PCA reads through 16 MiB POSIX read-only private mapping windows, unmapping
the previous window during each complete operator pass. It shares the count-store
record window primitive; arithmetic still uses Float64 log values. File size, bounds and values
are checked before arithmetic traversal. Scratch creation is exclusive, rejects
symlinks and uses mode 0600. Normal completion and thrown errors unlink scratch;
the existing publication owner cleans rejected staging. Abrupt process death
is not covered by the cleanup checks.

The report records cache SHA256, bytes, selected entry count, three source passes
and actual PCA entry visits. Before cache writing, byte and work bounds reject
oversized requests. Visits equal selected entries times
`2 * min(maximumBasis, selectedFeatures) + 3 * components`. PCA rank, residual
and orthogonality gates remain unchanged. The shared source scanner admits up to
64 GiB, two million cells and four billion source entries; the resident count
projection retains its five-million-entry default. Cache and entry-visit
allowances each default to two billion. Explicit cache allowances can reach
32 billion bytes (two billion 16-byte records); explicit work allowances can
reach 1,408 billion entry visits, matching the supported basis/component bounds.
PCA score storage and downstream row readers share a two-million-row ceiling;
quality artifacts share a 512 MiB ceiling across publication, verification and
snapshot readers. These are checked storage limits, not scale qualification.
The full HIRISA execution and comparison gates are recorded in
[the execution specification](../Benchmarks/HIRISA/PCA_EXECUTION.md) and
[comparison specification](../Benchmarks/HIRISA/PCA_COMPARISON.md).

## Original full experimental matrix qualification

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

## Windowed PCA follow-up

The v2 storage method `three-source-passes-windowed-selected-COO-v2` replaces the
whole-file map with the shared 16 MiB window reader. File format, source traversal
order, HVG selection, arithmetic accumulation order and PCA algorithm are retained.
Validation and every forward/transpose traversal use windows, including rewinding
from the final window to the first. Input vector dimensions reject explicitly.

A new numerical control crosses a window boundary, repeats both operators three
times, and checks exact products and visit accounting. This complements the
complete Baron and Hagai cache sizes, which both exceed one window. Source
reconstruction and resident metadata/score/report limits remain unchanged.

`VivoOmicsFileSnapshot` now owns shared POSIX snapshot copying/hashing for this
path and the persistent count store. Each caller keeps its existing byte cap.
This also removes the Foundation read-buffer retention found in the prior
count-store benchmark. Snapshot tests check exact multi-buffer copies, rejection
before creation when oversized, and refusal to overwrite an existing snapshot.

The final release passes 66 Swift tests in 18 suites, the count-store and streamed
H5AD regression checks, H5AD projection controls, and CSR/CSC/implicit-zero PCA
controls. Both complete real cohorts pass all eight workflow assertions and
retain every native PCA field exactly relative to the published product.
Independent Scanpy comparisons retain the original aligned score errors
(Baron `2.780e-12`, Hagai `1.515e-11`) and exact 2,000-gene selections.

Same-host publication observations on Apple M4 Pro / macOS 26.6:

| Cohort | Prior seconds | Windowed seconds | Prior resident MB | Windowed resident MB |
| --- | ---: | ---: | ---: | ---: |
| Baron | 6.72 | 7.68 | 302.42 | 264.60 |
| Hagai | 13.83 | 15.00 | 345.28 | 274.53 |

These are single-run observations (decimal MB): memory decreased and wall time
increased. Verification/repeat timings and memory are retained separately. This
change does not claim a throughput improvement. The snapshot helper and mapping
changes are measured together. [Full evidence](evidence/2026-09-09-windowed/README.md)
binds source, product, baseline, complete reports, references and retained failures.

## Complete Norman scale qualification

The windowed path now passes the complete Norman filtered release: 111,445 cells,
33,694 genes and 361,582,621 source nonzeros. The predeclared plan keeps 2,000
variable genes, 20 components, a 128-vector basis and the existing `1e-6` relative
residual gate. Only the explicit work allowance increases from the default two
billion to twenty billion entry visits; cache admission remains two billion bytes.

Native PCA uses 18,354,035 selected records (293,664,560 cache bytes), traversing
5,799,875,060 entries through 16 MiB windows. Independent Scanpy checks agree on
all feature statistics within the existing tolerances and exactly the same 2,000
genes. Aligned PCA score error is `1.0971e-12`; maximum native relative residual
is `3.0565e-12`. Full native source/report reconstruction passes.

Publication took 153.94 seconds with 1,376,321,536 maximum resident bytes;
verification took 158.47 seconds with 1,493,237,760 maximum resident bytes on the
Mac mini. These include QC, condition aggregates, normalization, feature
selection, PCA and report serialization/reconstruction. The source, report and
reference hashes and original logs are retained in
[Norman evidence](evidence/2026-09-09-norman/README.md).

The independent checker now releases the source AnnData owner before promoting
count values to FP64, avoids a redundant identity-row copy, shares immutable
sparse indices for its linear-space moment check, and releases the full matrix
once the PCA subset is materialized. Public Scanpy algorithms and numerical
tolerances are unchanged. Baron/Hagai result documents remain exactly identical
and the resident JSON input path also passes the existing Kang comparison.

This is complete-cohort exploratory PCA. Fitting preprocessing on every condition
does not qualify train-only preprocessing for held-out perturbation prediction.
The roughly 118 MiB native JSON report, aggregates and metadata remain resident;
report materialization is still a scale bottleneck. Million-cell, downstream
streaming graph/integration, parallel-kernel and Metal qualification remain open.
Native and reference hosts differ, so their timings are not a speed comparison.

The [standalone PCA bundle](PCA_BUNDLE.md) now separates this algorithm from
pseudobulk aggregation and stores scores/loadings as binary records. It shares
cell-quality accumulation and the same PCA engine with the legacy route.
