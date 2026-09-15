# Matched raw-count-to-graph timing and memory

This is a fresh three-repetition comparison on the same frozen **24,673-cell,
15,706-feature Kang H5AD**. Native NumiVivo runs the actual product commands
`singlecell-h5ad-pca` and `singlecell-pca-neighbors` from the raw sparse H5AD;
the CPU path uses the default FP64 owner and the Metal path uses the qualified
FP32 owner. Scanpy independently runs sparse raw counts → total-count
normalization to 10,000 → log1p → Seurat-flavor 2,000-feature HVG selection
with 20 mean bins → float64 ARPACK PCA (20 components, seed 7) → exact
sklearn-backed UMAP fuzzy neighbors (20 neighbors, one worker).

## Current PBMC3K receipt (2026-09-15)

The [current PBMC3K receipt](evidence/2026-09-15-pbmc3k-current/README.md)
executes the same raw-count-to-PCA-to-neighbor route on the pinned public
2,700-cell, 32,738-feature source. The native and independent Scanpy paths use
the same 2,000 selected features and 20 components; the final repetition has
zero missing or extra neighbors/edges and exact coordinates for all 77,916
fuzzy-graph edges. Native CPU and Metal binary graphs then pass verification,
and seeded Louvain replay at seeds 7, 19 and 42 gives identical label vectors
and pair relations (ten communities, zero disconnected communities).

The receipt retains the earlier failed namespace comparison: the first Scanpy
driver exposed display-name indices while the native mapping used `gene_ids`.
That failure is preserved, and the driver now selects the explicit mapped
feature-ID column. Median pipeline times are 2.0487 s CPU, 2.0197 s Metal and
6.1915 s Scanpy; these are one physical M4 Pro run with three repetitions and
separate output formats. This closes current numerical interoperability and
descriptive clustering for one unreported library. Donor replication,
biological preservation and outcome prediction remain unqualified.

## Measured result

| Path | Pipeline seconds (three runs) | Median | Peak RSS (three runs) | Median RSS |
| --- | --- | ---: | --- | ---: |
| Native CPU FP64 | 13.708, 13.957, 13.982 | **13.957 s** | 399.6, 399.0, 398.9 MiB | **399.0 MiB** |
| Native Metal FP32 | 10.795, 10.843, 10.936 | **10.843 s** | 410.9, 410.8, 432.3 MiB | **410.9 MiB** |
| Scanpy exact | 8.549, 8.178, 8.228 | **8.228 s** | 1,063.7, 1,082.1, 1,067.8 MiB | **1,067.8 MiB** |

For the native product path, Metal is **22.3% lower elapsed time (1.29x)** than
the default CPU path on this cohort. Its median RSS is 3.0% higher, so no
memory reduction is claimed. Scanpy is 41.0% lower than native CPU and 24.1%
lower than native Metal in this measured process boundary, while using 2.68x
the native CPU median RSS. These numbers are one physical M4 Pro desktop,
three fixed-order repetitions, and exercised caches; they are not confidence
intervals, cold-cache results, or universal backend speed claims.

## Numerical comparison

The final repetition passes exact cell-barcode and sample-ID identity checks,
ordered selected-feature checks for all 2,000 features, all 468,787 non-self
neighbor memberships, and all 706,016 fuzzy-graph edge coordinates for both
native backends versus Scanpy. Maximum common neighbor-distance differences are
7.66e-11 (CPU) and 1.48e-6 (Metal); maximum fuzzy-weight differences are
4.51e-6 (CPU) and 1.53e-5 (Metal). Weight and distance differences are retained,
not rounded away. The three Scanpy distance/connectivity matrices, PCA scores,
and identity files also match byte-for-byte after excluding the per-run timing
field in `report.json`.

The timing boundary includes process startup, source reads, arithmetic and
output writes. Native verification commands run after each timed pipeline and
are recorded as required checks, but are outside the elapsed interval. Native
publication snapshots the source and writes a provenance-bearing PCA and graph
bundle; Scanpy writes sparse matrices, scores and an identity report. The
algorithms and source cohort are matched, while these output formats and
bookkeeping costs remain different. Peak RSS is `/usr/bin/time -l` maximum
resident set size across the native PCA and graph commands, or the complete
Scanpy process.

This closes the stated end-to-end comparison for this cohort and stage. It does
not qualify model-fitting acceleration, million-cell scaling, approximate
neighbors, downstream clustering or embedding, clinical/biological outcome
prediction, or a production-default change. CPU FP64 remains the default.

## Reproduction and retained evidence

The source drivers are `run.py`, `scanpy.py` and `compare.py`. Run `run.py` only
with a new output directory; it retains partial output after a failure and
refuses to overwrite it. The authoritative retry is retained at
`/Users/n/numivivo-pca-end-to-end-20260914-retry` on the Mac mini, with protocol,
plans, all nine outputs and logs, native verification receipts, timing records
and `comparison.json`. The first driver-shadowing failure remains at
`/Users/n/numivivo-pca-end-to-end-20260914` and is not counted as a result.

The external inputs are the frozen source at
`/Users/n/numivivo-pca-borrowed-20260912/pca-1/original.h5ad`, its `plan.json`,
the qualified binary at
`/Users/n/numivivo-native-metal-knn-20260912/build/numivivo-omics`, and the
isolated Scanpy 1.12.4 environment at
`/Users/n/numivivo-metal-scverse-20260912/env`. `protocol.json` binds their
SHA-256 identities. This is a benchmark receipt, not evidence that NumiVivo
can predict a biological outcome from an unseen genome, cell state, treatment,
or patient.
