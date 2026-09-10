# Approximate PCA neighbors with measured recall

`singlecell-pca-neighbors` supports optional HNSW search over fitted or frozen-query
PCA bundles. Omitting `approximation` retains the exact Dispatch implementation.
Set `storage: "binary"` for the [streamed binary graph store](GRAPH_STORE.md);
omission preserves resident JSON output.
The native CPU implementation uses [hnswlib](https://github.com/nmslib/hnswlib)
v0.8.0, revision `3f3429661187e4c24a490a0f148fc6bc89042b3d`, with unmodified
Apache-2.0 headers, license and per-file hashes under `Sources/NumiVivoCore/ThirdParty/`.
The algorithm is described by [Malkov and Yashunin](https://arxiv.org/abs/1603.09320).
Python is used only for independent qualification.

```sh
numivivo singlecell-pca-neighbors pca-bundle --plan hnsw.json --output graph
numivivo singlecell-pca-neighbors-verify graph
```

```json
{
  "schemaVersion": 1,
  "inputKind": "fitted",
  "neighbors": {"neighbors": 15},
  "execution": {"workers": 1},
  "approximation": {
    "connections": 16,
    "constructionWidth": 200,
    "searchWidth": 128,
    "seed": 7,
    "maximumDistanceEvaluations": 500000000,
    "scoreCacheBytes": 33554432
  }
}
```

Use `inputKind: "query"` for a frozen query bundle. This builds neighbors among
query cells; it does not perform cross-reference annotation. HNSW requires one
worker explicitly; the exact mode's default of four workers is rejected here.
Exact-mode query/candidate block settings do not control HNSW's fixed 256-row tiles.
Integrated representations are rejected for these PCA inputs.

## Native ownership and limits

The C++ index stores UInt32 row identities. Its custom metric reads complete FP64
score records through a bounded LRU cache, computes squared Euclidean distances
without fused contraction, and returns Euclidean distances to Swift. Every decoded
record is checked for row, component and finite value. At least two score tiles
are retained so loading the second distance operand cannot invalidate the first.
The cache holds decoded values; temporary encoded reads use at most 256 KiB.
Cache accounting excludes container and allocator overhead.

Construction inserts rows serially in source order with the recorded seed.
Queries include self in slot zero, followed by retrieved neighbors sorted by
distance and row index. Search is approximate: global membership and cutoff-tie
selection can differ from exact search. No cross-toolchain bitwise guarantee is
made. The index is rebuilt for verification rather than serialized; reported
`indexStorageBytes` is hnswlib's serialized-size estimate and excludes allocator,
mutex and container overhead.

Limits are 1–64 PCs, k=2–128, connections=8–64, construction width from connections
through 512, search width from k through 1024, and 1–32 billion metric evaluations.
The Swift cache bound is 256 KiB–1 GiB. Defaults remain 500 million evaluations
and 32 MiB; larger allocations require explicit plan values. Shared C/Swift
constants admit up to two million rows. Streamed neighbor entries follow the
row and k bounds; symmetric graph bytes follow `rows × 2 × (k−1) × 16`. Resident JSON output remains capped at four million neighbor entries. Binary
output emits rows and constructs connectivity through a disk transpose and merge,
without resident neighbor/edge arrays. The HNSW index, identities, cell-scale
bookkeeping and part of input PCA reconstruction remain resident. Neither storage
mode alone establishes million-cell qualification. The complete HIRISA case
below supplies separate measured evidence.

`maximumDistanceEvaluations` counts all construction and query metric calls,
including repeated pairs. It is separate from exact-mode `maximumDistancePairs`.
Approximate `graph.distancePairs` and `execution.directedDistanceEvaluations` both
report that actual total; exact graphs retain their unique-pair convention.
`execution.hnsw` separates construction/query work, read bytes/calls, cache hits,
peak decoded cache bytes and estimated index bytes. Work-budget, I/O, malformed
record, nonfinite-distance and allocation failures reject publication. Cancellation
is checked between insertions and queries. Metric errors are propagated after
the upstream call returns so they do not unwind through borrowed visited lists.

The immutable lifecycle is shared with [exact neighbors](WINDOWED_NEIGHBORS.md):
reconstruct input, build graph, bind plan/input/graph/execution hashes, and rebuild
on verification. Executable-bound older PCA receipts must be refitted using the
current binary; changing receipt hashes is not a migration.

## Measured evidence

[Archived qualification](evidence/2026-09-09-hnsw/README.md) covers full Baron,
Hagai and Norman cohorts, each with 20 PCs and k=15. The predeclared gate was mean
strict recall ≥0.95 and fifth-percentile strict recall ≥0.85. Baron and Hagai use
every cell against independently qualified exact graphs. Norman uses a fixed
uniform sample of 2,048 queries, each searched exactly against all 111,445 cells.

| Cohort | Cells | Mean strict recall | Fifth percentile | Publish seconds | Peak resident bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Baron | 8,569 | 0.9999666572 | 1.0 | 6.85 | 222674944 |
| Hagai | 13,863 | 0.9999896951 | 1.0 | 15.11 | 300728320 |
| Norman | 111,445 | 0.9994070871 (sampled) | 1.0 (sampled) | 161.56 | 1357234176 |

These are single-run Mac mini M4 Pro lifecycle measurements including PCA input
reconstruction and output serialization, not isolated search timings or a same-host
scverse comparison. No small-cohort speed advantage is established. Norman used
467,109,415 metric calls, a 17,831,200-byte decoded score cache and an estimated
16,998,788-byte index. Its final graph contains 2,513,474 connectivity entries.
All score values fit within the default cache in these three cohorts; a separate
600-row, 64-component regression forces tile eviction and checks exact results.

Independent NumPy checks every returned distance; umap-learn checks the fuzzy
graph on the returned approximate neighbors. Maximum distance error is below
3.6e-15 and maximum fuzzy-weight error below 7.6e-6. Fuzzy topology matches that
reference exactly; this does not imply exact-neighbor topology. Minimum strict
recall is 12/14, 13/14 and 12/14 respectively (Norman sampled). No parameter tuning
or threshold relaxation was needed. All 75 Swift tests in 22 suites, both builds,
20 HNSW fixture commands and 22 query-input regression commands passed. A full
Baron exact-mode regression preserves the prior published graph bytes.

`run_hnsw.py` records native lifecycle commands. `check_hnsw_reference.py` validates
receipt bindings, distances, declared recall, bandwidths, kernel mass and graph
connectivity; it retains sampled exact neighbors and per-query recall.
`check_hnsw.py` covers fitted/query inputs, replay, cache invariance, small-fixture
exact results, resource rejection and rehashed tampering. These checks establish
numerical behavior and bounded recall evidence, not biological generalization,
downstream clustering/embedding quality, Metal execution or million-cell readiness.

## Complete HIRISA scale qualification

The [frozen HIRISA experiment](../Benchmarks/HIRISA/GRAPH_RESULTS.md) now passes
full native publication/reconstruction for 1,612,594 original cells, 20 PCs and
k=15. All returned distances and 34,707,084 fuzzy edge records pass independent
checks. The unchanged fixed 2,048-query panel has mean strict recall
0.9997209821428572 and fifth percentile 1.0 against all cells. Total metric work
is 6,963,750,537 evaluations; decoded score cache peaks at 258,015,040 bytes and
whole publication RSS at 4,360,978,432 bytes. This is a bounded complete-cohort
CPU graph qualification; index/metadata residency, downstream biological quality
and controlled performance remain separately measured work.
