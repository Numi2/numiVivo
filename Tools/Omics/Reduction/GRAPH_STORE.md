# File-backed PCA neighbor and fuzzy graph storage

Set `"storage": "binary"` in a `singlecell-pca-neighbors` plan to store neighbors
and connectivity as binary records. Exact Dispatch and approximate HNSW searches
both emit rows incrementally. Omitting `storage`, or using `"json"`, preserves the
resident JSON graph route and its four-million-neighbor-entry limit.

```json
{
  "schemaVersion": 1,
  "inputKind": "fitted",
  "neighbors": {"neighbors": 15},
  "execution": {"workers": 1},
  "approximation": {"connections": 16, "constructionWidth": 200, "searchWidth": 128, "seed": 7},
  "storage": "binary"
}
```

```sh
numivivo singlecell-pca-neighbors pca --plan plan.json --output graph
numivivo singlecell-pca-neighbors-verify graph
```

For exact search omit `approximation` and choose the Dispatch worker count and
pair-work budget explicitly. `inputKind: "query"` supports frozen query bundles.
These are graphs among query cells, not cross-reference annotation.

## Record contract

Every record is exactly 16 bytes, little endian: UInt32 row, UInt32 column, then
64 payload bits. Integers and floating-point values must be decoded according to
the file's role. The files are not sparse count matrices.

| File | Ordering / coordinates | Payload | Number of records |
| --- | --- | --- | ---: |
| `neighbors.bin` | Fixed-width cell-major; column is neighbor row ID; self first, then distance/index order | FP64 Euclidean distance | cells × k |
| `bandwidths.bin` | Cell-major, columns 0/1/2 | FP64 rho / sigma / absolute kernel-mass residual | cells × 3 |
| `offsets.bin` | Rows 0 through cells inclusive; column zero | UInt64 CSR entry offset | cells + 1 |
| `edges.bin` | Source row, then strictly increasing target row; no self edge | FP64 symmetric fuzzy-union weight | connectivity entries |

`graph.json` is a small `VivoPCAGraphStoreReport` with format identifier,
algorithm/options, dimensions, counts, components, isolation and all four file
fingerprints. Cell identities are stored once in `input/metadata.json`, whose
receipt is bound by the parent graph receipt. Consumers must branch on plan
storage/format; binary `graph.json` is not `VivoSingleCellNeighborGraph` JSON.

The report preserves existing work semantics: exact `distancePairs` counts unique
unordered pairs; approximate `distancePairs` counts all construction/query metric
evaluations including repeats. Approximate recall remains an independent check.

## Execution and memory

Search emits one row at a time into the shared 1 MiB record writer. HNSW has a
synchronous C row callback with explicit output-failure propagation. Exact search
emits rows after each bounded Dispatch batch; it retains only that batch's heaps
and results. No fixed-width all-cell neighbor array is allocated in binary mode.

Graph construction makes a sequential pass for global mean distance, then computes
bandwidths with the same function as resident graph construction. It writes positive
directed weights and counts incoming degrees. A prefix sum assigns disk locations;
source-major traversal scatters incoming records through checked positional writes.
Each incoming row is consequently ordered by source ID. A sequential row merge
combines outgoing and incoming weights using the existing fuzzy-union formula.
Output edges and offsets are written immediately, without graph-wide dictionaries,
edge arrays, queues or a graph-sized JSON encoding. Union-find counts components.
Private transpose files are removed on success and failure.

Readers map at most 16 MiB each. Only k outgoing edges are sorted at a time; even
an incoming hub is streamed without materializing its full row. Degree, prefix and
component arrays remain resident and scale with cell count. Input PCA metadata,
fitting state and the HNSW index retain their own resident costs. The reported
HNSW cache/index sizes do not measure total process memory.

Binary mode admits up to one million cells with k=2–128, subject to existing
source/PCA and search-work bounds. This permits up to 128 million neighbor records
and 254 million symmetric edge records; it does not certify successful million-cell
execution. Exact search remains quadratic and capped at 500 million unique pairs.
HNSW's metric budget remains separately bounded. Disk capacity is required for
retained input, output records and temporary directed/transpose files.

Publication reconstructs the PCA input before search, binds every output hash and
atomically moves the private directory. Verification reconstructs input and graph
and checks every binary file. Rehashing a changed binary file and its graph/parent
receipts cannot turn altered results into a valid reconstruction. Rebuilt
executables must refit PCA from original sources; old receipts are not relabeled.

## Qualification and remaining work

`check_graph_store.py` covers fitted/query inputs in exact and HNSW modes,
byte-identical replay, every binary numeric field versus JSON, overwrite/format
rejection and rehashed modifications to each binary file. Swift regressions cover
duplicates, zero-rho bandwidths, disconnected components and row-sink failures.
`run_graph_store.py` records full Baron, Hagai and Norman native lifecycles and a
same-executable JSON baseline. `check_graph_store_reference.py` compares every
record, FP64 bit, cell identity and PCA score against the previously independently
qualified full-cohort graphs. Real measurements are recorded in the accompanying
evidence archive and table below.

[File-backed clustering](FILE_CLUSTERING.md) consumes this binary store directly,
including windowed edges at every aggregation level. [File-backed embedding](FILE_EMBEDDING.md)
also consumes this store with a windowed mutable schedule. File-backed integration,
million-cell end-to-end qualification, biological stability and Metal/scverse
performance comparisons remain separate required work. This store removes graph
materialization from construction/publication; it does not complete the full
single-cell roadmap.

## Full-cohort measurements

The [archived qualification](evidence/2026-09-09-graph-store/README.md) checks every
binary record against the previously independently qualified full graphs. All
neighbor IDs/distances, bandwidths, CSR coordinates/weights, PCA score bytes and
cell identities are exact for Baron, Hagai and Norman, plus full-Baron exact
search. The same executable's JSON graph payloads also match the previous release
byte for byte. Norman retains the previous 2,048-query sampled recall scope;
its 1,671,675 neighbor and 2,513,474 connectivity records are all checked here.

| Cohort | Binary publish seconds | Binary peak resident bytes | JSON publish seconds | JSON peak resident bytes | Binary verify seconds | Verify peak resident bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baron, 8569 cells | 6.82 | 202817536 | 6.83 | 229736448 | 6.85 | 204668928 |
| Hagai, 13863 cells | 14.80 | 226689024 | 14.88 | 310018048 | 14.87 | 224559104 |
| Norman, 111445 cells | 157.97 | 477413376 | 161.43 | 1420509184 | 160.15 | 468582400 |

Norman's binary mode reduces the observed peak by 66.4% relative to JSON from the
same executable. Its four graph files total 74,094,880 bytes. These are single-run
native CPU lifecycle observations on M4 Pro/24 GiB/macOS 26.6, including input PCA
reconstruction and serialization. They do not establish a general speedup or a
same-host scverse comparison. Completed artifacts were copied for independent
checks while later cases ran; there was no competing native build or benchmark.

The original 75 Swift tests and two new graph-store tests pass (77 in 23 suites).
Fourteen full-cohort native commands, 22 binary graph fixture commands and 22
query-input regression commands pass. The fixture suites include expected
rejections; the archive preserves their actual statuses. Unrelated genomic CLI
routing changes arrived upstream during validation. The merged executable was
rebuilt and received fresh query-input and graph-store CLI checks; the full-cohort
measurements above remain explicitly bound to the earlier benchmark binary,
whose omics source bytes are unchanged in publication.
