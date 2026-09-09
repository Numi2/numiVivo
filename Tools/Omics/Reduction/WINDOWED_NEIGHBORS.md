# File-backed, parallel exact PCA neighbors

`singlecell-pca-neighbors` builds exact Euclidean kNN and the existing
UMAP-compatible fuzzy graph from a [fitted PCA bundle](PCA_BUNDLE.md) or a
[frozen query bundle](PCA_QUERY.md). It uses the shared neighbor heap, distance/index
tie order and fuzzy-graph calculation. It does not reconstruct a dense cells-by-genes
or cells-by-cells matrix.

```sh
numivivo singlecell-pca-neighbors pca-bundle --plan neighbors.json --output graph
numivivo singlecell-pca-neighbors-verify graph
```

A plan is explicit about the input representation:

```json
{
  "schemaVersion": 1,
  "inputKind": "fitted",
  "neighbors": {
    "neighbors": 15,
    "maximumDistancePairs": 200000000
  },
  "execution": {
    "workers": 4,
    "queryBlockRows": 128,
    "candidateBlockRows": 2048
  }
}
```

Use `inputKind: "query"` for frozen query scores. `neighbors` uses the existing
neighbor options, including local connectivity; integrated representation is
rejected because these inputs contain PCA scores. Self is slot zero, followed by
neighbors sorted by distance and source row index. Duplicate or tied points retain
that deterministic convention.

## Execution and memory

Each Dispatch worker owns an independent score reader, query tile, candidate tile
and bounded nearest-neighbor heaps. The reader maps at most 16 MiB of the immutable
score file. Each query block scans all candidate rows; score records are validated
against their exact row/component coordinates. Workers never share mutable heaps,
readers or floating-point accumulators. A locked result holder propagates success
or failure, and the caller merges blocks in row order after each batch joins.
Cancellation is checked between bounded parallel batches on the calling task.

Worker limits are 1-16, query tiles 1-512 rows and candidate tiles 1-8192 rows.
Final fixed-width kNN arrays are bounded to four million entries. Metadata, cell
identities, the final fuzzy graph and its JSON encoding remain resident; input PCA
reconstruction has its own resident fitting costs. Thus only score access and search
tiles are bounded independently of cohort size. Graph storage is not yet out of core.

Rows are independently owned to avoid cross-worker heap synchronization. This
computes both directed distances `(i,j)` and `(j,i)`, while the resident triangular
search computes each pair once. `graph.distancePairs` retains its existing unique
unordered-pair convention. `execution.json` explicitly records the doubled directed
evaluations, scalar distance terms and score-record reads. Serial and parallel
file-backed modes perform the same work; neither promises a speedup over the
triangular resident implementation.

Exact search remains quadratic. The existing maximum of 500 million unique pairs
still applies (50 million by default), so the complete 111,445-cell Norman cohort
and million-cell graphs are not supported by this exact route. Approximate search
with measured recall is required for that scale; increasing memory bounds does
not remove this computational limitation.

## Artifact lifecycle

The output retains a complete verified `input/` bundle, canonical `plan.json`,
`graph.json`, deterministic `execution.json` work accounting and `receipt.json`.
Graph construction first copies and reconstructs the input. CLI receipts bind the
exact executable and OS identity; older inputs must be rebuilt from their retained
sources using the current executable, not merely rehashed. Verification rebuilds
both the input and graph, then checks all receipt-bound files. Rehashed changes to
graph weights, execution counts or input scores are rejected. Publication requires
a new destination and cleans private staging after errors.

`graph.json` has the same `VivoSingleCellNeighborGraph` representation consumed by
the existing clustering and embedding implementations. This command does not run
those downstream algorithms. Query-only graphs describe proximity among the query
cells; they are not cross-reference label mapping or donor integration.

## Qualification tools

`check_window_neighbors.py` exercises fitted/query input, one/two/four workers,
uneven tiles, exact graph bytes, repeats and negative controls.
`run_window_neighbors.py` measures native publication with one and four workers
and reconstruction on complete Baron and Hagai cohorts and a held-out Baron donor.
`window_neighbors_reference.py` validates binary input/receipt bindings and feeds
all scores/edges to `check_neighbors.py`: independent SciPy exact membership/order
and distances, plus umap-learn fuzzy topology, bandwidths and weights.

The external Python checker materializes score arrays and bounded distance blocks;
it is qualification tooling, not the native execution path. End-to-end command
measurements include source/PCA reconstruction and graph serialization. They are
not isolated kernel benchmarks or biological validation.

## Complete-cohort results

The [archived qualification](evidence/2026-09-09-windowed-neighbors/README.md)
passes complete Baron (8,569 cells), Hagai (13,863 cells) and held-out Baron human3
(3,605 cells), each with 20 PCs and k=15. Serial and four-worker graph bytes are
identical; SciPy membership/order/distances match exactly. All fuzzy topology
matches umap-learn, with maximum weight error below 5.7e-6.

One/four-worker publication took 7.15/6.61 s, 15.92/14.66 s and 5.08/5.02 s
respectively on the Mac mini. Four-worker resident peaks were 233, 304 and 170 MB.
These single-run lifecycle timings include PCA reconstruction, which dominates
cost. They do not establish kernel-only speedups or a cross-platform advantage.
The gate passed 73 Swift tests in 21 suites, 13 successful real-data commands,
19 graph-fixture commands and 22 query-input regression commands. Initial compiler
and stale-input failures are retained alongside their successful follow-ups.
