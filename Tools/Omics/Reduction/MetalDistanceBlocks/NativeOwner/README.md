# Explicit Metal FP32 backend in the native neighbor owner

The product PCA neighbor command now supports `execution.backend: "metalFP32"`.
Omit that field for the unchanged CPU FP64 default. Metal requires one worker and
cannot be combined with HNSW. Both binary and resident JSON graph routes record
an explicit FP32 method and precision qualification in the graph and execution
report. Existing default encoding omits the new optional field.

```json
{
  "schemaVersion": 1,
  "inputKind": "fitted",
  "storage": "binary",
  "neighbors": {"representation": "pca", "neighbors": 20, "maximumDistancePairs": 100000000},
  "execution": {"backend": "metalFP32", "workers": 1, "queryBlockRows": 128, "candidateBlockRows": 8192}
}
```

Use `singlecell-pca-neighbors PCA --plan PLAN --output NEW_BUNDLE` and
`singlecell-pca-neighbors-verify BUNDLE`. Build the scoped actual CLI with
`Tools/Omics/H5AD/build.sh OUTPUT --with-cli`. Full app/package builds are not
qualified here.

## Bounded owner and numerical contract

Each command owns a queue, pipeline, one windowed score reader and three reused
shared buffers. Query/candidate buffers contain FP32 conversions of bounded FP64
read tiles. The distance buffer has queryRows × candidateRows FP32 elements
(4 MiB for the plan above), independent of total cohort size. Existing limits
remain queryRows 1–512, candidateRows 1–8192 and the unchanged pair-work budget.
No whole PCA score array or cells-by-genes matrix is allocated by this kernel.
Metadata, graph assembly and source reconstruction retain their separate bounds.

Candidate search is exhaustive. Distances use explicit FP32 arithmetic with fast
math disabled, widen to FP64 before square root, and enter the existing row/index
heap and fuzzy graph assembly. Every command must complete successfully before
its output is read or buffers reused. Cancellation is checked before each tile,
after GPU completion and before emitting rows. Conversion overflow, nonfinite
source values and nonfinite squared distances reject. There is no silent CPU
fallback. Physical Apple GPU admission matches the existing normalization owner.

## Executed checks

All 13,863 real Hagai cells and 277,260 neighbor entries pass actual owner
execution; indices and widened distances match the qualified FP32 reference
exactly. The owner driver takes 0.872 s for this stream and peaks at 37.4 MB for
the whole check process. This is not full CLI timing or a comparative benchmark.

Fresh actual CLI PCA publication passes, followed by CPU binary, Metal binary
and Metal JSON graph publication and source reconstruction. All ten CLI checks
pass, including invalid worker count, unknown backend and Metal/HNSW combination
rejection before publication. Single-run publication times are 9.41 s CPU binary,
8.52 s Metal binary and 8.59 s Metal JSON, including upstream verification and
output handling. These are not repeated speedup measurements.

All neighbor memberships and 395,338 edge coordinates match CPU FP64. One row
changes neighbor order. Metal binary and JSON neighbors, distances and weights
agree exactly. Maximum fuzzy-weight difference from CPU is 6.24e-6; the graph
retains one component and zero isolated cells. This is numerical graph evidence,
not embedding, clustering, integration or biological-preservation validation.

Focused controls pass for duplicate-distance ties and partial query/candidate
tiles, FP32 conversion overflow, squared-distance overflow and nonfinite input.
A real owner task cancels successfully after a 100 ms request delay. Work-budget
rejection and unchanged default execution keys pass. These checks do not prove
cross-device replay, concurrent source-mutation safety or exhaustive memory
safety. Million-cell execution remains limited by the existing pair budget and
other pipeline stages; it is not claimed here.

## Retained evidence

`verify.py` checks the archive, current owner source hashes and complete receipts.
It does not launch a fresh GPU run. The archive retains executed drivers, CLI
plans/logs, source/build hashes and numerical comparisons. Large binaries, input
snapshots and graph arrays remain hash-bound in its external manifest under
`/Users/n/numivivo-native-metal-knn-20260912`.

The first test-driver compilation incorrectly called Swift's `keys` property;
its source and error log are retained. The corrected driver passed. Product
source compilation passed with retained deprecation warnings. No failed run was
reclassified as a successful one. Downstream biological checks and broader
CPU/scverse performance comparison remain the next qualification gates.
