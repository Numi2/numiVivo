# Native Metal PCA-neighbor regression (2026-09-14)

The source revision `0ca4d1efb3e0fd0be95643ff6d54704705eb90a4` was built and
executed on the physical Apple M4 Pro Mac mini with
`NUMIVIVO_TEST_METAL=1 swift test --filter SingleCellWindowedNeighborTests`. All
three tests passed, including the new physical Metal regression.

The new test writes a bounded 257-by-4 sparse PCA-score fixture, runs
`VivoWindowedPCANeighbors` with `execution.backend = .metalFP32`, one worker,
and 23-by-31 query/candidate tiles, then compares the result with the resident
CPU FP64 oracle. Integer-valued coordinates make the membership and tie
structure exact while the run exercises FP32 conversion, Metal upload, tiled
dispatch, CSR assembly and widened distance/weight output. Neighbor indices,
CSR structure and `distancePairs` match exactly; distances and weights remain
within `1e-6`.

This is a physical backend and numerical contract regression. It does not qualify
cohort-scale timing, million-cell execution, graph biological preservation,
clustering, integration, or biological outcome prediction. The compressed log,
source hashes, machine metadata and command are bound by `execution.json` and
`manifest.json`.
