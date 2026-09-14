# Native Metal sparse feature-statistics qualification (2026-09-14)

The source revision `8e10488d` was built and executed on the physical Apple M4
host (`Mac16,12`) with
`NUMIVIVO_TEST_METAL=1 swift test --build-path /tmp/numivivo-metal-feature-stats-build --filter SingleCellMetalFeatureStatisticsTests`.
Both focused tests passed.

The new opt-in `metalFP32` reduction groups only sparse nonzero
log-normalized values by feature, reduces feature sums and squares on the GPU,
and widens the results before the existing implicit-zero correction and HVG
selection. The test compares every populated feature's count, mean and
nonzero centered moment with the CPU FP64 owner, then runs the reduction owner
through the Metal profile. No cells-by-features dense matrix is allocated.

This is a bounded physical dispatch and numerical-contract qualification. The
FP64 CPU path remains the default. The receipt does not establish general GPU
speedup, million-cell scaling, downstream model-fitting acceleration, biological
signal preservation or biological-outcome prediction. The compressed log,
source hashes, host metadata and command are bound by `execution.json` and
`manifest.json`.
