# Native Metal batched negative-binomial objective (2026-09-14)

The source revision `b976d7a2` was built and executed on the physical Apple M4 host (`Mac16,12`) with
`NUMIVIVO_TEST_METAL=1 swift test --build-path /tmp/numivivo-metal-nb-objective-build --filter SingleCellMetalNegativeBinomialTests`. All four focused tests passed.

The fixed-dispersion NB2 owner now has an opt-in `backend: .metalFP32` profile. It batches the mean-dependent
objective `y*log(mu) - (y + 1/a)*log(1+a*mu)` in bounded Metal tiles. Count-only log-gamma terms, coefficient
updates, diagnostics and the published `logLikelihood` remain on the exact FP64 CPU owner. The direct test
compares the widened GPU sum with an FP64 CPU oracle; the fit test verifies the opt-in record and exact report.

This is a bounded physical dispatch and numerical-contract qualification. CPU FP64 remains the default. The
FP32 count and mean bounds are explicit. The profile does not establish full GPU model-fitting acceleration,
cohort timing, FDR or interval calibration, million-cell scaling, biological preservation or biological-outcome
prediction. The compressed log, source hashes, host metadata and command are bound by `execution.json` and
`manifest.json`.
