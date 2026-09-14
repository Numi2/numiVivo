# Native Metal sparse PCA operators (2026-09-14)

The source revision `fdcd8f5b` was built and executed on the physical Apple M4
host (`Mac16,12`) with
`NUMIVIVO_TEST_METAL=1 swift test --build-path /tmp/numivivo-metal-pca-operators-build --filter SingleCellMetalPCAOperatorsTests`.
All three focused tests passed.

The resident Krylov PCA owner now has an opt-in `pcaOperatorsBackend =
metalFP32` profile. Projection uses the selected-expression CSR stream and
transpose uses a feature-major CSC stream, retaining only nonzero entries. GPU
outputs are widened to `Double` before the existing eigensolver and residual
checks. The direct test compares both operations with a CPU FP64 sparse oracle;
the integrated test exercises feature statistics plus the Metal projection and
transpose path on the paired synthetic reduction fixture.

This is a bounded physical dispatch and numerical-contract qualification. CPU
FP64 remains the default. The FP32 profile requires a declared residual
tolerance and does not establish broad GPU speedup, million-cell scaling,
model-fitting acceleration, biological-preservation or biological-outcome
prediction. The compressed log, source hashes, host metadata and command are
bound by `execution.json` and `manifest.json`.
