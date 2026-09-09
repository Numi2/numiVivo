# Complete Norman windowed PCA qualification

This run uses the unchanged native product from commit
`272f125b08d5766d4407c0c748c6bb361101a0b3`. The product and its 66-test qualification
are bound in the [windowed PCA source evidence](../2026-09-09-windowed/source-binding.json).
The native executable SHA-256 is
`71c4fd3c1e9279d98b09def064f3582daf16a246b63c53857a5e034091b1a799`.
The same isolated HDF5 2.2.0 runtime and Apple M4 Pro validation host are used.

Source provenance is inherited from the [complete Norman benchmark](../../../PerturbationPrediction/Norman/README.md).
The original H5AD SHA-256 is
`efde6f5301fe256725dce1d980f37bd96a13481a9a16135515897368e631affc`.
No cells or source genes were discarded for this scale check. `input-plan.json`
is the predeclared request; `plan.json` is the native canonical plan. Numerical
settings and tolerance remain at their original defaults. The work allowance
is explicit and does not change the algorithm or convergence threshold.

`publish.log` and `verify.log` contain the exact native receipts/status plus
`/usr/bin/time -l` measurements. Both commands exited successfully. `checks.json`
records commands, product/source identities, dimensions, sparse work and measured
memory. `report.json.gz` contains every original report byte; its decompressed
SHA-256 matches both the native receipt and independent reference result.
`reference.json` and `reference.log` retain the full Scanpy comparison and timing.
The reference allocation changes were regressed against both complete Baron/Hagai
H5AD inputs and the earlier resident Kang JSON dataset.

Retained artifacts:

- Native: `macmini:/Users/n/numivivo-norman-window-pca-20260909`.
- Native executable: `macmini:/Users/n/numivivo-window-pca-20260909/numivivo`.
- Reference: `/Users/home/numivivo-norman-window-pca-20260909`.
- Original local source: `/Users/home/numivivo-norman-20260909/original.h5ad`.

To reproduce with the verified source, pinned product and explicit HDF5 runtime:

```sh
numivivo singlecell-h5ad-pseudobulk original.h5ad \
  --plan input-plan.json --output new-bundle
numivivo singlecell-h5ad-pseudobulk-verify new-bundle
python3 Tools/Omics/Reduction/check_reference.py \
  --h5ad original.h5ad --report new-bundle/report.json --out new-reference.json
```

This is full-cohort numerical/storage qualification. It does not qualify unseen
perturbation prediction: preprocessing here sees all conditions. Metadata,
aggregates, PCA scores and report serialization remain resident. Neither a
million-cell run nor parallel/Metal execution is claimed. Native and reference
hosts differ, so the timing logs are not a CPU/scverse speed comparison.
