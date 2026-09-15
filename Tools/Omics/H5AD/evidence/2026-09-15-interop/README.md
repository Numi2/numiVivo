# Native AnnData interoperability run (2026-09-15)

The current repository revision passed the actual `Tools/Omics/H5AD/check_interop.py`
harness against a freshly built `h5ad-check` binary. The harness creates AnnData
fixtures with the pinned Python stack, exercises CSR, CSC and dense count layers,
then reopens native output in AnnData. It also checks exact UInt64 counts,
sparse canonicalization, raw's independent feature axis, categorical/nullable
metadata, explicit feature identity, malformed H5AD rejection and native
reimport.

`report.json` is the compact result receipt. `sources.sha256` binds the native
source set used by the scoped build. The generated H5AD fixtures are intentionally
not retained here because they are reproducible and the report is the durable
claim. This is interoperability evidence only; it is not a biological benchmark.

Reproduce from the repository root:

```sh
bash Tools/Omics/H5AD/build.sh /tmp/numivivo-h5ad-build-interop
python Tools/Omics/H5AD/check_interop.py \
  --binary /tmp/numivivo-h5ad-build-interop/h5ad-check \
  --out /tmp/numivivo-h5ad-interop-report
```

Pinned Python versions: AnnData 0.13.3.post0, h5py 3.16.0, NumPy 2.5.3,
pandas 3.0.5 and SciPy 1.18.1. Native HDF5: 2.2.0. Binary SHA-256:
`90aa101944410eaf1826e1a84d5cb0fa2e159f18f23e4dde3e1d6b07797caa8a`.
