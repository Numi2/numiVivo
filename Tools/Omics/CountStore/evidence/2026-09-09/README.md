# Count-store qualification evidence

This receipt set binds the modified source to the native release product. The
validation checkout started at `adc4b4d4ec0d3b932d08f0e15049cd52835f1598`;
`source-binding.json` records all modified source/test/build-file hashes and
confirms matching local and remote bytes. Documentation changes do not alter
that binary. `native-host.txt` records the actual build host and toolchain.

`initial/` is the first implementation's unchanged output, including its runner
which did not yet gate memory. Reconstruction passed, but the measured memory
failed the later 1 GiB requirement. `initial-memory-assessment.json` preserves
that distinction. `final-interop/` uses the fixed release, while
`initial-interop/` identifies the earlier product. `final-tests.log` and
`final-release.log.gz` retain the selected 65-test gate and full CLI release build.

The large source and record files are retained outside Git:

- Initial native artifacts: `macmini:/Users/n/numivivo-count-store-20260909`.
- Final native artifacts: `macmini:/Users/n/numivivo-count-store-final-20260909`.
- Full source on reference host: `/Users/home/numivivo-norman-20260909/original.h5ad`.
- Independent comparison: `/Users/home/numivivo-count-store-norman-gzip-check-20260909`.

The final root contains `original.h5ad`, `plan.json`, the qualified `numivivo`,
`store/` and `normalized/`. Source provenance and mapping are inherited from the
[complete Norman benchmark](../../../PerturbationPrediction/Norman/README.md).
Count and normalized payloads each contain 361,582,621 16-byte records. Their
hashes are receipt-bound; omitting those 11.57 GB from Git does not imply a
sampled comparison. Original and normalized metadata/statistics are archived
in gzip form where needed; their decompressed bytes retain receipt hashes.

To reproduce, build the full `numivivo` product from the bound source and supply
an HDF5 runtime through `NUMIVIVO_HDF5_LIBRARY`. Create a fresh native output
root with the verified source as `original.h5ad`, `norman-plan.json` as
`plan.json`, and the product as `numivivo`. Then run:

```sh
python3 Tools/Omics/CountStore/run_native.py --root NATIVE_ROOT
python3 Tools/Omics/CountStore/check_norman.py \
  --source LOCAL_SOURCE.h5ad --host macmini \
  --remote-root NATIVE_ROOT --out NEW_REFERENCE_DIRECTORY
```

The reference environment and numerical tolerances are recorded in the
independent result. Timings include source reconstruction and disk writes;
reference and native hosts differ, so no comparative speedup is claimed.

`reference-uncompressed-interrupted.*` retains the deliberately interrupted first
reference transfer. The final checker uses gzip over SSH; numerical comparisons
and hashes consume the complete decompressed payload, with gzip checksum and
stream termination validation. The interrupted run is not a qualification.
