# Windowed PCA and fixed-buffer snapshot qualification

The validation checkout starts at `e21e757d35a2da4532d41a73e78d12526a9976ab`.
`source-binding.json` records the changed native sources, tests and scoped source
list, with exact matching local/remote bytes. Full native work runs on
`macmini:/Users/n/numivivo-h5ad-20260909` and the prepared experimental cohorts
live under `macmini:/Users/n/numivivo-window-pca-20260909`.

The first test build failed because Swift could not type-check a compound test
expression in time. `initial-tests.log.gz` preserves that failure. The
intermediate window-only change passed 66 tests and built successfully; its
logs are retained separately. `tests.log.gz` is the final 66-test gate after
extracting the shared snapshot helper. Compressed logs preserve original bytes.

The final native runner compares both complete source cohorts against the
published whole-map product, then exercises reconstruction, exact repeats,
scratch cleanup and budget/convergence rejection through `check_streamed_cli.py`.
The comparison checks all PCA fields exactly, as well as metadata, QC, aggregate
counts, cache digest and sparse arithmetic work. The native baseline executable
is the retained count-store release at
`/Users/n/numivivo-count-store-final-20260909/numivivo`.

Input provenance and original real-data qualification remain in
[STREAMING.md](../../STREAMING.md) and the linked benchmark suite. Reproduce with
prepared `baron-input/` and `hagai-input/` directories containing `prepared.h5ad`
and `mapping.json`:

```sh
python3 Tools/Omics/Reduction/check_windowed_cli.py \
  --binary /path/to/new/numivivo --baseline /path/to/published/numivivo \
  --root /new/native-root
python3 Tools/Omics/Reduction/check_reference.py \
  --h5ad /path/to/prepared.h5ad --report /native-root/baron/bundle/report.json \
  --out /new/baron-reference.json
```

Both exact native regression and independent Scanpy comparison are required.
Same-host timing observations are single runs and do not establish general
speedup. Windowed mapping does not remove resident metadata, moments, Krylov
basis, scores, aggregates or reports. This is not million-cell, parallel-kernel,
GPU or biological-preservation qualification.
