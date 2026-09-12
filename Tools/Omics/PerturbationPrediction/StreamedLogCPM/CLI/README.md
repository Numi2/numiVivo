# Native streaming log-CPM command

```sh
numivivo-omics singlecell-logcpm-stream --plan plan.json \
  --stream-sha256 SHA256 --output new-report.json < records.bin
```

The plan contains featureIDs, groupIDs, rowGroups and full-axis UInt64 rowTotals.
Canonical records are little-endian row-u32, feature-u32, count-u64, strictly
ordered with positive counts. Optional `ordering` is `rowMajor` (default) or
`featureMajor`; the coordinate order must match the declared layout. Plan loading is capped at 256 MiB; the accumulator
retains row assignments/totals and group-by-feature sums, not a cell-by-gene
matrix. Therefore this is bounded streaming, not fully disk-backed metadata.
All source features must contribute to row totals before any panel projection.
Zero-count cells remain in group denominators.

The command requires an expected lowercase stream SHA256, validates all row
totals and the final digest before writing, and refuses an existing output.
The report contains normalized means and axes; it is not a provenance bundle.
Retain the plan, source receipts, expected digest and executing binary separately.

The complete scoped build with --with-cli passes. All values on the real
2,711-cell, 36,601-feature matrix match the prior verified output exactly.
Truncated, duplicate, wrong-total, out-of-axis, empty and wrong-hash streams
produce no report. An existing report is preserved. Zero-cell behavior passes.
The shared consumer also passes oversized-read, read-error and cancellation
checks, and exact real-matrix equality with record-splitting reads.

[receipt.json](receipt.json) binds the binary, sources and real matrix.
[check_cli_actual.py](check_cli_actual.py) and [Errors.swift](Errors.swift)
retain the executed checks and original host paths. First the build was invoked
without --with-cli, so the initial CLI check could not locate the binary; the
corrected isolated build and checks above passed. Existing Metal deprecation
warnings remain. Full-app build, performance and biological accuracy are outside
this qualification.

[Per-chunk buffer cleanup](Memory/README.md) reduces measured CLI peak RSS
from 109.4 MB to 25.6 MB on the same real matrix with identical results.

[Feature-major qualification](../FeatureMajor/README.md) now covers the complete
111,445-cell Norman store and all 237 conditions without count reordering.
