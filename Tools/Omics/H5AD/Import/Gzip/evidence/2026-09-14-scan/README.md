# Gzip-wrapped H5AD scan-path qualification

Date: 2026-09-14
Source commit: `8b5f4be7665b370c77dc186624cabd6cda329fcc`

The isolated native `SingleCellInterchangeTests` suite passed 7/7. A native
`singlecell-h5ad-write` created a 4-cell × 3-feature H5AD, Python added only an
external gzip wrapper, and native `singlecell-h5ad-pseudobulk` published and
verified the wrapped source. The result reported four source cells and four
canonical nonzero records. The bundle's `original.h5ad` bytes match the
compressed input exactly.

This qualifies the shared bounded decode boundary used by scan-based H5AD
workflows (streamed pseudobulk, count-store, PCA, program scoring and reference
queries). The expanded H5AD exists only in a private temporary directory and is
removed after the scan. Axis projection and annotation remain plain
HDF5-readable-source paths and are not covered by this receipt. This is
interoperability evidence, not biological outcome validation.
