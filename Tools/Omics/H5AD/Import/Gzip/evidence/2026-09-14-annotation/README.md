# Gzip-wrapped H5AD annotation qualification

Date: 2026-09-14
Source commit: `856a94526f4d4508a75e743574c4b09079c42e31`

The isolated native `SingleCellInterchangeTests` suite passed 7/7. A native
`singlecell-h5ad-write` created a 4-cell x 3-feature H5AD, Python added only an
external gzip wrapper, and native `singlecell-h5ad-annotate` published a new
plain H5AD with an explicit nullable `obs/native_score` column. The source
fingerprint in the annotation receipt matches the exact gzip bytes, and HDF5
reopening verifies the values `[0.5, NaN, -0.5, 0]` plus the immutable edit
journal under `uns/numivivo_edits`.

This qualifies the bounded private decode boundary for annotation. The decoded
view is removed after staging; the published output is a new plain H5AD. This is
interoperability and storage evidence, not biological annotation or outcome
validation.
