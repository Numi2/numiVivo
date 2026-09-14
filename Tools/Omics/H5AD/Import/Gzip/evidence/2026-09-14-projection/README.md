# Gzip-wrapped H5AD axis-projection qualification

Date: 2026-09-14
Source commit: `65b9e80ef536c848164a7c104205b091e5984a95`

The isolated native `SingleCellInterchangeTests` suite passed 7/7. A native
`singlecell-h5ad-write` created a 4-cell x 3-feature H5AD, Python added only an
external gzip wrapper, and native `singlecell-h5ad-project` published and
verified an explicit `[3,1]` cell and `[2,0]` feature projection. The result
reported source shape `[4,3]`, output shape `[2,2]`, and 33 element visits.
The bundle's `original.h5ad` bytes match the compressed input exactly.

This qualifies the shared bounded private decode boundary for complete-object
axis projection. The expanded H5AD exists only in a private temporary directory
and is removed after projection. Annotation remains a plain HDF5-readable source
path because it stages and edits the complete source object. This is
interoperability evidence, not biological outcome validation.
