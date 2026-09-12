# Kang source-identity mapping: partial exposure admission

The executed join verifies all eight frozen model donor IDs against sixteen
control/stimulated native sample IDs. Source H5AD barcode and sample order match
the native metadata across all 24,673 cells. Every donor/condition has source-
annotated B cells; the complete B-cell subset totals 2,651 cells. `mapping.json`
records per-condition counts and hashes the source H5AD, native metadata and
model donor list. No expression values or predictions are changed.

These sample IDs derive from the retained pertpy-hosted Kang benchmark annotations;
they are not original GEO/library accessions. The source description identifies
Kang 2018, doi:10.1038/nbt.4042, pertpy file 34464122. Batch is explicitly unreported.
Primary library accessions remain null. The existing primary exposure protocol
must not be treated as completed per-library admission on this evidence alone.

`map.py` reproduces the join using backed AnnData metadata reads and the retained
native JSON files, writing to a fresh directory. This does not reverify donor mean
expression, authenticate all upstream demultiplexing, or establish cross-study
participant independence. Primary library linkage is the remaining Kang step;
GSE181897 mapping is also still open.
