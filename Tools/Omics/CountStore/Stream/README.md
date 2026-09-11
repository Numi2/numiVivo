# Native canonical count streams

`singlecell-count-stream-pseudobulk` consumes raw records from standard input,
checks every record, and sums them by biological replicate, condition and optional
cell group. It supports remote source adapters without requiring a local copy of
their full count matrix. Source extraction remains the adapter's responsibility.

```sh
Tools/Omics/H5AD/build.sh /path/to/build --with-cli
/path/to/build/numivivo-omics singlecell-count-stream-pseudobulk \
  --plan axes.json --output new-bundle < records.bin
/path/to/build/numivivo-omics singlecell-count-stream-verify \
  new-bundle < records.bin
Tools/Omics/CountStore/Stream/test.sh /path/to/build
```

Each 16-byte record contains little-endian `UInt32` row index, `UInt32` feature
index and `UInt64` positive raw count. Rows increase; features are strictly
increasing within each row. Duplicate records, zero counts, invalid axes,
truncation, trailing bytes and integer overflow are errors. Empty cells remain
in the metadata and QC output. The producer must emit the complete declared
stream, close it, and check both its own and the native process's exit status.

The schema-1 plan contains `metadata` (`VivoSingleCellCountMetadata`),
`sourceDeclaration`, `rowNonzeros` and optional `rowTotals`. Cardinalities must
match every row. Supplied totals must describe the exact retained feature axis;
historical QC values from a different matrix must remain separate. Without
independent totals, native totals are computed from the records and should be
checked against an independent reader. No normalization or inferred annotation
is applied. Missing mitochondrial annotation leaves its fraction unavailable.

Input reads are at most 1 MiB. A read chunk and a parser buffer of up to 1 MiB
plus 15 bytes can coexist. A pool drains Foundation temporaries each iteration;
[paired allocator controls](Memory/README.md) verify this correction. Metadata,
per-cell QC and donor aggregates remain resident: this is **not** fully out-of-core cell
metadata. Bounds are two million cells, 100,000 features, four billion records,
five million aggregate nonzeros, a 256 MiB plan and a 512 MiB report. Serialized
size can impose a tighter bound. No Metal acceleration is claimed for this owner.

A new bundle contains the exact canonical plan, report and receipt. Publication
is transactional. The receipt binds the implementation, plan, report, stream
byte count and SHA-256 of every consumed record. Verification requires a complete
replay and recomputes QC and aggregates. Updating a forged report's receipt hash
does not bypass reconstruction. The bundle does not retain the full stream;
keep a reconstructible source or a separate stream archive. The source declaration
and an HTTP ETag are **not** a SHA-256 of unread upstream HDF5 bytes, and native
record validation is distinct from HDF5 decoding and biological validation.

The scoped tests compare fragmented streams with the resident count owners,
exercise corrupt and truncated inputs, validate row contracts and overflow,
and check transactional publication and forged-report rejection. These are
numerical/software controls. The [Parse IFN-beta admission](../../PerturbationPrediction/ParseIFNB/README.md)
tracks the separate complete experimental count execution.
