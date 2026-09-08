# Native single-cell counts to a reproducible report

This example is synthetic. It is not a measured biological dataset and does not
establish single-cell biological accuracy. The second library has twice the
counts of the first, while its library-size-normalized values are the same.
There are three features and three barcodes per library, including an empty cell.
The samples and biological-replicate IDs differ; the repeated barcode strings do
not refer to the same cells. Duplicate gene symbols retain their distinct IDs.

On the Apple machine, from the repository root:

```sh
swift build -c release
RUN="$(mktemp -d "${TMPDIR:-/tmp}/numivivo-singlecell.XXXXXX")"
.build/release/numivivo singlecell-run Examples/singlecell/manifest.json \
  --store "$RUN/artifacts" --output "$RUN/receipt.json"
.build/release/numivivo singlecell-verify "$RUN/receipt.json" \
  --store "$RUN/artifacts" --output "$RUN/verification.json"
.build/release/numivivo singlecell-export "$RUN/receipt.json" \
  --store "$RUN/artifacts" --output "$RUN/report.json"
```

The report contains the exact raw count dataset, per-cell quality metrics, a
separate optional normalized view, and raw pseudobulk counts with source-cell
indices. The expected per-cell totals are `5, 8, 0, 10, 16, 0`; the pseudobulk
count rows are `[2, 4, 7]` and `[4, 8, 14]`. Empty cells are retained, and their
mitochondrial fractions are missing rather than a fabricated zero percent.

`receipt.json` identifies the immutable input bundle, result and executing
implementation. Source bytes are archived in the existing artifact store.
Verification reads those bytes, not the current source paths, and reconstructs
the result. A later change to an input file creates a different input identity.
The same executable and operating-system identity are required for verification;
this is intentionally not a cross-implementation equivalence certificate.

## Supply another count dataset

Create a manifest with the same shape as `manifest.json`. Declare each library's
sample, biological replicate, optional donor, condition, batch, organism, source
origin and count unit. Supply mitochondrial feature IDs explicitly. Optional
`cellGroups` annotations map barcodes to supplied group labels; no cell types are
inferred. These metadata are user declarations, not independently verified facts.

The current importer accepts uncompressed UTF-8 `matrix.mtx`, three-column
`features.tsv` restricted to `Gene Expression`, and `barcodes.tsv`. Matrix Market
must use `coordinate integer general`, with features as rows and barcodes as
columns. There is no gzip decompression, HDF5/AnnData import, automatic modality
filtering, or gene-ID remapping. Source paths must remain inside the manifest
directory; absolute paths, `..`, symlinks and special files are rejected.

All libraries must share the exact ordered feature dictionary, annotation,
count unit and source evidence class. Do not mix measured and synthetic inputs
and then label the combined dataset measured. Technical libraries can share a
biological-replicate ID; inconsistent donor or organism declarations for that ID
are rejected. Pseudobulk rows group by replicate, condition and supplied cell
group. Donor IDs and contributing batches remain visible for later repeated-
measure or batch-aware analysis. Aggregation does not prove independent samples.

Only omit `normalizationTarget` when no normalized view is needed. Raw counts
always remain UInt64. Library-size normalization is not conversion to absolute
molecules, concentration, or a differential-expression result. There is no cell
filtering, doublet removal, batch correction, clustering, gene selection,
annotation inference or hypothesis testing in this profile.

## Qualification commands

```sh
bash Tools/Omics/run_portable_checks.sh
python3 Tools/Omics/check_cli.py --binary .build/release/numivivo \
  --out "$RUN/cli-checks"
```

The first command compiles and executes the actual Foundation-based Swift
implementation and safe filesystem reader. The second tests the built public
CLI, including artifact integrity, failed inputs, stable repeated runs, generic
workflow execution and cache validation. Python is test orchestration only; it
supplies no production calculation. The `Single-cell count workflows` action
retains the source, compiler, binary and execution logs.

See [the method and ownership contract](../../Documentation/Design/SINGLECELL_COUNTS.md)
and [the implementation audit](../../Documentation/Audit/2026-09-08_SINGLECELL_FOUNDATION.md).
