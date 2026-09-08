# Native single-cell counts to a reproducible report

This example is synthetic. It has three features and three barcodes per library,
including an empty cell. The second library has twice the raw counts of the
first but the same library-size-normalized values. Sample and replicate IDs
differ; repeated barcode strings do not refer to the same cells. Duplicate gene
symbols retain distinct IDs. This is not a measured biological dataset.

```sh
swift build -c release
RUN="$(mktemp -d "${TMPDIR:-/tmp}/numivivo-singlecell.XXXXXX")"
.build/release/numivivo singlecell-run Examples/singlecell/manifest.json \
  --store "$RUN/artifacts" --output "$RUN/receipt.json"
.build/release/numivivo singlecell-verify "$RUN/receipt.json" \
  --store "$RUN/artifacts" --output "$RUN/verification.json"
.build/release/numivivo singlecell-export "$RUN/receipt.json" \
  --store "$RUN/artifacts" --output "$RUN/report.json"
.build/release/numivivo singlecell-mex "$RUN/receipt.json" \
  --store "$RUN/artifacts" --output "$RUN/reexported-counts"
```

The report retains exact raw counts, per-cell quality metrics, a separate optional
normalized view, and raw pseudobulk counts with source-cell indices. Expected
cell totals are `5, 8, 0, 10, 16, 0`; pseudobulk rows are `[2, 4, 7]` and
`[4, 8, 14]`. Empty cells remain present in this count report, with missing rather
than fabricated zero mitochondrial fractions.

The receipt binds immutable input bytes, the result and the executing binary/OS.
Verification reconstructs the calculation from archived bytes, not current source
paths. Changed inputs get different identities. After rebuilding, regenerate a
receipt with the new executable; cross-version equivalence is not assumed.

## Supply another count dataset

Use the manifest shape shown here. Declare each library's sample, biological
replicate, optional donor, condition, batch, organism, evidence origin and count
unit. Supply mitochondrial feature IDs explicitly. Optional `cellGroups` maps
barcodes to supplied annotations; no cell type is inferred from counts.

Plain UTF-8 or gzip `matrix.mtx`, `features.tsv` and `barcodes.tsv` sources are
supported by the campaign. Gzip is decoded natively, with checksums and an
aggregate expansion budget. Original compressed bytes remain archived. Matrix
Market must use `coordinate integer general`, with features as rows and barcodes
as columns. Features require ID, name and `Gene Expression`. HDF5/AnnData,
automatic modality filtering and gene-ID remapping are not implemented.
Paths must remain relative to the manifest directory; `..`, symlinks and special
source files are rejected.

All libraries need the exact same ordered feature dictionary and annotations,
count unit and source evidence class. Technical libraries may share a biological-
replicate ID; inconsistent donor/organism metadata for that ID is rejected.
Pseudobulk aggregation preserves contributing donors, samples and batches, but
does not establish that declared replicates are statistically independent.

## Quality filtering and expression analysis

`singlecell-analyze` is a separate stage with an explicit analysis plan. It adds
recorded cell selection, feature summaries and optional donor-aware expression
contrasts without changing the original count result. This small count-only
fixture is not an appropriate differential-expression example: both libraries
are control samples, and it has insufficient replication/reference genes.
Use the [complete paired-donor example](../singlecell-cohort/README.md), generated
by `singlecell-example`, for the analysis commands and TSV/MEX exports.

## Qualification commands

```sh
bash Tools/Omics/run_portable_checks.sh
python3 Tools/Omics/check_cli.py --binary .build/release/numivivo \
  --out "$RUN/cli-checks"
python3 Tools/Omics/check_analysis_cli.py --binary .build/release/numivivo \
  --out "$RUN/analysis-cli-checks"
```

The Swift checks execute native numerical/source-reading code. Python only
orchestrates tests of the public executable and provides no production fit.
The `Single-cell count workflows` action retains source/compiler/binary identities
and execution logs; the existence of that action does not establish a passing run.

See [the count contract](../../Documentation/Design/SINGLECELL_COUNTS.md),
[the analysis contract](../../Documentation/Design/SINGLECELL_ANALYSIS.md), and
[the original scoped audit](../../Documentation/Audit/2026-09-08_SINGLECELL_FOUNDATION.md).
