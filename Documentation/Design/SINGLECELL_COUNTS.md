# Single-cell count-assay foundation

## Scope and ownership

The `Omics` module adds a count-assay payload and native transformations to the
existing scientific platform. It does not create another artifact database,
workflow scheduler, evidence-signature system or GPU simulation runtime.
`VivoMatrixMarketCounts` owns strict MEX import; `VivoSingleCellAnalysis` owns
concatenation, quality metrics, normalization and pseudobulk aggregation;
`VivoSingleCellCampaign` composes those authorities. Both the dedicated CLI and
`vivo.platform.singlecell` use that same campaign.

The source origin (`measured`, `synthetic`, `simulated`) distinguishes assay
inputs. It does not replace the platform's observed/derived/inferred evidence
classification, and a supplied `measured` label does not validate the experiment.
The original manifest and source-file bytes are the lineage record. Result
fingerprints establish content identity, not authorship or biological truth.

## Numerical and identity contracts

Raw integer counts are represented by canonical cell-by-feature CSR arrays:
row offsets, sorted unique feature indices and positive UInt64 counts. Zero
entries are implicit; empty cells remain explicit rows. All row sums and
pseudobulk additions reject UInt64 overflow. Import and raw aggregation never
pass through Double. Two feature IDs may have the same name; duplicate IDs are
rejected. Cells are identified by sample ID plus barcode, not barcode alone.

Quality metrics contain total counts, detected features and explicitly annotated
mitochondrial counts. A missing mitochondrial annotation or zero total yields
no mitochondrial fraction. No name-prefix heuristic substitutes for annotation.

When requested, normalization computes
`log1p((Double(count) / Double(cellTotal)) * targetSum)` in a separate sparse
floating-point view. It retains the original structure and row ordering. This
view is approximate for very large counts and cannot be supplied to an API
requiring raw counts. `targetSum` must be finite and positive.

Pseudobulk uses raw sums by supplied biological-replicate ID, condition and
optional group. Technical libraries may pool within this explicit key. Donor,
batch, sample and source-cell identities remain in the output; distinct
replicates do not become one simply because their donor ID matches. Conversely,
creating distinct replicate IDs does not prove statistical independence. This
module computes neither differential-expression statistics nor causal effects.

## File admission and immutable replay

Manifests reject unknown option names and unsupported schemas. Defaults admit
at most 64 libraries, 100,000 cells, 100,000 features, 2,000,000 nonzeros, 64 MiB of
aggregate source data and 16 KiB per line. Manifest bytes have a separate 2 MiB
limit and also count toward aggregate source bytes. A manifest can reduce, not
increase, the caller's ceiling. Text record limits bound line-array expansion;
the matrix profile allows at most `maximumNonzeros + 4096` lines, including
header, dimensions and comments. These are admission limits, not a measured RSS
bound or a demonstrated maximum useful dataset size.

`VivoSingleCellCampaignIO` uses the existing descriptor-relative
`VivoRootedFileStore`. The original manifest and admitted input bytes are
snapshotted before publication. Paths cannot escape the manifest directory;
symlinks, FIFOs and other nonregular source files are refused. This records the
bytes read; it does not establish that independently changing source files were
captured at a common experimental instant.

`VivoSingleCellArtifacts` stores the complete input bundle and an input- and
implementation-bound result through `VivoArtifactStore`. It publishes no mutable
reference. An interrupted operation can leave immutable orphan objects but cannot
return a success receipt before both objects exist. Receipts can be verified by
hash checks plus full native reconstruction. A correctly hashed but numerically
altered report must still fail. Verification does not read mutable source paths.

The CLI reuses `VivoWorkflowCLIImplementation` to bind the executing binary and
OS. Another binary or OS is rejected rather than silently reinterpreting a
receipt. A future cross-version migration needs a separate explicit comparison
contract. No bitwise cross-platform claim follows from this version.

## General workflow composition

The existing registry exposes `vivo.platform.singlecell`, version `1`, with input
port `input` of kind `vivo.singlecell-input-bundle-v1`, empty configuration and
output port `report` of kind `vivo.singlecell-report-v1`. Use the dedicated CLI's
stored input fingerprint as a `stored` workflow artifact; do not rewrite raw
counts through a generic floating-point JSON representation. The existing
scheduler supplies dependency handling, resource admission, immutable receipts
and validated cache reuse. The adapter supplies deterministic reconstruction of
its output. `Tools/Omics/check_cli.py` builds and executes this recipe rather than
merely checking the catalog entry.

## Explicit remaining work

This is bounded native CPU processing. There are no claimed Metal speedups,
out-of-core matrices, HDF5/AnnData or OME-Zarr adapters, compressed MEX ingestion,
cell filtering, doublet or ambient-RNA correction, donor-aware differential
expression, batch correction, embeddings, inferred annotation or perturbation
prediction. No biological benchmark is completed by the synthetic fixtures.
Those capabilities must extend the shared identities and source lineage rather
than silently change existing count semantics.

Format reference: [10x Genomics MEX output documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-mex-matrices).
