# Single-cell count-assay foundation

## Scope and ownership

The `Omics` module adds a count-assay payload and native transformations to the
existing scientific platform. It does not create another artifact database,
workflow scheduler, evidence-signature system or GPU simulation runtime.
`VivoOmicsSourceDecoder` owns bounded plain/gzip source decoding;
`VivoMatrixMarketCounts` owns strict MEX parsing; `VivoSingleCellAnalysis` owns
concatenation, quality metrics, normalization and pseudobulk aggregation;
`VivoSingleCellCampaign` composes these authorities. The dedicated count CLI and
`vivo.platform.singlecell` use that same campaign.

Explicit cell filtering and donor-aware expression are a separate, implemented
analysis stage described in [SINGLECELL_ANALYSIS.md](SINGLECELL_ANALYSIS.md).
Count ingestion does not silently run that analysis or delete low-quality cells.

The source origin (`measured`, `synthetic`, `simulated`) distinguishes assay
inputs. It does not replace the platform's observed/derived/inferred evidence
classification, and a supplied `measured` label does not validate the experiment.
The original manifest and source-file bytes, including compressed originals,
are the lineage record. Fingerprints establish identity, not biological truth.

## Numerical and identity contracts

Raw integer counts use canonical cell-by-feature CSR arrays: row offsets,
sorted unique feature indices and positive UInt64 counts. Zero entries are
implicit; empty cells remain explicit rows. Row sums and pseudobulk additions
reject UInt64 overflow. Import and raw aggregation never pass through Double.
Two feature IDs may have the same name; duplicate IDs are rejected. Cells are
identified by sample ID plus barcode, not barcode alone.

Quality metrics contain total counts, detected features and explicitly annotated
mitochondrial counts. Missing mitochondrial annotation or zero total yields no
mitochondrial fraction. No gene-name prefix heuristic substitutes for annotation.
When requested, normalization computes
`log1p((Double(count) / Double(cellTotal)) * targetSum)` in a separate sparse
floating-point view. It retains row/feature ordering and cannot substitute for
raw counts. `targetSum` must be finite and positive.

Pseudobulk uses raw sums by supplied biological-replicate ID, condition and
optional group. Technical libraries may pool within this explicit key. Donor,
batch, sample and source-cell identities remain in the output. Distinct replicate
IDs do not by themselves prove independence. The separate analysis stage applies
additional replication/design checks before statistical inference.

## File admission and immutable replay

Manifests reject unknown option names and unsupported schemas. Defaults admit
at most 64 libraries, 100,000 cells, 100,000 features, 2,000,000 nonzeros, and
16 KiB per line. Original bytes and expanded bytes each have an aggregate 64 MiB
allowance across the entire campaign. Manifest bytes have a separate 2 MiB limit
and count toward both aggregate allowances. Limits may be reduced, not silently
raised. The matrix profile permits at most `maximumNonzeros + 4096` lines,
including header, dimensions and comments. These are admission limits, not a
measured RSS bound or a demonstrated maximum useful dataset size.

Gzip is detected by magic bytes, not the filename. The native system-zlib decoder
validates checksums/trailers, supports bounded concatenated members, and rejects
truncation, trailing non-gzip data and excessive expansion. The strict MEX parser
still requires coordinate/integer/general matrices and Gene Expression features.
There is no automatic modality filtering or gene-ID reconciliation.

`VivoSingleCellCampaignIO` uses the existing descriptor-relative
`VivoRootedFileStore`. Sources are snapshotted before publication. Relative paths
cannot use `..`; symlinks, FIFOs and other nonregular source files are refused.
This records the bytes read, not an experimentally simultaneous snapshot of
independently changing files.

`VivoSingleCellArtifacts` stores the complete input bundle and an input- and
implementation-bound result through `VivoArtifactStore`. It changes no mutable
reference. Interrupted publication may leave immutable orphan objects but cannot
return a success receipt before both objects exist. Verification checks hashes
and reconstructs the native result from archived input, not mutable source paths.
A correctly hashed but numerically altered result still fails reconstruction.

The CLI reuses `VivoWorkflowCLIImplementation` to bind the executable and OS.
After rebuilding, regenerate count receipts with that executable before new
analysis. Cross-version equivalence is not inferred by relabelling old receipts.

## Workflow composition and interchange

`vivo.platform.singlecell`, version `1`, accepts `input` of kind
`vivo.singlecell-input-bundle-v1`, empty configuration and produces `report` of
kind `vivo.singlecell-report-v1`. Use the stored input fingerprint as a `stored`
workflow artifact rather than rewriting raw counts through generic JSON numbers.
The existing scheduler retains dependencies, admission, immutable receipts and
validated cache reuse. `vivo.platform.singlecell-analyze` adds the separate
analysis-plan input and shares the native cohort analysis authority.

`singlecell-mex` exports verified raw counts and metadata to a new directory.
Per-sample matrices retain count units, evidence class, ordered features, cell
groups and mitochondrial annotations. Exported-to-source row maps preserve
lineage when source rows were interleaved. The cohort commands also export
selected MEX libraries and TSV quality, feature, design and contrast tables.

## Remaining boundaries

This is bounded native CPU processing. No Metal speedup, out-of-core matrix,
HDF5/AnnData or OME-Zarr adapter, doublet/ambient correction, inferred annotation,
clustering or perturbation-prediction result is claimed. The implemented
log-expression inference is not a calibrated replacement for established
negative-binomial or weighted-expression pipelines. Synthetic fixtures establish
neither measured assay validity nor biological false-discovery calibration.

Format reference: [10x Genomics MEX output documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/cr-outputs-mex-matrices).
