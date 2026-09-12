# AnnData interoperability and the single-cell development sequence

NumiVivo reads local H5AD through the native HDF5 library; those import/export
and analysis commands do not require Python. The optional remote-source adapter
described below uses Python/HDF5 to produce records for a separate native stream
reader, and independent interoperability tests also use Python.
Install HDF5 (`brew install hdf5` on macOS), or set `NUMIVIVO_HDF5_LIBRARY` to
its shared library. Other workflows do not require HDF5. HDF5 calls are serialized
because installed libraries may not be thread safe.

## Current outcome and scale evidence

The [biological prediction assessment](BiologicalPrediction.md) states what can
be predicted from available data, including negative results and simple baselines.
The complete [HIRISA experiment](../Tools/Omics/Benchmarks/HIRISA/README.md) adds
1,612,594 cells, 131 libraries, 48 native/reference DE cases, 79 frozen donor
prediction folds, and complete PCA/graph qualification. Native seed-7 integration
and three Harmony references pass coarse preservation margins. The subsequent
[program diagnostic](../Tools/Omics/Benchmarks/HIRISA/INTEGRATION_PROGRAMS.md) fails
three sensitive comparisons in all four candidates, with 18/32 controls
insufficient. A separately declared [within-library fitting follow-up](../Tools/Omics/Benchmarks/HIRISA/PROGRAM_CALIBRATION.md)
meets the same margins with 28/32 sensitive controls, retaining four insufficient
controls and the original failures as development evidence. Full-cohort [native clustering publication/replay and independent
checks](../Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md) now pass;
broader biological preservation remains open.

## Remote count records

The [native count-stream reader](../Tools/Omics/CountStore/Stream/README.md) accepts
canonical little-endian row/feature/count records over standard input. It validates
strict ordering, positive counts, row cardinalities, optional exact row totals,
truncation and overflow. It computes cell QC and donor/condition/group aggregates,
binds every input byte by SHA-256, and requires stream replay to verify its bundle.
A source adapter must separately establish HDF5 decoding and source provenance.

The complete [Parse IFN-beta admission](../Tools/Omics/PerturbationPrediction/ParseIFNB/README.md)
passes on 725,031 cells and 1.37 billion records using version-pinned
HTTP ranges. This avoids a 227 GB local source snapshot. The matrix buffer is
bounded; metadata and QC are resident, and unread source bytes are not assigned
a fabricated whole-file fingerprint. No acceleration or biological-prediction
claim follows from count-stream validation.

A [temporary-memory correction](../Tools/Omics/CountStore/Stream/Memory/README.md)
keeps the native allocator control near 105 MB across 160 MB and 1.6 GB streams;
the original implementation grew to 1.71 GB. Exact counts/reports and four native
tests pass. These synthetic controls qualify the buffer-lifetime repair; the
complete Parse ingestion and all twelve native source replays now pass;
[offline artifact restoration](../Tools/Omics/PerturbationPrediction/ParseIFNB/COUNT_RESULTS.md) also passes. The [preparation binding](../Tools/Omics/PerturbationPrediction/ParseIFNB/README.md#immutable-preparation-binding)
also checks every donor projection against the immutable original source axes.

The [file-backed cell-axis/count path](../Tools/Omics/CountStore/CellAxis/README.md)
now implements disk-resident cell identities, QC and group membership with bounded
resident dictionaries and aggregate counts. Fourteen native tests pass. Complete
Parse cell-axis import and reopen preserve all 725,031 original identities, row
cardinalities and retained-matrix totals, with native peak RSS of 65,617,920 bytes
(62.6 MiB). [Complete paired count ingestion and source replay now pass](../Tools/Omics/CountStore/CellAxis/COUNT_RESULTS.md)
for every selected record. File-backed peak RSS is 165.6/181.8 MiB during
ingestion/replay, versus 1,005.0/1,024.7 MiB for the resident owner.
The [file-backed expression owner](../Tools/Omics/CountStore/Expression/README.md)
now consumes these aggregates using explicit source-receipt membership references.
On all 725,031 Parse cells, both native fits and reconstruction pass and all six
edgeR/limma/DESeq2 comparisons complete. The statistical result matches the prior
native owner exactly apart from membership representation. The [incremental report writer](../Tools/Omics/CountStore/Expression/Memory/README.md)
now reduces full-cohort native analysis peak RSS from 970.4 to 202.3 MiB with
all report bytes unchanged. Per-feature model state remains resident; this does
not qualify every pipeline stage or add an unseen biological prediction.
The import measurement excludes the source adapter and is not a whole-pipeline
memory or Metal performance result. Existing H5AD/PCA/graph/model bounds remain.

## Use

```
numivivo singlecell-h5ad-import experiment.h5ad --plan mapping.json --output imported
numivivo singlecell-run imported/manifest.json --store artifacts --output receipt.json
numivivo singlecell-analyze receipt.json --plan analysis.json --store artifacts --output analysis-receipt.json
numivivo singlecell-h5ad-write imported/dataset.json --output counts.h5ad
numivivo singlecell-h5ad-annotate imported/original.h5ad --plan annotations.json --output annotated.h5ad
numivivo singlecell-h5ad-project experiment.h5ad --plan projection.json --output projected
numivivo singlecell-h5ad-project-verify projected
```

[Axis projection](../Tools/Omics/H5AD/Projection/README.md) selects/reorders cells
and features across supported aligned AnnData slots, preserving raw's independent
feature axis and copying unstructured data. Its source-bound plan is distinct
from the count-import mapping below; it does not infer which matrix is raw counts.

The import mapping explicitly names the count array and experimental design.
AnnData does not standardize which layer contains raw counts or which columns
identify biological replicates. For example:

```json
{
  "schemaVersion": 1,
  "id": "experiment",
  "evidence": "measured",
  "sourceDescription": "Dataset accession, release and assay description",
  "countUnit": "umiCount",
  "matrixPath": "layers/counts",
  "sampleColumn": "sample",
  "groupColumn": "cell_type",
  "featureNameColumn": "gene_symbol",
  "mitochondrialFeatureIDs": [],
  "samples": [
    {
      "id": "donor1-control",
      "biologicalReplicateID": "donor1",
      "donorID": "donor1",
      "condition": "control",
      "batchID": "batch1",
      "organism": "NCBITaxon:9606"
    }
  ]
}
```

Supply every actual sample, with the study's actual replicate structure.
`barcodeColumn` is optional; otherwise the observation index supplies barcodes.
`groupColumn`, `featureIDColumn` and `featureNameColumn` are optional.
Use `featureIDColumn` for stable gene IDs when the AnnData index contains
nonunique gene symbols; the index remains the default display name. Missing
group annotations are preserved; missing selected design identities are rejected. Mitochondrial feature
identities are explicit, never inferred from names.

Preserve complete barcode identities when joining metadata. The
[Adamson GEO restoration](../Tools/Omics/PerturbationPrediction/Adamson/GEO_RESTORATION.md)
corrects 3,634 assignments caused by removing barcode suffixes before a metadata
join. Successful H5AD parsing and count reconstruction alone cannot establish
that deposited experimental labels are correct.

## Preservation contract

- Import publishes the **unchanged original.h5ad**, a native `dataset.json`, the
  mapping, a source/dataset/mapping fingerprint receipt with executable identity
  and HDF5 version, and a MEX manifest usable by the existing replayable workflow.
- The count projection reads X, raw/X or a named layer. Raw import uses
  raw/var and its independent feature axis, not the current var dictionary.
  It supports CSR, CSC, and dense
  arrays read through one-row HDF5 hyperslabs. It never allocates cells × genes.
  Sparse duplicates are summed with checked UInt64 arithmetic, indices sorted,
  and explicit zero entries removed. The original representation stays in the
  source file.
- Integer counts retain all UInt64 bits. Floating counts must be finite,
  nonnegative, integral, and no larger than 2^53. Selecting normalized fractional
  data fails instead of converting it to counts.
- Selected metadata columns support strings, categorical strings, and nullable
  strings, including unsigned categorical codes. Root, dataframe, and array
  encoding versions are checked against the
  [AnnData on-disk specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
- Source raw, additional layers, nullable/numeric metadata, category order,
  embeddings, graphs and uns are retained in the original file. Exporting that
  source returns the original object. The separate native annotation command
  adds or explicitly replaces derived fields in a new copy, preserving untouched
  datasets, dtypes, categories and nullable metadata.
- Writing `dataset.json` creates a **new count AnnData object**, with CSR X,
  observation design columns, nullable donor/group columns, feature names,
  mitochondrial annotations and native metadata in uns. Unique generated obs
  indices coexist with original barcode/sample columns. Use `barcodeColumn:
  "barcode"`, `sampleColumn: "sample"`, `groupColumn: "group"`, and
  `featureNameColumn: "name"` when importing native output again.
- Publication refuses existing destinations. Resident count import uses a
  64 MiB source limit; resident count processing is bounded at
  100,000 cells/features and 5 million nonzeros by default. Streaming processing
  and axis projection have separate routes and limits. HDF5 itself may retain
  decompression and variable-string buffers beyond the Swift sparse arrays.

## Editing an existing AnnData object

`singlecell-h5ad-annotate` uses a versioned `VivoH5ADAnnotationPlan` with the
SHA-256 of the exact source, required provenance text, and explicit `add` or
`replace` edits. It keeps the original cell/feature order. It cannot change X,
raw, observation/feature indices or delete fields. Each output embeds an immutable
plan/source/implementation record in `uns/numivivo_edits`; stdout includes the
output fingerprint. Neither free-form labels nor the provenance text establish
biological authority.

Supported destinations are columns in obs/var; numeric embeddings in obsm/varm;
sparse matrices in obsp/varp/layers; and typed values or nested dictionaries in
uns. Supported values include Float64, exact Int64/UInt64, booleans, UTF-8 strings,
ordered/unordered string categories with missing codes, nullable columns, and
CSR matrices. Missing floating results are written as ordinary NumPy NaNs; other
nullable types retain explicit masks. Supplied nonfinite floating values fail.

Source fingerprints prevent applying a result to another ordering or dataset.
Column and matrix dimensions must match the source axes. Embeddings have at most
256 dense components; graph/layer edits must be sparse. Edits are bounded to 256
fields, 2 million payload elements, 64 MiB of string content and 16 nesting levels.
Annotation edits require HDF5 1.12 or newer. Storage metadata is inspected
without reading matrix values; soft/external links, virtual/external datasets and
object references are rejected because their meaning may change during copying.
Existing destinations and conflicting add/replace modes fail. Failures and
cancellation before publication leave no output. Edited parent groups are copied
before mutation so hard-linked backups retain their old values.

Annotation accepts up to 16 GiB source files, 32 GiB outputs and two million
rows/features; the existing plan payload limits still apply. Snapshot and
publication use APFS copy-on-write when available, with bounded 1 MiB streaming
fallback. Publication keeps pinned destination directories, private temporaries
and atomic no-overwrite linking. Edited groups copy attributes and retain links
to untouched children, preserving aliases without copying all their datasets.
HDF5 metadata/decompression allocations remain separate from the copy buffer.
The [full real HIRISA annotation run](../Tools/Omics/H5AD/Annotation/Atlas/README.md)
adds exact source-row provenance to all 1,612,594 cells and verifies all 46 original
datasets (7,735,703,364 stored elements), their attributes and the unchanged source
hash. The complete annotated file also reopens in backed AnnData with both indices,
all 29 original obs/var columns and every added row value verified; sparse count
windows match while the matrix stays backed. Metadata remain resident. This is
annotation/storage qualification, not biological labels or predictions. The complete
[Adamson UPR ingestion check](../Tools/Omics/PerturbationPrediction/Adamson/README.md)
preserves the original 65,337 cells and 237,812,947 sparse entries and verifies
native aggregates against an independent SciPy reference.

A complete executable authoring/verification example is provided by
`Tools/Omics/H5AD/check_annotations.py`. Its JSON value syntax follows the typed
Swift enum, for example:

```json
{
  "path": "obs/native_score",
  "mode": "add",
  "value": {"nullableFloat64": {"values": [0.5, null, -0.5, 0]}}
}
```

Place edits in a plan with `schemaVersion: 1`, the source fingerprint's existing
`{"bytes": [32 byte values]}` representation, and `provenance`. The example above
requires exactly four source observations; it is not a command to annotate an
arbitrary dataset. Annotation execution itself is native Swift/HDF5, with no
Python subprocess.

## Public-file check

`Tools/Omics/H5AD/check_public_data.py` downloads the full public
[PBMC3K dataset used by Scanpy](https://scanpy.readthedocs.io/en/latest/generated/scanpy.datasets.pbmc3k.html)
from its pinned URL and verifies SHA-256 before use. It retains the download,
re-encodes the legacy H5AD through current AnnData, adds provenance annotations
natively and checks all original datasets/dtypes/attributes in the result. No
cells or features are dropped: 2,700 cells, 32,738 features and 2,286,884 nonzeros.
The required reference re-encoding is recorded explicitly; direct native legacy
import is not qualified.

With `--full-product --count-analysis`, this check also runs native count import,
count receipts, analysis and replay/export, comparing exact counts and QC to
Scanpy and log normalization to a 1e-10 absolute tolerance. The complete PBMC3K
matrix passed after raising the bounded count default to 5 million nonzeros.
This is not an out-of-core implementation, and the full replay/export peaked at
1.48 GB resident memory in the recorded debug build. PBMC3K supplies no
donor-aware DE or integration result.

The [experimental benchmark suite](../Tools/Omics/Benchmarks/README.md) records
paired-donor Kang and original-count Hagai NB comparisons, Haber count/QC
checks, full-source Baron streaming QC/pseudobulks, and rejected normalized Hagai
and fractional Muraro inputs. The broader benchmark and
cross-study calibration requirements remain open.

```
python Tools/Omics/H5AD/check_annotations.py --binary /path/to/numivivo --full-product --out /tmp/annotation-checks
python Tools/Omics/H5AD/check_public_data.py --binary /path/to/numivivo --full-product --out /tmp/public-file-check
```

## Streaming raw counts into pseudobulks

```
numivivo singlecell-h5ad-pseudobulk experiment.h5ad --plan stream-plan.json --output streamed
numivivo singlecell-h5ad-pseudobulk-verify streamed
```

The plan has `schemaVersion: 1`, `mapping` containing the H5AD import mapping
above, and `contrasts` containing the existing expression contrast definitions
(an empty list requests QC and aggregation only). The command retains every
source cell. It does not silently apply a resident-analysis filter. The native
NB cohort evaluator consumes the resulting pseudobulks with the original
sample, donor, cell and feature metadata.

The shared native reader now reads sparse indices/counts in slices of at most
65,536 elements. It canonicalizes duplicate coordinates within each CSR row or
CSC column and checks every UInt64 accumulation. Dense inputs use one row at
a time. The streaming route retains metadata, per-cell QC and grouped count
sums; it does not construct the complete resident sparse cell matrix. It copies
and hashes the H5AD source in 1 MiB blocks before analysis.

The [independent-donor null benchmark](../Tools/Omics/NegativeBinomial/IndependentNull/README.md)
now verifies the complete 1,199,390,735-byte Human Immune Health Atlas B/Plasma
source: all 160,632 cells, 233,853,565 raw entries and 108 donor aggregates match
independent count/QC calculations. Its old 1 GiB rejection exposed the source
admission mismatch; streamed PCA reference/query snapshots now share the same
64 GiB copy ceiling. That ceiling is capacity admission, not demonstrated scale.

The [complete HIRISA release](../Tools/Omics/Benchmarks/HIRISA/README.md) now
passes native publication/reconstruction on all 1,612,594 cells and 3.846 billion
entries. An independent streaming audit verifies every cell identity/QC value,
all 131 library aggregates and all sixteen inference-cohort memberships/counts.
Release publication/reconstruction take 280.00/285.23 seconds, peaking at
4.43/4.96 GB RSS. The 526.5 MB JSON report nearly reaches its 512 MiB limit;
this qualifies complete-source ingestion, not general million-cell PCA, graph,
integration or out-of-core execution. All 48 native and 48 reference DE runs are
complete; every native case passes independent conditional numerical checks.
The full inference archive retains each method's gene family, diagnostics and
warnings. All 79 frozen held-out prediction folds through the new aggregate
batch API pass native replay and independent NumPy reconstruction. Complete
results retain both gene families and weaker cases; ridge improves contrast
RMSE over the training-mean baseline in only four of sixteen contrasts. The
[subsequent preparation-transfer experiment](../Tools/Omics/Benchmarks/HIRISA/CONTEXT_TRANSFER.md)
now completes all 120 folds with native replay and independent reconstruction.
Full-cohort PCA, graph, clustering and seed-7 integration also have published
results; these numerical gates do not establish general biological prediction.

The output is an exchange bundle containing `original.h5ad`, `plan.json`,
`report.json` and `receipt.json`. Verification copies the source to a private
snapshot, checks source/plan/report/executable fingerprints and reconstructs
the report. Existing destinations are refused. This route uses the common
fingerprint types but does not yet publish its report into the artifact-store
DAG. It is not a substitute for general chunked transforms or sparse PCA.

Explicit limits are 64 GiB source bytes, 2 million cell identities, 100,000
features, 4 billion source sparse entries, 5 million aggregate nonzeros,
2 MiB encoded plan and 512 MiB encoded report. The pseudobulk CLI now uses
the same 2 MiB input-plan allowance; a complete Adamson mapping with 1,106
guide/GEM sample identities exposed the former 128 KiB CLI mismatch. These
are admission bounds; the complete HIRISA input demonstrates the scope above.
Metadata, group membership and JSON reports remain resident. APFS snapshot
cloning preserves separate content ownership and avoids allocating another full
source payload; other filesystems use bounded descriptor copying. Compression
depends on the installed HDF5 filters; gzip is qualified, while LZF was rejected
on the tested host because no native LZF filter was installed.

The [full Hagai streaming comparison](../Tools/Omics/Benchmarks/README.md#full-hagai-streaming-comparison)
retains all 32.85 million source nonzeros. The existing five-million resident
matrix limit is unchanged.

## Validation commands

```
bash Tools/Omics/H5AD/build.sh /tmp/numivivo-h5ad-build
python Tools/Omics/H5AD/check_interop.py --binary /tmp/numivivo-h5ad-build/h5ad-check --out /tmp/h5ad-checks
```

Use the versions pinned in `Tools/Omics/H5AD/requirements.txt`. The optional
public count/Scanpy check uses `Tools/Omics/Benchmarks/requirements.txt`. The test creates
AnnData files, imports/exports with the actual Swift owners, reopens results in
AnnData, checks exact counts, source bytes, missing annotations and metadata,
and checks malformed input rejection. Reports identify the executable hash and
reference versions. These are software interoperability fixtures, **not evidence
of biological accuracy or competitiveness**.

Verification on 2026-09-09 used AnnData 0.13.3.post0, h5py 3.16.0 and native
HDF5 2.2.0. All 24 interoperability cases and six full-product H5AD command
checks passed. The Mac mini package build and nine existing single-cell tests
passed; 120 portable checks, statistics/compression smoke checks, and the
existing cohort CLI's 13 assertions across 20 commands also passed. These
checks cover count exchange and existing workflow regression, not a whole-suite
scientific qualification. The test harness requires controlled error exits for
invalid inputs; crashes do not count as successful rejection.

The follow-up qualification passed 27 count-interoperability cases and 26
annotation cases, including controlled rejection of storage dependencies and
reference types. The six count CLI checks and nine existing single-cell tests
also passed. The full `numivivo` executable passed the annotation suite and
full-file PBMC3K preservation check with AnnData 0.13.3.post0, h5py 3.16.0 and
HDF5 2.2.0. This adds real-file preservation evidence; it does not satisfy the
multi-donor experimental benchmark requirement.

## Required development order and remaining evidence

The complete development objective remains open:

1. **AnnData/H5AD:** native count exchange, raw-axis import, explicit feature IDs
   and source-preserving annotation edits implemented. Native
   [axis projection](../Tools/Omics/H5AD/Projection/README.md) now filters/reorders
   unique cell/feature indices across supported aligned AnnData slots, retains
   raw's independent feature axis, and binds source/output/replay provenance.
   Repeated-axis checks now cover complete prepared Kang and Baron, all aligned
   slots and raw, with exact prior unique-selection bytes. Native legacy
   projection now preserves the complete original Kang file, including compound
   obs/var, structured PCA/UMAP and categories previously held in uns. Both full
   and repeated selections agree with AnnData and reconstruct exactly. Subsequent
   [numeric-category migration](../Tools/Omics/H5AD/Projection/README.md#numeric-legacy-categories-2026-09-12),
   [scalar embeddings](../Tools/Omics/H5AD/Projection/README.md#scalar-legacy-embeddings-2026-09-12)
   and [numeric indices](../Tools/Omics/H5AD/Projection/README.md#numeric-dataframe-indices-2026-09-12)
   now pass their declared native/AnnData conformance checks. The
   [explicit identity route](../Tools/Omics/H5AD/Projection/NUMERIC_INDEX_ROUTE.md)
   connects numeric-index files to annotation, count-store and pseudobulk without
   inferring biological IDs or silently coercing index values. Ragged embedding
   fields and arbitrary nested records remain unsupported; storage/work bounds
   still apply. The [complete original Kang analytical route](../Tools/Omics/H5AD/Projection/LEGACY_COUNT_ROUTE.md)
   now passes native annotation, streaming count import and pseudobulk:
   14,184,532 exact count records and 124 exact donor/condition/cell-type groups.
   A missing legacy uns dictionary tag was repaired as projection v2; the earlier
   failed annotation remains retained. Full PBMC3K count/QC
   and annotation preservation now pass; this single library is not a donor-DE benchmark.
   [Full-atlas annotation](../Tools/Omics/H5AD/Annotation/Atlas/README.md) now adds
   an exact source-row column to all 1,612,594 HIRISA observations and verifies
   every original HDF5 dataset element. The later
   [graph-community export](../Tools/Omics/H5AD/Annotation/Atlas/Clustering/README.md)
   binds all cell/sample identities and labels and reopens in backed AnnData;
   that export's count check is sampled, while all original metadata are checked.
   Annotation admits 16 GiB source/32 GiB output and two-million-element axes;
   existing total edit and string-byte budgets still apply. These are qualified
   interchange routes, not cell-type or outcome validation.
2. **Experimental benchmarks:** one eight-donor Kang B-cell contrast now passes
   exact Scanpy QC/pseudobulk checks and a descriptive PyDESeq2 comparison; Haber
   tuft-cell count/QC passes. A third treatment study,
   [Crowell cortex](../Tools/Omics/Benchmarks/Crowell/README.md), now adds all
   25,224 deposited nuclei from eight independent mice: seven population
   analyses/replays and 42 direct R fits. Full count/QC and offset checks pass;
   CPE correctly remains unavailable for insufficient animal replication.
   Eight strict conditional-refit p-value checks remain failed (maximum 1.84e-7);
   a stopping-tolerance diagnostic verifies stored Wald arithmetic and finds no
   changed BH<0.05 decisions. Large reference significance disagreements remain,
   especially excitatory neurons. Three donor-resolved treatment studies are
   now represented; broader reference sensitivity and calibration remain. Direct
   [R edgeR/limma-voom/DESeq2 comparisons](../Tools/Omics/Bioconductor/README.md)
   now cover the same Kang/Hagai counts and paired designs under fixed native
   and package normalization. Default native support-rank rejection withholds
   3,494 Kang and 247 Hagai eligible genes. The explicit experimental
   [active-donor NB policy](../Tools/Omics/NegativeBinomial/ActiveDonor/README.md)
   adds 3,429 Kang tests with independently checked gene-specific designs and
   uncertainty. It retains minimum replication and rank gates, leaving 65 Kang
   and 247 Hagai eligible genes untested. Correlation and numerical agreement
   do not establish selection-adjusted FDR or close broader calibration gaps.
   Full-scope original Hagai mouse count/QC, pseudobulks and a predeclared NB
   comparison now pass through the streaming route. Its three donor pairs are
   inferred explicitly from deposited sample prefixes, supported by the primary
   individual table. Cross-study FDR calibration remains open.
   A frozen [untreated-cell null benchmark](../Tools/Omics/NegativeBinomial/NullBenchmark/README.md)
   now completes ten splits each of Kang and Hagai, forty native fits/replays and
   120 direct R fits. Exact counts/designs and numerical checks pass. Native
   Hagai nevertheless reports twenty BH<0.05 sham calls across eight splits
   under each policy; Kang has none, with default support coverage still limited.
   These overlapping resamples expose a calibration concern and do not establish
   general FDR control, independent-donor replication or alternative-model power.
   A [stage audit](../Tools/Omics/NegativeBinomial/DispersionAudit/README.md)
   exactly reproduces all twenty DESeq2 references and completes eighty trend
   decompositions. Equal Hagai prior variance rules out a prior-strength
   explanation; gene-wise estimation and outlier admission remain the next
   concrete inference targets, without claiming the calibration issue fixed.
   The subsequent [profile audit](../Tools/Omics/NegativeBinomial/ProfileAudit/README.md)
   passes 2,305 valid native point checks and exposes missed better in-bounds
   objectives in the default DESeq2 references. Objective-checked Hagai reference
   sensitivity yields 25 sham calls instead of one; nine Kang parametric
   sensitivities fall back to local trends and remain failed method constraints.
   Matching the original reference does not justify a native numerical repair;
   statistical calibration and independent power evidence remain open.
3. **Negative-binomial DE:** an explicit native NB cohort model now runs through
   the count/analysis/replay/table CLI, with adjusted dispersion estimation,
   robust trend, prior shrinkage, offsets, shared paired/batch designs and
   diagnostics. The full Kang NB comparison and independent numerical checks
   pass; [evidence and limitations](../Tools/Omics/NegativeBinomial/README.md)
   remain distinct from multi-study calibration and production qualification.
   Old plans retain the log-linear baseline. Optional
   [count-likelihood contrast shrinkage](../Tools/Omics/NegativeBinomial/EffectShrinkage/README.md)
   now jointly refits nuisance coefficients under an explicit fixed normal prior,
   with conditional Laplace uncertainty and separate failure diagnostics. Original
   Wald/BH inference remains unchanged. An explicit
   [empirical weighted-quantile effect prior](../Tools/Omics/NegativeBinomial/EmpiricalPrior/README.md)
   now learns the contrast prior width from full-design available genes, records
   its contributing/excluded genes and fails independently of original inference.
   A frozen [held-out-animal count-risk study](../Tools/Omics/NegativeBinomial/HeldOutRisk/README.md)
   now evaluates all 112 overlapping Crowell training folds. Empirical shrinkage
   improves the primary depth-conditional NB log score by 0.0426 nats/gene over
   MLE, with 55/56 animal/population summaries improved; fixed shrinkage performs
   better in endothelial cells. All count isolation and conditional numerical
   checks pass. This is eight previously inspected animals, not cross-study
   effect-truth recovery or FDR/interval calibration.
   Optional [contrast likelihood-ratio inference](../Tools/Omics/NegativeBinomial/LikelihoodRatio/README.md)
   now passes 50 frozen native analyses and 443,705 independent gene/analysis
   checks, plus real H5AD publication/replay with byte-exact old Wald output.
   The null fit retains full-model dispersion and donor/support choices;
   omitted options preserve Wald behavior. Original edgeR stopping discrepancies
   remain recorded alongside a passing tighter-solver sensitivity. This does
   not close calibration: Hagai has 22 sham calls versus Wald's 20 under each
   policy, and active-donor Kang adds one across ten overlapping splits.
   The subsequent [QL stage study](../Tools/Omics/NegativeBinomial/QuasiLikelihood/README.md)
   completes 116 reference fits across 29 full-support analyses, with exact
   source checks and independently verified QL arithmetic. Native stable unit
   deviance now passes 7.88 million real fitted-count checks and 702 high-precision
   boundary cases. [Native conditional moments and adjusted residuals](../Tools/Omics/NegativeBinomial/QuasiLikelihood/Moments/README.md)
   now pass 25 focused Swift tests and independent checks of all 3,940,848
   available default-budget moments across 58 attempted real-data arms. Twenty
   gene/arm fits exhaust the declared work limit; a separately declared larger
   budget resolves nineteen, leaving one Hagai fit unavailable. The positive-mean
   scaling-underflow defect and its verified repair are retained. The subsequent
   [explicit adaptive evaluator](../Tools/Omics/NegativeBinomial/QuasiLikelihood/Moments/Adaptive/README.md)
   resolves all twenty original work-limit failures and passes all 58 arms and
   3,940,972 moments, including independent checks of the 124 recovered moments.
   Its 28 focused native tests pass; direct summation remains the default and
   the earlier failures remain recorded. The subsequent
   [native global QL scale/refit](../Tools/Omics/NegativeBinomial/QuasiLikelihood/GlobalScale/README.md)
   passes all 58 arms and 483,576 gene/arm fits, independent LOWESS/fixed-scale
   GLM comparisons and direct score checks. Its 31 focused native tests pass;
   no genes are omitted from either real-data scale update. Abundance covariates
   remain supplied reference inputs in that stage. The subsequent
   [native abundance integration](../Tools/Omics/NegativeBinomial/QuasiLikelihood/Abundance/README.md)
   derives those covariates directly, passing all 58 integrated arms,
   105 controlled cases and 35 focused native tests with original count fits
   unchanged. The final-arm SSH output stall, partial output and successful
   file-backed recovery remain retained. The subsequent
   [native robust QL moderation](../Tools/Omics/NegativeBinomial/QuasiLikelihood/Moderation/README.md)
   passes all 58 arms and 483,576 posterior rows, 265 numerical checks and
   43 focused native tests. All 111 prior profiles pass independent objective
   and trend checks; posterior error against a tightly optimized reference is
   at most 1.32e-6. Ordinary-default reference discrepancies remain recorded,
   including coarse optimization near the upper prior-DF boundary. The subsequent
   [native adjusted QL cohort test](../Tools/Omics/NegativeBinomial/QuasiLikelihood/Inference/README.md)
   completes the native chain on all 58 arms and 483,576 tests with exact upstream
   fits, residual DF, abundance and moderation. Independent constrained-fit,
   F-tail and BH checks pass, as do 50 native tests in nine suites. The existing
   pseudobulk expression owner also passes a real Kang check with 5,400 tests
   and all 15,706 original feature identities retained. The opt-in method leaves
   defaults unchanged and does not supply QL effect intervals. Modern adjusted
   QL explicitly disables the legacy-only Poisson bound. Legacy inference,
   varying-support borrowing and fresh calibration remain open. QL is not a
   demonstrated calibration fix: native-trend QL yields 39 Hagai sham calls
   versus Wald's 20, and three Kang calls versus two from the original adjusted
   reference and zero from Wald. These inspected families cannot qualify future
   method selection. Active-donor QL borrowing remains separate.
   A [fresh independent-donor experiment](../Tools/Omics/NegativeBinomial/IndependentNull/README.md)
   adds nine disjoint twelve-donor cohorts from the full public Human Immune
   Health Atlas B/Plasma release, with a protocol frozen before DE outcomes.
   Original Wald/LRT/QL and pinned edgeR/limma/DESeq2 families are retained.
   A separately declared all-cohort follow-up repairs an unrequested influence
   diagnostic dependency in singleton batches, retaining original failures and
   exact fit comparisons. One selected study and descriptive sham-call events
   do not establish universal FDR, power, interval coverage or a production default.
   Heavy-tailed priors, prior uncertainty, posterior coverage, effect-estimation
   risk and robust cross-study qualification remain.
   An explicit Gamma dispersion-trend option now passes controlled Kang/Hagai
   comparisons and independent numerical checks. The earlier robust log trend
   remains available with unchanged defaults; neither is promoted by correlation
   alone.
4. **Sparse feature selection/PCA/kNN/UMAP-compatible embeddings/clustering:**
   sparse Seurat-style HVG selection and centered matrix-free PCA now have exact
   selected-gene agreement and numerical PCA agreement with Scanpy on full Kang
   B-cell and PBMC3k matrices. Results use the existing analysis/replay artifact
   route; see [reduction evidence](../Tools/Omics/Reduction/README.md).
   Exact PCA kNN and UMAP-compatible fuzzy connectivity now have full real-data
   reference checks; see [neighbor graph](../Tools/Omics/Reduction/NEIGHBORS.md).
   Sparse multilevel Louvain now has independent objective, connectivity and
   three-seed heuristic reference comparisons on both datasets;
   [clustering limits](../Tools/Omics/Reduction/CLUSTERING.md) remain explicit.
   [File-backed clustering](../Tools/Omics/Reduction/FILE_CLUSTERING.md) shares
   that solver and stores original and aggregated edges in windowed CSR files.
   Native fixed-epoch UMAP-compatible optimization now has curve, gradient,
   schedule and real-data neighborhood-preservation checks;
   [embedding limits](../Tools/Omics/Reduction/EMBEDDING.md) remain explicit.
   Partition/embedding stability and biological validation remain. Never densify
   cells × genes; resident, pair-work and embedding-update limits still apply.
5. **Batch integration:** native single-covariate donor/batch correction now
   preserves original PCA and explicitly selects corrected downstream coordinates.
   [Kang qualification](../Tools/Omics/Reduction/INTEGRATION.md) compares donor
   mixing, cross-donor condition accuracy, and measured RNA program preservation
   with three pinned Harmony reference runs and a response-erasure control.
   [Full-cohort file integration](../Tools/Omics/Reduction/FILE_INTEGRATION.md) now
   matches the previous numerical solver exactly across three seeds on Kang and
   independent Hagai. Hagai passes its response-preservation margins; full Kang
   fails NK-cell recall preservation in both native and reference Harmony runs.
   A [native exact-tree performance trial](../Tools/Omics/Reduction/MNN_TREE.md)
   preserves all original MNN outputs on Hagai/Kang/Ding but is slower; the
   prototype remains archived and is not a production option.
   The subsequent [HNSW/local-kernel candidate](../Tools/Omics/Reduction/LocalMNN/README.md)
   completes all three original cohorts but fails type/program preservation on
   Kang and Ding despite over 99% exact-anchor recall; it is not promoted.
   The [tiled full-Gaussian option](../Tools/Omics/Reduction/GaussianKernel/README.md)
   preserves all anchors and original measured biological metrics on all three
   cohorts while accelerating correction. Its scalar underflow fallback is
   explicitly budgeted; scalar defaults and full-cohort limits remain intact.
   A [fixed matching-only intervention](../Tools/Omics/Reduction/FullGaussianMNN/README.md)
   reuses the earlier HNSW anchors with full Gaussian smoothing. Hagai and Ding
   pass the applicable measured gates, but Kang megakaryocyte recall loses 0.062218
   against the fixed 0.05 allowance. Four unavailable Kang classifier strata and
   partial Ding labels remain. Approximate matching is not promoted.
   The opt-in [native scale-aware MNN method](../Tools/Omics/Reduction/MNN_INTEGRATION.md)
   now passes the available preservation margins on full Kang, Hagai and Ding,
   with independent anchor, coordinate and neighbor agreement. It remains a
   development result on already inspected cohorts; untouched-study validation
   and Ding's partial source-label coverage remain open.
   The [full HIRISA annotation-retention diagnostic](../Tools/Omics/Benchmarks/HIRISA/ANNOTATION_RETENTION.md)
   now evaluates all 1,612,594 cells with five independent per-cell SVD checks.
   Native integration fails 37/146 supported, sensitive comparisons, including
   10/29 rare comparisons; all three Harmony candidates also fail. Another 325
   comparisons lack sufficient support or control sensitivity. These are
   annotation-recoverability measurements, not authoritative cell identities.
   Missing rare-type strata and the insensitive original Kang erasure control
   remain explicit. Full Baron is rejected for confounding. Multiple covariates,
   prospective mapping and general multi-donor competitiveness remain open.
6. **Annotation:** [standalone native binary program bundles](../Tools/Omics/Programs/BUNDLES.md)
   separate per-cell arrays from pseudobulk JSON and preserve explicit source
   IDs, exact-name matching, missing scores and native replay. [Full HIRISA
   publication, replay and independent score checks](../Tools/Omics/Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
   now pass. [Native fixed marker/program scoring](../Tools/Omics/Programs/README.md)
   now shares sparse arithmetic between resident and streamed H5AD routes, with
   exact IDs, signed weights, definition fingerprints, missing-gene coverage and
   null scores for empty libraries. Full Kang and Baron comparisons pass against
   independent sparse Scanpy-normalized products, with descriptive within-donor
   response observations. [Native frozen reference mapping](../Tools/Omics/ReferenceMapping/NATIVE.md)
   now fits training-only PCA and returns provenance-bound candidate kNN labels;
   all four Baron held-out donor folds match independent reference arithmetic.
   Optional [native balanced multinomial logistic mapping](../Tools/Omics/ReferenceMapping/Logistic/README.md)
   now fits training-only scale and class weights. All 8,569 candidate labels
   match the external classifier, and independent probability/objective checks
   pass. Macro-F1 improves over kNN in all four donors, while accuracy falls in
   two and acinar/Schwann failures persist. The default kNN bytes remain exact.
   An optional [explicit PCA feature panel and complete Kang–Ding transfer](../Tools/Omics/ReferenceMapping/CrossStudy/README.md)
   now permit different gene universes while retaining full source normalization.
   All 68,704 query cells pass independent sparse/PCA/probability checks, but both
   coarse-family targets fail (macro-F1 0.613/0.666); megakaryocyte recall is zero
   in both directions. Ding's partial label coverage and differing fine-label
   taxonomies remain explicit. The original default bytes remain exact.
   The [complete post-result diagnosis](../Tools/Omics/ReferenceMapping/CrossStudy/Diagnosis/README.md)
   verifies original Kang category labels and all query-cell transcript scores,
   then reconstructs both frozen gene-space classifiers. Same-named rare classes
   have sharply different measured RNA profiles, limiting source-label equivalence.
   Both transfer failures remain; no cells were relabelled or excluded.
   Calibrated annotation, novel-class rejection and independent multi-study
   biological qualification remain; labels are not authoritative.
7. **Perturbation prediction:** [real donor-held-out response baselines](../Tools/Omics/PerturbationPrediction/README.md)
   now evaluate all eight Kang and three Hagai donors using supplied controls and
   sealed treated outcomes. [Native donor-response prediction](../Tools/Omics/PerturbationPrediction/NATIVE.md)
   now fits and freezes all four baselines through streamed H5AD aggregation,
   with reconstruction and real-data numerical agreement on all eleven folds.
   The [complete Norman filtered release](../Tools/Omics/PerturbationPrediction/Norman/README.md)
   supplies 105 single-target and 131 paired-target conditions with explicit
   unseen-combination and unseen-target splits.
   [External composition baselines](../Tools/Omics/PerturbationPrediction/Norman/COMBINATIONS.md)
   now score every held-out pair with sealed outcomes, exact frozen replay and
   a target-shuffled control. [Native composition prediction](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_COMPOSITION.md)
   now selects verified condition counts, fits a frozen model and reproduces all
   786 query/method vectors exactly. Neither these baselines nor the donor models
   qualify genetic interactions.
   An [external unseen-target experiment](../Tools/Omics/PerturbationPrediction/Norman/UNSEEN_TARGETS.md)
   now scores all 105 target folds, with descriptor coverage for 102. Its
   co-response ridge model underperforms the mean-single baseline; this does not
   qualify a native unseen-target predictor. A [control-only descriptor follow-up](../Tools/Omics/PerturbationPrediction/Norman/CONTROL_DESCRIPTORS.md)
   covers 97 targets but also fails to beat the mean-response baseline.
   A [GO-informed target-kernel experiment](../Tools/Omics/PerturbationPrediction/Norman/GO_TRANSFER.md)
   covers 101 targets: fixed regularization modestly improves matched mean and
   shuffled baselines, while nested selection fails the primary all-gene comparison.
   The [native target-kernel owner](../Tools/Omics/PerturbationPrediction/Norman/NATIVE_TARGET_KERNEL.md)
   now fits reusable annotation models, rejects seen target IDs and known aliases,
   retains unsupported descriptors and reproduces the fixed GO result across all
   105 held-target folds. Its 517 available prediction vectors agree with frozen
   references within 2.45e-15; the 0.58% mean-RMSE gain remains a development result
   with substantial per-target failures.
   The subsequent [complete Replogle validation](../Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
   retains 30 targets across five technical gemgroups and completes all 150
   native folds. It meets the fixed all-gene criterion in all five groups with
   1.41–2.12% lower mean RMSE than the training mean, while 34/150 folds remain
   worse than no change. All 750 vectors pass independent checks and repeated
   scoring is identical. This separately collected K562/UPR experiment shares
   investigators and earlier target selection; it does not supply independent
   laboratory, temporal-knowledge or new-tissue validation.
   The [Kang–HIRISA cross-study comparison](../Tools/Omics/PerturbationPrediction/CrossStudyIFNB/README.md)
   now completes 26 cross/within donor folds over 11,884 exact shared symbols.
   An explicit native panel preserves full source-library denominators. Numerical
   checks pass for all 104 vectors, but ridge fails both directional primary
   comparisons and is worse than cross mean for all thirteen donors. These reused
   studies expose a combined context-transfer failure, not prospective validation.
   Optional [normal-model donor-response intervals](../Tools/Omics/PerturbationPrediction/Intervals/README.md)
   now pass all 26 native folds and 104 independent bound-array checks, retaining
   every prior point vector. Nominal 95% treated-expression coverage ranges from
   42.39% to 99.21% between the two transfer directions, with different widths and
   unavailable genes. This is implemented uncertainty output, not general
   calibration or Bayesian/mechanistic coupling.
   The [complete external GSE226572 duration-transfer test](../Tools/Omics/PerturbationPrediction/GSE226572/README.md)
   checks all 510.7 million raw stored values and 126,633 cells admitted by frozen
   initial QC. Native aggregation, fitting, all three donor predictions and
   independent checks pass. Across all 18 released treatment-time profiles,
   mean response improves RMSE by 2.00%, below its 5% target, with 8/18 cases worse
   than no change; ridge also fails. Nominal 95% coverage spans 85.37–95.45%.
   All five repeated scoring outputs are byte-exact. This whole-population,
   fixed-duration transfer is distinct from B-cell and learned temporal prediction.
   The subsequent [native duration owner](../Tools/Omics/PerturbationPrediction/Duration/README.md)
   adds `singlecell-duration-fit`, `singlecell-duration-predict` and reconstructing
   model/prediction verifiers. It uses explicit exposure hours, within-donor
   log-time interpolation, control-only query donors and full RNA denominators;
   extrapolation and identity mismatches are rejected. All three GSE226572 donor
   holdouts and 18 outcomes pass native/independent checks. Duration mean improves
   RMSE by 55.71% against no change and 24.83% against matched time-invariant mean;
   context ridge fails its secondary gate. This is development on inspected data,
   with only two training donors per fold, broad nominal intervals and two losses
   to the matched mean. Independent temporal/context validation remains open.
   The [complete joint-count donor study](../Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/README.md)
   now fits, predicts and scores all thirteen Kang/HIRISA donor omissions. Kang
   improves over training mean; HIRISA loses in all five donors. Later response
   shrinkage and transfer experiments retain those original failures.
   The [three-study held-out native experiment](../Tools/Omics/CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/ResponseShrinkage/StudyHeldOut/README.md)
   selects penalties using only training studies and balances study weights.
   All 75 donors, 11,800 genes and 885,000 predictions pass independent numerical
   checks. It improves slightly over balanced training mean in Kang/GSE181897,
   but fails both baselines in every HIRISA donor; no study reaches a 5% gain over
   both baselines. GSE181897 is explicitly reused for development here; its
   earlier external-test result is unchanged. These completed experiments do
   not establish calibrated uncertainty or reliable context transfer.
   Bayesian/mechanistic integration, reliable unseen-target gene prediction, unseen
   cell/tissue contexts and single-cell response distributions remain open.
8. **Multimodal:** the [native multi-assay core and 10x CITE-seq path](../Tools/Omics/Multimodal/README.md)
   now preserve independent feature spaces, assay row maps, exact counts, genomic
   intervals and spatial frames. All 5,247 public PBMC5k cells and both RNA/protein
   assays match h5py/SciPy, Scanpy and MuData; native H5MU export preserves missing
   assay rows. Native H5MU import now preserves explicit modality maps and selected CSR/CSC/dense
   count layers, with the full real CITE-seq dataset matching across native and
   MuData-written sources. The [complete real paired RNA/ATAC benchmark](../Tools/Omics/Multimodal/MULTIOME.md)
   now preserves all 2,711 nuclei, 36,601 genes and 98,319 peaks, with explicit
   cut-site units and exact agreement across native 10x and independent MuData
   imports. The [complete Visium spatial interchange benchmark](../Tools/Omics/Multimodal/SPATIAL.md)
   now checks all 4,039 spots, counts and pixel coordinates through native H5MU
   reimport, including standard spatial-array export. The [native Visium directory
   importer](../Tools/Omics/Multimodal/Visium/README.md) now reads original filtered
   RNA HDF5 and explicitly selected legacy/header position CSVs. Both formats
   preserve every count and coordinate in the full release and reproduce the
   earlier qualified dataset/H5MU bytes. Native replay, 13 Swift tests and 16
   lifecycle/regression commands pass. This is a bounded resident RNA/spot path;
   HD Parquet, Matrix Market directories and biological ATAC/spatial analysis,
   continuous measurements and joint multimodal analysis remain open.
9. **Out-of-core:** streamed H5AD normalization/HVG and explicitly memory-mapped
   selected-entry PCA now pass full Baron and Hagai comparisons against Scanpy;
   see [storage and qualification](../Tools/Omics/Reduction/STREAMING.md).
   Selected-entry PCA now uses 16 MiB mapping windows with a shared fixed-buffer
   snapshot reader, replacing the whole-cache map. Full Baron/Hagai native PCA
   outputs remain exact; independent Scanpy comparisons pass. The complete
   111,445-cell Norman matrix now also passes 2,000-feature/20-component PCA
   against Scanpy and native reconstruction, with 5.8 billion windowed entry
   visits. Its publication/verification peaks of 1.38/1.49 GB expose the remaining
   resident report and aggregate costs. A subsequent [standalone PCA bundle](../Tools/Omics/Reduction/PCA_BUNDLE.md)
   removes unused condition aggregation and JSON score/loading encoding. Complete
   Norman publication/reconstruction now peak at 443.2/443.4 MB; the same product's
   legacy route peaks at 1.386 GB. All scores, loadings, identities, QC and feature
   statistics remain exact, and independent Scanpy checks include stored centers.
   Full Norman count aggregation now exercises 361.6 million source entries
   through bounded slices, with 3.58 million aggregate nonzeros.
   A [persistent count store with windowed normalization](../Tools/Omics/CountStore/README.md)
   now retains all 361.6 million Norman entries as exact binary records. Native
   import, normalization and reconstruction verification each stay below 286 MB
   maximum resident memory on 111,445 cells. Every raw record matches SciPy and
   every normalized record matches Scanpy within 8.9e-16 absolute error.
   [Frozen query projection](../Tools/Omics/Reduction/PCA_QUERY.md) now applies
   training-only centers/loadings with query scores in 16 MiB file windows. All
   8,569 Baron cells pass donor-held-out projection checks against Scanpy/SciPy;
   existing label-reference reports remain exact. Training metadata, moments,
   basis and fitted scores still remain resident. [File-backed exact neighbors](../Tools/Omics/Reduction/WINDOWED_NEIGHBORS.md)
   now use independent Dispatch workers and bounded score tiles. Complete Baron,
   Hagai and held-out human3 graphs match SciPy membership/distances and umap-learn
   topology, with identical serial/parallel graph bytes. Final graph arrays remain
   resident and exact search remains quadratic. Optional [native HNSW](../Tools/Omics/Reduction/HNSW_NEIGHBORS.md)
   now builds the complete 111,445-cell Norman graph with a bounded score cache;
   strict mean recall is 99.9407% on 2,048 preselected exact-query checks. Complete
   Baron/Hagai recall exceeds 99.996%. A [binary graph store](../Tools/Omics/Reduction/GRAPH_STORE.md)
   now streams exact/HNSW neighbor rows and constructs CSR connectivity through a
   disk transpose and row merge. Graph arrays and their large JSON encoding are
   avoided. Complete Norman construction peaks at 477 MB versus 1.421 GB for
   same-executable JSON output, with every binary FP64 value exact to the qualified
   graph. Input PCA state, HNSW index and cell-scale bookkeeping remain resident.
   File-backed Louvain now consumes this graph and writes each aggregated level
   through stable disk scatter, preserving the previous solver's summation order.
   [File-backed embedding](../Tools/Omics/Reduction/FILE_EMBEDDING.md) now stores
   mutable edge schedules in a 16 MiB mapping window and reads only the requested
   PCA initialization columns into resident arrays. File-backed integration now
   keeps latent matrices in bounded mappings and is measured on full Kang/Hagai;
   [Complete HIRISA graph qualification](../Tools/Omics/Benchmarks/HIRISA/GRAPH_RESULTS.md)
   and [native seed-7 integration](../Tools/Omics/Benchmarks/HIRISA/FULL_INTEGRATION_RESULTS.md)
   now pass operational/numerical gates on all 1,612,594 cells. [Full-cohort
   clustering](../Tools/Omics/Benchmarks/HIRISA/FULL_CLUSTERING_RESULTS.md) and
   [native program scoring](../Tools/Omics/Benchmarks/HIRISA/NATIVE_PROGRAM_RESULTS.md)
   also pass publication, replay and independent checks. Broader biological
   preservation and native multi-seed integration remain open; resident
   metadata/bookkeeping limits remain explicit.
   The [complete Replogle UPR partition](../Tools/Omics/PerturbationPrediction/Replogle2020/README.md)
   also passes on 129,839,577 mixed-source entries, with RNA and guide counts
   separated. Projection uses charged work and an explicit output-byte allowance
   while retaining bounded transfers. Full and confident-cohort native RNA sums,
   QC and membership match independent AnnData/SciPy. The subsequent
   [fixed target-prediction experiment](../Tools/Omics/PerturbationPrediction/Replogle2020/RESULTS.md)
   completes all 150 folds and meets the primary comparison in all five technical
   groups, with modest gains and target-level failures. This adds expression
   evidence; it does not qualify full H5MU processing or new biological contexts.
10. **Metal:** the first explicit [sparse count-normalization backend](../Tools/Omics/CountStore/Metal/README.md)
    now passes all 14,184,532 original Kang records, every native replay and the
    actual product single-cell CLI. FP32 output has maximum absolute error
    1.24e-6 against FP64 Scanpy; raw counts/identities and default CPU bundles
    remain exact. M4 Pro end-to-end medians are 1.043 s CPU and 1.050 s Metal,
    establishing no speedup. CPU remains default. This is a bounded GPU
    arithmetic result, not million-cell GPU or biological qualification.
    The [profile-guided shared writer](../Tools/Omics/CountStore/Metal/Profile/README.md)
    preserves every bundle byte, with final medians 0.918 s CPU / 0.912 s Metal.
    CPU's slower first call leaves its three-run mean slightly worse; this does
    not establish general acceleration or a GPU advantage. Final release/ASAN
    record checks and thirteen actual CLI checks pass.
    PCA/kNN/model-fitting acceleration and matched end-to-end scverse speed and
    memory comparisons remain open.
11. **Other omics:** genomics/variants, bulk RNA, proteomics, metabolomics, spatial
    imaging, metabolic modeling and regulatory prediction on shared provenance
    remain outside this initial interoperability block.
12. **Cross-scale biology:** variant → regulation → RNA/cell state → protein and
    molecular mechanism → reaction/kinetics → cellular phenotype → tissue
    prediction requires executable, independently qualified links at each boundary.
    [AlphaGenome Atlas integration assessment](AlphaGenomeAtlas.md) identifies
    an external prediction source for the variant/regulation/RNA link. A separate
    [native genomic evidence layer and retrieval adapter](Design/ALPHAGENOME_ATLAS.md)
    already exists for bounded public-reference research. RNA/ATAC integration,
    live service qualification and downstream biological coupling remain open.
## Source-bound aggregation cohorts

[Cell selection](../Tools/Omics/H5AD/CELL_SELECTION.md) now applies explicit,
source-fingerprinted observation indices during native streamed pseudobulk
aggregation. It retains the full input and maps every selected report row back
to the original axis. The [Adamson author cohort](../Tools/Omics/PerturbationPrediction/Adamson/COHORT.md)
verifies 50,440 cells and 781,977,660 UMIs against independent sparse aggregation;
control definitions and predictive scoring remain pending.

### GSE181897 count handoff and external prediction (2026-09-11)

The [complete source audit and native handoff](../Tools/Omics/PerturbationPrediction/GSE181897/README.md)
validates 292,741,570 original RNA/antibody records and matches each modality's
library total for all 136,142 cells. External indexed-gzip extraction preserves
all 15,272 author B-lineage cells across 64 donors and all original condition
codes, retaining 20,303 RNA features. Native streamed aggregation processes
34,287,682 selected records and agrees exactly on all 3,818,645 nonzero values
in 379 observed groups. Reconstruction and repeated report bytes agree.

This count handoff does not qualify native direct gzip-H5AD import, full AnnData
slot preservation or multimodal outcome prediction. Historical primary author
code subsequently resolved the treatment codes. The [complete external prediction
test](../Tools/Omics/PerturbationPrediction/GSE181897/RESULTS.md) now checks all
124 native predictions across 62 paired donors and 11,800 shared genes.
HIRISA-trained mean passes the frozen 5% improvement target with a 5.62% gain;
Kang-trained mean fails. Independent checks and repeated scoring pass, but
HIRISA-trained nominal 95% treated-expression interval coverage averages only
35.54%. All incomplete pairs, absent genes, individual failures and insufficient
uncertainty remain reported.
