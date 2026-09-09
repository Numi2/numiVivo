# AnnData interoperability and the single-cell development sequence

NumiVivo now reads H5AD through the native HDF5 library. Python is only used by
an independent interoperability test, not by import, export, or analysis.
Install HDF5 (`brew install hdf5` on macOS), or set `NUMIVIVO_HDF5_LIBRARY` to
its shared library. Other workflows do not require HDF5. HDF5 calls are serialized
because installed libraries may not be thread safe.

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
- Publication refuses existing destinations. Resident count import and annotation
  editing use a 64 MiB source limit; resident count processing is bounded at
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

The output is an exchange bundle containing `original.h5ad`, `plan.json`,
`report.json` and `receipt.json`. Verification copies the source to a private
snapshot, checks source/plan/report/executable fingerprints and reconstructs
the report. Existing destinations are refused. This route uses the common
fingerprint types but does not yet publish its report into the artifact-store
DAG. It is not a substitute for general chunked transforms or sparse PCA.

Explicit limits are 1 GiB source bytes, 1 million cell identities, 100,000
features, 1 billion source sparse entries, 5 million aggregate nonzeros,
2 MiB encoded plan and 512 MiB encoded report. The CLI's common plan reader
has a tighter 128 KiB input limit. These are bounds, not demonstrated scale:
metadata, group membership and JSON reports remain resident. Compression
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
   General legacy encodings, ragged/structured aligned formats, duplicated axes
   and files beyond the explicit storage/work bounds remain. Full PBMC3K count/QC
   and annotation preservation now pass; this single library is not a donor-DE benchmark.
2. **Experimental benchmarks:** one eight-donor Kang B-cell contrast now passes
   exact Scanpy QC/pseudobulk checks and a descriptive PyDESeq2 comparison; Haber
   tuft-cell count/QC passes. Several independent donor-resolved studies, robust
   reference sensitivity and calibration remain. Direct
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
3. **Negative-binomial DE:** an explicit native NB cohort model now runs through
   the count/analysis/replay/table CLI, with adjusted dispersion estimation,
   robust trend, prior shrinkage, offsets, shared paired/batch designs and
   diagnostics. The full Kang NB comparison and independent numerical checks
   pass; [evidence and limitations](../Tools/Omics/NegativeBinomial/README.md)
   remain distinct from multi-study calibration and production qualification.
   Old plans retain the log-linear baseline. Effect shrinkage, fuller nuisance
   handling and robust cross-study qualification remain.
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
   Multiple independent donor studies, diverse cell types, rare states, multiple
   simultaneous covariates and prospective mapping remain unqualified. This does
   not yet establish competitive general multi-donor integration.
6. **Annotation:** [native fixed marker/program scoring](../Tools/Omics/Programs/README.md)
   now shares sparse arithmetic between resident and streamed H5AD routes, with
   exact IDs, signed weights, definition fingerprints, missing-gene coverage and
   null scores for empty libraries. Full Kang and Baron comparisons pass against
   independent sparse Scanpy-normalized products, with descriptive within-donor
   response observations. [Native frozen reference mapping](../Tools/Omics/ReferenceMapping/NATIVE.md)
   now fits training-only PCA and returns provenance-bound candidate kNN labels;
   all four Baron held-out donor folds match independent reference arithmetic.
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
   imports. Biological ATAC/joint-model and real spatial qualification, continuous
   measurements and joint multimodal analysis remain open.
9. **Out-of-core:** streamed H5AD normalization/HVG and explicitly memory-mapped
   selected-entry PCA now pass full Baron and Hagai comparisons against Scanpy;
   see [storage and qualification](../Tools/Omics/Reduction/STREAMING.md).
   Full Norman count aggregation now exercises 361.6 million source entries
   through bounded slices, with 3.58 million aggregate nonzeros.
   Metadata, moments, basis, scores and report remain resident. Streaming
   downstream graphs/integration, parallel kernels and million-cell qualification
   remain open.
10. **Metal:** only after stable algorithms; end-to-end CPU/scverse speed and
    memory comparisons remain for sparse transforms, PCA/kNN and model fitting.
11. **Other omics:** genomics/variants, bulk RNA, proteomics, metabolomics, spatial
    imaging, metabolic modeling and regulatory prediction on shared provenance
    remain outside this initial interoperability block.
12. **Cross-scale biology:** variant → regulation → RNA/cell state → protein and
    molecular mechanism → reaction/kinetics → cellular phenotype → tissue
    prediction requires executable, independently qualified links at each boundary.
    [AlphaGenome Atlas integration assessment](AlphaGenomeAtlas.md) identifies
    a candidate external prediction source for the variant/regulation/RNA link;
    no adapter or biological qualification is claimed.
