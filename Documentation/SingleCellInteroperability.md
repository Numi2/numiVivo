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
```

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
- Publication refuses existing destinations. The current source limit is 64 MiB;
  count processing remains bounded at 100,000 cells/features and 5 million
  nonzeros by default. This is not out-of-core execution. HDF5 itself may retain
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
one paired-donor Kang comparison, Haber count/QC checks and a rejected
non-count Hagai input. The multi-dataset donor-DE requirement remains open.

```
python Tools/Omics/H5AD/check_annotations.py --binary /path/to/numivivo --full-product --out /tmp/annotation-checks
python Tools/Omics/H5AD/check_public_data.py --binary /path/to/numivivo --full-product --out /tmp/public-file-check
```

## Validation

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
   and source-preserving annotation edits implemented. Axis-changing operations
   (cell/feature filtering/reordering with every aligned slot), direct legacy
   encoding support, and out-of-core projections remain. Full PBMC3K count/QC
   and annotation preservation now pass; this single library is not a donor-DE benchmark.
2. **Experimental benchmarks:** one eight-donor Kang B-cell contrast now passes
   exact Scanpy QC/pseudobulk checks and a descriptive PyDESeq2 comparison; Haber
   tuft-cell count/QC passes. Several independent donor-resolved studies, robust
   reference sensitivity and R edgeR/limma/DESeq2 comparisons remain.
3. **Negative-binomial DE:** an explicit native NB cohort model now runs through
   the count/analysis/replay/table CLI, with adjusted dispersion estimation,
   robust trend, prior shrinkage, offsets, shared paired/batch designs and
   diagnostics. The full Kang NB comparison and independent numerical checks
   pass; [evidence and limitations](../Tools/Omics/NegativeBinomial/README.md)
   remain distinct from multi-study calibration and production qualification.
   Old plans retain the log-linear baseline. Effect shrinkage, fuller nuisance
   handling and robust cross-study qualification remain.
4. **Sparse feature selection/PCA/kNN/UMAP-compatible embeddings/clustering:**
   implementation and reference comparisons remain; never densify cells × genes.
5. **Batch integration:** donor/batch-aware methods and biological-signal
   preservation evaluation remain.
6. **Annotation:** deterministic marker/program scoring with provenance, then
   learned reference mapping; automatic labels must not become authoritative.
7. **Perturbation prediction:** held-out perturbation/donor/context experiments,
   connection to Bayesian/mechanistic owners and real-data evaluation remain.
8. **Multimodal:** multiple feature spaces per cell for ATAC, CITE-seq, paired
   RNA/ATAC and spatial assays remain.
9. **Out-of-core:** chunked sparse storage, memory mapping, streaming transforms
   and native parallel kernels remain; current bounds do not qualify million cells.
10. **Metal:** only after stable algorithms; end-to-end CPU/scverse speed and
    memory comparisons remain for sparse transforms, PCA/kNN and model fitting.
11. **Other omics:** genomics/variants, bulk RNA, proteomics, metabolomics, spatial
    imaging, metabolic modeling and regulatory prediction on shared provenance
    remain outside this initial interoperability block.
12. **Cross-scale biology:** variant → regulation → RNA/cell state → protein and
    molecular mechanism → reaction/kinetics → cellular phenotype → tissue
    prediction requires executable, independently qualified links at each boundary.
