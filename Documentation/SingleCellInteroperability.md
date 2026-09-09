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
`groupColumn` and `featureNameColumn` are optional. Missing group annotations are
preserved; missing selected design identities are rejected. Mitochondrial feature
identities are explicit, never inferred from names.

## Preservation contract

- Import publishes the **unchanged original.h5ad**, a native `dataset.json`, the
  mapping, a source/dataset/mapping fingerprint receipt with executable identity
  and HDF5 version, and a MEX manifest usable by the existing replayable workflow.
- The count projection reads X or a named layer. It supports CSR, CSC, and dense
  arrays read through one-row HDF5 hyperslabs. It never allocates cells × genes.
  Sparse duplicates are summed with checked UInt64 arithmetic, indices sorted,
  and explicit zero entries removed. The original representation stays in the
  source file.
- Integer counts retain all UInt64 bits. Floating counts must be finite,
  nonnegative, integral, and no larger than 2^53. Selecting normalized fractional
  data fails instead of converting it to counts.
- Selected metadata columns support strings, categorical strings, and nullable
  strings. Root, dataframe, and array encoding versions are checked against the
  [AnnData on-disk specification](https://anndata.readthedocs.io/en/stable/fileformat-prose.html).
- Source raw, additional layers, nullable/numeric metadata, category order,
  embeddings, graphs and uns are retained in the original file. These are **not
  yet editable native AnnData fields**. Exporting that source returns the original
  object; it does not claim to include later filtering or analysis.
- Writing `dataset.json` creates a **new count AnnData object**, with CSR X,
  observation design columns, nullable donor/group columns, feature names,
  mitochondrial annotations and native metadata in uns. Unique generated obs
  indices coexist with original barcode/sample columns. Use `barcodeColumn:
  "barcode"`, `sampleColumn: "sample"`, `groupColumn: "group"`, and
  `featureNameColumn: "name"` when importing native output again.
- Publication refuses existing destinations. The current source limit is 64 MiB;
  count processing remains bounded at 100,000 cells/features and 2 million
  nonzeros by default. This is not out-of-core execution. HDF5 itself may retain
  decompression and variable-string buffers beyond the Swift sparse arrays.

## Validation

```
bash Tools/Omics/H5AD/build.sh /tmp/numivivo-h5ad-build
python Tools/Omics/H5AD/check_interop.py --binary /tmp/numivivo-h5ad-build/h5ad-check --out /tmp/h5ad-checks
```

Use the versions pinned in `Tools/Omics/H5AD/requirements.txt`. The test creates
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

## Required development order and remaining evidence

The complete development objective remains open:

1. **AnnData/H5AD:** native count exchange implemented; native editing/round-trip
   of general annotated objects, raw-axis selection, broader encoding coverage,
   and large real-file qualification remain.
2. **Experimental benchmarks:** several public datasets with known donors,
   perturbations and expected biology; same-data Scanpy and edgeR/limma/DESeq2
   comparisons remain to be implemented and run.
3. **Negative-binomial DE:** likelihood, dispersion estimation and shrinkage,
   offsets, paired/batch designs and diagnostics remain. Existing moderated
   log-linear DE is a baseline, not a production NB method.
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
