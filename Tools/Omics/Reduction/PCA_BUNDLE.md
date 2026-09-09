# Standalone H5AD PCA bundles

`singlecell-h5ad-pca` performs cell QC, all-gene normalization/HVG selection and
centered sparse PCA without computing condition-level pseudobulk aggregates.
It reuses `VivoH5ADReduction` and `VivoSingleCellReduction`; it introduces no
second PCA engine. Streamed PCA and pseudobulk now share their cell-quality
accumulator, including exact integer totals and missing mitochondrial fractions.

```sh
numivivo singlecell-h5ad-pca original.h5ad --plan pca-plan.json --output new-pca
numivivo singlecell-h5ad-pca-verify new-pca
```

The plan contains `schemaVersion: 1`, the existing explicit H5AD `mapping`, and
optional `reduction` settings from [streamed PCA](STREAMING.md). Projection centers
are retained by default and required; explicitly disabling them is rejected.
The canonical plan records this choice. There are no contrast or pseudobulk
settings in this standalone plan. Normalization still includes every source gene.

## Persistent representation

| Artifact | Meaning |
| --- | --- |
| `original.h5ad` | Exact immutable source snapshot |
| `plan.json` | Canonical mapping, normalization and PCA parameters |
| `metadata.json` | Complete source cell/feature identities and count provenance |
| `quality.json` | Per-cell exact counts, detected features and mitochondrial QC |
| `model.json` | All feature statistics, selection, projection centers, eigenvalue/residual diagnostics and storage/work accounting |
| `scores.bin` | Complete cell-major FP64 score matrix |
| `loadings.bin` | Complete selected-feature-major FP64 loading matrix |
| `receipt.json` | Separate fingerprints for source, plan, metadata, QC, model, scores, loadings and implementation |

Each matrix record is exactly 16 bytes: little-endian UInt32 row, UInt32 component,
then the Float64 value bits. Every row/component coordinate is emitted, including
zero and negative values. Rows follow metadata cell order for scores, and
`model.selectedFeatureIndices` for loadings; component indices start at zero.
The format identifier is `complete-row-major-u32-row-u32-component-f64-le/v1`.
These are transformed coordinates, never raw count matrices. The shared record
writer buffers 1 MiB; scores/loadings are not serialized as JSON matrices.

Verification reconstructs all artifacts from the retained H5AD and canonical
plan and checks every receipt-bound file. Rehashing altered score coordinates
or projection centers does not bypass reconstruction. Publication requires a
new destination and removes private staging on thrown failures. Private scratch
uses the shared 16 MiB record windows; the source snapshot uses fixed-buffer I/O.

The H5AD route retains its existing 1 GiB source, one-million-observation,
100,000-feature and one-billion-source-entry admission bounds. Selected scratch
is capped at two billion bytes; sparse entry visits default to two billion and
can be explicitly raised to twenty billion. These bounds are not qualifications
at all allowed sizes. Metadata, quality, feature statistics, fitted score/loading
arrays and PCA basis remain resident. This removes unused aggregation and large
JSON matrix encoding, not every resident structure.

## Independent checks

`check_pca_bundle.py` exercises native publication, reconstruction, exact repeats,
CSR/CSC equality, a zero-count cell, budget/center rejection, cleanup and rehashed
score/model tampering. `pca_bundle_reference.py` validates all matrix coordinates
and values and adapts a verified bundle to the existing independent Scanpy checker:

```sh
python3 Tools/Omics/Reduction/pca_bundle_reference.py \
  --bundle new-pca --out reference-input
python3 Tools/Omics/Reduction/check_reference.py \
  --h5ad original.h5ad --report reference-input/report.json --out reference.json
```

The reference adapter is external evidence tooling. It materializes arrays and
JSON for comparison; it is not the native execution path. The Scanpy checker also
validates stored projection centers against the selected normalized matrix means.

A fitted bundle can preserve preprocessing for later consumers, but fitting it
on every condition is not train-only preprocessing for a held-out perturbation
experiment. No learned annotation, integration, million-cell execution, parallel
kernel, GPU or biological-prediction qualification follows from this artifact alone.

## Complete Norman qualification

The [archived benchmark](evidence/2026-09-09-pca-bundle/README.md) covers all
111,445 cells, 33,694 genes and 361,582,621 source entries with 2,000 HVGs and
20 PCs. Standalone publication took 136.11 s at 443,203,584 maximum resident
bytes; reconstruction took 137.46 s at 443,432,960 bytes. The same product's
legacy pseudobulk/PCA route took 158.85 s at 1,386,496,000 bytes.

All scores, loadings, feature statistics, cell identities and QC exactly match
the previous full Norman result; the legacy report remains byte-for-byte
unchanged. Independent Scanpy checks pass, including retained centers, with
aligned relative score error 1.0971e-12. Scores occupy 35,662,400 binary bytes
and loadings 640,000 bytes. The gate also passed 69 Swift tests in 19 suites,
12 standalone CLI commands with 7 expected rejections, and existing streamed
H5AD/legacy CSR/CSC regression checks.

## Applying a training fit

[Windowed frozen query projection](PCA_QUERY.md) uses the retained centers and
binary loadings on disjoint new cells. Fits intended for this route must declare
`featureNamespace`; complete feature IDs, organism and count units must match.
