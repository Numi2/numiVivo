# Sparse feature selection and PCA

Enable `reduction` in the existing `singlecell-analyze` plan:

```json
{
  "schemaVersion": 1,
  "id": "sparse-pca",
  "normalizationTarget": 10000,
  "contrasts": [],
  "reduction": {
    "highlyVariableFeatures": 2000,
    "meanBins": 20,
    "components": 20,
    "maximumBasis": 128,
    "relativeResidualTolerance": 0.000001,
    "seed": 7
  }
}
```

The optional result is included in `singlecell-analysis-export` and reconstructed
by `singlecell-analysis-verify`. Plans without reduction retain their earlier
encoding and behavior. Scores carry `(sampleID, barcode)` identities in retained
cell order; loadings follow `selectedFeatureIndices` in source feature order.
Never attach scores to an H5AD object by row position without mapping these axes.

HVG selection follows Scanpy's Seurat dispersion method on library-normalized
counts recovered with `expm1` from the common log-normalized view. Moments include
implicit zeros and sample-variance correction; equal-width mean bins are right
closed. Undefined normalized dispersions are missing, not zero. All ties at the
nth finite score are retained. The common library-size denominator includes all
source genes, including genes subsequently excluded from PCA.

PCA centers the selected sparse log-normalized matrix without gene scaling.
Covariance products use sparse row traversal. The only dense data are the bounded
Krylov basis, its small projected eigensystem, scores and loadings; there is no
dense cells-by-genes matrix or full feature covariance. Twice reorthogonalized
Krylov iteration uses a deterministic seed. Component residuals and loading
orthogonality must pass before a result can be published. Insufficient basis or
numerical rank is an error; there is no silent fallback or tolerance relaxation.
The bounds are currently 10,000 selected features (including ties), 64 components
and 256 basis vectors, within the existing resident count limits. This is not yet
an out-of-core PCA implementation or a million-cell performance qualification.

## Verification

`check_cli.py` exercises the production count/artifact path, deterministic replay,
and controlled rejection of an insufficient basis. `check_reference.py` compares
all feature statistics and exact selected membership against Scanpy, then sparse
ARPACK PCA eigenvalues, variance ratios, subspace and aligned scores. It uses the
same retained cell identities and original raw imported counts. Reference code
never densifies the expression matrix.

```sh
python check_cli.py --binary /path/to/numivivo --imported /path/to/imported --out /new/native
python check_reference.py --dataset /path/to/imported/dataset.json --report /new/native/report.json --out /new/reference.json
```

The 2026-09-09 comparisons use all 2,651 Kang B cells × 15,706 genes and all 2,700
PBMC3k cells × 32,738 genes, each with 2,000 selected genes and 20 components.
Source provenance and preparation are retained in the existing experimental
benchmark suite. These establish numerical agreement on two real datasets;
they do not establish biological validity, batch integration, neighborhood graph,
clustering, or perturbation prediction. Those requirements remain open.

The preserved `pre-integration` evidence records the successful checks before
the independent Atlas contract commit reached main. Final publication evidence
is recorded separately after rebuilding the combined tree. Build warnings from
existing unrelated owners are retained in the build log.

Final combined-tree evidence: 23 Swift tests in five suites; six production
workflow assertions per dataset; 18 existing CLI assertions across 24 commands.
Both datasets retain exactly the same selected genes as Scanpy. Aligned relative
PCA score errors are 2.924e-12 (Kang) and 6.906e-12 (PBMC3k); maximum native
component residuals are 2.467e-12 and 7.573e-12 respectively. All means, sample
variances, bin assignments, dispersions and normalized dispersions are checked.
Scanpy stores normalized dispersions in float32, which determines their explicit
comparison tolerance. See `evidence/2026-09-09/source-state.json` for source and
binary identities, dataset reports for reference versions, and preserved stderr
for the intentional insufficient-basis rejections.

For the next stage, see [exact neighbors and fuzzy connectivity](NEIGHBORS.md).

[Multilevel Louvain clustering](CLUSTERING.md) consumes the verified sparse graph.

[Native UMAP-compatible coordinates](EMBEDDING.md) optimize the same sparse graph.

[Donor/batch integration](INTEGRATION.md) retains original PCA and requires an
explicit corrected representation for downstream graphs. Its real-data checks
measure donor mixing and response preservation separately.
