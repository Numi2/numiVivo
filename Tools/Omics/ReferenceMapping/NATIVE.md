# Native frozen reference mapping

`singlecell-reference-fit` trains sparse HVG/PCA from a labelled training H5AD and
archives the exact source, mapping, label provenance, feature namespace, complete
gene universe, normalization target, selected features, training means, loadings,
scores and labels. `singlecell-reference-map` projects an independently supplied
query and applies exact uniform kNN voting using the existing bounded neighbor
heap. No dense cells-by-genes or query-by-reference distance matrix is built.

```sh
numivivo singlecell-reference-fit training.h5ad --plan fit.json --output reference
numivivo singlecell-reference-verify reference
numivivo singlecell-reference-map query.h5ad --plan query.json --reference reference --output mapped
numivivo singlecell-reference-map-verify mapped
```

Fit plan:

```json
{
  "schemaVersion": 1,
  "mapping": { "...": "the complete existing H5AD mapping, with groupColumn naming source labels" },
  "reduction": {
    "normalizationTarget": 10000,
    "pca": { "highlyVariableFeatures": 2000, "meanBins": 20, "components": 20, "maximumBasis": 128, "seed": 7 }
  },
  "featureNamespace": "HUMAN_GENE_SYMBOL",
  "labelProvenance": "Exact source and version of the supplied reference labels",
  "neighbors": 15
}
```

Query plan uses `schemaVersion`, its own complete `mapping`, and the same
`featureNamespace`. It must omit `mapping.groupColumn`: the inference code does
not read query labels. Count unit and organism must match the reference. All
source genes must be present with exact IDs, although their order may differ.
Missing or additional genes are rejected because they would change the library
normalization denominator; no alias guessing or zero imputation occurs.

The training fit enables optional `retainProjectionCenters` on the existing PCA
owner. Ordinary reduction plans omit the new option and keep their old output
encoding. Query normalization uses its own full library sum, but feature
selection, centering and loadings remain frozen from training. Empty query
libraries have null scores and no candidate label; empty reference libraries are
rejected. Reference/query `(sampleID, barcode)` overlaps are rejected. Shared
donors are reported explicitly; absence of overlap alone does not prove an
independent biological evaluation.

Neighbor ordering is squared Euclidean distance followed by reference row index.
Class order is sorted source labels; tied votes select the first class. Each
mapped row retains scores, reference neighbor indices/distances, integer class
votes and a `candidateLabel`. Dividing votes by k yields uncalibrated proportions.
This does not provide unknown-class rejection, calibrated uncertainty or
authoritative cell identities. The rare-cell failures in the external benchmark
remain relevant; matching its numerical output does not remove them.

A fit bundle contains `training/` (the native streamed analysis), `plan.json`,
`model.json`, and `receipt.json`. A map bundle includes a complete archived
`reference/`, the original query H5AD, plan, report and receipt. Verification
reconstructs the reference from its original counts and supplied labels, then
reconstructs the query result. Hashes bind source, plan, result and implementation.
Existing output directories are refused and failed staging directories removed.
Reference copying is restricted to named, size-bounded files.

Bounds: the existing H5AD reader's 1 GiB source/100 million nonzero limits apply;
reference labels are limited to 100,000 training cells and 256 classes, k to 128,
PCA components to 64. Query scores, votes and neighbor arrays have explicit size
bounds. `maximumDistanceOperations` and `maximumProjectionUpdates` default to
2 billion each (maximum 20 billion). Reported distance operations count actual
nonempty query rows; the allocation check conservatively budgets all rows.
Reference metadata/PCA scores and query scores remain resident. This is not
million-cell, approximate-neighbor, GPU or cross-study performance qualification.

`check_native.py` uses the predeclared full four-fold Baron benchmark to compare
native centers, frozen query scores, neighbor indices/distances, votes and labels
with the independent Python reference. It also exercises source reconstruction,
exact replay, label-free queries with reversed feature order, empty libraries,
contract mismatch/overlap/budget/overwrite rejection, and altered model rejection.
The external balanced-logistic baseline remains external; native multinomial
training and broader biological validation are still open.

## Product qualification, 2026-09-09

The release build passed and all 46 single-cell Swift tests in 13 suites passed.
The full Baron four-fold product workflow covers all 8,569 query cells. All
reference fits and query mappings reconstruct; repeated query reports are
byte-identical. Frozen score maximum absolute errors by donor are 2.132e-14,
2.843e-14, 2.221e-14 and 2.088e-14. Neighbor indices, candidate labels and vote
proportions match the independent scikit-learn reference exactly. The original
streamed PCA report is byte-identical when center export is omitted.

See [product evidence](evidence/2026-09-09-native/summary.json), the
[complete fit example](evidence/2026-09-09-native/human1/fit.json) and
[label-free query mapping](evidence/2026-09-09-native/human1/query.json).
The source linkage manifest binds all eight native inputs to the already archived
H5AD projection receipts and original Baron snapshot. Reports/models are gzip
compressed in the repository; local product bundles retain complete originals.
Initial fixture-basis and scoped operation-accounting failures are retained.
