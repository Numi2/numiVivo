# Native donor-response prediction

The native path fits all four predeclared baselines from a streamed H5AD training
cohort and predicts known-perturbation responses from a separate control-only
cohort. It reuses native H5AD pseudobulk aggregation and the existing bounded
symmetric eigensolver. No dense cells-by-genes matrix is created: only the small
donor-by-gene matrices and donor kernel are dense.

```sh
numivivo singlecell-perturbation-fit training.h5ad --plan fit.json --output model
numivivo singlecell-perturbation-verify model
numivivo singlecell-perturbation-predict control.h5ad --plan query.json --reference model --output predictions
numivivo singlecell-perturbation-prediction-verify predictions
```

A fit plan has schemaVersion=1, the complete existing H5AD `mapping`, an explicit
`featureNamespace`, `perturbationID`, `controlCondition`, `treatmentCondition`, and
source `provenance`. The mapping must contain only the declared two conditions.
Every training donor must have exactly one aggregated control/treatment pair with
the same biological replicate ID. Technical libraries combine through the
existing explicit replicate mapping. One organism and one cell-group identity
(including explicitly unannotated data) are supported per model. There is no
implicit selection of cells or arbitrary collapse of distinct cell groups.

A query plan contains schemaVersion=1, its own complete `mapping`, the same
`featureNamespace` and `perturbationID`. All declared samples must be controls.
Each query donor has one aggregated control profile and must be absent from
training. Organism, count unit, cell group and the complete feature-ID universe
must match; source gene order may differ. Empty libraries are rejected. No
missing-gene imputation, automatic alias mapping or treated-outcome access is
performed. Condition identity remains explicit caller-supplied assay metadata.

The versioned baseline method fixes normalization at natural-log(1+CPM), context
selection at >=10 total training counts and expression in >=2 training donor
pairs, population control scaling and alpha=1. Exactly constant control columns
use their common value as the center and a scale of one, avoiding division by a
roundoff-only standard deviation. Nonconstant columns use the training population
mean/SD. The ridge system solves `(K+I) dual = centered response` and rejects a
maximum scaled residual above 1e-9. No validation-donor hyperparameter tuning or
query refitting occurs.

Every result contains no-change, mean-response, median-response and context-ridge
estimates, each with the unclipped change, nonnegative predicted treated profile,
applied change and implied CPM sum. Feature IDs and control profiles are retained.
These are uncalibrated log-expression point estimates; they do not necessarily
form a closed CPM composition and are not raw count libraries. The previously
observed weaker folds and composition deviations remain visible. There is no
Bayesian interval, heterogeneous single-cell response distribution, causal
identification, unseen-perturbation prediction or mechanistic interpretation.

A model bundle archives the complete native streamed training bundle, the fit
plan, fitted model and hash receipt. A prediction bundle archives that reference,
the complete query H5AD, plan, predictions and receipt. Verification reconstructs
training from original counts and then reconstructs prediction. Source, plan,
model/result and implementation hashes are bound; recomputing a receipt after
altering a fitted model does not bypass reconstruction. Existing outputs are
refused and failed staging directories removed.

Bounds: existing H5AD streaming source limits apply (1 GiB and 1 billion nonzeros).
Training supports 2-64 donors, at most 2 million donor-gene values and 100 million
fit work units (`donors² × genes`). Pseudobulk materialization is limited to 128
groups and 4 million group-gene values. Encoded model/prediction documents are
limited to 64 MiB. This does not remove the outstanding million-cell,
multimodal, parallel-kernel or GPU qualification requirements.

`check_native.py` constructs each of the eleven held-out folds through native
H5AD projection, then compares fitted context transforms, dual coefficients,
all four prediction vectors and CPM sums with the predeclared external reference.
It also checks reconstruction, repeated prediction, reordered source features,
contract mismatch, donor overlap, treated query, missing gene, overwrite and
rehashed-model tampering. It needs the exact source bundles and saved external
benchmark outputs, whose source fingerprints are checked before projection.

## Product qualification, 2026-09-09

The release build and 48 single-cell tests in 14 suites passed. A subsequent
constant-control regression plus the two existing perturbation tests passed,
covering 49 distinct tests in total without changing product sources.

All eight Kang donor folds (15,706 genes) and three Hagai folds (22,048 genes)
passed source-bound fit/prediction reconstruction and exact prediction replay.
No-change, mean and median predictions match the external reference exactly.
Maximum ridge prediction errors are 1.044e-13 for Kang and 4.113e-13 for Hagai;
maximum scaled solve residuals are below 9.487e-13. These numerical checks retain
the external benchmark's biological limitations and weaker donor results.

[Kang evidence](evidence/2026-09-09-native/kang/summary.json) and
[Hagai evidence](evidence/2026-09-09-native/hagai/summary.json) retain every fold.
The source-link manifest identifies compressed originals; each fold archives its
native projection plans/receipts plus compressed model and prediction reports.
The exact source H5AD can be decompressed and projected to regenerate training
and control inputs. Verification of an archived receipt requires its qualified
implementation identity; a rebuilt runtime must regenerate evidence instead of
rewriting the old implementation hash. Fit/query examples are preserved beside
the reports, along with all expected-rejection logs.
