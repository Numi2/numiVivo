# GSE181897 external IFN-beta prediction protocol

Declared 2026-09-11 before reading this cohort's count matrix, fitting models or
inspecting prediction scores. Dataset discovery and primary GEO experimental
metadata have been read. Earlier Kang, HIRISA and their transfer/interval results
have been inspected; model selection here is therefore informed by those
failures. GSE181897 is a new test collection for this project, not an assertion
of independent investigators or verified absence of participant overlap.

## Admission and complete cohort

Use the original GEO `GSE181897_concat.4.raw.h5ad.gz`, explicitly described by
its submitter as raw integer RNA and protein counts with individual covariates.
Retain the source bytes and verify gzip integrity and decoded identity. Inspect
metadata before count aggregation. Source treatment is 500 IU/mL IFN-beta for
9 hours, with a same-incubation unstimulated control. GEO reports 64 individuals,
five stimuli and twelve pooled sequencing reactions. Pools are technical units,
not independent biological donors. Preserve author donor, treatment and cell
identity fields; do not infer these from outcome expression or barcode patterns.
If essential mappings cannot be established, retain the admission failure and
resolve it before prediction.

The endpoint follows the prior B-cell study question. Include every released,
author-annotated B cell in control or IFN-beta, for every donor with both
conditions. List each unmatched donor, ambiguous label and excluded category;
do not select donors or genes by prediction error, response strength or cell
count. Require at least one cell and positive RNA library count in each paired
stratum; report sparse strata. Do not redefine B cells from the treated outcome
with a newly fitted classifier. Author treated-cell labels make this an
annotation-conditioned endpoint, not prospective cell identification.

Separate RNA from surface-protein features using primary feature metadata.
Verify finite, nonnegative integer counts over the complete source matrix in
bounded sparse chunks; no dense cells-by-features array. Preserve all source
RNA features in each donor's library denominator. The prediction panel is the
intersection of the previous 11,884-symbol Kang-HIRISA panel with exact unique
source RNA symbols in this cohort. No alias guessing, zero padding or post-score
feature filtering. Freeze the full mapping and unsupported features. Feature
metadata may determine this intersection before counts are aggregated.

## Frozen predictors and endpoints

Train separately on all eight original Kang donor pairs and all five original
HIRISA enriched-B-cell donor pairs. Reuse their source-qualified count bundles
and full original normalization denominators. Fit the existing native no-change,
mean, median and context-ridge methods, with unchanged ridge settings, on the
metadata-defined panel. Request the existing mean-response nominal 95% Student-t
future-donor intervals. There is no new training, calibration, hyperparameter
choice, variance floor or uncertainty correction using GSE181897. All eligible
GSE181897 donors are final query/test donors. Query inputs contain control counts
only; treated counts enter scoring after native predictions are frozen.

Report both training origins separately, for every eligible donor. The primary
point-prediction criterion is a mean-response reduction of at least 5% in
equal-donor all-panel RMSE versus no change in each origin separately. This is a
prospective useful-improvement threshold for this experiment, not an established
clinical criterion. Retain the prior ridge comparison as secondary: ridge must
beat both no change and the training mean. Do not choose a winning origin or
method after scores are available. RMSE is on treated log1p(CPM), equivalently
the response after applying the same observed query-control offset; keep the
existing nonnegative clipping and report its impact.

For uncertainty report pointwise coverage, lower/upper misses, mean width and
available/total genes per donor, in both raw response and transformed treated
spaces. Genes, cells and reused training sets are not independent replicates.
Report equal-donor means and donor ranges, including cases worse than no change.
The nominal 95% intervals are not called calibrated merely because their mean
coverage approaches 95%; they also assume independent normal donor responses
and a fixed observed query control. No simultaneous-gene or ridge intervals are
claimed. No result qualifies phenotype, tissue function or clinical benefit.

## Execution and evidence

Freeze source identities, cohort admission, feature mapping, all plans and
query controls before native fitting. Execute both complete query cohorts,
verify native models and predictions, and freeze outputs before scoring.
Independently reconstruct all point estimates and interval bounds, including
full-library normalization, training-only scaling and variance. Repeat scoring
for identical result bytes. Keep all failures, unsupported donors/features,
source/executable versions and complete restorable results.

Use indexed reads of the original gzip H5AD to avoid a redundant 3 GB decoded
copy on a nearly full disk. This is external source preparation; it is not a
claim that the native product directly reads gzip H5AD. Native H5AD pseudobulk
transport and native model execution must be reported separately from original
single-cell ingestion. Any follow-up tuning after this test requires new
validation before a stronger generalization claim.
