# Frozen HIRISA preparation-transfer experiment

Declared after the enriched-population donor-response results and integration
diagnostics were inspected, at source revision `5d2d7bfd`. No preparation-transfer
models or PBMC response scores have been computed for this specification.
It is a new held-out endpoint within a previously inspected study, not an
independent-study replication or a prospectively collected experiment.

## Question and complete cohort

Can a response model trained on enriched immune cells predict the average IFNα
RNA response in the corresponding author-labelled lineage cultured within whole
PBMCs? Also test the reverse direction. Both directions withhold the query donor
from every training preparation. Use all five original donors, all matched
IFNα/unstimulated libraries and all eligible cells, retaining the complete
1,612,594-cell source and a partition ledger for every source cell.

The original [experimental methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/methods/)
describe enriched and whole-PBMC cultures and the same 21-hour IFNα treatment.
Use the exact GEO metadata for donor, experiment batch and pool matching.
`culture_IFNa` maps explicitly to IFNa and `culture_no_stim` to none for the PBMC
arm. Fresh samples are not matched cultured controls. Enriched controls retain
the original same-batch, same-pool matching; do not pool the separate T-cell
control batches. Preserve the original condition labels in source provenance.

## Fixed annotation mapping

Use exact deposited `celltype.l1` labels in both preparations:

| Lineage ID | Exact author label | Matching enriched preparation |
| --- | --- | --- |
| B | B | Bcell |
| Mono | Mono | Monocyte |
| NK | NK | NK |
| CD4-T | CD4 T | Tcell |
| CD8-T | CD8 T | Tcell |
| other-T | other T | Tcell |

Do not merge T subtypes, use confidence-score cutoffs, relabel cells, select DE
genes or train a new annotator. Author `DC` and `other` cells and labels outside
their matching enriched preparation remain in the exclusion ledger. All source
libraries and cells remain preserved. The author's [analysis methods](https://apps.allenimmunology.org/aifi/resources/ifn-response/analysis/)
describe reference label transfer and subsequent curation. These deposited labels
are predictions, not independently verified identities. Treated-cell labels
define the evaluation stratum and may depend on treatment; this is conditional
population-average expression prediction, not prospective cell identification.

## Pool matching and donor independence

First pair every cultured-PBMC treatment library with its control by original
donor, experiment batch and pool. Combine raw counts across those matched pool
pairs only within the same donor, preparation, lineage and biological condition.
Retain every contributing accession and batch/pool identity. Never turn pools
or technical libraries into additional donors. One donor consequently contributes
one control and one treated aggregate in each preparation and lineage.

Select all cells satisfying the fixed library and label mapping. A usable donor
pair requires at least ten cells in each aggregate and positive count libraries.
Membership and cell-count eligibility are fixed before reading response values;
zero-count libraries become explicit execution failures. Keep all five requested
held-out folds per lineage and direction, including failed/unavailable folds.
Use every eligible other donor for training, requiring at least two training
pairs under the existing native owner. Never replace a failed fold with a
different donor, pool, lineage or preparation.

## Model and comparisons

Use the unchanged native `paired-donor-log1p-CPM-response-baselines-alpha1-v1`
owner and [previously qualified numerical settings](PREDICTION_EXECUTION.md):
full-feature library totals, log1p(CPM), training-only feature selection and
control scaling, equal donor weights, ridge alpha 1, zero clipping and no CPM
reclosure. Retain noChange, meanResponse, medianResponse and contextRidge.
Neither PCA nor integrated embeddings enter fitting or prediction.

Freeze sixty cross-preparation folds: six lineages × five held-out donors × two
directions. For each, add a matched within-preparation reference trained on the
other donors in the query preparation, yielding sixty reference folds. Query
controls and scoring targets are identical between each pair of models. The
reference quantifies a preparation-transfer penalty; it does not select a model
or a preferred direction. The same donor is excluded from both training sets.

Report all twelve lineage/direction contrasts. Primary endpoint is equally
weighted donor-fold mean all-18,082-gene response RMSE, accompanied by every
fold result and coverage. A cross-preparation ridge advantage requires lower
RMSE than both noChange and its cross-preparation training-mean baseline in
that contrast. Missing folds prevent an unqualified complete-contrast claim.
Report cross-minus-within errors for every learned method without tuning to
that comparison. Retain the training-selected feature family as secondary;
it cannot rescue an all-gene failure. Include MAE, response correlation, clipped
predictions and implied CPM totals as in the prior experiment. No calibrated
intervals or significance claims follow from five donors.

## Execution and verification gates

Freeze exact source/protocol/annotation membership, matched libraries, aggregate
keys, folds and model-owner hashes before preparation-transfer fitting. Store
the training, query-control and scoring-target roles separately. Freeze every
prediction before scoring treated outcomes. The fitter receives only selected
training pairs; the predictor receives only held-out controls.

Native aggregation must reconstruct original sparse counts for every selected
cell and gene. Check original membership and integer aggregates independently
against source HDF5/SciPy; no dense cell-by-gene representation. Check donor and
source-cell exclusion, original feature identities, all native models and
predictions against an independent NumPy implementation, and native replay.
Mutating excluded/held-out treated counts must not change training or query
inputs. Preserve failures, unsupported strata and exact runtime identities.

Preparation is also associated with experiment batch and culture composition.
This design measures transfer across those combined contexts; it cannot isolate
cell-cell effects, technical batch, enrichment or other causal mechanisms.
It does not qualify unseen perturbations, tissues, species, patients, individual
cell distributions, clinical outcomes, Bayesian coupling or Metal performance.
