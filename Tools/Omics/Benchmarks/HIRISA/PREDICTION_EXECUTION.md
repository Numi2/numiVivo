# Frozen HIRISA donor-held-out prediction execution

2026-09-10, before any HIRISA response-model fitting. Source preparation, paired
DE fitting and some DE concordance results were already inspected. Those results
do not select the model, genes, donors or treatments below. The original
PROTOCOL.md already prespecified leave-one-donor-out prediction for all sixteen
enriched contrasts. This specification fixes the existing native owner's
settings and scoring without changing that cohort.

Use all **79 folds**: five held-out donors in fifteen contrasts and four in the
predeclared Bcell IFNg contrast. The latter trains on three paired donors; all
other folds train on four. Keep the original accession, donor, matched batch,
condition and complete gene identities. The training input contains only the
eligible training donor pairs. The query contains only the held-out donor's
matched control library. Its treated library is a scoring target and must not
be an input to fitting, feature selection, normalization, query projection,
regularization or prediction. Preserve every fold failure.

Use the existing `VivoPerturbation` owner, method
`paired-donor-log1p-CPM-response-baselines-alpha1-v1`, without outcome-driven
changes. Sum original integer counts per library and normalize each library by
its own total to log1p(CPM). Training donors have equal weight. Context features
require at least ten training counts and expression in at least two training
donor pairs. Centers and population standard deviations use training controls
only; constant controls use scale one. Divide standardized contexts by the
square root of the selected feature count. The dual ridge system adds exactly
one to its diagonal. Do not tune alpha, select a DE signature, use author cell
labels to filter the source, or choose between methods from held-out outcomes.

Retain every existing estimate: noChange, meanResponse, medianResponse and
contextRidge. Predicted treated log1p(CPM) is clipped at zero, exactly as in the
owner. Report unclipped and applied responses and the implied CPM sum; do not
silently renormalize predicted expression or describe it as raw counts.

Freeze model and query outputs before opening the corresponding treated target
for scoring. Check all native coefficients, selected features, context moments,
ridge solve and point predictions against an independent NumPy reconstruction
using the same training/query counts. The target cannot be used for this
reconstruction. Compare each method on two explicit feature families: the full
18,082-gene source and that fold's training-selected context genes. Within each
family, report treated-expression RMSE/MAE and response RMSE/MAE/Pearson
correlation against the actual held-out treated-minus-control log1p(CPM).
Correlation is undefined for a constant prediction and must remain null.
Report fold-level results and equally weighted donor-fold summaries within each
contrast. No post-hoc gene selection, successful-fold-only summary, pooled-cell
replicates or selective treatment reporting.

The fold manifest binds source/design/protocol identities and the current model
owner source hashes. Transport must retain original cell membership and exact
library aggregates. Native raw-ingestion verification is required before the
inference and prediction chains can claim native end-to-end linkage. Reusing
verified aggregates must not give the model access to held-out treated rows.
The existing publication route repeatedly scans training source snapshots; any
reusable aggregate route needs its own count, membership, donor-exclusion and
replay qualification before use.

This is prediction of a known perturbation in an unseen donor, not prediction
of a new perturbation identity or biological preparation. The PBMC transfer
mapping remains a separate, unfitted specification. No calibrated uncertainty,
causal mechanism or general biological validity follows from these folds.
