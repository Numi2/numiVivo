# Kang–HIRISA IFN-beta transfer protocol

Declared 2026-09-11 before cross-study model fitting or scoring. Both studies and
within-study outcomes have already been inspected. This is a new transfer
experiment on reused studies, not untouched external or prospective validation.
No predictor tuning, winner selection or calibrated uncertainty claim is planned.

Use all eight Kang B-cell donor pairs and all five HIRISA enriched-Bcell IFNb
pairs from the previously source-qualified count aggregates. Retain every source
gene for each library's normalization denominator and all cells in these admitted
cohorts. Do not equate enriched B-cell preparations with pure annotated B cells.

Train on all five HIRISA donors, predict the observed controls of all eight Kang
donors; train on all eight Kang donors, predict the observed controls of all five
HIRISA donors. Also run every within-study leave-one-donor-out reference on the
same output panel. Test donors' treated counts enter scoring only. Inputs,
training memberships, feature mapping and plans are frozen before fitting;
all predictions are frozen before scores are computed.

Kang uses IFN-beta-treated PBMCs from lupus donors for six hours. HIRISA uses
healthy-donor enriched B cells, IFN-beta 100 units/mL for 21 hours and a fixed-RNA
probe assay. Timing, preparation, health status, library chemistry and study vary
together. Interpret transfer as a combined context shift, not an isolated causal
effect or interchangeable assay measurement. Source author annotations condition
the outcome strata; prospective identification of treated cells is not qualified.

The shared panel contains only exact, unique source gene-symbol matches. Preserve
HIRISA's original Ensembl IDs and symbols in the mapping, and Kang's original
symbol IDs. Do not guess aliases, collapse ambiguous names, strip suffixes or
pad absent genes with zero. Exclude unmatched/ambiguous features only from the
output/context panel, never from either full-source library denominator. Retain
all exclusions and quantify coverage. This is symbol-matched transfer, not a
verified common reference genome/annotation release.

Use the unchanged four numerical baselines: no-change, training mean, training
median, and training-control-context ridge with alpha=1. Context features are
chosen only within the fixed shared panel using >=10 training counts and >=2
expressing donor pairs. Standardize using training control population moments,
scale by sqrt(context width), use a training mean-response intercept. Normalize
with natural log1p(CPM) before selecting the panel. Clip predicted treated values
below zero; retain unclipped changes and panel implied-CPM subtotal without
reclosing. Treat the query control as observed, not an inferred absolute RNA level.

Primary endpoint: equally weighted donor mean response RMSE across every shared
feature, separately for each transfer direction. The fixed ridge comparison
passes a direction only if it is strictly better than both no-change and the
cross-study training-mean response. Report relative gains without interpreting
an arbitrarily small difference as practically sufficient. Report all per-donor
RMSE/MAE/correlations, worse-than-baseline counts, and matched within-study errors.
Do not pool the two directions or count overlapping folds as independent studies.
The RNA endpoint does not establish disease prediction or calibrated intervals.

Native execution must perform every fit/prediction and reconstruct its bundle.
An independent NumPy/scikit-learn calculation must verify normalization, selected
features, model coefficients and all prediction vectors. Repeated scoring must
be identical. Test that missing panel genes reject, extra nonpanel counts affect
the denominator, and existing strict full-universe predictions retain their
behavior. Keep historical biological results bound to their original executable.
