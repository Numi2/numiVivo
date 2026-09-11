# Exposure-duration model development: GSE226572

Frozen before the first duration-model fit. The fixed Kang-to-GSE226572 experiment
has already been scored and inspected. This is retrospective model development,
not a new independent external validation. Preserve its failed 5% gate.

## Cohort and identity

Reuse the exact complete source-verified GSE226572 aggregates from the
[external experiment](../GSE226572/README.md): all 24 raw libraries, 126,633 cells
admitted by initial author QC, 36,601 measured genes, three literal donors,
two pooled culture-matched zero-hour controls per donor and six treated times
per donor. These are whole-population RNA outcomes, not authoritative cell types.
All cultures lasted 36 hours; IFN-beta1a exposure was staggered, 1,000 U/mL.
Use the unchanged 12,993 unique shared-symbol response panel from that experiment.
All 36,601 features remain in every library-normalization denominator.
No new QC, barcode selection, relabeling, variable-gene optimization or tuning.

## Fixed model and split

Run exactly three leave-one-donor-out folds. Fit the other two donors' control
and all six treated aggregates; present only the held-out donor's pooled control
to prediction. Predict all six of that donor's actual exposure hours.

For each training donor, form treatment-minus-control log1p(CPM) vectors.
Include a zero-hour, zero-response anchor. Interpolate each donor's curve linearly
in log1p(hours), returning exact vectors at observed knots. Refuse extrapolation
past any training donor's last time. Average the interpolated responses with equal
donor weight; also retain their median. Context ridge uses alpha=1, training-only
control centering/scaling, expressed-in-two-donors and total-count-at-least-10
features, and the existing normalized linear kernel. Do not tune alpha, time
transformation, panel, normalization or feature thresholds after scoring.

Clip predicted treated logRNA at zero. Retain unclipped and applied responses.
Nominal 95% future-donor Student-t bounds use the two interpolated donor responses
(df=1), sample variance and sqrt(1+1/n). Exactly constant responses have unavailable
bounds. These assumptions do not establish calibration; interpolation error,
control uncertainty, dose transfer and causal or mechanistic kinetics are absent.

## Comparisons and declared gates

Primary metric: per-outcome RMSE across every panel feature, averaged equally
across the six times within each donor, then equally across all three donors.
Compare durationMean with noChange and a matched time-invariant training mean.
The latter averages the six observed nonzero responses within each training donor,
then averages the two donor means, and is applied unchanged at every query time.
This matched baseline uses the same training donors and features, and ignores time.

Primary development gate: durationMean reduces primary RMSE by at least 5% against
BOTH noChange and matched time-invariant mean. Secondary: durationContextRidge
must beat BOTH durationMean and noChange. Report every method, donor and time,
including negative results; do not choose the best time or method post hoc.
The prior eight-Kang-donor fixed-response results are a separate training-origin
comparison and do not isolate the effect of the time model.

## Execution and independent arithmetic

Fit and predict with the actual Swift owners and CLI. Verify every model and
prediction by reconstruction from frozen H5AD snapshots. Repeat predictions and
compare report bytes. Freeze all native predictions before this experiment's
scoring. Independently reconstruct full-source normalization, per-donor curves,
log-time interpolation, context selection/scaling, NumPy linear solve, clipped
points and SciPy Student-t bounds. Check every coordinate at absolute tolerance
1e-8 / relative tolerance 1e-8 (intervals 1e-8). Retain all raw predictions,
source dependencies, commands, exit codes, executable/source hashes and metrics.
Repeat scoring and require exact result bytes. No subsampling or cherry-picking.

Passing would qualify a bounded held-out-donor RNA development result in this
study. It would not demonstrate unseen perturbations, external tissue/lab transfer,
cell identity, protein/immune outcomes, phenotype, treatment benefit or clinical use.
