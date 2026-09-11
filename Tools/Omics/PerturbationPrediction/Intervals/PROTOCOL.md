# Donor-response predictive-interval protocol

Declared 2026-09-11 before fitting or examining interval coverage. The Kang,
HIRISA and cross-study point-prediction results have already been inspected.
This is native uncertainty-method development and empirical coverage assessment
on reused studies, not new independent biological validation.

Implement one explicitly requested pointwise interval for a future donor's
response under an independent, identically distributed normal donor-response
model with unknown mean and variance. Center it on the training mean response,
not context ridge. For n paired training donors, use sample response variance
s² with n−1 degrees of freedom and half-width
`t((1+coverage)/2, n−1) × s × sqrt(1+1/n)`.
The extra future-observation variance distinguishes prediction from confidence
in the mean. Reuse the native Student-t numerical owner. See
[NIST prediction limits](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/predlimi.htm),
with one future observation.

Zero or exactly constant training-response variance yields an unavailable
interval, not a zero-width assertion of certainty. Retain those genes and their
missingness. Do not fit a variance floor, shrinkage prior, coverage correction,
feature filter or threshold from test outcomes. A new biological donor is the
prediction unit; genes and cells are not independent donors.

Use the already frozen Kang–HIRISA shared-gene panel, source-normalization
conventions and all 26 cross/within folds. Add only nominal coverage 0.95 to
native fit plans, freeze those plans, then complete every native fit/prediction
and replay before interval scoring. Keep all four point predictions unchanged.
Preserve unclipped response limits and separately transformed nonnegative
predicted-treated limits. The observed query control is treated as fixed;
measurement error and query-specific heteroscedasticity are not separately
modeled. Clipping can create a boundary point mass, so report coverage both
before and after that transformation.

For each donor, report available/total genes, fraction covered among available
genes, misses below/above, and mean interval width in log1p(CPM) units. Report
clipped and unclipped coverage separately. Summarize by equal donor weighting
within query study and cross/within mode; retain every donor and unavailable
feature. Do not pool genes into a binomial confidence interval, treat overlapping
folds as independent studies, or call nominal 95% empirically calibrated merely
because arithmetic/replay checks pass. Report deviations from 95% without tuning.
No calibrated uncertainty or production promotion is an acceptance outcome of
this small, reused-cohort experiment.

Independently verify sample variances, Student-t quantiles and all four interval
bound arrays with NumPy/SciPy, plus unchanged point predictions against the
previous frozen native results. Repeated scoring must be identical. Test invalid
coverage, constant-response unavailability, future-donor versus mean uncertainty,
query clipping and unchanged default serialization. Keep failed attempts and
all executable/source identities. This does not attach uncertainty to the
kinetic target-engagement model or qualify a cross-scale Bayesian coupling.
