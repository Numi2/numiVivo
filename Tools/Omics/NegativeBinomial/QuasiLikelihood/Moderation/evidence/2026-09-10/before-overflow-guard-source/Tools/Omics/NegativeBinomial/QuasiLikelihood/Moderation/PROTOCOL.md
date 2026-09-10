# Native robust QL moderation

Freeze before new native results. Base 6cb43383fce3de26b9a48361370c3cfea29c73bf.
Use all 58 original Kang, Hagai and Crowell full-support arms. Extract exact
native abundance, final adjusted residual DF and quasi-dispersion from the
published native-abundance stage. Preserve input and upstream artifact hashes;
do not rerun count fitting or replace these data with synthetic inputs.

Implement the modern unequal-DF scaled-F prior method: precision-weighted
local linear smoothing of bias-corrected log quasi-dispersion, bounded profile
likelihood estimation of prior DF, probability-based robust reweighting,
outlier-specific prior DF and posterior quasi-dispersion. The existing equal-DF
untrended variance prior is not this model. Original reference comparisons use
limma 3.68.5 fitFDistUnequalDF1 and squeezeVar with the modern method selected.

Required numerical owners are weighted LOWESS (weight-based neighborhoods,
tricube regression, reference npts/delta geometry and constant/small-family
handling), log-minus-digamma, trigamma, stable log-gamma increments and log F
tails. Implement from mathematical definitions and primary documentation;
reference package source remains external. Native existing unweighted LOWESS
may be reused for equal precision weights, with zero robustness updates for
the prior trend. Weighted smoothing must not silently ignore prior weights.

Declare all bounds, eligible/floored observations, weight clipping, profile
boundaries, work exhaustion and unavailable fits. Keep every feature indexed.
Do not fabricate posterior values for an invalid input family. Use stable
logarithmic probabilities and preserve zero-probability handling separately
from unavailable numerical evaluation. Native family computation is CPU FP64.

Qualify controlled numerical boundaries against independent high-precision
calculations and pinned R functions. Require weighted smoother agreement within
2e-7 relative to max(1,abs(reference)), log-tail absolute error at most 2e-7,
and relative log-minus-digamma/trigamma error at most 2e-10. For profile fitting,
check native minima independently against the same-input objective, with a
per-informative-weight objective gap at most 2e-8. Include exact endpoints of
the declared shape interval [1,4999] (prior DF [2,9998]); distinguish boundary
solutions from interior fits. Reference default optimize stopping is coarse
near the upper boundary, so additionally use a tightly converged reference
optimizer with explicit endpoint evaluation. Report default-reference
hyperparameter differences descriptively rather than treating a default
stopping point as the exact optimum.

On all real arms, independently check the actual native trend inputs, tail
probabilities, robust weights, posterior arithmetic and final prior outputs.
Posterior arithmetic must agree within 2e-12 relative to max(1,abs(reference)).
Compare same-input tightened-reference posterior estimates within 2e-5 and
retain any disagreement for profile/smoother/tail arbitration; do not relax
the tolerance or drop a family to obtain a pass. Retain original default
reference differences, all attempted/failed runs and explicit boundary cases.

Use compact stage inputs and outputs, with native output persisted remotely
before compressed IPv4 transfer. Check host workload ownership before native
runs, measure a pilot, and checkpoint completed outputs by exact hashes.
This stage does not supply the constrained cohort QL hypothesis test, Poisson
bound or varying-support borrowing contract, and does not establish FDR,
coverage, power, biological truth or production calibration.
