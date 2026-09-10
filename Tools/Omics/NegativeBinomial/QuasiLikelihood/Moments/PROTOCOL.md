# Native NB deviance moments and adjusted residuals

Freeze before new native/reference moment results. The target is native QL
residual adjustment, followed by global scale estimation, robust unequal-DF
prior fitting and cohort inference. This stage does not claim the latter steps.

Implement direct, mode-centered NB/Poisson probability recurrence, accumulating
mass and the first two raw moments of the previously qualified unit deviance.
Bound omitted left and right probability and deviance moments using geometric
ratio bounds. Normalize retained weights and propagate truncation bounds to
mean and variance. Default relative tolerance 1e-10; no sampled or unbounded
implicit tail cutoff, fitted Chebyshev tables or asymptotic fallback. Return
explicit failure if at most 1,000,000 evaluated support points cannot establish
the requested bounds. Inputs outside representable arithmetic fail explicitly.
Bounds cover omitted probability, not floating-point interval rounding.

For observations with mean mu, NB trend dispersion phi and supplied average
QL scale s, evaluate NB deviance moments at mean mu/s and dispersion phi.
Let m=E[D], v=Var[D], A=2m/v and K=2m^2/v. Fit-information weights are
mu/(1+mu*phi/s); native weighted QR gives leverages h. Adjust the observed
unit deviance at dispersion phi/s by A, and residual DF by (1-h)*K. Preserve
all moment/tail/convergence diagnostics. Match the reference's declared
residual-leverage threshold 1e-4: below it both contributions are zero.
Do not fit a new mean, dispersion, scale or prior in this conditional stage.
Observation weights and heterogeneous active-donor prior borrowing are outside
this initial contract and remain required future work.

Freeze a Cartesian moment grid: mean {1e-8,1e-4,0.01,0.1,0.5,1,3,10,30,100,1000}
x dispersion {0,1e-8,0.001,0.1,0.7,1,4}. Independently evaluate moments using
SciPy NB/Poisson PMFs and high precision checks for low means. Compare native
means, variances, A and K with relative tolerance 2e-7 and retain all failures.
Include invalid/degenerate inputs, deliberate work exhaustion, increasing
precision, reparameterized weighted designs and unit-leverage observations in
native tests.

Apply the native residual adjustment to both modern QL arms of all 29 source
cases from the published QL stage study, with the original saved counts, means,
designs, trend dispersions and average QL scales. Record every attempted gene,
all errors and support work. Compare native leverages to reference within 1e-9
absolute. Compare adjusted deviances and DF descriptively, retaining differences
from the reference approximations; numerical equality is not assumed. Use
independent direct-moment checks to arbitrate disagreements, never retune toward
package output. No new inference, BH decisions or calibration claims are made.

Execution may be checkpointed by case, with hashes and terminal receipts, but
full-family qualification requires all 58 arms. Measure the first case before
scheduling the full computation; if direct evaluation is too costly, retain
that measured limitation and use this engine to qualify a future accelerated
moment evaluator, rather than claiming full-family coverage from a subset.
