# Fixed-dispersion contrast likelihood-ratio qualification

Freeze before evaluating LRT results. Use all forty original untreated-cell
Kang/Hagai null analyses (ten splits, two support policies), plus all ten
previous empirical-prior treatment analyses spanning Kang, Hagai and Crowell.
These data have been inspected before; this is method sensitivity, not an
untouched calibration experiment. No count, donor, filter, normalization,
dispersion, influence or support policy is selected using the new LRT scores.

Recompute each original request with the current production cohort owner and
require exact original feature results, design and dispersion diagnostics.
Then change only `negativeBinomialOptions.testMethod` to `likelihoodRatio`.
The same unpenalized full-model fits must be preserved. Refit nuisance
coefficients subject to the requested linear contrast equaling zero, at the
same gene-specific final dispersion. No extra pseudocount or effect prior
enters the test. Retain null coefficients, means, convergence, signed raw
statistic and failures. Reject negative LR below -1e-7; clamp only smaller
roundoff to zero. Use chi-square(1) tails and BH over available LRT genes.
Boundary, replication, support and influence exclusions remain in force.
Null convergence failures withhold a test; no silent Wald fallback.

Independently check all admitted genes with NumPy/SciPy: source identities and
counts; exact inherited design/dispersion; zero constrained contrast and
null-space score equations (scaled score <= 1.1e-7); relative null mean
reconstruction <= 1e-10; LR statistic using long-double likelihood differences
within 2e-6; chi-square tail and BH within 2e-7. Independently fit the same
full/null models with pinned edgeR glmFit/glmLRT, prior.count=0 and fixed
native dispersions/offsets, for every available gene of the ten treatment
analyses. Accept LR differences <= 2e-5 and probabilities <= 2e-6; retain
package convergence/messages and every mismatch. No new dispersion fit is
presented as an isolated LRT comparison.

Report all original and new method-specific tested families, exclusions,
BH<0.05 calls and changes. Overlapping null splits do not constitute new donor
replication. Treatment calls are not effect truth or power. This asymptotic
fixed-dispersion LRT does not incorporate dispersion uncertainty and is not
quasi-likelihood or evidence of general FDR control. No automatic default
promotion follows from fewer sham calls or agreement with a reference.

Qualify the production Swift build and focused tests on the physical Mac mini.
Exercise H5AD publication/replay for a real paired contrast and verify an old
Wald plan produces the same report with the new binary. The compact evaluation
harness compiles the actual production owners and is only measurement tooling.
