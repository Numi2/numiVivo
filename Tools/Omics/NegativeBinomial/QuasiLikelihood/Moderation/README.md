# Native robust QL variance moderation

`VivoOmicsQLModeration.fit` estimates an abundance-dependent, unequal-residual-DF
scaled-F prior and posterior quasi-dispersions. Its `fit(from:)` entry point
consumes a complete [native abundance/global-scale fit](../Abundance/README.md),
with no supplied reference abundance or prior parameters. This completes the
native prior/posterior stage; the constrained cohort QL hypothesis test, Poisson
bound, varying-support borrowing and calibration remain open. Cohort defaults
are unchanged.

For residual quasi-dispersion `v`, residual shape `a = residualDF/2` and prior
shape `d`, the prior uses corrected log variance `log(v) + log(a) - digamma(a)`
and precision `1/trigamma(a)`. Native precision-weighted local linear smoothing
estimates its abundance trend. A bounded profile likelihood estimates `d` on
`[1,4999]`, including both exact parameter endpoints. Stable log-gamma increments
avoid subtracting two large log-gamma values in that objective.

The first fitted F distribution supplies two-sided probabilities. BH values
above 0.3 become unit weights; the remaining probabilities reweight the prior
profile. Upper-tail outliers receive individual prior DF through probability
mixing and a monotone cumulative adjustment. Posterior variance is

```text
posterior = (residualDF * v + priorDF * priorScale) / (residualDF + priorDF)
```

This is the modern probability-weighted unequal-DF method documented by
[limma](https://bioc.r-universe.dev/limma/doc/manual.html), rather than the older
Winsorized equal-DF estimator. Reference source remains in the installed package;
the native implementation owns the mathematical calculations. The qualification
script changes only the optimizer binding in an in-memory reference function,
retaining the installed function and its ordinary-default output separately.

## Explicit domain and diagnostics

The family accepts 3–100,000 rows, finite nonnegative variances up to `1e100`,
residual DF in `[0,1e6]`, and optional finite abundance coordinates. It requires
three positive informative variances. Values below `1e-12` times the positive
informative median are floored in the prior likelihood; the posterior uses the
original variance and DF. Rows with residual DF below 0.01 retain zero likelihood
weight in both prior fits. This deliberately preserves their exclusion through
robust refitting; the pinned reference's recursive call can lose its initial
low-DF exclusion. All 58 real families have DF above 1.02, so this distinction
does not affect their reference comparisons. The low-DF native contract has
separate focused coverage and is not claimed as exact reference parity.

Trend weights are normalized by their upper quartile and clipped to
`[1e-8,100]`. An all-zero raw-weight smoother or zero precision upper quartile
is unavailable. Weighted neighborhoods accumulate precision mass, with tricube
local regression and interpolation between selected anchors. Equal weights
reuse native LOWESS with zero robust updates; small families use a weighted
global line. Returned diagnostics retain clipping indices, work, anchors,
constant fits, floored/excluded indices, both profiles, evaluated objective
points, boundaries, screening probabilities, robust weights and posterior rows.

The default profile budget is 128 evaluations per profile, with a shared
100-million feature-evaluation budget. Each smoother has a separate
100-million neighborhood-visit budget. Profile minimization stops at parameter
bracket width `1e-10`. Exhausted work, invalid shapes, unavailable tails and
nonfinite results throw explicit errors. A failed family has no fabricated
posterior. `fit(from:)` rejects incomplete upstream families.

Special-function shapes are bounded to `[1e-8,1e6]`; F log-statistics to
`[-1000,1000]`. Log F tails use an incomplete-beta continued fraction with
explicit convergence checks. An underflowed probability remains distinct from
an unavailable log-tail calculation. These are native CPU FP64 calculations.

## Qualification on real experimental data

All **58** original Kang, Hagai and Crowell arms pass, retaining all **483,576**
gene/arm rows from the published abundance stage. No count fit is rerun or
substituted. The prior performs 111 profiles, including 53 robust refits;
20 profiles select the upper boundary and 91 select interior solutions.

| Independent comparison | Maximum error | Frozen limit |
| --- | ---: | ---: |
| Precision-weighted trend | 1.59e-13 | 2e-7 |
| Profile objective value per informative weight | 7.50e-12 | 2e-8 |
| Objective gap to tightly optimized same-input minimum | 8.78e-14 | 2e-8 |
| Native log upper-tail probability | 1.70e-12 | 2e-7 |
| Robust screening/refit weights | 1.16e-13 | 2e-7 |
| Posterior arithmetic | 4.43e-16 | 2e-12 |
| Posterior against tightened reference with endpoints | 1.32e-6 | 2e-5 |

Trend, weight and posterior errors divide by `max(1, abs(reference))`; log-tail
errors are absolute. The reference pins limma 3.68.5 and edgeR 4.10.5. Default R
optimization can stop before the upper endpoint: ordinary-default posterior
discrepancy reaches **0.0962** under the same error normalization, and common
prior DF discrepancy reaches **0.2291**. These differences are retained, not
described as exact default-reference agreement. Independent objective checks
and tightly optimized reference endpoints justify the native parameter choice;
they do not establish better biological calibration.

All **265** numerical cases pass: 25 weighted smoothers, 15 shape cases checked
at 80-digit precision, and 225 log-tail cases. Maximum errors are 2.66e-12 for
smoothing, 4.45e-16 relative for shape functions and 4.66e-10 absolute for log
tails. All **43** focused count-model tests in eight suites pass on the physical
M4 Pro with Swift 6.3.3, including lower/upper endpoint and arithmetic-overflow
assertions on the final sources. Tests cover robustness, original-value posterior
arithmetic, low-DF exclusion, work failures and native abundance integration.

The native executable SHA256 is
`48f4f13238ea8b2f350ccec383903ef827c8b3db7479cc00b9e82ea2f32946bf`.
The largest final process resident size is 29,343,744 bytes; profiles visit
47,414,224 feature terms across all arms. These measurements concern the compact
moderation stage only and are not end-to-end speed or million-cell evidence.

An initial special-function build failed before producing an executable; its
source and errors are retained. The first R pilot checker failed while closing
an already-closed file connection, before comparisons. The repaired reader
reused the same native output. All native results were spooled on the Mac mini
and transferred as compressed files over IPv4.

A final boundary audit found that a small-family weighted line could return
nonfinite values after intermediate overflow. A retained native before/after
probe reproduces that result and the repaired explicit error. All 58 real-data
arms, 265 numerical checks and 43 tests were rerun on the final binary; every
real-data fit and feature identity remains exactly equal to the earlier candidate.
Both candidate families and their source identities are retained.

## Evidence and reproduction

The [summary](evidence/2026-09-10/summary.json) and
[manifest](evidence/2026-09-10/manifest.json) retain complete compact inputs,
native outputs, reference diagnostics, hashes, tests and failed attempts.
Upstream abundance/count matrices and native executables are externally
referenced. The archive verifier validates stored members, not remote file
availability.

```sh
bash build.sh NATIVE_RUN/build
bash build_numerics.sh NATIVE_RUN/numerics-build
python prepare.py --abundance-root ABUNDANCE_RUN --ql-root QL_RUN --out RUN/inputs
python check_numerics.py --out RUN/numerics-attempt-2 --binary NATIVE_RUN/numerics-build/numerics
python run.py --inputs RUN/inputs --out RUN/family --remote-out NATIVE_RUN/family --binary NATIVE_RUN/build/nb-ql-moderation --limit 1
python run.py --inputs RUN/inputs --out RUN/family --remote-out NATIVE_RUN/family --binary NATIVE_RUN/build/nb-ql-moderation
python summarize.py --root RUN
python archive.py --root RUN --out ARCHIVE
python archive.py --out ARCHIVE --verify
```

The published final summary selects `family-final`, `numerics-final-check` and
`swift-tests-final.log` using the corresponding `summarize.py` options. Its
`--comparison-family-directory family` check verifies unchanged real-data fits
against the retained candidate before the overflow guard.

Use new directories for failed-attempt recovery, retaining earlier artifacts.
Completed receipts are reused only with exact input, binary, protocol, checker
and output hashes. Collect native source/build/test logs and `execution.json`
before summary/archive generation. The frozen [protocol](PROTOCOL.md) fixes all
arms and limits. Remaining work includes constrained QL cohort hypotheses,
Poisson bounds and fresh sham/held-out calibration. The parent study's 39 Hagai
sham calls versus Wald's 20 remain unresolved evidence.
