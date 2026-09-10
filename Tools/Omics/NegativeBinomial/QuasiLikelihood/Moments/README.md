# Native deviance moments and adjusted QL residuals

`VivoOmicsNBResidualAdjustment` now owns conditional NB/Poisson deviance moments
and adjusted QL residuals in Swift. It uses the previously qualified native
unit-deviance calculation. It does not estimate the global QL scale, fit a
robust prior, or perform a QL cohort test. Existing Wald/LRT plans remain unchanged.

## Method

`moments(mean:dispersion:relativeTolerance:maximumTerms:)` starts at the count
distribution's mode and accumulates relative probability weights and the first
two deviance moments in both directions. Compensated sums and final weight
normalization avoid underflow from starting with P(Y=0) at large means. No
reference-package coefficient tables, empirical approximation or Gamma fallback
enters the native calculation.

Success requires explicit bounds on omitted probability, mean and variance.
The NB right tail uses a geometric bound on successive probability ratios and
a linear upper envelope for deviance; the Poisson envelope is quadratic. The
left tail uses its decreasing reverse probability ratios and zero-count
deviance as an upper bound. First/second-moment tail bounds are propagated
through normalization and variance subtraction. These bound truncation, not
floating-point interval rounding. The default tolerance is 1e-10 and the default
limit is one million evaluated counts per moment. Work exhaustion is an error,
with no substituted asymptotic estimate. The kernel retains a fixed number of
scalars while summing support points; returned observation diagnostics and the
measurement harness still consume memory.

For supplied fitted means mu, NB trend dispersion phi and average QL scale s,
`adjustedResiduals` evaluates moments at mean mu/s and dispersion phi. With
m=E[D] and v=Var[D], it uses A=2m/v and K=2m²/v. The observed deviance at
variance dispersion phi/s is multiplied by A. Native weighted QR of the
supplied design gives leverages h; adjusted DF is (1-h)K. Residual leverage
below the declared 1e-4 threshold contributes zero to both sums. This follows
the moment-matching principle described by
[Chen et al., edgeR v4](https://doi.org/10.1093/nar/gkaf018), evaluated here by
direct probability summation. No new means, trends, average scales or priors are
fitted by this conditional operation.

Zero-mean standalone distributions return zero information. Conditional
scaling that underflows a positive mean to zero fails explicitly. This was a
real boundary defect in the first implementation: the failed probe and old
owner source are retained. The repair also checks transformed dispersion and
aggregate representability. All 58 real-data arm outputs are byte-identical
between the original and repaired binaries. Thirty-seven independent check
receipts were reused only after verifying exact output equality; the remaining
checks ran on the repaired output.

## Measured evidence, 2026-09-10

The frozen [protocol](PROTOCOL.md) uses the two modern QL reference arms for
all 29 full-support Kang, Hagai and Crowell analyses from the
[parent stage study](../README.md). Counts, means, design, trend dispersions and
average QL scales are bound to those original input/reference hashes. Active-
donor prior borrowing and observation weights are outside this stage.

- All 25 focused Swift tests in four suites pass on the physical M4 Pro.
- The complete 77-point NB/Poisson grid passes independent PMF and 80-digit
  checks, with maximum relative error 7.97e-11.
- All 58 real-data arms were attempted. Fifty-four have complete default-budget
  output. Twenty gene/arm fits in the other four exhaust the declared limit.
- All 3,940,848 available default-budget moments pass independent, bounded-batch
  PMF summation. Maximum relative moment/scale/DF error is 6.08e-10; maximum
  PMF mass discrepancy is 2.85e-9. No high-precision rescue was needed.
- Native/reference leverage error is at most 3.89e-15. Relative differences
  from edgeR's approximate adjusted deviances and DF reach 0.001212 and
  0.000626 respectively; these descriptive differences are retained rather
  than treated as exact package agreement.

The completed default moments required 2,135,439,271 count-support evaluations.
Native batch-process peak RSS ranged from 159.4 to 279.4 MB. These measurements
include the JSON harness and shared-host execution; they establish neither
performance superiority nor million-cell performance.

## Work-limit failures remain

The primary failures are nine Hagai treatment genes under each trend and one
Crowell excitatory-neuron gene under each trend. A separately declared
[ten-million-point sensitivity](WORK_LIMIT_SENSITIVITY.md) preserves all twenty
failures and reruns their exact inputs through the same native owner with only
its public work bound increased.

Nineteen fits then complete, adding 118 independently checked moments. Maximum
relative moment error is 1.57e-9, and the largest completed moment uses
9,384,448 support evaluations. Hagai feature index 19389 under the edgeR trend
still exhausts ten million evaluations. Its entire six-observation residual
result remains unavailable. The sensitivity does not change the default limit
or replace failed primary receipts. This demonstrates the need for an
accelerated evaluator qualified against the direct moment engine before
full-family native QL integration.

The study is deliberately `completed-with-resource-limit-failures`. It is not
a native QL cohort method, posterior calibration, power evidence or FDR-control
qualification. The earlier QL sham-discovery concern is unchanged.

## Reproduce

```sh
bash build.sh NATIVE_RUN/build
python check_grid.py --out RUN --binary NATIVE_RUN/build/nb-moments
python run.py --ql-root PRIOR_QL_RUN --out RUN --binary NATIVE_RUN/build/nb-moments --jobs 3
python check_real.py --ql-root PRIOR_QL_RUN --root RUN
bash build_sensitivity.sh NATIVE_RUN/sensitivity-build
python check_sensitivity.py --ql-root PRIOR_QL_RUN --initial-root RUN --out SENSITIVITY_RUN --binary NATIVE_RUN/sensitivity-build/nb-moments-work-sensitivity
```

`run.py` checkpoints only terminal case/arm receipts and verifies existing
output hashes before reuse. `--case` permits measurement of the declared first
case; full-family coverage still requires all 58 arms. `check_real.py
--available` checks currently terminal arms and explicitly reports incomplete
coverage. This qualification's repaired rerun lives under `final/` and preserves
the initial run; `compare_requalification.py` checks every output before copying
any unchanged arithmetic receipt. `summarize.py` and `archive.py` consume that
retained layout.

The final native executable SHA256 is
`935909b17d818c1c10727c4a169e058393b3552b57f0ce9dc68a62c4df64195c`.
The sensitivity executable SHA256 is
`d91b824de8cfa99ad3f4a30361c8b845a5eea0153e45a6d134c67db41989fc8f`.
Both compile the same final production owner on Swift 6.3.3. Computation is CPU
FP64, with no Metal claim.

The [manifest](evidence/2026-09-10/manifest.json) binds complete case receipts,
final per-gene checks, grid results, independent checks, the underflow failure
and repair, work sensitivity, source hashes and native logs. Complete large
native arrays and prior QL fitted matrices remain externally retained at exact
listed paths and hashes; they are not all included in Git. The
[summary](evidence/2026-09-10/summary.json) keeps successful arithmetic and
unavailable work distinct.
