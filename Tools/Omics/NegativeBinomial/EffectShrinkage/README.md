# Explicit count-likelihood effect shrinkage

An optional `negativeBinomialOptions.effectPriorStandardDeviationLog2` adds a
zero-centered normal prior on the requested log2 contrast. The caller supplies
a finite SD in `[0.01, 100]`; omission retains the existing inference path.
For example:

```json
"negativeBinomialOptions": {
  "trend": "gammaParametric",
  "effectPriorStandardDeviationLog2": 1.0
}
```

At the final fixed NB dispersion, the native optimizer maximizes the **raw
count log likelihood minus `(c' beta)^2 / (2 s^2)`**, with `s = log(2) * SD_log2`.
It jointly refits the nuisance coefficients and uses observed-information
Newton steps, an augmented QR penalty row, line search and a scaled score
convergence check. Only the requested contrast receives a prior; equivalent
full-rank design encodings preserve the estimate. The prior row does not add
biological observations, donors or residual degrees of freedom. The original
512-observation bound and positive-count support/replication gates remain.

Only genes admitted to the original Wald test are eligible. A failed or
unconverged shrinkage attempt is recorded independently and cannot remove a
Wald test or change the BH family. Original feature effects, standard errors,
intervals, p-values and adjusted p-values remain unchanged.

Reports add `negativeBinomial.effectShrinkage` with eligible, converged and
failed counts. Each eligible feature diagnostic contains `effectShrinkageFit`
and, on failure, `effectShrinkageError`. The raw fit's coefficients, effect,
prior SD and posterior SD use **natural logs**. Divide effect and SD by `log(2)`
for log2 units. `converged: false` is an unavailable estimate, even when raw
iteration diagnostics exist. The SD is the contrast projection of the inverse
observed negative log-posterior Hessian: a **conditional Laplace approximation**,
not a replacement Wald SE. It treats dispersion and the supplied prior as fixed.
No calibrated posterior interval or new significance test is produced.

This option supplies a fixed, contrast-specific normal prior. The separate
[empirical-prior option](../EmpiricalPrior/README.md) estimates its width using
weighted quantiles; neither uses apeglm's heavy-tailed shrinkage. They do not resolve the retained
null-calibration findings in [NullBenchmark](../NullBenchmark/README.md) and
[ProfileAudit](../ProfileAudit/README.md). Prior uncertainty, posterior
coverage, effect-estimation risk, power and general biological calibration
remain open.

## Reproduction

`prepare.py` freezes the measured Kang and Hagai treatment plans under both
original and active-donor policies, preserving original reports and receipts.
The qualification prior is fixed at **1 log2 unit before running the fits**.
`run_native.py` publishes and reconstructs all four native H5AD bundles, then
stores lossless report/plan/receipt archives with exact source hashes. Original
source H5ADs remain externally retained; `restore_bundle.py --archive ...
--source <hash-matching-original.h5ad> --out ...` restores a native bundle.

`check.py --root ...` independently solves the penalized score equations with
SciPy from the unpenalized coefficients and checks the observed-information
Laplace SD, count likelihood, penalized objective, means, effects and stationarity
for **every tested gene**. It compares original counts, design, inference,
trend and unshrunk diagnostics exactly. Solver status messages are retained;
qualification requires the declared numerical residual tolerances regardless
of the solver's flag. The checker records all per-gene discrepancies and fails
if any declared tolerance or monotonic-shrinkage check fails.

Swift regression tests additionally cover equivalent design encodings, offset
shifts, nuisance refitting, a weak-prior limit, deficient support, nonconvergence,
invalid priors, the 512-observation limit and active-donor propagation. These
numerical fixtures are not measured biological evidence.

## Native qualification, 2026-09-10

The final physical M4 Pro CLI SHA256 is
`2d9e524cbb1cdd552fc9f52a1c5237cc904f9faec6eb6bc30b7f49d479ee5433`.
All 17 focused Swift tests passed. Four complete native treatment runs and four
replays passed; a separately restored Kang archive also replayed successfully.

| Measured scope | Shrunk/tested genes | Active-support fits | Original BH <= 0.05 calls |
| --- | ---: | ---: | ---: |
| Kang, original donor policy | 5,400 | 0 | 850 |
| Kang, active-donor policy | 8,829 | 3,429 | 815 |
| Hagai, original donor policy | 12,426 | 0 | 4,827 |
| Hagai, active-donor policy | 12,426 | 0 | 4,827 |

All **39,081 gene/contrast fits** passed the unchanged independent numerical
tolerances, with zero native failures and zero unsuccessful SciPy solver flags.
These are four analyses of two measured studies, with overlapping genes and
observations; they are not four independent biological benchmarks. Counts, QC,
designs, dispersion trends, original fit diagnostics and every original feature
inference field remained exactly equal to their archived baselines. The BH
counts above describe those unchanged original tests, not shrinkage significance.

The maximum native/reference natural-log effect difference was `6.638e-8`;
maximum posterior SD difference `4.394e-9`; maximum relative mean difference
`2.395e-7`; maximum log-likelihood difference `7.167e-7`. All independently
recomputed scaled scores were below `1e-7`. Full per-gene checks, source hashes,
plans, receipts, logs and compressed reports are in
[evidence/2026-09-10](evidence/2026-09-10/manifest.json).

The first implementation used ambiguous `counts.map(Double.init)`, which Swift
resolved to the UInt64 **bit-pattern** initializer. Its 3,617 failed shrinkage
attempts, failed 16-test run (five issues), source snapshot, reproduction and
unqualified-report inventory are retained in the evidence. None of that first
attempt's shrinkage is qualified, including fits flagged converged. Its four
original Wald inference outputs remained exact. Explicit `Double(value)`
conversion fixed the defect; a 70-digit scalar count-likelihood regression now
covers the conversion and posterior curvature. All four measured runs and
independent checks were repeated with the corrected executable. Build logs also
retain pre-existing warnings from unrelated Metal and workflow sources.
