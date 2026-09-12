# Native nonlinear context kernel: transfer gate remains failed

This candidate replaces linear control-profile similarity with a centered RBF
kernel. Squared distance is averaged across all 11,600 fixed features; bandwidth
is the median squared pair distance within each training fold, recomputed for
inner selection. A zero median uses bandwidth 1. Study weights remain equal in
total, with a weighted kernel center and unpenalized response intercept. The native
solver uses donor-space Cholesky and keeps the original nonnegative output clipping.

The protocol was written before native fitting. All 75 donors remain: GSE181897
62, HIRISA 5, Kang 8. Each outer fold excludes the entire query study's treated
outcomes. Training-study holdouts select penalty .01/.1/1/10/100/infinity, with
exact ties preferring stronger regularization. No bandwidth grid, response-scale
tuning, new feature selection or Parse query/target evaluation is performed.

| Held-out study | Penalty | Gain over training mean | Gain over no change | Fixed gate |
| --- | ---: | ---: | ---: | --- |
| GSE181897 | infinity (training mean) | 0.000% | 4.633% | FAIL |
| HIRISA | 0.1 | 14.720% | 16.492% | PASS |
| Kang | 1 | -0.070% | 7.172% | FAIL |

The unchanged gate requires at least 5% lower equal-donor mean RMSE against both
baselines in every study. **Not promoted.** Nonlinear control similarity does not
resolve this development transfer failure. HIRISA's gain is also smaller than the
prior linear candidate's gain. This experiment provides no reason to replace the
product default or claim reliable context prediction.

The standalone native Swift source compiled and completed all three fits. Prediction
hashes were frozen before scoring. Independent NumPy weighted-kernel eigensystem
calculations reproduce every inner loss, selected penalty, response weight, baseline
and all 870,000 predicted values. Maximum prediction difference is 8.04e-14. Every
225 donor/model RMSE comparison also passes scalar compensated-sum verification.
These checks qualify numerical computation, not biological accuracy or a full
product build. The three studies are reused development cohorts; Parse was already
inspected during prior work and is not an untouched future validation set.

`ContextKernel.swift`, `run.py`, `verify_score.py`, protocol, execution/freeze receipts,
logs and every donor score are retained here. `manifest.json` binds these artifacts,
the native executable, inputs and outputs. Large runtime files remain at
`/Users/n/numivivo-context-rbf-20260912`. Input construction verifies the frozen
source-cohort dependencies and omits each query's treated profile from its native
input. The source arrays remain dependencies under the prior context evaluation
workspace, not a new raw-count data preparation. To reproduce, use a fresh directory,
compile with `swiftc -O -parse-as-library ContextKernel.swift -o context-kernel`,
then run `run.py` and `verify_score.py` with NumPy and one BLAS thread. Preserve the
specified source dependencies and adapt only the output workspace paths.
