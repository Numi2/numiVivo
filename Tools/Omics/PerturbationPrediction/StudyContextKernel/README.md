# Native multigene context kernel: one study improves, transfer remains unqualified

A native Swift model uses the full untreated expression profile to predict a
shared weighting of training donors' response vectors. It improves HIRISA over
both baselines, unlike the preceding gene-wise model, but loses to the
study-balanced mean response on Kang and GSE181897. It is **not promoted**.

| Held-out study | Donors | Candidate RMSE | Mean-response RMSE | No-change RMSE | Gain over mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| HIRISA | 5 | 0.265097 | 0.336026 | 0.342458 | 21.11% |
| Kang | 8 | 1.138815 | 1.135852 | 1.224197 | −0.26% |
| GSE181897 | 62 | 1.073059 | 1.069742 | 1.121616 | −0.31% |

All five HIRISA donors improve over both baselines, with 22.59% gain over no
change. Six of eight Kang donors and 55 of 62 GSE181897 donors lose to the
mean-response baseline; one GSE181897 donor also loses to no change. Only HIRISA
meets the unchanged 5% gain requirement over both baselines. No pooled or
winner-by-study model is substituted for the full result.

## Method and exclusions

Hold each entire study out, using the exact prior 11,800 shared features and all
75 donor pairs. Inputs contain training controls and treated observations but
only query controls. No query treated value enters native fitting, selection
or prediction. Each training study has equal total weight, split equally among
its donors. Inner leave-one-training-study-out selection minimizes equal-study
mean treated log1p CPM MSE. Penalties are .01, .1, 1, 10, 100 and infinity;
exact ties prefer stronger regularization. All three folds select one.

For training control profiles X, center each gene by the weighted training
mean. Let Z_i = sqrt(w_i) (X_i - mean), K = Z Z^T / G. Solve
(K + lambda I) a = Z (q - mean) / G. With e_i = sqrt(w_i) a_i, use response
weights beta_i = w_i + e_i - w_i sum(e). Predict max(0, q + beta^T (T-X)).
Infinity returns the balanced mean response. The same coefficients apply to
all output genes, allowing control genes to inform other genes' responses.
Weights sum to one but can be negative: this is linear extrapolation, not a
convex donor mixture. There is no gene variance scaling or feature selection.

The donor-space kernel avoids a genes-by-genes matrix. The standalone driver
bounds training and query donors to 128 and each donor-by-gene payload to two
million entries; it is not a cell-level out-of-core solver or product CLI.
This is a log-expression point model, not a negative-binomial count posterior.
Dose, time, cell composition and other protocol differences are not modeled
explicitly; control expression alone does not identify all relevant context.

## Verification and reproduction

All 885,000 native predictions, training-study losses and response weights agree
with a NumPy eigensystem implementation using independent matrix products and
solves. Native fitting uses scalar kernel accumulation and Cholesky. All 225
donor/baseline RMSE values are checked with scalar compensated summation.
Four numerical/invariance checks cover cross-gene prediction, gene permutation,
duplicated feature-space normalization and query-batch independence. Seven
rejection checks cover leaked query outcomes, overlap, invalid axes and values.

Compile ContextKernel.swift with Swift 6 optimized mode. The evidence archive
contains source, preparation, tests, all native inputs/outputs, source cohort
arrays, prediction freeze, verification and scoring scripts. In a fresh directory,
restore the archive, retain its input identities, and adjust machine paths.
Run test.py, run.py, verify.py, score.py and verify_scores.py with NumPy available;
run.py requires a new native directory and test.py a new tests directory.
The source script and binary hashes were frozen before native prediction.
Use verify_archive.py to check the retained evidence without refitting.

The original public H5AD/count-source audit is inherited from the prior
[three-study experiment](../../CountObservation/Joint/Adaptive/Full/DonorExclusion/Prediction/ResponseShrinkage/StudyHeldOut/README.md);
source cohort and metadata hashes are checked again, not regenerated from H5AD
in this experiment. These are previously inspected development studies.
Participant overlap beyond study identifiers remains unverified. No fresh
external, uncertainty, clinical or general biological-outcome qualification follows.

Local storage failures occurred before prediction execution. Native runs moved
to the Mac mini with the same binary and frozen inputs. No failed biological
result was dropped or prediction run restarted.

The [exposure metadata admission review](Exposure/README.md) identifies a primary
GSE181897 protocol reference and distinguishes it from donor-level exposure mapping.
The current model inputs contain no explicit dose or duration predictors.
