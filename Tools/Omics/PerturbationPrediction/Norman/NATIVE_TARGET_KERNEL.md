# Native target-kernel prediction

`VivoTargetKernel` fits a reusable model from control and single-target condition
counts plus independently supplied target annotations. Queries supply an unseen
target identity and its annotation terms. The native owner uses term-set Jaccard
similarity and kernel ridge with an unpenalized intercept to predict log1p-CPM
expression in the exact declared training context.

This implements the fixed-regularization candidate from [GO_TRANSFER.md](GO_TRANSFER.md).
Norman remains a reused development benchmark. The earlier co-response,
control-correlation and nested-GO failures remain part of the evidence. The owner
accepts other declared term namespaces; that interface does not establish their
predictive value.

## Input and model contract

The input training JSON is the existing [composition training artifact](NATIVE_COMPOSITION.md),
with control first and then selected target-condition rows. Use
`singlecell-composition-prepare` to select counts from a freshly reconstructed
native pseudobulk bundle. The fitter receives only this selected training file.
Target-specific held-out responses, paired responses and scoring outcomes are
not fields of the fit or query plans.

```sh
numivivo singlecell-composition-prepare source-bundle --plan selection.json --output training.json
numivivo singlecell-target-kernel-fit training.json --plan descriptors.json --output model
numivivo singlecell-target-kernel-verify model
numivivo singlecell-target-kernel-predict model --plan queries.json --output prediction
numivivo singlecell-target-kernel-prediction-verify prediction
```

Each descriptor has `targetID`, optional `featureID`, `status` and sorted unique
`terms`. `available` requires a feature identity and nonempty terms; `noData`
requires an identity and empty terms; `unresolvedIdentity` requires no identity
and no terms. Training descriptors must cover exactly the selected target IDs.
Feature identities cannot be shared by two training targets.

A fit plan contains schema version 1, the exact composition `context`,
`descriptorNamespace`, annotation `source` fingerprint, `provenance`, `targets`,
`regularization` (default 1) and `maximumWork` (default 200 million). A query plan
contains schema version 1, matching context/namespace, an annotation source
fingerprint, unique query IDs with descriptors, and its work limit. A source
fingerprint binds supplied annotation provenance; it does not authenticate an
annotation service or require that future queries use the same capture file.
Unknown input settings reject.

Prediction rejects any target ID already present in training and any query
feature ID matching a known training feature, including unsupported training
descriptors. It cannot resolve an undisclosed alias or authenticate an independently
supplied identity. Missing annotations and unresolved identities remain explicit.
An otherwise available query with no shared terms returns `noSharedTerms` and no
target-specific weights or expression. Generic baselines remain available.

The model freezes original feature IDs, complete control/target CPM denominators,
sorted target IDs, descriptors, supported training rows and the inverse of
`K + lambda I`. A native Cholesky solve checks the inverse residual. With
`C = inverse(K + lambda I)`, `u = C 1`, and query similarity row `k`, the weights are
`k C + (1 - k C 1) u^T / (1^T u)`. The intercept preserves a constant response.
Negative predicted log-expression is clipped to zero; clipping counts and implied
CPM sums are retained without renormalization.

Every report includes no-change, all-training-target mean and supported-target
mean baselines. Supported queries also retain a fixed negative control that
rotates supported training responses by one position relative to descriptors.
This is one matched shuffle, not a permutation significance test. Regularization
is fixed by the supplied plan; this owner does not select it using outcomes.

## Artifact integrity and resource bounds

Fit bundles contain training, canonical plan, model and receipt. Prediction
bundles copy the complete model reference and retain query, report and receipt.
Both verifiers snapshot input files privately, check fingerprints and the current
implementation identity, reconstruct the model from training and descriptors,
and compare canonical derived bytes. The existing product identity includes the
executable hash and operating-system version; copying the same binary to another
OS version does not authorize its original receipts. Prediction verification also recomputes the
report. Hashing an edited derived artifact again does not bypass reconstruction.
New destinations are required and incomplete staged outputs are removed.

There are 2–256 training targets, at least two supported training descriptors,
1–100,000 features and at most ten million condition-feature values. Existing
sparse-training limits apply, including five million nonzeros, positive libraries
and exact integer conversion through a total of at most 2^53. Descriptors have at
most 4,096 terms each and 100,000 terms per plan. Queries are limited to 64 and
bounded output/work estimates are checked before materialization. Lambda must
be finite and between 0.001 and 1,000. Training/plan/model/report file limits are
64/2/128/512 MiB. These are resident condition matrices and JSON artifacts; the
owner does not create a dense cell-by-gene matrix or qualify million-cell fitting.

## Qualification procedure

`run_native_target_kernel.py` runs all 105 Norman held-target folds with the actual
native executable. It prepares all single-condition counts once from a freshly
verified source, then selects control and the other 104 single targets before
invoking each fit. All paired conditions are excluded. Mutating the held row must
leave every selected training input unchanged. Each fold runs fit, model verify,
predict and prediction verify. The first fold also checks exact repeated
prediction and reordered-training model bytes.

`check_native_target_kernel.py` checks every selected count against the frozen
independent source, every native gene vector against the previously frozen GO
and fixed-shuffle predictions, and every supported query against an independent
centered-kernel scikit-learn solve. It separately scores sealed responses on the
all-gene and frozen training-selected top-1,000 panels. Clipping classification
differences near floating-point zero are retained explicitly.

`check_target_kernel_integrity.py` exercises the product CLI with synthetic
numerical controls and invalid inputs, including rehashed model, descriptor,
query and report edits. Swift tests cover analytical weights and matched shuffles,
unsupported descriptors, identities, resource limits and artifact reconstruction.
These controls establish software behavior; they are not biological validation.

## Results on 2026-09-10

On the physical M4 Pro Mac mini, the release product, scoped H5AD build and all
95 single-cell tests in 28 suites passed. The focused suite passed 10 tests,
including the five new target-kernel tests. CLI qualification passed 26 commands,
including 22 expected rejections. Rehashed model/descriptor edits failed model
reconstruction; rehashed query/report edits failed prediction reconstruction.
The synthetic centered-kernel comparison differed by at most 3.56e-15.

Fresh native aggregation covered all 111,445 cells and reproduced the existing
full-source report hash exactly. Every one of 105 held-target fits and predictions
passed reconstruction, for 423 commands including preparation and replay checks.
All prepared counts matched the independent reference exactly. Across 33,694
genes, all 517 available baseline/query vectors matched frozen reference outputs
within 2.45e-15. All 101 supported queries also matched independent centered-kernel
predictions within 1.78e-15; the maximum query-weight difference was 3.27e-16.
Comparison and scoring files were byte-identical when repeated.

For the same 101 supported targets, mean per-target response RMSE was:

| Method | All genes | Training-selected top 1,000 | Worse than no change, all genes |
| --- | ---: | ---: | ---: |
| No change | 0.133391 | 0.122163 | 0 |
| Mean of all training singles | 0.127227 | 0.105720 | 23 |
| Mean of supported training singles | 0.127239 | 0.105737 | 22 |
| Fixed GO kernel, lambda 1 | 0.126490 | 0.102073 | 29 |
| Matched fixed shuffle | 0.127379 | 0.104376 | 26 |

These reproduce the earlier fixed-GO development results: about 0.58% lower
all-gene RMSE than the all-single mean. The kernel remains worse than that mean
for 42 of 101 targets and worse than no change for 29. It is a modest aggregate
gain with retained failures, not reliable general unseen-target prediction.
C19orf26, C3orf72 and KIAA1804 remain unresolved; IER5L has no usable GO terms.
All three generic baselines are available even for these four targets; no
target-specific expression is fabricated for them.

There were 162 query/method clipping-count differences against the independent
matrix-product calculation. The archived `check_clipping.py` reconstructs native
ordered arithmetic and matches every native expression value and clipping count
exactly. Every sign disagreement is at a magnitude of at most 1.57e-17. Per-gene
values are retained in `clipping-detail.json`; the differences are not silently
discarded. Implied CPM-sum diagnostics differ by at most 1.41e-8 from independent
summation.

The retained first model was also restored on the other host. Its original
receipt correctly rejected macOS build 25G5028f because it was qualified on build
25G72. Refitting with the unchanged executable on the local OS produced
byte-identical model/query/report values and new platform-bound receipts. The
original rejection and receipts are retained; cross-OS receipt replay is not
claimed.

The 105-fold driver took 714.8 seconds including fresh source reconstruction,
423 CLI commands and evidence handling. Median fit/predict command times were
1.35/1.49 seconds; maximum fit/predict resident memory was 577,617,920 bytes.
These include JSON and reconstruction costs and are not a speed comparison with
the Python experiment or a GPU qualification. The first direct invocation of the
non-executable scoped build script returned permission denied; invoking the
tracked script with `bash` succeeded, and both logs remain in the evidence.

The qualified executable SHA-256 is
`8ac67c1100f58c9191650cc956bc0db82dff07359d7d56be749cbebb51ec5e24`.
Exact evidence and readable summaries are under
[`evidence/2026-09-10-native-target-kernel`](evidence/2026-09-10-native-target-kernel).

## Retention and interpretation

All per-fold prediction reports, plans, queries, receipts, model hashes and command
logs are retained. One complete fitted model is retained, together with the full
prepared single-condition training file. The other 104 model files are not
retained in full: reconstruct each by selecting its recorded training target list
with the driver's `subset` function, serializing with `write`, and running the
native fit using its saved plan. Check the exact model hash in that fold's manifest.
This avoids storing many almost identical copies of condition counts and CPM.

Annotation coverage and scientific interpretation follow the pinned GO capture:
101 supported targets, three unresolved source-symbol identities and one gene
without usable GO terms. Current annotations do not establish temporal independence
from Norman. Same-study held-target prediction does not qualify unseen donors or
tissues, cell-level distributions, calibrated uncertainty, genetic interactions,
causality, reaction kinetics, Bayesian/mechanistic coupling or Metal execution.
