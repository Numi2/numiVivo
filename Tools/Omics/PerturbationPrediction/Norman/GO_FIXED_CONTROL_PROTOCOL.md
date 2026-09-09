# Post-result matched fixed-lambda control

Declared 2026-09-10 after the primary GO_TRANSFER_PROTOCOL.md results: fixed
lambda=1 has lower mean RMSE than both means, but nested lambda selection fails
the primary all-gene comparison. The original protocol included a shuffled
nested model, which does not match the fixed model's regularization choice.
This audit addresses that comparator gap; it is not a new untouched benchmark.

Keep the frozen original fixed-model query weights, all 105 folds, original
source counts, source annotations, support mask and evaluation panels. For each
supported held target, rotate sorted supported training responses by one position,
then apply exactly the fixed lambda=1 weights. Independently refit the centered
kernel with scikit-learn at lambda=1 to verify the shuffled prediction. Neither
the weights nor the held outcomes can change any selection or model parameter.

Freeze this complete control archive before scoring. Compare its full all-gene
and training-top-1,000 per-target metrics against the original fixed, nested,
mean and no-change results. Repeat prediction and scoring exactly. Retain every
failure and coverage exclusion. A gain supports this fixed candidate only on
this inspected development cohort, not native, unseen-study, donor/tissue,
causal or single-cell distribution qualification.
