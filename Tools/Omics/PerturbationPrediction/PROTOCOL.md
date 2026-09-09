# Held-out donor response protocol

Predeclared 2026-09-09 before fitting these benchmarks.

Use the complete independently verified pseudobulk count axes from Kang B cells
(8 paired human donors, ctrl/stim) and Hagai fibroblasts (3 paired mouse donors,
unstimulated/LPS6). Run every leave-one-donor-out fold separately within each study.
The held-out donor's control library is supplied at prediction time. Its treated
library is sealed until scoring. No cross-species training or disease inference.

Transform each library independently to natural-log(1 + counts per million),
using all source genes in the denominator. Do not reuse DE size factors or a
feature panel selected using held-out treated counts. Donors receive equal weight;
cells are not treated as independent biological replicates.

Predeclare four response baselines: zero change, training-donor mean change,
training-donor median change, and linear-kernel ridge prediction of response from
control context. For ridge, select genes with >=10 total training counts and
nonzero counts in at least two training donor pairs; standardize control genes
using training population means/SDs (zero SD becomes one), divide by sqrt(feature
count), and use alpha=1 with a training response-mean intercept. No tuning on
held-out donors. Fit all response genes; the selected panel is only for context.
Verify scikit-learn KernelRidge against an independent NumPy linear solve.

Predict treated logCPM by adding the estimated change to the held-out control,
then clip below zero. Score the resulting applied change and preserve both the
unclipped estimate and final prediction. Evaluate full source genes, the
training-expressed context panel, the top 200 absolute training-mean response
genes (ID tie break), and pre-existing study marker panels. Report response RMSE,
MAE, Pearson correlation when defined, sign agreement on nonzero observed effects,
and explained response sum of squares relative to the no-change baseline.
Also report expression-space metrics to expose inflated correlations from
baseline expression. Preserve per-gene predictions, training membership, panels,
transforms, model parameters and all fold metrics. No successful fold selection.

Negative controls: remove the held-out treated sample before model construction;
mutating that sample must not change training/model/prediction artifacts. Predictions
must reproduce from a frozen model and the supplied control alone. Repeated runs
must match arrays and metrics exactly. Retain failures rather than tune them away.

This establishes deterministic external baseline evidence and a future native
contract. It does not qualify an implemented native perturbation predictor,
single-cell response distributions, causal identification, unseen perturbation
identity, or unseen tissue/cell type. There is no Bayesian interval or mechanistic
claim. Hagai has only two training donors per fold; biological generalization
claims and uncertainty calibration are consequently very limited.
