# Replogle UPR: fixed unseen-target prediction results

**The fixed native GO predictor meets the predeclared primary criterion in all
five original gemgroups.** Across all 33,694 RNA genes, equal-target mean response
RMSE improves by **1.41–2.12%** over the training-mean response. All 150 held-guide
folds completed, with no failed or omitted predictions. This is a modest measured
gain in a separately collected K562 UPR experiment; it does not establish reliable
prediction for every target or new biological contexts.

The subsequent [complete target-consistency diagnostic](TargetConsistency/README.md)
finds that SCYL1, SRP68 and SRP72 lose to training mean in all five technical
groups; SRP72 also loses to no change in every group. This does not change the
original aggregate result or establish independent biological replication.

## Complete primary results

The [protocol](PROTOCOL.md) preceded expression-matrix acquisition. It fixes the
native direct-GO Jaccard kernel, lambda 1, unpenalized intercept, full-source-gene
log1p CPM, clipping at zero, and no CPM reclosure. Each fit receives its own
gemgroup's pooled control and the other 29 target conditions. The held target's
expression is excluded from the fitter and query. All 30 targets have usable
descriptors, so the all-training and supported-training means coincide.

Lower all-gene response RMSE is better; values are natural-log(1+CPM) units and
are averaged equally over all 30 held targets within each gemgroup.

| Original technical group | No change | Training mean (both) | Fixed GO | Fixed shuffle | GO gain over mean | GO worse than no change | GO worse than mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| gemgroup-1 | 0.182897 | 0.171417 | **0.167786** | 0.171998 | 2.12% | 10/30 | 6/30 |
| gemgroup-2 | 0.211268 | 0.197815 | **0.194765** | 0.199233 | 1.54% | 8/30 | 7/30 |
| gemgroup-3 | 0.210073 | 0.198920 | **0.195679** | 0.199727 | 1.63% | 7/30 | 7/30 |
| gemgroup-4 | 0.189801 | 0.179391 | **0.175925** | 0.180121 | 1.93% | 7/30 | 9/30 |
| gemgroup-5 | 0.193951 | 0.180638 | **0.178085** | 0.181375 | 1.41% | 2/30 | 8/30 |

Every gemgroup beats all three required comparators. Nevertheless, GO is worse
than no change in **34/150** individual folds and worse than the mean in
**37/150**. These are counts of technical-group/target outcomes, not 150
independent biological replicates. Every per-target result, failure list, RMSE,
MAE, correlation and sign metric is retained; no favorable group or target was
selected after scoring.

The secondary panel uses the top 1,000 genes by mean CPM across the control and
training targets only, with feature-ID tie breaking. Its GO/mean response RMSE
is 0.115894/0.130663, 0.133169/0.149635, 0.124241/0.142007,
0.108216/0.124455 and 0.112073/0.127221 in gemgroups 1–5, respectively.
All five also beat the shuffle on this secondary panel. The primary full-gene
criterion determines the result; panel-specific gains cannot replace it.

## Inputs, identity and separation

The [complete original cohort](README.md) and [native preparation](TRAINING_PREPARATION.md)
retain all 40,997 source cells, all RNA and guide assay features separately, and
the exact 32,829-cell confident selection. Five native training inputs have 31
rows each: one pooled control and all 30 guides. Both controls pool only within
their original gemgroup. The 150 count exclusions were frozen before descriptor
resolution, and all held-count mutation checks passed.

[Nominal gene identity evidence](IDENTITIES.md) combines original guide labels,
FBA's independent reproduction of the guide/gene sequence table, exact source
Ensembl identities, reference sequence core matches, and one explicit HGNC
previous-symbol record. It is not a newly retrieved original supplement or
validation of full guide sequences, off-target specificity or knockdown efficacy.
All 30 exact Ensembl identities have usable direct GO terms from MyGene build
`20260906`. No captured direct annotation cites either benchmark paper, but
current knowledge remains potentially informed by these studies.

All 150 descriptor plans, query inputs and selected-training hashes were frozen
before the first native fit. All 150 native predictions were then frozen before
the independent checker/scorer read held responses. No hyperparameter tuning,
feature-denominator change or outcome-based identity selection occurred. The
fixed one-position response shuffle uses the same supported training targets.

## Numerical and execution evidence

The actual scoped single-cell product ran on the physical M4 Pro using the same
executable qualified for ingestion and training preparation, SHA-256
`36df30cd11afffb9bca0f88369af2dca3260759ddd9945c2c02d36073b1574db`.
There were **601 successful native commands**: 150 fits, 150 model replays,
150 predictions, 150 prediction replays and one identical repeated prediction.
No new Swift implementation, full-application qualification or GPU prediction
is claimed in this experiment.

The independent checker reconstructs all five training matrices from the
original-count SciPy reference and compares all **750 prediction vectors**.
It uses both direct NumPy weights and an independently centered scikit-learn
kernel-ridge fit. All comparisons pass the retained tolerances. Maximum vector
error against the direct reference is **1.78e-15** and maximum weight error is
**3.33e-16**. Implied CPM sums agree within tolerance and remain unrenormalized.

There are 246 vector-level differences in clipped-feature counts, bounded by
the number of reference features within 1e-11 of zero; these rounding-sensitive
diagnostics are retained separately. They are not presented as exact cross-library
clip-count agreement. Response metrics subtract the same verified native control
vector from prediction and truth, making the no-change response exactly zero.

Two complete checker/scorer executions produce identical scientific artifacts.
The independent reference environment uses NumPy 2.5.3, scikit-learn 1.9.0 and
SciPy 1.18.1. No dense cell-by-gene matrix is created; fitting and reference
calculations use sparse counts and small condition-by-gene matrices.

## Scope of the result

This is independent collection of data for testing the unchanged algorithm,
not transfer of Norman-trained coefficients. Replogle and earlier studies share
investigators and K562 experimental systems; UPR targets were selected using
earlier Adamson results. Five gemgroups are technical/capture conditions, not
independent donors or laboratories. Results retain original IDs without assuming
experiment-number-to-platform mappings.

The result supports a small average improvement in predicting the measured RNA
response of a withheld nominal target within these contexts. It does not
establish prospective target selection, calibrated uncertainty, individual-cell
distributions, causal mechanisms, RNA-to-phenotype transfer, new tissues or
clinical outcomes. Adamson's separate experimental-role gate remains unresolved.

## Reproduction and retained evidence

`run_target_kernel.py` freezes inputs and orchestrates the native owners;
`check_target_kernel.py` verifies the frozen outputs and scores every fold.
The [archive](evidence/2026-09-11-target-kernel/manifest.json) retains exact source
captures, plans, all native command logs and metadata, numerical checks, complete
metrics and the scoring-repetition receipt. Bundled source files preserve their
original bytes and hashes. Large prediction vectors, training files and the
first complete model remain external under recorded stored and decoded hashes
in `macmini:/Users/n/numivivo-replogle2020-20260911/target-kernel`, with a verified
local copy. Other models reconstruct from frozen inputs and their retained hashes.

The native input-freeze SHA-256 is
`951526564ad9ed0b4d97ce275f4f12753b0da661124b4bf27bc628cb9b46ffb3`;
the prediction-freeze SHA-256 is
`735e659449a6870bafc719e2613df8a5176ed814a16c4503f310e3e64619cb70`.
Retain the earlier source and preparation archives under their original identities.
Destinations for new runs must not already exist.

```sh
python3 Tools/Omics/PerturbationPrediction/Replogle2020/run_target_kernel.py \
  /path/to/qualified-study --descriptors /path/to/identity-author-analysis/descriptors \
  --out /path/to/new-native-predictions
python3 Tools/Omics/PerturbationPrediction/Replogle2020/check_target_kernel.py \
  /path/to/qualified-study --native /path/to/new-native-predictions \
  --descriptors /path/to/identity-author-analysis/descriptors --out /path/to/new-checks
python3 Tools/Omics/PerturbationPrediction/Replogle2020/verify_target_archive.py \
  Tools/Omics/PerturbationPrediction/Replogle2020/evidence/2026-09-11-target-kernel
```
