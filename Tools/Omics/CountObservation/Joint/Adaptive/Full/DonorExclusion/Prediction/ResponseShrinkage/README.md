# Training-selected response shrinkage

This development experiment follows inspected joint-count and response-transport failures. It is not a fresh external test. It models log1p pseudobulk RNA response from the untreated pseudobulk profile, retaining the original 13 donor exclusions.

For each gene, let x = log1p(control CPM), y = log1p(treated CPM), and d = y - x. Fit an unpenalized response intercept and a ridge-regularized control slope:

    b = sum((x - mean(x)) * (d - mean(d))) /
        (sum((x - mean(x))^2) + lambda * (trainingDonors - 1))
    prediction = max(0, queryX + mean(d) + b * (queryX - mean(x)))

Zero denominator gives zero slope. Infinite penalty is the training-mean response limit. Penalties are 0, 0.01, 0.1, 1, 10, 100 and infinity. Select one penalty per outer fold by mean squared error across all available genes and inner held-out donors; ties prefer stronger shrinkage. Then refit on all outer training donors. Query treated outcomes are read only by the separate scorer, after predictions are frozen.

The input shards were exported using outer-training count/dispersion availability. Inner splits do not refit those availability filters: this is not fully nested dispersion qualification. The fitted values use observed pseudobulk CPM directly; this candidate is neither a new NB observation likelihood nor the original Bayesian posterior, and supplies no calibrated uncertainty. Outer scoring uses exactly the original conditionalPrediction gene population. Training feature availability can therefore exceed the scored population.

`freeze.json` binds selection code before execution. `selection.json` records every training input hash, all candidate inner losses, donor identities and fitted model hashes. `predict.py` checks exact training donors and outer exclusion. `verify_models.py` compares the first available gene in every shard with an independent augmented least-squares solution; it does not independently reproduce every CV loss. `verify_predictions.py` recomputes every prediction using an affine expression and independently verifies query input hashes and donor exclusion. `score.py` compares the frozen original outcomes and baselines without replacing their scores.

Reproduction uses the retained `numivivo-joint-fold-fits-20260912` and `numivivo-joint-fold-prediction-20260912` study directories. Adapt the explicit paths in scripts if restoring elsewhere. Run select.py, verify_models.py, predict.py, verify_predictions.py, then score.py in that order using the scientific Python environment. The predictions directory must not already exist.

## Complete development result

| Dataset | Candidate RMSE | Original joint RMSE | Training-mean RMSE | Gain over mean |
| --- | ---: | ---: | ---: | ---: |
| Kang, 8 donors / 65,895 genes across folds | 1.049574 | 0.999211 | 1.311326 | 19.96% |
| HIRISA, 5 donors / 65,447 genes across folds | 0.088835 | 0.158447 | 0.090109 | 1.41% |

Every donor improves against training mean. HIRISA improves over the previous
joint and transport candidates but still misses the existing 5% improvement
gate. Kang passes that aggregate gate but is less accurate than the original
joint predictor. Do not promote this as a qualified replacement or combine
model winners by dataset after scoring. Training-only regularization helped;
fresh independent evaluation and calibrated uncertainty remain necessary.

All Kang training folds selected lambda 0.1. HIRISA selected lambda 1 in four
folds and 0.1 in fold 03. Selection read no outer query-treated outcomes, but
these development cohorts had already been inspected before this experiment.
All 131,342 predictions passed independent numerical recomputation. Independent
augmented least-squares checks matched 3,383 genes (one available gene per input
shard), with maximum coefficient difference below 1e-14. Full CV losses were not
independently reimplemented. The original count-model failure statuses and
original evaluation files were retained unchanged.

`evidence.tar.gz` retains scripts, the pre-execution selection freeze, all 13
fitted models and predictions, input/model hash bindings, verification, scores
and logs. Raw training/query shards remain in their original retained studies.
