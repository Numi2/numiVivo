# Held-out combination prediction protocol

Predeclared 2026-09-09 before scoring the paired-condition outcomes.

Use the full, qualified Norman filtered release and its already published split:
one control, 105 single-target training conditions, all 131 paired-target query
conditions. Both constituent singles are observed for every query. Keep every
source feature and condition. This is pooled condition-mean prediction, not
donor generalization, unseen-target prediction, or single-cell distribution
prediction. The 105 unseen-target folds remain a separate task.

Normalize each condition independently to CPM with the complete 33,694-feature
library denominator. Let c be control CPM, s[t] single-target CPM, and d[t] =
log1p(s[t]) - log1p(c). Give each single condition equal weight. Freeze these
training-only quantities before making any paired outcome available to scoring.
The predictor's inputs are a training-only archive and query target identities;
it must not accept a paired outcome archive.

For query targets a and b, predeclare:

1. `noChange`: log1p(c).
2. `meanSingleResponse`: log1p(c) + mean over all 105 d[t].
3. `additiveLogResponse`: log1p(c) + d[a] + d[b].
4. `meanConstituentLogResponse`: log1p(c) + (d[a] + d[b])/2.
5. `additiveCPMResponse`: log1p(max(c + (s[a]-c) + (s[b]-c), 0)).
6. `shuffledAdditiveLogResponse`: additive log prediction after rotating the
   sorted training-target names by one position. This fixed identity-negative
   control is not a permutation test or a fitted model.

Clip negative predicted log expression to zero. Preserve model quantities and
unclipped responses through deterministic reconstruction from the frozen model;
save final per-gene predictions. Do not reclose predicted CPM compositions or
interpret them as count predictions. Report implied CPM sums and clipping counts.
No hyperparameters or scales may be tuned on the paired outcomes.

Predeclare gene panels, all selected without paired outcomes:

* All source genes.
* Training-expressed: >=10 total counts and expression in >=2 control/single
  libraries.
* Top 1,000 training-mean CPM genes, with Ensembl-ID tie breaking.
* Query-specific top 200 absolute additive single-target log responses, with
  Ensembl-ID tie breaking. The query identities select already observed singles;
  paired counts do not enter selection.

For each query, baseline and panel report response RMSE, MAE, Pearson when
defined, sign agreement on nonzero observed responses, and explained response
sum of squares relative to zero change. Report expression-space correlations
separately. Aggregate metrics weight query conditions equally; preserve all
per-query results and worse-than-no-change cases. Report comparisons against the
shuffled-target baseline without inventing biological-replicate confidence.

Checks: paired outcomes are physically separated from fitting/prediction inputs;
training-only model reload reproduces predictions; reversing target order gives
identical predictions; unknown/repeated targets and nonfinite/negative counts
are rejected. Mutating every sealed paired outcome leaves training, model and
prediction arrays unchanged. A second run must reproduce every array and metric.
Compare additive-log predictions with an independently fitted scikit-learn
linear least-squares model on the single-target identity design. This verifies
the arithmetic, not independent biological truth.

The original [Norman study](https://doi.org/10.1126/science.aax4438) motivates the
combinatorial task. [Ahlmann-Eltze et al.](https://www.nature.com/articles/s41592-025-02772-6)
motivate an additive reference; their feature processing, conditions and splits
differ, so this is not a numerical paper reproduction or a comparison against
their deep-learning runs. [Systema](https://www.nature.com/articles/s41587-025-02777-8)
motivates inspecting systematic response separately from target-specific benefit.
No GEARS, CPA, foundation-model, Bayesian uncertainty or mechanistic qualification
is claimed by this benchmark. An additive model cannot establish genetic
interaction prediction merely by predicting an unseen pair.
