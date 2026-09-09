# Unseen-target co-response baseline: insufficient predictive transfer

The fixed co-response descriptor model does **not** improve on the mean-single
baseline on the full Norman held-target benchmark. This negative result changes
the next implementation decision: do not promote it as the native unseen-target
solution. Better descriptors or models need new predeclared evaluations.

The [protocol](UNSEEN_TARGETS_PROTOCOL.md) was written before scoring. Each of 105
folds removes every condition involving its held target. Only control and the
other 104 singles are used; all paired counts are unused. The query is the held
single, across all 33,694 genes. Of 105 targets, 102 have unique exact source-symbol
matches. C19orf26, C3orf72 and KIAA1804 remain unsupported for this descriptor;
the no-change and mean-single baselines still score those three.

Descriptors contain a target gene's response to other available perturbations.
This is target-specific information measured without applying the held
intervention. It is not a perturbation identity lookup. Training descriptors do
contain self-intervention entries, unlike the query: this transfer-distribution
mismatch is a possible limitation, not an established explanation of failure.
No descriptor or hyperparameter was revised after observing results.

## Results

Mean per-target log1p-CPM response RMSE on the same 102 supported targets:

| Method | All genes | Training top 1,000 | Worse than no-change, all genes |
| --- | ---: | ---: | ---: |
| No change | 0.133473 | 0.122522 | 0 |
| Mean training-single response | 0.127156 | 0.105677 | 23 |
| Co-response ridge | 0.132229 | 0.109411 | 44 |
| Shuffled co-response ridge | 0.133835 | 0.112642 | 45 |

Ridge beats the mean baseline on only 18/102 all-gene targets (43/102 on the
training top-1,000 panel). It beats its matched shuffled control on 65/102 and
67/102 respectively. This fixed rotation is not a permutation significance test.
Per-target scores preserve all failures, correlations, sign agreement, and MAE.
Predictions are clipped at zero log-expression without CPM reclosure; diagnostics
retain clipped feature counts and implied CPM sums. Unsupported numeric rows are
explicit NaNs and never silently scored as predictions.

## Verification and scope

Both complete runs reproduced every stored array and receipt exactly. All 105
fold input matrices remain identical after changing every excluded condition's
counts; an AHR full refit also reproduces predictions exactly after that mutation.
Condition reordering leaves selected counts identical. Five malformed-count
cases reject. The independent [scikit-learn SVD Ridge implementation](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html)
agrees within 1.49e-14 absolute response error on all supported folds, for both
ordinary and shuffled outputs. The script uses the native-qualified condition
reference and creates no dense cell-by-gene matrix.

Prediction and scoring are separate commands. The prediction orchestrator reads
the full pinned reference to select folds and check exclusion; the fitter only
receives selected control/training arrays. This is function-level isolation,
not process-level sealed storage. Each prediction archive is finalized and hashed
before the separate scoring command reads held outcomes. The two score runs are
also exactly identical. The prediction bundle freezes outputs and fold metadata;
it is not a reusable native model bundle.

The repeated prediction run took 15.38 seconds, peak RSS 1,143,783,424 bytes on
the local Apple M4, using one BLAS thread. This includes reference load, all fits,
independent checks, isolation checks and archive compression; it is not a native
speedup or million-cell qualification. Runtime/package details and logs are in
`evidence/2026-09-09-unseen-targets`.

These are pooled K562 CRISPRa condition means, without independent biological
replication. Pair responses involving unseen targets, new cell types/donors,
Bayesian uncertainty, causal mechanisms and reaction kinetics remain unqualified.
No Swift runtime behavior changed in this experiment.

## Reproduce

Use the benchmark Python environment and the existing byte-pinned `reference.npz`
from `prepare.py`. Invoke from this directory:

```sh
python unseen_targets.py predict --reference /path/reference.npz --out first
python unseen_targets.py score --reference /path/reference.npz --predictions first --out scores
python unseen_targets.py predict --reference /path/reference.npz --out repeat
python check_unseen_targets.py --reference /path/reference.npz --first first --repeat repeat --output checks.json
python unseen_targets.py score --reference /path/reference.npz --predictions repeat --out repeat-scores
```

All output destinations must be new. The evidence retains first predictions,
fold membership, checks, both run logs, score tables and provenance. The pinned
input reference remains in the existing Norman qualification; it is not duplicated.
