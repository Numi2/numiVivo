# Shrinkage transfer failure: observed range is insufficient as an admission gate

This post-hoc diagnostic retains every frozen prediction from the native
26-fold transfer experiment. It does not change the model, select a new penalty,
exclude genes, calibrate an OOD score, or promote a passing subset.

| Transfer direction | Query genes outside observed training range | RMS slope correction, log1p CPM | Candidate / mean RMSE inside range | Outside range |
| --- | ---: | ---: | ---: | ---: |
| HIRISA → Kang, 8 donors | 90.37–92.09% | 0.118–0.159 | 0.998–1.003 | 0.993–1.001 |
| Kang → HIRISA, 5 donors | 71.42–73.23% | 1.653–1.712 | 1.820–1.876 | 3.467–3.617 |

The direction with greater observed-range coverage fails much more severely.
Kang → HIRISA also fails inside the training range for every donor. Range
membership alone therefore does not establish a useful accuracy gate here.
These are per-gene marginal ranges from few donors, not a multivariate
applicability model or a biological threshold.

The candidate adds `slope * (queryControl - trainingControlMean)` to the mean
response before clipping. Its correction is much larger in Kang → HIRISA.
For the post-clipping prediction shift delta and mean-baseline residual r, the
paired squared-error difference is exactly:

    candidateError² - meanBaselineError² = delta² + 2 * delta * r

All per-gene identities pass. In each of the five failing transfer donors, the
cross term is negative: the shift aligns with reducing baseline error in this
aggregate algebraic sense. Its squared magnitude nevertheless overwhelms that
benefit. The squared-shift term is 111.4–113.6% of the net error increase, offset
by a negative cross term. This is an exact decomposition of these predictions;
it does not prove a biological mechanism or justify rescaling from inspected
query outcomes.

## Consequence for the next model experiment

Do not install a range filter as a fix, promote the inside-range subset, or
select different penalties from the failed transfer outcomes. The next
selection design needs held-out study/context structure in its training data,
not only within-study donor validation. A third condition-qualified study can
support such development; its assay, population, dose/time and feature contract
must be explicit. Already inspected cohorts remain development evidence, and a
subsequent untouched outcome cohort is still needed for a fresh claim. Simply
expanding the penalty grid on these same donor holdouts would not address that
design gap.

## Evidence and reproduction

`diagnose.py` reads only the original frozen training/query input and fitted
output files; it binds their hashes, records every gene's observed range and
prediction correction, and recomputes the unchanged predictions.
`score_strata.py` separately reads the frozen outcome counts, preserving all
library denominators and the exact 11,884-gene panel, and decomposes paired
errors in all 26 folds. `verify.py` independently reconstructs every support
record using scalar operations: all 308,984 records pass.

Run the three scripts in that order with NumPy available. They reference the
retained original `numivivo-shrinkage-transfer-20260912` and
`numivivo-cross-study-ifnb-20260911/inputs` study directories; adapt the explicit
paths if restoring elsewhere. The archive includes scripts, reports, logs and
hashes. Per-gene support files remain in the local study with hashes recorded in
`external-artifacts.json`; they can be regenerated from the prior published
native input/output archive. No large original atlas or count file is duplicated.
