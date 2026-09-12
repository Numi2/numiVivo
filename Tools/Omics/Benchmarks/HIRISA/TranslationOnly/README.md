# Full-HIRISA translation-only ablation: preservation still fails

Removing donor-specific rotations reduces annotation failures to **22/146**
eligible comparisons, including **5/29 rare** comparisons. The previous rigid
projection fails 31/146 (9 rare), and native correction fails 37/146 (10 rare).
Translation resolves 24 native failures but introduces nine others. Mean eligible
recall is 0.649932, below the rigid projection's 0.650730; fewer threshold failures
do not imply better mean recall. No method is promoted.

The experiment retains all **1,612,594 cells, 20 PCs and 114 annotation folds**.
For each donor, add the difference between its frozen native-corrected mean
and original PCA mean to every original PCA row. No rotation, scaling, labels
or hyperparameter fitting enter this transformation. Target means inherit the
prior correction's assumptions and transductive use of the full cohort.

| Unchanged diagnostic | Translation-only result |
| --- | ---: |
| All / eligible annotation comparisons | 471 / 146 |
| B intermediate failures | 5 / 8 eligible |
| CD14 Mono failures | 0 / 8 eligible |
| Coarse response contrasts passing | 15 / 16 |
| Sensitive mixed-program comparisons passing | 14 / 14 |
| Sensitive within-library program comparisons passing | 28 / 28 |
| Mean conditional donor scatter fraction | 0.072873 |

All 18 insensitive mixed-program and four insensitive within-library comparisons
remain in the evidence. Complete program flags remain false; sensitive-subset
passes are not complete qualification. The B-cell IFNa coarse contrast fails its
individual-fold preservation rule: maximum balanced-accuracy loss is 0.103212,
even though within-donor response vectors are preserved to numerical tolerance.
Donor scatter is closer to original PCA (0.079429) than native correction
(0.020288) or rigid projection (0.028972); this diagnostic cannot distinguish
technical batch from donor biology.

This ablation separates the effect of rotation from donor-mean translation.
CD14 Mono annotation failures resolve, while B intermediate failures and a
response classification failure remain. It does not show that translation
provides adequate integration or that a label-specific correction would generalize.
The next method must improve alignment without trading away these complete gates.

## Execution and reproduction

The protocol and input hashes were frozen before fitting; output hashes were
frozen before scoring. All coordinates were checked against the separately
ordered centered expression, maximum error 3.55e-15, with centered-norm error
8.98e-16. This algebraic cross-check is not an independent scientific reference.
The existing SVD annotation oracle passes all 114 folds, original PCA confusion
reproduces exactly, and the within-library program oracle passes. This was an
actual Python development run on inspected HIRISA data, not native-product or
prospective validation. No product default changes.

`python3 verify.py` checks every archive member hash, full-cohort gate counts and
output binding. [summary.json](summary.json) and [manifest.json](manifest.json)
record outcomes, environment and dependencies. The archive retains the frozen
protocol, fit/evaluation scripts, models, all scored comparisons and oracle
results, prior native results, logs and terminal status. The 516,030,080-byte
coordinate payload is retained separately on the Mac mini and hash-bound in the
manifest. Original HIRISA source, PCA and reference inputs remain external
study dependencies; this small result archive is not a standalone raw-data replay.

For fresh execution, restore the manifest-bound study inputs, extract the
archive scripts into a fresh directory, and run `fit.py STUDY OUT`, then
`evaluate.py --study STUDY --root OUT --python PYTHON`. Retain the protocol and
adjust archived absolute dependency paths for the host. Do not overwrite this
completed evaluation or reuse it as an untouched validation cohort.
