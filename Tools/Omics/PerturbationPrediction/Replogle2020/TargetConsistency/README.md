# Replogle target consistency across technical groups

The original fixed GO predictor improves average all-gene response RMSE in each
of five gemgroups, but three targets are worse than the training-mean baseline
in **all five**: sgSCYL1, sgSRP68 and sgSRP72. sgSRP72 also loses to no change in
all five. This retrospective diagnostic retains all 30 targets and 150 folds.
It does not fit a new model or constitute independent biological replication.

| Target losing to training mean in every group | Mean GO RMSE | Mean baseline RMSE | Relative excess error |
| --- | ---: | ---: | ---: |
| sgSCYL1 | 0.185188 | 0.182529 | 1.46% |
| sgSRP68 | 0.192492 | 0.186869 | 3.01% |
| sgSRP72 | 0.171239 | 0.159803 | 7.16% |

The table averages RMSE equally over the five technical groups for each target;
relative excess is the ratio of those averages. It is not a confidence interval.
SRP72 mean error also exceeds no change by 5.37%.

| Number of groups with worse GO error | Targets worse than training mean | Targets worse than no change | Targets worse than fixed shuffle |
| ---: | ---: | ---: | ---: |
| 0 | 19 | 17 | 21 |
| 1 | 1 | 3 | 4 |
| 2 | 2 | 3 | 0 |
| 3 | 3 | 4 | 0 |
| 4 | 2 | 2 | 1 |
| 5 | 3 | 1 | 4 |

SCYL1, SEC61A1, SRP72 and TELO2 lose to the fixed shuffle in every group.
These comparisons do not establish a mechanism, guide specificity or a causal
reason for failure. Five capture groups share the experimental system. Do not
interpret repeated technical behavior as five independent biological studies.
The original modest aggregate PASS and all failed individual folds remain valid.

## Verification and next experiment

`diagnose.py --output NEW_JSON` uses the published source scoring archive through
Git, verifies stored/decoded/member hashes, requires all 750 full-gene method
scores, and reproduces every original group-level failure list. It reports every
target and all three comparator histograms. No expression vector is changed or
newly scored. The first diagnostic attempt tried to parse archived log text as
JSON and failed before output; the corrected reader verifies all member hashes
but parses only JSON members.

The next useful model experiment is nested training-target-only regularization
selection, retaining all 30 outer held-target folds and all five groups. The
outer target must be excluded from every inner fit and selection loss. Compare
the resulting native predictions with the fixed lambda-1 model and unchanged
baselines, freeze predictions before outer scoring, and retain all failures.
This is now development reuse of Replogle, so even a gain would require another
untouched study before an independent generalization claim. Do not exclude these
repeatedly failing targets or select their hyperparameters using their outcomes.
