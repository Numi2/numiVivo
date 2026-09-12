# Nested training-target regularization: complete, unchanged predictions

All 150 outer target/gemgroup folds select **lambda 1**, the original fixed value,
from 0.01, 0.1, 1, 10 and 100. Selection minimizes equal-inner-target mean
full-33,694-gene RMSE after clipping, with all 29 inner targets retained. Exact
loss ties choose the larger lambda. Every outer target is removed from the
count matrix before the selection calculation, with the selected matrix hash
checked against the original frozen exclusion record.

The efficient leave-one-target-out weight calculation agrees with explicit
28-target augmented-system solves for every inner fold and grid value; maximum
absolute weight error is 8.88e-16. This is numerical qualification of selection,
not a fresh predictive result. All five original technical groups and all 30
targets remain included. The original fixed-model predictions and diagnostic
were previously inspected, so this is development reuse, not independent study
validation. Identical selected regularization offers no expected repair of the
identified target failures.

The complete native run passes all 601 commands: 150 fits, model replays,
predictions and prediction replays, plus one repeated prediction. All **150 full
reports are byte-identical** to the original fixed-lambda reports. Independent
NumPy and centered scikit-learn reconstruction passes all 750 vectors, with
maximum prediction error 1.78e-15 and weight error 3.33e-16. All five original
aggregate criteria still pass; 37/150 individual folds still lose to training
mean and 34/150 to no change. There is **no predictive improvement**.

All predictions were frozen before the new scorer read held responses. The run
uses the original qualified study executable, not a new full application build.
`python3 verify.py` checks the retained selection losses, source hashes, complete
execution record and original failure counts. The evidence archive includes the
executed checker and dependencies, command logs, plans, receipts, full scores,
report-comparison hashes and setup failures. Large inputs, models and prediction
arrays remain externally hash-bound in the recorded remote workspace. Retained
evidence verification does not execute a new prediction experiment.

This closes the proposed regularization search. Changing lambda within this
predeclared grid does not repair the target failures; repeated tuning on these
same outcomes cannot establish generalization. Further model changes need a
specific hypothesis and training-only selection, followed by a new untouched
study for an independent claim.

`protocol.json` and `select_regularization.py` were written before this selection
execution; `selection.json.gz` binds protocol, selector, original input freeze,
all selected-count hashes and all inner losses. Scripts currently use the
recorded remote study paths. The native driver imports the study's original
`prepare_training.py` and also requires a sibling copy to hash its provenance.

Retained setup failures: `attempt-1/` contains the initial selector named
`select.py`, which shadowed a Python standard-library module and failed import;
`native-attempt-1/` and its log preserve the driver stopping before native fitting
because the sibling provenance helper was absent. The renamed selector and
copied unchanged helper resolve these setup failures without changing the
selection rule or data.
