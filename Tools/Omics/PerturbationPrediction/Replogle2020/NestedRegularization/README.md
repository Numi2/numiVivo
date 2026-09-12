# Nested training-target regularization: selection complete, native run pending

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

The native driver freezes all selected plans before issuing actual fit, model
replay, prediction and prediction replay commands. Its complete 150-fold run is
in progress at `/Users/n/numivivo-replogle-nested-20260912/native`; no new outer
score or completed-run claim is made here. Full native predictions must be frozen
before the independent scorer reads held responses. The current run uses the
original qualified study executable, not a new full application build.

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
