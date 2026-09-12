# Parse context prediction: failed transfer test

The frozen context candidate fails on all twelve Parse donors: mean donor RMSE
is 0.375746, versus 0.310710 for training-mean response and 0.284361 for no change.
It worsens error by 20.93% and 32.14%, respectively, failing the fixed requirement
of at least 5% improvement against both baselines. It is worse than training mean
in all twelve donors and worse than no change in eight. No donor was excluded.

This evaluates 139,200 native predictions across the exact 11,600-feature panel.
Training uses all 75 donors from GSE181897, HIRISA and Kang; penalty selection
holds out each training study and selects 1. Parse queries contain only PBS
means from the frozen 72,446 B-cell membership. Full-source per-cell log1p CPM
normalization precedes panel projection; no aliases or renormalization are used.

The protocol was frozen before prediction. Parse treated values had previously
been accessed for engineering verification, which is explicitly disclosed;
prediction hashes were frozen before outcome scoring. No outcome-driven model,
feature, cell-selection or acceptance-rule changes were made. This is an external
study transfer stress test. Participant independence and matched exposure remain
unverified, so it is not protocol-matched independent validation.

Native predictions, response weights and all training-only inner losses match
an independent NumPy eigensystem calculation; maximum prediction difference is
6.44e-15. All 36 donor/baseline errors also match scalar compensated summation
within 5.56e-17. Exact source feature symbols and PBS query values were checked
again before scoring. Numerical correctness does not repair the biological failure.
The candidate remains unpromoted, also retaining its two failed development-study
transfer gates. Parse is now an inspected evaluation cohort, not an untouched
future test set.

## Retained evidence

[Protocol](protocol.json), [all donor scores](scores.json), and
[numerical verification](verification.json) summarize the execution.
[evidence.tar.gz](evidence.tar.gz) retains all 17 study files, including the
native source/binary, exact training arrays and metadata, panel contract,
PBS-only prediction input, native output, freeze receipts, scoring targets and
execution/verification/scoring scripts. [manifest.json](manifest.json) binds
all members; every member passed checks on both hosts.

The archive is 11,616,098 bytes with SHA256
`0098361925c553d9520e5386d1a9322d6aecb00396f9e593211178eb2399fd1b`.
Original count-source and normalization evidence remain separate dependencies
linked from [the normalization replay](../../StreamedLogCPM/ParseReplay/README.md).

To recheck numerical predictions, extract into a new directory and run
`python verify.py` with NumPy installed. The original run.py and score.py retain
execution-host source paths and refuse to overwrite their completed outputs;
fresh source replay requires restoring those dependencies and using a separate
study directory. Do not overwrite this retained failed evaluation.

The next modeling decision should address the observed cross-study failure
using development cohorts and explicit exposure/population information. Any
revision using Parse outcomes is development reuse and needs a new untouched
cohort before an independent transfer claim.

The [retrospective error decomposition](Diagnosis/README.md) identifies excessive
response magnitude in all twelve donors together with imperfect directional
alignment. It does not fit or validate a correction.

The [frozen weight audit](WeightAudit/README.md) confirms equal study base weights
and nonnegative response weights for every query. A nonnegative-weight constraint
would not repair this failure.
