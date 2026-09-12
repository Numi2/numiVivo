# Training-only response shrinkage: transfer remains failed

The native candidate jointly selects response scale (1, .75, .5, .25, 0) and
penalty (.01, .1, 1, 10, 100, infinity) using training-study holdouts.
Exact ties prefer smaller scale, then stronger penalty. Response scaling precedes
nonnegative treated-expression clipping; both original baselines are unchanged.

| Held-out study | Donors | Scale | Gain over training mean | Gain over no change | Gate |
| --- | ---: | ---: | ---: | ---: | --- |
| HIRISA | 5 | 0.75 | 30.114% | 31.566% | PASS |
| Kang | 8 | 1 | -0.263% | 6.993% | FAIL |
| GSE181897 | 62 | 0.5 | -0.422% | 4.231% | FAIL |

All folds select penalty 1. All 75 donors and 11,600 frozen features remain.
The fixed gate requires at least 5% lower equal-donor mean RMSE against both
baselines in every study. Two studies fail; the candidate is not promoted.
No Parse input or outcome enters this execution. The design was motivated by
its inspected failure, and these three cohorts are reused development data.
This is not independent biological validation.

The standalone native Swift candidate compiled and executed on the Mac mini.
Independent NumPy eigensystem calculations verify all inner losses, selection,
response weights, predictions and baselines. Every donor score also matches
scalar compensated summation within 1e-12. No product default was changed;
a full-product build is outside this qualification.

[verified-scores.json](verified-scores.json) contains all donor scores and numerical
comparisons. [protocol.json](protocol.json) records the frozen design.
[evidence.tar.gz](evidence.tar.gz) retains all 17 study files, including source,
binary, inputs, outputs, freeze receipts, driver and verifier/scorer.
[manifest.json](manifest.json) binds every member. All hashes passed on both hosts.
The archive is 19,997,005 bytes, SHA256
`f073ebedcf7e1c79abba357293cc738c9af1eb9af835b7da429ed76eb9e0bbcd`.

Scoring uses separately bound cohort arrays from the prior exact-panel study.
Scripts retain execution-host paths and reject overwriting completed outputs.
Restore those dependencies in a new study directory for reproduction; this
archive alone is not a standalone raw-count replay.

Shrinkage improves HIRISA but does not resolve study transfer. Further development
needs exposure/population-aware training and an untouched evaluation cohort
before another independent transfer claim.
