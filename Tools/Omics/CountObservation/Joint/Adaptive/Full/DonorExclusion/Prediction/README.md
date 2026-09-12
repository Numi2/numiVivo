# Control-only joint prediction in omitted donors

The first two complete Kang donor omissions improve RNA point-prediction error
on their eligible genes. These are reused development donors, not an independent
biological validation cohort. This retained September 12 snapshot covers two of
13 folds; the ongoing run is not fully qualified. No uncertainty-calibration
claim follows from these fits.

| Omitted donor | Predicted / source genes | Joint RMSE | No-change RMSE | Training-mean RMSE | RMSE gain vs no change / training mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kang:patient_101 | 8,353 / 15,706 | 0.931587 | 1.346774 | 1.252205 | 30.83% / 25.60% |
| Kang:patient_1015 | 7,790 / 15,706 | 0.592070 | 0.747017 | 0.751154 | 20.74% / 21.18% |

RMSE uses `log1p(CPM)`: the prediction is the log of one plus the conditional
mean treated rate, and the outcome is measured treated pseudobulk CPM using
full measured-library depths. No change uses measured query-control CPM. The
training-mean baseline adds the average training-donor log response to query
control, clipped at zero in log space. Every model is compared on identical
eligible genes. The complete-run development threshold was frozen at 5% pooled
RMSE gain against both baselines, with every worse donor reported separately;
these two folds do not establish an all-fold pass.

The folds retain respectively 7,351 and 7,911 unavailable dispersions, plus two
and five leaf-budget limits. None is filled with a successful prediction. There
are 1,103 and 1,230 latent point-mass predictions. Fitted point masses do not prove
zero biological or parameter uncertainty. Historical external treated-interval
undercoverage remains unresolved.

## Execution and checks

`Run.swift` reads a frozen adaptive model and the omitted donor's control-count
histogram. It checks donor identities and calls the native adaptive predictor;
no query-treated counts or full-cohort fitted weights enter prediction. A fixed
one-cell 10,000-read plan supplies a standardized count-moment scenario, separate
from observed treated libraries and latent RNA-rate scoring. Prediction is
streamed in 64-gene shards using the existing qualified native library.

The individual-cell SciPy reference reconstructs input histograms and checks all
posterior weights, rate moments, planned-count moments, source fingerprints and
unavailable states. Across the two folds, **9,887,511** comparisons pass with
maximum scaled discrepancy **1.28e-12**. Wrong features, training-donor queries,
wrong model donors and duplicate count bins are each rejected before prediction.
After all prediction bytes are hashed and checked, a separate scoring process
opens outcomes. A second scorer uses SciPy CSC sums and scalar reductions to
check every outcome, baseline and aggregate metric. An initial scorer mistook a
root feature dataset for a donor group; the structural fix and failed-stage note
are retained, without changing the scientific endpoint.

`advance.py` advances verified complete fits through prediction, numerical checks,
scoring and scoring checks. It uses a separate process from fitting and checks
the fitting owner's actual command while waiting. Partial output is retained on
failure. The execution copies are frozen under the study root; do not replace
those scripts or binaries during an active run.

## Reproduction

Build a new executable with `build.sh QUALIFIED_BUILD OWNER_SOURCE NEW_BINARY`.
`runtime-freeze.json` must bind that executable's hash and path. Freeze the
scoring plan before outcome scoring. `run.py FIT_STUDY PREDICTION_STUDY TAG`
requires a fully independently checked fit. Then run `verify.py`, `score.py`
and `verify_scores.py` with the same three arguments. `advance.py` applies that
sequence to all declared folds as they finish.

The [retained two-fold evidence](evidence/2026-09-12-kang-two-folds/retention.json)
contains complete fit and prediction shards, manifests, numeric checks, scores,
recipes and the exact prediction executable. The archive deduplicates bytes and
uses parts below GitHub's file-size limit. Source caches remain in the published
full-cohort archives and are identified by their original hashes.

Restore with `restore.py EVIDENCE_DIRECTORY NEW_DIRECTORY FULL_CACHE_DIRECTORY`.
Then run the restored `recipes/verify.py` and `recipes/verify_scores.py` using
`NEW_DIRECTORY/fits`, `NEW_DIRECTORY/prediction` and each retained fold tag.
Restoration validates every object and logical file. The cache links and restored
files are shared hard links and must be treated as read-only.

## Complete-run summary

`summarize.py FIT_STUDY PREDICTION_STUDY [OUTPUT]` validates score-report hashes,
the frozen endpoint plan, prediction shard hashes and independent scoring
receipts before pooling errors. It reports pooled and per-origin RMSE, MAE,
excluded cases, point masses and every worse donor. Incomplete groups retain a
null development pass result, even if their partial gains exceed 5%. Missing
folds remain listed against the original 13-fold protocol. No complete-run or
independent-biological-validation claim is inferred from completed subsets.

The summary checker reconstructs aggregate metrics directly from every scored
gene row and checks that an empty run cannot pass, an altered score is rejected,
and an altered endpoint plan is rejected. The retained software-validation
snapshot covers 24,490 gene-folds from the first three locally scored Kang folds;
it is a partial progress summary, not a replacement for full-fold evidence.


## Five-fold progress snapshot

The [retained reports](evidence/2026-09-12-five-fold-progress/manifest.json)
cover three Kang and two HIRISA omissions. Both HIRISA donors lose to the
training-mean response: RMSE is 57.5% worse for HIRISA-00 and 74.4% worse for
HIRISA-01. Their pooled joint RMSE is 0.143048 versus 0.086439 for training mean
and 0.332588 for no change. The three Kang folds improve pooled RMSE by 26.1%
versus no change and 22.5% versus training mean.

Pooling all five gives gains of 28.1% and 21.6%, respectively, but masks the
HIRISA losses. All groups retain a null full-run criterion because eight of the
thirteen declared folds are pending in this snapshot. The new directory retains
reports and hash bindings; complete raw shards remain in the live study pending
full archival publication. It does not replace the earlier complete two-fold
archive or establish independent biological validation.
