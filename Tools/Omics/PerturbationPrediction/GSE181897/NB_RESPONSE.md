# GSE181897 native negative-binomial response run

**Status on 2026-09-14: the complete native NB2 run executes and verifies, but
it does not meet the frozen 5% transfer target for either training origin.** It
adds a count-based conditional RNA estimate; it does not establish prediction
of a general biological outcome.

This is a separate, opt-in run on the already frozen GSE181897 inputs. The
training models use 8 paired Kang donors or 5 paired HIRISA donors. Each query
contains control RNA only. The 62 eligible donor treated rows are used only for
the final held-out comparison, not for fitting. The source aggregate retains
all 379 observed groups, 20,303 RNA features and the author-resolved
`source-code:B` IFN-beta / `source-code:C` control roles.

## Native result

The model uses `responseModel: "negativeBinomial"` with the native NB2
treatment effect, per-feature dispersion, standard error and status arrays. The
existing no-change, mean-response, median-response and context-ridge estimates
remain in every prediction report. RMSE is the equal-donor mean over all 11,800
panel features in natural-log(1+CPM) units.

| Training origin | NB tested / panel | No-change RMSE | Mean-response RMSE | NB2 RMSE | NB2 gain vs no-change | NB2 better than no-change | NB2 better than mean |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Kang (8 donors) | 4,695 / 11,800 | 1.121616 | 1.131581 | 1.117699 | 0.35% | 46 / 62 | 44 / 62 |
| HIRISA (5 donors) | 11,616 / 11,800 | 1.121616 | 1.058543 | 1.091140 | 2.72% | 62 / 62 | 0 / 62 |

The NB estimate is also scored on its identified features only: mean RMSE is
0.976306 for Kang and 1.097403 for HIRISA. Those denominators are different
from the full-panel comparison and are reported to expose feature support, not
to replace the preregistered all-panel gate. HIRISA's existing mean-response
baseline remains the strongest of these estimates; the NB2 route does not
improve it.

The complete run contains 124 predictions (62 per origin), 32 four-donor query
bundles, 620 estimate vectors and 32/32 native prediction verifications. The
fresh b7 source aggregate also passes native verification. The compact
[checks](evidence/2026-09-14-nb2/checks.json),
[per-donor scores](evidence/2026-09-14-nb2/scores.json),
[summary](evidence/2026-09-14-nb2/summary.json),
[execution manifest](evidence/2026-09-14-nb2/manifest.json) and
[query log](evidence/2026-09-14-nb2/query-run.log) retain the evidence without
committing the large H5AD, model or prediction bundles.

The four legacy baseline means reproduced in this score are unchanged from the
earlier frozen comparison to the displayed precision; only the new
`negativeBinomialEffect` row is additional.

## Support and batch handling

Kang reports 4,695 `tested`, 4,114 `filteredLowExpression` and 2,991
`rankDeficientSupport` features. HIRISA reports 11,616 `tested`, 89
`filteredLowExpression` and 95 `rankDeficientSupport` features. The model
stores `adjustForBatch=false` for both origins: Kang has one shared unreported
pool, while HIRISA's pool is nested within donor. Adding donor and nested-batch
columns would consume the paired design's residual degrees of freedom. The
published b7 guard omits that redundant adjustment and rejects other
non-identifiable designs instead of silently fitting a rank-deficient model.

The earlier f287 HIRISA attempt failed with
`design shape, finite entries or residual degrees of freedom`. That log remains
as [preserved pre-guard evidence](evidence/2026-09-14-nb2/pre-guard-hirisa-fit.log);
the successful b7 run does not overwrite it. Kang's NB effects are byte-identical
before and after the guard apart from the explicit batch-adjustment request
flag.

## Reproduction and interpretation

The native run used revision `b7c3fc942b341629786e47bf1739faea85f3525c` on the
physical Mac mini. The release binary SHA-256 is
`63ca5e565235b45a642e52fb70541e8fe1567814a824d38231cdd53f8b748dbb`; the native
HDF5 library SHA-256 is
`a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`; and the
fresh source aggregate report SHA-256 is
`4fc2f93b1165c70a32bac68e8295ff0979a67a78ec48fe87c0da1f43df41fbee`.
The reproducible scoring recipe is
[score_nb_response.py](score_nb_response.py); it validates source controls,
library totals, feature identity, donor coverage and every estimate before
writing scores.

The result is evidence for a bounded conditional molecular-response estimate in
a reused development comparison. It does not provide calibrated NB uncertainty,
an unseen-perturbation model, causal or mechanistic identification, or a link
from RNA to phenotype, disease progression, tissue function or treatment
benefit. A stronger claim requires a new frozen biological endpoint and
independent held-out validation.
