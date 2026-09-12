# Response transport candidate: rejected for promotion

This development experiment follows inspection of the completed Kang/HIRISA joint-count results. It is not independent validation and does not replace their frozen benchmark scores.

For each learned atom (control c, treated t), transport its log1p response to query control pseudobulk CPM q:

    transportedCPM = max(0, (1 + q) * (1 + t) / (1 + c) - 1)
    prediction = log1p(sum(controlOnlyPosteriorWeight * transportedCPM))

Weights are the original frozen control-only conditional probabilities. Each atom is clipped before averaging in CPM. No new parameter was selected or fitted. This is a deterministic transport heuristic, not the original model posterior, a refitted count likelihood, or calibrated uncertainty. Query control is treated as fixed. The existing training-only fits and their dispersion restrictions remain prerequisites.

The candidate and endpoint were frozen before this candidate's predictions or scoring. The original outcomes had already been inspected; the freeze does not make this an untouched test. Prediction code reads only existing input and posterior shards. Scoring is a separate step against the same conditionalPrediction population as the original study. All 13 donors are retained.

| Dataset | Transport RMSE | Original joint RMSE | Training-mean RMSE | Transport gain over mean |
| --- | ---: | ---: | ---: | ---: |
| Kang, 8 donors | 1.222553 | 0.999211 | 1.311326 | 6.77% |
| HIRISA, 5 donors | 0.100816 | 0.158447 | 0.090109 | -11.88% |

HIRISA improves relative to the original joint predictor but loses to the training-mean baseline in every donor. Kang loses accuracy relative to the original joint predictor. Do not promote this candidate as a replacement or combine winners by dataset after seeing these outcomes. A revised model needs training-only selection and fresh independent evaluation. These results narrow the next task: coordinate transport alone is insufficient; investigate estimation and regularization of conditional donor responses.

`test.py` checks response identities, clipping, averaging and rejection of invalid inputs. `verify.py` independently evaluates the rate-ratio expression and binds original input/output hashes to the frozen verification reports. `score-attempt1.json` preserves an initial population-assertion failure; the correction explicitly selects the frozen eligible status and does not alter predictions.

`evidence.tar.gz` retains scripts, freeze, all candidate predictions, input/posterior hash bindings, numerical verification, scores and logs. It does not duplicate the original raw fit and posterior shards. Reproduction requires those retained original study directories. Run `test.py`, then `predict.py SOURCE NEW_DEST`, and the separate verification/scoring scripts with their study paths adjusted to the restored locations.
