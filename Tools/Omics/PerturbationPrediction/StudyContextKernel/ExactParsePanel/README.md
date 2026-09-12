# Exact Parse feature-panel development evaluation

The separately versioned `context-parse-exact11600-v1` panel uses all 11,600
exact symbol matches between the original 11,800-feature context model and the
Parse feature axis, in original model order. The rule and coordinates were frozen
before fitting this revision. No aliases, zero filling, response-based feature
selection, donor changes or library renormalization were applied. The original
model remains unchanged. This panel is compatible by exact symbol; that does
not prove sequence or assay equivalence.

All 75 donor predictions were executed with the retained optimized native Swift
context-kernel binary. Each outer fold excludes the entire query study's treated
outcomes; penalty selection uses only the other studies. The three studies were
previously inspected, so this remains development reuse, not untouched external
validation. No Parse response was fitted or scored.

| Held-out study | Donors | Gain over training mean | Gain over no change | >=5% over both |
| --- | ---: | ---: | ---: | --- |
| HIRISA | 5 | 21.154% | 22.793% | PASS |
| Kang | 8 | -0.263% | 6.993% | FAIL |
| GSE181897 | 62 | -0.309% | 4.338% | FAIL |

**Not promoted.** Exact panel compatibility does not repair the observed transfer
failure. All folds select penalty 1 using the unchanged training-only procedure.
No donor failure or study is hidden by a pooled mean. The endpoint is equal-donor
mean RMSE of average log1p CPM RNA profiles, with full-source RNA denominators.
This does not establish protein, cellular phenotype, tissue or clinical outcomes.

Independent NumPy eigensystem checks verify every inner loss, response weight
and all 870,000 predicted values; maximum prediction difference is 3.56e-15.
A separate scalar check verifies all 225 donor/baseline scores; maximum
difference is 4.16e-14. Exact feature coordinates were independently reconstructed
from metadata. These checks establish computation, not biological accuracy.

## Evidence and reproduction

[Summary](summary.json), [archive](evidence.tar.gz) and [member manifest](manifest.json)
retain the panel freeze, projected cohorts, all fold inputs, native binary/source,
predictions, numerical checks and donor-level scores. Every archived member was
independently checked against its size and SHA256. Execution used an existing
optimized standalone CPU binary, not a fresh full-product build.

After extracting into a new directory, run `python verify.py`, `python score.py`
and `python verify_scores.py` with NumPy installed. For a fresh native run, make
a separate copy without `native/` and `prediction-freeze.json`, then run
`python run.py`. That driver freezes predictions before numerical verification
and scoring. To regenerate preparation, restore the original cohort dependencies
identified by SHA256 in `panel-freeze.json`, adjust only their local paths, and
run `prepare.py` in a new directory. Do not overwrite the retained execution.

The [Parse primary dose/reagent](../../ParseIFNB/BCellCounts/NEXT_STEPS.md#primary-dose-and-reagent-resolved)
is now known. Cross-study exposure comparability and participant overlap remain
unverified.
This revision supplies an explicit panel contract and its failed development
evaluation; it does not authorize describing a future Parse result as protocol-
matched independent validation.
