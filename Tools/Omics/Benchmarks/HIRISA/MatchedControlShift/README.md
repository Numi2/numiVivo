# Matched-control donor shifts: failed full-cohort correction

Equal-weight matched-control shifts fail **34/146** eligible annotation
comparisons, including **11/29 rare** comparisons. Mean eligible recall is
0.647907. The preceding native-target translation fails 22/146 (five rare);
rigid projection fails 31/146 (nine rare). This candidate is not promoted.

Shift estimation uses **403,449 control cells** across five preparations and
five donors. Controls are `none` for Bcell, Monocyte, NK and Tcell and
`culture_no_stim` for PBMC; Fresh PBMC is not the cultured intervention control.
For each preparation, compute each donor's original PCA mean and the equal-donor
grand mean. Average grand-mean-minus-donor-mean equally over preparations and
apply that one shift to all cells from the donor. Every donor/control group has
observations. No treated coordinate, native corrected coordinate, author cell
label, rotation or outcome-selected hyperparameter enters shift estimation.

The original PCA uses the full cohort. Observed query-donor controls enter the
shift, so this is transductive integration, not unseen-donor prediction. Equal
weights remove raw cell-yield weighting across preparations and donors; they do
not remove cell-composition differences inside each control group or distinguish
technical nuisance from real donor biology.

| Complete unchanged evaluation | Result |
| --- | ---: |
| Cells / components / annotation folds | 1,612,594 / 20 / 114 |
| All / eligible annotation comparisons | 471 / 146 |
| Native annotation failures resolved / new failures | 21 / 18 |
| Coarse response contrasts passing | 15 / 16 |
| Sensitive mixed-program comparisons passing | 14 / 14 |
| Sensitive within-library comparisons passing | 28 / 28 |
| Mean conditional donor scatter fraction | 0.083336 |

Donor scatter exceeds original PCA's 0.079429. All insensitive and insufficient
comparisons remain retained; complete program flags remain false. Within-donor
response vectors are preserved to numerical tolerance, yet response
classification and annotation preservation still fail. Thus eliminating treated
values from shift estimation does not make a global additive donor correction
adequate for this cohort. The evidence motivates a method that can address local
composition/alignment, with complete preservation gates, rather than further
selection among global shifts based on inspected outcomes.

## Verification and reproduction

The protocol and input hashes precede fitting, and output hashes precede
scoring. All transformed coordinates pass the algebraic translation check.
A separate indexed-accumulation script exactly reproduces the 25 control counts,
means and donor shifts. All 114 annotation SVD oracle folds, original baseline
reconstruction and the within-library program oracle pass. These are numerical
checks; the biological failures stand. This is an executed Python development
experiment on inspected data, not native product or independent qualification.

`python3 verify.py` checks every archive member hash, gate accounting, output
binding and the control-mean oracle receipt. [summary.json](summary.json) and
[manifest.json](manifest.json) record results and provenance. The archive retains
protocol, fitting/evaluation/oracle sources, model values, all comparisons, logs
and completion status. Internal evaluator filenames use `translation` because
the same unchanged translation harness is reused. The full 516,030,080-byte
coordinate payload is retained on the Mac mini and bound by the manifest;
original HIRISA source, PCA and program references are external dependencies.

Restore those dependencies and extract the scripts into a new output directory
to run `fit.py STUDY OUT`, followed by `evaluate.py --study STUDY --root OUT
--python PYTHON` and `check_control_means.py`. Adjust the latter script's retained
absolute study path for the host. Keep the frozen protocol; do not overwrite this
completed run. Fresh independent validation remains necessary before promotion.
