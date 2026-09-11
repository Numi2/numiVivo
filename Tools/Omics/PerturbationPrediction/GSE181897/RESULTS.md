# Complete external GSE181897 B-cell prediction test

**HIRISA-trained mean passes the frozen primary target; Kang-trained mean fails.**
All 62 eligible external donors were predicted from control RNA only using each
training origin. HIRISA mean reduces equal-donor RMSE by 5.62% and beats no-change
in every donor. Its nominal 95% intervals nevertheless severely under-cover.
This is a bounded external RNA point-prediction result, with model-dependent
transfer and failed uncertainty. It does not establish general biological outcomes.

## Complete point-prediction results

The unchanged [protocol](PROTOCOL.md) requires at least 5% lower equal-donor,
all-panel RMSE versus no-change for the mean response, assessed separately for
both training origins. Ridge's secondary criterion is improvement over both
no-change and mean response. Neither origin was selected after scoring.
RMSE uses natural-log(1+CPM) over all 11,800 panel features.

| Training origin | No change | Mean response | Median response | Context ridge | Mean gain vs no-change | Mean primary | Ridge secondary |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| Kang, 8 donors | 1.121616 | 1.131581 | 1.121135 | 1.110158 | −0.89% | FAIL | PASS |
| HIRISA, 5 donors | 1.121616 | 1.058543 | 1.058579 | 1.079193 | 5.62% | PASS | FAIL |

Kang mean is worse than no-change in **41/62** donors; Kang ridge is worse than
no-change in 21/62 and worse than mean in 2/62. HIRISA mean and ridge both beat
no-change in all 62 donors, but ridge is worse than mean in **62/62**. Thus the
experiment does not pass both origins' primary requirements. Every donor's RMSE,
MAE, clipping count and source stratum size is retained in
[scores.json.gz](evidence/2026-09-11-prediction/scores.json.gz), alongside all
independent comparisons and interval assessments.

## Predictive intervals

The frozen normal-donor model uses Student-t critical values and the
future-observation factor `sqrt(1+1/n)`. Intervals concern mean-response prediction
for a future donor, conditional on its measured control; they are not confidence
intervals for the mean improvement or intervals for ridge.

| Training origin | Available / total genes | Mean treated coverage | Donor coverage range | Mean treated width | Mean unclipped-response coverage |
| --- | ---: | ---: | --- | ---: | ---: |
| Kang | 10,496 / 11,800 | 92.39% | 79.14–97.47% | 4.4861 | 92.39% |
| HIRISA | 11,799 / 11,800 | 35.54% | 30.95–42.84% | 0.3097 | 31.01% |

Widths use log1p(CPM) units. Constant-response features have unavailable
intervals: 1,304 for Kang and one for HIRISA. The coverage denominator includes
only available features; all missing features and below/above misses are retained.
Kang intervals are much wider and some donors under-cover. HIRISA's narrow
intervals fail badly in this external context despite its useful point estimates.
Neither nominal 95% nor aggregate improvement supplies calibrated uncertainty.

A subsequent [complete cell-sampling diagnosis](Uncertainty/README.md) finds
substantial observed sampling variability and poor zero-count coverage. Adding
variance estimated from observed control **and treated** cells raises HIRISA
treated coverage to 94.25%, but this is an outcome-informed diagnostic, not a
predictor or independent validation. It leaves the results above unchanged and
identifies the need for a count-based observation model.

## Primary treatment identity

The first audit found incomplete current author code and an external curator's
uncertain control label. The original count receipts remain unchanged with
unqualified condition meanings. Historical primary code now resolves that gap:

- [Author IFN comparison](https://github.com/yelabucsf/clue/blob/5e1f3e2a3ec9ddbf4d594e7b0c0a2f5f2a5fe6f0/prod/ifn_compare/ifn.ipynb), notebook cells 11/12/16: the ordered B/G comparison is labelled IFN-beta/IFN-gamma.
- [Author all-stimulation preprocessing](https://github.com/yelabucsf/clue/blob/5e1f3e2a3ec9ddbf4d594e7b0c0a2f5f2a5fe6f0/prod/deseq/fix_read_counts.ipynb), cells 12/13: every stimulation is compared with C as control. An older production notebook uses the same mapping; metadata lists A/B/G/R/P as stimuli, excluding C.

The primary commit is `5e1f3e2a3ec9ddbf4d594e7b0c0a2f5f2a5fe6f0`, tree
`e23ff74d73ee97c68b8cdca444e7500e30ecfcff`, before the repository reset removed
production notebooks from the current tree. The audit verifies each complete
download against its Git blob and records SHA256, URL and zero-based notebook
cell numbers. Notebook code was inspected; author notebook computations were
not executed or used as prediction outputs. No RNA pattern or prediction score
determined the roles. This establishes the analysis-code condition mapping,
not a new measurement of reagent exposure or validation of cell identities.

## Frozen cohort and inputs

The original source and complete count qualification are in [README](README.md).
The source describes 500 IU/mL IFN-beta for nine hours and matched-incubation
controls. All 62 donors with positive B/C RNA libraries enter the test, comprising
4,938 author-annotated B-lineage cells, with 4–122 cells per donor/condition.
Donors 5 and 23 lack B and remain excluded with their original counts reported.
No new QC, abundance threshold or outcome-dependent donor filter was introduced.
Treated author labels condition this endpoint; this is not prospective cell-type
identification. Pools are technical units, not extra biological replicates.

The panel keeps all 11,800 exact matches to the previous 11,884-symbol panel and
retains the 84 absent features. RNA normalization uses all 20,303 original query
RNA features and each training study's full original library denominator.
All eight Kang and five HIRISA training donor pairs are reused under unchanged
model settings. Training/query H5AD rows are explicitly labelled donor aggregates;
they are not represented as individual cells. Queries contain only controls.

Source context differs: Kang uses lupus-donor PBMC B cells after six hours;
HIRISA uses healthy-donor enriched B cells after 21 hours at 100 units/mL and a
different assay. No temporal or dose model was added. Participant overlap is
unverified and some investigators are shared. The result is a separate data
collection, not independent-laboratory validation or proof of biological equivalence.

## Execution, checks and retained failures

The actual qualified native executable ran on the physical Mac mini. Its SHA256
is `0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b`;
HDF5 library SHA256 is `a00ffbf8ab94ad81f67231a1ae01df748689e1c35a3615a57f88f4710b8d213e`.
This reuses that historical implementation; it is not a new full-package build
qualification. Inputs froze before fitting and all predictions froze before
scoring. Four-donor publications bounded temporary storage while retaining all
62 donors for each origin.

All 124 native predictions reconstruct. Independent NumPy/SciPy/scikit-learn
checks validate both fits, 496 estimate vectors, 248 interval assessments and all
379 original native source aggregates. Maximum point and bound differences are
`9.82e-14` and `3.56e-15`. All six score/check files repeat byte-exactly.

The initial run stopped verifying Kang chunk 14 because temporary reconstruction
exhausted disk space. Its published prediction and failed log were preserved;
after removing exact committed working-copy duplicates, only that verification
was retried and remaining chunks continued. Earlier completed predictions were
not refitted. A metadata API tree-alias assertion and a first local checker call
before artifact transfer finished are also retained; neither read outcome scores.

The complete 366 logical native files occupy 671,133,531 bytes when expanded.
They are retained as 110 independently decoded, hash-verified objects totaling
154,174,098 compressed bytes. One complete four-donor prediction bundle from
each origin was additionally restored and verified natively. This is separate
from the original native reconstruction of all predictions.

## Retention and reproduction

The [compact evidence manifest](evidence/2026-09-11-prediction/manifest.json)
retains executed scripts, source-role audit, complete plans, all scores/checks,
both execution attempts, logs and artifact identities. Full native objects and
H5AD inputs remain at **both** `/Users/home/numivivo-gse181897-20260911` and
`macmini:/Users/n/numivivo-gse181897-20260911`. The manifest lists every external
file and hash. The compact repository archive alone cannot restore these files.
Original author notebooks remain source-hash-bound locally and retrievable from
the pinned author commit; they are not republished in the evidence archive.

The recipes below use the already qualified study layout and scientific Python
environment. Existing completed outputs must be preserved; use a separate copy
with new output destinations for a fresh run. `resume_prediction.py` is the exact
one-time disk-failure recovery recipe, not a generic resume command.

```sh
python prepare_prediction.py
# On the physical Mac mini, after transferring frozen prediction-inputs/:
python3 run_prediction.py
# After complete native results have been copied back:
python score_prediction.py --out new-scores
python restore_prediction.py --native prediction-execution \
  --bundle prediction-HIRISA-00 --out restored-HIRISA-00
numivivo singlecell-perturbation-prediction-verify restored-HIRISA-00
```

The source protocol remains unchanged. GSE181897 is now scored development
knowledge for subsequent changes: improved uncertainty or model selection using
these results needs a new validation experiment before stronger claims.
