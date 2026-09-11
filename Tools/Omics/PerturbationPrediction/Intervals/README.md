# Native donor-response predictive intervals

**Native arithmetic and replay: PASS. General predictive calibration: not
established.** Optional nominal 95% intervals now execute on all 26 frozen
Kang–HIRISA cross/within donor folds. The existing four point predictions remain
exactly unchanged. This is an uncertainty-method assessment on previously
inspected studies, not new independent biological validation.

## Empirical results

All 11,884 source-matched genes remain in the panel. Genes with unavailable
intervals are retained in the denominator of availability, but cannot contribute
a coverage observation. Coverage and width below are averaged equally across
donors among each donor's available genes. They are not confidence intervals
computed by treating genes as independent donors.

| Query study / mode | Donors | Available genes | Raw response coverage | Treated coverage after clipping | Mean treated width |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kang / within | 8 | 87.65% | 92.94% | 92.94% | 4.230723 |
| Kang / cross from HIRISA | 8 | 99.98% | 33.08% | 42.39% | 0.279829 |
| HIRISA / within | 5 | 99.97% | 94.69% | 94.70% | 0.467843 |
| HIRISA / cross from Kang | 5 | 88.57% | 99.21% | 99.21% | 5.152851 |

Widths use natural log1p(CPM) units. Kang within-study donor coverage ranges from
**81.45% to 98.92%** after clipping; its mean near 95% hides substantial variation.
HIRISA within-study coverage ranges from 93.29% to 96.57%. Cross-to-Kang coverage
ranges from 36.73% to 51.02%, well below nominal coverage. Cross-to-HIRISA intervals
cover broadly but are much wider than its within-study intervals; high coverage
alone does not establish useful uncertainty. The complete donor scores retain
below/above misses, availability and both unclipped and clipped widths.

These studies differ in health status, stimulation duration, preparation and
assay chemistry. Normal, exchangeable donor responses are an explicit model
assumption, not an experimental finding. The [preceding point-prediction test](../CrossStudyIFNB/README.md)
and its failed cross-ridge comparisons remain unchanged. There is no tuning of
variance floors, features, coverage levels or correction factors using these
coverage results. Few biological donors, correlated genes and reused studies
prevent promoting these nominal intervals as generally calibrated.

## Contract and method

Set `donorResponseIntervalCoverage` in a native perturbation fit plan to request
an interval, for example `0.95`. Allowed finite values are 0.5 through 0.99.
Omitting it preserves the historical model and prediction serialization.
This option works with the default full feature universe or an explicit
`responseFeatureIDs` panel; full-library normalization remains unchanged.

For each gene, the native owner computes training-donor sample response variance
with n−1 degrees of freedom. Its interval is centered on the training **mean
response**, with half-width `t × s × sqrt(1+1/n)`. This is the
[NIST normal-model prediction limit](https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/predlimi.htm)
for one future observation, using the existing native Student-t quantile owner.
It includes future-donor variation as well as uncertainty in the estimated mean.
It is not the narrower confidence interval for the population mean.

`donorResponseVariances` stores optional per-gene sample variances in the model.
Exactly constant or zero-variance responses produce null entries, not zero-width
certainty. Each query's `meanResponsePredictiveInterval` retains the nominal
coverage, donor count, degrees of freedom, Student critical value, unavailable
feature indices, and four nullable bound arrays:

- `unclippedResponseLower` / `unclippedResponseUpper`.
- `predictedTreatedLower` / `predictedTreatedUpper`, computed by adding the
  observed query control and clipping below zero.

Clipping can create a point mass at zero, so transformed and raw coverage must
remain distinct. Query control measurement uncertainty is not modeled separately.
The intervals apply to the mean-response baseline only; they are not context-ridge
intervals, simultaneous coverage across genes, or Bayesian kinetic-model output.
There is no RNA-to-phenotype or clinical outcome qualification.

## Verification and identities

The [protocol](PROTOCOL.md) was declared before interval fitting/coverage scoring,
but after prior point-prediction outcomes were inspected. All 26 native fits,
model verifications, predictions and prediction verifications completed, plus
one exact repeated prediction: **105 successful native commands**.

Seven Swift tests pass, covering the existing panel behavior plus predictive
versus mean uncertainty, clipping, invalid coverage, constant-response
unavailability and omitted-field compatibility. Eight actual CLI regression
commands retain their expected success/rejection outcomes; the historical
18,082-gene HIRISA numerical model and four predictions remain exact when the
option is omitted. A compressed complete bundle was restored and both its model
and prediction were verified again with the qualified native executable.

NumPy 2.5.3 / SciPy 1.18.1 independently verify all sample variances, availability
masks, quantiles and **104 interval-bound arrays**. Maximum bound difference is
4.264e-14, maximum variance difference 3.553e-15, and maximum Student critical-value
difference 1.022e-14. All **104 prior point vectors** match exactly. Two complete
scoring runs produce identical checks, comparisons, scores and summaries.
The check also verifies every archived native file's bytes and hashes.

Physical execution: M4 Pro Mac mini, Swift 6.3.3, macOS 26.6. This qualifies the
scoped CPU Omics owner/CLI, not a full-application build or GPU performance.

- Native executable SHA256: `0f08cc423a10a2910962b25039c90731f2ceff013ec32e9a6df09d547a87ae4b`.
- Interval input freeze: `5c188899cfea5a7875333f62c7b54924ebd2a7911f77697acde9915f7470ed1d`.
- Interval prediction freeze: `f86389a2b206a24fdbd17e3f916ab8b25497b92b0ddaf015689cc128ab0878ff`.
- Independent scorer: `1ddef0beff146a25078967676808b595993ce24c6268f62c6e706cc45d663de0`.

## Reproduce and restore

Use the frozen [CrossStudyIFNB inputs](../CrossStudyIFNB/README.md). From the
repository root, build and run:

```sh
bash Tools/Omics/H5AD/build.sh /new/runtime --with-cli
bash Tools/Omics/PerturbationPrediction/CrossStudyIFNB/test.sh /new/runtime
NUMIVIVO_HDF5_LIBRARY=/path/to/libhdf5.dylib \
python Tools/Omics/PerturbationPrediction/Intervals/run.py \
  --inputs /path/to/cross-study/inputs --binary /new/runtime/numivivo-omics \
  --out /new/intervals
python Tools/Omics/PerturbationPrediction/Intervals/score.py \
  --inputs /path/to/cross-study/inputs --native /new/intervals \
  --previous-native /path/to/cross-study/native-complete --out /new/scores
```

The runner retains complete native outputs in 26 verified archives totaling
191,456,140 bytes rather than 524,424,073 uncompressed bytes. Only each finished
invocation's scratch is removed after all native checks, archive-member byte
verification and an open-handle check. All source inputs remain retained; failed
invocations leave their scratch for diagnosis. Restore one complete bundle with:

```sh
python Tools/Omics/PerturbationPrediction/Intervals/restore.py \
  --native /path/to/intervals/native --fold cross-HIRISA-to-Kang-patient_101 \
  --out /new/restored-fold
```

The output contains ordinary `model` and `prediction` directories for the native
verify commands. Restoration rejects mismatched archives and unsafe member paths.
Use the original qualified executable to replay historical receipts.

The [compact evidence archive](evidence/2026-09-11) retains exact metadata, plans,
freezes, scores, verification logs and source code. Full native archives remain
in `/Users/n/numivivo-donor-intervals-20260911/native/bundles` on `macmini` and
`/Users/home/numivivo-donor-intervals-20260911/native/bundles` locally. Their member
hashes and compressed identities are retained in the prediction freeze. The
compact repository archive alone does not contain the complete experimental data.
