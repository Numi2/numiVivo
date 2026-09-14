# Cross-study IFN-beta NB2 follow-up

**Status: native execution and independent scoring PASS; biological transfer
comparison FAIL in both directions.** This is an opt-in negative-binomial
follow-up to the historical log-linear Kang-HIRISA experiment. It asks whether
the count-based response model improves held-out donor transfer for the same
known IFN-beta perturbation. It does not test an unseen perturbation.

The frozen inputs contain eight paired Kang donors, five paired HIRISA donors,
and 11,884 exact, unique gene-symbol matches. Each source library keeps all of
its measured features in its normalization denominator; the shared panel is
only the output and model-selection space. Cross-study folds train on every
donor in the other study. Within-study folds leave the query donor out and are
reported as matched references. Held-out treated rows are read only by the
independent scorer after model fitting.

## Cross-study result

Lower response RMSE is better. Values are equally weighted donor means over the
complete shared panel in natural log1p(CPM) units. NB2 is the native
`negativeBinomialEffect` estimate. The fixed transfer comparison requires an
estimate to beat both no-change and the cross-study training mean.

| Training → query | Donors | No change | Cross mean | Cross ridge | NB2 | NB2 vs no change | NB2 better no change / mean | Decision |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| HIRISA → Kang | 8 | 1.220670 | 1.135572 | 1.178674 | 1.184337 | +2.98% | 8 / 0 | FAIL: worse than cross mean |
| Kang → HIRISA | 5 | 0.341591 | 0.521068 | 0.596829 | 0.343493 | −0.56% | 2 / 5 | FAIL: worse than no change |

The cross-study NB2 estimate therefore does not qualify reliable context
transfer. HIRISA-trained NB2 helps the Kang no-change baseline for all eight
donors but remains worse than the simple HIRISA training mean. Kang-trained NB2
is slightly worse than no change for HIRISA even though it beats the cross-study
mean for all five donors. These are reused, previously inspected studies with
health status, exposure duration, preparation and assay chemistry changing
together.

The matched within-study references are useful execution controls, not new
transfer evidence:

| Query study / mode | No change | Mean | Ridge | NB2 | NB2 better no change / mean |
| --- | ---: | ---: | ---: | ---: | ---: |
| Kang / within | 1.220670 | 1.156071 | 1.132625 | 1.187743 | 8 / 3 |
| HIRISA / within | 0.341591 | 0.106039 | 0.107549 | 0.149177 | 5 / 0 |

## Evidence and limits

All 26 folds completed 105 native commands: fit, model verification, prediction
and prediction verification for each fold, plus one exact repeated prediction.
Each report contains five estimates, giving 130 frozen prediction vectors. The
independent pure-Python scorer reconstructs every source library from the
frozen NPZ counts, verifies all native training sparse rows and model moments,
checks the prediction transformations, and scores all 26 held-out treated rows.
The largest reconstructed model difference is `4.30e-13`.

The opt-in NB2 fits use paired-donor count effects with per-feature dispersion,
standard-error and status diagnostics. Cross-study models identify 11,675
features for HIRISA → Kang and 4,696 for Kang → HIRISA; unavailable features
remain explicit no-change fallbacks. No batches are adjusted: HIRISA batches are
donor-nested and Kang reports one uninformative `unreported` batch. This is a
model identifiability decision, not evidence that the assays share a batch.

The native executable was built from revision `8e534a11d6e85f0b27980d44db348ae62c96f9e3`
on the physical M4 Pro Mac mini. The execution binary SHA-256 is
`d69a9df8978860e345008827bf53a1a4ce0109b08a1f337f757a38e1a6807895`.
The input freeze is `3fde88150541b518f6bbd0e69c87f6d7fa35bafb7877b23703f9163cbbce8c5f`,
the native prediction freeze is
`3eb80ac012a19fb53546faf835d2971bece4b15c098e501fe107ce6c7685d39f`,
the runner is `87bcecda728eff6e68d7d29c1136b38b3f109af341b690d3e18dd0d163659e63`,
and the independent scorer is
`e608903b9fcb7bb5514c9818ccd47480dcca01ed3a70e3e157dae5c224c7cc73`.
The scored checks, summary and per-fold rows are respectively
`0c6f6781e25707386de63fdbd83788f2a70d0f1fd0a6aad25f6a6e072cf17fcd`,
`bb228b3137a74a41fa634859bae0bacd975bd654d866730e05d83b67f2dbe1e2` and
`a8fcf5e6242765b345e2b574c93ee453090ffb82e679a76537c08ed44a7ef9f3`.
Compact score, freeze and command evidence is retained in
[evidence/2026-09-14-nb2](evidence/2026-09-14-nb2).

This result supports a bounded conditional molecular RNA estimate for a known
IFN-beta response. It does not establish causal mechanism, calibrated
uncertainty, unseen-study generalization, single-cell response distributions,
phenotype, disease progression, treatment benefit or any other general
biological outcome.

## Reproduce

Build the native CLI from the repository revision, then run the frozen folds
sequentially and score them after prediction is complete:

```sh
bash Tools/Omics/H5AD/build.sh /new/runtime --with-cli
NUMIVIVO_HDF5_LIBRARY=/path/to/libhdf5.dylib python Tools/Omics/PerturbationPrediction/CrossStudyIFNB/run_nb2.py --inputs /frozen/inputs --binary /new/runtime/numivivo-omics --out /new/native
python Tools/Omics/PerturbationPrediction/CrossStudyIFNB/score_nb2.py --inputs /frozen/inputs --native /new/native --out /new/score
```

The scorer uses only the Python standard library. The large H5AD, source-count
and per-fold model payloads remain on the Mac mini under
`/Users/n/numivivo-cross-study-ifnb-20260911`; the compact evidence directory
contains the hashes and summary needed to identify that external payload.
