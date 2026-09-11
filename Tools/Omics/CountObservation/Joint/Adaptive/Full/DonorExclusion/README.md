# Training-only calibration for donor-excluded joint fits

All **13 donor omissions** now have native, source-bound count-calibration
parameters for every source gene: eight Kang folds with 15,706 genes and five
HIRISA folds with 18,082 genes. This produces **432,116 gene/condition records**.
All **3,803,077** independent numerical comparisons pass; the maximum scaled
discrepancy is `2.87e-13`. These are inputs for subsequent joint likelihood fits,
not completed donor-held-out predictions or biological validation.

The [earlier paired-moment diagnosis](../../../../Paired/README.md) already
recomputed marginal calibration for each omission, but retained covariance
diagnostics rather than the dispersion parameters needed by the joint fitter.
This execution exposes the complete calibration models using the same qualified
original-cell moments. It does not reread or select a new cell or gene panel.

`prepare.py` creates one input per omission containing only the remaining donors'
complete moments. `Calibrate.swift` rejects an input containing its declared
excluded donor, pairs conditions by identity, invokes the native
`VivoCountObservationCalibration.fit` for each condition, and fingerprints the
exact training input bytes. The excluded donor's numeric observations are absent
from the input, rather than merely ignored by the optimizer. Reversing the
training group order reproduces every fitted model field exactly, apart from the
expected input/source fingerprint change. An input containing its excluded donor
is rejected before fitting.

`verify.py` binds each input to the retained original moment groups minus exactly
one donor. It independently reconstructs every dispersion numerator, denominator,
raw/admitted dispersion, boundary flag, latent/measurement variance and Gamma
prior parameter, including all unavailable statuses. Gamma-prior unavailability
alone does not make a valid cell dispersion unavailable to the joint model.

## Eligibility changes before joint fitting

All changes below refer to availability of **both** cell dispersions, rather than
the separate Gamma-prior admission rule. The full-cohort eligible counts are
8,414 Kang genes and 13,267 HIRISA genes.

| Origin | Donor omissions | Previously eligible genes becoming unavailable per omission | Previously unavailable genes becoming eligible per omission |
| --- | ---: | ---: | ---: |
| Kang | 8 | 10–670 | 9–51 |
| HIRISA | 5 | 109–532 | 165–295 |

These changes demonstrate why full-cohort dispersions must not be substituted
into a donor-excluded fit. They do not estimate parameter uncertainty or establish
interval coverage. Dispersion remains a plug-in moment estimate under the
existing cell/depth assumptions, with the original kernel-domain limits retained.
No omitted donor is replaced, and no unavailable gene is silently filled in.

## Reproduction and next boundary

The full HIRISA evidence archive includes this study under
`donor-exclusion/study`, its recipes under `donor-exclusion/recipes`, and both
original moment inputs. `protocol.json` binds all 13 inputs; `runtime-freeze.json`
binds the actual M4 Pro executable and the unchanged native library used by the
full-cohort fitter. The native calibration stage took 33.57 seconds across the
sequential SSH runs; this is not end-to-end joint-fitting performance.

Build the calibration harness with `bash DonorExclusion/build.sh SCOPED_BUILD`
after the standard scoped H5AD build. `prepare.py STUDY PAIRED_MOMENTS_STUDY`
creates the 13 input files; `run.py STUDY` uses the pinned executable named by
`runtime-freeze.json` over SSH to `macmini`.

To reconstruct the independent check after restoring that archive:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 \
  RESTORED/donor-exclusion/recipes/verify.py \
  RESTORED/donor-exclusion/study RESTORED/study
```

The next stage must fit a separate joint model from each fold's training donors
and these fold-specific dispersions, then produce control-only predictions before
scoring omitted outcomes. Existing Kang/HIRISA observations are development data;
this does not replace a frozen test in an independent biological context.
