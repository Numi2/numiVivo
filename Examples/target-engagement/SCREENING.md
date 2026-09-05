# Metal-screened kinetic inference

This workflow uses a reusable Metal approximation to screen proposed kinetic parameters, then applies the existing FP64 likelihood correction before accepting survivors. It does not substitute an approximate GPU posterior for the declared statistical model.

**Execution status:** the new screen, changed sampler and example regressions have not been built, typechecked or run during this increment. No GPU result or speedup is asserted below. The existing synthetic inputs are software fixtures, not measured BTK or zanubrutinib data. Earlier portable check results do not cover these changed sources.

## Check the screen on Apple hardware

From the repository root, with a supported Apple SDK/toolchain:

```sh
swift run numivivo engagement-screen-check \
  Examples/target-engagement/synthetic-inference.json \
  --screen-policy Examples/target-engagement/screen-policy.json \
  --output /tmp/numivivo-screen-check.json \
  --store /tmp/numivivo-screen-artifacts
```

The command evaluates repeated batches, reversed batches, each candidate alone, and renamed candidate ordinals. A candidate's score should be unchanged by those arrangements. It also reports the difference from a diagonal FP64 score and from the authoritative FP64 likelihood. These are finite probes, not universal validation or performance measurements.

## Fit with FP64-corrected GPU screening

```sh
swift run numivivo engagement-fit-screened \
  Examples/target-engagement/synthetic-inference.json \
  --screen-policy Examples/target-engagement/screen-policy.json \
  --store /tmp/numivivo-screen-artifacts \
  --checkpoint /tmp/numivivo-screen-checkpoint.json \
  --output /tmp/numivivo-screened-posterior.json
```

The artifact store is required for fitting so the screen implementation description and actual-device probe report are retained. The fitter reruns the probes and refuses to start if batch-independence checks fail. It does not silently fall back to an unscreened fit.

All initial particle likelihoods, tempering weights and accepted stored likelihoods remain FP64. A rejected screen proposal avoids a new FP64 evaluation; a screen survivor still needs the exact acceptance correction. The result separately records screen evaluations, screened-out proposals and FP64 evaluations. Fewer exact calls do not necessarily mean lower total runtime or better exploration.

The default screen uses one retained arena per calibration case, with shared buffers for the whole candidate cohort. Timestep selection uses the prior bounds, not the current batch. The supplied policy limits retained buffers to 256 MiB and explicitly bounds steps, command chunking and screening evaluations. It is not a whole-process memory cap.

## Resume and predict

```sh
swift run numivivo engagement-fit-screened \
  Examples/target-engagement/synthetic-inference.json \
  --screen-policy Examples/target-engagement/screen-policy.json \
  --store /tmp/numivivo-screen-artifacts \
  --resume /tmp/numivivo-screen-checkpoint.json \
  --checkpoint /tmp/numivivo-screen-checkpoint.json \
  --output /tmp/numivivo-screened-resumed.json

swift run numivivo engagement-predict \
  /tmp/numivivo-screened-posterior.json \
  --output /tmp/numivivo-screened-predictions.json
```

Existing result files require explicit `--force` to replace. Using the same resume/checkpoint file authorizes continuing that saved checkpoint. Resume checks the training problem, priors, numerical policy and screen identity. A changed executable, device description, shaders or screen configuration invalidates that screen-bound restart. The same posterior result remains compatible with the existing FP64 prediction and experiment-design commands because its target likelihood has not changed.

A cancelled or failed partial stage does not replace the last completed checkpoint. Work after the checkpoint may be repeated during resume. Artifact publication and checkpoint-file replacement are separate writes, not one cross-filesystem transaction.

## Current scope

The screen supports exact-valued Gaussian observations, including inferred noise scale, residual SD and bias. When the authoritative assay model is correlated, the screen uses diagonal marginals and the second-stage FP64 correction retains the full covariance. Censored-data screening is rejected. The normal `engagement-fit` workflow remains available for the existing reference-supported assays.

Unsupported FP32 clocks, nonzero subnormal parameters, excessive fixed-step budgets and numerical failures are errors, not zero likelihoods. Failed candidates are not silently omitted from the inferred distribution. The screen does not add donor hierarchy, finite-drug physiology coupling or a cellular efficacy model.

## Regression executable

```sh
swift run --package-path Examples/target-engagement screened-posterior-checks
swift run --package-path Examples/target-engagement screened-posterior-checks --metal
```

The first invocation still builds the Apple package but tests delayed-acceptance logic with deterministic Swift screening functions. It checks imperfect, perfect and constant approximations, exact stored likelihoods, ledgers, checkpoint continuation, missing/changed screens, failures, budgets and legacy optional fields. The second additionally runs actual GPU probes, a synthetic screened fit and the existing predictive workflow.

Neither invocation was executed during this development increment. Retain the actual revision, toolchain, binary/device identity and logs when running them. Qualification also requires repeated screened/unscreened posterior comparisons and timing/mixing measurements rather than only a successful executable exit.

See `Documentation/Design/METAL_SCREENED_INFERENCE.md` for the acceptance equations, implementation contracts and remaining qualification requirements.
