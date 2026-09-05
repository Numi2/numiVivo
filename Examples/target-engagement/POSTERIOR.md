# Fit kinetic uncertainty, then predict held-out conditions

This is a worked software example, not a measured BTK/zanubrutinib dataset. Every parameter and observation in `synthetic-inference.json` is labelled assumed. The observation values were generated analytically, not by executing the new inference or forward runtime. No posterior result, Apple build, or new regression pass is claimed by this example's presence.

## Analytic fixture

The model is reversible binding with kon = 100000 M^-1 s^-1 and koff = 0.2 s^-1, zero covalent conversion, zero target turnover, and initially free target. A constant unbound exposure L is applied from t=0 through t=10 seconds, then removed. The exact occupied fraction is:

```text
a = kon * L
occupied(t <= 10) = a/(a+koff) * [1 - exp(-(a+koff)*t)]
occupied(t > 10)  = occupied(10) * exp(-koff*(t-10))
```

The two calibration exposures are 1e-6 M and 2e-6 M. The validation exposure is 1.5e-6 M; the test exposure is 0.5e-6 M. Conditions have separate leakage-group labels. Observations are noise-free analytic values at 1, 2, 5, 10 and 20 seconds. SD = 0.01 is an assumed Gaussian observation model, not a measured assay precision. Initial fractions and the exposure discontinuity are explicit.

The fitter estimates association and dissociation jointly, sharing each across the four same-context conditions. Independent log-uniform priors are declared over [30000,300000] M^-1 s^-1 and [0.03,1] s^-1. The base model retains the generating values for provenance, but each likelihood call overwrites the two fitted fields with its candidate values. Neither the base generating values nor held-out outcomes initialize or optimize the posterior particles.

## Run on the Apple development environment

From the repository root:

```sh
swift run numivivo engagement-fit \
  Examples/target-engagement/synthetic-inference.json \
  --checkpoint /tmp/numivivo-posterior-checkpoint.json \
  --output /tmp/numivivo-posterior.json \
  --store /tmp/numivivo-posterior-artifacts

swift run numivivo engagement-predict \
  /tmp/numivivo-posterior.json \
  --output /tmp/numivivo-posterior-predictions.json \
  --store /tmp/numivivo-posterior-artifacts

swift run numivivo engagement-sensitivity \
  /tmp/numivivo-posterior.json \
  --output /tmp/numivivo-posterior-sensitivity.json
```

Existing result files are rejected unless `--force` is supplied. A requested checkpoint file is updated only at complete tempering-stage boundaries. The artifact store receives immutable checkpoints and results; its normal identity and I/O rules apply.

Resume uses the same inference problem and sampler configuration:

```sh
swift run numivivo engagement-fit \
  Examples/target-engagement/synthetic-inference.json \
  --resume /tmp/numivivo-posterior-checkpoint.json \
  --checkpoint /tmp/numivivo-posterior-checkpoint.json \
  --output /tmp/numivivo-posterior-resumed.json
```

The matching resume/checkpoint path explicitly authorizes continuing that checkpoint file. Training, prior, sampler or numerical-policy changes invalidate the saved plan identity. A failed run retains its last completed stage; an incomplete population cannot be passed to the predictive command as a posterior. A checkpoint-store or file failure can occur after an in-memory stage has committed, so inspect the returned failure and available artifacts rather than assuming rollback.

The current posterior likelihood uses the native FP64 reference solver. `--backend metal` is intentionally not a fitting option: an approximate GPU likelihood must have separate numerical qualification and identity before it can replace this distribution. The forward `engagement-run --backend metal` workflow remains separate.

## Read the results

The fit result includes the problem, prior definitions, final joint particles, their log likelihoods, RNG checkpoint state, tempering stages, incremental-weight ESS, mutation acceptance and initial-ancestor diversity. Beta=1 means the declared tempering path completed. It does not demonstrate convergence or globally adequate posterior exploration.

The predictive report contains one result per case, preserving calibration/validation/test labels. Latent credible intervals express conditional parameter uncertainty. Measurement-predictive intervals additionally include the assumed SD. These are pointwise intervals. They do not include omitted mechanisms, unknown exposure error or patient variability. Missing measurement SD leaves predictive fields null. Any failed posterior particle suppresses the affected case's intervals rather than being dropped.

The sensitivity report evaluates local derivatives at the posterior mean in prior coordinates. It compares h and h/2 estimates, forms the local information matrix, and reports weak parameter combinations. Numerical rank is withheld when derivative checks fail. Full local rank is not a global identifiability or biological-validity certificate.

## Regression executable

```sh
swift run --package-path Examples/target-engagement posterior-checks
```

The supplied regression source checks correlated Gaussian sampling, stage-boundary resume, changed-prior rejection, failed/omitted evaluations, held-out-data independence, the analytic kinetic fixture, predictive intervals and local information. It was not run during this development increment. Preserve the actual toolchain, revision and logs when running it; no results from the older kinetic checks should be attributed to this new code.

Detailed contracts and limitations are in `Documentation/Design/KINETIC_POSTERIOR_INFERENCE.md` at the repository root.
