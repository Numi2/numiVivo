# Conditional kinetic posterior inference

Status: source implementation. This increment was not built, typechecked, executed, statistically qualified, or validated against biological measurements during authoring. The existing 44-check kinetic reference record predates this increment and does not qualify the new posterior code. Apple builds and regression execution remain a user-side handoff.

## What is implemented

The new path fits context-qualified kinetic parameters to the existing target-engagement study representation, propagates joint posterior particles to calibration and held-out cases, and computes local likelihood-information diagnostics. It reuses the existing FP64 target-engagement forward solver, study schema, SplitMix64 generator, canonical artifact hashing and pinned-root document I/O. It does not replace `VivoAdaptiveEnsembleOptimizer`, relabel its elite archive as a posterior, or introduce another physiology engine.

The command surface is `engagement-fit`, `engagement-predict`, `engagement-sensitivity` and `posterior-help`. Fitting is a deterministic-FP64-likelihood workflow. The generic sampler accepts a batched evaluator so a qualified GPU forward backend can be integrated later with a distinct likelihood identity. The current CLI does not silently substitute approximate Metal or surrogate evaluations into an FP64 posterior.

## Probability model and prior semantics

`VivoPosteriorParameter` supports two independent bounded prior families:

- `uniformPhysical`: x = lower + u (upper - lower).
- `uniformLogPhysical`: x = exp(log(lower) + u [log(upper) - log(lower)]).

Here u is uniform in (0,1). These are different probability distributions, not interchangeable optimization transformations. The sampler operates in prior coordinates u, targeting L(x(u))^beta. The physical prior is defined by the pushforward, so an extra physical-space Jacobian is neither missing nor to be added again. Arbitrary marginal priors and correlated prior distributions are not supported by this version. Correlations induced by the likelihood are represented by the joint particle population.

Kinetic binding fields include association, dissociation, inactivation, target turnover, competitor association/dissociation/concentration, and a multiplicative exposure scale. Units must match the field. Sharing is explicit through case identifiers. Overlapping writers and shared rates across unequal chemical/biological contexts are rejected. Exposure-scale priors must remain inside each model's declared unbound-concentration domain.

Baseline target abundance is deliberately not an inferable field for this occupancy-only reservoir model. Normalized target-fraction dynamics do not depend on that abundance. Likewise, association rate and an unknown exposure multiplier can be confounded; the software must not claim that a posterior interval proves those mechanisms independently identified.

## Training and held-out separation

The prepared likelihood contains only cases labelled `calibration`, sorted by identifier, the bindings projected onto those cases, the explicit priors, the forward numerical policy and a likelihood-method identifier. Held-out observation values are not included in that identity or read by likelihood evaluation. The full result retains the original study, including validation and test cases, for later evaluation. Changing training data invalidates a checkpoint's plan identity; changing held-out observations alone does not alter the fitted likelihood.

The existing study validator rejects a declared leakage group split across calibration/validation/test roles. Authors still must choose defensible groups. The code cannot discover undisclosed shared donors, assay batches or other dependencies from labels alone.

The kinetic likelihood uses independent Gaussian measurement errors with known positive SD. Missing SD is rejected for fitting, not replaced with zero or a default. The full normalized Gaussian log density is accumulated with compensated summation. Assay correlation, uncertain noise scales, censoring, model-discrepancy terms and population random effects require additional explicit observation models. Their uncertainty is not included in current intervals.

Every candidate is applied to a copy of the original kinetic model. Original evidence and priors remain in the result. Temporary candidate rates are marked as proposals/assumptions, not measured or independently validated fitted values. No scalar parameter overwrites the original source experiment in place.

## Tempered sequential Monte Carlo

`VivoTemperedPosteriorSampler` implements the following bounded computation:

1. Draw particles independently in uniform prior coordinates and evaluate their deterministic log likelihoods.
2. Choose the next inverse temperature beta by bisection on the incremental-weight effective sample size. A temperature step below the declared minimum stops the run; the sampler never forces beta to one merely to report completion.
3. Form stable, normalized incremental likelihood weights.
4. Estimate a regularized full covariance in prior coordinates and freeze it for the upcoming mutation stage.
5. Systematically resample, retaining initial-ancestor identifiers.
6. Apply a fixed number of Metropolis sweeps. Proposals mix a correlated Gaussian random walk with independent draws from the prior coordinates. Out-of-prior random walks are rejected; coordinate-wise reflection is not used because it would require different treatment for correlated proposals.
7. Commit the complete population, stage diagnostics, RNG state and evaluation ordinal before advancing.

Completed stages have equal particle weights because each contains resampling followed by invariant mutation. Metropolis acceptance is determined by the tempered log-likelihood difference in the uniform prior coordinates. The proposal kernels are fixed during a stage. Adaptive SMC is a finite-particle approximation, not an exact posterior certificate.

The pre-resampling ESS measures weight concentration. It is not an independent-posterior-sample count, an MCMC chain ESS, or proof of convergence. Ancestor diversity and mutation acceptance are recorded separately. Beta=1 indicates completion of the declared annealing path, not adequate exploration of every posterior mode. Repeated seeds, larger populations, additional mutation sweeps and scientific validation remain required qualification work. No marginal-likelihood estimate is emitted.

## Failure, cancellation and restart

The current sampler contract accepts finite deterministic log likelihoods. It does not implement zero-probability support through arbitrary numeric penalties, stochastic likelihood estimators or pseudo-marginal state. A batch must return each exact candidate once; returned coordinates and physical values are checked against the request even when completion order differs.

A numerical or evaluator-contract failure stops inference and retains the last complete stage. Failed candidates are recorded; they are not removed to condition the posterior on successful simulations. A resource limit also leaves the previous committed checkpoint intact. The output explicitly distinguishes a failed/tempered population from a completed posterior, and downstream prediction rejects incomplete results.

An actor reservation prevents concurrent calls to `run` from sharing scratch state across suspension. `checkpoint()` can observe only the last committed population. Stage work uses local RNG and particle copies. A cancelled or failed partial stage does not replace accepted state. Checkpoints bind training likelihood, priors, sampler configuration and forward numerical policy. Resume reproduces the committed random stream and candidate ordering under the same implementation and numerical environment; this behavior is covered by supplied but unexecuted regressions.

A progress callback runs after stage commit. A failing checkpoint file write does not imply rollback of evaluator side effects or of the in-memory committed stage. Immutable-store publication and an ordinary checkpoint file are separate writes, not a cross-filesystem crash transaction. The CLI makes this state visible through its failure result. Hashes establish content identity, not authorship or proof of numerical correctness.

## Posterior prediction

The predictor evaluates every joint particle, preserving rate/exposure correlations. It never reconstructs a population from independent marginal intervals. For every declared study observation it reports the latent mean, SD and pointwise equal-tail credible interval. Where assay SD is known, it also computes quantiles of the Gaussian mixture over particle predictions by deterministic CDF bisection, plus log predictive density and a pointwise coverage indicator.

Measurement-predictive intervals include the specified assay noise; latent credible intervals do not. Unknown SD leaves predictive fields null. Correlation entries for collapsed prior coordinates are also null rather than invented zero correlations. Tiny-scale parameter SD is computed with scaled deviations to avoid squaring small physical values into zero.

A failed particle suppresses the affected case's interval report rather than renormalizing surviving particles. Requested cases and their failures remain visible. Calibration cases retain an in-sample label. Synthetic observations retain their origin. Intervals are conditional and pointwise; they are not simultaneous confidence bands, patient-specific uncertainty guarantees or clinical safety claims.

## Local information diagnostics

The sensitivity command evaluates the calibrated observation model at the posterior mean in prior coordinates. It estimates each derivative at h and h/2, whitens by the declared measurement SD, forms J'J, and computes its symmetric eigensystem using small FP64 algebra. Weak eigen-directions identify locally weak combinations of parameters on the chosen prior-coordinate scale.

Numerical rank is withheld if finite-difference consistency fails. A full local rank does not establish global or structural identifiability, posterior convergence, causal mechanism correctness, or transferable biological prediction. Bounds can require one-sided finite differences. Changing priors changes the coordinate scale of this diagnostic.

## Apple-native execution boundary

All new inference and forward-likelihood implementation is native Swift with FP64 numerical/statistical state. The target adapter evaluates bounded concurrent candidate tasks; completion order is canonicalized before sampler mutation. The small dense covariance and information matrices have at most 32 dimensions. Native Metal target-engagement execution remains available through the existing forward workflow, but the new posterior CLI uses the declared reference likelihood until a distinct GPU likelihood has numerical error qualification. There is no CUDA/Python dependency or unsupported Apple-GPU FP64 claim.

The existing Metal arena has parameter-major, environment-strided storage and configuration-bound checkpoints. Those are the intended integration points for future parameter-varying GPU likelihood cohorts. Constructing one independent GPU runtime per posterior particle is not the intended production design. A GPU backend must bind precision, integration policy, device-sensitive semantics and failure behavior into its likelihood identity; CPU and GPU results must not be mixed under one unchanged posterior checkpoint.

## Files and commands

New Kit files are under `Sources/NumiVivoKit/Calibration/`:

- `VivoPosteriorModel.swift`: prior, plan, candidates, particles, stage ledger and checkpoint contracts.
- `VivoPosteriorNumerics.swift`: weighting, temperature selection, covariance, resampling and predictive-mixture numerics.
- `VivoTemperedPosteriorSampler.swift`: actor-owned batched SMC execution.
- `VivoTargetPosteriorProblem.swift`: typed kinetic binding, calibration-only likelihood and fitting integration.
- `VivoTargetPosteriorPrediction.swift`: joint prediction and pointwise interval evaluation.
- `VivoTargetPosteriorSensitivity.swift`: local information and finite-difference diagnostics.

The CLI router is `Sources/NumiVivoCLI/VivoPosteriorCLICommands.swift`. An analytic, explicitly synthetic four-condition example is `Examples/target-engagement/synthetic-inference.json`. Its workflow is documented in `Examples/target-engagement/POSTERIOR.md`. The `posterior-checks` executable covers correlated Gaussian sampling, restart identity and deterministic continuation, numerical failure retention, omitted-candidate rejection, training/held-out separation, the synthetic kinetic fit, predictive intervals and local sensitivities.

No new test-results log or passed-check count was recorded in this increment. Execute and retain results on the actual source revision before treating these paths as qualified.

## Primary methodological references

- Del Moral, Doucet and Jasra (2006), Sequential Monte Carlo Samplers, JRSS B 68:411-436. https://doi.org/10.1111/j.1467-9868.2006.00553.x
- Stan User's Guide, reparameterization and change of variables. https://mc-stan.org/docs/stan-users-guide/reparameterization.html
- PyMC's SMC documentation describes the same prior-to-posterior tempered distribution family; it is a methodological reference, not a runtime dependency. https://www.pymc.io/projects/docs/en/v5.17.0/api/generated/pymc.smc.sample_smc.html
