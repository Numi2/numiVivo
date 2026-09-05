# Metal-screened inference with an FP64 target

Status: committed source implementation. This increment has not been built, typechecked, executed on a GPU, benchmarked, or statistically qualified during authoring. The earlier 94-check portable manifest describes older exact source snapshots; it does not qualify the changed sampler or new GPU path. No new passed-test record is produced here.

## One statistical target

The existing tempered SMC fitter can use a deterministic Metal screen before its authoritative FP64 likelihood evaluation. This is delayed acceptance, not replacement of the statistical target by an approximate GPU likelihood. The priors, training likelihood, assay model, tempering weights and stored particle likelihoods remain those of the existing FP64 implementation.

Let l(u) be the authoritative log likelihood and s(u) the finite deterministic screen in uniform prior coordinates. At inverse temperature beta, a symmetric proposal u to u' uses two independent random acceptance draws:

```text
alpha1 = min(1, exp(beta * (s(u') - s(u))))
alpha2 = min(1, exp(beta * ((l(u') - l(u)) - (s(u') - s(u)))))
```

Only stage-one survivors require a new FP64 evaluation. The second ratio corrects the approximation. Initial particles, all incremental weights, temperature selection and accepted stored likelihoods use l. The existing fixed-covariance Gaussian random walk and independent uniform-prior proposal are symmetric in these coordinates. Out-of-prior proposals remain rejected.

A poor screen can reduce acceptance, worsen mixing or make execution slower. Correct acceptance ratios are not a claim of adequate finite-particle exploration, numerical qualification or biological accuracy. A screen must be an immutable function of one candidate and its context, independent of batch companions, lane position, ordinal and prior calls. No adaptive refitting is performed.

Method references: Christen and Fox (2005), *Markov chain Monte Carlo Using an Approximation*, JCGS 14(4):795-810, DOI 10.1198/106186005X76983; Del Moral, Doucet and Jasra (2006), DOI 10.1111/j.1467-9868.2006.00553.x. These are methodological references, not runtime dependencies.

## Fixed numerical plan

`VivoMetalTargetCasePlan` compiles one ordinary F1 ProgramPack and immutable schedule per calibration case. A conservative rate bound uses declared prior bounds and the full exposure schedule. It selects a fixed dyadic timestep independent of the current candidate batch. Choosing the step from the current batch maximum would make scores depend on companion candidates and violate the screen contract.

The scheduler represents exposure discontinuities and observation times explicitly and rejects unsupported FP32 clock resolution. Nonzero subnormal inputs are rejected rather than silently flushed to zero. Parameter bounds cover initialization and prior support with outward FP32 rounding; these are labelled numerical support, not confidence intervals. Step and operation counts are bounded per case and in aggregate. Overly stiff or long requests fail the declared execution budget rather than being labelled biologically impossible.

The supported graph has free, reversible, covalent and optional competitor target states. Arbitrary rules, stochastic species and temporal state are not part of this offline screening profile.

## Reused Metal execution

The implementation reuses `VivoMetalArena`, the target source compiler and the existing `nvivo_f1_heun_predict` and `nvivo_f1_heun_correct` kernels. There is one retained arena per calibration case, not per particle. No new reaction equations, circulation engine or authoritative biological-state owner are introduced.

Parameter-major/environment-strided uploads, external-input values, assay parameters, score accumulators and failure flags persist across proposal batches. Padded capacity duplicates lane zero for computation only; it is not additional statistical evidence. Bounded command chunks reuse pipelines and buffers. Full trajectories are not read back for likelihood evaluation.

The new shader adds only external-exposure updates and observation scoring. Each lane checks its target fractions and accumulates independent Gaussian log scores using compensated FP32 summation. The CPU receives scores/status and combines cases in FP64. Standard ProgramPack monitors and state validation remain active. No candidate is silently removed after a numerical failure.

The workspace actor reserves shared buffers across asynchronous execution. Cancellation releases them only after submitted GPU commands complete. Hardware command failures poison the workspace. Retained-buffer limits cover owned arenas and buffers, not all process, driver or pipeline allocations. No measured speedup is claimed.

## Assay and approximation scope

Version one screens exact-valued observations. It includes candidate-specific noise scaling, additional residual SD and bias. Missing reported uncertainty is not assigned zero: an explicit positive residual model must supply it.

Within-case assay correlation is omitted only by the approximate Metal score; the full correlated likelihood remains in FP64 acceptance correction. Probe reports compare with both a diagonal FP64 calculation and the authoritative likelihood, separating deliberate assay approximation from integration/precision differences. Censored-data GPU screening is rejected. The existing unscreened workflows retain their prior censoring support.

Hierarchical donor variation, general correlated censoring and full biological model discrepancy are not implemented by this increment.

## Identity and restart

`VivoPosteriorScreeningPolicy` binds an optional screen fingerprint and evaluation budget into the execution plan, separately from the authoritative likelihood fingerprint. Absent optional fields retain the legacy encoding shape. Screened checkpoints cannot resume with a missing or changed screen.

The Metal identity includes the executing binary SHA-256, shader-source hashes, OS, device registry/name, numerical settings, fixed case programs/schedules and workspace description. Running bundle resources are assumed immutable. Hashes establish identity, not trust in arbitrary code or proof of scientific correctness.

Completed-stage ledgers distinguish authoritative evaluations, screen evaluations, screened-out proposals, out-of-prior proposals and final acceptance. Resampled current-particle screen scores are recomputed at each stage; GPU scratch buffers are not checkpoint state. Exact evaluation ordinals remain sequential for actual FP64 evaluations.

Partial-stage failure or cancellation retains the last completed population and RNG state. Resume can repeat work attempted after that boundary; completed-stage budgets do not measure every failed process attempt. Probe/setup work is separately reported. Artifact publication and checkpoint-file updates are individual writes, not a cross-filesystem crash transaction.

`VivoTargetPosteriorRecord` validates the authoritative problem separately from screening metadata and then validates the actual screen-bound ledger. Existing FP64 prediction and experiment-design workflows can therefore consume screened results without changing their statistical target.

## Commands and checks

`engagement-screen-check` performs actual-device repeated-batch, permutation, single-candidate/padding and changed-ordinal probes. It compares scores bitwise across batch arrangements and reports FP64 differences. Finite probes cannot establish correctness across the whole prior.

`engagement-fit-screened` requires an artifact store, persists the screen description and probes, and refuses to start when determinism probes fail. There is no silent CPU fallback. The original `engagement-fit` remains the unscreened FP64 route. See `Examples/target-engagement/SCREENING.md` for invocations.

The new `screened-posterior-checks` executable covers imperfect, perfect and constant screens; exact-likelihood retention; separate evaluation ledgers; serialization/resume; changed-screen rejection; failure and budget behavior; and legacy optional fields. `--metal` adds actual GPU probes, a synthetic screened fit and existing prediction. These are supplied regression sources, not recorded passes.

## Source locations

Under `Sources/NumiVivoKit/Calibration/`, this increment extends `VivoPosteriorModel.swift`, `VivoTemperedPosteriorSampler.swift` and `VivoTargetPosteriorProblem.swift`, and adds `VivoMetalTargetScreenPlan.swift`, `VivoMetalTargetLikelihoodScreen.swift` and `VivoMetalTargetScreenChecks.swift`.

`Sources/NumiVivoShaders/Resources/NumiVivoTargetLikelihood.metal` is explicitly included by the root package manifest. `Sources/NumiVivoCLI/VivoScreenedInferenceCLICommands.swift` exposes the workflow through the existing command router and artifact I/O.

## Qualification handoff

Still unexecuted here: Swift typechecking, full package/ProgramPack compilation, shader compilation, actual-device probes, screened statistical regressions, end-to-end posterior comparisons and performance measurements. The prior 94-check result must not be cited as covering these changed sources.

Qualification should retain the exact revision and binary/device/OS identity, compare repeated screened and unscreened runs, refine the screen timestep, and measure total time and mixing as well as avoided FP64 calls. This increment does not establish clinical BTK predictions, cellular efficacy, tissue response or patient-specific treatment outcomes.
