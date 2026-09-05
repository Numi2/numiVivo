# Context-qualified kinetics and target engagement

Status: implemented source with a targeted portable FP64 regression run. Full Apple package compilation, native ProgramPack integration execution, Metal/Core ML behavior, and experimental predictive validation were not performed in this development pass. This document describes this increment, not the validation status of concurrent molecular-engine work elsewhere in the repository.

## Executable data flow

`VivoTargetEngagementExperiment` binds `VivoCovalentKineticPack`, a `VivoUnboundExposureTrace`, initial target fractions, and requested sampling times. It can be evaluated with the portable reference or lowered to the existing VivoProgram language and compiled into an ordinary F1 ProgramPack. The Metal cohort runner uses `VivoTransactionalMolecularRuntime`; it does not install a second reaction engine.

`VivoPhysiologyExposureAdapter` projects actual physiology snapshots into a trace with explicit units and concentration meaning. It returns the exact evidence bytes and fingerprint for the shared artifact store. `VivoTargetEngagementCompiler.physiologyBridge` instead prepares a live one-to-one link for the existing molecular/physiology coordinator. The live link requires an analyte that already represents unbound drug; a total concentration is not silently treated as unbound.

`VivoTransitionStateDerivation` supplies a narrowly scoped chemistry-to-kinetics boundary: a qualified activation free energy can become a conditional unimolecular conversion rate. This component does not generate the free energy, run QM/MM, estimate association kinetics, or establish that a chosen reaction coordinate is adequate.

## Model equations and assumptions

Let F, R, C, and Q be fractions of baseline target abundance that are free, reversibly drug-bound, covalently drug-bound, and competitor-bound. Let L(t) be prescribed free-drug concentration in M. The optional competitor concentration X is constant and unbound. Define a = kon L, c = konX X, b = koff, q = koffX, k = kinact, and d = target turnover. Time is in seconds.

```text
dF/dt = -(a + c + d) F + b R + q Q + d
dR/dt = a F - (b + k + d) R
dC/dt = k R - d C
dQ/dt = c F - (q + d) Q
```

The constant source d corresponds to synthesis d times baseline target concentration. Thus d(F+R+C+Q)/dt = d[1-(F+R+C+Q)]. Initial fractions must sum to one within the numerical contract. Drug occupancy is (R+C)/(F+R+C+Q); covalent occupancy is C divided by that total.

This is a conditional reservoir model. Ligand depletion, state-dependent turnover, saturable target synthesis, multiple competitors, metabolism, transport, and cellular efficacy are not inferred by these equations. They require their own represented mechanisms. The existing physiology runtime remains the exposure authority when used; the adapter does not reproduce its dynamics.

Context equality includes compound, target, variant, site, chemical state, host context, temperature, pH, and ionic strength. Version one deliberately rejects undeclared context transfer rather than extrapolating rates. Identity strings must be meaningful, but equality and hashes alone do not validate their scientific correctness.

## Numerical implementations

The FP64 reference augments the four-state vector with a constant coordinate. It evaluates a nonnegative uniformization series with scaling and squaring for each constant-exposure interval. The small-interval series uses 32 terms, while squaring and propagation have explicit work bounds. Samples split intervals exactly; exposure discontinuities do not reset target state. Positivity, finite values, and target balance are checked without clipping or renormalization. Initial data are also checked against the requested numerical tolerance, including a time-zero-only request.

The reported balance residual is not an error bound on all observables. Analytic regression comparisons establish selected numerical behavior, not universal accuracy or biological validity. A small-system FP64 run is appropriate even on Apple hardware; forcing every operation onto a GPU is not a requirement.

The Metal implementation lowers target fractions as dimensionless occupancy species and keeps free drug in molar units. It batches up to 4096 common-kinetics experiments over the existing parameter/topology tables. The actual immutable exposure schedule bounds the association rate used to choose a dyadic RK2 step. Output samples and exposure tables have bounded aggregate capacities. Exact FP32-supported time boundaries are required; unsupported clocks fail rather than shifting an intervention. There is no silent reference fallback.

The chosen RK2 step is a conservative stability/positivity restriction, not a demonstrated global-error tolerance. Apple qualification must compare against the reference and perform timestep refinement. Candidate publications are checked before release. Individual trajectories are not silently removed from an ensemble following failure.

## Evidence and uncertainty

Every kinetic parameter has a canonical unit, origin, uncertainty declaration, and source locator. Measured, fitted, and calculated inputs require an immutable evidence digest; assumed inputs may be explicit without a digest. A digest must refer to persisted evidence, not be invented to satisfy the validator. The current validator checks the digest representation; source availability, authenticity, and scientific applicability are separate responsibilities.

`unknown` is not zero variance. A bounded parameter range is a sensitivity assumption, not a confidence interval. A declared log-normal marginal is not a fitted joint posterior. The deterministic reference and Metal outputs are nominal trajectories, not predictive intervals.

Transition-state conversion uses:

```text
k = kappa * (kB T / h) * exp(-activationGibbsFreeEnergy / (R T))
```

It requires a pre-reactive bound-complex reference, a nonnegative activation Gibbs free energy, and a classical transmission probability in (0,1]. It rejects electronic-energy-only input, separated-reactant reference states, and unrepresentable rates. Barrier standard deviation can be propagated to conditional log-rate standard deviation, but this omits uncertainty in mechanism, sampling adequacy, transmission, and context transfer. The downstream parameter therefore retains unknown total uncertainty. Assumptions in a derivation remain visible in the derived parameter's origin.

The physical relation and standard-state distinction follow IUPAC Gold Book definitions T06470 and G02631:

- https://goldbook.iupac.org/terms/view/T06470
- https://goldbook.iupac.org/terms/view/G02631

## Held-out studies

`VivoTargetEngagementStudyEvaluator` evaluates fixed model parameters. It supports calibration, validation, and test partitions and rejects a declared leakage group split across roles. It does not use held-out observations to adjust parameters. The dataset author remains responsible for grouping shared donors, conditions, assay batches, and other dependence correctly.

Observation types are drug occupancy, covalent occupancy, free-target fraction, and competitor occupancy. Values are fractions rather than percentages; noisy observations are not clipped to [0,1]. Observation times must explicitly occur in the experiment. Standardized residuals use supplied measurement noise only and must not be labelled predictive-interval coverage.

All requested cases remain in the report. A numerical failure suppresses aggregate RMSE/MAE for the affected partition instead of improving a score by discarding failed cases. Measured observation counts are reported separately; synthetic fixture agreement is not experimental validation.

## Source integration

| Location | Responsibility |
|---|---|
| `Kinetics/VivoCovalentKinetics.swift` | Context, parameter evidence, exposure, target fractions, experiment contracts |
| `Kinetics/VivoTargetEngagementProgramSource.swift` | Ordinary typed F1 VivoProgram generation |
| `Kinetics/VivoKineticsDocumentIO.swift` | Public bounded I/O through existing pinned-root file implementation |
| `Pharmacology/VivoTargetEngagementReference.swift` | Portable FP64 nominal propagation |
| `Pharmacology/VivoTargetEngagementIntegration.swift` | Native compilation, physiology adapters, fingerprinted reference run records |
| `Pharmacology/VivoTargetEngagementMetalRunner.swift` | Existing-runtime Metal exposure cohorts |
| `Chemistry/VivoTransitionStateRate.swift` | Typed conditional free-energy conversion |
| `Chemistry/VivoTransitionStateDerivation.swift` | Persistent derivation evidence and explicit kinact replacement |
| `Calibration/VivoTargetEngagementStudy.swift` | Held-out occupancy evaluation |
| `NumiVivoCLI/VivoTargetEngagementCLICommands.swift` | Integrated command surface |

Kit-relative paths above are under `Sources/NumiVivoKit/`; the CLI is under `Sources/`. No independent artifact store, campaign scheduler, or biological reaction engine was introduced.

## Compatibility changes

Core ML surrogate construction now requires exactly one uncertainty source. Missing uncertainty no longer produces a zero array. Outputs must have the declared batch/output shape. Fixed uncertainty still needs appropriate calibration evidence; supplying a vector does not establish its quality.

`VivoSurrogateAuthorityGate` now validates its contract fingerprint, serializes predictions across backend suspension, requires initial authoritative evaluation, and retains pending refresh until a matching generation/context/input-bound request is acknowledged. The old no-argument `recordAuthoritativeEvaluation()` is unavailable. Callers must supply the request, authoritative result fingerprint, bounded output values, and successful numerical status. This is a trusted participant protocol, not proof that an untrusted caller actually performed the computation. Checkpointing/reconstructing a gate cannot grant authority: a new gate starts requiring refresh.

Legacy ProgramPack F2 requires `maximumSubsteps=1`. Adaptive redrawing after a realized rejection is not accepted as a statistically qualified method. This restriction prevents the internal multi-attempt path; caller-driven retries and exact/hybrid stochastic qualification still require an explicit policy and statistical checks.

The calibration compiler now derives parameter indices from actual metadata table order, rather than accessing a nonexistent `ParameterMetadata.index`. Strategy validation also rejects integer overflow when computing its parameter-count-dependent minimum population.

## Verification and remaining work

Executed: the 44-assertion portable harness with Swift 6.2.1 on x86_64 Linux, optimization enabled and warnings treated as errors. It covers selected analytic kinetics, turnover, competition, pulse/washout, stiffness, sampling-boundary invariance, invalid inputs, generated JSON structure, conditional rate conversion, held-out grouping, and failure-preserving metrics. Exact source identities are in `Tools/Kinetics/portable-results.json`.

Not executed here: the full Swift/C++ package, generated ProgramPack compilation, Metal shader execution, Core ML integration, filesystem integration through the full package, or the supplied Apple integration executable. Syntax parsing of new integration files is not typechecking or execution.

Remaining work for the larger biological objective includes data-backed target models and exposure calibration, uncertainty inference and identifiability, validated chemistry-derived kinetics, finite drug mass balance where needed, independently checked selectivity panels, and a calibrated cellular-response endpoint. This increment does not claim tissue/patient prediction or drug safety. Its immediate result is an executable, provenance-aware target-engagement workflow that those later capabilities can use.
