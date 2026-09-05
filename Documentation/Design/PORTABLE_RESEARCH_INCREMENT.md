# Portable research increment: measured-data boundaries and executed inference

This increment develops the portable portions of the next NumiVivo revision: posterior execution checks, richer assay likelihoods, a finite-drug reaction operator, single-measurement experimental design, and a published aggregate-data comparison boundary. It does not claim completion of GPU inference, hierarchical donor modelling, whole-body pharmacology, chemistry-derived protein rates, or a calibrated cellular endpoint.

## Execution status

The selected native Swift numerical implementation was compiled and exercised on Linux x86_64 with Swift 6.2.1. The retained run passed **94 assertions**: 72 inference/assay checks and 22 domain checks. Production numerical sources were optimized with whole-module optimization; the small test drivers were compiled separately without optimization. Both used Swift 6 mode with warnings as errors. Final test-driver updates reused the already compiled numerical library after verifying that all 15 production source files, the RNG excerpt and support adapter were unchanged. This is stated in the retained manifest, rather than represented as a fresh full-tree build.

The harness uses Foundation and a test-only artifact-identity adapter. SHA-256 uses system OpenSSL on Linux. Production CryptoKit, ProgramPack compilation, artifact storage, full CLI, Apple SDK integration and Metal execution were not validated by this harness. No clinical prediction was produced. `Tools/Posterior/Results/manifest.json` binds the exact numerical source and report bytes; a later changed source is not qualified by these records.

The previous posterior document's original 'not executed during authoring' statement is historical. This increment supersedes it only for the source snapshots and operations in the new portable manifest. It does not convert the older Apple integration example or the whole package into a validated release.

## Assay-aware inference and prediction

`VivoGaussianAssay.swift` introduces exact, left-censored, right-censored and interval-censored Gaussian observations. Intervals describe individual observations, not cohort ranges. Stable log-CDF and interval calculations avoid probability floors in extreme tails and loss of significance for narrow intervals. Tests compare with mpmath at 100 decimal digits and independent NumPy covariance calculations.

The optional `VivoTargetAssayModel` adds noise scaling, additional residual SD, bias and exponential within-case temporal correlation. The posterior binder can estimate those quantities jointly with kinetics. Missing reported SD can be used only when an explicit positive residual-SD model supplies the uncertainty. That residual SD is an assumption or inferred quantity, not a fabricated laboratory measurement. Parameter-dependent covariance contributes both its quadratic residual term and log determinant to the likelihood.

Training-only identities include the relevant assay model. Held-out values do not enter fitting. Legacy independent-Gaussian requests retain their v1 likelihood identity and original scalar accumulation order when no training assay extension is active. The extended path receives its own method identity. It does not reuse a saved posterior under changed likelihood semantics.

The predictor propagates complete joint kinetic/assay particles. Its schema-2 result separates latent occupancy intervals, measurement-predictive intervals, exact-observation densities and censoring-event probabilities. It includes the mixture of the full case likelihood for correlated observations; sums of marginal scores are not labelled as a joint score. Unknown measurement uncertainty and undefined correlations remain null. Censored placeholders never receive an exact-observation coverage test.

Correlated censoring is explicitly unsupported. An independent product of marginal censoring probabilities would not be the correlated likelihood. The old mean-only J'J sensitivity path likewise rejects extended training assays rather than omit information carried by covariance or censoring. Hierarchical donor effects and a general structural model-discrepancy process remain future work; an additive residual term is not a claim that all discrepancy has been represented.

## What the inference checks establish

Repeated-seed checks exercise a correlated Gaussian posterior, and stage-boundary serialization/resume reproduces the same particle and RNG history in this environment. Failed and omitted evaluator results stop inference rather than removing problematic parameter regions. End-to-end synthetic fits exercise unknown reported SD with inferred residual noise, correlated errors, censored prediction and held-out separation.

Twenty-four prior-predictive trials compare with an analytic truncated-normal posterior. The nominal 90% intervals covered the generating parameter in 21 of 24 trials. Ranks and coverage are retained as diagnostics, not promoted to a general simulation-based calibration certificate. SMC particles can be dependent; small fixed-seed experiments cannot establish universal posterior exploration or uncertainty coverage.

For one synthetic kinetic/noise example, posterior means were 98306.53226280556 M^-1 s^-1, 0.19259433341375246 s^-1 and residual SD 0.006043685735284298 in fraction units. The generating kinetic values were 100000 and 0.2, with fixed synthetic residuals. These are not BTK measurements or clinical parameter estimates.

## Finite-drug reaction operator

`VivoFiniteDrugReactions` operates on a local fixed-volume state containing free drug, free target, reversible complex, covalent complex, inactive metabolite, eliminated drug equivalents, and target-synthesis/removal ledgers. Association consumes one free drug and one free target; dissociation restores both. Bound drug removed with target turnover becomes the explicitly declared inactive metabolite. Free-drug clearance/metabolism and metabolite clearance preserve the total drug-equivalent ledger.

Each suboperator uses an exact positive reaction update. Symmetric composition with deterministic step doubling controls local splitting error. It does not use stochastic retries, concentration clipping or mass repair. Drug and target material balances are checked separately from the local-error estimator. The caller owns accepted state, dosing, transport and coupled commit/rollback. This is an extension point for the existing physiology coordinator, not an alternate circulation simulator or a completed PBPK implementation.

Checks cover equal-reactant association, finite ligand limitation, synthesis/turnover, zero-duration identity, work-budget failure and comparison against independently written eight-state ODE equations solved with SciPy Radau and DOP853. For the retained fixture, the native end-state maximum absolute difference from Radau was 5.931103484153953e-14 M. The two reference integrators agreed within 1e-16 M. These are fixture-specific results, not universal solver error bounds. Active-metabolite binding, state-specific turnover and spatial transport are not introduced by this operator.

## Experiment selection

`VivoTargetExperimentDesigner` ranks proposed single Gaussian measurements using expected information gain conditional on the finite joint posterior. A candidate supplies no observed outcome. Its template identifies the context and parameter-sharing rules explicitly, and new exposure/time inputs remain inside the declared applicability domain. The computation includes inferred assay uncertainty and reports conditional Monte Carlo standard error, two half-sample estimates and information per declared relative cost.

Tests check zero information for indistinguishable components, the log(2) limit for well-separated two-component predictions, decreasing information under high noise, and invariance to candidate ordering. The objective concerns the complete fitted parameter vector, including nuisance noise parameters. It is not automatically information about a single biological endpoint, a globally optimal multimeasurement plan or a treatment recommendation. No laboratory or calendar action occurs.

## Published-data benchmark boundary

`Examples/btk-aggregate-benchmark/observations.json` transcribes published zanubrutinib nodal occupancy summaries from Tam et al., Blood 2019, DOI 10.1182/blood.2019001160. It preserves medians, ranges, cohort sizes, visit labels and rounded percentages as aggregate facts. No patient-level observations, measurement SDs, exposure histories or kinetic parameters are invented from those summaries. Direct source-page access presented a challenge; the transcription used indexed primary-article text and identifies that retrieval boundary.

`VivoAggregateOccupancyEvaluator` compares those summaries with explicitly supplied complete virtual cohorts of the appropriate size. A virtual individual's parameter posterior cannot be relabelled as a population. Population-model and visit/time-mapping evidence are required. The comparison is descriptive; it is not an individual-observation likelihood, a fitted clinical BTK model or a safety/efficacy conclusion. Portable tests use labelled synthetic cohort predictions only. Missing patient-level and exposure inputs remain recorded in the manifest.

## Application integration

The new `VivoResearchCLICommands` is routed alongside the existing chemistry, posterior, engagement, MD and structure commands. It adds `finite-drug-run`, `engagement-design`, `occupancy-benchmark-inspect`, `occupancy-benchmark-compare` and `research-help`, reusing the existing bounded document I/O and artifact store. CLI source was syntax-parsed, but the full Apple executable was not built or run here.

Examples and invocation instructions are in `Examples/target-engagement/PORTABLE_RESEARCH.md` and `Examples/btk-aggregate-benchmark/README.md`. The portable runner is `bash Tools/Posterior/run_portable_checks.sh <new-output-directory>`. Reference generation with mpmath/NumPy/SciPy is optional validation tooling, not a production dependency.

Remaining priorities include parameter-varying Metal posterior execution with separate numerical qualification, covariance/censoring-aware information diagnostics, hierarchical biological variation, evidence-backed exposure/kinetics for actual cohorts, and a measured cellular-response endpoint. This increment supplies executable, checked components for those tasks without claiming the missing biological evidence already exists.
