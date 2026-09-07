# QM/MM dividing-surface transmission and recrossing

NumiVivo supports a classical dynamical transmission coefficient built on the same fingerprinted Born–Oppenheimer Hamiltonian used for a qualified QM/MM PMF. The current estimator is a paired signed reactive-flux calculation: each statistically qualified near-surface state launches exact `+v` and `-v` branches and the accepted coefficient is a stable late-time plateau, not an eventual product-commitment ratio.

The supported chain is:

```text
qualified QM/MM PMF
        |
        v
independent harmonic window at the declared dividing surface
        |
        v
finite-band reweighting + correlation-aware surface states
        |
        v
exact paired +v / -v fresh unbiased NVE launches
        |
        v
signed product-indicator history hP(+v,t) - hP(-v,t)
        |
        v
late-time transmission plateau kappa
        |
        v
surface-derived positive-flux normalization + qualified PMF
        |
        v
context-qualified rate / replica qualification / kinetic pack
```

## Surface ensemble

`VivoQMMMSurfaceEnsembleRequest` is a typed immutable request. It contains the qualified PMF, exact NVT dynamics configuration, initial state, independent random seed, a harmonic window centered at the PMF dividing surface, bounded equilibration/production schedules, minimum accepted-step separation, finite surface band, Gaussian surface-kernel bandwidth, requested state count, minimum effective sample count and autocorrelation limit.

`VivoQMMMSurfaceEnsemble.run` composes the dividing-surface umbrella with the supplied `VivoMDCandidateForceProvider` and executes through `VivoMDMetalRuntime`. The complete production coordinate trace determines a conservative autocorrelation stride. Eligible states are separated by at least the greater of that statistical stride and the requested accepted-step separation.

The selected states are reweighted by the inverse umbrella factor and a Gaussian delta-kernel around the dividing surface. Their normalized statistical weights sum to one. The artifact records effective surface samples, the canonical positive coordinate flux `0.5 <|dot(xi)|>`, its sampling error, exact coordinates/velocities/cell/clock, a reconstructable state-payload fingerprint and the biased checkpoint fingerprint.

The biased checkpoint fingerprint is execution provenance only. It is not an unbiased Hamiltonian restart artifact.

## Exact paired unbiased NVE shooting

`VivoQMMMDynamicalTransmission.run` verifies that the PMF, surface ensemble, retained system and base BO provider identify the same Hamiltonian. The shooting configuration must be force-equivalent to the surface dynamics except that it is NVE with no thermostat.

For each retained surface state, the reaction-coordinate velocity is recomputed from the mapped coordinate Jacobian and exact particle velocities and compared with the stored value. Its magnitude defines the positive flux speed. NumiVivo then constructs the product-directed velocity field and its exact time reverse.

Both branches start from new `VivoClassicalInitialState` objects containing the exact stored positions, cell and time. Each starts a fresh `VivoMDMetalRuntime` using the **unbiased base provider**. The biased umbrella checkpoint itself is never restored into either branch.

Both branches are propagated through the entire declared shooting horizon. At every observation step the artifact records the reaction coordinate and whether it lies in the declared product or reactant commitment basin. The final branch outcome is based on the requested number of consecutive terminal commitment observations. Candidate rejection is a numerical failure, not a recrossing event.

## Paired signed reactive-flux estimator

For state `i`, let `w_i` be its finite-band surface statistical weight and `v_i` its positive reaction-coordinate speed. Its flux weight is

```text
W_i = w_i |v_i|
```

At each observation time the estimator is

```text
kappa(t) = sum_i W_i [hP_i(+v,t) - hP_i(-v,t)] / sum_i W_i
```

where `hP` is one when that branch is in the product commitment basin and zero otherwise. The reported `transmissionCoefficient` is the mean of the final configured number of `kappa(t)` observations.

The result retains the entire coefficient history, plateau range, effective flux-weighted sample count, unresolved paired-flux fraction, coefficient standard error, paired branch histories and commitment counts. Acceptance requires a positive physical coefficient no greater than one, a plateau range below the configured threshold, sufficient effective flux samples, bounded unresolved flux and bounded sampling error.

A nonphysical, unstable or insufficiently sampled plateau remains inspectable evidence with `converged=false`; it cannot be applied to a kinetic rate.

This is a **classical recrossing correction**. It does not represent tunnelling, alternate mechanisms, missing chemical-state populations or errors in the underlying electronic/force-field model.

## Rate application and flux normalization

`VivoQMMMDynamicalTransmission.applying` accepts only a reconstructable converged result whose PMF evidence exactly matches the target `VivoQMMMFreeEnergyRateRequest`. It sets the transmission model to `.calculated`, binds kinetic evidence to the exact paired-shooting evidence fingerprint, and installs `VivoQMMMSurfaceFluxNormalization` derived from the same qualified surface ensemble.

The rate layer therefore uses sampled positive dividing-surface flux when computed transmission evidence is available rather than mixing an independently estimated one-dimensional kinetic prefactor with the shooting ensemble. The PMF still supplies the equilibrium surface-to-reactant probability ratio and remains bound to the same reaction coordinate and Hamiltonian.

PMF profile-height uncertainty is not propagated into the flux-normalized rate as though it were an Eyring barrier uncertainty. PMF sampling, surface-flux sampling, transmission sampling, replica dispersion and model sensitivity remain distinct evidence channels.

## Workflow artifacts and CLI

Executable force providers contain runtime closures and are not serialized into the workflow DAG. Surface collection and paired NVE shooting remain Swift execution operations with a real `VivoMDCandidateForceProvider` and Metal runtime.

Once paired branch evidence exists, deterministic artifact operations are available:

| Operation | Purpose |
|---|---|
| `vivo.platform.qmmm-transmission-analyze` | Reconstruct the v3 paired reactive-flux result from a typed transmission request, exact retained-system/provider fingerprints and paired branch evidence. |
| `vivo.platform.qmmm-apply-transmission` | Validate the result and produce the exact PMF rate request carrying calculated transmission and sampled surface-flux normalization. |

The serializable input types are `VivoQMMMDynamicalTransmissionAnalysisRequest` and `VivoQMMMComputedTransmissionApplicationRequest`. The analysis request stores `VivoQMMMTransmissionPair` values, not executable provider state.

Equivalent CLI entry points are:

```sh
numivivo qmmm-transmission-analyze transmission-evidence.json --output transmission.json
numivivo qmmm-transmission-apply transmission-application.json --output calculated-rate-request.json
```

## Qualification expectations

Software qualification has two independent layers:

1. Deterministic estimator tests reconstruct finite-band surface evidence and paired branch histories, verify a known signed-flux plateau, reject unstable plateaus and Hamiltonian transplantation, and reconstruct rate application.
2. A real Apple Metal integration test executes the biased surface collector and then launches exact paired fresh unbiased NVE branches through the existing runtime.

Scientific production qualification additionally requires independent surface and PMF trajectories, coordinate/window and surface-band sensitivity, shooting-horizon and commitment-basin sensitivity, electronic-model and QM-region sensitivity, alternate chemical states/pathways, and external chemical or experimental validation. Synthetic software fixtures do not establish those claims.
