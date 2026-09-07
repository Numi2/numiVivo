# QM/MM dividing-surface transmission and recrossing

NumiVivo supports a classical dynamical transmission coefficient built on the same fingerprinted Born–Oppenheimer Hamiltonian used for a qualified QM/MM PMF. This closes the software path from a sampled dividing surface to a computed recrossing correction without interpreting biased umbrella time as reaction time.

The supported chain is:

```text
qualified QM/MM PMF
        |
        v
independent harmonic window at the declared dividing surface
        |
        v
step-separated canonical surface states
        |
        v
fresh unbiased NVE launches from exact stored positions/velocities
        |
        v
flux-weighted product/reactant commitment outcomes
        |
        v
computed classical transmission coefficient kappa
        |
        v
flux-normalized QM/MM rate request
        |
        v
replica qualification / kinetic-pack application
```

## Surface ensemble

`VivoQMMMSurfaceEnsembleRequest` is a typed immutable request. It contains the qualified PMF, the exact NVT dynamics configuration used by that PMF, an initial state, an independent random seed, a harmonic window centered at the PMF dividing surface, bounded equilibration/production schedules, accepted-step separation, the surface tolerance and the requested number of states.

`VivoQMMMSurfaceEnsemble.run` composes the dividing-surface umbrella with the supplied `VivoMDCandidateForceProvider` and executes through `VivoMDMetalRuntime`. It records accepted positions, velocities, cell, clock, reaction coordinate and the source checkpoint fingerprint. Surface evidence is bound to the qualified PMF evidence, retained classical system, base BO provider, biased surface provider and exact execution fingerprint.

The stored checkpoint fingerprint is provenance for the biased trajectory. It is not used to restore biased Hamiltonian state for shooting.

## Unbiased NVE shooting

`VivoQMMMDynamicalTransmission.run` validates that the PMF, surface ensemble, retained system and base BO provider identify the same Hamiltonian. The shooting configuration must be force-equivalent to the surface dynamics except that it is NVE with no thermostat.

For every surface state, NumiVivo evaluates the reaction-coordinate velocity. Negative product-directed velocities are time reversed. A new `VivoClassicalInitialState` is then created from the exact stored physical positions, cell, time and selected velocities, and a fresh `VivoMDMetalRuntime` is started with the **unbiased base provider**. The biased checkpoint itself is never restored into the NVE trajectory.

Each trajectory is propagated until it satisfies the declared product commitment range, returns to the reactant commitment range, or reaches the finite shooting horizon. Commitment can require multiple consecutive observations. Candidate rejection is a numerical failure rather than an implied recrossing event.

## Flux-weighted estimator

Surface states are weighted by the magnitude of their product-directed reaction-coordinate velocity. The computed coefficient is

```text
kappa = sum(product-committed flux weights) / sum(all surface flux weights)
```

The result also records the effective flux-weighted sample count, unresolved flux fraction, coefficient standard error, counts for product commitment/reactant recrossing/unresolved outcomes and immutable trajectory evidence.

A converged result requires the declared minimum effective flux samples, maximum unresolved flux fraction and maximum coefficient standard error. A zero product coefficient is retained as inspectable evidence but does not qualify rate application.

This is a **classical recrossing correction**. It does not represent quantum tunnelling, alternate reaction pathways, missing chemical states or errors in the electronic/force-field model.

## Rate application

`VivoQMMMDynamicalTransmission.applying` accepts only a reconstructable converged result whose PMF evidence exactly matches the target `VivoQMMMFreeEnergyRateRequest`. It changes the rate request from assumed/external transmission to `.calculated`, stores the computed coefficient and binds the kinetic evidence to the transmission evidence fingerprint. The existing flux-normalized rate calculation is then used unchanged.

PMF profile-height uncertainty is not propagated into this flux-normalized rate as though it were an Eyring barrier uncertainty. Replica dispersion, finite PMF sampling, transmission sampling and model sensitivity remain separate evidence channels.

## Workflow artifacts

Executable force providers contain runtime closures and therefore are not serialized into the workflow DAG. Surface collection and NVE shooting remain Swift execution operations with an actual `VivoMDCandidateForceProvider` and Metal runtime.

Once trajectories exist, two deterministic artifact nodes are available:

| Operation | Purpose |
|---|---|
| `vivo.platform.qmmm-transmission-analyze` | Reconstruct a `VivoQMMMDynamicalTransmissionResult` from a typed transmission request, exact retained-system/provider fingerprints and trajectory evidence. |
| `vivo.platform.qmmm-apply-transmission` | Validate that result and produce a new `VivoQMMMFreeEnergyRateRequest` containing the calculated transmission coefficient. |

The serializable input types are `VivoQMMMDynamicalTransmissionAnalysisRequest` and `VivoQMMMComputedTransmissionApplicationRequest`. They do not claim to serialize or reproduce the BO force-provider implementation.

## Qualification expectations

Software qualification should include both layers:

1. Deterministic estimator tests that reconstruct typed surface evidence and verify flux weighting, convergence gates, PMF/provider transplantation rejection and rate application.
2. A real Apple Metal integration test that runs the biased surface collector and hands the stored physical state to fresh unbiased NVE shooting through the existing runtime.

Scientific production qualification additionally requires independent surface and PMF trajectories, coordinate/window sensitivity, shooting-horizon and commitment-basin sensitivity, electronic-model and QM-region sensitivity, alternate chemical states/pathways, and external chemical or experimental validation.
