# QM/MM activation free energies and context-qualified rates

This workflow sits on top of NumiVivo's existing `VivoMDCandidateForceProvider` contract. It does not introduce another molecular-dynamics engine, periodic electrostatics implementation, polarization solver or adaptive-region solver.

The purpose is narrow: turn a prepared and already-qualified chemical transformation into a sampled activation PMF under the same QM/MM Born–Oppenheimer Hamiltonian used for dynamics, then produce a context-bound unimolecular rate without treating biased trajectory time as reaction time.

## Execution chain

```text
prepared molecular structure + native classical system
                         |
                         v
mapped reaction connectivity / qualified saddle
                         |
                         v
QMMM Hamiltonian plan + electronic configuration
                         |
                         v
VivoMDCandidateForceProvider
(fixed / periodic / polarizable / adaptive implementations remain authoritative)
                         |
                         v
analytic conservative umbrella composition
                         |
                         v
Langevin NVT production windows
                         |
                         v
exact accepted-state block checkpoints + scalar CV traces
                         |
                         v
unbinned MBAR + overlap/autocorrelation diagnostics
                         |
                         v
qualified PMF bound to exact system/provider/dynamics/reaction/environment
                         |
                         v
flux-normalized conventional TST
                         |
                         v
existing VivoTransitionStateDerivation / kinetic context
```

## Reaction coordinates

The initial production contract supports:

- one interatomic distance; and
- a difference of two interatomic distances.

Atom indices address the canonical prepared molecular structure. They are resolved once to physical classical particles, so virtual-site particle ordering cannot change the reaction coordinate. Periodic displacements use the same exact minimum-image geometry as the retained Hamiltonian.

For these two coordinate classes, the coordinate mass metric is geometry independent. NumiVivo therefore computes

```text
g_xi = sum_i |d xi / d r_i|^2 / m_i
```

from the exact physical particle masses rather than introducing a fitted kinetic prefactor.

## Conservative umbrella composition

`VivoQMMMUmbrellaBias.provider` wraps an existing complete Born–Oppenheimer force provider. For window center `xi0` and force constant `k`, the added potential is

```text
U_bias = 1/2 k (xi - xi0)^2
```

and the added physical-particle force is its analytic derivative through the reaction-coordinate Jacobian.

The wrapped provider still owns the original QM/MM electronic contribution. Retained classical forces remain owned by `VivoMDMetalRuntime`. Periodic Ewald/PME, induced dipoles, adaptive partitions, link-atom redistribution and dependent-site redistribution are therefore not reimplemented by the free-energy layer.

The current umbrella wrapper is NVT/NVE only. A barostat changes a periodic minimum-image reaction coordinate through the trial cell; that affine derivative is not approximated or silently omitted.

## Durable production

`VivoQMMMFreeEnergyArchiveRunner` executes windows sequentially to bound Apple unified-memory pressure. Production is split into configured blocks. A durable checkpoint is published only after a complete block and contains:

- the exact accepted-state MD checkpoint;
- current window and production-step position;
- scalar reaction-coordinate samples;
- sampled potential energies; and
- the immutable execution request identity.

Resume verifies the retained classical-system fingerprint, the composed umbrella-provider execution fingerprint and the accepted MD clock before continuing. A rejected or cancelled block does not replace the last durable block.

## MBAR analysis and acceptance

`VivoQMMMFreeEnergy.analyze` performs the following operations:

1. conservative autocorrelation thinning for each window;
2. unbinned multistate Bennett acceptance-ratio fixed-point reconstruction;
3. adjacent-window phase-space overlap checks;
4. kernel reconstruction of a one-dimensional PMF;
5. local effective-sample diagnostics in the reactant basin and at the dividing surface; and
6. a conditional finite-sample uncertainty estimate.

A profile that fails overlap, effective-sample, autocorrelation or MBAR residual checks is returned as nonconverged. It cannot be promoted by the rate API.

The stored `activationFreeEnergyKJPerMol` is a PMF profile-height diagnostic: dividing-surface PMF minus the minimum within the declared reactant range. It is useful for inspecting the profile but is **not directly used as the kinetic activation free energy**.

## Flux-normalized rate

For a one-dimensional coordinate `xi`, the conventional classical transition-state flux is calculated from the converged PMF as

```text
k_TST = [rho(xi*) / integral_R rho(xi) dxi]
        * sqrt(R T g_xi / (2 pi))
```

where the bracketed term has units of inverse nanometres and the positive-velocity factor is in nanometres per picosecond under NumiVivo's kJ/mol, Da, nm and ps units. The result is converted to inverse seconds.

This avoids assigning a kinetic meaning to umbrella elapsed time and avoids treating PMF peak-minus-basin-minimum as an Eyring barrier. For compatibility with the existing kinetic infrastructure, NumiVivo constructs the unique nonnegative Eyring-equivalent activation free energy whose `kBT/h` expression reproduces this direct PMF flux exactly. The explicitly supplied classical transmission probability is then applied through the existing `VivoTransitionStateDerivation` path.

A flux greater than the Eyring thermal prefactor implies a negative/barrierless equivalent activation free energy and is rejected by this contract rather than forced into the existing bound-complex TST model.

## Reaction and environment qualification

A rate-producing PMF is not accepted from a bare profile JSON. `VivoQMMMFreeEnergyQualification` binds it to:

- canonical structure fingerprint;
- retained classical-system fingerprint;
- Born–Oppenheimer provider fingerprint;
- MD numerical configuration fingerprint;
- chemical-state identifier;
- typed environment (`explicitSolution` or `proteinEnvironment`);
- environment/host identifier;
- deterministically revalidated mapped reaction connectivity;
- declared reactant and product endpoint identities; and
- exact reactive atom mapping.

The qualifying API reruns `VivoReactionConnectivity.validate`, so a stored `converged=true` flag is not treated as authority. The reaction-coordinate atoms must occur in the qualified saddle's structure mapping.

This is the path that permits `proteinEnvironment`. The older local-RRHO `VivoPreparedReactionRate` remains intentionally unable to relabel an isolated/continuum stationary-point calculation as a protein activation free energy.

## CLI analysis and rate conversion

Offline retained traces can be reconstructed with:

```sh
numivivo qmmm-free-energy-analyze analysis-request.json \
  --output activation-pmf.json
```

A context-qualified PMF can then be converted through the kinetic contract with:

```sh
numivivo qmmm-free-energy-rate rate-request.json \
  --output rate.json
```

Optionally, `--kinetics kinetic-pack.json` writes a separate kinetic pack with only the bound-complex inactivation parameter replaced. Association, dissociation and target-turnover parameters are retained.

Actual QM/MM umbrella dynamics are exposed through the Swift execution API because the authoritative BO provider is an executable, fingerprinted closure assembled from the existing fixed/polarizable/adaptive QM/MM implementations. Offline JSON must not pretend to serialize such a provider.

## What this does not establish

Passing the numerical gates establishes convergence only for the declared Hamiltonian, reaction coordinate, windows and sampled state. It does not establish:

- correct protonation or tautomer populations;
- complete reactive-conformer populations;
- absence of an alternative mechanism;
- electronic-structure or force-field accuracy;
- calibrated adaptive-region promotion free energies;
- negligible finite-size error;
- dynamical recrossing or tunnelling when transmission is assumed; or
- experimental predictive accuracy.

A realistic protein-rate qualification therefore still needs a chemically justified prepared system, mapped reaction, window design, convergence campaign, sensitivity checks and external chemical validation. The software path now keeps those questions separate instead of filling missing terms with electrostatic-solvent convergence or a local harmonic barrier.
