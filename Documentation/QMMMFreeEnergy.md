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
VivoQMMMFreeEnergyForceFactory
(fixed / PME / polarizable / adaptive assembly)
                         |
                         v
VivoMDCandidateForceProvider
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
independent replica agreement
                         |
                         v
flux-normalized conventional TST
                         |
                         v
existing kinetic context / covalent kinetic pack
```

## Reaction coordinates

The initial production contract supports one interatomic distance and a difference of two interatomic distances. Atom indices address the canonical prepared molecular structure. They are resolved once to physical classical particles, so virtual-site particle ordering cannot change the reaction coordinate. Periodic displacements use the same exact minimum-image geometry as the retained Hamiltonian.

For these coordinate classes the coordinate mass metric is geometry independent:

```text
g_xi = sum_i |d xi / d r_i|^2 / m_i
```

NumiVivo computes it from the exact physical particle masses rather than introducing a fitted kinetic prefactor.

## Unified force setup

`VivoQMMMFreeEnergyForceFactory` is the single programmatic assembly surface for production sampling. Fixed-region mode accepts an existing `VivoQMMMHamiltonianPlan` and selects direct analytic HF or reciprocal Metal PME from the electronic configuration. Adaptive mode delegates to `VivoAdaptiveQMMMDynamics` and returns its retained permanent-classical baseline plus the complete multipartition correction provider.

The factory verifies that the returned provider owns the exact retained system that the MD runtime will integrate, and checks finite versus periodic geometry and PME requirements before sampling starts. Polarization, adaptive partitioning, link atoms and virtual sites remain implemented by their existing authoritative components.

## Conservative umbrella composition

`VivoQMMMUmbrellaBias.provider` wraps an existing complete Born–Oppenheimer force provider. For window center `xi0` and force constant `k`:

```text
U_bias = 1/2 k (xi - xi0)^2
```

The added physical-particle force is the analytic derivative through the reaction-coordinate Jacobian. Retained classical forces remain owned by `VivoMDMetalRuntime`. Periodic Ewald/PME, induced dipoles, adaptive partitions, link-atom redistribution and dependent-site redistribution are not reimplemented by the free-energy layer.

The current umbrella wrapper is NVT/NVE only. A barostat changes a periodic minimum-image reaction coordinate through the trial cell; that affine derivative is not approximated or silently omitted.

## Durable production

`VivoQMMMFreeEnergyArchiveRunner` executes windows sequentially to bound Apple unified-memory pressure. Production is split into configured blocks. A durable checkpoint is published only after a complete block and contains the exact accepted-state MD checkpoint, current window and production-step position, scalar reaction-coordinate samples, sampled potential energies and immutable execution-request identity.

Resume verifies the retained classical-system fingerprint, composed umbrella-provider execution fingerprint and accepted MD clock before continuing. A rejected or cancelled block does not replace the last durable block.

## MBAR analysis and acceptance

`VivoQMMMFreeEnergy.analyze` performs conservative autocorrelation thinning, unbinned MBAR reconstruction, adjacent-window overlap checks, kernel reconstruction of the one-dimensional PMF, local effective-sample diagnostics and a conditional finite-sample uncertainty estimate.

A profile that fails overlap, effective-sample, autocorrelation or MBAR residual checks is nonconverged and cannot be promoted by the rate API. The stored `activationFreeEnergyKJPerMol` is only a PMF profile-height diagnostic: dividing-surface PMF minus the minimum within the declared reactant range. It is not directly used as the kinetic activation free energy.

## Flux-normalized rate

For a one-dimensional coordinate `xi`, conventional classical transition-state flux is calculated as

```text
k_TST = [rho(xi*) / integral_R rho(xi) dxi]
        * sqrt(R T g_xi / (2 pi))
```

The bracketed term has units of inverse nanometres and the positive-velocity factor is in nanometres per picosecond under NumiVivo's kJ/mol, Da, nm and ps units. The result is converted to inverse seconds.

This avoids assigning kinetic meaning to umbrella elapsed time and avoids treating PMF peak-minus-basin-minimum as an Eyring barrier. For compatibility with existing kinetic infrastructure, NumiVivo constructs the unique nonnegative Eyring-equivalent activation free energy whose `kBT/h` expression reproduces this direct PMF flux. The explicitly supplied classical transmission probability is then applied through `VivoTransitionStateDerivation`.

A flux greater than the Eyring thermal prefactor implies a negative/barrierless equivalent activation free energy and is rejected by this contract rather than forced into the bound-complex TST model.

## Reaction and environment qualification

A rate-producing PMF is not accepted from a bare profile JSON. `VivoQMMMFreeEnergyQualification` binds it to canonical structure fingerprint, retained classical-system fingerprint, Born–Oppenheimer provider fingerprint, MD numerical configuration fingerprint, chemical-state identity, typed environment, host identity, deterministically revalidated mapped reaction connectivity, endpoint identities and exact reactive atom mapping.

The qualifying API reruns `VivoReactionConnectivity.validate`, so a stored `converged=true` flag is not authority. Reaction-coordinate atoms must occur in the qualified saddle's structure mapping. This is the path that permits `proteinEnvironment`; the local-RRHO prepared-rate path remains unable to relabel an isolated calculation as a protein activation free energy.

## Independent PMF replicas

A single numerically converged PMF is useful evidence but is not the production replication gate. `VivoQMMMReplicatedFreeEnergyRate` accepts two or more independently qualified rate requests and requires all replicas to share the exact kinetic context, environment, reaction coordinate, analysis configuration, structure/system/provider/dynamics fingerprints, reaction-connectivity binding and transmission model. Stochastic window seeds must be disjoint across all replicas.

Replica agreement is evaluated directly in `ln(k)` because rate is exponentially sensitive to free energy. A secondary PMF profile-barrier range is retained as a diagnostic. The result records the geometric-mean rate, sample standard deviation across independent replica log rates, mean within-replica conditional variance, their quadrature combination, explicit acceptance issues and an immutable evidence fingerprint.

A replicated result can exist with `converged=false` so the disagreement remains inspectable, but it cannot replace the inactivation parameter of a `VivoCovalentKineticPack` until the configured agreement gates pass.

## CLI analysis and rate conversion

```sh
numivivo qmmm-free-energy-analyze analysis-request.json \
  --output activation-pmf.json

numivivo qmmm-free-energy-rate rate-request.json \
  --output rate.json

numivivo qmmm-free-energy-replicated-rate replicas.json \
  --output replicated-rate.json
```

Both rate commands accept `--kinetics kinetic-pack.json`. The single-PMF command writes a pack using that qualified rate. The replicated command writes a pack only when independent-replica agreement passes. Association, dissociation and target-turnover parameters are retained.

Actual QM/MM umbrella dynamics are exposed through the Swift execution API because the authoritative BO provider is an executable, fingerprinted closure. Offline JSON must not pretend to serialize such a provider.

## What this does not establish

Passing the numerical gates establishes convergence only for the declared Hamiltonian, reaction coordinate, windows and sampled state. It does not establish correct protonation or tautomer populations, complete reactive-conformer populations, absence of alternative mechanisms, electronic-structure or force-field accuracy, calibrated adaptive-region promotion free energies, negligible finite-size error, negligible recrossing or tunnelling when transmission is assumed, or experimental predictive accuracy.

A realistic protein-rate qualification therefore still needs a chemically justified prepared system, mapped reaction, window design, multiple independent production replicas, convergence and sensitivity campaigns, and external chemical validation. The software path keeps those questions explicit instead of filling missing terms with continuum-solvent convergence or a local harmonic barrier.
