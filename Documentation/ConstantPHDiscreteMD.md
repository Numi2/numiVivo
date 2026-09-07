# Discrete constant-pH MD/MC

`VivoConstantPHDiscreteMD` is NumiVivo's first native dynamical bridge between explicitly declared chemical states and pH-dependent sampling. It samples a semigrand ensemble by alternating fixed-state molecular dynamics with discrete Metropolis-Hastings state proposals.

This is a runtime API, not a serialized workflow calculation. Each chemical state owns an executable complete Hamiltonian callback and a propagation callback. Those closures may wrap the existing classical, QM/MM, PME, polarizable or adaptive force-provider machinery, but executable force providers are not represented as JSON.

## Statistical model

For state `i`, let

- `U_i(R)` be the complete physical potential energy returned by its executable Hamiltonian at conformation `R`;
- `B_i` be its evidence-bound reference semigrand correction at `referencePH`;
- `n_i` be its relative bound-proton count;
- `T` be temperature;
- `pH` be the target pH.

The state weight is proportional to

```text
exp[-( U_i(R) + B_i + n_i R T ln(10) (pH-referencePH) ) / (R T)]
```

Only differences in proton count and reference correction matter. A one-unit pH increase changes the relative reservoir weight of two states differing by one bound proton by exactly a factor of ten.

The proposal graph is explicit and reciprocal. A move is selected uniformly from the current state's neighbors, so the Metropolis-Hastings acceptance includes

```text
ln q(reverse)/q(forward) = ln(degree(current)/degree(proposed))
```

This prevents proposal-graph degree from biasing the stationary distribution.

## MD/MC cycle

For every attempt:

1. Propagate `mdStepsPerAttempt` under the **current state's complete Hamiltonian**.
2. Retain the accepted coordinates, velocities, cell and clock as one state-independent physical snapshot.
3. Select one reciprocal graph neighbor.
4. Evaluate the current and proposed complete Hamiltonians on the **exact same snapshot**.
5. Add the explicit reference-bias and proton-reservoir contributions.
6. Apply the Hastings-corrected Metropolis criterion.
7. If accepted, change only the chemical-state identity. Coordinates and velocities remain unchanged.
8. Commit state identity, physical snapshot and RNG state atomically.

Kinetic energy cancels from the state move because v1 requires all executable states to share the same physical particle and mass manifold.

## Identity and restart

Each executable state declares both a common `physicalManifoldFingerprint` and its own `hamiltonianFingerprint`. A checkpoint binds:

- the complete sampling configuration;
- all executable Hamiltonian identities;
- the common physical manifold;
- current chemical state;
- exact physical coordinates, velocities, cell and MD clock;
- RNG state;
- visit and acceptance accounting;
- the complete proposal/acceptance trace.

A checkpoint cannot be resumed under another state graph, pH, Hamiltonian set or physical manifold.

## Evidence requirements

`VivoConstantPHStateDefinition` requires explicit origin/evidence for every reference-state correction. The sampler does not turn an assumed pKa or guessed protonation preference into calculated evidence merely because the Monte Carlo algebra is correct.

The complete physical energy callbacks must use one coherent energy convention. State-dependent classical parameters, QM/MM electron sectors, boundary-charge models, polarization models and periodic corrections must be represented by the state's executable Hamiltonian identity. The sampler does not repair mismatched energy references.

## Deliberate v1 boundary

The v1 implementation requires one fixed physical particle/mass manifold across states. It therefore does **not** yet support proposals that literally create/delete atoms or change masses. Production protonation states must currently be represented on a common latent physical topology by the caller's Hamiltonian implementation.

It also does not:

- enumerate protonation states or tautomers;
- predict pKa values;
- calculate the reference semigrand corrections;
- decide which titratable sites are chemically relevant;
- prove that the supplied state graph is complete;
- provide pH replica exchange yet;
- infer reaction pathways or proton-transfer mechanisms.

These are separate scientific capabilities rather than hidden defaults.

## Relation to static state thermodynamics

`VivoQMMMChemicalStateThermodynamics` performs deterministic pH reweighting when relative state free energies are already known. `VivoConstantPHDiscreteMD` instead samples conformations and state identity together through alternating MD/MC. The static calculator remains useful for independent thermodynamic checks and externally supplied state free energies.

For a realistic protein workflow, the next development stages are:

```text
explicit state catalog
      -> discrete constant-pH MD/MC
      -> pH replica exchange / exchange diagnostics
      -> state populations + uncertainty
      -> state-specific QM/MM reaction rates
      -> rapid-equilibrium or explicit chemical-state kinetic network
```

The reference statistical formulation follows the discrete constant-pH semigrand MD/MC family: fixed-state MD is periodically interrupted by protonation-state Monte Carlo moves, with pH represented through the proton chemical potential.