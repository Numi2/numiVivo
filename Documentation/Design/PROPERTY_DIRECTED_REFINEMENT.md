# Property-directed electronic-space refinement

Status: bounded native numerical implementation with selected-CI and orbital-information support, content-addressed workflow operations, CLI routing, anchor-Hamiltonian export, and regression fixtures. Source implementation and portable regression agreement do not establish the full Apple package, GPU performance, chemical accuracy, or biological applicability.

## Scientific scope

The controller selects electronic detail using changes in an explicitly requested fixed-geometry profile. Correlation diagnostics explain represented electronic dependencies; they do not replace the physical comparisons. Inputs are prepared Hamiltonians in declared orbital frames, an ordered path with physical cross-geometry overlaps, mandatory active orbitals, complete transport groups, a predeclared candidate-block universe, discovery points, held-out confirmation points, and numerical/resource limits.

The first and last points are endpoints. The requested interior barrier point is a comparison point, not an optimized or characterized transition state. A selected space is accepted only within the declared candidate universe. Unexamined orbitals are listed. Finite electronic sensitivity is neither a reference-error bound nor an activation Gibbs free energy or kinetic rate.

## Existing numerical authorities

`VivoDirectCI.diagonalize` supplies the same Davidson implementation for explicitly selected determinant spaces and complete small sectors. `VivoDirectHamiltonian.connections` supplies the common Slater–Condon action. `VivoSelectedCI` grows a bounded determinant space without enumerating the complete sector. Its connected scan accumulates coherent amplitudes from all selected determinants before forming the omitted residual or ranking new determinants.

Projected residual, external residual and their full norm remain distinct. PT2 is unavailable when nonzero coupling has a near-zero or nonnegative denominator; an arbitrary denominator floor is not substituted into a reported physical correction. Near-degenerate candidates are promoted. Truncated external-candidate storage throws rather than reporting an incomplete residual as complete. Selection capacity, iteration capacity, and global operator-work exhaustion are explicit failure/termination states. Residual convergence does not independently certify ground-root identity, total spin, basis accuracy or molecular validity.

`VivoSelectiveOrbitalInformation` requests exact one-/two-spatial-orbital marginals from the existing CI wavefunction trace. It does not reconstruct a general two-spatial-orbital density operator from ordinary one-/two-particle RDMs. Occupation, double occupation, single entropy, pair entropy and mutual information retain separate definitions. An incomplete pair query must not be interpreted as a sparse matrix with zero-valued omitted edges.

## Refinement sequence

`VivoPropertyDirectedSpace` transports complete orbital subspaces through `VivoOrbitalPathTransport`, preserving the existing physical state-overlap calculations. It never multiplies neighboring overlap matrices to invent an unchecked nonadjacent physical overlap.

Every expansion uses the same discovery geometries and nested active spaces. Occupied orbitals promoted into the active space leave the doubly occupied core; they are not silently made empty. Whole declared orbital groups are promoted together. The score uses the maximum normalized forward, reverse, reaction-energy and relative-profile change divided by charged computational work. Absolute-energy changes remain visible but are not substituted for barrier sensitivity. Cost is a declared operator/auxiliary work accounting measure, not a wall-clock speed prediction.

All predeclared candidates are considered. Individually small changes are insufficient: the controller must also compare the complete candidate union. When discovery suggests adequate sensitivity, it evaluates all path points and all adjacent physical state-overlap edges, including the held-out points. A failed holdout terminates; the same points are not reused for tuning while still described as held out. If a candidate union cannot fit the allowed space, no low-sensitivity conclusion is issued.

Current solvers for this loop are a single tightly converged direct-CI root or a tightly converged selected-CI root. Existing CASSCF, multistate CASSCF and ECC-DMET are not replaced, but they are not silently inserted as alternative refinement backends.

## Workflow and immutable handoff

`VivoCorrelationRefinementOperations` supplies native operations for selected CI, selective orbital information and property refinement. `VivoNativeSolverDispatch` recognizes the versioned selected-CI request explicitly and rejects unknown solver schemas/fields rather than retrying an older method. A request declares whether a bounded selected-CI probe may be retained while unconverged. The default accepted-solver policy rejects it.

The shared workflow binds inputs, configuration, numerical budget, implementation fingerprint and output kinds. Its cache path validates outputs again. Property-refinement validation currently replays the bounded deterministic calculation; this is deliberately conservative and not an inexpensive cached certificate. Execution and validation each have their own invocation of the declared numerical budget; logged solver work does not include OS allocation overhead or claim measured hardware performance.

`chemistry-refine` persists incomplete exploratory reports as such and returns exit code 2. Receipt existence authenticates data, not scientific success. `chemistry-export-space` accepts only a sensitivity-established, reconstructible report and an explicitly declared path point. It materializes the chosen projected Hamiltonian and attaches request, evidence, parent-Hamiltonian, orbital-frame and partition fingerprints to provenance. The existing `chemistry-solve` interface consumes this Hamiltonian without another electronic solver stack.

Export does not extrapolate to new geometries, change a running trajectory, or qualify a new production Hamiltonian. Any MD/QM/MM production-model change must remain separately versioned and independently equilibrated/qualified. Electronic anchor sensitivity must not bypass the existing PMF, transmission, replicated-rate or chemical-qualification requirements.

## Entry points

Library: `VivoSelectedCI`, `VivoSelectiveOrbitalInformation`, `VivoPropertyDirectedSpace`, `VivoCorrelationRefinementOperations.materializeHamiltonian`.

The general workflow catalog registers the same operations through `VivoPlatformRefinementOperations`. `workflow-template property-refinement` composes refinement, anchor export, selected CI, Hamiltonian-validated state extraction and selective orbital information. Failed confirmation blocks the downstream calculation; an unconverged probe cannot enter the accepted-state extraction node.

CLI: `chemistry-solve`, `chemistry-correlations`, `chemistry-refine`, `chemistry-export-space`, plus templates listed by `chemistry-help`. See the [complete algebraic example](../../Examples/property-directed-refinement/README.md).

## Qualification boundary and remaining development

The portable harness exercises actual production numerical source with deterministic algebraic fixtures and negative tests. It is not a mocked solver, but it is also not the complete Apple package or a protein reaction. Apple artifact/CLI integration has a separate test driver and full-module tests. Retain exact source hashes, compiler/platform information and failures when running either.

The following remain separate development and scientific-qualification work: automatic chemical orbital seeding from prepared molecules; scalable broad-system probes beyond the finite orbital/ERI representation; multistate and optimized-orbital refinement; adaptive fragment/bath/QM-region changes through the common workflow; cost-calibrated allocation between electronic refinements and actual sampling/transmission trajectories; overlap-qualified free-energy corrections; downstream kinetic-observable sensitivities; and measured Metal/TensorOps acceleration. No speedup, enzyme-barrier accuracy, or experimental agreement is asserted here.

## Methodological sources

The implementation combines established building blocks; their combination is a development hypothesis, not a claim of unique scientific priority.

- Selected CI and perturbative selection: Holmes, Tubman and Umrigar, *Heat-bath Configuration Interaction*, arXiv:1606.07453. The native implementation uses its documented deterministic residual/contribution ranking; it is not presented as a reproduction of every HCI screening algorithm.
- Orbital-information-driven representation optimization: *Quantum Information-Assisted Complete Active Space Optimization*, arXiv:2309.01676. Its formal construction does not make an arbitrary approximate mutual-information cut a certified reaction-rate bound.
- Existing NumiVivo conventions: [ECC-DMET integration](ECC_DMET_INTEGRATION.md), [native correlated execution](NATIVE_CORRELATED_EXECUTION.md), and [artifacts and provenance](ARTIFACTS_AND_PROVENANCE.md).
