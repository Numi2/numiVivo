# Mapped reaction connectivity

## Scope

`VivoReactionConnectivity` extends local saddle characterization with explicit
endpoint assignment and sampled displacement/step refinement. The implementation
must be qualified by the `Mapped reaction connectivity` workflow; this document
contains no asserted passing run, measured error or paper-reproduction result.

The endpoint references are independently characterized minima. A dissociated
endpoint is a list of separately characterized molecules with an explicit map
from each molecule's atoms to the saddle atom identifiers. A small total gradient
at widely separated fragments is not treated as a single-molecule RRHO minimum.

No chemical identity is inferred from a generic bond-length cutoff. No atom
permutation is introduced to improve an endpoint fit. Element, isotope, spin
population and Gaussian basis primitives are checked against the mapped saddle
atoms. Both endpoints must account for every saddle atom exactly once.

## Physical-model consistency

The saddle and every component use the same electronic solver, SCF settings,
solvent configuration and correlated-solvent settings. Gas FCI and equilibrium
FCI solvent models remain different from RHF-reference-frozen solvent models.
The existing nuclear qualification rejects inappropriate removal of rigid modes
in frozen external fields; that restriction is retained.

Independent anchored ECC component surfaces are not accepted by this version.
An atom map does not establish a shared orbital gauge or a consistent frozen-core
energy across different fragment surfaces. Those calculations require an explicit
common reaction model rather than relaxing the compatibility check.

## Endpoint acceptance

For each final branch geometry, all supplied endpoint candidates are evaluated.
Exactly one candidate must satisfy every criterion:

* Each mapped component fits its characterized reference within the declared
  mass-weighted proper-rotation RMSD. Quaternion alignment permits rotations and
  translation, not reflection or atom permutation.
* The final gradient satisfies the configured endpoint threshold.
* The total electronic energy agrees with the sum of component reference
  electronic energies within the declared bound.
* Different endpoint components are separated by at least the declared minimum
  intercomponent distance.

The energy test matters even when components are far apart: distance by itself
does not establish negligible interaction or relaxed intramolecular structure.
A unique match is required. Two admissible matches are an ambiguity, not grounds
to pick the first candidate. The reverse and forward branches must reach different
mapped endpoint identifiers. Reference minima and the saddle are reconstructed
through the existing nuclear validator before descent evidence is accepted.

## Independent numerical refinements

Four calculations start from the same validated saddle:

1. The declared initial unstable-mode displacement and integration step.
2. Half the displacement at the original integration step.
3. The original displacement at half the integration step.
4. Half both quantities.

This varies the two numerical choices independently instead of comparing only
a diagonal refinement where two errors could cancel. Each run uses the production
mass-weighted descent implementation and the same physical arc-length limit.
Failure to reach endpoint gradient thresholds is not declared convergence.

Endpoint labels determine branch correspondence; the arbitrary sign of an
unstable eigenvector does not determine which chemical endpoint is called forward.
Every pair of trials is compared in both directions, giving twelve comparisons.

Comparisons interpolate energies and mapped pair distances onto explicitly
sampled common physical arc lengths. They do not extrapolate one trajectory to
cover the other. Both the energy difference and maximum mapped pair-distance
difference must remain below their declared tolerances. The common interval must
have a minimum declared length. These are sampled error checks, not a uniform
continuous-curve error bound or a proof that no unobserved intermediate event exists.
Proper-rotation chirality is checked for the endpoints; intermediate refinement
observables are the explicitly recorded mapped pair distances and energies.

The result retains all four trajectories, both endpoint assignments per trial,
the twelve comparison records and the descent electronic-evaluation count.
`converged` is false when the recorded refinements do not satisfy the contract.
A failed refinement is diagnostic evidence, not a successful connection.

## Numerical authority and workflow

`VivoReactionConnectivityWorkflow.operation(implementationFingerprint:)` plugs
into `VivoChemistryWorkflow`. The task consumes one canonical request artifact
containing the exact endpoint records, atom maps, saddle and numerical settings.
It publishes a connection only when refinement succeeds. Output validation reruns
the physical calculation, so edited endpoint labels or success flags cannot
substitute for computed evidence. The operation rejects a nested endpoint budget
that differs from the task's numerical resource contract.

`VivoReactionConnectivityWorkflow.plan` constructs the corresponding one-node
artifact plan. This is a library/workflow entry point. It does not add a new
`reaction-run` enum case or claim that the existing command parser already accepts
this request schema.

Each of the four descent runs retains the existing nuclear solve budget. A
conservative bound on their sum is checked before execution. Independent endpoint
and saddle reconstruction have their own request-level budgets. The recorded
descent count must not be presented as total execution cost including validation.

## Executable qualification

The Swift package tests are in
`Tests/NumiVivoIntegrationTests/ReactionConnectivityTests.swift`:

```bash
swift test -c release --jobs 3 -Xswiftc -enable-testing \
  --filter ReactionConnectivityTests
```

The suite checks proper-rotation alignment, rejection of reflected/permuted mapped
structures, linear and atomic limits, malformed refinement settings, and an actual
H3/H2/H finite-basis connection with independent displacement and step refinements.
The integration calculation writes its exact request and full result when
`NUMIVIVO_TEST_ARTIFACTS` is set. Numerical success is determined by the run, not by
the presence of the fixture or this documentation.

The CI workflow retains the source commit, source archive, compiler/platform
identity, source hashes, test log and numerical outputs even on failure. It uses
the actual production nuclear and electronic implementations; there is no
analytical surrogate in place of the electronic energy surface.

## What this does not establish

A passing H3/STO-3G connection is not an acrylamide-methanethiolate or BTK result.
Mapped reaction inputs, the appropriate electronic/solvent model and independent
chemical reference results remain necessary for that reproduction.

Local connectivity and harmonic free-energy calculations are separate claims.
This operation does not produce a kinetic parameter, tunneling correction,
transmission coefficient, conformational ensemble free energy or experimentally
validated rate. It does not turn the bounded CI implementation into a scalable
whole-protein solver or certify every biology/Metal runtime in the application.
