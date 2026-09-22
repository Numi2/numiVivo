# Protein-material mechanical and thermal stress testing

## Current scope

This adds one native research workflow under `NumiVivoKit/ProteinMaterials`, not
another molecular simulator. Its first use is conditional mechanical/thermal
screening of supplied, prepared structures. It is inspired by Nature Chemistry,
DOI [10.1038/s41557-025-01998-3](https://doi.org/10.1038/s41557-025-01998-3).
It is not a reproduction of that paper or of its experimental outcomes.

Implemented source includes mass-weighted COM restraints and analytical physical
forces; exact discrete protocol-work accounting; noncovalent D-H-A geometry and
time-weighted persistence; native-contact retention, contact-distance error and
radius of gyration; whole-molecule periodic reconstruction; fixed-cell NVT stage
execution through the existing Metal MD owner; immutable checkpoint/metric journals
with verification and accepted-prefix resume; and complete-replica candidate
screening under matching declared numerical protocols.

**Local measured evidence:** 15 numerical XCTest cases passed in a standalone
Swift 6.2.1 Linux package containing the actual numerical source and existing
`VivoVector3D` declaration. The recorded source hashes are in
`Tools/ProteinMaterials/evidence/2026-09-22/portable.json`. Those tests do not build
the complete Apple package or run Metal. Additional native integration tests and
an Apple CI workflow are supplied; their existence is not a passing qualification.
Use the actual workflow outcome for the exact commit.

## Public commands

```sh
numivivo protein-stress-help
numivivo protein-stress-example --output fixture.json
numivivo protein-stress-validate fixture.json --output validation.json
numivivo protein-stress-run fixture.json --store artifacts --output receipt.json
# Use the actual immutable journal hash printed by the run, not a receipt hash.
numivivo protein-stress-verify fixture.json --store artifacts --journal JOURNAL_SHA256 --output verification.json
numivivo protein-stress-run fixture.json --store artifacts --resume JOURNAL_SHA256 --output resumed.json
numivivo protein-stress-campaign campaign.json --store artifacts --output campaign-report.json
```

The example is a two-particle harmonic fixture, not a protein. Inputs and exports
reuse `VivoMDAtomicFileExport` bounds and no-clobber atomic publication. Run and
campaign commands return 75 for an incomplete/rejected result, 1 for a CLI/input
error, and 0 only for completion. `protein-stress-validate` reports input admission,
not a simulated or experimentally validated result.

`VivoProteinStressRequest` contains the source `VivoClassicalSystem`, its original
`VivoMDConfiguration`, exact `VivoMDCheckpoint`, source description, replica ID,
explicit random seed, stages, sampling interval, physical-particle selection,
optional pull definition, explicit D-H-A triplets, geometric criteria and explicit
reference contacts. Use existing molecular preparation/MD workflows to obtain a
parameterized, relaxed, equilibrated source. The new runner never fabricates
velocities from a coordinate trajectory or declares an illustrative equilibration
duration sufficient.

## Loading and thermal protocol

A stage holds a fixed harmonic reference and Langevin target temperature for an
explicit number of integration steps. At a stage boundary the Hamiltonian and/or
thermostat parameters change explicitly while exact positions, velocities, cell,
physical time and accepted RNG step are retained. The bias force provider remains
a deterministic function of geometry and immutable configuration. It has no
mutable clock and does not weaken `VivoMDCandidateForceProvider` semantics.

For mass-weighted centers x_A and x_B, q is either their radial separation or the
projection of x_B-x_A onto a declared normalized laboratory direction. The bias is

    U(x; r) = 0.5 * k * (q(x) - r)^2

Its analytic force is distributed with the same normalized mass weights as the
coordinate. Groups are nonempty, unique, disjoint and massive. Coordinates are in
nm, stiffness in kJ mol^-1 nm^-2, forces in kJ mol^-1 nm^-1 internally and pN in
reports. The exact SI conversion is 1e24 / 6.02214076e23.

External protocol work at a fixed-coordinate parameter switch is U_new(x)-U_old(x).
The initial insertion of a nonzero bias is included. Temperature changes do not
rescale velocities instantaneously; the subsequent existing thermostat handles
dynamics. Protocol work is not heat, equilibrium folding free energy, or a
Jarzynski free-energy estimate. No such estimator is implemented here.

This is **piecewise-constant staged pulling**, not a silently approximated smooth
moving restraint. Reconstructing the runtime at each stage makes the authority
boundary explicit but has a cost. Refine the stage interval and MD timestep and
compare to an independent implementation before interpreting physical results.
There is no measured GPU speedup or real-protein throughput claim.

Reference: the official [GROMACS pulling documentation](https://manual.gromacs.org/current/reference-manual/special/pulling.html)
describes pull-coordinate geometries and nonequilibrium work. This implementation
must be compared at the same masses, coordinate, force constant, boundary
convention and protocol, not just the same initial structure.

## Periodic geometry and valid profile

`VivoMDBarostatPlan` supplies the existing deterministic molecular membership,
ordered spanning tree and closure edges. `VivoProteinWholeMoleculeLayout` walks
that tree using the existing exact `VivoPeriodicCell.minimumImage` operation,
checks non-tree edges, and rejects periodic winding or inconsistent connectivity.
The two group centers are never independently minimum-imaged; that would replace
a large molecular extension with another periodic image. All pull and analysis
atoms must belong to one connected molecule. Covalent edge images must remain
resolvable; this is not a bond-breaking or topology-changing simulation.

The first profile is fixed-cell Langevin NVT with ordinary FP32 positions. Periodic
biased runs must satisfy the existing provider's PME/electrostatic contract. LJ
requires the existing smooth switching convention. NPT bias moves, compensated
external-force coordinates, induced-dipole ownership and Drude particles are
rejected, not downgraded. Use existing NPT preparation before the stress run.
Multiple disconnected assembly members, intermolecular bond-network metrics and
arbitrary chemical typing are not claimed.

## Structural observables

D-H-A identities and donor/acceptor assignments are supplied explicitly; no
protonation or acceptor chemistry is inferred. The default geometric classification
is donor-acceptor <=0.35 nm, hydrogen-acceptor <=0.25 nm, and D-H-A >=150 degrees
(straight is 180 degrees). These are declared geometric thresholds, not a unique
physical definition of bonding or extra force-field terms. No springs are added
to stabilize hydrogen bonds.

Persistence uses explicit observation times and left-constant time weighting.
First losses and reformations are observed events; unsampled transitions and
continuous bond lifetimes are unresolved. Stage-entry and source frames can share
a time. Analyze persistence within a stage, or choose one declared boundary frame
per unique time before calling the strict persistence function.

Native-contact distances and selected mass-weighted radius of gyration are
translation/rotation invariant. Empty contacts produce missing retention and
contact-distance error, not a perfect score. Contact-distance RMS error is not
aligned-coordinate RMSD. No secondary-structure classifier or automatic
"unfolded" threshold is supplied. The largest sampled tensile force is labelled
as such, not assumed to be an unfolding force.

## Accepted-state journal and failure handling

The existing `VivoArtifactStore` owns request bytes, exact checkpoints, journal
entries, receipts and campaign reports. Each immutable journal entry links its
predecessor and records source/stage/sample/rejection identity, accepted progress,
checkpoint hash, cumulative protocol work and reconstructible structural metrics.
A unique run reference avoids collisions between concurrent runs of one request.

A numerical rejection publishes the last accepted state and a rejection
certificate. A rejected journal cannot be resumed as an eligible candidate. There
is no automatic retry, new random seed, smaller timestep or replacement replica.
A cancellation, device failure or storage failure returns a distinct incomplete
receipt and the most recent immutable durable prefix. Unsaved numerical progress
is not promoted. Resume reconstructs the original sampling grid and never
reinitializes accepted velocities or RNG step. A finished prefix remains finished.

Verification bounds the entire traversal, verifies artifact kinds and hashes,
checks exact source and stage identities, replays the declared time increments,
reconstructs each metric, verifies every zero-time transition and recomputes
protocol work. It does not independently recompute all MD forces or prove the
physical correctness of a trajectory. Content hashes establish integrity, not
trusted authorship or experimental accuracy.

## Candidate screening

A campaign supplies 2–32 distinct seed/ID replicas per candidate. Replicas within a
candidate share the same physical system and analysis definition. Across
candidates, admitted numerical settings, stage schedule, bias stiffness/direction,
sampling and geometric thresholds must match. The caller must document attachment
and structural-selection correspondence. Matching settings do not establish
force-field comparability or equivalent equilibrated initial ensembles.

Every declared replica runs sequentially through the same MD owner. Results retain
failed/rejected/cancelled runs. A candidate with an incomplete replica has no
successful-subset aggregate score. Complete candidates receive descriptive
between-replica means, sample standard deviations and standard errors, where
available. Repeated velocity seeds alone do not establish independent structural
exploration. A conditional Pareto set considers mean sampled peak tensile force
and final retained-contact fraction, only when both are present. These are not
calibrated uncertainty bounds or a ranking of experimentally proven materials.

## Completion gates still open

1. Actual complete-package Apple build, all native tests, public CLI execution,
   rejection restoration and prefix replay on a named physical GPU.
2. Independent force/energy and protocol-work comparison, including periodic
   groups and virtual solvent, stage/timestep refinement and temperature response.
3. Prepared real-protein reference panel with complete topology/parameters,
   calibrated sampling, held-out physical observations, and explicit failed cases.
   Paper-specific force fields and preparation cannot be silently substituted.
4. Automatic donor/acceptor typing and validated unfolding/event detection; native
   contact scoring is not by itself a thermal melting model.
5. Backbone/sequence generation with installed model weights, constraint handling,
   iterative proposal/evaluation and an independently scored improvement. This
   change invokes no RFdiffusion, ProteinMPNN or folding-model inference.
6. Calibrated molecular-to-material constitutive transfer, explicit network and
   crosslink architecture, concentration and solvent conditions, and NumiLab/Matter
   validation. No hydrogel stiffness or constitutive model is inferred from one
   protein's sampled peak force.

The numerical methods and workflow are inside NumiVivo. Macroscopic material
execution remains with NumiLab/Matter; this change neither duplicates that solver
nor gates molecular screening on completion of NumiBrain or the human simulator.
