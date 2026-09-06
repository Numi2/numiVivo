# Reusable MD → QM/MM interface

This interface prepares a finite, electrostatically embedded QM/MM sample from
an existing MD snapshot or accepted checkpoint. It is general platform code:
there are no paper-specific atoms, reactions, targets or reference energies.
The authoritative inputs remain `VivoMolecularStructureDocument`,
`VivoClassicalSystem`, and the sampled MD state.

## Entry points

```swift
let frame = try VivoQMMMCompiler.prepareMD(
    document: document,
    system: classicalSystem,
    checkpoint: checkpoint, // snapshot: also supported
    request: baseRegion,     // explicit QM atom IDs and electron counts, including caps
    solventPolicy: .init(radiusNM: 0.5),
    configuration: .init(maximumImagePairs: 50_000_000),
    budget: budget
)

// A backend must return energy plus forces on every nucleus AND embedding charge.
let electronic = try electronicBackend(frame.region.electronicSystem)
let forces = try VivoQMMMForceMapper.assemble(
    document: document, system: classicalSystem,
    frame: frame, electronic: electronic, budget: budget
)

// Required before using a new backend/partition as a force-capable interface.
let check = try VivoQMMMForceChecks.check(
    document: document, system: classicalSystem, frame: frame,
    electronic: electronic, budget: budget,
    electronicEnergy: electronicEnergyBackend
)
guard check.passed, check.completePhysicalAtomCheck else {
    throw VivoChemistryError.convergence("QM/MM force qualification failed or was partial")
}
```

The two backend functions above are supplied by the caller; they are not newly
invented solver names. The force backend returns
`VivoQMMMElectronicForceResult`; the independent energy callback returns Hartree.
`VivoQMMMElectronicFiniteDifferences.evaluate` is a bounded reference adapter for
energy-only backends. It differentiates all nuclear and embedding-charge centers
and labels its derivative method explicitly. It is not an analytic-gradient
implementation. Electronic response, basis/Pulay derivatives and solver
convergence remain the backend's responsibility.

The original `prepare` and selection-only `promoteSolventRegion` APIs remain.
The latter does **not** unwrap its input. Use `prepareMD` to connect selection,
reconstruction, promotion, link placement and identity mapping in one operation.
Do not set `coordinatesAreUnwrappedFiniteCluster` on wrapped MD coordinates to
bypass preparation. The new path sets it only after reconstruction succeeds.

## Coordinate and molecular contract

The sampled cell is authoritative, including `nil` for a finite sample. The
imported structure cell cannot silently replace a changed NPT cell. Preparation
is host-side FP64 work at an explicit sample boundary; it does not mutate the
GPU MD state, velocities, clock or integrator.

A shared topology graph combines chemical bonds, classical bonds, distance
constraints and supported virtual-site parent dependencies. It establishes a
physical atom↔particle bijection and complete connected molecular components.
Residue labels select solvent identities but never define molecular closure.
Bond/constraint traversal uses the existing bounded closest-lattice-image MD
helper, including skewed triclinic cells. Inconsistent cycles, periodic winding
and ill-conditioned cells outside that helper's search capacity are errors.
Missing connectivity is not repaired by distance guessing.

The component containing the smallest base QM atom is anchored at that atom's
sampled position. Other QM components and each remaining molecule receive a
single whole-molecule lattice translation chosen deterministically near the
base QM region. Within a connected molecule, physical atoms are unwrapped by
connectivity. Solvent heavy/all-atom distance rules are shared with imaging.
There is no second independent minimum-image selection after placement:
promotion uses Euclidean distances in the **exact finite coordinates** passed
to the electronic compiler.

Each promoted solvent molecule is complete, even when represented by multiple
residue records. A molecule is never shortened to meet an atom/residue cap.
Partial base ownership of a solvent molecule, incomplete solvent connectivity,
or a solvent label connected to non-solvent atoms without explicit whole-molecule
selection is rejected. Promotion assumes closed-shell added molecules and uses
formal charges to update both electron sectors. Even electron count does not
establish a physical spin state. Per-residue `minimumDistanceNM` records the
minimum distance of its enclosing selected molecule; residue electron counts
sum to the molecular addition.

All original particle slots survive: physical atoms, MM boundary atoms, and
supported linear virtual sites. Sites are checked against unwrapped parents
modulo the sampled cell and then reconstructed from their parent weights in the
final geometry. Stale sites are not silently repaired. Supported negative
weights are retained. QM-owned site charges are excluded from embedding; their
existing cross-partition LJ parameters are retained. Sites with mixed QM/MM
parents, unmodeled site rules, nested/nonlinear site constructions, Drude
particles, and unresolved Z1 site ownership fail explicitly.

## Mapping and artifacts

`VivoMDQMMMPreparedFrame` retains the source step/time/configuration, optional
accepted-checkpoint payload digest, finite cluster, promotion decision, final
request and prepared region. Its cluster contains:

- Original `atomToParticle` and per-particle atom/molecule identities, with no reindexing.
- Integer lattice images **added** to physical source coordinates; virtual-site images are `nil` because sites are reconstructed.
- Complete molecular membership, site construction rules, sampled cell and source/system fingerprints.

For a physical particle, `rFinite = rSource + image[0]*a + image[1]*b + image[2]*c`.
The anchor remains in the source frame; no artificial recentering is introduced.
`validate` reconstructs the region from the finite geometry and checks mappings.
`validateSource(..., snapshot:)` additionally audits lattice translations against
the original immutable sample. A digest or structurally valid mapping alone is
not a replacement for retaining that source. Workflow output validation replays
the complete preparation from the original input payloads.

Nuclear `structureAtomIndex` and charge `classicalParticleIndex` mappings remain
in `VivoElectronicSystem`. Artificial hydrogen nuclei have explicit link records,
not invented structure atom IDs. Electronic force results are bound to the
canonical electronic-system digest and require every center in the same order.
Missing MM reactions, mismatched coordinates, nonfinite values or unsupported
energy conventions are errors, not zero-filled arrays.

## Forces and energy accounting

Internal electronic quantities are Hartree and Hartree/Bohr. Forces mean
**minus** energy derivatives. MD-facing forces are also available as
`particleForcesKJPerMolPerNM`, using the shared atomic-unit constants.

For the compiler's fixed-distance cap,

```text
u = (rM - rQ) / |rM - rQ|
rL = rQ + d*u
JM = (d/r) * (I - u*uᵀ)
JQ = I - JM
FM += JMᵀ * FL
FQ += JQᵀ * FL
```

A constant-ratio split is incorrect for this geometry. The implemented Jacobian
preserves cap force and torque; several caps on one atom accumulate. Radial cap
force goes to the QM endpoint. Direct real-nucleus forces and all MM embedding
reactions are added by source identity. Cross-partition LJ returns equal and
opposite raw particle forces, including actual virtual-site centers. Artificial
link atoms have no LJ parameters. Original mixing rules, pair tables and LJ
exception scales/overrides are retained.

Only after all these contributions are combined are linear site forces
redistributed once as `Fparent += weight * Fsite`. Virtual slots in the final
result are zero. The legacy LJ-only QM view is already redistributed; it must
not be added to the new full-particle raw LJ array. Derivative method metadata,
source/cluster/region identities, net force and centroid-referenced net torque
are retained in the assembled result and/or check report.

The assembled energy is the existing finite-cluster Z1 electronic convention
plus cross QM–MM LJ. It is **not the complete MD Hamiltonian**: MM self energy,
retained classical bonded terms and their forces are not included here. Do not
add this result blindly to a full classical force field, which would double
count replaced QM and QM–MM interactions. A dynamical caller needs an explicit
retained/replaced classical-term partition before integrating these forces.
The interface does not implement periodic Ewald/PME QM/MM, cell stress/virial,
adaptive-region dynamics, free-energy corrections or polarizable embedding.
A complete molecular finite cluster is not a periodic electronic calculation.

## Force checks

`VivoQMMMForceChecks` compares mapped physical forces with total finite-cluster
energy differences at steps `h` and `h/2`. Each physical displacement rebuilds
virtual sites, link nuclei and charge coordinates; membership and lattice images
remain fixed. Links are never displaced independently. The check combines
componentwise absolute/relative error, step-halving stability, agreement with
the supplied baseline energy, total force and total torque. An explicitly
sampled atom subset is labeled partial through `completePhysicalAtomCheck`.
Energy-call limits and pair/allocation budgets fail closed. Reports cannot
establish the chemical accuracy of the underlying electronic model.

## Workflow and CLI integration

The normal `VivoPlatformOperations.registry` exposes:

| Operation | Inputs | Outputs |
|---|---|---|
| `vivo.platform.md-qmmm` | original structure document, classical system, accepted checkpoint | prepared `frame`, electronic `system`, complete cluster `mapping` |
| `vivo.platform.qmmm-forces` | original structure document, classical system, prepared frame, electronic-center force result | mapped interaction `forces` |

The preparation configuration is `VivoWorkflowMDQMMMConfiguration` (explicit
base request, optional solvent policy, cluster limits). Force assembly has empty
configuration. The preparation `system` output has the standard
`vivo.electronic-system` kind and feeds existing integral/electronic workflow
nodes without a file-format conversion. The force input kind is
`vivo.qmmm-electronic-forces`; availability of an energy solver alone does not
claim an analytic force producer. Existing `workflow-plan`/`workflow-run` and
artifact-store validation apply; no special paper-oriented command is required.
The isolated `md-snapshot` adapter remains unchanged and still retains
periodicity rather than silently converting periodic input to isolated QM.

## Regression coverage and execution

`MDQMMMInterfaceTests` covers wrapped water and virtual sites, permuted particle
order, sampled versus imported cells, skewed closest images, nonperiodic
samples, whole molecules spanning residues, capacity failures, stale sites,
periodic winding, checkpoint lineage/source audits, missing MM reactions,
wrong forces, multiple fixed-distance caps, and two-step energy/force agreement.
`PlatformQMMMWorkflowTests` executes both registered operations and checks output
reconstruction and tamper rejection. The existing `SolventPromotionTests`
retains selection-only API regression coverage.

```sh
swift test -c release --filter 'MDQMMMInterfaceTests|PlatformQMMMWorkflowTests|SolventPromotionTests'
```

The general platform Apple workflow builds the unchanged full package and runs
these suites alongside existing platform/Metal tests. It requires each named
suite to execute and preserves the source, compiler, binary digest, logs and
native fixture artifacts. These tests include synthetic center potentials to
isolate coordinate/force mapping. They are not reproductions of an experimental
system and are not chemical accuracy benchmarks. Portable checks with support
types do not qualify native artifact hashing, the full package, Metal execution
or the workflow registry; consult the actual native run before claiming those.
