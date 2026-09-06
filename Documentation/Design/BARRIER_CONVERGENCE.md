# Fixed-geometry barrier convergence

`VivoBarrierConvergence` runs a declared sequence of reduced electronic problems
against full-space FCI at exactly the same mapped nuclear geometries. It extends
the existing native reaction workflow without adding another runtime or store.
Production execution is Swift/FP64 with the existing Accelerate contractions on
Apple; Python/PySCF is used only by the independent conformance adapter.

## What is and is not assessed

The report measures forward, reverse and reaction electronic-energy differences
using the declared interior evaluation point and the first/last input structures.
It also measures the entire relative electronic profile and absolute total-energy
errors. A large common offset cannot hide profile error; good total energies do
not by themselves qualify a barrier difference.

Calling a point `barrierPointIdentifier` is not stationary-point characterization.
The campaign does not optimize nuclei, reconstruct an IRC, assign endpoints,
apply RRHO corrections, or create a rate. The existing nuclear/connectivity tools
must qualify those claims separately on the same model and input geometry.
First/last supplied structures are not automatically isolated reactants or
products at infinite separation. No paper result is inferred from these inputs.

All nuclei retain explicit atom mapping, basis primitives, electron populations,
and frozen external point charges. Nuclear coordinates use Bohr. Changing basis,
composition or spin within one profile is rejected. The complete FCI reference
sector is capacity-checked before expensive AO integral generation.

## Common physical orbital frame

Default preparation diagonalizes the one-electron core Hamiltonian in a complete
orthonormal AO frame. This is not an RHF claim. Adjacent Gaussian cross overlaps
are integrated natively, and group-wise Procrustes transport proceeds in both
directions from the declared anchor. An explicit complete anchor coefficient
matrix may replace the default. Physical subspace loss fails the calculation;
equal orbital dimensions are insufficient. Exact reference-CI overlaps are
checked at every path edge to detect sampled ground-state discontinuities.

The optional `ensembleOrbitals` policy forms a weighted full-reference 1-RDM in
that transported common frame. Its weights are declared inputs. One ordered
natural-orbital rotation is applied to every Hamiltonian, reference wavefunction,
AO coefficient frame and cross overlap. It is not separately optimized at each
geometry and is not fitted to reference barrier errors. The occupation spectrum,
particle trace and applied rotation are retained. A reduced CAS boundary that
splits a near-degenerate occupation subspace is rejected.

This reference-assisted policy requires full CI reference data. It does **not**
demonstrate a computational speedup, a scalable active-space selector, or
prediction without the reference. It is an explicit conformance/development tool
for finding which reduced models retain a given energy profile.

## Approximation levels

A request uses one family throughout its hierarchy:

- CASCI: strictly nested active spaces, occupied core and empty external orbitals.
  Core orbitals may be promoted into the active space, not demoted into the empty
  external space. True CI residuals are checked after every solve.
- ECC: the actual integrated `VivoECCReactionPath`, including physical bath
  continuity and configured inner matching, with increasing declared cluster
  sizes. CASCI is not used as an undocumented substitute for failed embedding.
  The current ECC implementation still requires its supported closed-shell
  sector. Inner reference/matching/objective and fragment identities are bound.

Resource or numerical failures remain explicit per-level outcomes. Later
full-bath success does not remove or conceal a failed earlier level. Malformed
requests fail immediately rather than becoming scientific observations.

The point-evaluation contract charges reference and CAS point solves, and reserves
both ECC execution and reconstruction limits before entry. Nested matching and
operator contractions keep their own limits. This is not a claim that the report
counts every matrix product or every nested solver iteration. Cache validation
is a separate bounded reconstruction run.

## Acceptance and artifact semantics

The default accuracy targets are 0.001 Hartree for maximum barrier/reaction error,
maximum relative-profile error and change between consecutive levels. The input
must declare at least two stable genuinely reduced levels. A full-space endpoint
never counts as one of those reduced levels.

`reducedAccuracyEstablished` requires that the final consecutive reduced window
meets both reference-error criteria and successive-level stability, and that all
later larger levels also agree. `reducedAccuracyNotEstablished` means the levels
ran but failed those scientific criteria. `incompleteLevelExecution` preserves
one or more numerical/resource/subspace failures. Accepted reduced identity is
absent for both negative outcomes.

A successful `reaction-run` can export a **negative assessment**. That is a valid
numerical observation, not permission to use an inaccurate model. Consumers
must inspect `assessment` and `acceptedReducedLevelIdentifier`. No kinetic rate
is installed. The existing artifact store binds input, implementation, operation
and resource identities and verifies SHA-256. Output validation re-executes the
campaign and compares reference values, orbital frames, individual outcomes,
error metrics and final assessment. It cannot promote a forged success flag or
silently remove an unfavorable level.

## Reproducible molecular inputs

Two committed examples describe nine collinear H + H2 exchange-scan structures
in explicit hydrogen 6-31G Cartesian Gaussian functions (six spatial orbitals,
three electrons). The nuclear scan is declared, not optimized. It is **not** the
previously characterized H3/STO-3G saddle and does not inherit that qualification.
The CAS levels contain three, four, five and six active spatial orbitals.

```bash
swift build -c release
.build/release/numivivo reaction-template h3-barrier-convergence --output plain.json
.build/release/numivivo reaction-template h3-barrier-convergence-ensemble --output ensemble.json
.build/release/numivivo reaction-run ensemble.json --store .numivivo/chemistry-artifacts --output ensemble.result.json
```

Equivalent committed inputs are `Examples/reaction/barrier-h3-631g.json` and
`Examples/reaction/barrier-h3-631g-ensemble.json`. The latter adds uniformly
weighted common-reference natural orbitals; every other scientific criterion
remains unchanged. Neither example represents acrylamide-methanethiolate or BTK.

## Independent comparisons and regression tests

`ReferenceAdapters/PySCF/barrier_convergence_oracle.py export` reads only declared
inputs and reconstructs normalized Gaussian integrals, physical overlaps, orbital
transport, reference FCI and every CASCI level with pinned PySCF 2.8.0. The adapter
checks CAS, not ECC via a CAS replacement. Its output is bound to the exact input
bytes and hashed independently. It never reads native energies during export.

`check` uses only the Python standard library and the actual production executable.
It verifies template/input agreement, all point and level energies, barrier
references, orbital/wavefunction overlaps, optional natural occupations, negative
accuracy assessments, cache identity and invalid/corrupted input rejection.

`BarrierConvergenceTests` also exercises the actual ECC hierarchy, retains
fractional inactive-population failures, rejects forged reports and tests the
positive assessment logic on an explicitly algebraic repeated-geometry control.
That positive control is not chemical evidence.

The read-only `Barrier convergence conformance` workflow builds the complete
Apple release product and all test targets, runs the full Swift suite (including
existing real Metal and mapped connectivity checks), and compares both production
campaigns to independently generated references. It never rewrites source during
validation. Reports, logs, fixtures, exact source and compiler identities are
retained even when a job fails. Consult the measured audit rather than inferring
success from this implementation description.

## Method references

- CovAngelo's path-consistency motivation: https://arxiv.org/abs/2604.10487
- Independent full-CI reference implementation: https://pyscf.org/user/ci.html
- CAS reference conventions: https://pyscf.org/user/mcscf.html

These describe methods/reference interfaces, not evidence that a particular
NumiVivo campaign is accurate. The executable result supplies that assessment.
