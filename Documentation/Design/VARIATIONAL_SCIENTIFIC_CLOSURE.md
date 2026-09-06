# Variational scientific closure

This document defines the bounded scientific capabilities added by the shared
residual barrier campaign and the coherent global fragment embedding. It is a
method and contract description; measured outcomes belong in
`Documentation/Audit/scientific-closure-observations.json`.

## Scope

The implementation has two separate variational families:

- `VivoResidualBarrierCampaign` enriches one nested global Fock space with
  Hamiltonian-action residual directions collected from every supplied
  geometry.
- `VivoVariationalEmbedding` unions overlapping fragment Fock spaces at the
  wavefunction level, removes duplicate directions, and solves one normalized
  CI state whose density is used by the smooth C-PCM feedback loop.

Neither family replaces the original orbital CAS campaign or the original
ECC-DMET energy formulation. The new global density is a distinct
variational CI/C-PCM formulation. It must not be described as democratic
ECC-DMET/PCM self-consistency.

## Residual-enriched reduced space

For each supplied geometry, the campaign constructs the existing native
electronic Hamiltonian and first solves the declared seed space. Given a
normalized projected CI state `c`, omitted Hamiltonian action is represented by

`r = (I - P P^T) H c`,

where `P` is the current orthonormal Fock-space basis. Independent residual
directions from all geometries are combined, orthonormalized, and appended to
one shared space. Every geometry is then resolved in that same enlarged space.
The sequence is nested and no fitted energy correction, perturbative shift, or
full-reference eigenvector is inserted into the approximation.

The eigenproblem is reduced, but the determinant sector and spatial orbital
modes are not removed. The implementation retains the complete determinant
vector for the declared sector and all six spatial orbital modes in the H3
6-31G example. This is not a smaller-qubit representation, a demonstrated
speedup, or a scalable protein-level FCI method.

The unchanged acceptance contract is:

- maximum barrier/reaction error: `0.001` Hartree;
- maximum relative-profile error: `0.001` Hartree;
- maximum change between consecutive levels: `0.001` Hartree;
- at least two consecutive genuinely reduced levels must satisfy the contract.

The accepted window is the final consecutive reduced window. A level may meet
the reference-error tests and still fail the successive-change test. In the
current H3 residual campaign, dimensions 27 and 36 meet the reference-error
tests, but the 27-level change from dimension 18 is `0.005267459604509472`
Hartree and therefore does not pass the stability threshold. The accepted
window is dimensions 54 and 63. The predecessor orbital-only campaign remains
a separate negative assessment and is not relabeled.

Every result retains the baseline assessment, all unsuccessful levels, physical
state-overlap checks, projected residuals, full residuals, determinant-sector
size, operation budget and the accepted reduced-level identity. Validation
re-executes the bounded campaign and rejects altered energies, vectors,
assessment flags or removed failures.

## Coherent global density and C-PCM

Each fragment contributes a Fock subspace. The implementation forms the
orthonormal union, removes duplicate columns, and diagonalizes one global CI
Hamiltonian in that union. The density is reconstructed from the normalized
global CI state, including cross terms between amplitudes from overlapping
fragment spaces. It is not a sum of independently normalized impurity density
matrices.

The result checks electron count, coefficient normalization, occupation bounds,
energy reconstruction and projected/full residuals. The external residual is
reported explicitly: small projected residual does not imply full-sector
convergence.

For correlated solvent feedback, the global density drives the smooth C-PCM
reaction field. Iteration acceptance requires density, reaction-potential and
energy convergence, plus stationarity of the returned correlated state in the
returned field. The electron-solvent contribution is included once through the
declared polarization-energy convention.

Existing ECC fragment/bath projectors can be imported only when inactive
occupied columns are declared explicitly. Fractional occupations are not
rounded into invented inactive determinants. The imported ECC frame is a source
of a variational space, not proof that the original ECC-DMET energy functional
has become globally self-consistent.

## Workflow and provenance

The reaction workflow exposes `residualBarrier`, `globalEmbedding`,
`connectivity` and `reproductionPreflight`. Production artifacts bind the
request, implementation identity, result, method label, density, projectors and
acceptance state. Cache reuse reconstructs and validates the numerical result;
changing stored scientific values is rejected.

The paper preflight checks declared content hashes, mapped structures, electron
populations, method/reference identities, and the required BTK snapshot,
topology and environment roles. A negative readiness report is an intentional
result. It does not authenticate author provenance and it does not perform a
paper reproduction.

## Required interpretation

The supplied H3/6-31G geometries are a fixed electronic scan, not an optimized
saddle or a paper reaction path. The global LiH example establishes a coherent
variational density and solvent feedback within its selected space while its
external residual remains nonzero. Neither result supplies a Gibbs barrier,
endpoint connectivity, transmission coefficient, kinetic rate, large-system
scaling result, or acrylamide-methanethiolate/BTK reproduction.
