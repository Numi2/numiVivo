# Shared-orbital ECC electronic paths

The numerical path is implemented by `VivoECCReactionPath`, molecular preparation
by `VivoMolecularECCPath`, and persistent execution by
`VivoMolecularECCPathWorkflow`. It extends the integrated ECC implementation; it
does not replace it with the older standalone QIO objective evaluator.

## Scientific contract

The input is an **ordered electronic profile**, not a claimed transition-state
search. Every geometry has the same ordered atoms, Gaussian basis recipe, alpha
and beta populations and frozen external point charges. The molecular route
requires closed-shell RHF preparation and records explicit atom identifiers.
Nuclear coordinates are in Bohr. The path coordinate is ordered, finite and has
an explicit unit. Atom identities are supplied and checked, not inferred from
coordinate proximity or an energy match. Bond connectivity can change along a
chemical path; atomic composition and electron number cannot silently change.

Preparation constructs a physical electronic Hamiltonian at every geometry and
integrates the adjacent Gaussian AO overlap matrices. With the RHF coefficients
C_i, the physical orbital overlap is C_i^T S_(i,i+1) C_(i+1). Equal orbital counts
are not evidence that this overlap is an identity matrix. Overlap spectra must
be contractions between orthonormal finite orbital spaces.

Within the declared transport groups, sequential orthogonal Procrustes alignment
sets a common orbital convention. One shared rotation U = exp(-K) then acts on
all aligned Hamiltonians. Its nonredundant parameters are restricted to the
existing ECC locality groups. Empty locality groups fix the shared rotation;
they do not implicitly allow unrestricted mixing.

This is a finite-basis implementation of the simultaneous path-orbital direction
proposed in the paper's outlook. The implemented objective and constraints are
explicit NumiVivo choices; this is not a claim to reproduce unpublished code.

## Nested execution

Every shared-rotation evaluation performs, at each geometry:

```text
physical Hamiltonian in aligned/shared frame
  -> correlated reference and orbital information
  -> fragment/bath/inactive-environment construction
  -> impurity solution
  -> complete selected-moment and number feedback when configured
  -> physical energy, particle closure and QIO objective
```

Single-fragment mode continues to omit property matching. Partition mode retains
its separate correlation potential at each geometry. The orbital rotation is
shared, not the correlation potentials or the Hamiltonians. Inner calculations
start from the deterministic initial potential rather than inheriting a
history-dependent potential from a different finite-difference trial.

Independent pointwise orbital optimization is disabled inside this outer loop.
The path objective is the weighted mean of each point's fragment-averaged ECC
QIO objective. Weights must be positive and sum to one. A bounded finite-difference
shared gradient and residual-decreasing line search update U. Acceptance requires
shared stationarity, settled point energies and converged inner calculations.
These are first-order local conditions, not a global-optimum guarantee.

The configured point-evaluation budget counts full point ECC evaluations,
including finite-difference and rejected trial evaluations. Each point also
retains the existing inner matching and numerical resource limits. Conservative
path allocations are checked before retaining all geometry results. This is not
an unbounded trajectory-scale implementation.

## Continuity is stronger than equal dimensions

The declared number of bath orbitals is fixed for every fragment across every
geometry and trial. Ranked environment candidate truncation is disabled in this
version; unsupported truncation is rejected rather than silently introducing
candidate discontinuities.

Adjacent fragment and bath projectors are compared with the **physical** orbital
overlap, not a coordinate-space identity. Their minimum singular values must
satisfy configured lower bounds. Same-geometry optimization trials are checked
against the accepted frame as well. A same-size but discontinuous subspace is
not accepted merely because its matrix dimensions match.

The correlated reference vectors also have exact determinant overlaps computed
using the existing fermionic orbital-frame machinery. A low overlap rejects a
reference-state discontinuity even when the bath dimensions have not changed.
This is a sampled ground-reference continuity guard, not an excited-state root
selector or proof that no crossing lies between sampled geometries. Refine the
path or introduce an explicit multistate model when the guard fails; do not
label a failed continuity test a successful calculation.

Stored results distinguish transport rotations, shared rotation, continuity
metrics, per-geometry ECC results, point-evaluation count and convergence state.
An empty bath has no bath singular-value diagnostic; it is represented as absent,
not as a fabricated measurement.

## Results and cache validation

The output contains the exact request fingerprint, all snapshot fingerprints,
path coordinates, relative electronic energies referenced to the first point,
per-point physical energy decompositions and a mandatory energy-meaning label:

`electronic-profile; no transition-state characterization, thermal correction,
standard-state correction or kinetic export`

The validator rebuilds physical orbital transport from Gaussian cross overlaps,
checks every inner ECC output, reconstructs the complete path, and recomputes the
shared orbital gradient. A collection of converged fixed-frame point outputs is
not enough to certify a different shared rotation. Continuity and relative energy
claims are checked numerically on fresh and reused outputs.

AO integrals and RHF references remain individual nodes using the original
operation identities. The final shared-path node depends on all of them. Changing
one geometry invalidates its two preparation nodes and the shared result, not the
other geometry nodes. Changing the objective weights invalidates the shared
result but leaves the ten preparation nodes of the five-point example reusable.
The original electronic commands use the same binary/backend fingerprint, so
an isolated-point calculation can reuse preparation produced by a path request.

Cache reuse protects reproducibility and detects content corruption. It is not
an authentication system for untrusted executables, and matching a reference
fixture does not establish the validity of an unrelated chemical model.

## Commands

On macOS, build the actual scoped production numerical module and command classes:

```bash
bash Tools/NativeChemistry/build_scoped_cli.sh /tmp/numivivo-path-cli
/tmp/numivivo-path-cli/numivivo-chemistry chemistry-path-template h2-stretch \
  --output h2.path.json
/tmp/numivivo-path-cli/numivivo-chemistry chemistry-path-run h2.path.json \
  --store .numivivo/chemistry-artifacts --output h2.profile.json
```

The complete product's main command router also dispatches these commands, but
its unrelated biology/Metal targets are not exercised by the scoped builder.
`chemistry-path-help` describes the command contract. With `--output`, the sibling
`<output>.receipt.json` records each task, cache reuse and result identity. The
result identity addresses the stored output envelope, not the exported payload.
Request-file aliases and output paths inside the store are rejected before an
overwrite. Unknown or duplicate options are errors, not alternate method choices.

## Reproducible checks

```bash
bash Tools/check_ecc_path.sh /tmp/numivivo-path-conformance
```

This builds the real Apple numerical module and command classes, runs 26 native
path checks, performs 20 mandatory independent PySCF comparisons, and runs 21 production CLI/store checks. PySCF 2.8.0 is an isolated conformance dependency only.
`VIVO_PATH_REFERENCE_DIR` can specify already-generated hashed reference fixtures;
`VIVO_ORACLE_PYTHON` can specify the oracle interpreter. Neither is required by
normal native path execution. Use an empty output directory for integration
checks; the runner does not remove previous artifacts.

The five H2 geometries span 1.2 to 2.0 Bohr in STO-3G, with a nonidentity initial
shared rotation. The calculation actually optimizes that rotation; it is not a
repeated fixed-orbital single-point example. Physical cross overlaps are compared
independently with PySCF. Full-bath ECC energies are compared to independent FCI
energies at every geometry. Signed SCF orbital rephasing must preserve the result.
Negative tests include missing/nonphysical overlaps, wrong bath size, same-size
subspace loss, reference discontinuity, budget exhaustion, forged energies and
shared frames, changed atom mapping, spin population and external charges.

The full-bath H2 limit qualifies the integration and numerical conventions. It
is **not** evidence of reduced-fragment reaction-barrier accuracy, a nontrivial
multi-fragment reaction path, or the paper's acrylamide-methanethiolate/BTK result.
The present path command does not combine smooth C-PCM with ECC, optimize nuclear
positions, locate saddle points, calculate thermal/standard-state corrections,
or install a kinetic rate. Those require their own coherent execution and
qualification rather than relabeling this electronic profile.

The `Shared ECC path conformance` Actions workflow retains exact repository
source, compiler/platform identity, hashes, reports and failed logs. The audit at
`Documentation/Audit/ecc-path-conformance-observations.json` records what actually
ran and distinguishes the Apple and portable local checks.

## Recorded Apple result

Code commit `1937512c3371e70a611f12221605abddd31a7343` passed all 67 new
checks in Actions run `33956068649`. The runner used Apple Swift 6.1.2 targeting
arm64 macOS 15 with the existing FP64 Accelerate matrix path. Both independent
oracle generation and Apple execution jobs succeeded. The prior native chemistry
conformance workflow also passed on the same code commit (run `33956068641`).

The five-point shared objective fell from -1.38213864366309 to
-4.036234335459631 and the final shared-gradient norm was 7.0568662025038975e-6.
The successful optimization reports 375 complete point evaluations; validation
and separate regressions perform additional evaluations. Lower QIO objective is not itself
proof of improved reaction-energy accuracy. Maximum total-energy difference
against the independent H2 FCI reference was 6.60e-10 Hartree, and the maximum
relative-profile difference was 7.89e-10 Hartree. The full-bath native ECC/FCI
comparison differed by at most 2.67e-15 Hartree. All compiled source identities
in the retained report were verified against the exact archived code commit.

## Primary method reference

CovAngelo, outlook on simultaneous orbital rotations across reaction-path points:
https://arxiv.org/html/2604.10487v1
