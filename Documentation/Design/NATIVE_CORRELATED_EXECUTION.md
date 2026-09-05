# Native correlated electronic execution

Status: source integrated into `NumiVivoKit` and the `numivivo` command router.
The scoped FP64 numerical runner passed on Linux. The complete Apple package,
Accelerate execution, artifact-store execution and command-line workflow were
not run in that environment. This document does not certify the whole project.

## One execution path

`VivoElectronicWorkflowRequest` describes an electronic system, Gaussian basis,
method and resource budget. Coordinates are in **Bohr**, energies in **Hartree**.
This is an electronic boundary; `VivoMolecularStructure` retains its existing
nanometre convention. No input is converted by guessing its numerical magnitude.

The correlated molecular path is:

```text
system + basis artifacts
  -> AO integrals
  -> converged RHF
  -> real spatial embedded Hamiltonian
  -> selected MP2 / CI / CCSD / CAS solver
  -> validated, content-addressed result
```

These are separate `VivoChemistryWorkflow` DAG operations. Changing the final
solver does not change the integral or RHF task identity. Reuse still requires
the same implementation, configuration, input hashes, output contracts and
resource settings. Output validation is performed on cache hits as well as new
execution. The workflow does not silently replace a failed requested method.

The application identifies its implementation using the executing binary's
SHA-256, OS version, architecture and FP64 backend class. A human-readable version
label alone is insufficient. These identities establish reproducibility and
integrity, not authorship, scientific correctness or trust in arbitrary binaries.
The cache is an application-owned local store, not an authentication system.

## Commands

From a complete checkout on a supported Apple development machine:

```bash
swift build -c release
.build/release/numivivo chemistry-template h2-ccsd --output h2.json
.build/release/numivivo chemistry-run h2.json \
  --store .numivivo/chemistry-artifacts --output h2.result.json
```

The second invocation with unchanged inputs can reuse verified stages:

```bash
.build/release/numivivo chemistry-run h2.json \
  --store .numivivo/chemistry-artifacts --output h2.repeated.json
```

These full-package commands are supplied for Apple-side execution; they were not
executed by the Linux development runner. `chemistry-help` describes all flags.

`chemistry-template` supports `h2-ccsd`, `h2-casci`, `h2-casscf`, `h2-lda`,
`h2-lda-cpcm`, `h2-cpcm-rhf`, `solver-ccsd` and `solver-fci`. The H2 example is
STO-3G at a bond length of 1.4 Bohr. Its all-orbital CAS example has no nonredundant
orbital rotations; the four-orbital regression below exercises actual CASSCF
orbital optimization.

For a prepared `VivoEmbeddedHamiltonian`:

```bash
.build/release/numivivo chemistry-template solver-fci --output solver.json
.build/release/numivivo chemistry-solve hamiltonian.json \
  --solver solver.json --output many-body.result.json
```

`chemistry-solve` also accepts `--budget budget.json`. Both execution commands
accept `--store`; its default is `.numivivo/chemistry-artifacts`. With `--output`,
a sibling `<output>.receipt.json` records task, input, implementation and result
fingerprints, with per-stage reuse flags. The result fingerprint addresses the
stored workflow output envelope, not the raw exported payload bytes. Output and
receipt paths may not alias
an input request, solver or budget file. Unknown/duplicate options are rejected.

The molecular correlated preparation route currently requires closed-shell RHF.
Direct embedded-Hamiltonian CC/CI calculations can specify other fixed alpha and
beta populations. This is not an unrestricted-orbital AO-to-CC preparation path.

## Projective CCSD and response properties

`VivoCoupledCluster` constructs the fixed-alpha/fixed-beta determinant sector and
uses it to evaluate `exp(-T) H exp(T)`. It does **not** diagonalize that matrix to
obtain the CCSD energy. `T` contains only single and double excitations relative
to the selected reference. Exponentials terminate at their nilpotence order.

Excitations use interleaved spin-orbital indices and explicit fermionic phases.
Each operator is phase-normalized so that its action on the reference yields a
positive canonical determinant. The serialized amplitude vector therefore has
an explicit excitation list; it is not an undocumented conventional `t2` tensor.

The projected residual and analytic commutator Jacobian drive bounded Newton
steps with a residual-decreasing line search. Iteration limits, singular
Jacobians and failed line searches are reported explicitly. A converged root is
not proof that CCSD found the global ground state or that a single-reference
approximation is appropriate for that molecule.

The optional Lambda solve supplies biorthogonal left/right states and unrelaxed
one- and two-particle response density matrices. The implementation checks the
Lambda residual, left/right overlap, particle trace and 2-RDM contraction.
Normalized right-wavefunction diagnostics remain separate from projective CC
energy and response properties.

**CC response matrices are not positive density operators for orbital entropy.**
`VivoManyBodySolverResult.informationState` is intentionally nil for CCSD and MP2.
Explicit CI/CAS states are available for the existing information analysis.

This implementation stores a full determinant-sector Hamiltonian and is a
bounded small-system conformance engine. It is not a tensor-contraction CCSD
implementation suitable for large molecular bases. Determinant, memory and
aggregate operator-application limits reject oversized work before unbounded
allocation. No GPU acceleration of this determinant-sector solver is claimed.

## CASCI and CASSCF

`VivoActiveSpace` records doubly occupied core, active and explicitly frozen
orbital columns. CASCI preserves this basis. CASSCF optimizes nonredundant
core/active/external rotations, excluding frozen columns and rotations within
one subspace.

The orbital gradient is analytic: it reconstructs the full-space density from
the active CI cumulant, occupied core and empty external orbitals. Its dense
contractions use the common FP64 matrix path. Transported limited-memory BFGS,
curvature checks, bounded rotations and an Armijo search control optimization.
Every candidate transforms the original integrals to avoid repeatedly rotating
already-rounded four-index tensors.

This is ground-state, fixed-alpha/fixed-beta CAS optimization, not state-averaged,
root-followed or explicitly total-spin-adapted CASSCF. A near-degenerate active
root can reject and request an explicit root policy. Active-space CI remains
bounded. Neither DMRG nor a scalable multireference solver is implied.

## Fragment number matching

`VivoFragmentNumberMatcher` solves a shared chemical potential against the sum
of fragment populations. Each cluster receives `-mu*N_fragment`, with a checked
orthogonal fragment projector. Applying a uniform shift to all cluster orbitals
cannot change a fixed-sector total electron count and is not used as a substitute.

Bracketing, monotonicity checks and bisection are bounded. Acceptance requires
the population residual, not merely a small chemical-potential interval. Plateaus,
unreachable targets and discontinuities do not become successful matches.
Physical cluster energies remove the artificial chemical-potential term.
Overlapping cluster energies must not be summed as a DMET total energy.

This is a multi-fragment building block. The single-fragment formulation in the
CovAngelo paper does not require this property-matching loop. This addition is
**not** a complete ECC-DMET correlation-potential self-consistency cycle and does
not itself qualify path-consistent QIO reaction energies.

## LDA and C-PCM consolidation

The restricted PZ81 LDA SCF now optionally includes the C-PCM reaction field.
Geometry-dependent AO grid values are cached for density/XC contractions;
`VivoCPCMOperator` caches surface AO integrals and nuclear potentials. The same
operator is used by restricted HF and LDA rather than rebuilding it on every
SCF iteration. Dense contractions use Accelerate FP64 on Apple and the existing
portable FP64 implementation otherwise.

The surface equation remains matrix-free PCG with a memory/work contract,
positive-density/electron-count checks, a recomputed **true** equation residual
and the existing Gauss-law diagnostic. Final canonicalization rechecks energy,
density, SCF residual and LDA grid electron count. Initial iteration traces no
longer serialize infinity as an absent previous-density comparison.

The cavity is the existing point-charge tessellation, not a newly validated
high-order smooth cavity. The implemented functional is spin-unpolarized LDA,
not B3LYP or a range-separated functional from the paper. C-PCM plus external
QM/MM point charges still rejects without a combined-cavity/outlying-charge
policy. Analytic nuclear solvent derivatives are not supplied here.

Electronic energy plus electrostatic polarization is **not a complete Gibbs
free energy**. Thermal, zero-point, standard-state and non-electrostatic solvent
contributions are absent. Such output must not be relabeled as an activation
Gibbs energy and installed as a biological rate.

## Scoped reproducible checks

```bash
bash Tools/check_correlated_chemistry.sh /tmp/numivivo-correlated-check
```

This compiles selected production numerical sources with Swift 6 and runs 29
checks without Python, CUDA, provider SDKs or artifact-store test substitutes.
It writes compiler/platform records, exact input-source SHA-256 hashes, numerical
results and reusable synthetic Hamiltonian/solver JSON files into the requested
directory. On Apple the shared dense algebra can use Accelerate; the recorded
development observation used portable FP64 on Linux, not that Apple branch.

Checks cover H2 CCSD/FCI agreement; a four-electron case in which they differ;
Lambda RDM energies and an independent finite-difference response derivative;
five analytic CAS orbital derivatives; nontrivial CASSCF optimization; budget,
reference, method and nonconvergence rejection; fragment matching and fixed-N
invariance; finite JSON traces; and coupled LDA/C-PCM convergence/accounting.

The numerical observations in `Documentation/Audit/native-correlated-observations.json`
refer to the selected local source snapshots listed there, not a complete build
of a repository commit. Source hashes identify exactly the snapshots that ran; they do not claim every
pre-existing file was byte-identical to current `main`. The full CLI/store integration was separately Swift-6
typechecked against inspected interface declarations; those interface stand-ins
were not executed and are not shipped in this repository.

No independent PySCF/ORCA conformance campaign, Metal execution, full Apple build,
acrylamide-methanethiolate reaction or BTK reproduction is claimed. The next
scientific qualification must include external small-molecule comparisons and
the actual paper geometries/protocol, not only these algebraic checks.

## Method references

- CovAngelo methodology and single-fragment formulation: arXiv:2604.10487,
  https://arxiv.org/abs/2604.10487.
- Public coupled-cluster reference interface: https://pyscf.org/user/cc.html.
- Public multiconfigurational reference interface: https://pyscf.org/user/mcscf.html.

These are methodology/conformance references. Their inclusion is not a claim
that an external package comparison was executed in this development increment.
