# Integrated ECC-DMET execution and external conformance

Status: native finite-basis implementation, integrated with the electronic CLI
and content-addressed workflow. This document specifies the implemented model,
not a claim of reaction-barrier accuracy or full CovAngelo reproduction.

## Two distinct modes

`VivoECCDMETConfiguration.mode` is authoritative.

- `singleFragment` follows the paper's distinction: one embedded fragment,
  no property-matching potential and no artificial global-number iteration.
- `selfConsistentPartition` is the multi-fragment extension. Fragments must form
  a complete, disjoint partition of the supplied orthonormal orbital basis.
  A finite, explicitly declared set of fragment-local one-/two-body moments and
  the total fragment electron count are matched self-consistently.

Section 2.3.4 of arXiv:2604.10487 describes selected one-/two-particle potentials
and explicitly omits property matching from its implemented single-fragment
scheme. NumiVivo does not rename a chemical-potential-only loop as ECC-DMET.
The multi-fragment reconstruction convention below is a specified NumiVivo
extension, not a claim that unpublished proprietary code has been reproduced.

## Execution and numerical authority

`VivoECCDMET.solve` accepts a real spatial `VivoEmbeddedHamiltonian`, immutable
configuration, optional restart orbital rotation/potential, and resource budget.
The system uses chemist ERIs, interleaved alpha/beta spin ordering and Hartree
energies. Molecular preparation still converts from the existing structure units
at its boundary; this solver never guesses coordinate units.

Each orbital frame performs the complete sequence:

```text
physical Hamiltonian in the current orbital frame
  -> correlated reference with optional auxiliary moment potentials
  -> normalized reference state, RDMs and orbital information
  -> ranked candidates and rank-checked correlated-density bath
  -> inactive-environment density and projected physical cluster Hamiltonian
  -> impurity solution and optional global fragment-number matching
  -> selected moment residuals and inner feedback
  -> QIO objective, physical energy and particle closure
```

The current correlated reference is FCI or CISD; impurity states are FCI.
Direct FCI uses the existing matrix-free Davidson solver, not a stored square
determinant Hamiltonian. Explicit correlated reference vectors/RDMs are still
bounded and can be expensive. This is not DMRG or a scalable whole-protein
correlated-reference calculation. CC response RDMs are not accepted as positive
density operators for the entropy objective.

### Bath construction and inactive environment

The environment-fragment one-RDM block supplies left singular vectors through
the existing `VivoBathBuilder`. Its spectral bath cannot exceed the rank of that
block. Mutual information, orbital entropy and a declared cumulant score rank
candidate environment orbitals. Candidate truncation records and bounds omitted
one-RDM coupling. A declared minimum bath size can add independent ranked
vectors; the result reports these as **supplemental** rather than as SVD rank.

`VivoCIOrbitalFrame` evaluates the exact alpha/beta exterior-power transformation
of the reference vector into the complete fragment/bath/environment frame.
It retains interleaved-spin fermionic signs, checks normalization and budgets,
and performs no eigensolve. This enables physical reference RDM and information
analysis in exactly the same orbital frame as the projected Hamiltonian.

For cluster columns C, Q = I - C C^T and the inactive spin-summed density is
D_env = Q D_ref Q. The environment contribution is J[D_env] - K[D_env]/2.
The current Hamiltonian boundary requires a spin-unpolarized environment;
spin-polarized reference densities reject rather than losing their distinction.
The physical cluster one-body matrix is C^T (h + V_env) C, and its two-body
integrals are transformed in the same frame. Its scalar is zero: the physical
constant is restored once in the global energy expression.

### Self-consistency, QIO and failure states

In partition mode, real Hartree multipliers modify the **reference** Hamiltonian
using the declared fragment-local operators. The one-body potential is traceless
to remove its uniform gauge. Selected operator moments from the correlated
reference are compared to impurity moments; a bounded finite-difference
Levenberg step and line search update the multipliers. Every trial rebuilds the
reference, baths, inactive density, impurity Hamiltonians and states.

A separate shared chemical potential shifts only the fragment number projectors
inside impurity Hamiltonians. It does not uniformly shift the whole fixed-N
cluster. Number matching accepts the population residual, not merely a small
chemical-potential bracket. Auxiliary shifts are removed from reported physical
cluster energies.

The outer QIO loop optimizes rotations U = exp(-K) within declared locality
classes. An empty class list fixes all orbitals; it does not enable unrestricted
rotations. Each finite-difference derivative and line-search trial runs the inner
feedback anew. Bath dimension changes invalidate a derivative or reject a trial;
the code does not compare objectives from silently different active dimensions.

Success requires the configured orbital-gradient criterion, inner selected-moment
and number criteria, physical electron closure, and settled physical energy.
Iteration/evaluation limits, stalled searches and inactive-particle inconsistency
are explicit failures. Finite-difference first-order stationarity is not a global
minimum or positive-curvature certificate. The method does not round fractional
inactive occupations into a convenient fixed integer impurity population.

### Physical energy convention

A projected correlated reference can have fractional cluster particle count, so
its expectation is contracted directly; it is not mislabeled as an integer-N
impurity state. `VivoECCEnergyContribution` stores both reference and impurity
contributions, and their difference.

Single-fragment output uses:

    E = <reference|H_physical|reference>
        + <impurity|H_cluster|impurity>
        - <reference|H_cluster|reference>.

This retains the inactive reference contribution and its interactions rather
than dropping them or counting a core twice.

Partition mode defines a reference-restored democratic correction. With w_p = 1
for a fragment orbital and zero for its bath, the one-body energy operator is
(w_p+w_q)/2 times (h_bare + V_env/2); the two-body integral weight is
(w_p+w_q+w_r+w_s)/4. The sum of **high-minus-reference** weighted cluster
contributions is added to the physical reference energy. Overlapping full
impurity energies must not be summed as a total energy. Auxiliary correlation
potentials and chemical-potential energies are not part of the physical scalar.

This convention has tested exact full-bath and occupied separable-core limits.
A complete global N-representable stitched two-RDM, a variational upper bound for
arbitrary partitions, and chemical accuracy of truncated baths are not implied.

## Workflow and CLI integration

`chemistry-run` dispatches the original and advanced request schemas explicitly.
`chemistry-solve` dispatches legacy methods, direct CI, tensor CCSD, multistate
CASSCF and the typed ECC solver request. Unknown schemas, ambiguous methods and
malformed advanced requests fail; none are retried as a lower-fidelity method.

```bash
numivivo chemistry-template h2-ecc-dmet --output h2.ecc.json
numivivo chemistry-run h2.ecc.json \
  --store .numivivo/chemistry-artifacts --output h2.ecc.result.json

numivivo chemistry-template solver-ecc-dmet --output ecc-solver.json
numivivo chemistry-solve hamiltonian.json \
  --solver ecc-solver.json --output embedded-result.json
```

The H2 template is an exact full-bath, fixed-orbital example. Nontrivial orbital
optimization, nonzero initial one-/two-body feedback and a reduced cluster with
occupied inactive environment are exercised by `ECCIntegrationChecks.swift`.
A successful H2 template alone is not evidence of nontrivial QIO optimization.

Advanced molecular requests also cover density-fitted MP2 with a supplied
auxiliary basis, tensor CCSD, direct CI, multistate CASSCF and smooth C-PCM. They
share the existing artifact store, DAG scheduling and primitive hash format.
Changing only the high-level solver preserves integral/RHF/Hamiltonian reuse.
Results exported with `--output` retain a sibling `.receipt.json`.

The ECC output validator reconstructs the final reference, bath projectors,
impurity states, RDMs, moment/number conditions and physical energy decomposition.
It independently recomputes the outer orbital gradient. It rejects a cached
potential that requires further feedback to converge; it never silently repairs
that potential and then certifies the stale result.

The shared `VivoFingerprint` and `VivoArtifactValidationError` declarations now
live in `Artifacts/VivoArtifactPrimitives.swift`, independent of the program
compiler and biological host model. Valid JSON fingerprint encoding is unchanged.
Constructing an invalid standalone fingerprint now reports an artifact-validation
error instead of a ProgramPack header error; callers matching that error case
should adjust. The generic JSON canonical-data helper now reports a Foundation
`EncodingError` instead of a Metal-runtime configuration error for unsupported
root shapes. Its successful canonical bytes and accepted root shapes are
unchanged. No successful decoding/identity convention changed.

## Executable external conformance

```bash
bash Tools/check_external_chemistry.sh /tmp/numivivo-conformance
```

This creates an isolated pinned PySCF environment unless one is supplied with
`VIVO_ORACLE_PYTHON`. Python is an oracle dependency, not a runtime dependency of
NumiVivo numerics or the electronic CLI. `VIVO_REFERENCE_DIR` can select existing
hashed fixtures. The complete three-molecule fixture set and all expected methods
are mandatory; missing comparisons are not silently skipped.

The runner freezes the actual production Swift compilation inputs and records
source hashes, compiler/platform identity and reference manifest. It executes:

- 45 molecular comparisons for H2, LiH and water: overlap/core integrals, HF,
  MP2, tensor CCSD, three FCI roots, fitted ERIs/rank/HF/MP2, smooth solvent for H2
  and water, and state-averaged LiH CASSCF energies.
- 25 native ECC tests: exact frame transforms, bath/particle constraints, actual
  one-/two-moment feedback, inactive-energy restoration, nontrivial QIO updates,
  and rejection of forged convergence claims.
- 9 smooth-solvent diagnostics and negative tests, plus 29 earlier regressions.
- 25 independent ECC reconstruction comparisons when an oracle interpreter is
  present; CI runs this reconstruction as a separate mandatory job.

`compare_ecc_result.py` independently builds PySCF FCI actions, restricts CISD
states, computes density matrices and NumPy bath SVDs, and reconstructs physical
cluster and global energies. The supplied orbital frame/potential defines the
fixed point being checked; native energies are not inputs to an oracle fit.

Smooth-solvent checks cover Lebedev moments, the variational reaction-potential
identity, translational invariance, true surface-equation residual, electron
count and allocation limits, and rejection of an outlying point charge without
an enclosing-cavity policy. Agreement with PySCF is for matched model and cavity
settings; it does not provide nonelectrostatic solvation terms, analytic nuclear
solvent gradients, thermal free energies or a reaction-rate validation.

### Actual Apple integration, without interface substitutes

```bash
bash Tools/NativeChemistry/build_scoped_cli.sh /tmp/numivivo-cli
python3 Tools/NativeChemistry/check_cli_integration.py \
  /tmp/numivivo-cli/numivivo-chemistry \
  /tmp/numivivo-conformance/reference /tmp/numivivo-cli-checks
```

The macOS integration build uses the real production electronic command class,
planners, numerics, CryptoKit and artifact store. A tiny test executable entrypoint
replaces only the full program's unrelated command router. Static linkage keeps
the executable fingerprint bound to the numerical implementation. No generated
interface mirrors, fake stores or numerical substitutes are present.

The `Native chemistry conformance` Actions workflow makes oracle generation,
Apple numerics/CLI and independent ECC reconstruction separate visible jobs.
Artifacts retain failing logs as well as successful reports. The scoped build is
not the complete biology/Metal product; report its actual runner architecture
and scope rather than inferring a full Apple build from its name.

See `Documentation/Audit/ecc-conformance-observations.json` for retained execution
provenance and outcomes. That audit distinguishes portable local observations
from GitHub Actions observations. No paper reaction or BTK reproduction is claimed.

## Primary methodology references

- CovAngelo, Section 2.3 and Appendix B: https://arxiv.org/html/2604.10487v1
- PySCF PCM implementation: https://pyscf.org/_modules/pyscf/solvent/pcm.html
- Pinned external reference software: PySCF 2.8.0, generated by
  `ReferenceAdapters/PySCF/export_native_reference.py`.
