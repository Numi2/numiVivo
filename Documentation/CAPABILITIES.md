# Capability map and current boundaries

This is a source-backed map of the implementation, not a feature certification. The public package combines mature contracts, new numerical source and experimental methods. A successful import, a function named after a method, or a committed example does not establish numerical agreement or biological applicability.

The navigation and implementation references below describe the reviewed source beginning at `b77fb7885ce2c484a10f2509798cb4a05a294047`, with subsequent integration corrections. Consult the actual revision when a method changes. This map deliberately does not assign an unmeasured speedup or maximum biological scale.

## Molecular preparation and dynamics

| Capability | Implemented source | Important boundary |
|---|---|---|
| Structure identity | Atoms, bonds, residues, chains, conformers, selections, mapping and periodic cells. | Archival representation does not imply every geometry is executable by every backend. |
| Interchange | PDB, mmCIF, SDF/MOL V2000, MOL2 and a strict SMILES subset. | Full stereochemistry, V3000, query semantics and other unrepresented chemistry are not general supported import claims. |
| Force fields | Unit-explicit particles and bonded terms, nonbonded exceptions, direct type-pair tables and AMBER topology/restart import. | Native general atom typing, AM1-BCC generation and protonation prediction remain distinct from importing prepared parameters. |
| Virtual sites and constraints | Linear virtual sites, parent-force redistribution and bounded distance-constraint projection; executable AMBER preparation helpers. | General virtual-site frames, polarizable/Drude execution and every AMBER extension are not supported. Constraint-solver accuracy still requires independent checks. |
| MD integrators | NVE and Langevin-middle NVT source; molecular-center Monte Carlo NPT. | Ensemble correctness, conservation, equilibration and speed are not established by implementation alone. |
| Nonbonded execution | Bounded neighbor construction, LJ, cutoff/reaction-field electrostatics and a PME mesh/FFT path. | PME settings are planning controls, not certified force-error bounds. Compare exception conventions, reciprocal accuracy, cutoff treatment and pressure sensitivity with independent references. |
| Periodic geometry | Orthogonal-cell execution preflight; general triclinic archival cells. | The current nearest-image MD profile rejects skew cells. A truncated-octahedron import must not be silently converted into an orthogonal simulation. |
| Protocols and storage | Minimize → NVT → NPT → production stages; explicit reconfiguration; immutable restart cursors; bounded binary trajectory chunks. | Supplied durations are illustrative. Trajectories are not checkpoints. Interrupted minimization restarts from accepted geometry rather than serialized optimizer history. |

Implementation: [MD source](../Sources/NumiVivoKit/MD), [structure foundation](Design/MOLECULAR_FOUNDATION_WAVE_A.md), [protocol contract](Design/MD_PROTOCOL_WORKFLOW.md), [trajectory format](Design/MD_TRAJECTORY_ARCHIVE.md).

## Electronic structure, embedding and reaction research

| Capability | Implemented source | Important boundary |
|---|---|---|
| Gaussian electronic structure | Native all-electron Cartesian Gaussian integrals, matrix algebra and restricted/unrestricted Hartree–Fock. | Dense FP64 work with explicit budgets; not a general density-fitting, ECP or spherical-basis implementation. |
| Density-functional work | Restricted LDA and its numerical-grid/functional components. | Not the full functional, dispersion and gradient coverage of established quantum-chemistry packages. |
| Correlation | Embedded-Hamiltonian records, restricted MP2 and small-system configuration interaction used by the native example. | Dense tensor storage and determinant growth limit the intended problem size. Do not infer production CCSD, CASSCF or DMRG coverage from the roadmap. |
| QM/MM and solvent | Electrostatic embedding, boundary-link/Z1 and LJ coupling machinery; C-PCM reaction-field and restricted self-consistent solvent work. | A continuum electrostatic polarization term is not a complete activation Gibbs free energy. Boundary and energy accounting need independent verification. |
| Orbital information and QIO | Correlated orbital-information machinery, subspace transport and shared-path orbital optimization under bounded CI and evaluation budgets. | This is not a reproduced ECC-DMET protein barrier calculation. Near-degenerate states require treatment beyond the current single-state path contract. |
| Quantum-algorithm research | Fermion/Pauli mapping, state-vector and variational-method source, plus resource-estimation components. | These are experimental algorithm modules, not a general QPU service or established quantum advantage. Do not assume Metal execution merely from the platform target. |
| Reaction paths and networks | Dedicated geometry, reaction-path and reaction-network modules. | A path representation or constrained scan does not certify a transition state, connected mechanism, free energy or reaction rate. |

Start with [native chemistry](../Examples/native-chemistry/README.md). Inspect [QM](../Sources/NumiVivoKit/QM), [QMEnv](../Sources/NumiVivoKit/QMEnv), [ManyBody](../Sources/NumiVivoKit/ManyBody), [Embedding](../Sources/NumiVivoKit/Embedding), [Quantum](../Sources/NumiVivoKit/Quantum), [ReactionPath](../Sources/NumiVivoKit/ReactionPath) and [ReactionAtlas](../Sources/NumiVivoKit/ReactionAtlas).

## Reaction kinetics, physiology and biological control

| Capability | Implemented source | Important boundary |
|---|---|---|
| Typed programs | Semantic validation, dimensional checks, executable pack references and runtime monitoring. | Earlier language illustrations may contain future features; accepted backend contracts are authoritative. A mathematical program is not a validated genetic construct. |
| Deterministic/stochastic kinetics | ProgramPack deterministic and bounded stochastic paths; a separate exact-SSA/tau/RK2 runtime with UInt32 counts and FP32 continuous state. | Mixed numerical authorities are currently separated by actual dependency components. This is not arbitrary connected-network hybrid dynamics. |
| Temporal and spatial behavior | Supported temporal rule/monitor operations and concentration finite-volume transport. | General delayed/refractory state, temporal rate/gate operators, membrane-interface physics, discrete hopping and live authority migration are not completed by the existing fidelity labels. |
| Target engagement | Exposure-driven free, reversible, covalent and optional competitor occupancy, shared turnover, FP64 reference and Metal cohort integration. | Exposure is externally maintained; it is not a finite drug pool. Occupancy is not cell survival, efficacy, toxicity or patient-specific response. |
| Energy-to-rate evidence | Conditional conversion from an explicitly identified activation Gibbs free energy relative to a bound pre-reactive complex. | An electronic barrier alone is rejected. Other kinetic parameters and total uncertainty are not automatically inferred. |
| Physiology and coupling | Compartment/partition models, molecular–physiology exchange, shared-state contracts and participant infrastructure. | Source compilation or logical transaction ownership is not whole-organism validation or crash-atomic distributed execution. |
| Calibration and studies | Likelihood, estimation and campaign components; target-engagement held-out study evaluation. | A study evaluator is not itself a parameter fitter or posterior sampler. Imported or synthetic observations must retain their evidence class. |

See [target engagement](Design/TARGET_ENGAGEMENT.md), [hybrid execution](Design/EXECUTABLE_HYBRID_RUNTIME.md), [ProgramPack limits](Design/PROGRAM_PACK_METAL_BACKEND.md), and [multicellular design](Design/MULTICELLULAR_AND_PARTITION.md).

## Evidence, ownership and execution

The shared artifact store provides content identity and checked object/reference persistence. MD protocols reuse it for state, observations and trajectory prefixes. Checkpoints bind represented state to model/configuration and, where implemented, a separately versioned numerical profile. A restored seed and accepted-step index identify a random namespace; floating-point scheduling, mesh accumulation and changed numerical implementations can still change trajectories.

Surrogate declarations carry authority and uncertainty contracts. A learned prediction must not silently become mechanistic evidence, and missing uncertainty is not zero uncertainty. Model validation, numerical validation, biological validation and treatment authorization are separate concepts.

The current full-package target is an Apple toolchain with Swift 6, C++23 and Metal. Selected Foundation/C++ components have portable harnesses. Historical portable results are documented in [the audit](../AUDIT.md) and the relevant example guides; they do not qualify the current full Apple build or any untested GPU module.

## Next development priorities

Consolidate shared semantics and execution interfaces before expanding duplicate runtimes. Qualify force laws, electrostatics, constrained ensembles, stochastic behavior, checkpoint continuation and failure recovery on Apple hardware. Expand chemistry coverage against independent references, then establish a carefully qualified connection from molecular energetics to dynamic observations.

The intended endpoint is a reproducible molecular-to-biological research workflow. It is not a claim that the repository already replaces every application discussed in the CovAngelo paper.
