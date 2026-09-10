# Capability map and current boundaries

This is a source-backed map of the implementation, not a feature certification. The public package combines mature contracts, new numerical source and experimental methods. A successful import, a function named after a method, or a committed example does not establish numerical agreement or biological applicability.

The navigation and implementation references below describe the reviewed source beginning at `b77fb7885ce2c484a10f2509798cb4a05a294047`, with subsequent integration corrections. Consult the actual revision when a method changes. This map deliberately does not assign an unmeasured speedup or maximum biological scale.

## Single-cell analysis and biological prediction

The single-cell entries were reviewed against published evidence available at
`aa85e1a0` on 2026-09-10. Historical receipts retain their named runtime identities.

| Capability | Implemented and measured scope | Remaining boundary |
| --- | --- | --- |
| H5AD and multi-assay interchange | Native count import/export, explicit metadata mapping and source preservation; separate RNA/protein/ATAC spaces and H5MU routes. | Supported encodings and bounds are explicit; interchange does not validate experimental labels or joint biological inference. |
| Count inference | Native negative-binomial fitting, dispersion/shrinkage, offsets and paired/batch designs; real-data numerical and R-reference comparisons. | Null-benchmark failures, varying support and broader FDR/interval calibration remain. Method concordance is not ground truth. |
| Reduction and integration | Sparse HVG/PCA, neighbors, Louvain, UMAP-compatible optimization and donor correction. Full HIRISA PCA/graph and seed-7 integration have numerical/replay evidence. | Full HIRISA clustering qualification remains open. Coarse integration margins pass, but marker/program, rare-cell and native multi-seed preservation do not yet close. |
| Known-treatment donor prediction | Native control-context ridge plus three simple baselines. All 79 HIRISA folds complete; ridge beats no-change in 14/16 contrasts and mean response in 4/16. | Requires the new donor's control profile and matching context. Two contrasts fail against no-change; no calibrated intervals or general phenotype claim. |
| Target/composition prediction | Native composition models reproduce all 131 Norman held-out pairs; native GO kernel evaluates 105 held-target folds, with 101 supported descriptors. | GO gain over mean is only 0.58% on reused data, with 29/101 worse than no-change. Reliable independent-study, unseen-context and mechanistic prediction remain unqualified. |
| Annotation | Sparse fixed marker/program scores and provenance-bound candidate reference labels. | Labels remain nonauthoritative; annotation calibration and novel-class rejection remain separate. |

Use the [prediction assessment](BiologicalPrediction.md) for input requirements,
positive and negative results and next evidence gates, and the
[single-cell roadmap](SingleCellInteroperability.md) for commands and detailed
owner-specific limits. Count/latent storage is bounded or streamed in specific
stages; cell-scale metadata and several model arrays still reside in memory.
No general Metal speedup for this pipeline is established.

## Molecular preparation and dynamics

| Capability | Implemented source | Important boundary |
|---|---|---|
| Structure identity | Atoms, bonds, residues, chains, conformers, selections, mapping and periodic cells. | Archival representation does not imply every geometry is executable by every backend. |
| Interchange | PDB, mmCIF, SDF/MOL V2000, MOL2 and a strict SMILES subset. | Full stereochemistry, V3000, query semantics and other unrepresented chemistry are not general supported import claims. |
| Force fields | Unit-explicit particles and bonded terms, nonbonded exceptions, direct type-pair tables and AMBER topology/restart import. | Native general atom typing, AM1-BCC generation and protonation prediction remain distinct from importing prepared parameters. |
| Virtual sites and constraints | Linear, out-of-plane and local-coordinate virtual-site definitions with parent-force redistribution; bounded distance-constraint projection and executable AMBER preparation helpers. | Each represented site model needs its own qualification. Unrepresented site definitions, Drude particles and unsupported AMBER extensions remain rejected. Compensated mode currently requires physical atoms. |
| MD integrators | NVE and Langevin-middle NVT source; molecular-center Monte Carlo NPT. | Ensemble correctness, conservation, equilibration and speed are not established by implementation alone. |
| Compensated coordinates | Opt-in high/low FP32 positions, RATTLE drift and exact-word checkpoint storage in numerical profile v6. | Fixed-cell classical NVE/NVT with physical atoms; excludes minimization, NPT, virtual/Drude sites, external force providers and FP32 trajectory archives. Exact stored words do not guarantee general bitwise PME replay. |
| Nonbonded execution | Bounded neighbor construction, LJ, cutoff/reaction-field electrostatics and a PME mesh/FFT path. | PME settings are planning controls, not certified force-error bounds. Compare exception conventions, reciprocal accuracy, cutoff treatment and pressure sensitivity with independent references. |
| Periodic geometry | Orthogonal-cell execution preflight; general triclinic archival cells. | The current nearest-image MD profile rejects skew cells. A truncated-octahedron import must not be silently converted into an orthogonal simulation. |
| Adaptive sampling and state handoff | Seeded replica blocks, explicit observable diagnostics, bounded accepted-prefix selection, exact checkpoint/snapshot export, immutable verification and mapped geometry seeds for fresh reaction qualification. | A selected prefix can remain unconverged; geometry transfer does not establish QM ensemble sampling. [Full-route evidence](Audit/2026-09-08_PREPARED_REACTION_CAMPAIGN.md), [export interface](Design/MOLECULAR_SAMPLING_EXPORT.md) and [reaction handoff](Design/PREPARED_REACTION_CAMPAIGN.md). |
| Protocols and storage | Minimize → NVT → NPT → production stages; explicit reconfiguration; immutable restart cursors; bounded binary trajectory chunks. | Supplied durations are illustrative. Trajectories are not checkpoints. Interrupted minimization restarts from accepted geometry rather than serialized optimizer history. |

Implementation: [MD source](../Sources/NumiVivoKit/MD), [structure foundation](Design/MOLECULAR_FOUNDATION_WAVE_A.md), [protocol contract](Design/MD_PROTOCOL_WORKFLOW.md), [trajectory format](Design/MD_TRAJECTORY_ARCHIVE.md).

The [v3 native qualification audit](Audit/2026-09-08_MD_NUMERICAL_QUALIFICATION.md) supplies independent charged-PME refinement, finite-time Langevin velocity statistics, fixed-geometry constrained thermalization and ideal molecular NPT transition evidence for its named fixtures. It also records numerical-contract rejection for old continuation state and explicit-grid resource admission. These measured cases do not qualify general interacting or constrained equilibrium ensembles.

The [frontier MD audit](Audit/2026-09-08_FRONTIER_MD_REFERENCE_PANEL.md) adds seven
prepared realistic systems, independent energies/forces, strict constraints and
a separately identified smooth-model NVE refinement panel. It records the
torsion-sign, coordinate-precision and RATTLE repairs, failed original inputs,
sharp-cutoff limitations and nondeterministic PME accumulation separately.

## Electronic structure, embedding and reaction research

| Capability | Implemented source | Important boundary |
|---|---|---|
| Gaussian electronic structure | Native all-electron Cartesian Gaussian integrals, restricted/unrestricted Hartree–Fock, Coulomb density-fitting factors, factorized restricted Hartree–Fock and analytical HF nuclear derivatives. | Explicit FP64 work and memory budgets; supported factorized paths need their own reference qualification. ECP and general spherical-basis coverage are not established. |
| Density-functional work | Restricted LDA and its numerical-grid/functional components. | Not the full functional, dispersion and gradient coverage of established quantum-chemistry packages. |
| Correlation | Embedded-Hamiltonian records, dense and factorized restricted MP2, and small-system configuration interaction used by the native example. | Factorization avoids selected dense tensors; determinant growth and retained intermediates still constrain scale. Do not infer production CCSD, CASSCF or DMRG coverage from the roadmap. |
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
| Conditional reaction encounters | Fully reconstructed molecular TST rates feed a tagged mapped reactant's first-event probability, flux and local sensitivities under explicit maintained reservoirs. | Molecularity, standard-state units, component identity, temperature and transmission remain bound. No finite-pool depletion, competing-channel multiplicity, reverse kinetics, occupancy or physical uncertainty is inferred. [Finite native/CLI campaign](Audit/2026-09-08_PREPARED_REACTION_CAMPAIGN.md). |
| Nuclear statistics and global uncertainty | One-dimensional Wigner/Eckart factors, fixed-cell ring-polymer equilibrium sampling, correlated-chain/bead convergence evidence and nonlinear propagation of supplied joint draws. | These are not multidimensional instantons, real-time quantum recrossing or inference of missing uncertainty. A prepared system remains unqualified until its own convergence and model-applicability evidence passes. |
| Reactive surrogate acceleration | Group-held-out energy/force delta fitting, Metal inference, authoritative endpoint correction and grouped external-domain coverage evidence. | Corrected equilibrium proposals are not physical MD. Coverage does not make learned forces authoritative or establish speedup, transport or residence time. |
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
