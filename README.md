# NumiVivo

**From atoms to biological behavior. Built for Apple silicon.**

NumiVivo is a research platform for molecular simulation and programmable biology. It brings molecular structures, GPU molecular dynamics, native electronic structure, reaction networks, and reproducible experiments into one Swift, C++ and Metal codebase.

The ambition is to follow a molecular change across scales: how a structure moves, where electrons interact, which reactions become possible, and how those reactions affect a cell or physiological model. Every transition should preserve the quantities, assumptions and evidence needed to understand the result.

[Get started](#get-started) · [Explore the capabilities](#explore-the-capabilities) · [Examples](#choose-an-experiment) · [Documentation](Documentation/README.md) · [Implementation status](Documentation/CAPABILITIES.md)

> **Research software, under active development.** The capabilities below describe source implementations and their intended workflows—not a fully qualified release. Apple package integration, GPU numerical behavior and performance require qualification. The [capability map](Documentation/CAPABILITIES.md) separates implemented methods, current restrictions and planned work.

## One scientific question, several scales

NumiVivo is being developed to connect a molecular hypothesis to a dynamic biological model:

```text
Structures and force fields
        ↓
Molecular dynamics and sampled configurations
        ↓
Electronic structure, QM/MM and orbital embedding
        ↓
Reaction energetics and explicitly qualified kinetic evidence
        ↓
Reaction networks, target engagement and physiological dynamics
        ↓
Versioned experiments, observations and reproducible results
```

This is the integration direction, not a claim that every arrow is already an automated, validated calculation. In particular, an electronic-energy difference is not automatically an activation free energy, reaction rate, binding affinity or biological outcome.

## Explore the capabilities

| Area | What the source provides | Explore |
|---|---|---|
| **Molecular foundations** | Canonical atoms, bonds, residues, conformers and periodic cells; selections and atom mapping; PDB, mmCIF, SDF/MOL V2000, MOL2 and a strict SMILES subset; topology preparation and unit-explicit force-field records. | [Structure and force-field foundation](Documentation/Design/MOLECULAR_FOUNDATION_WAVE_A.md) |
| **Apple-native molecular dynamics** | Metal force kernels, neighbor construction, minimization, NVE and Langevin NVT, molecular-center Monte Carlo NPT, PME electrostatics, distance constraints and linear virtual sites. Preparation protocols retain stage identity, restart state and bounded trajectory output. | [MD protocols](Documentation/Design/MD_PROTOCOL_WORKFLOW.md) |
| **Native electronic structure** | FP64 Cartesian Gaussian integrals, restricted/unrestricted Hartree–Fock, restricted LDA, embedded Hamiltonians, MP2 and small-system configuration interaction. The chemistry path runs without Python callbacks or a CUDA runtime. | [Native chemistry example](Examples/native-chemistry/README.md) |
| **Embedding and reaction research** | QM/MM electrostatic and boundary-link machinery, C-PCM reaction-field work, correlated orbital information, orbital-subspace alignment and a bounded path-consistent QIO optimizer. These are experimental methods, not a reproduced protein-reaction result. | [Embedding source](Sources/NumiVivoKit/Embedding) · [QM environment source](Sources/NumiVivoKit/QMEnv) |
| **Programmable reaction dynamics** | Typed molecular programs and compiled ProgramPacks; deterministic kinetics; discrete stochastic simulation; exact-SSA/tau-leap/RK2 execution across dependency-separated components; temporal rules, monitors and concentration transport within declared backend limits. | [Hybrid runtime](Examples/hybrid-reaction-runtime/README.md) · [ProgramPack backend](Documentation/Design/PROGRAM_PACK_METAL_BACKEND.md) |
| **Target engagement and physiology** | Exposure-driven reversible binding, covalent conversion, competition and turnover; a native FP64 reference and an existing-runtime Metal path; physiological exchange and molecular–physiology coupling contracts. | [Target-engagement example](Examples/target-engagement/README.md) |
| **Reproducible experiments** | Content-addressed artifacts and tasks, evidence references, configuration-bound checkpoints, staged protocols, compact trajectory chunks, observation records and explicit failure results. | [Artifacts and provenance](Documentation/Design/ARTIFACTS_AND_PROVENANCE.md) · [Trajectory storage](Documentation/Design/MD_TRAJECTORY_ARCHIVE.md) |

The repository also contains experimental quantum-algorithm, reaction-path, reaction-network, population, calibration and surrogate components. Their presence is not a claim of complete Qiskit, ORCA, GROMACS or physiological-modeling equivalence. See the [source-backed capability map](Documentation/CAPABILITIES.md) before choosing a backend.

## Built around Apple silicon

**Native execution rather than a Python simulation loop.** Swift manages scientific objects, experiments and concurrent operations. C++23 implements compilation, validation and portable reference components. Metal executes molecular-dynamics and kinetic kernels. Precision-sensitive chemistry retains FP64 native CPU calculations rather than forcing every calculation onto the GPU.

**State stays where the calculation needs it.** The GPU runtimes use persistent buffers, compiled tables and explicit capacities. Private simulation state is separated from host-visible commands, diagnostics and requested observations. Full coordinate readback is an explicit sampling or checkpoint operation.

**A failed candidate must not become accepted state.** Runtime transactions distinguish proposed and accepted evolution. Checkpoints bind the represented state to its model and numerical configuration. The MD numerical profile is versioned separately so a restart cannot silently substitute a different algorithm.

**Reproducibility has a defined scope.** Counter-based random namespaces retain seeds and accepted-step identity. They do not guarantee bitwise-identical trajectories across devices or compiler versions; PME accumulation and floating-point execution have additional reproducibility limits. Integrity hashes establish which bytes were used, not whether the scientific model is correct.

## Get started

The primary development target is an **Apple-silicon Mac with macOS 15 or newer**, Swift 6 and an Apple SDK providing Metal and a C++23-capable toolchain. The package also declares iOS 18 support; it is not a released iOS application. The native package has no Python or CUDA runtime dependency.

Run from a terminal on the Apple machine:

```sh
git clone https://github.com/Numi2/numiVivo.git
cd numiVivo
swift build -c debug
.build/debug/numivivo --help
```

These are build and execution instructions, not a recorded successful build of the current revision. Keep the commit SHA with compiler diagnostics and numerical results. Portable checks do not qualify the complete Apple package.

### First experiment: exposure and target occupancy

Use a supplied synthetic fixture—no downloaded dataset or prepared protein is required:

```sh
RUN_DIR="$(mktemp -d "${TMPDIR:-/tmp}/numivivo.XXXXXX")"

.build/debug/numivivo engagement-validate \
  Examples/target-engagement/synthetic-pulse.json

.build/debug/numivivo engagement-run \
  Examples/target-engagement/synthetic-pulse.json \
  --backend reference \
  --output "$RUN_DIR/engagement.json"
```

Inspect target-state fractions over the supplied exposure schedule. The inputs and observations are **synthetic mathematical fixtures**, not measured pharmacology or treatment recommendations. The [example guide](Examples/target-engagement/README.md) adds the Metal backend, compiled kinetic source, conditional free-energy conversion and held-out observation evaluation.

## Choose an experiment

| Start with | What it exercises | Input |
|---|---|---|
| [Native chemistry](Examples/native-chemistry/README.md) | Gaussian integrals → Hartree–Fock → embedded Hamiltonian → MP2/FCI → orbital information. | Included H₂/STO-3G example. |
| [Hybrid reactions](Examples/hybrid-reaction-runtime/README.md) | Exact SSA, tau-leaping and RK2 on separate reaction components; checkpoint and resume. | Included synthetic model and counts. |
| [Target engagement](Examples/target-engagement/README.md) | Exposure, binding, covalent conversion, competition and turnover. | Included synthetic experiments. |
| [Prepared-system MD](Examples/md-preparation/README.md) | AMBER import → minimization → NVT → NPT → production samples. | Your prepared topology and restart. |
| [Digital tissue homeostasis](Examples/digital-tissue-homeostasis/README.md) | Molecular-control, host-context and coupling concepts. | Included research fixtures; check the declared runtime support. |

### A prepared molecular system, one MD protocol

For an appropriately prepared **orthogonal periodic system**, the existing AMBER bridge provides the starting point. Once `system.json` and `initial.json` have been imported:

```sh
.build/debug/numivivo md-protocol-template system.json > protocol.json

# Review the durations, temperature, pressure and minimization gate first.
.build/debug/numivivo md-protocol-validate protocol.json \
  --system system.json --state initial.json

.build/debug/numivivo md-protocol-run protocol.json \
  --system system.json --state initial.json \
  --store ./md-artifacts > receipt.json
```

The receipt identifies durable restart state. Coordinate samples are stored as bounded binary chunks, not an ever-growing in-memory JSON trajectory. Later stages preserve state unless their velocity initialization is explicitly changed. A failed required minimization or rejected dynamics candidate blocks the protocol rather than silently changing the experiment.

The template is illustrative: its duration does not establish equilibration. [Import, sampling and resume details →](Documentation/Design/MD_PROTOCOL_WORKFLOW.md)

## Molecular biology inside NumiLab

The wider goal is a coupled experimental environment in which molecular reactions respond to cell state, tissue organization, physical transport and nervous-system activity.

| System | Intended responsibility |
|---|---|
| **NumiVivo** | Molecular representations, chemistry, reaction dynamics and molecular experiments. |
| **NumiTissue** | Cell populations, development and tissue organization. |
| **NumanX** | Geometry, mechanics, fluids and physical coupling. |
| **NumiBrain** | Neural, autonomic and embodied control. |

NumiVivo contains coupling contracts and participant infrastructure. Fully qualified cross-repository biological simulations remain an integration objective, not something established by this table.

## Development direction

The next milestone is **an integrated, numerically qualified research workflow**, not simply more named solvers.

The priorities are to consolidate model semantics and execution ownership; qualify the Apple MD and kinetic backends; expand electronic-structure and embedding methods against independent references; and connect reaction energetics to observations and kinetics without losing thermodynamic meaning. Membrane interfaces, general triclinic dynamics, delayed/refractory state, live fidelity migration and larger correlated calculations require additional work.

The chemistry program is informed by [CovAngelo, arXiv:2604.10487](https://arxiv.org/abs/2604.10487). Its methods motivate native QM/MM, compact embedding and consistent reaction-path treatment. This repository does not claim to reproduce that paper's results or contain its authors' unpublished implementation.

## Documentation and contribution

[Documentation index](Documentation/README.md) · [Capability map and limits](Documentation/CAPABILITIES.md) · [Execution audit](AUDIT.md) · [Contribution guide](CONTRIBUTING.md) · [Security policy](SECURITY.md)

```text
Sources/NumiVivoCore/      C++ compiler, validation, pack format and reference logic
Sources/NumiVivoKit/       Swift scientific modules, runtimes, artifacts and workflows
Sources/NumiVivoShaders/   Metal kernels and explicit shader-module loading
Sources/NumiVivoCLI/       numivivo command-line interface
Examples/                 Model fixtures and executable examples
Schemas/                  Versioned data contracts
Documentation/            Methods, architecture, audits and implementation limits
Tools/                    Development tools and portable qualification harnesses
```

Report reproducible bugs through the repository's Issues tab, including the commit, command, platform and relevant diagnostics. Keep confidential research data out of public reports. Contributions should include the affected contract, an example and an explicit validation boundary; see [CONTRIBUTING.md](CONTRIBUTING.md).

## Scientific boundary and license

NumiVivo models biological systems; it does not validate a therapy, authorize an experiment or establish safety in a living organism. An executable model is not proof of a realizable biological construct. Numerical checks, calibration evidence and biological validation are distinct.

Licensed under [Apache License 2.0](LICENSE). See [NOTICE](NOTICE). External datasets, force-field parameters, model weights and third-party materials retain their own terms. For research use, cite the exact repository revision and the methods and source data used by the calculation.
