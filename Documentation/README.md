# NumiVivo documentation

Start with the [project overview](../README.md) for the scientific scope and first experiment. Use the [capability map](CAPABILITIES.md) before relying on a numerical method, import format or cross-system connection.

## Choose a workflow

| Task | Entry point |
|---|---|
| Build the Apple package and run a supplied experiment | [Get started](../README.md#get-started) |
| Calculate a small native electronic-structure example | [H₂/STO-3G chemistry](../Examples/native-chemistry/README.md) |
| Compare discrete and continuous reaction components | [Hybrid reaction runtime](../Examples/hybrid-reaction-runtime/README.md) |
| Model exposure-driven target occupancy | [Target engagement](../Examples/target-engagement/README.md) |
| Import a prepared molecular system and run MD stages | [Prepared-system MD](../Examples/md-preparation/README.md) |
| Inspect or resume an MD run | [Protocol execution and restart](Design/MD_PROTOCOL_WORKFLOW.md) |
| Read and verify stored coordinate samples | [Trajectory archive format](Design/MD_TRAJECTORY_ARCHIVE.md) |
| Explore molecular-control and host-context concepts | [Digital tissue homeostasis](../Examples/digital-tissue-homeostasis/README.md) |

Examples are computational fixtures, not experimental protocols. Their README files identify required inputs and whether a command is an Apple-side qualification task. Synthetic values are not biological measurements.

## Follow the implementation

**Structure and molecular dynamics.** Read [Molecular foundation](Design/MOLECULAR_FOUNDATION_WAVE_A.md), [MD protocols](Design/MD_PROTOCOL_WORKFLOW.md), [trajectory storage](Design/MD_TRAJECTORY_ARCHIVE.md), and the [Wave B consolidation audit](Audit/WAVE_B_CONSOLIDATION.md). The current MD protocol command path is `VivoMDProtocolCLICommands` → `VivoMDProtocolRunner` → `VivoMDMetalRuntime`; trajectory persistence uses the shared `VivoArtifactStore`.

**Reaction execution.** Read the [ProgramPack backend](Design/PROGRAM_PACK_METAL_BACKEND.md), [executable hybrid runtime](Design/EXECUTABLE_HYBRID_RUNTIME.md), and [exact/hybrid design](Design/EXACT_AND_HYBRID_STOCHASTIC.md). The executable contracts take precedence over older high-level fidelity descriptions.

**Electronic structure and embedding.** Start with the [native example](../Examples/native-chemistry/README.md), then inspect the [QM](../Sources/NumiVivoKit/QM), [QM environment](../Sources/NumiVivoKit/QMEnv), [many-body](../Sources/NumiVivoKit/ManyBody) and [embedding](../Sources/NumiVivoKit/Embedding) modules. These are bounded research implementations. Consult the [capability map](CAPABILITIES.md) for distinctions between small-system methods and the larger intended chemistry workflow.

**Physiology, populations and kinetics.** Read [Target engagement](Design/TARGET_ENGAGEMENT.md), [multicellular and partition models](Design/MULTICELLULAR_AND_PARTITION.md), and the source in [Physiology](../Sources/NumiVivoKit/Physiology), [Population](../Sources/NumiVivoKit/Population) and [Coupling](../Sources/NumiVivoKit/Coupling).

**Experiments and evidence.** Read [Artifacts and provenance](Design/ARTIFACTS_AND_PROVENANCE.md), [surrogate authority](Design/SURROGATE_AUTHORITY.md), and [artifact-store hardening](Audit/ARTIFACT_STORE_HARDENING.md). Experiment definitions, source identity, numerical configuration and evidence are separate data objects; a checksum is not a scientific validation result.

## Status and qualification

[Capability map](CAPABILITIES.md) · [Audit entry point](../AUDIT.md) · [Execution audit](Audit/2026-09-04_EXECUTION_AUDIT.md) · [Wave B audit](Audit/WAVE_B_CONSOLIDATION.md)

Historical audit reports apply to the source identities they record. They do not establish that subsequent changes compiled or passed the same checks. Source inspection, portable assertions, Apple package compilation, GPU numerical comparisons and biological validation are different evidence levels.

A useful qualification record includes the exact commit, input/configuration fingerprints, toolchain, device, command, observed result and stated acceptance criterion. Avoid reporting an entire module as validated because one example completed.

## Contribute or report a problem

Use [CONTRIBUTING.md](../CONTRIBUTING.md) for source ownership and change requirements. Use [SECURITY.md](../SECURITY.md) for the security-reporting boundary. Include small nonconfidential inputs and reproducible commands when reporting a numerical or integration failure.
