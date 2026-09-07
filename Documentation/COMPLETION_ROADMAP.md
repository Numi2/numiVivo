# NumiVivo completion roadmap

The product objective is a leading Apple-native tool for reproducible molecular-to-biological research. The scope is the connected workflow described in the README, capability map and engineering principles. Passing one build, one molecular fixture or one numerical suite does not complete that objective.

This roadmap orders the work without replacing the intended architecture with a smaller product. Unsupported methods must remain explicit, and scientific claims must retain their actual evidence class. An unavailable author dataset can block a particular reproduction; it does not block improving the general platform.

## Acceptance requirements

| Requirement | Evidence needed for completion | Current work or remaining gap |
| --- | --- | --- |
| Reliable execution and recovery | Clean full Apple build; verified output/cache identities; cancellation and concurrent caller tests; injected storage-failure recovery; uninterrupted versus resumed results across complete protocols. | September 2026 work closes shared-task cancellation and MD resume/publication defects. Native verification is recorded separately; these fixes do not qualify every runtime. |
| Correct molecular and stochastic dynamics | Independent force/energy checks, conservation and constrained-ensemble checks, PME convergence, exact stochastic distributions, hybrid authority conservation, transactional rejection, and restart invariance over the supported numerical profiles. | Seven dedicated hybrid native conformance tests passed on M4 Pro; their exact scope is recorded in the reliability audit. MD electrostatics, constraints, thermostats, barostats and force-provider combinations require a maintained coverage matrix, rather than a passing harmonic example alone. |
| Qualified chemistry and sampled uncertainty | Independent electronic references, mapped path/connectivity evidence, declared free-energy and nuclear approximations, convergence studies, surrogate held-out errors and exact acceptance accounting; uncertainties remain correlated and missing sources remain visible. | The reaction gate is repaired at `9f9da42`. Revision `59e1868` adds bead-convergence, correlated-sampling and surrogate-coverage assessments; their declared scope is in the nuclear/surrogate contract. No protein-scale kinetic or general chemical accuracy claim follows from small-system gates. |
| Connected research workflows | Complete prepared-system → sampling → electronic/path calculation → conditional kinetic evidence → observable workflow, with units, mapping, context, parameter provenance and assumptions retained at every edge; both successful and missing-input cases. | Existing adapters and fixtures cover portions of the route. A representative full-route campaign, reproducible from published inputs and commands, is still needed. |
| Complete declared runtime semantics | Source, native tests and counterexamples for the accepted language/backend contract; persistent state includes everything that can affect continuation; coupled participants commit the same logical boundary. | ProgramPack delays/refractory behavior, broader spatial interfaces, live authority migration and stronger multi-participant execution remain separate from the well-mixed hybrid backend. Follow the linked design contracts rather than interpreting fidelity names as implementation evidence. |
| Demonstrated useful scale and performance | Published systems, source/binary/device identities, correctness tolerances, preparation and warmup policy, elapsed time, throughput, memory, and accepted-work accounting; compare like-for-like numerical models against named independent references where available. | No broad speedup or maximum biological scale is established. Measure complete workflows and bottlenecks before promoting an acceleration claim. |
| A usable, maintainable research product | Documented complete examples, discoverable CLI/API operations, actionable bounded errors, readable progress and result inspection, automated native regression gates, install/release checks, and retained evidence for failures as well as successes. | Workflow activity now has a concrete ownership model. CLI examples, operation coverage, long-run inspection and release qualification must stay aligned with the actual package. |

## Execution order

1. Completed the cancellation, MD recovery and hybrid conformance milestone with [native Mac mini evidence](Audit/2026-09-07_COMPLETION_RELIABILITY.md): a clean complete build, 224 baseline tests, 22 final-source regressions and 57 CLI checks. This milestone does not complete the acceptance table.
2. Completed the [streaming archive recovery milestone](Audit/2026-09-07_ARCHIVE_STREAMING.md): a real 100,001-chunk archive validates, resumes and extends through production code while preserving its prefix. Fresh CLI processes verify the archive with memory independent of its length for the measured fixture. The explicit materializing index retains its bound. Turn new integrity or recovery failures into owning-runtime repairs; do not relax gates to obtain a pass.
3. Expand the MD force/ensemble/continuation matrix and complete nuclear/surrogate qualification in parallel, respecting physical GPU ownership.
4. Run a prepared, mapped full-route research campaign using explicit published inputs and report every remaining approximation.
5. Use measured correctness and resource bottlenecks to prioritize scale improvements and broader runtime semantics. Keep permanent interfaces and one numerical authority per represented state.
6. Audit the entire acceptance table against current source, native runs, artifacts and release behavior before claiming product completion.

## Source contracts

- [Capability map](CAPABILITIES.md) and [engineering principles](PRINCIPLES.md).
- [General workflows](Design/GENERAL_WORKFLOW_PLATFORM.md), [MD protocols](Design/MD_PROTOCOL_WORKFLOW.md), and [trajectory archives](Design/MD_TRAJECTORY_ARCHIVE.md).
- [Hybrid execution](Design/EXECUTABLE_HYBRID_RUNTIME.md), [ProgramPack backend](Design/PROGRAM_PACK_METAL_BACKEND.md), and [transactional runtime](Design/TRANSACTIONAL_RUNTIME.md).
- [Prepared molecular workflows](PreparedMolecularWorkflows.md), [QM/MM free energy](QMMMFreeEnergy.md), and [nuclear, uncertainty and surrogate methods](Design/NUCLEAR_UNCERTAINTY_SURROGATE.md).

Measured validation records must identify the tested source and scope. Historical `VALIDATION_STATUS.md` results remain evidence for their named revisions, not an automatic certificate for later code.
