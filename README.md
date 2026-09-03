# NumiVivo

**Apple-native software for specifying, compiling, simulating, verifying, and operating molecular programs in living systems.**

NumiVivo treats engineered biological behavior as an executable, versioned program. A program describes what a cell or tissue observes, how it combines signals over time, what state it retains, what molecular actions it may take, and the constraints under which it must stop.

```text
molecular inputs
    ↓
sensing → state → decision → biological action
    ↑                            ↓
    └──── feedback, limits, and verification
```

NumiVivo is being built as the molecular-programming layer of NumiLab:

```text
NumiLab
├── NumanX      transport, mechanics, fluids, contact, and mechanotransduction
├── NumiTissue  cells, development, differentiation, and tissue organization
├── NumiBrain   neural, autonomic, and organism-level control
└── NumiVivo    intracellular reactions, synthetic circuits, molecular memory,
                biological control programs, and experimental execution
```

## Core objective

A researcher should be able to state a desired biological behavior, compile it into explicit molecular mechanisms, simulate it across cell and tissue contexts, inspect uncertainty and failure modes, and export a reproducible design for controlled experimental validation.

NumiVivo is not a sequence editor wrapped around a simulator. It is an integrated system with five authoritative stages:

1. **Specify** — a typed language for inputs, logic, state, timing, actions, constraints, monitoring, and termination.
2. **Compile** — lower behavior into receptors, promoters, RNA regulators, protein switches, enzymes, transport processes, and reaction networks.
3. **Simulate** — execute deterministic, stochastic, spatial, multicellular, and coupled-physics models at explicit fidelity levels.
4. **Verify** — evaluate reachability, limits, hazards, context sensitivity, uncertainty, and bounded shutdown behavior.
5. **Operate** — run versioned experiments, interventions, measurements, replay, provenance, and transactional coupling to NumiLab.

## Apple-silicon execution model

NumiVivo is designed around Apple silicon rather than ported from a CPU-first or CUDA-first runtime.

- Swift coordinates experiments, artifacts, command submission, and structured concurrency.
- C++23 owns compilation, graph analysis, immutable contracts, and portable reference numerics.
- Metal owns production reaction, regulation, stochastic, diffusion, transport, reduction, and verification kernels.
- Persistent simulation state remains GPU-resident across steps.
- Private heaps hold authoritative hot state; shared memory is limited to compact commands and diagnostics.
- Programs are compiled into fixed tables, offsets, capacities, execution cohorts, and fingerprints before dispatch.
- Structure-of-arrays layouts provide coalesced access across large populations of cells.
- Sparse active sets and event queues avoid work on biologically inactive mechanisms.
- Counter-based random streams make stochastic execution reproducible across dispatch order changes.
- Each step runs against shadow state and is committed only after numerical and safety checks pass.
- Learned surrogates may run through MLX or Metal 4 tensor operations, but never replace authoritative constraints or state transitions.

## Fidelity contract

Every result declares the model class that produced it.

| Level | Scope | Primary use |
|---|---|---|
| **F0** | Boolean, threshold, temporal, and finite-state molecular logic | architecture search and static verification |
| **F1** | Deterministic well-mixed reaction networks | calibrated dynamics and parameter screening |
| **F2** | Stochastic kinetics with discrete molecular counts | noise, switching, rare-event, and population distributions |
| **F3** | Spatial reaction–diffusion with membranes and compartments | gradients, localization, transport, and morphology |
| **F4** | Multicellular tissue programs coupled to mechanics, perfusion, immunity, and organism state | integrated in-vivo hypotheses |

Promotion between levels is error-controlled and recorded. Agreement with a lower-fidelity model is not treated as biological validation.

## Program example

```yaml
apiVersion: numivivo.org/v1alpha1
kind: VivoProgram
metadata:
  name: local-inflammatory-controller
spec:
  target:
    cellType: engineered-chondrocyte

  inputs:
    - id: il1b
      source: extracellular
      unit: nM
    - id: tissue-damage
      source: host
      unit: normalized
    - id: cell-stress
      source: intracellular
      unit: normalized

  state:
    - id: inflammatory-memory
      type: leaky-integrator
      halfLife: 6 h

  rules:
    - id: activate-therapy
      when: sustained(il1b > 0.8 nM and tissue-damage > 0.55, for: 20 min)
      then:
        express: il1ra
        rate: proportional(il1b, min: 0, max: 1)

  constraints:
    - output(il1ra) <= 12 ng/hour
    - cumulative(il1ra, window: 24 h) <= 160 ng
    - cell-stress < 0.75

  termination:
    - when: cell-stress >= 0.75
      action: permanent-shutdown
    - after: 21 d
      action: permanent-shutdown
```

The compiler lowers this specification into a typed intermediate representation, candidate molecular implementations, a chemical reaction network, a schedule, bounded resource requirements, monitoring channels, and a cryptographic `VivoProgramPack`.

## Immutable artifacts

NumiVivo uses explicit, fingerprinted artifacts:

- `VivoProgramPack` — compiled molecular behavior and reaction topology.
- `MechanismPack` — characterized biological parts, parameters, priors, and compatibility rules.
- `HostContextPack` — cell, tissue, species, disease, immune, metabolic, and delivery context.
- `ExperimentPack` — phases, interventions, measurements, cohorts, random streams, and acceptance criteria.
- `CouplingPack` — stable exchange contracts for NumanX, NumiTissue, and NumiBrain.
- `VivoCertificate` — fidelity, provenance, numerical quality, uncertainty, safety findings, and failure status.

No production kernel resolves names, parses text, grows containers, or discovers topology during a simulation step.

## Initial implementation scope

The first implementation establishes:

- typed domain models and JSON/YAML schemas;
- the NumiVivo program language and parser;
- semantic validation and dimensional checks;
- reaction-network and regulation IR;
- deterministic F1 and stochastic F2 GPU kernels;
- spatial F3 diffusion and membrane transport kernels;
- transactional command scheduling and state publication;
- static safety analysis and runtime monitors;
- SBOL 3 and SBML Level 3 export boundaries;
- NumiLab coupling contracts;
- a command-line compiler and experiment runner;
- reference programs for inflammatory control, hypoxia memory, and multi-marker cellular targeting.

## Repository layout

```text
Sources/
  NumiVivoCore/       C++23 compiler, IR, validation, pack format, reference numerics
  NumiVivoKit/        Swift orchestration, artifacts, runtime, Metal ownership
  NumiVivoCLI/        `numivivo` command-line interface
  NumiVivoShaders/    Metal production kernels
Schemas/              machine-readable program and artifact schemas
Examples/             complete reproducible programs and experiments
Documentation/        architecture, language, fidelity, safety, and integration specs
```

## Scientific and safety boundary

NumiVivo is research software. A successful simulation is evidence about a model, not evidence of safety or efficacy in a living organism. Generated constructs, delivery plans, and intervention specifications require independent biological review, appropriate containment, ethics oversight, and applicable regulatory authorization before physical use.

The software is designed to make assumptions and failure conditions explicit. It must not hide uncertainty behind a single score or convert an unverified simulation into an experimental instruction.

## Standards

NumiVivo is designed to interoperate with:

- SBOL 3 for machine-readable synthetic biological designs;
- SBML Level 3 for reaction-network exchange;
- SED-ML and COMBINE archives for simulation experiments;
- ontology-backed identifiers for molecular entities, processes, units, and evidence;
- FAIR provenance and content-addressed artifacts.

## Status

Active foundational development. The public interfaces and artifact formats will remain versioned while the implementation matures.

## License

Apache License 2.0. See `LICENSE`.
