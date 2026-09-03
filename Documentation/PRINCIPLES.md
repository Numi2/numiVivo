# NumiVivo engineering principles

## Biological ground truth is authoritative

The simulator exists to make biological hypotheses explicit and computable. Internal consistency, numerical convergence, and agreement between Numi components are necessary quality checks; none of them supersedes experimental evidence.

Each parameter and mechanism carries an evidence class:

- `observed` — directly measured in the declared context;
- `derived` — calculated from declared observations;
- `calibrated` — fitted against declared data;
- `inferred` — estimated from a statistical or learned model;
- `assumed` — selected by the model author;
- `hypothetical` — introduced to test an unverified mechanism.

The runtime propagates these classes into result certificates.

## Context is part of the program

Molecular behavior depends on cell type, tissue state, species, developmental stage, disease state, delivery mode, metabolism, immune activity, mechanics, and time. A program without a compatible `HostContextPack` is incomplete and cannot be promoted beyond F0.

## Compilation removes ambiguity

Human-readable programs are not executed directly. Compilation resolves:

- identifiers and ontology references;
- dimensions and unit conversions;
- mechanism implementations;
- reaction topology and stoichiometry;
- program-state layout;
- monitor and termination semantics;
- cell and compartment cohorts;
- bounded capacities and dispatch dimensions;
- random-stream namespaces;
- coupling channels;
- artifact fingerprints.

The production runtime consumes immutable numeric tables.

## Transactional execution

For logical step `k`, the runtime reads committed state `S_k` and writes only to shadow state `S*`. It evaluates numerical validity, physical bounds, program constraints, runtime monitors, and coupling acceptance before publishing `S_{k+1}`.

```text
S_k → shadow execution → validation → commit S_{k+1}
                                  └→ reject and retain S_k
```

A failure is typed and local to the smallest safe transaction domain. Partial molecular state is never published.

## Fidelity is a contract

Fidelity levels declare represented phenomena, numerical methods, unresolved effects, and valid conclusions. Promotion requires a stated error measure and a certificate. Demotion is allowed only when the requested observable remains within its declared tolerance.

## Bounded execution

Programs compile to explicit limits for cells, species, reactions, regulations, events, monitors, spatial voxels, delayed transitions, and output records. A runtime may reject or substep a request; it may not grow unbounded data structures on the GPU timeline.

## Reproducible stochasticity

Random values are generated from counter-based streams. Stream identity is independent of thread scheduling and includes:

```text
program fingerprint
experiment fingerprint
replicate seed
cell or voxel identity
mechanism identity
logical step
sample lane
```

## Apple-native data ownership

Persistent authoritative state uses private Metal resources. Shared memory is reserved for compact command packets, diagnostics, and publication boundaries. CPU and GPU share an address space on Apple silicon, but shared accessibility is not used as a reason to expose the entire hot state to CPU mutation.

## Safety is not a score

Safety analysis produces findings, assumptions, unresolved hazards, counterexamples, monitor coverage, and termination guarantees. A scalar ranking may summarize but cannot replace these records. Missing evidence remains missing evidence.

## No hidden fallback

Unsupported kinetics, units, mechanisms, standards features, or platform capabilities produce explicit diagnostics. The compiler and runtime do not silently substitute a simpler biological model.
