# Contributing to NumiVivo

NumiVivo combines molecular biotechnology, numerical simulation, compiler design, and Apple-native GPU systems. Contributions must preserve scientific traceability and deterministic artifact semantics.

## Design requirements

Every change must satisfy the following constraints:

- Biological assumptions are explicit and attached to evidence or a declared modelling hypothesis.
- Fidelity is named. F0, F1, F2, F3, and F4 results are not interchangeable.
- Authoritative state uses FP32 or a more accurate representation. Reduced precision is restricted to derived caches, learned approximations, and preconditioners with an explicit error boundary.
- Runtime topology, capacities, names, units, and semantic references are compiled before production dispatch.
- Production simulation does not depend on Python callbacks or per-step CPU interpretation.
- GPU steps are transactional. Failed validation cannot publish partial state.
- Random streams are reproducible from artifact fingerprint, experiment seed, cell identifier, mechanism identifier, and logical step.
- New file formats are versioned, little-endian, content-addressed, and reject unknown required features.
- Safety findings cannot be downgraded silently.
- A simulation result must retain the exact program, context, mechanism library, runtime configuration, and source fingerprints that produced it.

## Source organization

- `NumiVivoCore` contains C++23 compilation, graph analysis, immutable formats, and portable reference logic.
- `NumiVivoKit` contains Swift orchestration, artifact ownership, experiment scheduling, diagnostics, and public APIs.
- `NumiVivoShaders` contains Metal kernels and only the minimum Swift code required to expose packaged shader resources.
- `NumiVivoCLI` contains the command-line surface. Domain logic does not belong in the CLI.

## Change procedure

1. Define the affected semantic contract.
2. Update the schema or format version when serialized meaning changes.
3. Implement the compiler/runtime change.
4. Add diagnostics for invalid or unsupported states.
5. Update documentation and one complete example.
6. Record any numerical, biological, or platform limitation that remains.

## Commit conventions

Use focused conventional commits:

```text
feat(core): lower cooperative binding into reaction IR
feat(metal): add bounded stochastic reaction dispatch
fix(pack): reject overlapping section ranges
docs(safety): define irreversible shutdown evidence
```

Generated artifacts, local experimental data, and compiled Metal libraries are not committed unless they are an intentional release fixture.

## Scientific claims

Code comments and documentation must distinguish:

- experimentally observed facts;
- values imported from a named dataset;
- calibrated parameters;
- inferred parameters;
- design assumptions;
- unvalidated hypotheses.

Do not describe model agreement as biological validation. Do not include experimental instructions that exceed the repository's computational research scope.
