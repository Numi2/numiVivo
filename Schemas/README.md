# NumiVivo schemas

The files under `Schemas/v1alpha1` define the human-readable source and manifest contracts currently implemented by NumiVivo.

## Authority

JSON Schema performs structural validation. The C++ or Swift semantic compiler remains authoritative for:

- unit compatibility and conversion;
- cross-reference resolution;
- graph topology;
- evidence and context requirements;
- numerical limits;
- safety findings;
- generated fingerprints;
- campaign ledgers;
- checkpoint signatures;
- device capability and memory limits.

Passing a schema does not imply that a model is executable, biologically valid, safe, or physically realizable.

## Source versus generated artifacts

Source documents may omit generated `fingerprint` fields when the corresponding semantic compiler creates them. Generated manifests and packs must include and verify their fingerprints.

Examples:

- `population-model.json` is source and may omit its model fingerprint.
- `campaign-definition.json` is source and has no campaign-manifest fingerprint.
- `campaign-manifest.json` is generated and must include candidate, job, ledger, and manifest fingerprints.
- checkpoint manifests are embedded inside the binary checkpoint envelope and are authenticated with the payload.

Do not insert arbitrary fingerprint strings merely to satisfy a generated-artifact schema. Compile the source document and retain the produced artifact.

## Current schemas

| File | Contract |
|---|---|
| `vivo-program.schema.json` | Molecular behavior, reaction, rule, constraint, and termination source |
| `mechanism-synthesis.schema.json` | Sequence-free mechanism problem or part library |
| `exact-ssa-model.schema.json` | Exact discrete stochastic source model |
| `hybrid-stochastic-plan.schema.json` | Generated reaction-authority partition |
| `population-model.schema.json` | Spatial phenotype-population source model |
| `physiological-partition.schema.json` | Reversible analyte partition source model |
| `numilab-coupling.schema.json` | Transactional participant-channel contract |
| `surrogate-contract.schema.json` | Bounded surrogate authority contract |
| `campaign.schema.json` | Campaign definition or generated manifest |
| `checkpoint-manifest.schema.json` | Manifest embedded in a checkpoint envelope |

## Versioning

A breaking change to field meaning, canonical ordering, default behavior, units, update semantics, or fingerprint construction requires a new schema or artifact major version.

Adding an optional field requires:

1. a documented default;
2. semantic compiler support;
3. canonical fingerprint behavior;
4. backward-reader behavior;
5. an example update.

Unknown required features are rejected. Readers do not silently interpret a newer major version as the current one.

## Limits

Schema collection limits mirror format capacities where practical, but the semantic compiler and runtime may impose stricter limits based on:

- configured resource policy;
- host address space;
- ProgramPack ABI widths;
- Metal buffer limits;
- device working-set size;
- algorithmic complexity;
- safety policy.

A schema maximum is not a promise that the maximum model fits a particular device.

## External validation

The schemas target JSON Schema draft 2020-12. External validation has not yet been performed as part of the current development pass. The first verification gate must check every source example against the schemas and then against the semantic compiler.
