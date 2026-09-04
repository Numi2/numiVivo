# Artifacts and provenance

## Objective

Every NumiVivo result must be traceable to the exact source, compiled program, mechanism assumptions, host context, experiment definition, coupling contract, runtime configuration, device profile, random streams, and committed state that produced it.

Provenance is represented by immutable content identities and hash-linked execution records. File names and directory locations are not identities.

## Artifact classes

### VivoProgram source

The human-readable behavior specification. It contains target context, inputs, species, state, parameters, reactions, rules, constraints, and termination semantics.

### ProgramPack

The compiled binary program consumed by production runtimes. It contains fixed-layout strings, species, parameters, stoichiometry, reactions, expression bytecode, actions, rules, monitors, cohorts, sparse incidence, and a runtime contract.

### MechanismPack

The selected abstract mechanisms, evidence records, parameter priors, compatibility assumptions, and performance envelopes that realize a program behavior in the model.

### HostContextPack

The organism, tissue, cell, disease, delivery, immune, metabolic, transport, and evidence context against which a program is interpreted.

### ExperimentPack

The scheduled interventions, measurements, stop conditions, candidate parameters, replicate identities, seeds, checkpoint cadence, and acceptance criteria.

### CouplingPack

The participant, channel, unit, shape, bounds, update semantics, cadence, and transactional ownership contract for NumanX, NumiTissue, NumiBrain, NumiVivo, physiology, environment, and experiment participants.

### CampaignManifest

The fully expanded candidate and replicate schedule. It fixes job ordering, seeds, candidate identities, and a hash-chain ledger before execution begins.

### VivoCertificate

A result record containing fidelity, runtime disposition, numerical status, safety findings, uncertainty, device identity, artifact identities, and committed time.

### Checkpoint envelope

A binary package containing a committed runtime state and the manifest needed to validate and resume it.

## Fingerprints

SHA-256 is used as the current content-identity function. A fingerprint is encoded as 64 lower-case hexadecimal characters.

Fingerprints are applied at several levels:

- original source bytes;
- normalized typed artifacts;
- individual binary sections;
- complete ProgramPack content;
- campaign candidates;
- campaign jobs;
- campaign manifest;
- checkpoint authenticated content;
- adaptive arena layout;
- mechanism candidate assignment.

A fingerprint binds bytes or a precisely defined canonical encoding. It is not a semantic equivalence claim.

## Canonical JSON

Typed JSON artifacts use UTF-8 encoding with sorted object keys and unescaped slashes for internal fingerprint generation. Arrays preserve their declared order unless the artifact contract explicitly sorts them before encoding.

Semantically unordered sets are normalized before hashing. Examples include:

- candidate parameter assignments sorted by parameter identifier;
- mechanism slot assignments sorted by slot identifier;
- reaction and species membership sorted numerically;
- metadata encoded with sorted keys.

A format must never change its canonicalization rules without a version change.

## ProgramPack integrity

ProgramPack has:

- magic and major/minor version;
- compiler ABI version;
- fidelity and feature flags;
- section descriptors;
- source fingerprint;
- content fingerprint;
- per-section fingerprints;
- explicit offsets, sizes, strides, counts, and alignments.

The reader rejects:

- unsupported major versions;
- unknown required sections;
- duplicate section types;
- malformed stride or count;
- invalid alignment;
- integer overflow;
- ranges outside the file;
- header overlap;
- section overlap;
- per-section fingerprint mismatch;
- complete-content fingerprint mismatch;
- missing required sections.

## Campaign ledger

Campaign compilation creates a deterministic candidate list and then a deterministic replicate list.

For job `j`, the ledger digest binds:

```text
ledger format tag
previous job digest
job index
candidate fingerprint
replicate index
derived seed
```

The first job uses a 64-zero previous digest. `ledgerHead` is the final job digest.

The validated artifact loader verifies:

- contiguous candidate and job indices;
- unique candidate fingerprints;
- finite and unique parameter assignments;
- candidate fingerprints;
- candidate reference from every job;
- previous-digest chain;
- every job digest;
- final ledger head;
- complete manifest fingerprint.

A campaign runner must refuse a manifest whose ledger does not verify.

## Checkpoint envelope

The current checkpoint format contains:

```text
fixed header
canonical JSON manifest
binary payload sections
optional Ed25519 signature
```

The header authenticates the manifest and payload together with SHA-256. An optional Ed25519 signature is computed over that digest.

Every section descriptor contains:

- identifier;
- encoding;
- payload offset and length;
- element count and stride;
- independent fingerprint.

The decoder validates total length, manifest length, payload length, signature length, section count, section ranges, overlap, layout, and fingerprints before exposing any section.

## Checkpoint lineage

A checkpoint may name a parent checkpoint fingerprint. This creates an explicit lineage:

```text
checkpoint_0 -> checkpoint_1 -> checkpoint_2
```

A branch is represented by two children naming the same parent. It is not represented by overwriting a previous checkpoint.

Checkpoint metadata also binds:

- runtime kind;
- owning artifact;
- source artifact;
- host context;
- experiment;
- coupling contract;
- logical step and time;
- random-stream version;
- signing-key identity.

## Required checkpoint state

A checkpoint is complete only when it captures all state needed for deterministic continuation. Depending on runtime, this includes:

- every authoritative state buffer;
- temporal state;
- rule refractory state;
- random counters;
- exact-event counters;
- delayed-event queues and due times;
- population fields;
- physiology state;
- adaptive-fidelity representation and topology;
- experiment intervention cursor;
- measurement and checkpoint cursors;
- coupling ledger identity;
- shutdown lifecycle state.

A state-only array is a snapshot, not necessarily a resumable checkpoint.

## Signatures

A signature validates that the holder of a named private key signed the authenticated checkpoint digest. It does not establish scientific validity, biological safety, regulatory approval, or trustworthy source data.

Signing-key rotation creates a new key identifier. Verification policy determines whether unsigned checkpoints are accepted.

## Bounded loading

`VivoValidatedArtifactLoader` checks a JSON artifact before typed decoding:

- regular-file requirement;
- optional symbolic-link rejection;
- maximum byte count;
- maximum nesting depth;
- maximum node count;
- maximum string and key bytes;
- maximum array size;
- maximum object-member count.

Typed validators then enforce schema version, identities, finite values, model bounds, generated fingerprints, and ledger integrity.

The loader does not execute embedded content, follow URLs, load dynamic libraries, invoke a shell, or resolve external code.

## Artifact store

The content-addressed store uses the artifact fingerprint as the storage identity. Writes are atomic:

1. Encode the complete artifact.
2. Calculate and verify its fingerprint.
3. Write to a temporary file in the destination filesystem.
4. Flush and close the file.
5. Rename into the content-addressed path.
6. Verify the stored bytes before returning success.

An existing object with the same identity must contain identical bytes. It is never overwritten with different content.

## Run bundles

A run bundle should contain:

- all source documents;
- all compiled packs;
- campaign manifest;
- device and runtime profile;
- step certificates;
- event stream;
- measurement results;
- checkpoints;
- safety report;
- uncertainty report;
- export artifacts;
- one top-level manifest that fingerprints every member.

The bundle must remain interpretable without access to mutable external databases. External evidence references can remain links, but the values actually used by the run must be included or content-addressed.

## Redaction and sensitive data

Artifacts must not include patient-identifying data, raw clinical records, personal genomic data, access tokens, private keys, or proprietary sequence libraries unless an explicitly authorized storage and privacy system is used.

Diagnostics should use bounded identifiers and indices. They must not dump arbitrary source documents or state buffers by default.

## Scientific interpretation

Integrity proves that bytes and execution lineage are consistent. It does not prove that:

- assumptions are correct;
- parameters are identifiable;
- the represented biology is complete;
- a model extrapolates to a new host;
- a result is safe or effective.

Those limitations remain part of the certificate and cannot be removed by signing an artifact.
