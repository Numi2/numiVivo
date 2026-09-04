# Transactional runtime semantics

## Purpose

Molecular, population, physiological, and coupled simulations must not expose partially advanced state. NumiVivo therefore treats every logical step as an atomic transaction over all state that can influence future behavior.

The minimum transaction domain includes:

- continuous or discrete molecular state;
- program memory and temporal operators;
- random-stream counters;
- delayed-event queues;
- spatial population fields;
- physiological compartment state;
- coupling inputs and candidate publications;
- monitor and termination state;
- experiment cursor and scheduled-intervention cursor;
- checkpoint lineage metadata.

## State model

For committed step `k`, let:

```text
S_k = {
  molecular state,
  temporal state,
  stochastic counters,
  delayed queues,
  spatial fields,
  physiological state,
  experiment state
}
```

A step proposal creates candidate state `S*` from `S_k` and an immutable command `C_k`:

```text
S* = Advance(S_k, C_k)
```

Candidate execution cannot alter `S_k`. A validation function produces a typed decision:

```text
Validate(S_k, S*, C_k) -> Commit | Retry(dt') | Reject | Shutdown
```

Only `Commit` publishes `S*` as `S_{k+1}`.

## Command immutability

A command fixes:

- transaction identifier;
- logical step and substep indices;
- source and artifact fingerprints;
- requested and candidate time step;
- absolute target time;
- fidelity and numerical mode;
- random-stream namespace;
- input and coupling ranges;
- intervention and measurement boundaries;
- publication requests;
- bounded capacities.

Kernels may read the command but cannot modify it. Retries create a new command with the same logical-step identity and a different substep index and time step. Random draws are keyed so a retry does not accidentally reuse an unrelated stream.

## Single-runtime step

The authoritative sequence is:

1. Clear the candidate status and event counters.
2. Copy or reference committed state as the immutable transaction base.
3. Apply bounded external inputs to the base representation.
4. Execute reaction, stochastic, population, transport, or physiology stages.
5. Mature delayed events due at or before the candidate boundary.
6. Execute behavioral rules and state transitions.
7. Evaluate constraints and termination conditions.
8. Validate finite values, bounds, integer invariants, capacity, and conservation checks.
9. Produce candidate publications and events.
10. Complete the command buffer.
11. Inspect runtime status.
12. Commit by swapping state ownership, or reject without publication.

Candidate readback is not interpreted as committed output before step acceptance.

## Retry semantics

A retry is permitted only when:

- the reported failure is classified as step-size-dependent;
- the caller permits adaptive reduction;
- the reduced step remains above the configured minimum;
- the attempt count remains below the configured maximum;
- no permanent failure, capacity overflow, or terminal shutdown was reported.

The default reduction is a factor of two. A runtime may use a more specific accepted bound when it is derived from a declared reaction, diffusion, advection, delayed-event, or coupling limit.

Retries always start from the same committed state. Candidate state, candidate queues, and candidate publications from the rejected attempt are discarded.

## Delayed events

A delayed reaction contains two separate transitions:

1. Consume reactants and schedule a signed event extent.
2. At the due boundary, apply the product-side state change.

The committed delayed queue and candidate delayed queue are separate. A rejected transaction cannot append to the committed queue or mature committed events. Queue overflow is a hard rejection because dropping a delayed event would change model semantics.

## Stochastic state

Random-stream counters are part of the checkpoint and migration contract. The key namespace includes sufficient identity to make event generation independent of dispatch order.

Exact SSA and tau-leap authorities cannot independently update reactions that share dynamic species. The hybrid planner first partitions the reaction hypergraph into connected components, then assigns one authority to each component.

When a cohort changes representation:

- discrete to continuous preserves the represented count subject to FP32 exactness checks;
- continuous to discrete requires a non-negative value sufficiently close to an integer;
- random counters remain stable when the migration contract requests preservation;
- delayed queues are preserved only when event semantics are unchanged.

## Coupled transaction

A coupled step contains several participants and a coordinator. Each participant implements:

```text
prepare(transaction) -> candidate state and publications
stageCommit(transactionID)
releaseCommit(transactionID)
rollback(transactionID)
```

The coordinator executes:

1. Validate the channel contract and transaction identity.
2. Ask every participant to prepare.
3. Intersect participant step bounds.
4. Repeat prepare at a smaller common step when permitted.
5. Validate required input and output channels.
6. Validate publication shape, units, bounds, and fingerprints.
7. Ask every participant to stage commit.
8. Release every staged commit.
9. Append one coordinator-ledger record.

A prepare or stage failure triggers rollback for every participant that may hold candidate state. Release must be idempotent because a process failure can occur after some participants have published.

## Fixed-point coupling

Strongly coupled molecular and physiological systems may require candidate iteration. A fixed-point transaction repeats candidate exchange without committing:

```text
physiology candidate -> molecular exposure
molecular candidate -> secretion or uptake
repeat until residual <= tolerance
```

Iteration has explicit maximum count, norm, relaxation, and failure policy. Non-convergence results in retry or rejection. It does not publish the last iterate as a successful state.

## Safety responses

Runtime monitors map to one of these responses:

- `record` — preserve the finding without changing state acceptance;
- `clamp` — accepted only for monitors with explicitly implemented clamp semantics;
- `reject-step` — discard candidate state;
- `substep` — retry from committed state with a smaller step;
- `reversible-shutdown` — stop execution while preserving a restart path;
- `permanent-shutdown` — terminal runtime state for the current artifact lineage.

Unsupported response semantics are treated as rejection, not ignored.

## Events and publications

Events are append-only records with bounded capacity. An event record includes the transaction, logical step, lane or spatial index, mechanism or monitor identifier, severity, and candidate time.

Publications are selected numeric projections of candidate state. Their descriptors define:

- producer and consumer;
- unit;
- tensor shape;
- update semantics;
- lower and upper bounds;
- whether the channel is required;
- schema version.

A publication that does not match its contract prevents joint commit.

## Checkpoints

A checkpoint is valid only at a committed boundary. The checkpoint manifest records:

- runtime kind;
- artifact, source, host, experiment, and coupling fingerprints;
- parent checkpoint identity;
- logical step and time;
- random-stream version;
- section layouts and fingerprints;
- optional signing-key identity.

The binary envelope authenticates the manifest and payload together. Each section also has an independent fingerprint. Section ranges must be bounded and non-overlapping.

## Recovery

Recovery loads the entire committed transaction domain before execution resumes. Missing random counters, delayed queues, topology state, intervention cursors, or coupling lineage make a checkpoint incomplete for deterministic continuation.

A runtime must reject:

- unsupported checkpoint major versions;
- fingerprint mismatch;
- invalid or missing required signature;
- artifact-context mismatch;
- malformed section layout;
- a random-stream version it cannot reproduce;
- a representation it cannot migrate explicitly.

## Invariants

The following conditions hold at every externally observable boundary:

1. State is either the previous committed state or the complete next committed state.
2. A rejected step leaves authoritative state unchanged.
3. A committed step has exactly one certificate and ledger identity.
4. Candidate publications become visible only with their owning commit.
5. All participants in a coupled step advance to the same logical time.
6. Random and delayed-event state advance atomically with molecular state.
7. A checkpoint refers only to committed state.
8. Safety findings cannot be removed by a numerical retry.
9. Surrogate output cannot bypass monitors or become safety authority.
10. Unknown required features cause explicit rejection.
