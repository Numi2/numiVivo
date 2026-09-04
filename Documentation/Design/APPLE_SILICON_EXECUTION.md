# Apple-silicon execution architecture

## Scope

NumiVivo is implemented as an Apple-native molecular and tissue simulation system. This document defines the production execution contract. It is not a performance claim and does not substitute for device-specific measurement.

The runtime has three authorities:

- **C++23** owns parsing, semantic compilation, reaction-graph analysis, mechanism synthesis, immutable pack construction, and portable reference logic.
- **Swift 6** owns artifact lifetimes, structured concurrency, experiment scheduling, device selection, coupling transactions, checkpoints, and public APIs.
- **Metal** owns production state evolution, stochastic events, transport, spatial population dynamics, reductions, monitoring, and state migration.

Python may be used by external research workflows, but it is not part of a production step and cannot mutate authoritative runtime state.

## Resource ownership

### Authoritative state

Persistent state is allocated with private Metal storage. Each transactional runtime maintains at least two copies:

```text
committed state
candidate state
```

Integrators that require intermediate values also maintain stage state and derivative or event scratch. Delayed-reaction queues, random counters, temporal state, population fields, physiology state, and adaptive-fidelity topology are included in the transaction domain when their values can affect the next commit.

The CPU does not directly mutate private state. Inputs are uploaded through bounded staging buffers, and observations are copied into explicit readback buffers. Shared storage is restricted to:

- compact command records;
- status and failure records;
- sparse coupling updates;
- publication requests;
- event headers;
- explicit snapshots and certificates.

### Immutable tables

ProgramPack, mechanism, topology, stoichiometry, rate-law, cohort, transport, and monitor tables are compiled before execution. Production kernels receive numeric offsets and bounded counts. They do not resolve names, parse documents, perform ontology lookup, grow containers, or discover model topology.

Immutable tables should be shared across campaign jobs that have the same content fingerprint. Mutable state remains isolated per job or environment.

## Data layout

The primary molecular layout is species-major structure-of-arrays:

```text
index = species × laneCount + lane
```

This supports contiguous access when a kernel evaluates one species across many cells, voxels, environments, or campaign candidates. Reaction-centered kernels use precompiled stoichiometry ranges, while species-centered gather kernels use sparse incidence offsets.

Spatial population state uses:

```text
index = phenotype × voxelCount + voxel
```

Physiological partition state uses:

```text
index = analyte × compartmentCount + compartment
```

Layouts are explicit in runtime contracts and checkpoints. A consumer must not infer a layout from a file name or array length.

## Dispatch organization

Compilation groups reactions by rate law, stiffness class, gate presence, delay semantics, stochastic eligibility, and spatial requirements. Cohort descriptors contain bounded reaction ranges and preferred dispatch widths.

The runtime dispatches homogeneous cohorts to reduce control-flow divergence. Heterogeneous biological programs are therefore transformed into several numerical cohorts rather than executed through a single branch-heavy kernel.

Threadgroup width is selected from the compiled pipeline's execution width and maximum thread count. It is not hard-coded to one Apple GPU generation.

## Numerical authority

### F0 — logic

The runtime evaluates typed conditions, temporal operators, state transitions, outputs, monitors, and termination rules. No continuous reaction dynamics are implied.

### F1 — deterministic reactions

Well-mixed continuous state uses FP32 authoritative values and a two-stage explicit integrator. Compiler and runtime bounds determine whether a requested step is admissible. Negative or non-finite candidates are rejected or retried with a smaller step.

### F2 — stochastic reactions

Discrete state uses UInt32 molecular counts. Counter-based random streams are derived from artifact, experiment, replicate, lane, mechanism, logical-step, and sample identities. Dispatch order does not define stream identity.

Two stochastic paths are present:

- bounded tau-leaping for event-dense cohorts;
- exact direct-method SSA for low-count or critical cohorts.

Shared-species reaction components remain under one numerical authority per transaction.

### F3 — spatial transport

Reaction and transport are split through explicit transaction stages. Diffusion and upwind advection use precomputed geometry and boundary modes. The scheduler intersects reaction, diffusion, advection, intervention, measurement, checkpoint, delayed-event, and coupling step bounds.

### F4 — coupled tissue execution

NumiVivo coordinates with NumanX, NumiTissue, NumiBrain, physiological distribution, and external experiment participants through a two-phase protocol:

```text
prepare candidates
validate all candidates
stage commit for all participants
release publication for all participants
```

Any prepare or stage failure triggers rollback of every prepared participant. Candidate outputs are not observable as committed state before joint release.

## Adaptive fidelity and allocation

`VivoAdaptiveFidelityPlanner` chooses numerical authority, representation, chunk count, migration contract, and step limits. `VivoAdaptiveStateArena` applies that plan to real private Metal heaps.

For every cohort and chunk, the arena allocates:

- a bounded state region containing the declared number of state copies;
- random-counter storage;
- delayed-event storage;
- topology storage;
- auxiliary numerical scratch.

Heap size is computed from Metal size/alignment requirements, not by summing unaligned logical sizes. Allocation is rejected when the aligned heap exceeds the declared working-set budget.

Reconfiguration follows a replacement transaction:

1. Validate the new plan and runtime layouts.
2. Allocate and zero a complete replacement heap.
3. Convert or copy every preserved state region.
4. Preserve random counters and delayed queues only when declared.
5. Rebuild spatial or tissue topology when required.
6. Read the migration status.
7. Swap arena ownership only if every migration succeeded.

Continuous-to-discrete conversion rejects negative, non-finite, overflowing, or materially fractional counts. Discrete-to-continuous conversion reports integer values that cannot be represented exactly in binary32.

## Campaign execution

Campaign jobs are content-addressed and grouped by immutable topology. The batch planner uses the device working-set profile and reserves memory for the operating system, command resources, and other runtimes.

A batch may combine jobs only when:

- its immutable topology fingerprint matches;
- its numerical authority is compatible with policy;
- the aligned resident estimate remains inside the usable budget;
- the batch job limit is not exceeded.

Execution waves bound the sum of simultaneously resident batches. The campaign manifest fixes candidate ordering, replicate ordering, seeds, and the hash-chain ledger before execution begins.

## Surrogates

Surrogate execution is acceleration, not authority. The gate checks:

- model, training-data, mechanism, and context fingerprints;
- input shape and declared domain;
- normalized extrapolation;
- uncertainty output;
- finite values and output bounds;
- mandatory authoritative refresh interval.

Rejected predictions fall back to the authoritative mechanism. Accepted predictions still pass normal runtime monitors and cannot satisfy a safety certificate by themselves.

## Failure publication

Every production dispatch writes a compact status record. Failures are typed, include the first offending index when available, and retain the committed state.

The runtime may return:

- committed;
- committed with a reduced step;
- rejected;
- reversible shutdown;
- permanent shutdown.

A command-buffer completion status is necessary but not sufficient for commit. Numerical, biological, safety, capacity, coupling, and provenance checks must also pass.

## Current limitations

The source tree has not yet received a package-wide compiler or shader verification pass. Therefore API spelling, platform availability, Swift actor isolation, ABI layout, and Metal language compatibility remain unverified until the designated Apple-silicon integration environment performs that work.

Performance, accuracy, and memory claims must be produced by versioned benchmark campaigns on named hardware. This document defines intended execution behavior only.
