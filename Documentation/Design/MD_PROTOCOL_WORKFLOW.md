# Staged molecular dynamics protocols

## One execution and persistence authority

VivoMDProtocolRunner drives the existing VivoMDMetalRuntime; it does not implement another integrator, force field, random stream or checkpoint arena. Large sampled coordinate arrays use VivoMDTrajectoryArchiveWriter over VivoArtifactStore. Numerical checkpoints, stage reports, transitions, observations and restart cursors are verified objects in that same store.

The supplied template is minimization -> NVT -> NPT -> NPT production. Its step counts and target conditions are illustrative inputs, not the CovAngelo paper's reproduced parameters or evidence of equilibration. The template requires minimizer convergence before proceeding. Protocols can declare other finite sequences of the supported numerical profiles.

## Stage entry

Every stage has one full MD configuration and explicit velocity policy. A source initial-state file has no velocities, so the first stage must say zero or maxwellBoltzmann. Later stages preserve velocities by default. Maxwell-Boltzmann initialization needs an explicit independent seed and target temperature; it projects sampled velocities to the constraint tangent space without forcing a chosen realized kinetic energy. Initialization is persisted before stepping, so resumption of a running stage does not repeat it.

Changing configuration is not a checkpoint-identity bypass. VivoMDStageTransfer emits an immutable transition record binding source checkpoint/configuration, destination configuration, declared velocity operation and prepared destination state. Global accepted step, physical time, cell and positions are retained. A post-initialization entry checkpoint completes the transition. The old runtime is released before allocating the new arena.

## Execution and sampling

Stage-local completed-step counters control sampling and persistence schedules and survive restart. Sampled positions and optional observables are obtained under one runtime reservation. Observation-only ticks are stored even when no positions are sampled. Position samples are compact FP32 binary chunks; observation links contain explicit potential, kinetic and total energy, temperature and assumed constrained degrees of freedom. No sampling interval is interpreted as a scientific independence claim.

Each checkpoint interval persists an accepted MD checkpoint, flushes the trajectory prefix, then writes an immutable protocol cursor. Only after those objects exist does the named reference point to that cursor. A write failure may leave unreachable objects but cannot invalidate the prior durable cursor. There is no ever-growing in-memory trajectory array. Completed-stage references are bounded by the 128-stage plan limit. Observation and chunk indices are linked content-addressed sequences rather than repeatedly rewritten lists.

A stage report binds entry/exit checkpoints, optional configuration transition, trajectory manifest, observation tail, progress and any minimization/rejection certificate. Completed trajectories are sealed immutable views. The writer is not irreversibly sealed before the report is stored, so a failed finalization can still publish a consistent running prefix.

## Restart and failures

A restart verifies numerical profile, plan identity, current/entry checkpoint identity, accepted-step progress, prior stage reports, and active trajectory/observation boundaries. It forks a new UUID-named checkpoint reference; the original run's prefix is unchanged. The active stage has already been initialized. A finished stage transitions once into the next stage; a blocked numerical or convergence gate is not silently bypassed.

No adaptive timestep or stochastic retry occurs inside the runner. A rejected MD candidate or unsuccessful required minimization blocks the run and persists its report. Cooperative task cancellation saves the most recent exportable accepted boundary. If the device runtime is poisoned or storage fails, the receipt points to the last durable cursor and includes the failure instead of claiming newer state was saved. In-flight minimization is restartable from its saved accepted geometry, but line-search history is not serialized; this is explicitly a restarted minimizer, not exact optimizer continuation.

A process crash can lose work since the last durable cursor; the checkpoint interval controls that window. A named reference is an atomic filesystem pointer, not a distributed multi-writer lease. UUID references prevent accidental same-name writers in this runner. Checkpoint hashes establish integrity, not authorship or trajectory accuracy. The numerical profile is source-implemented, not GPU-qualified.

## User-side commands

```sh
# After building numivivo on Apple silicon:
numivivo md-protocol-template system.json > protocol.json
# Review/edit the explicit durations, conditions and minimization gate.
numivivo md-protocol-validate protocol.json --system system.json --state initial.json
numivivo md-protocol-run protocol.json --system system.json --state initial.json \
  --store ./md-artifacts > receipt.json

# Use the checkpointReference value printed in receipt.json or stderr.
numivivo md-protocol-resume protocol.json --system system.json --store ./md-artifacts \
  --reference md-RUN-UUID-checkpoint > resumed-receipt.json
numivivo md-protocol-inspect --store ./md-artifacts --reference md-RUN-UUID-checkpoint

# trajectoryManifest hashes are stored in stage reports and protocol cursors.
numivivo md-trajectory-inspect --store ./md-artifacts --manifest SHA256 --verify
```

`--checkpoint SHA256` can replace `--reference` for cursor inspection/resume. JSON on stdout is small metadata except the template; scientific state remains in the rooted store. Legacy md-run/md-minimize remain available for compatibility, but md-protocol-run is the bounded-output path for multi-stage work.

## Deferred gates

This pass did not run a Swift package build, Metal compile, MD simulation, checkpoint trajectory comparison, statistical ensemble test or throughput benchmark. The portable/native numerical contracts, force agreement, PME grid convergence, constrained trajectory statistics, NPT detailed balance and resource-pressure behavior require independent qualification. The current MD nearest-image profile is restricted to orthogonal cells; general triclinic structure storage must not be mistaken for unrestricted triclinic dynamics.
