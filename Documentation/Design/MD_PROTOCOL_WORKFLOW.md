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

A stage report binds entry/exit checkpoints, optional configuration transition, trajectory manifest, observation tail, progress and any minimization/rejection certificate. Completed trajectories are sealed immutable views. The writer remains retryable until the terminal cursor and its reference have been published.

As soon as a minimizer or dynamics stage returns a terminal outcome, the runner retains that exact outcome and its certificates before any fallible checkpoint or report operation. Once exported, its accepted checkpoint is retained too. Finalization constructs a separate terminal cursor and publishes the checkpoint, sealed trajectory view and report before publishing that cursor and updating its reference. Only then does the runner release the outcome and its runtime. A report-write failure cannot turn an unsuccessful minimizer's new geometry into a running cursor with a fresh iteration budget, or cause a rejected dynamics candidate to be tried again by the same runner.

While a terminal outcome is pending, every persistence or `run()` retry attempts to publish that same outcome before any numerical work. If storage still fails, the receipt reports the previous confirmed durable cursor and includes a persistence diagnostic. A newly published reference can point only to a complete terminal cursor; it is acknowledged in the receipt only after reference publication succeeds. The unpublished accepted geometry may be lost if the process exits; resume then starts from the older durable boundary. The runner does not claim to have saved a new numerical state without its corresponding terminal gate result.

## Restart and failures

A restart verifies numerical profile, plan identity, initial-state identity, current/entry checkpoint identity, accepted-step progress, prior stage reports, and active trajectory/observation boundaries before allocating its device runtime. The verification walks the completed stage sequence: every report has the declared schema, stage, entry/exit checkpoints, clocks and output references; every transition is recomputed from the preceding accepted exit and the declared destination velocity policy. Physical coordinates, time and cell must survive stage entry. Dependent sites are reconstructed by the runtime; Maxwell-Boltzmann velocities remain the explicitly declared initialization operation, not a host-side distribution qualification.

Dynamics clock validation compares the declared step count and timestep with the persisted time using a conservative FP64 positive repeated-addition error bound. This accepts the runtime's accumulated rounding without replaying billions of additions and rejects clocks unrelated to the accepted step count. It does not claim bit-exact replay from a time formula.

Successful minimization reports must include a bound certificate with consistent iteration accounting and a force/convergence result that satisfies the declared gate. A rehashed report cannot turn an unsuccessful required minimization into a successful prior stage merely by setting `successful: true`. Successful dynamics reports require all requested steps, and blocked dynamics reports require a rejected candidate certificate bound to the accepted exit.

The minimizer re-evaluates accepted-state energy after rejected trials, so a resume check cannot demand bitwise monotonic readback from unordered PME charge accumulation. Energies must remain finite, and the declared force/convergence gate still controls success. A rejected Double step-scale shrink may underflow to zero; its blocked certificate remains resumable as blocked. Host regressions cover one FP32 ULP of re-evaluation difference and the exact scale-underflow arithmetic without claiming measured PME noise.

Resume forks a new UUID-named checkpoint reference; the original run's prefix is unchanged. The active stage has already been initialized. A finished stage transitions once into the next stage; a blocked numerical or convergence gate is not silently bypassed.

A publication failure can occur after the numerical step commits but before its scheduled sample or observation is linked. The v1 cursor remains sufficient to recover that boundary: trajectory frame counts/ranges and the observation tail identify which current outputs exist. A running prefix may omit only an output due at its current accepted checkpoint. Resume publishes each missing output from that checkpoint before taking another step or finishing the stage. An already appended frame is not duplicated when only its observation needs repair. Older gaps, excess outputs, and incomplete finished-stage output schedules are rejected because the saved checkpoint cannot reconstruct an earlier state. The same reconciliation applies when the caller retries `run()` on an existing runner after a transient write failure.

Trajectory restart validation checks the link index and final coordinate payload. Observation validation binds the tail, ordinal, schedule, physical clock and finite observables; it does not replay the simulation or exhaustively reread earlier observation payloads. These integrity checks do not establish trajectory accuracy or scientific equilibration.

The trajectory check uses streaming `restart` validation without a default length cutoff or an in-memory array of all chunk hashes. Its completed validation binding is passed directly to the active writer, so protocol resume does not traverse the same index twice. Long validation remains cancellable at each read. Persistence after cancellation deliberately uses the writer's already-validated flushed manifest through a bounded store read, preserving the ability to publish a pending terminal outcome without entering a new cancellable inspection traversal.

No adaptive timestep or stochastic retry occurs inside the runner. A rejected MD candidate or unsuccessful required minimization blocks the run and persists its report. Cooperative task cancellation saves the most recent exportable accepted boundary. If the device runtime is poisoned or storage fails, the receipt points to the last durable cursor and includes the failure instead of claiming newer state was saved. In-flight minimization is restartable from its saved accepted geometry, but line-search history is not serialized; this is explicitly a restarted minimizer, not exact optimizer continuation.

A process crash can lose work since the last durable cursor; the checkpoint interval controls that window. A named reference is an atomic filesystem pointer, not a distributed multi-writer lease. UUID references prevent accidental same-name writers in this runner. Checkpoint hashes establish integrity, not authorship or trajectory accuracy. Native recovery execution is qualified only to the scope of the named tests below.

`MDProtocolResumeTests` covers host-side gate/transition/clock tampering after hashes are recomputed, accumulated FP64 clock rounding, recoverable current-output boundaries, and rejection of older missing samples. Its native regressions resume a real Metal accepted checkpoint with missing/both/partial publications, compare every trajectory frame and final checkpoint with uninterrupted execution, check exact observation-step accounting, and cover recovery at the last step. A separate native sequence starts a minimizer, transfers into Maxwell-Boltzmann-initialized NVT, completes both stages and resumes the finished protocol; an unsuccessful required native minimization remains blocked after resume. Running that suite on the target hardware is required evidence; the presence of these tests alone is not a passing result.

The suite also forces an actual observation write failure by placing an empty directory at the precomputed first observation object's otherwise absent content address inside its temporary store. It checks the failed receipt and accepted durable checkpoint, removes that exact test-owned sentinel, and verifies both same-runner retry and forked resume against uninterrupted coordinates, velocities and observations. No production storage failure hook is involved.

A second real I/O regression blocks the exact report object for a required, unsuccessful one-iteration minimization that has changed accepted geometry. Repeated same-runner failures must leave the old durable entry unchanged. A fork from that entry repeats only the original budget. After removing the test-owned report sentinel, both pending runners must publish the original unsuccessful certificate and exact original exit checkpoint, and further resume must remain blocked.

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

The [7 September 2026 reliability audit](../Audit/2026-09-07_COMPLETION_RELIABILITY.md) records a clean full build and nine passing MD recovery regressions, including real Metal simulation, complete checkpoint/trajectory comparisons and storage failures. The broader portable/native numerical contracts, force agreement, PME grid convergence, constrained trajectory statistics, NPT detailed balance, throughput and resource-pressure behavior still require independent qualification. The current MD nearest-image profile is restricted to orthogonal cells; general triclinic structure storage must not be mistaken for unrestricted triclinic dynamics.
