# Wave B: staged execution and bounded trajectory output

## Status

This continuation adds source implementations for a staged MD protocol runner,
explicit stage-transfer records, thermal velocity initialization, protocol-level
checkpoint resumption, streamed binary trajectory archives, and CLI commands.
It also removes a read/write race in virtual-site force redistribution.

These additions reuse `VivoMDMetalRuntime`; they do not add another integrator,
force-field interpretation, constraint solver, or barostat. A successful protocol
result means its configured stages completed. It does not certify equilibrium,
ensemble sampling, force accuracy, or biological validity.

No Apple package build, Metal compile, MD simulation, PME force comparison,
constraint qualification, NVT/NPT distribution study or benchmark was performed
as part of this source implementation. The existing Wave B numerical backend
still requires those checks.

## Public interfaces

- `VivoMDProtocol`, `VivoMDProtocolStage`: explicit minimization and dynamics stages.
- `VivoMDProtocolRunner.run`: schedules the existing runtime and emits typed events.
- `VivoMDStageTransfer`: records a configuration-changing boundary without calling
  checkpoint restore under an incompatible configuration.
- `VivoMDVelocityInitializer`: one-time thermal initialization at fixed positions.
- `VivoMDProtocolCheckpoint`: exact protocol/stage progress plus a runtime checkpoint.
- `VivoMDTrajectoryArchiveWriter` and `VivoMDTrajectoryArchiveReader`: bounded,
  checksummed FP32 trajectory export and sequential verification.
- `VivoMDAtomicFileExport`: bounded regular-file reads and atomic explicit exports.

The CLI adds `md-protocol-template`, `md-protocol-check`, `md-protocol-run`,
`md-trajectory-inspect`, and `md-protocol-help`. Existing `md-run`, `md-minimize`,
AMBER, structure, and reaction-runtime commands retain their distinct routes.

## Stage semantics

The convenience factory produces minimization, NVT equilibration, NPT
 equilibration, and NPT production stages. The caller supplies step counts.
These counts are a schedule, not a declaration that a particular system has
converged or equilibrated. Required minimization convergence is enabled by default;
a failed convergence criterion stops the protocol instead of silently proceeding.

A stage boundary carries positions, velocities, periodic cell and physical time.
Its new configuration has its own fingerprint. Dynamics step numbering starts
at zero in a newly constructed stage runtime. Thermostatted stages must have
distinct seeds, so restarting local step numbering does not repeat an earlier
stage's random namespace. The convenience factory derives these seeds from the
base seed and stores them explicitly in the resulting protocol.

The factory thermalizes once after minimization. Later stages preserve velocities.
The initialization recipe is not reapplied during checkpoint resumption. Custom
protocols that begin with dynamics must provide either an initialization recipe
or an explicit initial velocity array through the SDK.

Each stage exclusively owns its runtime. It releases that runtime before the
next stage is constructed, rather than keeping every stage's GPU arena alive.
The small stage history does not retain the complete trajectory. Full state is
copied to the host only for requested snapshots and checkpoint/stage boundaries.

## Thermal initialization

Velocities are sampled with variance `R*T/m` in nm/ps units. A separate seeded
stream is used for preparation. Optional center-of-mass momentum removal applies
only at initialization; it is not a continuous COM-removal operation.

At fixed coordinates, distance-constraint velocity rows are projected in the
mass metric until a normalized residual reaches the requested threshold. This
is a one-time Double-precision preparation operation. Virtual-site velocities
are derived from the physical parent velocities. There is no final rescaling of
kinetic energy to force an exact instantaneous temperature.

Preparation-time periodic displacement uses a bounded closest-lattice-image
search. This fixes image selection for this preparation routine only. It does
not change the minimum-image implementation in the MD force/constraint shaders.
Highly ill-conditioned periodic cells exceed the preparation search budget and
are rejected rather than assigned an unverified nearest image.

## Resumption

A protocol checkpoint identifies the exact protocol, classical system, current
stage, configuration, stage-local accepted step, stage start time, completed
history and current runtime state. It cannot be used to change the model,
protocol or thermostat settings without an explicit new-stage operation.

Checkpoints are emitted at stage entry, configured accepted-step intervals and
stage boundaries. On a recoverable execution error the runner attempts to
publish the most recent accepted-boundary state. If the runtime cannot provide
that state, the preceding successfully published checkpoint remains the recovery
point. A failed output write does not imply rolling back already accepted MD.

A restored completed stage is not rerun. Resuming after minimization does not
repeat thermal initialization already represented by a dynamics checkpoint.
Within a dynamics stage, accepted step and seed come from the exact runtime
checkpoint. The numerical reproducibility guarantees are those of that runtime;
this orchestration layer does not establish cross-device bitwise equality.

## Output ownership and limits

Protocol-run CLI output paths are preflighted before GPU construction. Inputs
cannot also be output paths. A report, checkpoint and trajectory have distinct
destinations. First publication honors `--force`; subsequent interval checkpoints
replace the file owned by that run atomically.

The runner emits observations independently of trajectory sampling. An
observation is not discarded merely because a coordinate frame was not due at
the same step. Observation metadata has a declared maximum count. Coordinate
frames can be streamed without an in-memory array of all sampled coordinates.

A resumed invocation writes a new trajectory segment. A boundary frame may
appear in both segments; frame ordinals are local to each archive, and each
frame carries its stage identifier, configuration identity, local step and
physical time. The code does not append to a previously finalized archive.

## Binary archive v1

All fixed-width integers are little endian. All positions and optional velocities
are packed triples of FP32 scalars, in nm and nm/ps respectively. The writer
requires lossless FP32 representability; it does not silently quantize arbitrary
Double-precision inputs. Runtime snapshots are already FP32-authoritative.

Header:

```
8 bytes  ASCII NVMDTRJ1
4 bytes  UInt32 canonical JSON header length
N bytes  canonical JSON header
32 bytes SHA-256(header JSON)
```

Each frame:

```
4 bytes  ASCII FRM1
8 bytes  UInt64 payload length
payload:
    4 bytes  UInt32 canonical JSON frame-metadata length
    N bytes  frame metadata: ordinal, stage, config, local step, time, cell
    12 * particleCount bytes  positions
    12 * particleCount bytes  velocities, only when declared by the header
32 bytes SHA-256(tag || payload length || payload)
```

The record chain starts at the header digest. Each frame updates it as
`SHA-256(previousChain || frameDigest)`.

Footer:

```
4 bytes  ASCII END1
8 bytes  UInt64 frame count
32 bytes final record-chain digest
32 bytes SHA-256(tag || frame count || chain)
```

A reader reaches successful EOF only after verifying that footer and confirming
there are no trailing bytes. It enforces header/metadata/frame bounds before
allocation, verifies each frame checksum, checks record ordinals and finite
coordinates, and rejects truncation. Returning one verified frame is not a claim
that the entire archive has been consumed or verified. `verifiedReceipt()` is
available only after footer verification.

The writer uses an exclusive temporary file and syncs complete content before
atomic publication. No-clobber publication cannot replace an existing directory
entry. Explicit overwrite replaces that entry rather than following a final
symlink. An unfinished archive is not published under the destination filename.
The receipt includes a SHA-256 digest of all published bytes. These hashes detect
corruption and identify content; they are not digital signatures or authorization.

A trajectory is not a checkpoint. It may omit velocities and contains no promise
of complete thermostat/barostat/random-state restoration. Use the protocol or
runtime checkpoint for resumption.

## Source-audit correction

`nvivo_md_redistribute_virtual_force` formerly wrote `+= 0` for particles with no
parent-incidence entries. A virtual-site thread therefore wrote a force record
that parent threads were simultaneously reading. The kernel now returns before
any such write. Source site records remain read-only during the gather, while
only physical parent forces are updated. Site energy remains assigned once and
is not duplicated onto parents.

## Remaining numerical/source review boundaries

The following are not established by this protocol layer:

- General triclinic nearest-image correctness in every MD kernel; preparation's
  closest-image helper must not be mistaken for a backend-wide replacement.
- A total PME force-error bound implied by `pmeTolerance`. Its present real-space
  planning criterion and mesh spacing require independent force/energy checks.
- Constraint-tangent minimization convergence, complete long-range exception
  treatment, cell-dependent PME discretization across NPT proposals, and robust
  molecular unwrapping for large periodic components.
- Atomic metadata capture for every standalone runtime client that might invoke
  concurrent reads/mutations. The protocol avoids that competition by exclusive
  runtime ownership, but is not a replacement for auditing the public runtime.
- Physical time versus the FP32 integration timestep over very long trajectories,
  as well as restart equivalence across numerical implementation versions.

These remain explicit qualification and implementation tasks. Do not interpret
this continuation as a completed Wave B numerical certification.
