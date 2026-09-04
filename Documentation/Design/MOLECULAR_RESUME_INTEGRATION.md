# Molecular transaction and resume integration

## Public source APIs

`VivoTransactionalMolecularRuntime` now exposes `checkpoint()`, `resumeCheckpoint()`, `restore(VivoMolecularCheckpoint)`, and `restore(VivoMolecularResumeCheckpoint)`. This supersedes the earlier implementation note in EXECUTABLE_HYBRID_RUNTIME.md that these APIs remained unexposed.

A raw v2 checkpoint is supported for nonspatial configurations. The new configuration-bound resume envelope stores the full runtime configuration alongside v2 state, including grid dimensions, spacing, boundary mode, environment count, seed and numerical settings. Restoring the envelope requires exact configuration equality. Neither route changes the immutable model or seed.

Capture is permitted only at a ready, accepted boundary. A stopped or failed actor cannot be resumed through a checkpoint mutation; the caller must construct an appropriately configured new runtime. Models with delayed reactions are rejected by these checkpoint APIs because the current ProgramPack arena does not own a delayed-event queue. The checkpoint must not imply preservation of unrepresented state.

Restore validates model/content identities, shapes, seed, fidelity, finite state, species bounds, parameter bounds, count integrality and nonnegative transport coefficients. Zero-parameter models export an empty semantic parameter payload and recreate the arena's private dummy scalar only during import. Restoration allocates and uploads a complete replacement arena; authoritative buffers and the clock are changed only after successful upload. Temporary memory therefore includes both arenas and transfer staging; an allocation failure leaves the previous accepted arena owned by the actor.

## Transaction changes

The actor reserves itself before its first GPU suspension. A second prepare, checkpoint, snapshot or state mutation cannot reuse resources while an operation is in flight. A permanent stop received during GPU execution prevents the completed candidate from becoming publishable. Cancellation is observed after GPU completion, before pending publication. Command failures stop the actor.

Command constants are copied into each encoder through `setBytes`; no encoded substep refers to a shared CPU uniform buffer that is subsequently overwritten. Spatial half-reactions receive different random substep namespaces. Dispatch counts are bounded by the UInt32 shader grid. Internal time accumulates in Double, and the legacy FP32 absolute-time interface fails explicitly when a requested step is below its representable resolution instead of claiming an advancing timestamp.

Coupling inputs are restricted to externally owned/input species and reject duplicate destinations. Transport and velocity mutations validate finite data. Checkpoint import checks species and parameter bounds before allocation. The status field used for shutdown diagnostics is the actual `firstMonitor` ABI field.

## Verification and remaining kernel boundary

This is source integration, not execution evidence. The original ProgramPack shader catalog still has naming and argument-layout inconsistencies with its legacy Metal resources. These public APIs and transaction changes do not resolve that separate backend mismatch. The executable hybrid path uses its own independently loaded Metal source and ABI and is not affected by the legacy binding convention.

No build, shader compile, runtime invocation, checkpoint round trip, numerical test or performance benchmark was run in this pass. Deferred validation must include Swift concurrency diagnostics, both shader paths, raw and configuration-bound checkpoint round trips, cancellation during preparation, reentrant mutation attempts, and restore allocation failure.
