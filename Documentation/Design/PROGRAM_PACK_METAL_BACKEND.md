# ProgramPack Metal execution backend

## Source integration status

The previous primary molecular catalog named kernels that did not match the bundled argument-buffer implementation. The active catalog now loads `NumiVivoProgramPackRuntime.metal`, whose twelve entry points match `VivoTransactionalMolecularRuntime` buffer bindings and the current ProgramPack v1 record layouts. `VivoRuntime` delegates to that engine instead of owning a divergent second arena and encoder implementation.

This supersedes the earlier backend-mismatch notes in EXECUTABLE_HYBRID_RUNTIME.md and MOLECULAR_RESUME_INTEGRATION.md. It is source integration, not evidence that an Apple compiler or GPU has accepted the implementation. No build, shader compilation, simulation, numerical test or benchmark was run in this development pass.

## Supported execution

F0 copies staged input state, executes ordered rules, evaluates monitors, validates and publishes. F1 uses two-stage Heun integration over sparse stoichiometric incidence. Its rates include zero order, mass action, Hill activation/repression, Michaelis-Menten, reversible mass action, passive/saturable exchange, degradation and pure custom expression bytecode. Parameters use a private parameter-major FP32 matrix with one column per environment. Time-dependent pure rate expressions use start/end stage time.

F2 uses the declared nonnegative rate-law values as intensities and requires count-valued reaction state. Reversible laws sample forward and reverse intensities separately. The bounded Poisson implementation sums independent low-mean Poisson samples, rejects exhausted draws, and requires a smaller step when either directional mean exceeds 1024. Negative results are rejected, not clipped. Count state is bounded to exactly representable FP32 integers through 2^24. The dedicated hybrid runtime is the distinct UInt32 implementation and supports explicit combinatorial same-reactant propensities; the two APIs must not be confused.

F3/F4 use reaction/finite-volume-transport/reaction splitting for concentration state. Each internal face uses paired diffusion and upwind-advection fluxes, scaled by the minimum adjacent volume fraction. Neighbor lookups remain inside the same environment. No-flux, periodic and absorbing boundary modes are represented. A local outgoing-flux bound requests subdivision before positivity can be lost through excessive transport. Externally owned species are excluded using ProgramPack metadata, even when a global velocity field is present. F4 still needs an external coupled participant for tissue mechanics; the mode alone is not a complete tissue solver.

Rules are evaluated serially within each lane in compiler order; independent lanes run in parallel. Set/add, expression/degradation actions, event requests and shutdown actions are implemented. Temporal rule/monitor operators retain Float2 state in the transaction's temporal buffer. Constraint monitors and termination monitors have distinct trigger polarity. A requested clamp has no general expression-to-state inverse and therefore rejects rather than claiming a clamp occurred.

## Explicitly unsupported contracts

Host preflight rejects delayed reactions, nonzero rule refractory periods, temporal operators inside concurrently evaluated reaction rates/gates, custom bytecode stochastic propensities, non-count F2 reaction state, and count-valued spatial transport. Nonzero membrane permeability requires an explicit membrane-interface model and is rejected by transport configuration. These cases need additional authoritative state or numerical semantics; no omitted queue, timer, count-to-concentration conversion or membrane geometry is fabricated.

## ABI and resources

ProgramPack has a 128-byte header and 72-byte descriptors. The shader reads the validated directory and follows private-buffer section offsets. Species, reactions, expressions, actions, rules, monitors, incidence and stoichiometry use the same strides as the compiler and native pack inspector. The command/status/event ABIs are 128/64/32 bytes respectively, with shader static assertions and existing Swift layout checks.

Every pass receives a copied command through `setBytes`. The transport pass additionally binds ProgramPack at index 8 for ownership metadata. There are no CPU writes into a single shared uniform buffer between pending half-steps. Rule/monitor dispatch is one thread per lane rather than one simultaneous writer per rule and lane.

The pipeline catalog loads molecular and physiology modules separately and caches compiled libraries by device and module identity. Partial final threadgroups are permitted. The package copies individual shader sources unchanged to the resource bundle root. This avoids a platform-dependent source-processing rule replacing the exact runtime source or concatenating unrelated legacy ABIs. Historical shader sources remain in the repository but are not selected by the active molecular catalog.

## Checkpoint access through the public runtime

```swift
let pack = try VivoProgramPack(data: Data(contentsOf: packURL))
let configuration = VivoRuntimeConfiguration(fidelity: .deterministic)
let runtime = try await VivoRuntime.make(pack: pack, configuration: configuration)
let result = try await runtime.step()
let saved = try await runtime.resumeCheckpoint()
let replacement = try await VivoRuntime.make(pack: pack, configuration: configuration)
try await replacement.restore(saved)
```

This is an API example, not a recorded run. Configuration-bound resume stores the spatial layout, numerical policy, seed and complete accepted state represented by this backend. Restore is a replacement-arena operation and cannot resume a stopped actor by bypassing its lifecycle. Full snapshots/checkpoints are explicit operations, not per-step state downloads.

## Required downstream verification

On an Apple-silicon environment, verify shader entry points and buffer reflection, rate-law units and signs, Heun convergence, count conservation and stochastic moments, paired finite-volume flux balance, temporal trigger/retry semantics, lifecycle reentrancy, resource-copy lookup, and checkpoint continuation. No speedup, maximum biological scale, clinical validity or cross-device bitwise identity is established by these source changes.

Reference: Apple Swift package copy-resource rule, https://developer.apple.com/documentation/packagedescription/resource/copy(_:)
