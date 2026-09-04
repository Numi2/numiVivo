# Execution-focused development — 4 September 2026

## Delivered source

- Executable disjoint exact-SSA, tau-leap and RK2 GPU runtime over `VivoExactSSAModel`, with true UInt32 discrete state and FP32 continuous state.
- Dependency-checked authority plans, default-exact plan builder, sparse GPU incidence, exact-event continuation across bounded dispatches, compact publications, joint candidate commit/discard, full accepted-boundary checkpoint and replacement-arena restore.
- `plan-hybrid`, `run-hybrid`, `hybrid-help` CLI commands, checkpoint resume, result metadata and a mixed-authority four-lane synthetic example at `Examples/hybrid-reaction-runtime`.
- Public molecular v2 checkpoint/restore and configuration-bound spatial resume APIs on both `VivoTransactionalMolecularRuntime` and `VivoRuntime`.
- In-flight actor reservations, stopped/cancelled candidate protection, copied per-pass constants, separate split-step random namespaces and a Double accumulated molecular clock.
- A new ProgramPack Metal backend aligned with current compiler records, named entry points and actor buffer bindings; native integrity/semantic preflight; ordered per-lane rules; bounded stochastic state; finite-volume concentration transport and ownership-aware publications.
- One authoritative molecular runtime behind the public facade; separate shader-module loading and per-device library reuse; explicit copied shader resources.
- Digest-chain generic encoding helper moved to file scope with unchanged hash field names and checked index growth.

## Remaining development boundaries

The dedicated hybrid path is a well-mixed kinetic runtime, not a general VivoProgram rules/delays/spatial interpreter. The ProgramPack path rejects delayed reactions, rule refractory periods, temporal rate/gate operators, count-valued spatial transport and membrane permeability without an interface model. General live authority migration, delay/refractory state integration, production adapters into the other NumiLab repositories, and a fused shared-GPU multi-participant commit protocol remain separate work.

This development pass did not run a Swift/C++ build, Metal compile, simulation, test suite, benchmark, or CI job. Static source/ABI inspection and successful Git commits do not establish compiler acceptance, numerical agreement or measured performance. Start from the exact committed source when performing those checks on Apple silicon.

See `Design/EXECUTABLE_HYBRID_RUNTIME.md`, `Design/MOLECULAR_RESUME_INTEGRATION.md`, and the superseding backend integration note `Design/PROGRAM_PACK_METAL_BACKEND.md` for precise scope and interfaces.
