# Executable hybrid Metal runtime

## Implemented path and scope

`VivoExactSSAModel` + `VivoHybridStochasticPlan` -> `VivoHybridExecutionCompiler` -> `VivoHybridReactionRuntime` -> `NumiVivoHybridExecution.metal`.

This is a real GPU execution path, not a scheduler that invokes full-table solvers independently. It implements mixed exact SSA, tau-leap and RK2 execution over disjoint connected components of a sequence-free, well-mixed kinetic model. Six source laws are supported: zero order, first order, second order with distinct inputs, second order with identical inputs, Hill activation and Hill repression. For identical inputs the discrete propensity is `rate*n*(n-1)/2`; the RK2 continuous limit is `rate*x*x/2`.

The compiler verifies the source fingerprint, complete/non-overlapping reaction ownership, exact agreement of uncoupled species, unique component identities, and ownership of every stoichiometric and propensity-only dependency. A supplied planner descriptor cannot hide a cross-component dependency. The executable fingerprint is recomputed from actual model identity, authority masks, and exact reaction groups. It does not trust the plan's descriptive fingerprint.

`VivoHybridExplicitPlanBuilder` defaults unspecified connected components to exact SSA. Approximate authorities require explicit overrides. A work budget cannot silently change the numerical authority. The existing statistical planner may also supply plans, but its heuristics are not treated as a scientific accuracy guarantee.

## GPU layout

Private, tracked Metal heap resources retain FP32 current/stage/candidate/derivative arrays, UInt32 current/candidate count arrays, immutable reaction and sparse incidence records, authority masks, reaction fluxes, sampled firings, and exact continuation records. Species-major indexing keeps adjacent independent lanes contiguous. Deterministic flux evaluation runs per reaction/lane; derivative gather runs per species/lane. State updates need no floating-point atomics.

Exact cohorts dispatch one owner per component/lane. Its selected reaction list is explicit. Other authorities cannot read or write those dynamic species. Constant species are copied unchanged. RK2 evaluates each reaction twice; tau sampling and incidence gather address only tau-owned records. Pass constants are copied with `setBytes`, so one shared CPU buffer is not overwritten between encoded stages. Threadgroup width derives from each compiled pipeline, and partial final groups are permitted.

One command queue orders all authorities and validations. Host interaction per step consists of a 32-byte status and requested compact publications. Full state is read only for explicit snapshots/checkpoints. Pipeline sets are cached by Metal device registry ID. Allocation planning includes heap alignment, buffer-size limits, shared output capacity, live allocation, and startup staging. Replacement restore accounts for the old arena remaining resident.

## Transactions, cancellation and SSA continuation

`prepareStep` reserves the actor before the first await. An in-flight reservation prevents another prepare, snapshot, restore or mutation from touching shared GPU resources. Successful preparation retains a candidate. `commitPreparedStep` flips both authoritative buffer pairs and updates the Double clock/UInt64 accepted-step index in one non-suspending actor operation. `discardPreparedStep` publishes nothing.

Exact SSA continuation preserves elapsed candidate time, event count and random cursor across GPU work chunks. It never redraws a candidate because one dispatch reached its event budget. Total work exhaustion produces `exactWorkBudgetExceeded`, not an automatically accepted shorter horizon. Completed lanes do not run again during continuation. Count changes for a selected event are validated before any event component is written.

Random values come from Philox-style counter rounds with separate SSA/tau domains and accepted-step/lane/reaction-or-cohort namespaces. Uniform values are generated on an FP32 midpoint grid strictly inside (0,1). Event-time nonprogress, invalid propensity, UInt32 underflow/overflow, nonfinite continuous state and bounded Poisson-sampler exhaustion reject the joint candidate. Tau sampling never substitutes a Gaussian result after sampler exhaustion and never clips negative state into validity.

## Checkpoints and reproducibility boundary

Accepted-boundary checkpoints store numerical ABI, model/executable-plan identities, seed, accepted step, Double time, lane/species shape, FP32 state and UInt32 counts. Restore checks identity, seed, version, shape, finite/nonnegative continuous data and zeroed inactive representations. It uploads into a replacement arena and switches ownership only after successful completion. A allocation/upload failure does not partially overwrite the accepted arena.

The accepted step and seed determine random stream reconstruction. In-flight SSA candidates are not exported. Changing exact work-chunk size is intended not to change the accepted trajectory, provided all work completes; that property still requires numerical verification. Cross-device or cross-compiler bitwise trajectory equality is not promised. `exactSSA` identifies the direct-event algorithm; rates, logarithms and local waiting times use FP32 and are not exact real arithmetic. Tau-leaping and RK2 remain approximations. The plan's maximum-step value is advisory, not an error estimator.

## Explicit remaining boundaries

This path does not execute general VivoProgram expression bytecode, temporal rules, delayed reactions, monitors, spatial transport or live authority migration. A spatial plan is rejected rather than treated as well mixed. The hybrid implementation loads its own shader library and owns its own checkpoint contract. The separate ProgramPack backend and molecular v2 restore integration are described in [PROGRAM_PACK_METAL_BACKEND.md](PROGRAM_PACK_METAL_BACKEND.md); those changes do not qualify the hybrid path.

The original implementation pass did not execute a package build, shader compilation, simulation, statistical test, checkpoint round trip or performance benchmark. Source/API/ABI review is not execution evidence.

## Native conformance suite

`Tests/NumiVivoIntegrationTests/HybridRuntimeTests.swift` exercises this dedicated backend on an actual production-selected Metal device. It fails when suitable hardware is absent. The suite covers:

- joint exact-SSA, tau-leap and RK2 execution over disjoint components, count conservation, RK2 agreement with the analytical Heun recurrence, inactive representation invariants, partial final SIMD groups, and publications preserving all UInt32 count bits;
- bitwise equality of accepted snapshots and checkpoints across exact event dispatch chunk sizes of 1 and 128;
- mixed-authority checkpoint continuation, refusal of incompatible identities/seed/ABI or malformed/noncanonical payloads, unchanged accepted state after refused restore, and prepared-candidate reservation/discard;
- a guaranteed UInt32-overflow rejection and a guaranteed exact-work-budget exhaustion, verifying that neither commits the independent continuous component or returns accepted publications;
- Poisson birth moments for exact SSA and low-mean tau inversion, a larger-mean tau PTRS case, and the Binomial survivor distribution of exact first-order decay.

Statistical checks use 8,192 independent lanes and fixed seeds. Mean and unbiased sample-variance tolerances are declared as six analytical standard errors before execution; the small-mean Poisson cases also check the zero-event probability. These checks do not qualify arbitrary tau-leap error, every supported propensity law, rare-event tails, connected-network hybrid methods or cross-device identity.

Run with `swift test -c release --jobs 3 -Xswiftc -enable-testing --filter HybridRuntimeTests` on an available Apple-silicon GPU. Setting `NUMIVIVO_TEST_ARTIFACTS` writes the accepted-state, continuation and moment reports outside the repository. The presence of the suite is not a passing result; retain its command output, exact revision and physical-device evidence when reporting qualification. No performance claim follows from this suite.

The [7 September 2026 reliability audit](../Audit/2026-09-07_COMPLETION_RELIABILITY.md) records all seven tests passing on the physical M4 Pro, including the distribution observations and their predeclared tolerances. Cross-device equality and the remaining method boundaries above are unchanged.

## API references consulted

- Apple, `MTLComputeCommandEncoder.setBytes(_:length:index:)`: https://developer.apple.com/documentation/metal/mtlcomputecommandencoder/setbytes(_:length:index:)
- Apple, `MTLComputeCommandEncoder`: https://developer.apple.com/documentation/metal/mtlcomputecommandencoder
- Swift language guide, concurrency and actor isolation: https://docs.swift.org/swift-book/documentation/the-swift-programming-language/concurrency/
