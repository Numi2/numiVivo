# NumiVivo execution and artifact audit

Audit date: 4 September 2026. Baseline: `a392672481c9f60b112b5b6f291bf684c0690623`.

## Decision

Do not treat the repository's earlier architecture descriptions as evidence that its execution paths are validated. The audit found defects in executable-pack validation, numerical units, actor reentrancy, persisted identity, and artifact filesystem access. The highest-impact inspected paths have been repaired or made to reject unsupported behavior explicitly. Portable checks pass for the components described below. Full Apple compilation and GPU behavior remain unverified.

The next project milestone should be a correctness-qualified execution release, not another expansion of the biological feature list. Freeze the compiler/runtime contracts while closing the remaining source issues and running the Apple gates. No claim of clinical safety, biological fidelity, benchmark leadership, or complete whole-repository defect coverage is made.

## Scope and method

This was a targeted, in-depth audit of the NumiVivo repository's execution-critical path: native ProgramPack serialization and inspection; source units and lowering; Swift pack construction; molecular/physiology transaction interfaces; physiological dosing, checkpoints and prepared tables; molecular–physiology bridge compilation and coordination; immutable object storage and named references. The underlying reaction kernels and hybrid runtime were reviewed as context, but were not executed or statistically validated.

Repository contents were read through the GitHub connector. A partial local source workspace was used for portable compilation and regression checks because a complete network checkout was not available. Native dependencies used by those checks were copied from the repository and their Git blob identities verified. The check harness uses the actual pack inspector, serializer, unit registry, JSON parser, SHA-256 implementation, filesystem helper and coupling-iteration helper—not substitutes for those implementations.

This is not a completed audit of every calibration optimizer, tissue/population kernel, standards exporter, experiment adapter, CLI branch, external NumiLab repository, or biological mechanism. Those areas must not inherit a “passed audit” label from these results.

## Findings and implemented repairs

### A01 — Executable references survived valid hashes

The original native inspector checked section hashes and directory ranges, but did not establish that the indices and bytecode inside those sections were safe for the GPU. A file can have internally consistent hashes and still contain an out-of-range reaction, action, species, parameter, or expression reference.

`PackValidation.cpp` now checks all executable v1 tables, record strides and minimum alignment, bounded counts, string boundaries and UTF-8, numeric ranges, parameter references, reaction arity, action ownership, monitor responses, and cohort coverage. Sparse incidence is reconstructed from reaction stoichiometry and compared against the actual CSR rows. VM programs are checked for valid opcodes, stack depth, termination and program boundaries. Mutable temporal slots and temporal program references have unique control ownership. Repeated table references have explicit work budgets to prevent bounded-size files from causing unbounded validation work.

Serialization invokes the same inspector before returning a pack. The v1 content hash still excludes the header; the reader therefore binds source/fidelity header fields to the hashed compilation manifest and feature flags to the hashed runtime contract. This preserves the existing wire identity without pretending an integrity hash authenticates an author.

### A02 — Swift direct construction bypassed inspection

`VivoProgramPack(data:)` previously performed weaker checks than the native inspector. It now invokes native hash and semantic validation before making a pack available to runtime callers. Data slices are normalized before wire offsets become collection indices. Fingerprint decoding preserves the 32-byte invariant with a bounded container read instead of allocating an arbitrary array before checking its length.

### A03 — Source values and state behavior could change during lowering

Source narrowing could saturate excessive doubles, turn tiny nonzero values into zero, or round discrete counts before a pack-level check could recover the original value. `SourceSemantics.cpp` now checks finite FP32 representability before lowering and verifies source counts, counters and count-valued inputs are nonnegative integers within the ProgramPack FP32 count range. The separate UInt32 hybrid runtime remains the correct path for larger integer populations.

Several declared state types were not actually lowered into their required behavior. The source preflight now rejects unimplemented leaky-integrator, latch, timer, finite-state and permanent-memory declarations rather than emitting ordinary scalars and losing their meaning. Scalar/counter declarations cannot carry ignored half-life or finite-state labels. This is containment of an implementation gap, not implementation of those missing behaviors.

### A04 — Kinetic numeric units were not normalized

A rate declared per minute could reach a seconds-based runtime unchanged. Mixed concentration scales could share one reaction extent even though their numeric values used different units. The mandatory source preflight checks dimensions, numeric bases, reference categories, rate/action output units and conversion-factor representability. Built-in rate parameters must explicitly use the matching species numeric basis and seconds. Noncanonical values are rejected with a conversion diagnostic; this pass does not silently rescale shared parameter tables.

Built-in kinetic units are also checked when importing cached packs by reconstructing the labeled reaction model. General custom bytecode in v1 has lost literal unit annotations, so its imported dimensional semantics cannot be completely reconstructed. New authored custom expressions are checked before lowering. A typed/versioned IR is still required to close the arbitrary-import gap.

### A05 — Physiology preparation was reentrant

The physiology actor did not reserve shared buffers until preparation completed. A second actor invocation could enter while the first awaited GPU completion and overwrite command/input state. Preparation, readback and checkpoint capture now reserve the actor before their first suspension. GPU command failures stop the actor. Command constants are copied into each encoder. No checkpoint can combine an earlier state readback with a later dose cursor.

### A06 — Restore and dose boundaries were insufficiently constrained

Physiology restore now validates table identity/shape, analyte bounds and dose-cursor consistency, uploads a replacement arena, and changes ownership only after successful upload. A near-term dose boundary is never crossed by rounding a shortened interval upward to the configured minimum step. The runtime reports the incompatible minimum instead. Infusion sources are applied as symmetric pre/post half-kicks rather than entirely before clearance. Numerical convergence of the complete scheme still needs verification.

`VivoPhysiologyModelValidator` checks decoded prepared CSR, clearance, compartment, analyte, initial-state and dose tables before allocation. These structural/numeric checks do not prove physiological mass balance or biological calibration.

### A07 — Bridge rate bounds used state dimensions

Rate transfer offsets/minimum/maximum were converted as concentrations rather than concentration per time. The compiler now applies the correct dimension and scale for each update mode. It uses enumerated species indices instead of a nonexistent `SpeciesMetadata.index` field. It rejects nonfinite/unrepresentable transfers, invalid ownership, duplicate destinations and conflicting replacement writers. Mapping cardinality is checked before allocating default index arrays or Cartesian products; the total expansion has a bounded budget.

### A08 — Convergence depended on units

The old residual denominator contained an unqualified numeric `1`, making the stopping decision depend on whether a channel used, for example, a larger or smaller concentration unit. The new helper uses a component-relative residual invariant to nonzero multiplicative unit changes. A bounded Aitken relaxation update may accelerate fixed-point iteration; it never extrapolates beyond the proposed state. Twelve portable checks cover the helper. They do not establish convergence for arbitrary nonlinear or stochastic biological systems.

### A09 — Joint commit was described more strongly than implemented

Two independently throwing actor commits are not a crash-atomic distributed transaction. The coordinator now exclusively owns participants, reserves the entire operation, blocks observations during preparation/release, validates both candidates before release, and does not observe cancellation between release calls. An unexpected release failure permanently quarantines the coordinator rather than pretending a committed participant was rolled back.

This provides logical atomicity to the coordinator's API under its ownership contract, not process-crash atomicity. A durable commit journal or single fused publication authority would be needed for the latter. The coordinator compares actual horizons and Double clocks, rather than accepting mismatched intervals under a tolerance scaled by one second. It does not automatically reduce a stochastic horizon after a realized failure.

### A10 — Source fingerprints did not identify executable state

The source model fingerprint alone did not bind prepared CSR coefficients, clearance overrides, bridge transforms or numerical configuration. The coordinator now computes an executable identity from the full prepared models and configurations. Joint checkpoints store both participant states, accepted generation and feedback history under that identity. Joint restore constructs replacement participants and swaps ownership only after both restores succeed.

Standalone physiology now has a preferred `VivoPhysiologyResumeCheckpoint` format with the same full-execution binding. The older raw checkpoint API remains for compatibility and does not acquire this stronger guarantee automatically; new durable resume workflows should use the envelope or joint checkpoint.

### A11 — Artifact metadata and filesystem paths were insufficiently bound

A valid object hash could be accompanied by forged caller-supplied descriptor metadata. References now have to match the persisted descriptor and verified object. Immutable metadata cannot be silently rewritten for identical bytes. Persisted dates are returned in their canonical serialized precision.

Storage now operates through a pinned directory descriptor. Descendant directory/file opens reject symlinks. Reads reject special files, use nonblocking opens before checking regular-file type, and enforce allocation/read bounds. Immutable entries use no-clobber publication; named references use atomic name replacement. Exact UTF-8 names are hashed into v2 reference filenames to avoid case-folding and Unicode filesystem aliases. Read-only v1 fallback and deletion tombstones preserve compatibility without resurrecting deleted names.

This is not a sandbox against administrators, hostile mounts, or another process with the same privileges. The caller-selected root is trusted. Hashes and embedded signing keys still require a separate source trust policy. See `ARTIFACT_STORE_HARDENING.md` for the filesystem boundary and portability details.

### A12 — Safety-report language exceeded the analysis

The safety report double-counted some irreversible actions and described distinct named signals as independent. Counting is corrected. The report now distinguishes reference coverage and syntactic shutdown presence from demonstrated biological independence, actual output bounds, reachable termination, or physical safety. Constant-false control conditions are reported. This remains a heuristic risk report, not a formal safety proof.

## Verification actually performed

| Check | Executed result | What it does not establish |
|---|---:|---|
| Native pack/semantic fixture and rehashed mutations | 55 passed | Full compiler pipeline, GPU execution, exhaustive fuzzing |
| AddressSanitizer and UndefinedBehaviorSanitizer on that native run | No findings in that run | Absence of all memory/undefined-behavior defects |
| Rooted filesystem helper, real POSIX I/O | 20 passed | APFS/iOS behavior, crash durability, all process races |
| Coupling iteration helper | 12 passed | Nonlinear/stochastic coupled-system convergence |
| Selected Swift syntax parsing | 7 files passed | Swift–Metal–C++ integration typechecking |
| Targeted C++ syntax/type checks | 4 changed translation units passed | Complete package linking or Apple SDK compatibility |

Total targeted executed assertions: **87**. Toolchains: Clang 17.0.0 and Swift 6.2.1 on `x86_64-unknown-linux-gnu`. The native harness recomputes hashes after most payload corruptions, so malformed references cannot pass merely because the integrity check was the only tested defense. The local source identity manifest and execution log are included next to this report.

Reproduce without building the NumiVivo package, invoking Metal, or using the network:

```sh
bash Tools/Audit/run_portable_checks.sh
```

The default native harness enables ASan and UBSan. `NVIVO_AUDIT_SANITIZERS=0` is an explicit portability escape for toolchains without sanitizer support; a run with that option is not sanitizer-checked. `CXX` and `SWIFTC` select installed toolchains. Scratch executables are removed at exit.

## Compatibility consequences

Recompile source programs under the new checks. A previously accepted artifact may now be rejected because its source meaning was not actually executable, its units were mismatched, its counts would round, its table semantics were invalid, or its resource footprint exceeded supported bounds. Do not turn off validation to preserve a misleading result. Convert values and units explicitly, or express unsupported state behavior using represented reactions/rules. Keep the rejected artifact and diagnostics as provenance.

Use full-execution-bound physiology/joint checkpoints for new workflows. Treat old raw physiology checkpoints as legacy data whose prepared model/configuration must be established independently. Existing v1 artifact-store references remain readable through the checked compatibility path; new writes use v2 names. Conflicting metadata for identical raw bytes now requires a separate provenance envelope.

## Remaining source work—not Apple-only tasks

1. **Typed IR and strict source decoding.** Preserve dimensions, reference categories, state-kind semantics and compiler-contract version in the serialized format. Reject unknown/multiple expression operators and unknown required source fields before interpretation. Audit integer conversion edge cases and locale-dependent numeric parsing. This audit did not complete a strict decoder rewrite.
2. **Stochastic rejection and numerical qualification.** The dedicated hybrid SSA continuation must be statistically qualified. The standalone legacy molecular adaptive retry path can select new draws after rejection; it is not an unbiased adaptive stochastic solver. The joint coordinator disables that retry path, but the legacy standalone path still needs an explicit policy or a correct stochastic-bridge/post-leap method. No inference should rely on selecting only successful stochastic runs.
3. **Global allocation authority.** Per-arena sizing and sequential construction do not provide a process-wide reservation across all live runtimes, replacement arenas, checkpoints and coupled NumiLab components. Add one shared budget authority and bounded streaming checkpoint I/O before claiming maximum scale.
4. **Unimplemented dynamics.** Delay/refractory queues, membrane geometry, count-valued spatial hopping, richer state lowering and live hybrid-authority migration remain unsupported where the runtime rejects them. Rejection is the current honest contract, not feature completion.
5. **Broader modules and trust.** Calibration/result provenance coverage, population/tissue kernels, standards round trips, parser fuzzing, signature trust anchors and production adapters in other NumiLab repositories still need dedicated audits. This report does not certify them.

## Apple-silicon handoff

Perform these gates in order on the same committed revision. Preserve toolchain/device metadata and raw logs with the result artifacts.

**Gate 1: package and native interface.** Run the portable checks first, then `swift build -c debug`. Resolve compiler/linker diagnostics across all targets; selected-source parsing is not a substitute. Confirm that the C++ public module exports the intended C ABI and does not depend on accidental header inclusion. Do not interpret this audit as evidence that the whole package currently builds.

**Gate 2: shader ABI and transactions.** With Metal API validation enabled, instantiate every pipeline selected by the molecular, hybrid and physiology paths. Check copied resource lookup, reflected buffer indices, scalar/record sizes, threadgroup limits and memory hazards. Exercise simultaneous prepare/snapshot/restore requests, cancellation at every suspension, dose boundaries, allocation failure, GPU command failure, and a failure during the second participant release. No mixed accepted state may be returned.

**Gate 3: numerical and resume qualification.** Use analytically checkable decay/exchange systems, closed-system mass balance, spatial flux conservation, stochastic moments and rare-event distributions, exact-SSA work-chunk invariance, and step refinement. Compare uninterrupted and checkpoint-resumed accepted trajectories under a documented reproducibility tolerance. Test both raw compatibility imports and new execution-bound resume formats. Validate long-horizon clocks and dose schedules.

Only after those gates pass should throughput, lane capacity, TensorOps/Metal 4 specialization or additional biological mechanisms become the acceptance target. No speedup or biological validity is established by the portable checks.

## References

The repairs are grounded primarily in the repository implementations and the checked source identities. Relevant platform contracts:

- Swift SE-0306, actor reentrancy: https://github.com/swiftlang/swift-evolution/blob/main/proposals/0306-actors.md
- Apple Metal, command encoder: https://developer.apple.com/documentation/metal/mtlcomputecommandencoder
- Apple open(2), O_NOFOLLOW/O_EXCL/nonblocking behavior: https://developer.apple.com/library/archive/documentation/System/Conceptual/ManPages_iPhoneOS/man2/open.2.html

These references describe platform semantics; they are not evidence that the complete NumiVivo implementation has been verified.
