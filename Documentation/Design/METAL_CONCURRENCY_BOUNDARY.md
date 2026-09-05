# Legacy Metal SDK concurrency boundary

NumiVivo remains in Swift 6 language mode. The tested macOS 15 SDK imports Metal
Objective-C device, command-queue and resource protocols without complete Sendable
annotations. Metal imports use `@preconcurrency` as an SDK migration boundary;
NumiVivo's own actors, tasks and value types retain strict concurrency checking.
This is not a declaration that arbitrary concurrent buffer mutation is safe.

Device allocation and command-queue submission are shareable services. Pipeline
states are treated as immutable after construction. Command buffers and encoders
remain local to one encoding operation. Mutable arena state belongs to its runtime;
CPU writes, GPU submissions and readbacks must retain runtime reservations until
command completion. No live buffer pointer is returned as a public state snapshot.
Initial resources are populated before their owning runtime is returned.

ProgramPack interventions now run synchronously on the transactional actor: a
checkpoint/read/modify/restore cannot interleave with a step while suspended. A
reversible shutdown changes the same lifecycle read by the transaction engine,
not an independent facade flag. Configuration-bound restoration preserves spatial
geometry and all temporal, velocity, volume and random-state fields.

The import annotation does not prove correct GPU synchronization. Actual Metal
execution and concurrent-operation rejection are separate qualification tests.
Legacy population and partition paths likewise need per-entrypoint runtime tests;
a passing compile is not evidence for every possible concurrent use.

Primary references:
- https://www.swift.org/migration/documentation/swift-6-concurrency-migration-guide/commonproblems/#Preconcurrency-Import
- https://developer.apple.com/documentation/metal/mtlcommandqueue
- https://developer.apple.com/documentation/metal/setting-up-a-command-structure

When the deployment SDK annotates the needed protocols, remove redundant import
annotations as diagnosed by the compiler, without weakening the runtime ownership
and GPU completion rules.
