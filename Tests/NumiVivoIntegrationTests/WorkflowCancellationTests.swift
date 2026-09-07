import Foundation
import Testing
@testable import NumiVivoKit

/// Numerical work deliberately remains suspended after cancellation so tests
/// can prove that callers detach and replacements do not overlap that work.
private actor CancellationWorkProbe {
    private var releases: [Int: CheckedContinuation<Void, Never>] = [:]
    private var releasingAll = false
    var calls = 0
    var active = 0
    var peak = 0
    var observedCancellation = 0

    func execute() async -> [String: Data] {
        calls += 1; active += 1; peak = max(peak, active)
        let call = calls
        if !releasingAll { await withCheckedContinuation { releases[call] = $0 } }
        if Task.isCancelled { observedCancellation += 1 }
        active -= 1
        return ["value": Data("7".utf8)]
    }
    func release(_ call: Int) { releases.removeValue(forKey: call)?.resume() }
    func releaseAll() {
        releasingAll = true
        let pending = Array(releases.values)
        releases.removeAll()
        for continuation in pending { continuation.resume() }
    }
}

private actor CancellationCompletion<Value: Sendable> {
    private(set) var result: Result<Value, Error>?
    var finished: Bool { result != nil }
    func finish(_ result: Result<Value, Error>) { self.result = result }
}

private struct CancellationObservedRun<Value: Sendable>: Sendable {
    let task: Task<Void, Never>
    let completion: CancellationCompletion<Value>
    func cancel() { task.cancel() }
}

/// Observe completion without awaiting a possibly stuck task.value. A timeout
/// can therefore unwind the test and release every suspended numerical call.
private actor CancellationRuns {
    private var cancellations: [@Sendable () -> Void] = []
    private var completions: [@Sendable () async -> Bool] = []

    func start<Value: Sendable>(_ operation: @escaping @Sendable () async throws -> Value) -> CancellationObservedRun<Value> {
        let completion = CancellationCompletion<Value>()
        let task = Task {
            do { await completion.finish(.success(try await operation())) }
            catch { await completion.finish(.failure(error)) }
        }
        cancellations.append { task.cancel() }
        completions.append { await completion.finished }
        return .init(task: task, completion: completion)
    }
    func cancelAll() { for cancel in cancellations { cancel() } }
    func allFinished() async -> Bool {
        for finished in completions { if !(await finished()) { return false } }
        return true
    }
}

private struct CancellationContext: Sendable {
    let store: VivoArtifactStore
    let workflow: VivoChemistryWorkflow
    let probe: CancellationWorkProbe
    let operation: VivoChemistryOperation
    let request: VivoChemistryTask
    let executor: VivoWorkflowExecutor
    let recipe: VivoWorkflowRecipe
    let runs: CancellationRuns
}

@Suite(.serialized) struct WorkflowCancellationTests: Sendable {
    private func directory() throws -> URL {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-cancel-\(UUID().uuidString)")
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        return root
    }
    private func operation(_ probe: CancellationWorkProbe) throws -> VivoChemistryOperation {
        .init(identifier: "test.cancellation", version: "1",
              implementationFingerprint: try VivoCanonicalJSON.fingerprint(Data("cancel-tests".utf8)),
              outputs: [.init(name: "value", kind: "test.number")],
              execute: { _, _, _ in await probe.execute() },
              validateOutputs: { _, _, output, _ in
                  guard output == ["value": Data("7".utf8)] else { throw VivoChemistryError.invalid("test output") }
              })
    }
    private func task(_ operation: VivoChemistryOperation) -> VivoChemistryTask {
        .init(operation: operation.identifier, version: operation.version,
              implementationFingerprint: operation.implementationFingerprint, inputs: [],
              configuration: .object([:]), outputs: operation.outputs)
    }
    private func until(_ condition: @Sendable () async -> Bool) async throws {
        let deadline = ContinuousClock.now.advanced(by: .seconds(10))
        while !(await condition()) {
            guard ContinuousClock.now < deadline else { throw VivoChemistryError.invalid("cancellation test handshake timed out") }
            try await Task.sleep(for: .milliseconds(1))
        }
    }
    private func completed<Value: Sendable>(_ run: CancellationObservedRun<Value>) async throws -> Result<Value, Error> {
        try await until { await run.completion.finished }
        return await run.completion.result!
    }
    private func value<Value: Sendable>(_ run: CancellationObservedRun<Value>) async throws -> Value {
        try await completed(run).get()
    }
    private func cancelled<Value: Sendable>(_ run: CancellationObservedRun<Value>) async throws {
        switch try await completed(run) {
        case .success: Issue.record("cancelled caller returned success")
        case .failure(is CancellationError): break
        case .failure(let error): Issue.record("unexpected cancellation error: \(error)")
        }
    }
    private func cleanup(_ context: CancellationContext) async {
        await context.runs.cancelAll()
        await context.probe.releaseAll()
        do {
            try await until {
                let callersDone = await context.runs.allFinished()
                let coreDone = await context.workflow.activity().isEmpty
                let executorDone = await context.executor.activity().isEmpty
                return callersDone && coreDone && executorDone
            }
        } catch { Issue.record("cancellation test cleanup did not drain: \(error)") }
    }
    private func withContext(_ body: @Sendable (CancellationContext) async throws -> Void) async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root), probe = CancellationWorkProbe()
        let operation = try operation(probe)
        let registry = try VivoWorkflowRegistry(implementationFingerprint: operation.implementationFingerprint,
            definitions: [.init(operation: operation, inputKinds: [:], summary: "gated test", validationScope: "constant test output", validateConfiguration: { _ in })])
        let recipe = VivoWorkflowRecipe(identifier: "cancel-recipe", artifacts: [], nodes: [
            .init(identifier: "source", operation: operation.identifier, inputs: [:])
        ], outputs: [.init(name: "answer", node: "source", port: "value")])
        let context = CancellationContext(store: store, workflow: VivoChemistryWorkflow(store: store), probe: probe,
            operation: operation, request: task(operation), executor: VivoWorkflowExecutor(store: store, registry: registry),
            recipe: recipe, runs: CancellationRuns())
        do { try await body(context) }
        catch { await cleanup(context); throw error }
        await cleanup(context)
    }

    @Test func cancellingOneOfTwoCallersPreservesTheSharedExecution() async throws {
        try await withContext { context in
            let workflow = context.workflow, probe = context.probe, operation = context.operation, request = context.request
            let first = await context.runs.start { try await workflow.run(request, using: operation) }
            let second = await context.runs.start { try await workflow.run(request, using: operation) }
            try await until { await workflow.activity().first?.waitingCallers == 2 }
            try await until { await probe.calls == 1 }
            first.cancel()
            try await cancelled(first) // Must finish before numerical work is released.
            #expect(await workflow.activity().first?.waitingCallers == 1)
            #expect(await workflow.activity().first?.cancellationRequested == false)
            await probe.release(1)
            let result = try await value(second)
            #expect(!result.reused)
            #expect(await probe.calls == 1)
            #expect(await probe.observedCancellation == 0)
            #expect(await workflow.activity().isEmpty)
            let cached = await context.runs.start { try await workflow.run(request, using: operation) }
            #expect(try await value(cached).reused)
        }
    }

    @Test func abandonedGenerationDrainsBeforeAReplacementAndCannotPublishItsOutput() async throws {
        try await withContext { context in
            let store = context.store, workflow = context.workflow, probe = context.probe
            let operation = context.operation, request = context.request
            let abandoned = await context.runs.start { try await workflow.run(request, using: operation) }
            try await until { await probe.calls == 1 }
            abandoned.cancel()
            try await cancelled(abandoned)
            #expect(await workflow.activity().first?.cancellationRequested == true)
            #expect(await workflow.activity().first?.waitingCallers == 0)
            let replacement = await context.runs.start { try await workflow.run(request, using: operation) }
            try await until { await workflow.activity().first?.waitingCallers == 1 }
            #expect(await probe.calls == 1)
            await probe.release(1)
            try await until { await probe.calls == 2 }
            #expect(await probe.observedCancellation == 1)
            #expect(await probe.peak == 1)
            #expect(try await store.list(kind: "chemistry-task-receipt").isEmpty)
            #expect(try await store.list(kind: "chemistry-output").isEmpty)
            await probe.release(2)
            let result = try await value(replacement)
            #expect(!result.reused)
            let cached = await context.runs.start { try await workflow.run(request, using: operation) }
            #expect(try await value(cached).reused)
            #expect(await probe.calls == 2)
            #expect(await workflow.activity().isEmpty)
        }
    }

    @Test func callerQueuedBehindDrainingWorkCanAlsoCancel() async throws {
        try await withContext { context in
            let store = context.store, workflow = context.workflow, probe = context.probe
            let operation = context.operation, request = context.request
            let first = await context.runs.start { try await workflow.run(request, using: operation) }
            try await until { await probe.calls == 1 }
            first.cancel(); try await cancelled(first)
            let queued = await context.runs.start { try await workflow.run(request, using: operation) }
            try await until { await workflow.activity().first?.waitingCallers == 1 }
            queued.cancel(); try await cancelled(queued)
            await probe.release(1)
            try await until { await workflow.activity().isEmpty }
            #expect(await probe.calls == 1)
            #expect(try await store.list().isEmpty)
        }
    }

    @Test func alreadyCancelledCallerDoesNotStartWorkOrCreateCacheEntries() async throws {
        try await withContext { context in
            let workflow = context.workflow, operation = context.operation, request = context.request
            let run = await context.runs.start {
                withUnsafeCurrentTask { $0?.cancel() }
                return try await workflow.run(request, using: operation)
            }
            try await cancelled(run)
            #expect(await context.probe.calls == 0)
            #expect(await workflow.activity().isEmpty)
            #expect(try await context.store.list().isEmpty)
        }
    }

    @Test func cancellingRecipeDuringNumericalWorkDoesNotPublishAReport() async throws {
        try await withContext { context in
            let executor = context.executor, recipe = context.recipe, probe = context.probe
            let run = await context.runs.start { try await executor.run(recipe) }
            try await until { await probe.calls == 1 }
            run.cancel()
            try await cancelled(run)
            await probe.release(1)
            // Numerical return precedes workflow retirement; wait for the latter
            // so a late, erroneous cache publication cannot escape this assertion.
            try await until { await executor.activity().isEmpty }
            #expect(try await context.store.list(kind: "vivo.workflow-run").isEmpty)
            #expect(try await context.store.list(kind: "chemistry-task-receipt").isEmpty)
        }
    }

    @Test func restartingAnExecutorWaitsForOldWorkAndRejectsItsAbandonedOutput() async throws {
        try await withContext { context in
            let executor = context.executor, recipe = context.recipe, probe = context.probe, store = context.store
            let abandoned = await context.runs.start { try await executor.run(recipe) }
            try await until { await probe.calls == 1 }
            abandoned.cancel()
            try await cancelled(abandoned)
            #expect(await executor.activity().first?.cancellationRequested == true)

            let replacement = await context.runs.start { try await executor.run(recipe) }
            try await until { await executor.activity().first?.waitingCallers == 1 }
            #expect(await executor.activity().first?.cancellationRequested == true)
            #expect(await probe.calls == 1)
            #expect(await probe.active == 1)
            #expect(try await store.list(kind: "vivo.workflow-run").isEmpty)

            await probe.release(1)
            try await until { await probe.calls == 2 }
            #expect(await probe.observedCancellation == 1)
            #expect(await probe.peak == 1)
            #expect(try await store.list(kind: "chemistry-output").isEmpty)
            #expect(try await store.list(kind: "chemistry-task-receipt").isEmpty)

            await probe.release(2)
            let result = try await value(replacement)
            #expect(result.report.allTasksSucceeded)
            guard let outcome = result.report.nodes.first?.outcome, case .succeeded(_, _, _, let reused) = outcome else {
                Issue.record("replacement recipe did not return its successful node"); return
            }
            #expect(!reused)
            #expect(await executor.activity().isEmpty)
            let cached = await context.runs.start { try await executor.run(recipe) }
            let cachedResult = try await value(cached)
            guard let cachedOutcome = cachedResult.report.nodes.first?.outcome, case .succeeded(_, _, _, let reusedCache) = cachedOutcome else {
                Issue.record("cached recipe did not return its successful node"); return
            }
            #expect(reusedCache)
            #expect(await probe.calls == 2)
            #expect(await probe.peak == 1)
        }
    }
}
