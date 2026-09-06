import Foundation
import Testing
@testable import NumiVivoKit

private actor WorkflowProbe {
    var calls = 0
    var active = 0
    var peak = 0
    func evaluate(_ cfg: VivoJSONValue, _ data: [String: Data]) async throws -> [String: Data] {
        calls += 1; active += 1; peak = max(peak, active)
        defer { active -= 1 }
        try await Task.sleep(for: .milliseconds(40))
        if cfg == .object(["fail": .boolean(true)]) { throw VivoChemistryError.convergence("deliberate test failure") }
        return ["value": data["value"]!]
    }
}
@Suite(.serialized) struct PlatformWorkflowTests {
    private func identity(_ name: String = "platform-test-implementation") throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(name.utf8))
    }
    private func registry(_ probe: WorkflowProbe, id: VivoFingerprint) throws -> VivoWorkflowRegistry {
        let operation = VivoChemistryOperation(identifier: "test.copy", version: "1", implementationFingerprint: id,
            outputs: [.init(name: "value", kind: "test.value")], execute: { cfg, input, _ in try await probe.evaluate(cfg, input) },
            validateOutputs: { _, input, output, _ in
                guard input["value"] == output["value"] else { throw VivoChemistryError.invalid("copy test integrity") }
            })
        return try .init(implementationFingerprint: id, definitions: [.init(operation: operation, inputKinds: ["value": "test.value"],
            summary: "Scheduler test operation, not a scientific kernel", validationScope: "test byte equality", validateConfiguration: { cfg in
                guard cfg == .object([:]) || cfg == .object(["fail": .boolean(true)]) else { throw VivoChemistryError.invalid("test settings") }
            })])
    }
    private func recipe(_ policy: VivoWorkflowSchedulingPolicy = .init()) -> VivoWorkflowRecipe {
        .init(identifier: "scheduler-test", artifacts: [.init(identifier: "source", source: .json(kind: "test.value", payload: .integer(7)))], nodes: [
            .init(identifier: "a", operation: "test.copy", inputs: ["value": .artifact(identifier: "source")]),
            .init(identifier: "b", operation: "test.copy", inputs: ["value": .artifact(identifier: "source")]),
            .init(identifier: "c", operation: "test.copy", inputs: ["value": .output(node: "a", port: "value")])
        ], outputs: [.init(name: "answer", node: "c", port: "value")], policy: policy)
    }
    private func directory() throws -> URL {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-platform-\(UUID().uuidString)", isDirectory: true)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: true); return url
    }
    private func rejects(_ body: () throws -> Void) {
        do { try body(); Issue.record("expected recipe rejection") } catch {}
    }
    private func record<T: Encodable>(_ value: T, _ name: String) throws {
        guard let path = ProcessInfo.processInfo.environment["NUMIVIVO_TEST_ARTIFACTS"] else { return }
        let directory = URL(fileURLWithPath: path, isDirectory: true)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true)
        try VivoCanonicalJSON.encode(value).write(to: directory.appendingPathComponent(name + ".json"), options: .atomic)
    }
    @Test func preflightRejectsCyclesTypesMissingPortsAndBudgets() throws {
        let registry = try registry(WorkflowProbe(), id: identity()), valid = recipe()
        let plan = try VivoWorkflowPlanner.compile(valid, registry: registry)
        #expect(plan.topologicalOrder == ["a", "b", "c"])
        var changed = valid
        changed.nodes[0].inputs["value"] = .output(node: "c", port: "value")
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.nodes[0].inputs.removeValue(forKey: "value")
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.artifacts[0].source = .json(kind: "wrong.kind", payload: .integer(7))
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.nodes[0].version = "999"
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.nodes[0].configuration = .object(["unrecognized": .boolean(true)])
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.policy.maximumReservedBytes = 128 * 1024 * 1024
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.nodes[0].resources = .init(budget: .init(maximumBytes: Int.max))
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        changed = valid; changed.outputs[0].port = "missing"
        rejects { _ = try VivoWorkflowPlanner.compile(changed, registry: registry) }
        try record(plan, "platform-static-plan")
    }
    @Test func orderingAndImplementationBindingsAreExplicit() throws {
        let registry = try registry(WorkflowProbe(), id: identity()), first = recipe()
        var second = first; second.nodes.reverse(); second.artifacts.reverse(); second.outputs.reverse()
        #expect(try first.fingerprint() == second.fingerprint())
        #expect(try VivoWorkflowPlanner.compile(first, registry: registry) == VivoWorkflowPlanner.compile(second, registry: registry))
        second.nodes[0].configuration = .object(["fail": .boolean(true)])
        #expect(try first.fingerprint() != second.fingerprint())
        let other = try self.registry(WorkflowProbe(), id: identity("different-build"))
        #expect(try VivoWorkflowPlanner.compile(first, registry: registry).implementationFingerprint != VivoWorkflowPlanner.compile(first, registry: other).implementationFingerprint)
    }
    @Test func deduplicationResumeAndReportVerification() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root), probe = WorkflowProbe()
        let registry = try registry(probe, id: identity()), executor = VivoWorkflowExecutor(store: store, registry: registry)
        let first = try await executor.run(recipe())
        #expect(first.report.allTasksSucceeded)
        #expect(await probe.calls == 2) // a/b share a task; c consumes a provenance envelope.
        #expect(first.report.maximumAdmittedTasks == 2)
        let second = try await executor.run(recipe())
        #expect(second.report.nodes.allSatisfy { $0.outcome.reused })
        #expect(await probe.calls == 2)
        let verified = try await executor.verify(first.artifact.fingerprint)
        #expect(verified.report.allTasksSucceeded)
        #expect(verified.report.exports == first.report.exports)
        try record(first.report, "platform-deduplicated-run")
    }
    @Test func independentBranchesRunConcurrentlyButRespectAdmission() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        var request = recipe()
        request.artifacts.append(.init(identifier: "other", source: .json(kind: "test.value", payload: .integer(8))))
        request.nodes[1].inputs["value"] = .artifact(identifier: "other")
        let probe = WorkflowProbe(), registry = try registry(probe, id: identity())
        let result = try await VivoWorkflowExecutor(store: VivoArtifactStore(rootURL: root), registry: registry).run(request)
        #expect(await probe.peak == 2)
        #expect(result.report.maximumAdmittedTasks == 2)
        #expect(result.report.maximumReservedBytes <= request.policy.maximumReservedBytes)
        let otherRoot = try directory(); defer { try? FileManager.default.removeItem(at: otherRoot) }
        request.policy.maximumReservedBytes = try request.nodes[0].reservationBytes()
        let serialProbe = WorkflowProbe()
        let serial = try await VivoWorkflowExecutor(store: VivoArtifactStore(rootURL: otherRoot), registry: self.registry(serialProbe, id: identity())).run(request)
        #expect(serial.report.maximumAdmittedTasks == 1)
        #expect(await serialProbe.peak == 1)
        #expect(serial.report.allTasksSucceeded)
    }
    @Test func failedBranchesRemainVisibleAndDoNotInvalidateIndependentResults() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let probe = WorkflowProbe(), store = try VivoArtifactStore(rootURL: root)
        let executor = try VivoWorkflowExecutor(store: store, registry: registry(probe, id: identity()))
        var request = recipe(); request.nodes[0].configuration = .object(["fail": .boolean(true)])
        let failed = try await executor.run(request)
        #expect(!failed.report.allTasksSucceeded)
        guard case .failed = failed.report.nodes[0].outcome else { Issue.record("failure removed"); return }
        guard case .blocked(let dependencies) = failed.report.nodes[2].outcome else { Issue.record("blocked child executed"); return }
        #expect(dependencies == ["a"])
        #expect(failed.report.nodes[1].outcome.outputs != nil)
        #expect(failed.report.exports[0].artifact == nil)
        let before = await probe.calls
        let repeated = try await executor.run(request)
        #expect(repeated.report.nodes[1].outcome.reused)
        #expect(await probe.calls == before + 1)
        request.nodes[0].configuration = .object([:])
        let repaired = try await executor.run(request)
        #expect(repaired.report.allTasksSucceeded)
        #expect(repaired.report.nodes[1].outcome.reused)
        try record(failed.report, "platform-partial-failure")
    }
    @Test func invalidGraphsAndCorruptInputsCannotExecuteTasks() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let probe = WorkflowProbe(), store = try VivoArtifactStore(rootURL: root)
        let executor = try VivoWorkflowExecutor(store: store, registry: registry(probe, id: identity()))
        var request = recipe(); request.nodes[0].inputs["value"] = .output(node: "c", port: "value")
        do { _ = try await executor.run(request); Issue.record("cycle accepted") } catch {}
        #expect(await probe.calls == 0)
        #expect(try await store.list().isEmpty)
        let input = try await store.put(data: Data("7".utf8), kind: "test.value", mediaType: "application/json")
        try Data("8".utf8).write(to: root.appendingPathComponent(input.objectPath))
        request = recipe(); request.artifacts[0].source = .stored(kind: "test.value", fingerprint: input.fingerprint)
        do { _ = try await executor.run(request); Issue.record("corrupt input executed") } catch {}
        #expect(await probe.calls == 0)
    }
    @Test func realMolecularOperationsShareAnUpstreamHamiltonian() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let store = try VivoArtifactStore(rootURL: root), request = try VivoWorkflowTemplates.molecularAnalysis()
        let executor = VivoWorkflowExecutor(store: store, registry: registry), workflow = VivoChemistryWorkflow(store: store)
        let result = try await executor.run(request)
        #expect(result.report.allTasksSucceeded && result.report.nodes.count == 6)
        let exports = Dictionary(uniqueKeysWithValues: result.report.exports.map { ($0.name, $0.artifact!) })
        let h = try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self, from: await workflow.payload(artifact: exports["hamiltonian"]!.artifact, expectedKind: "vivo.embedded-hamiltonian"))
        let ci = try VivoCanonicalJSON.decode(VivoManyBodySolverResult.self, from: await workflow.payload(artifact: exports["fci"]!.artifact, expectedKind: "vivo.many-body-result"))
        #expect(abs(ci.energyHartree - (try VivoDirectCI.solve(h).roots[0].energyHartree)) < 1e-8)
        let cached = try await executor.run(request)
        #expect(cached.report.nodes.allSatisfy { $0.outcome.reused })
        try record(result.report, "platform-molecular-analysis")
    }
    @Test func realMetalSegmentsContinueWithoutResettingAcceptedState() async throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: identity())
        let store = try VivoArtifactStore(rootURL: root), request = try VivoWorkflowTemplates.molecularDynamics()
        let executor = VivoWorkflowExecutor(store: store, registry: registry), workflow = VivoChemistryWorkflow(store: store)
        let result = try await executor.run(request)
        #expect(result.report.allTasksSucceeded && result.report.nodes.count == 2)
        #expect(result.report.maximumAdmittedMetalTasks == 1)
        guard let export = result.report.exports.first(where: { $0.name == "checkpoint" })?.artifact else { Issue.record("missing MD checkpoint"); return }
        let checkpoint = try VivoCanonicalJSON.decode(VivoMDCheckpoint.self, from: await workflow.payload(artifact: export.artifact, expectedKind: "vivo.md-checkpoint"))
        #expect(checkpoint.acceptedStep == 20)
        #expect(abs(checkpoint.timePS - 0.004) < 1e-12)
        let cached = try await executor.run(request)
        #expect(cached.report.nodes.allSatisfy { $0.outcome.reused })
        try record(result.report, "platform-metal-segments")
    }
}
