import Foundation

public enum VivoWorkflowNodeOutcome: Codable, Sendable, Equatable {
    case succeeded(task: VivoFingerprint, receipt: VivoFingerprint, outputs: [VivoChemistryOutputReceipt], reused: Bool)
    case failed(category: String, message: String)
    case blocked(dependencies: [String])
    public var outputs: [VivoChemistryOutputReceipt]? {
        if case .succeeded(_, _, let outputs, _) = self { return outputs }; return nil
    }
    public var reused: Bool {
        if case .succeeded(_, _, _, let reused) = self { return reused }; return false
    }
}
public struct VivoWorkflowNodeRecord: Codable, Sendable, Equatable {
    public let identifier: String
    public let outcome: VivoWorkflowNodeOutcome
}
public struct VivoWorkflowExportRecord: Codable, Sendable, Equatable {
    public let name: String
    public let node: String
    public let port: String
    /// Nil when its producer failed or was blocked. No fabricated empty output.
    public let artifact: VivoChemistryOutputReceipt?
}
public struct VivoWorkflowRunReport: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/workflow-run/v1"
    public let schema: String
    public let recipeFingerprint: VivoFingerprint
    public let recipeArtifact: VivoFingerprint
    public let planArtifact: VivoFingerprint
    public let implementationFingerprint: VivoFingerprint
    public let nodes: [VivoWorkflowNodeRecord]
    public let exports: [VivoWorkflowExportRecord]
    public let maximumAdmittedTasks: Int
    public let maximumAdmittedMetalTasks: Int
    public let maximumReservedBytes: Int
    public let allTasksSucceeded: Bool
    public let meaning: String
}
public struct VivoWorkflowExecution: Sendable {
    public let report: VivoWorkflowRunReport
    public let artifact: VivoStoredArtifact
}

/// Orchestration only. Every accepted numerical output is executed, validated
/// and content-addressed by the existing VivoChemistryWorkflow authority.
/// Per-invocation admission limits include configured numerical + input + output
/// reservations; they do not measure RSS or constrain unregistered applications.
public struct VivoWorkflowExecutor: Sendable {
    private let store: VivoArtifactStore
    private let registry: VivoWorkflowRegistry
    public init(store: VivoArtifactStore, registry: VivoWorkflowRegistry) { self.store = store; self.registry = registry }

    public func run(_ recipe: VivoWorkflowRecipe) async throws -> VivoWorkflowExecution {
        let plan = try VivoWorkflowPlanner.compile(recipe, registry: registry)
        try Task.checkCancellation()
        let workflow = VivoChemistryWorkflow(store: store)
        var inputs: [String: VivoChemistryTaskInput] = [:]
        // Check *all* external content before importing inline inputs or running
        // independent branches. A corrupt ancestor is not an ordinary cache miss.
        for artifact in recipe.artifacts {
            if case .stored(let kind, let fingerprint) = artifact.source {
                let descriptor = try await store.descriptor(for: fingerprint)
                guard descriptor.byteCount <= UInt64(recipe.policy.maximumReservedBytes) else {
                    throw VivoChemistryError.resourceLimit("external workflow artifact exceeds admission cap")
                }
                _ = try await workflow.payload(artifact: fingerprint, expectedKind: kind)
                inputs[artifact.identifier] = .init(name: artifact.identifier, artifact: fingerprint, kind: kind)
            }
        }
        for artifact in recipe.artifacts.sorted(by: { $0.identifier < $1.identifier }) {
            if case .json(let kind, let payload) = artifact.source {
                let value = try await store.put(data: VivoCanonicalJSON.encode(payload), kind: kind, mediaType: "application/json")
                inputs[artifact.identifier] = .init(name: artifact.identifier, artifact: value.fingerprint, kind: kind)
            }
        }
        let recipeArtifact = try await store.put(data: VivoCanonicalJSON.encode(recipe), kind: "vivo.workflow-recipe", mediaType: "application/json")
        let planArtifact = try await store.put(data: VivoCanonicalJSON.encode(plan), kind: "vivo.workflow-plan", mediaType: "application/json")
        let nodes = Dictionary(uniqueKeysWithValues: recipe.nodes.map { ($0.identifier, $0) })
        let planned = Dictionary(uniqueKeysWithValues: plan.nodes.map { ($0.identifier, $0) })
        var records: [String: VivoWorkflowNodeOutcome] = [:]
        var peakTasks = 0, peakMetal = 0, peakBytes = 0
        while records.count < nodes.count {
            try Task.checkCancellation()
            // Preserve blocked descendants, while unrelated branches continue.
            var marked = true
            while marked {
                marked = false
                for id in plan.topologicalOrder where records[id] == nil {
                    let failed = planned[id]!.dependencies.filter { records[$0] != nil && records[$0]!.outputs == nil }
                    if !failed.isEmpty { records[id] = .blocked(dependencies: failed); marked = true }
                }
            }
            let ready = plan.topologicalOrder.filter { id in
                records[id] == nil && planned[id]!.dependencies.allSatisfy { records[$0]?.outputs != nil }
            }
            if ready.isEmpty {
                guard records.count == nodes.count else { throw VivoChemistryError.invalid("validated workflow scheduler made no progress") }
                break
            }
            var batch: [String] = [], reserved = 0, metal = 0
            for id in ready {
                let p = planned[id]!, isMetal = p.numericalBackend.hasPrefix("metal")
                if batch.count >= recipe.policy.maximumConcurrentTasks { break }
                if isMetal && metal >= recipe.policy.maximumConcurrentMetalTasks { continue }
                guard p.reservedBytes <= recipe.policy.maximumReservedBytes - reserved else { continue }
                batch.append(id); reserved += p.reservedBytes; if isMetal { metal += 1 }
            }
            guard !batch.isEmpty else { throw VivoChemistryError.resourceLimit("no ready node fits the declared workflow admission policy") }
            peakTasks = max(peakTasks, batch.count); peakMetal = max(peakMetal, metal); peakBytes = max(peakBytes, reserved)
            // Resolve ports before task creation. No child task mutates shared
            // dictionaries. Completed task receipts, not closure state, feed edges.
            var jobs: [(String, VivoChemistryTask, VivoChemistryOperation)] = []
            for id in batch {
                let node = nodes[id]!, operation = try registry.definition(node.operation).operation
                let resolved = try node.inputs.keys.sorted().map { port -> VivoChemistryTaskInput in
                    switch node.inputs[port]! {
                    case .artifact(let identifier):
                        guard let input = inputs[identifier] else { throw VivoChemistryError.invalid("unresolved workflow input") }
                        return .init(name: port, artifact: input.artifact, kind: input.kind)
                    case .output(let producer, let output):
                        guard let input = records[producer]?.outputs?.first(where: { $0.name == output }) else {
                            throw VivoChemistryError.invalid("unresolved successful workflow output")
                        }
                        return .init(name: port, artifact: input.artifact, kind: input.kind)
                    }
                }
                let task = VivoChemistryTask(operation: operation.identifier, version: operation.version,
                    implementationFingerprint: operation.implementationFingerprint, inputs: resolved,
                    configuration: node.configuration, outputs: operation.outputs, resources: node.resources)
                jobs.append((id, task, operation))
            }
            let completed = try await withThrowingTaskGroup(of: VivoWorkflowNodeRecord.self) { group in
                for (id, task, operation) in jobs {
                    group.addTask {
                        try Task.checkCancellation()
                        do {
                            let result = try await workflow.run(task, using: operation)
                            try Task.checkCancellation()
                            return .init(identifier: id, outcome: .succeeded(task: result.taskFingerprint,
                                receipt: result.receiptFingerprint, outputs: result.outputs, reused: result.reused))
                        } catch is CancellationError { throw CancellationError() }
                        catch {
                            let category: String
                            switch error {
                            case VivoChemistryError.convergence: category = "numerical-convergence"
                            case VivoChemistryError.resourceLimit: category = "resource-limit"
                            case is VivoArtifactStoreError: category = "artifact-integrity-or-io"
                            default: category = "operation-failed"
                            }
                            return .init(identifier: id, outcome: .failed(category: category,
                                message: String(String(describing: error).prefix(2048))))
                        }
                    }
                }
                var values: [VivoWorkflowNodeRecord] = []
                for try await value in group { values.append(value) }
                return values
            }
            for value in completed { records[value.identifier] = value.outcome }
        }
        try Task.checkCancellation()
        let report = VivoWorkflowRunReport(schema: VivoWorkflowRunReport.schema, recipeFingerprint: plan.recipeFingerprint,
            recipeArtifact: recipeArtifact.fingerprint, planArtifact: planArtifact.fingerprint,
            implementationFingerprint: registry.implementationFingerprint,
            nodes: records.keys.sorted().map { .init(identifier: $0, outcome: records[$0]!) },
            exports: recipe.outputs.sorted(by: { $0.name < $1.name }).map { output in
                .init(name: output.name, node: output.node, port: output.port,
                      artifact: records[output.node]?.outputs?.first(where: { $0.name == output.port }))
            }, maximumAdmittedTasks: peakTasks, maximumAdmittedMetalTasks: peakMetal, maximumReservedBytes: peakBytes,
            allTasksSucceeded: records.values.allSatisfy { $0.outputs != nil },
            meaning: "execution and declared per-operation output validation; failed and blocked work retained; reuse revalidates existing task outputs; no change to scientific evidence classifications and no inferred cross-scale rate or efficacy")
        let artifact = try await store.put(data: VivoCanonicalJSON.encode(report), kind: "vivo.workflow-run", mediaType: "application/json")
        return .init(report: report, artifact: artifact)
    }

    /// Verify a stored run's manifest, reconstruct each task's input identity,
    /// and rerun the recipe under the current implementation. This may execute
    /// previously failed tasks; it never accepts a user-edited success flag.
    public func verify(_ artifact: VivoFingerprint) async throws -> VivoWorkflowExecution {
        let descriptor = try await store.descriptor(for: artifact)
        guard descriptor.kind == "vivo.workflow-run" else { throw VivoChemistryError.invalid("expected workflow run artifact") }
        let saved = try VivoCanonicalJSON.decode(VivoWorkflowRunReport.self, from: await store.data(for: artifact, verify: true))
        guard saved.schema == VivoWorkflowRunReport.schema, saved.implementationFingerprint == registry.implementationFingerprint else {
            throw VivoChemistryError.invalid("workflow run implementation/schema differs; run its recipe as a new execution instead")
        }
        let recipeDescriptor = try await store.descriptor(for: saved.recipeArtifact)
        guard recipeDescriptor.kind == "vivo.workflow-recipe" else { throw VivoChemistryError.invalid("run recipe artifact kind") }
        let recipe = try VivoCanonicalJSON.decode(VivoWorkflowRecipe.self, from: await store.data(for: saved.recipeArtifact, verify: true))
        let plan = try VivoWorkflowPlanner.compile(recipe, registry: registry)
        let planDescriptor = try await store.descriptor(for: saved.planArtifact)
        guard planDescriptor.kind == "vivo.workflow-plan", saved.recipeFingerprint == plan.recipeFingerprint,
              try VivoCanonicalJSON.decode(VivoWorkflowPlan.self, from: await store.data(for: saved.planArtifact, verify: true)) == plan else {
            throw VivoChemistryError.invalid("run recipe/plan binding")
        }
        let rebuilt = try await run(recipe)
        guard saved.nodes.count == rebuilt.report.nodes.count, saved.exports == rebuilt.report.exports,
              saved.allTasksSucceeded == rebuilt.report.allTasksSucceeded, saved.meaning == rebuilt.report.meaning,
              saved.maximumAdmittedTasks == rebuilt.report.maximumAdmittedTasks,
              saved.maximumAdmittedMetalTasks == rebuilt.report.maximumAdmittedMetalTasks,
              saved.maximumReservedBytes == rebuilt.report.maximumReservedBytes else {
            throw VivoChemistryError.invalid("run manifest differs on task reconstruction")
        }
        for (a, b) in zip(saved.nodes, rebuilt.report.nodes) {
            guard a.identifier == b.identifier else { throw VivoChemistryError.invalid("run node identity changed") }
            switch (a.outcome, b.outcome) {
            case (.succeeded(let ta, let ra, let oa, _), .succeeded(let tb, let rb, let ob, _)):
                guard ta == tb, ra == rb, oa == ob else { throw VivoChemistryError.invalid("run task, receipt or outputs changed") }
            default:
                guard a.outcome == b.outcome else { throw VivoChemistryError.invalid("run failure/blocked status changed; new execution is not the saved report") }
            }
        }
        return rebuilt
    }
}
