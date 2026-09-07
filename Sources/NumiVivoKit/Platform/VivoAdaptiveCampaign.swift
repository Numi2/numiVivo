import Foundation
import Dispatch

public struct VivoCampaignArtifactBinding: Codable, Sendable, Equatable {
    public let artifactIdentifier: String
    public let actionIdentifier: String
    public let exportName: String
    public init(artifactIdentifier: String,actionIdentifier: String,exportName: String) {
        self.artifactIdentifier = artifactIdentifier; self.actionIdentifier = actionIdentifier; self.exportName = exportName
    }
}
public struct VivoAdaptiveCampaignAction: Codable, Sendable, Equatable {
    public var proposal: VivoRefinementActionProposal
    public var recipe: VivoWorkflowRecipe
    public var metricExport: String
    public var bindings: [VivoCampaignArtifactBinding]
    public init(proposal: VivoRefinementActionProposal,recipe: VivoWorkflowRecipe,metricExport: String = "metric",
                bindings: [VivoCampaignArtifactBinding] = []) {
        self.proposal = proposal; self.recipe = recipe; self.metricExport = metricExport; self.bindings = bindings
    }
}
public struct VivoAdaptiveCampaignRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/adaptive-refinement-campaign/v1"
    public var schema: String
    public var identifier: String
    public var criteria: [VivoRefinementCriterion]
    public var actions: [VivoAdaptiveCampaignAction]
    public var maximumActions: Int
    /// Admission units, not a wall-clock promise or a count of GPU instructions.
    /// Individual native recipes retain their own numerical and memory budgets.
    public var maximumDeclaredWorkUnits: Int
    public init(identifier: String,criteria: [VivoRefinementCriterion],actions: [VivoAdaptiveCampaignAction],
                maximumActions: Int = 32,maximumDeclaredWorkUnits: Int = 1_000_000_000) {
        schema = Self.schema; self.identifier = identifier; self.criteria = criteria; self.actions = actions
        self.maximumActions = maximumActions; self.maximumDeclaredWorkUnits = maximumDeclaredWorkUnits
    }
}
public struct VivoAdaptiveCampaignPlan: Codable, Sendable, Equatable {
    public let requestFingerprint: VivoFingerprint
    public let implementationFingerprint: VivoFingerprint
    /// Bound producer outputs use kind-checked symbolic fingerprints here.
    /// Actual recipes are resolved and fingerprinted again before execution.
    public let actionPlans: [String: VivoWorkflowPlan]
    public let interpretation: String
}
public enum VivoAdaptiveCampaignPlanner {
    public static func compile(_ request: VivoAdaptiveCampaignRequest,registry: VivoWorkflowRegistry) throws -> VivoAdaptiveCampaignPlan {
        guard request.schema == VivoAdaptiveCampaignRequest.schema, !request.identifier.isEmpty, request.identifier.utf8.count <= 512,
              (1...128).contains(request.maximumActions), request.maximumDeclaredWorkUnits > 0 else {
            throw VivoChemistryError.invalid("adaptive campaign identity or admission budget")
        }
        try VivoAdaptiveRefinementPolicy.validate(criteria: request.criteria,proposals: request.actions.map(\.proposal))
        let actions = Dictionary(uniqueKeysWithValues: request.actions.map { ($0.proposal.identifier,$0) })
        var plans: [String: VivoWorkflowPlan] = [:]
        while plans.count < actions.count {
            let ready = actions.keys.sorted().filter { plans[$0] == nil && actions[$0]!.proposal.prerequisites.allSatisfy { plans[$0] != nil } }
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("campaign dependency cycle") }
            for identifier in ready {
                let action = actions[identifier]!
                var recipe = action.recipe
                // Recipes admit only immutable inline JSON or stored artifacts.
                guard Set(action.bindings.map(\.artifactIdentifier)).count == action.bindings.count else {
                    throw VivoChemistryError.invalid("duplicate campaign artifact binding")
                }
                for binding in action.bindings {
                    guard action.proposal.prerequisites.contains(binding.actionIdentifier),
                          let producer = actions[binding.actionIdentifier], let plan = plans[binding.actionIdentifier],
                          let output = producer.recipe.outputs.first(where: { $0.name == binding.exportName }),
                          let node = plan.nodes.first(where: { $0.identifier == output.node }),
                          let port = try registry.definition(node.operation).operation.outputs.first(where: { $0.name == output.port }),
                          !recipe.artifacts.contains(where: { $0.identifier == binding.artifactIdentifier }) else {
                        throw VivoChemistryError.invalid("campaign binding must name a declared prerequisite export and a new artifact input")
                    }
                    let symbolic = try VivoCanonicalJSON.fingerprint(Data("campaign-symbolic:\(binding.actionIdentifier):\(binding.exportName)".utf8))
                    recipe.artifacts.append(.init(identifier: binding.artifactIdentifier,source: .stored(kind: port.kind,fingerprint: symbolic)))
                }
                let plan = try VivoWorkflowPlanner.compile(recipe,registry: registry)
                guard let exported = recipe.outputs.first(where: { $0.name == action.metricExport }),
                      let node = recipe.nodes.first(where: { $0.identifier == exported.node }),
                      VivoPlatformAdaptiveOperations.metricOperations.contains(node.operation),
                      let metric = try registry.definition(node.operation).operation.outputs.first(where: { $0.name == exported.port }),
                      metric.kind == VivoPlatformAdaptiveOperations.metricKind else {
                    throw VivoChemistryError.invalid("campaign metric must be produced by a native reconstructing evidence adapter")
                }
                // A confirmation result may not become a tuning input for its
                // own criterion, even indirectly through an action dependency.
                if action.proposal.role == .discovery {
                    var ancestors = Set(action.proposal.prerequisites), frontier = Array(ancestors)
                    while let ancestor = frontier.popLast() {
                        let parent = actions[ancestor]!
                        if parent.proposal.role == .confirmation && parent.proposal.criterionIdentifier == action.proposal.criterionIdentifier {
                            throw VivoChemistryError.invalid("held-out confirmation feeds discovery for the same criterion")
                        }
                        for dependency in parent.proposal.prerequisites where ancestors.insert(dependency).inserted { frontier.append(dependency) }
                    }
                }
                plans[identifier] = plan
            }
        }
        return .init(requestFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request)),
            implementationFingerprint: registry.implementationFingerprint,actionPlans: plans,
            interpretation: "Predeclared native workflow actions and kind-checked dependencies. Symbolic plans do not execute work or qualify a model. Observed fresh execution latency calibrates ordering only; each metric retains its own physical unit, target binding and native acceptance checks.")
    }
}
public struct VivoAdaptiveCampaignStep: Codable, Sendable, Equatable {
    public let decision: VivoRefinementDecision
    public let record: VivoRefinementActionRecord
    public let workflowRunArtifact: VivoFingerprint?
    public let resolvedRecipeFingerprint: VivoFingerprint?
}
public struct VivoAdaptiveCampaignCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/adaptive-refinement-checkpoint/v1"
    public let schema: String
    public let request: VivoAdaptiveCampaignRequest
    public let implementationFingerprint: VivoFingerprint
    public let previousCheckpoint: VivoFingerprint?
    public let steps: [VivoAdaptiveCampaignStep]
}
public enum VivoAdaptiveCampaignTermination: String, Codable, Sendable {
    case criteriaSatisfied, confirmationFailed, actionLimit, declaredWorkLimit, actionsExhausted, executionIncomplete
}
public struct VivoAdaptiveCampaignReport: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/adaptive-refinement-report/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let implementationFingerprint: VivoFingerprint
    public let checkpointArtifact: VivoFingerprint
    public let steps: [VivoAdaptiveCampaignStep]
    public let finalDecision: VivoRefinementDecision
    public let chargedDeclaredWorkUnits: Int
    public let unexecutedActionIdentifiers: [String]
    public let termination: VivoAdaptiveCampaignTermination
    public let interpretation: String
}
public struct VivoAdaptiveCampaignExecution: Sendable {
    public let report: VivoAdaptiveCampaignReport
    public let artifact: VivoStoredArtifact
}

/// Executes the selected complete recipes through VivoWorkflowExecutor. No
/// production Hamiltonian is mutated and no asynchronous work outlives a call.
/// Checkpoints preserve timing observations for deterministic policy replay.
public actor VivoAdaptiveCampaignExecutor {
    private let store: VivoArtifactStore
    private let registry: VivoWorkflowRegistry
    private let executor: VivoWorkflowExecutor
    private let payloads: VivoChemistryWorkflow
    private var activeInvocation = false
    public init(store: VivoArtifactStore,registry: VivoWorkflowRegistry) {
        self.store = store; self.registry = registry
        executor = VivoWorkflowExecutor(store: store,registry: registry); payloads = VivoChemistryWorkflow(store: store)
    }
    private func metric(_ execution: VivoWorkflowExecution,action: VivoAdaptiveCampaignAction) async throws -> VivoRefinementMetricEvidence {
        guard execution.report.allTasksSucceeded,
              let receipt = execution.report.exports.first(where: { $0.name == action.metricExport })?.artifact,
              receipt.kind == VivoPlatformAdaptiveOperations.metricKind else {
            throw VivoChemistryError.convergence("adaptive action has failed/blocked work or no native metric export")
        }
        let data = try await payloads.payload(artifact: receipt.artifact,expectedKind: receipt.kind)
        let result = try VivoCanonicalJSON.decode(VivoRefinementMetricEvidence.self,from: data)
        try result.validate(); return result
    }
    private func resolved(_ action: VivoAdaptiveCampaignAction,executions: [String: VivoWorkflowExecution]) throws -> VivoWorkflowRecipe {
        var recipe = action.recipe
        for binding in action.bindings {
            guard let receipt = executions[binding.actionIdentifier]?.report.exports.first(where: { $0.name == binding.exportName })?.artifact,
                  !recipe.artifacts.contains(where: { $0.identifier == binding.artifactIdentifier }) else {
                throw VivoChemistryError.invalid("adaptive prerequisite output is unavailable or collides with an input")
            }
            recipe.artifacts.append(.init(identifier: binding.artifactIdentifier,source: .stored(kind: receipt.kind,fingerprint: receipt.artifact)))
        }
        return recipe
    }
    public func resume(_ checkpointID: VivoFingerprint) async throws -> VivoAdaptiveCampaignExecution {
        let descriptor = try await store.descriptor(for: checkpointID)
        guard descriptor.kind == "vivo.adaptive-refinement-checkpoint" else { throw VivoChemistryError.invalid("expected campaign checkpoint") }
        let checkpoint = try VivoCanonicalJSON.decode(VivoAdaptiveCampaignCheckpoint.self,from: await store.data(for: checkpointID,verify: true))
        return try await run(checkpoint.request,resuming: checkpointID)
    }
    public func run(_ request: VivoAdaptiveCampaignRequest,resuming checkpointID: VivoFingerprint? = nil) async throws -> VivoAdaptiveCampaignExecution {
        guard !activeInvocation else { throw VivoChemistryError.invalid("campaign executor already has an active invocation") }
        activeInvocation = true; defer { activeInvocation = false }
        let plan = try VivoAdaptiveCampaignPlanner.compile(request,registry: registry)
        let actions = Dictionary(uniqueKeysWithValues: request.actions.map { ($0.proposal.identifier,$0) })
        var steps: [VivoAdaptiveCampaignStep] = [], executions: [String: VivoWorkflowExecution] = [:]
        var charged = 0, previousCheckpoint = checkpointID
        if let checkpointID {
            let descriptor = try await store.descriptor(for: checkpointID)
            guard descriptor.kind == "vivo.adaptive-refinement-checkpoint" else { throw VivoChemistryError.invalid("campaign checkpoint kind") }
            let saved = try VivoCanonicalJSON.decode(VivoAdaptiveCampaignCheckpoint.self,from: await store.data(for: checkpointID,verify: true))
            guard saved.schema == VivoAdaptiveCampaignCheckpoint.schema, saved.request == request,
                  saved.implementationFingerprint == registry.implementationFingerprint, saved.steps.count <= request.maximumActions else {
                throw VivoChemistryError.invalid("campaign checkpoint request/implementation binding")
            }
            for step in saved.steps {
                let decision = try VivoAdaptiveRefinementPolicy.decide(criteria: request.criteria,proposals: request.actions.map(\.proposal),
                    records: steps.map(\.record),remainingDeclaredWorkUnits: request.maximumDeclaredWorkUnits-charged)
                guard step.decision == decision, decision.selectedActionIdentifier == step.record.actionIdentifier,
                      let action = actions[step.record.actionIdentifier] else { throw VivoChemistryError.invalid("campaign decision history does not reconstruct") }
                if let run = step.workflowRunArtifact {
                    let execution = try await executor.verify(run)
                    let recipe = try resolved(action,executions: executions)
                    let rebuiltPlan = try VivoWorkflowPlanner.compile(recipe,registry: registry)
                    guard execution.report.recipeFingerprint == rebuiltPlan.recipeFingerprint,
                          step.resolvedRecipeFingerprint == rebuiltPlan.recipeFingerprint else {
                        throw VivoChemistryError.invalid("campaign recipe changed between selection and execution")
                    }
                    executions[step.record.actionIdentifier] = execution
                    if step.record.failure == nil {
                        guard step.record.metric == (try await metric(execution,action: action)) else {
                            throw VivoChemistryError.invalid("campaign metric differs from reconstructed native output")
                        }
                    } else {
                        guard !execution.report.allTasksSucceeded || step.record.metric == nil else {
                            throw VivoChemistryError.invalid("campaign failure record is inconsistent")
                        }
                    }
                } else {
                    guard step.record.failure != nil, step.record.metric == nil else { throw VivoChemistryError.invalid("campaign success without a workflow run") }
                }
                charged += action.proposal.declaredWorkUnits; steps.append(step)
            }
        }
        func checkpoint(_ records: [VivoAdaptiveCampaignStep],previous: VivoFingerprint?) async throws -> VivoFingerprint {
            let value = VivoAdaptiveCampaignCheckpoint(schema: VivoAdaptiveCampaignCheckpoint.schema,request: request,
                implementationFingerprint: registry.implementationFingerprint,previousCheckpoint: previous,steps: records)
            let descriptor = try await store.put(data: VivoCanonicalJSON.encode(value),kind: "vivo.adaptive-refinement-checkpoint",mediaType: "application/json")
            _ = try await store.setReference("adaptive-campaign-"+plan.requestFingerprint.hex,to: descriptor)
            return descriptor.fingerprint
        }
        previousCheckpoint = try await checkpoint(steps,previous: previousCheckpoint)
        var decision: VivoRefinementDecision
        var termination: VivoAdaptiveCampaignTermination
        while true {
            try Task.checkCancellation()
            decision = try VivoAdaptiveRefinementPolicy.decide(criteria: request.criteria,proposals: request.actions.map(\.proposal),
                records: steps.map(\.record),remainingDeclaredWorkUnits: request.maximumDeclaredWorkUnits-charged)
            if decision.allCriteriaSatisfied { termination = .criteriaSatisfied; break }
            if decision.criteria.contains(where: \.confirmationFailed) { termination = .confirmationFailed; break }
            if steps.count >= request.maximumActions { termination = .actionLimit; break }
            guard let selected = decision.selectedActionIdentifier else {
                let withoutBudget = try VivoAdaptiveRefinementPolicy.decide(criteria: request.criteria,proposals: request.actions.map(\.proposal),
                    records: steps.map(\.record),remainingDeclaredWorkUnits: Int.max)
                termination = withoutBudget.selectedActionIdentifier != nil ? .declaredWorkLimit
                    : (steps.contains(where: { $0.record.failure != nil }) ? .executionIncomplete : .actionsExhausted)
                break
            }
            let action = actions[selected]!, start = DispatchTime.now().uptimeNanoseconds
            charged += action.proposal.declaredWorkUnits
            var runID: VivoFingerprint?, recipeID: VivoFingerprint?, observation: VivoRefinementMetricEvidence?, failure: String?
            var uncached = false
            do {
                let recipe = try resolved(action,executions: executions)
                let localPlan = try VivoWorkflowPlanner.compile(recipe,registry: registry)
                recipeID = localPlan.recipeFingerprint
                let execution = try await executor.run(recipe)
                runID = execution.artifact.fingerprint; executions[selected] = execution
                uncached = execution.report.nodes.allSatisfy { node in
                    if case .succeeded(_,_,_,let reused) = node.outcome { return !reused }; return false
                }
                observation = try await metric(execution,action: action)
            } catch is CancellationError { throw CancellationError() }
            catch { failure = String(String(describing: error).prefix(2048)) }
            let elapsed = Double(DispatchTime.now().uptimeNanoseconds-start)/1e9
            let record = VivoRefinementActionRecord(actionIdentifier: selected,metric: observation,elapsedSeconds: elapsed,
                entirelyUncachedExecution: uncached && failure == nil,failure: failure)
            steps.append(.init(decision: decision,record: record,workflowRunArtifact: runID,resolvedRecipeFingerprint: recipeID))
            previousCheckpoint = try await checkpoint(steps,previous: previousCheckpoint)
        }
        let done = Set(steps.map { $0.record.actionIdentifier })
        let report = VivoAdaptiveCampaignReport(schema: VivoAdaptiveCampaignReport.schema,requestFingerprint: plan.requestFingerprint,
            implementationFingerprint: registry.implementationFingerprint,checkpointArtifact: previousCheckpoint!,steps: steps,
            finalDecision: decision,chargedDeclaredWorkUnits: charged,unexecutedActionIdentifiers: actions.keys.filter { !done.contains($0) }.sorted(),
            termination: termination,
            interpretation: "Native recipe execution under a predeclared criterion/action pool. Unit- and observable-bound metrics retain native qualification; statistical uncertainty, electronic sensitivity and local kinetic uncertainty are not combined into one invented error bar. Failed confirmation stops tuning. Fresh measured latency calibrates ordering only; cached execution does not calibrate fresh cost. Declared admission work is not elapsed time. A satisfied campaign does not establish unexamined model adequacy, chemical accuracy, rate validity or performance superiority.")
        let artifact = try await store.put(data: VivoCanonicalJSON.encode(report),kind: "vivo.adaptive-refinement-report",mediaType: "application/json")
        return .init(report: report,artifact: artifact)
    }
}
