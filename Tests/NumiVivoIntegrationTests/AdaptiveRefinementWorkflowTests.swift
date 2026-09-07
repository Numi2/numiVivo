import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct AdaptiveRefinementWorkflowTests {
    private func registry() throws -> VivoWorkflowRegistry {
        try VivoPlatformOperations.registry(implementationFingerprint: VivoCanonicalJSON.fingerprint(Data("adaptive-integration-v1".utf8)))
    }
    private func root() -> URL { FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-adaptive-\(UUID().uuidString)") }
    @Test func bothNativeCampaignsPlan() throws {
        for name in ["equilibrium-correction","electronic-crosscheck"] {
            let request = try VivoAdaptiveCampaignExamples.make(name)
            let plan = try VivoAdaptiveCampaignPlanner.compile(request,registry: registry())
            #expect(plan.actionPlans.count == 2)
            #expect(plan.actionPlans.values.allSatisfy { !$0.nodes.isEmpty })
        }
    }
    @Test func realEquilibriumCampaignExecutesCachesAndResumes() async throws {
        let directory = root(); defer { try? FileManager.default.removeItem(at: directory) }
        let store = try VivoArtifactStore(rootURL: directory), registry = try registry()
        let executor = VivoAdaptiveCampaignExecutor(store: store,registry: registry)
        let request = try VivoAdaptiveCampaignExamples.make("equilibrium-correction")
        let first = try await executor.run(request)
        #expect(first.report.termination == .criteriaSatisfied)
        #expect(first.report.steps.count == 2)
        #expect(first.report.steps.allSatisfy { $0.record.failure == nil && $0.record.metric?.nativeChecksPassed == true })
        let replay = try await executor.resume(first.report.checkpointArtifact)
        #expect(replay.report.termination == .criteriaSatisfied)
        #expect(replay.report.steps == first.report.steps)
        let freshRequestReplay = try await executor.run(request)
        #expect(freshRequestReplay.report.termination == .criteriaSatisfied)
        #expect(freshRequestReplay.report.steps.allSatisfy { !$0.record.entirelyUncachedExecution })
        #expect(first.report.steps.map(\.record.metric) == freshRequestReplay.report.steps.map(\.record.metric))
    }
    @Test func electronicCampaignUsesExistingRefinementBackend() async throws {
        let directory = root(); defer { try? FileManager.default.removeItem(at: directory) }
        let store = try VivoArtifactStore(rootURL: directory)
        let request = try VivoAdaptiveCampaignExamples.make("electronic-crosscheck")
        #expect(!request.criteria[0].requiresDisjointConfirmationSources)
        let execution = try await VivoAdaptiveCampaignExecutor(store: store,registry: registry()).run(request)
        #expect(execution.report.termination == .criteriaSatisfied)
        #expect(execution.report.steps.map(\.record.actionIdentifier) == ["direct","selected"])
    }
    @Test func wrongObservableCannotSatisfyCampaign() async throws {
        let directory = root(); defer { try? FileManager.default.removeItem(at: directory) }
        var request = try VivoAdaptiveCampaignExamples.make("equilibrium-correction")
        request.criteria[0].observableFingerprint = String(repeating: "f",count: 64)
        let executor = try VivoAdaptiveCampaignExecutor(store: VivoArtifactStore(rootURL: directory),registry: registry())
        await #expect(throws: VivoChemistryError.self) { try await executor.run(request) }
    }
    @Test func unsampledTargetCannotPassNativeMetric() async throws {
        let registry = try registry()
        let operation = try registry.definition("vivo.platform.mbar-target-refinement").operation
        let sampled = VivoAdaptiveCampaignExamples.equilibriumRequest(seed: 81)
        let request = VivoMBARTargetRefinementRequest(sampledStateIdentifiers: sampled.sampledStateIdentifiers,
            targetIdentifiers: sampled.targetIdentifiers,matchingSampledStateIndices: [0,nil],
            binEdges: sampled.binEdges,coordinateUnit: sampled.coordinateUnit,samples: sampled.samples,configuration: sampled.configuration)
        let input = ["request":try VivoCanonicalJSON.encode(request)]
        let outputs = try await operation.execute(.object([:]),input,.init())
        try operation.validateOutputs(.object([:]),input,outputs,.init())
        let metric = try VivoCanonicalJSON.decode(VivoRefinementMetricEvidence.self,from: #require(outputs["metric"]))
        #expect(!metric.nativeChecksPassed)
    }
    @Test func untrustedMetricProducerAndHoldoutFeedbackReject() throws {
        var request = try VivoAdaptiveCampaignExamples.make("equilibrium-correction")
        request.actions[0].metricExport = "result"
        #expect(throws: VivoChemistryError.self) { try VivoAdaptiveCampaignPlanner.compile(request,registry: registry()) }
        request = try VivoAdaptiveCampaignExamples.make("equilibrium-correction")
        var tune = request.actions[0]
        tune.proposal.identifier = "tune-on-confirmation"
        tune.proposal.prerequisites = ["confirmation"]
        request.actions.append(tune)
        #expect(throws: VivoChemistryError.self) { try VivoAdaptiveCampaignPlanner.compile(request,registry: registry()) }
    }
    @Test func exhaustedAdmissionBudgetNeverMeansSatisfied() async throws {
        let directory = root(); defer { try? FileManager.default.removeItem(at: directory) }
        var request = try VivoAdaptiveCampaignExamples.make("equilibrium-correction")
        request.maximumDeclaredWorkUnits = 1
        let result = try await VivoAdaptiveCampaignExecutor(store: VivoArtifactStore(rootURL: directory),registry: registry()).run(request)
        #expect(result.report.termination == .declaredWorkLimit)
        #expect(result.report.steps.isEmpty && !result.report.finalDecision.allCriteriaSatisfied)
    }
}
