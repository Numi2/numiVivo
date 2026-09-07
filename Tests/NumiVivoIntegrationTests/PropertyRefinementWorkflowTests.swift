import Foundation
import Testing
@testable import NumiVivoKit

/// Apple/full-module integration checks. Portable numerical assertions are in
/// Tools/NativeChemistry/PropertyRefinementChecks.swift; they do not replace this suite.
@Suite(.serialized) struct PropertyRefinementWorkflowTests {
    private func implementation() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data("property-refinement-test-implementation-v1".utf8))
    }
    private func json<T: Encodable>(_ value: T) throws -> VivoJSONValue {
        try VivoCanonicalJSON.decode(VivoJSONValue.self,from: VivoCanonicalJSON.encode(value))
    }
    @Test func selectedSolverDispatchIsExplicit() throws {
        let request = VivoSelectedCISolverRequest()
        let decoded = try VivoNativeSolverDispatch.decode(VivoCanonicalJSON.encode(request),implementationFingerprint: implementation())
        #expect(decoded.operation.identifier == "vivo.native.selected-ci")
        let restored = try VivoCanonicalJSON.decode(VivoSelectedCISolverRequest.self,from: VivoCanonicalJSON.encode(decoded.configuration))
        #expect(restored == request)
        guard case .object(var bad) = try json(request) else { Issue.record("request must be an object"); return }
        bad["unrecognizedSelectionRule"] = .boolean(true)
        #expect(throws: VivoChemistryError.self) {
            try VivoNativeSolverDispatch.decode(VivoCanonicalJSON.encode(VivoJSONValue.object(bad)),implementationFingerprint: implementation())
        }
    }
    @Test func boundedProbeCannotBeRelabeledAsAcceptedEigenpair() async throws {
        let h = try VivoPropertyRefinementExamples.algebraicPath().points[0].hamiltonian
        let cfg = VivoSelectedCIConfiguration(maximumDeterminants: 1,selectionBatchSize: 1)
        let probe = VivoSelectedCISolverRequest(configuration: cfg,requireConverged: false)
        let op = try VivoCorrelationRefinementOperations.selectedCI(implementationFingerprint: implementation())
        let inputs = ["hamiltonian": try VivoCanonicalJSON.encode(h)]
        let output = try await op.execute(json(probe),inputs,.init())
        let result = try VivoCanonicalJSON.decode(VivoSelectedCIResult.self,from: #require(output["electronic"]))
        #expect(!result.converged)
        try op.validateOutputs(json(probe),inputs,output,.init())
        let accepted = VivoSelectedCISolverRequest(configuration: cfg,requireConverged: true)
        #expect(throws: VivoChemistryError.self) { try op.validateOutputs(json(accepted),inputs,output,.init()) }
        let changedSeeds = VivoSelectedCISolverRequest(configuration: cfg,initialDeterminants: [4],requireConverged: false)
        #expect(throws: VivoChemistryError.self) { try op.validateOutputs(json(changedSeeds),inputs,output,.init()) }
    }
    @Test func orbitalInformationIsBoundToActualInputState() async throws {
        let op = try VivoCorrelationRefinementOperations.orbitalInformation(implementationFingerprint: implementation())
        let state = VivoCIState(orbitalCount: 2,alphaElectrons: 1,betaElectrons: 0,
                               determinants: [1,4],coefficients: [1/sqrt(2),1/sqrt(2)])
        let cfg = try json(VivoOrbitalInformationSelection(orbitals: [0,1],pairs: [.init(0,1)]))
        let inputs = ["state": try VivoCanonicalJSON.encode(state)]
        let outputs = try await op.execute(cfg,inputs,.init())
        try op.validateOutputs(cfg,inputs,outputs,.init())
        let different = VivoCIState(orbitalCount: 2,alphaElectrons: 1,betaElectrons: 0,determinants: [1],coefficients: [1])
        #expect(throws: VivoChemistryError.self) {
            try op.validateOutputs(cfg,["state": VivoCanonicalJSON.encode(different)],outputs,.init())
        }
    }
    @Test func refinementCachesOnlyRequestBoundReconstructibleEvidence() async throws {
        let request = try VivoPropertyRefinementExamples.algebraicPath()
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-refinement-\(UUID().uuidString)")
        defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root)
        let input = try await store.put(data: VivoCanonicalJSON.encode(request),kind: "vivo.property-directed-space-request",mediaType: "application/json")
        let id = try implementation(), op = VivoCorrelationRefinementOperations.propertyDirectedSpace(implementationFingerprint: id)
        let workflow = VivoChemistryWorkflow(store: store)
        let task = VivoChemistryTask(operation: op.identifier,version: op.version,implementationFingerprint: id,
            inputs: [.init(name: "request",artifact: input.fingerprint,kind: "vivo.property-directed-space-request")],
            configuration: .object([:]),outputs: op.outputs,resources: .init(budget: request.budget))
        let first = try await workflow.run(task,using: op), second = try await workflow.run(task,using: op)
        #expect(!first.reused && second.reused)
        #expect(first.outputs == second.outputs)
        let artifact = try #require(first.outputs.first)
        let data = try await workflow.payload(artifact: artifact.artifact,expectedKind: artifact.kind)
        let report = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceResult.self,from: data)
        #expect(report.sensitivityEstablishedWithinDeclaredPool)
        var changedBudget = request.budget; changedBudget.maximumOperatorApplications += 1
        #expect(throws: VivoChemistryError.self) {
            try op.validateOutputs(.object([:]),["request": VivoCanonicalJSON.encode(request)],["refinement": data],changedBudget)
        }
    }
    @Test func anchorExportCarriesEvidenceAndReproducesChosenEnergy() throws {
        let request = try VivoPropertyRefinementExamples.algebraicPath()
        let result = try VivoPropertyDirectedSpace.run(request)
        let h = try VivoCorrelationRefinementOperations.materializeHamiltonian(result,
            selection: .init(pointIdentifier: "barrier-point"),budget: request.budget)
        #expect(h.orbitalIdentifiers == ["required","reaction-sensitive"])
        #expect(h.provenance["refinementEvidenceSHA256"]?.count == 64)
        #expect(h.provenance["refinementPartitionSHA256"]?.count == 64)
        let ci = try VivoDirectCI.solve(h,configuration: .init(residualTolerance: 1e-12))
        let expected = try #require(result.confirmation?.chosen.points.first { $0.pointIdentifier == "barrier-point" })
        #expect(abs(ci.roots[0].energyHartree-expected.energyHartree) < 1e-10)
        #expect(throws: VivoChemistryError.self) {
            try VivoCorrelationRefinementOperations.materializeHamiltonian(result,
                selection: .init(pointIdentifier: "unseen-geometry"),budget: request.budget)
        }
    }
    @Test func exportedHamiltonianValidatorRejectsTamperedEnergy() async throws {
        let request = try VivoPropertyRefinementExamples.algebraicPath(), report = try VivoPropertyDirectedSpace.run(request)
        let op = try VivoCorrelationRefinementOperations.refinedHamiltonian(implementationFingerprint: implementation())
        let cfg = try json(VivoRefinedHamiltonianSelection(pointIdentifier: "reactant"))
        let inputs = ["refinement": try VivoCanonicalJSON.encode(report)]
        let outputs = try await op.execute(cfg,inputs,request.budget)
        var altered = try JSONSerialization.jsonObject(with: #require(outputs["hamiltonian"])) as! [String:Any]
        altered["constantEnergyHartree"] = 10.0
        let bad = try JSONSerialization.data(withJSONObject: altered)
        #expect(throws: VivoChemistryError.self) {
            try op.validateOutputs(cfg,inputs,["hamiltonian": bad],request.budget)
        }
    }
    @Test func platformRecipePlansTypedStagesAndRejectsWrongVersions() throws {
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: implementation())
        let recipe = try VivoWorkflowTemplates.make("property-refinement")
        let plan = try VivoWorkflowPlanner.compile(recipe,registry: registry)
        #expect(plan.topologicalOrder == ["refine","anchor","selected","state","information"])
        var wrong = recipe; wrong.nodes[2].version = "1"
        #expect(throws: VivoChemistryError.self) { try VivoWorkflowPlanner.compile(wrong,registry: registry) }
        wrong = recipe; wrong.nodes[4].inputs["state"] = .output(node: "selected",port: "electronic")
        #expect(throws: VivoChemistryError.self) { try VivoWorkflowPlanner.compile(wrong,registry: registry) }
    }
    @Test func platformRecipeExecutesAndResumesNativeStages() async throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-refinement-recipe-\(UUID().uuidString)")
        defer { try? FileManager.default.removeItem(at: root) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: implementation())
        let executor = try VivoWorkflowExecutor(store: VivoArtifactStore(rootURL: root),registry: registry)
        let recipe = try VivoWorkflowTemplates.propertyRefinement()
        let first = try await executor.run(recipe)
        #expect(first.report.allTasksSucceeded && first.report.nodes.count == 5)
        let repeated = try await executor.run(recipe)
        #expect(repeated.report.allTasksSucceeded)
        #expect(repeated.report.nodes.allSatisfy { $0.outcome.reused })
        #expect(first.report.exports == repeated.report.exports)
    }
    @Test func unconfirmedRecipeCannotExportOrRunDownstreamSolver() async throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-refinement-blocked-\(UUID().uuidString)")
        defer { try? FileManager.default.removeItem(at: root) }
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: implementation())
        let executor = try VivoWorkflowExecutor(store: VivoArtifactStore(rootURL: root),registry: registry)
        var recipe = try VivoWorkflowTemplates.propertyRefinement()
        let original = try VivoPropertyRefinementExamples.algebraicPath()
        var object = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(original)) as! [String:Any]
        object["maximumPointEvaluations"] = 1
        let request = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceRequest.self,
            from: JSONSerialization.data(withJSONObject: object))
        recipe.artifacts[0].source = .json(kind: "vivo.property-directed-space-request",payload: try json(request))
        let execution = try await executor.run(recipe)
        #expect(!execution.report.allTasksSucceeded)
        let records = Dictionary(uniqueKeysWithValues: execution.report.nodes.map { ($0.identifier,$0) })
        #expect(records["refine"]?.outcome.outputs != nil)
        guard let anchor = records["anchor"], case .failed = anchor.outcome else {
            Issue.record("unconfirmed report must fail at anchor export"); return
        }
        guard let solver = records["selected"], case .blocked = solver.outcome else {
            Issue.record("downstream solver must be blocked, not silently downgraded"); return
        }
    }

}
