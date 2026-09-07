import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MolecularMultistateRefinementWorkflowTests {
    @Test func molecularRecipeUsesVersionedNativeStages() throws {
        let id = try VivoCanonicalJSON.fingerprint(Data("molecular-refinement-v2-tests".utf8))
        let registry = try VivoPlatformOperations.registry(implementationFingerprint: id)
        let recipe = try VivoWorkflowTemplates.make("molecular-property-refinement")
        let plan = try VivoWorkflowPlanner.compile(recipe,registry: registry)
        #expect(plan.topologicalOrder == ["prepare","refine","anchor"])
        #expect(recipe.nodes[1].version == "2" && recipe.nodes[2].version == "2")
    }
    @Test func molecularPrimitivesReuseExistingAOAndHFIdentities() throws {
        let request = try VivoPropertyRefinementExamples.molecularHydrogenStretch()
        let id = try VivoCanonicalJSON.fingerprint(Data("molecular-refinement-v2-tests".utf8))
        let plan = try VivoMolecularSpacePreparationWorkflow.plan(request,requestArtifact: id,
            systemArtifacts: Array(repeating: id,count: request.snapshots.count),basisArtifact: id,implementationFingerprint: id)
        #expect(plan.nodes.count == 2*request.snapshots.count+1)
        #expect(plan.nodes[0].operation.identifier == "vivo.native.ao-integrals")
        #expect(plan.nodes[1].operation.identifier == "vivo.native.hf-from-integrals")
        #expect(plan.resultOutput == "request")
    }
    @Test func molecularArtifactRejectsForgedSelection() async throws {
        let request = try VivoPropertyRefinementExamples.molecularHydrogenStretch()
        let id = try VivoCanonicalJSON.fingerprint(Data("molecular-refinement-v2-tests".utf8))
        let op = VivoMolecularSpacePreparationWorkflow.operation(implementationFingerprint: id,suppliedPrimitives: false)
        let inputs = ["request":try VivoCanonicalJSON.encode(request)]
        let outputs = try await op.execute(.object([:]),inputs,request.budget)
        try op.validateOutputs(.object([:]),inputs,outputs,request.budget)
        var bad = outputs
        var json = try JSONSerialization.jsonObject(with: #require(outputs["request"])) as! [String:Any]
        json["mandatoryActiveOrbitals"] = [0]
        bad["request"] = try JSONSerialization.data(withJSONObject: json)
        #expect(throws: VivoChemistryError.self) { try op.validateOutputs(.object([:]),inputs,bad,request.budget) }
    }
    @Test func optimizedAnchorPreservesEveryRequestedRootEnergy() throws {
        let original = try VivoPropertyRefinementExamples.multistateCancellation()
        var json = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(original)) as! [String:Any]
        var cfg = VivoMultiStateCASSCFConfiguration(weights: [0.5,0.5],rootLabels: ["lower","upper"],
            optimization: .init(gradientTolerance: 1e-8,energyToleranceHartree: 1e-12),followRoots: false)
        cfg.davidson.residualTolerance = 1e-12
        let solver = VivoSpaceRefinementSolver.stateAveragedCASSCF(configuration: cfg,
            states: .init(labels: ["lower","upper"],groups: [[0],[1]]))
        json["solver"] = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(solver))
        let request = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceRequest.self,from: JSONSerialization.data(withJSONObject: json))
        let result = try VivoPropertyDirectedSpace.run(request)
        let h = try VivoCorrelationRefinementOperations.materializeHamiltonian(result,
            selection: .init(pointIdentifier: "barrier-point"),budget: request.budget)
        let roots = try VivoDirectCI.solve(h,configuration: .init(roots: 2,residualTolerance: 1e-12),budget: request.budget)
        let point = try #require(result.confirmation?.chosen.points.first { $0.pointIdentifier == "barrier-point" })
        let expected = try #require(point.roots)
        for i in expected.indices { #expect(abs(roots.roots[i].energyHartree-expected[i].energyHartree) < 1e-9) }
        #expect(h.provenance["refinementOrbitalFrameSHA256"] != nil)
    }
}
