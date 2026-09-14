import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct VivoCrossScaleEvidenceTests {
    private func fingerprint(_ value: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(value.utf8))
    }

    private func nodes() throws -> [VivoCrossScaleNode] {
        let source = try fingerprint("source")
        let model = try fingerprint("model")
        return VivoCrossScaleStage.allCases.enumerated().map { index, stage in
            .init(id: "node-\(index)", stage: stage, quantity: stage.rawValue,
                  unit: stage == .variant ? "allele" : "research-unit",
                  contextFingerprint: source,
                  evidenceClass: index == 0 ? .measured : .predicted,
                  sourceFingerprint: source, modelFingerprint: index == 0 ? nil : model)
        }
    }

    private func links(status: VivoCrossScaleLinkStatus) throws -> [VivoCrossScaleLink] {
        let source = try fingerprint("link-source")
        let model = try fingerprint("link-model")
        let validation = try fingerprint("link-validation")
        return VivoCrossScaleStage.allCases.dropLast().enumerated().map { index, stage in
            .init(id: "link-\(index)", fromNodeID: "node-\(index)", toNodeID: "node-\(index + 1)",
                  relation: "maps-\(stage.rawValue)",
                  sourceFingerprints: status == .unavailable ? [] : [source],
                  modelFingerprint: status == .unavailable ? nil : model,
                  validationFingerprint: status == .qualified ? validation : nil,
                  heldOutObservations: status == .qualified ? 3 : 0,
                  status: status, limitations: ["research fixture only"])
        }
    }

    private func graph(status: VivoCrossScaleLinkStatus) throws -> VivoCrossScaleEvidenceGraph {
        .init(id: "cross-scale-fixture", sourceDescription: "Synthetic structural fixture",
              nodes: try nodes(), links: try links(status: status), gaps: [])
    }

    @Test func qualifiedPathRequiresEveryAdjacentHeldOutLink() throws {
        let qualified = try graph(status: .qualified)
        let assessment = try qualified.assess()
        #expect(assessment.readiness == .qualifiedResearchPath)
        #expect(assessment.missingBoundaries.isEmpty)
        #expect(assessment.hypothesisBoundaries.isEmpty)
        #expect(assessment.qualifiedBoundaries.count == 7)
        #expect(assessment.canReportResearchOutcome)
        try qualified.requireQualifiedResearchPath()
        #expect(try VivoCanonicalJSON.decode(VivoCrossScaleEvidenceGraph.self,
                                              from: VivoCanonicalJSON.encode(qualified)) == qualified)
    }

    @Test func hypothesisAndMissingLinksCannotAuthorizeOutcome() throws {
        let hypothesis = try graph(status: .hypothesis)
        let hypothesisAssessment = try hypothesis.assess()
        #expect(hypothesisAssessment.readiness == .hypothesisOnly)
        #expect(hypothesisAssessment.hypothesisBoundaries.count == 7)
        #expect(!hypothesisAssessment.canReportResearchOutcome)
        #expect(throws: (any Error).self) { try hypothesis.requireQualifiedResearchPath() }

        let incomplete = VivoCrossScaleEvidenceGraph(id: "incomplete", sourceDescription: "Synthetic structural fixture",
                                                     nodes: try nodes(), links: [], gaps: [])
        let incompleteAssessment = try incomplete.assess()
        #expect(incompleteAssessment.readiness == .incomplete)
        #expect(incompleteAssessment.missingBoundaries.count == 7)
        #expect(!incompleteAssessment.canReportResearchOutcome)
    }

    @Test func invalidStageOrderAndQualifiedEvidenceAreRejected() throws {
        let source = try fingerprint("source")
        let model = try fingerprint("model")
        let validation = try fingerprint("validation")
        let baseNodes = try nodes()
        let skipped = VivoCrossScaleLink(id: "bad-skip", fromNodeID: "node-0", toNodeID: "node-2",
                                         relation: "skips", sourceFingerprints: [source], modelFingerprint: model,
                                         validationFingerprint: validation, heldOutObservations: 1,
                                         status: .qualified, limitations: ["fixture"])
        #expect(throws: (any Error).self) {
            try VivoCrossScaleEvidenceGraph(id: "bad-order", sourceDescription: "Synthetic structural fixture",
                                            nodes: baseNodes, links: [skipped], gaps: []).validate()
        }

        let missingReceipt = VivoCrossScaleLink(id: "bad-receipt", fromNodeID: "node-0", toNodeID: "node-1",
                                                relation: "maps", sourceFingerprints: [source], modelFingerprint: model,
                                                validationFingerprint: nil, heldOutObservations: 0,
                                                status: .qualified, limitations: ["fixture"])
        #expect(throws: (any Error).self) {
            try VivoCrossScaleEvidenceGraph(id: "bad-evidence", sourceDescription: "Synthetic structural fixture",
                                            nodes: baseNodes, links: [missingReceipt], gaps: []).validate()
        }
    }
}
