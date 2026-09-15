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

        let duplicateGap = VivoCrossScaleGap(from: .variant, to: .regulation,
                                              reason: "contradictory fixture gap",
                                              requiredEvidence: "held-out regulatory observation")
        #expect(throws: (any Error).self) {
            try VivoCrossScaleEvidenceGraph(id: "bad-gap", sourceDescription: "Synthetic structural fixture",
                                            nodes: baseNodes, links: [try links(status: .qualified)[0]],
                                            gaps: [duplicateGap]).validate()
        }
    }

    @Test func atlasCaptureProjectsHypothesisAndRetainsDownstreamGaps() throws {
        let digest = String(repeating: "a", count: 64)
        let binding = VivoGenomicCaseBinding(parentReportSHA256: digest, caseFingerprint: digest,
            dataClass: "synthetic", assembly: "GRCh38", referenceSHA256: digest,
            hlaAlleles: ["HLA-A*02:01"], tumorRNASampleID: "rna")
        let candidateID = try VivoAtlasEvidence.digest(Data("candidate".utf8))
        let input = VivoAtlasVariantInput(candidateID: candidateID, chromosome: "chr1",
            start0: 9, stop0: 10, reference: "A", alternate: "C")
        let request = try VivoAtlasEvidence.plan(binding: binding, inputs: [input],
            requestedScorers: ["SYNTHETIC_RNA"], ontologyTerms: ["CL:0000000"])
        let requestData = try VivoCanonicalJSON.encode(request)
        let query = try #require(request.variants.first)
        let capture = VivoAtlasCapture(schema: "numivivo.org/atlas-capture/v1",
            requestSHA256: try VivoAtlasEvidence.digest(requestData), mode: "fixture",
            provider: "google-deepmind/alphagenome-atlas", clientCommit: VivoAtlasEvidence.clientCommit,
            clientVersion: "0.0.0-synthetic", adapterSHA256: digest,
            retrievedAt: "2026-09-08T00:00:00Z", serviceVersion: nil,
            serviceVersionStatus: "not-exposed-by-pinned-api",
            termsURI: "https://deepmind.google.com/science/alphagenome/terms",
            permittedUse: "noncommercialResearch", trainingPermitted: false,
            referenceCheck: "synthetic-not-biological", referenceSHA256: digest,
            outcomes: [.init(variantKey: query.key, status: .available, diagnostic: nil)],
            tables: [.init(scorer: "SYNTHETIC_RNA", isSigned: true,
                observations: [.init(variantKey: query.key, metadata: ["gene": .string("G1")])],
                tracks: [["ontology": .string("CL:0000000")]], rawScores: [[0.4]], quantiles: [[0.7]])])

        let projection = try VivoCrossScaleAtlasBuilder.build(request: request,
            capture: capture, requestData: requestData)
        try projection.validate()
        #expect(projection.coverage.requestedVariants == 1)
        #expect(projection.coverage.availableVariants == 1)
        #expect(projection.coverage.scoreCells == 1)
        #expect(projection.graph.links.first?.status == .hypothesis)
        #expect(projection.graph.nodes[1].evidenceClass == .predicted)
        let assessment = try projection.graph.assess()
        #expect(assessment.hypothesisBoundaries == ["variant → regulation"])
        #expect(assessment.missingBoundaries.count == 6)
        #expect(!assessment.canReportResearchOutcome)
        #expect(try VivoCanonicalJSON.decode(VivoCrossScaleAtlasProjection.self,
            from: VivoCanonicalJSON.encode(projection)) == projection)

        var unavailable = capture
        unavailable.outcomes[0].status = .noData
        unavailable.tables = []
        let missing = try VivoCrossScaleAtlasBuilder.build(request: request,
            capture: unavailable, requestData: requestData)
        #expect(missing.coverage.availableVariants == 0)
        #expect(missing.coverage.noDataVariants == 1)
        #expect(missing.graph.links.first?.status == .unavailable)
        #expect(try missing.graph.assess().readiness == .incomplete)
    }
}
