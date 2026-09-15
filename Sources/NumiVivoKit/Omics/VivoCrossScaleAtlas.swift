import Foundation

/// Coverage and provenance for the native Atlas-to-cross-scale projection.
/// The model fingerprint identifies the provider/client identity recorded by
/// the capture; it is not an invented immutable service model version.
public struct VivoCrossScaleAtlasCoverage: Codable, Sendable, Equatable {
    public let requestFingerprint: VivoFingerprint
    public let captureFingerprint: VivoFingerprint
    public let modelFingerprint: VivoFingerprint
    public let requestedVariants: Int
    public let availableVariants: Int
    public let noDataVariants: Int
    public let failedVariants: Int
    public let scorerNames: [String]
    public let scoreCells: Int

    public init(requestFingerprint: VivoFingerprint, captureFingerprint: VivoFingerprint,
                modelFingerprint: VivoFingerprint, requestedVariants: Int,
                availableVariants: Int, noDataVariants: Int, failedVariants: Int,
                scorerNames: [String], scoreCells: Int) {
        self.requestFingerprint = requestFingerprint; self.captureFingerprint = captureFingerprint
        self.modelFingerprint = modelFingerprint; self.requestedVariants = requestedVariants
        self.availableVariants = availableVariants; self.noDataVariants = noDataVariants
        self.failedVariants = failedVariants; self.scorerNames = scorerNames
        self.scoreCells = scoreCells
    }
}

/// A cross-scale graph plus the exact Atlas coverage that produced it.
/// `graph` remains assessable by the normal cross-scale admission gate.
public struct VivoCrossScaleAtlasProjection: Codable, Sendable, Equatable {
    public let graph: VivoCrossScaleEvidenceGraph
    public let coverage: VivoCrossScaleAtlasCoverage

    public init(graph: VivoCrossScaleEvidenceGraph, coverage: VivoCrossScaleAtlasCoverage) {
        self.graph = graph; self.coverage = coverage
    }

    public func validate() throws {
        try graph.validate()
        guard coverage.requestedVariants >= 0,
              coverage.availableVariants >= 0,
              coverage.noDataVariants >= 0,
              coverage.failedVariants >= 0,
              coverage.availableVariants + coverage.noDataVariants + coverage.failedVariants == coverage.requestedVariants,
              coverage.scoreCells >= 0,
              coverage.scorerNames.count <= VivoAtlasEvidence.maximumVariants * 16,
              Set(coverage.scorerNames).count == coverage.scorerNames.count,
              coverage.scorerNames.allSatisfy(vivoOmicsID) else {
            throw VivoOmicsError.invalid("Atlas cross-scale coverage")
        }
    }
}

/// Converts a verified Atlas capture into the first edge of the ordered
/// cross-scale graph. Available Atlas scores become a hypothesis link from a
/// measured variant node to a predicted regulation node. Every later edge is
/// explicitly unavailable until its own model, validation artifact and
/// held-out observations are supplied.
public enum VivoCrossScaleAtlasBuilder {
    public static let implementation = "numivivo.org/cross-scale/atlas-hypothesis/v1"

    public static func build(request: VivoAtlasRequest, capture: VivoAtlasCapture,
                             requestData: Data) throws -> VivoCrossScaleAtlasProjection {
        try VivoAtlasEvidence.validate(capture: capture, request: request, requestData: requestData)
        let requestFingerprint = try VivoCanonicalJSON.fingerprint(requestData)
        let captureData = try VivoCanonicalJSON.encode(capture)
        let captureFingerprint = try VivoCanonicalJSON.fingerprint(captureData)
        let modelIdentity = "\(capture.provider)|\(capture.clientCommit)|\(capture.serviceVersionStatus)"
        let modelFingerprint = try VivoCanonicalJSON.fingerprint(Data(modelIdentity.utf8))
        let contextFingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request.binding))

        let statusByVariant = Dictionary(uniqueKeysWithValues: capture.outcomes.map { ($0.variantKey, $0.status) })
        let available = capture.outcomes.filter { $0.status == .available }.count
        let noData = capture.outcomes.filter { $0.status == .noData }.count
        let failed = capture.outcomes.filter { $0.status == .failed }.count
        var scoreCells = 0
        for table in capture.tables {
            let product = table.observations.count.multipliedReportingOverflow(by: table.tracks.count)
            guard !product.overflow, product.partialValue <= Int.max - scoreCells else {
                throw VivoOmicsError.limit("Atlas score-cell coverage")
            }
            scoreCells += product.partialValue
        }
        let coverage = VivoCrossScaleAtlasCoverage(
            requestFingerprint: requestFingerprint, captureFingerprint: captureFingerprint,
            modelFingerprint: modelFingerprint, requestedVariants: capture.outcomes.count,
            availableVariants: available, noDataVariants: noData, failedVariants: failed,
            scorerNames: Array(Set(capture.tables.map(\.scorer))).sorted(), scoreCells: scoreCells)

        let graphID = "atlas-cross-scale-\(String(requestFingerprint.hex.prefix(16)))"
        let stages = VivoCrossScaleStage.allCases
        let nodes = stages.enumerated().map { index, stage in
            let hasPrediction = index == 1 && available > 0
            let evidenceClass: VivoCrossScaleEvidenceClass = index == 0 ? .measured : (hasPrediction ? .predicted : .unavailable)
            let source: VivoFingerprint? = index == 0 ? requestFingerprint : (hasPrediction ? captureFingerprint : nil)
            let model: VivoFingerprint? = hasPrediction ? modelFingerprint : nil
            let quantity: String
            let unit: String
            switch stage {
            case .variant:
                quantity = "declared genomic variant"; unit = "allele"
            case .regulation:
                quantity = "Atlas molecular-effect score"; unit = "provider-score"
            default:
                quantity = "unavailable downstream quantity"; unit = "unavailable"
            }
            return VivoCrossScaleNode(id: "\(graphID)-\(stage.rawValue)", stage: stage,
                quantity: quantity, unit: unit, contextFingerprint: contextFingerprint,
                evidenceClass: evidenceClass, sourceFingerprint: source, modelFingerprint: model)
        }

        var links: [VivoCrossScaleLink] = []
        for (index, stage) in stages.dropLast().enumerated() {
            guard let next = stage.next else { continue }
            let firstEdge = index == 0
            let atlasAvailable = firstEdge && available > 0
            let status: VivoCrossScaleLinkStatus = atlasAvailable ? .hypothesis : .unavailable
            let sourceFingerprints = atlasAvailable ? [requestFingerprint, captureFingerprint] : []
            let modelFingerprintForLink = atlasAvailable ? modelFingerprint : nil
            let limitations: [String]
            if atlasAvailable {
                limitations = VivoAtlasEvidence.limitations + [
                    "This edge is an external molecular-effect hypothesis; no allele-resolved held-out regulatory measurement is attached.",
                    "The Atlas capture does not establish the downstream RNA/cell-state boundary."
                ]
            } else if firstEdge {
                limitations = ["The Atlas capture has no available variant score; no prediction is substituted for noData or failed queries."]
            } else {
                limitations = ["No independently validated model, validation artifact and held-out observations are supplied for this boundary."]
            }
            links.append(VivoCrossScaleLink(id: "\(graphID)-\(stage.rawValue)-to-\(next.rawValue)",
                fromNodeID: nodes[index].id, toNodeID: nodes[index + 1].id,
                relation: firstEdge ? "atlas-molecular-effect-hypothesis" : "downstream-boundary",
                sourceFingerprints: sourceFingerprints, modelFingerprint: modelFingerprintForLink,
                validationFingerprint: nil, heldOutObservations: 0, status: status,
                limitations: limitations))
        }
        let graph = VivoCrossScaleEvidenceGraph(id: graphID,
            sourceDescription: "Verified AlphaGenome Atlas capture projected into an explicit cross-scale hypothesis graph",
            nodes: nodes, links: links, gaps: [])
        let result = VivoCrossScaleAtlasProjection(graph: graph, coverage: coverage)
        try result.validate()
        // Keep the status lookup in the construction path so a future change
        // cannot accidentally treat an unrequested variant as available.
        guard statusByVariant.count == capture.outcomes.count else {
            throw VivoOmicsError.invalid("Atlas outcome identity")
        }
        return result
    }
}
