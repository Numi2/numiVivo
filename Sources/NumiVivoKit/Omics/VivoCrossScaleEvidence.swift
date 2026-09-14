import Foundation

/// The ordered stages in the long-term variant-to-tissue research path.
/// Stage order is part of the contract: a graph cannot silently skip a
/// biological boundary or feed a downstream result back into an upstream one.
public enum VivoCrossScaleStage: String, Codable, CaseIterable, Sendable {
    case variant
    case regulation
    case rnaCellState
    case protein
    case molecularMechanism
    case kinetics
    case cellularPhenotype
    case tissuePrediction

    public var order: Int { Self.allCases.firstIndex(of: self)! }
    public var next: VivoCrossScaleStage? {
        let i = order + 1
        return i < Self.allCases.count ? Self.allCases[i] : nil
    }
    public var displayName: String {
        switch self {
        case .variant: "variant"
        case .regulation: "regulation"
        case .rnaCellState: "RNA/cell state"
        case .protein: "protein"
        case .molecularMechanism: "molecular mechanism"
        case .kinetics: "reaction/kinetics"
        case .cellularPhenotype: "cellular phenotype"
        case .tissuePrediction: "tissue prediction"
        }
    }
}

/// How a node or link is supported.  These values describe evidence class;
/// they do not turn a prediction into a measurement or a clinical claim.
public enum VivoCrossScaleEvidenceClass: String, Codable, Sendable {
    case measured
    case predicted
    case simulated
    case unavailable
}

/// Admission state for a transition. `qualified` is intentionally strict:
/// it requires an identified model, a validation artifact and at least one
/// held-out observation. A hypothesis can be carried through the graph but
/// cannot authorize an end-to-end outcome report.
public enum VivoCrossScaleLinkStatus: String, Codable, Sendable {
    case qualified
    case hypothesis
    case unavailable
}

public struct VivoCrossScaleNode: Codable, Sendable, Equatable {
    public let id: String
    public let stage: VivoCrossScaleStage
    public let quantity: String
    public let unit: String
    public let contextFingerprint: VivoFingerprint
    public let evidenceClass: VivoCrossScaleEvidenceClass
    public let sourceFingerprint: VivoFingerprint?
    public let modelFingerprint: VivoFingerprint?

    public init(id: String, stage: VivoCrossScaleStage, quantity: String, unit: String,
                contextFingerprint: VivoFingerprint, evidenceClass: VivoCrossScaleEvidenceClass,
                sourceFingerprint: VivoFingerprint?, modelFingerprint: VivoFingerprint?) {
        self.id = id; self.stage = stage; self.quantity = quantity; self.unit = unit
        self.contextFingerprint = contextFingerprint; self.evidenceClass = evidenceClass
        self.sourceFingerprint = sourceFingerprint; self.modelFingerprint = modelFingerprint
    }
}

public struct VivoCrossScaleLink: Codable, Sendable, Equatable {
    public let id: String
    public let fromNodeID: String
    public let toNodeID: String
    public let relation: String
    public let sourceFingerprints: [VivoFingerprint]
    public let modelFingerprint: VivoFingerprint?
    public let validationFingerprint: VivoFingerprint?
    public let heldOutObservations: Int
    public let status: VivoCrossScaleLinkStatus
    public let limitations: [String]

    public init(id: String, fromNodeID: String, toNodeID: String, relation: String,
                sourceFingerprints: [VivoFingerprint], modelFingerprint: VivoFingerprint?,
                validationFingerprint: VivoFingerprint?, heldOutObservations: Int,
                status: VivoCrossScaleLinkStatus, limitations: [String]) {
        self.id = id; self.fromNodeID = fromNodeID; self.toNodeID = toNodeID; self.relation = relation
        self.sourceFingerprints = sourceFingerprints; self.modelFingerprint = modelFingerprint
        self.validationFingerprint = validationFingerprint; self.heldOutObservations = heldOutObservations
        self.status = status; self.limitations = limitations
    }
}

/// A machine-readable record of a missing boundary. Gaps are retained instead
/// of being inferred from an absent file or silently treated as zero evidence.
public struct VivoCrossScaleGap: Codable, Sendable, Equatable {
    public let from: VivoCrossScaleStage
    public let to: VivoCrossScaleStage
    public let reason: String
    public let requiredEvidence: String

    public init(from: VivoCrossScaleStage, to: VivoCrossScaleStage, reason: String, requiredEvidence: String) {
        self.from = from; self.to = to; self.reason = reason; self.requiredEvidence = requiredEvidence
    }
}

public enum VivoCrossScaleReadiness: String, Codable, Sendable {
    case incomplete
    case hypothesisOnly
    case qualifiedResearchPath
}

public struct VivoCrossScaleAssessment: Codable, Sendable, Equatable {
    public let readiness: VivoCrossScaleReadiness
    public let missingBoundaries: [String]
    public let hypothesisBoundaries: [String]
    public let qualifiedBoundaries: [String]
    /// True means the graph satisfies the research-path gate only. It never
    /// means clinical, treatment or patient-specific authorization.
    public let canReportResearchOutcome: Bool
}

public struct VivoCrossScaleEvidenceGraph: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let sourceDescription: String
    public let nodes: [VivoCrossScaleNode]
    public let links: [VivoCrossScaleLink]
    public let gaps: [VivoCrossScaleGap]

    public init(schemaVersion: Int = 1, id: String, sourceDescription: String,
                nodes: [VivoCrossScaleNode], links: [VivoCrossScaleLink], gaps: [VivoCrossScaleGap]) {
        self.schemaVersion = schemaVersion; self.id = id; self.sourceDescription = sourceDescription
        self.nodes = nodes; self.links = links; self.gaps = gaps
    }

    public static let maximumNodes = 1_024
    public static let maximumLinks = 2_048
    public static let maximumGaps = 256

    private func label(_ value: String) -> Bool {
        !value.isEmpty && value.utf8.count <= 256 &&
        !value.contains("\0") && !value.contains("\t") &&
        !value.contains("\n") && !value.contains("\r")
    }

    /// Validates identity, evidence class, adjacent stage order and the strict
    /// fields required before a link can be called qualified.
    public func validate() throws {
        guard schemaVersion == 1, vivoOmicsID(id), !sourceDescription.isEmpty,
              sourceDescription.utf8.count <= 16_384,
              !nodes.isEmpty, nodes.count <= Self.maximumNodes,
              links.count <= Self.maximumLinks, gaps.count <= Self.maximumGaps else {
            throw VivoOmicsError.invalid("cross-scale graph schema, identity or bounds")
        }
        guard Set(nodes.map(\.id)).count == nodes.count,
              Set(links.map(\.id)).count == links.count else {
            throw VivoOmicsError.invalid("cross-scale node/link identifiers must be unique")
        }
        let nodeByID = Dictionary(uniqueKeysWithValues: nodes.map { ($0.id, $0) })
        for node in nodes {
            guard vivoOmicsID(node.id), label(node.quantity), label(node.unit) else {
                throw VivoOmicsError.invalid("cross-scale node identity or quantity")
            }
            switch node.evidenceClass {
            case .measured:
                guard node.sourceFingerprint != nil, node.modelFingerprint == nil else {
                    throw VivoOmicsError.invalid("measured node requires a source and no predictive model")
                }
            case .predicted, .simulated:
                guard node.sourceFingerprint != nil, node.modelFingerprint != nil else {
                    throw VivoOmicsError.invalid("predicted or simulated node requires source and model fingerprints")
                }
            case .unavailable:
                guard node.sourceFingerprint == nil, node.modelFingerprint == nil else {
                    throw VivoOmicsError.invalid("unavailable node cannot carry source or model evidence")
                }
            }
        }
        var stagePairs = Set<String>()
        for link in links {
            guard vivoOmicsID(link.id), vivoOmicsID(link.fromNodeID), vivoOmicsID(link.toNodeID),
                  vivoOmicsID(link.relation), !link.limitations.isEmpty,
                  link.limitations.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 4_096 }),
                  link.heldOutObservations >= 0,
                  link.sourceFingerprints.count <= 64,
                  Set(link.sourceFingerprints).count == link.sourceFingerprints.count else {
                throw VivoOmicsError.invalid("cross-scale link identity, limitations or evidence bounds")
            }
            guard let from = nodeByID[link.fromNodeID], let to = nodeByID[link.toNodeID], from.id != to.id,
                  to.stage.order == from.stage.order + 1 else {
                throw VivoOmicsError.invalid("cross-scale links must connect adjacent ordered stages")
            }
            let pair = "\(from.stage.rawValue)->\(to.stage.rawValue)"
            guard stagePairs.insert(pair).inserted else {
                throw VivoOmicsError.invalid("cross-scale stage boundaries must have one link")
            }
            switch link.status {
            case .qualified:
                guard !link.sourceFingerprints.isEmpty, link.modelFingerprint != nil,
                      link.validationFingerprint != nil, link.heldOutObservations > 0 else {
                    throw VivoOmicsError.invalid("qualified link requires source, model, validation and held-out observations")
                }
            case .hypothesis:
                guard !link.sourceFingerprints.isEmpty else {
                    throw VivoOmicsError.invalid("hypothesis link requires source evidence")
                }
            case .unavailable:
                guard link.sourceFingerprints.isEmpty, link.modelFingerprint == nil,
                      link.validationFingerprint == nil, link.heldOutObservations == 0 else {
                    throw VivoOmicsError.invalid("unavailable link cannot carry predictive or validation evidence")
                }
            }
        }
        var gapPairs = Set<String>()
        for gap in gaps {
            guard let to = gap.from.next, to == gap.to, !gap.reason.isEmpty, !gap.requiredEvidence.isEmpty,
                  gap.reason.utf8.count <= 4_096, gap.requiredEvidence.utf8.count <= 4_096 else {
                throw VivoOmicsError.invalid("cross-scale gap must name one adjacent boundary")
            }
            guard gapPairs.insert("\(gap.from.rawValue)->\(gap.to.rawValue)").inserted else {
                throw VivoOmicsError.invalid("cross-scale gaps must be unique")
            }
        }
    }

    /// Summarizes the seven required adjacent boundaries without treating
    /// missing links as failures of a model or as biological zeros.
    public func assess() throws -> VivoCrossScaleAssessment {
        try validate()
        let nodeByID = Dictionary(uniqueKeysWithValues: nodes.map { ($0.id, $0) })
        var missing: [String] = [], hypotheses: [String] = [], qualified: [String] = []
        for stage in VivoCrossScaleStage.allCases.dropLast() {
            guard let next = stage.next else { continue }
            let label = "\(stage.displayName) → \(next.displayName)"
            if let link = links.first(where: {
                nodeByID[$0.fromNodeID]?.stage == stage && nodeByID[$0.toNodeID]?.stage == next
            }) {
                switch link.status {
                case .qualified: qualified.append(label)
                case .hypothesis: hypotheses.append(label)
                case .unavailable: missing.append(label)
                }
            } else {
                missing.append(label)
            }
        }
        let readiness: VivoCrossScaleReadiness
        if !missing.isEmpty { readiness = .incomplete }
        else if !hypotheses.isEmpty { readiness = .hypothesisOnly }
        else { readiness = .qualifiedResearchPath }
        return .init(readiness: readiness, missingBoundaries: missing,
                     hypothesisBoundaries: hypotheses, qualifiedBoundaries: qualified,
                     canReportResearchOutcome: readiness == .qualifiedResearchPath)
    }

    /// The explicit admission gate for an end-to-end research outcome. A
    /// hypothesis-only or incomplete graph returns a bounded error naming the
    /// first missing boundary; callers must not substitute an outcome claim.
    public func requireQualifiedResearchPath() throws {
        let assessment = try assess()
        guard assessment.canReportResearchOutcome else {
            let boundary = (assessment.missingBoundaries + assessment.hypothesisBoundaries).first ?? "unknown boundary"
            throw VivoOmicsError.invalid("cross-scale research outcome path is not qualified at \(boundary)")
        }
    }
}
