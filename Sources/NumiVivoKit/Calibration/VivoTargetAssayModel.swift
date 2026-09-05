import Foundation

/// One observation model per case. An independent case boundary is not an
/// automatic assertion that experiments share no donors or batch effects.
public struct VivoTargetAssayModel: Codable, Equatable, Sendable {
    public let caseIdentifier: String
    public let noise: VivoGaussianAssayNoise
    public let support: [String: VivoAssaySupport]
    public let evidence: VivoKineticEvidence
    public init(caseIdentifier: String, noise: VivoGaussianAssayNoise = .init(),
                support: [String: VivoAssaySupport] = [:], evidence: VivoKineticEvidence) {
        self.caseIdentifier = caseIdentifier; self.noise = noise; self.support = support; self.evidence = evidence
    }
    public func validate(for item: VivoTargetEngagementStudyCase) throws {
        try noise.validate(); try evidence.validate(origin: .assumed)
        let ids = Set(item.observations.map(\.identifier))
        guard caseIdentifier == item.identifier, support.keys.allSatisfy(ids.contains) else {
            throw VivoPosteriorError.invalid("assay support refers to a foreign case or observation")
        }
        for value in support.values { try value.validate() }
    }
}

extension VivoKineticFitField {
    public var isAssayParameter: Bool {
        switch self {
        case .assayNoiseScale, .assayNoiseFloor, .assayBias, .assayCorrelationFraction, .assayCorrelationTime: return true
        default: return false
        }
    }
}

extension VivoTargetPosteriorProblem {
    public var usesExtendedAssays: Bool { !(assays ?? []).isEmpty || bindings.contains { $0.field.isAssayParameter } }
    public var usesExtendedTrainingAssays: Bool {
        let ids = Set(study.cases.filter { $0.partition == .calibration }.map(\.identifier))
        return (assays ?? []).contains { ids.contains($0.caseIdentifier) }
            || bindings.contains { $0.field.isAssayParameter && $0.caseIdentifiers.contains(where: ids.contains) }
    }
    public func assaySupport(caseIdentifier: String, observationIdentifier: String) -> VivoAssaySupport {
        assays?.first { $0.caseIdentifier == caseIdentifier }?.support[observationIdentifier] ?? .exact
    }
    public func assayNoise(caseIdentifier: String, values: [Double]) throws -> VivoGaussianAssayNoise {
        guard values.count == bindings.count, study.cases.contains(where: { $0.identifier == caseIdentifier }) else {
            throw VivoPosteriorError.invalid("assay parameter layout or case identity")
        }
        let base = assays?.first { $0.caseIdentifier == caseIdentifier }?.noise ?? .init()
        var scale = base.scale, floor = base.additionalStandardDeviation, bias = base.bias
        var rho = base.correlatedFraction, tau = base.correlationTimeSeconds
        for (binding, value) in zip(bindings, values) {
            guard value.isFinite, value >= binding.parameter.lower, value <= binding.parameter.upper else {
                throw VivoPosteriorError.invalid("assay candidate outside its prior")
            }
            guard binding.caseIdentifiers.contains(caseIdentifier) else { continue }
            switch binding.field {
            case .assayNoiseScale: scale = value
            case .assayNoiseFloor: floor = value
            case .assayBias: bias = value
            case .assayCorrelationFraction: rho = value
            case .assayCorrelationTime: tau = value
            default: break
            }
        }
        let result = VivoGaussianAssayNoise(scale: scale, additionalStandardDeviation: floor, bias: bias,
                                           correlatedFraction: rho, correlationTimeSeconds: tau)
        try result.validate()
        return result
    }
    internal func validateAssays() throws {
        let byCase = Dictionary(uniqueKeysWithValues: study.cases.map { ($0.identifier, $0) })
        guard (assays ?? []).count <= study.cases.count,
              Set((assays ?? []).map(\.caseIdentifier)).count == (assays ?? []).count else {
            throw VivoPosteriorError.invalid("duplicate assay models or assay capacity")
        }
        for assay in assays ?? [] {
            guard let item = byCase[assay.caseIdentifier] else { throw VivoPosteriorError.invalid("unresolved assay case") }
            try assay.validate(for: item)
        }
        let lower = bindings.map { $0.parameter.lower }, upper = bindings.map { $0.parameter.upper }
        for item in study.cases {
            let low = try assayNoise(caseIdentifier: item.identifier, values: lower)
            let high = try assayNoise(caseIdentifier: item.identifier, values: upper)
            if high.correlatedFraction > 0 {
                guard item.observations.count <= 256,
                      item.observations.allSatisfy({ assaySupport(caseIdentifier: item.identifier, observationIdentifier: $0.identifier) == .exact }) else {
                    throw VivoPosteriorError.invalid("correlated censoring or oversized correlated assay is unsupported")
                }
            }
            if item.partition == .calibration {
                for observation in item.observations { _ = try low.standardDeviation(reported: observation.standardDeviation) }
            }
        }
    }
}
