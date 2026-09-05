import Foundation

public struct VivoTransitionStateDerivation: Sendable {
    public let request: VivoTransitionStateRateRequest
    public let estimate: VivoTransitionStateRateEstimate
    public let parameter: VivoKineticParameter
    public let evidenceFingerprint: VivoFingerprint
    public let evidenceData: Data

    private struct Evidence: Codable {
        let schemaVersion: UInt32
        let numericalMethod: String
        let request: VivoTransitionStateRateRequest
        let estimate: VivoTransitionStateRateEstimate
    }
    public static func calculate(_ request: VivoTransitionStateRateRequest) throws -> Self {
        let estimate = try VivoTransitionStateRateEstimator.estimate(request)
        let data = try VivoCanonicalJSON.encode(Evidence(schemaVersion: 1,
            numericalMethod: "numivivo.conditional-classical-tst.v1", request: request, estimate: estimate))
        let fingerprint = try VivoCanonicalJSON.fingerprint(data)
        // Assumed upstream inputs must not become high-confidence evidence just
        // because the deterministic conversion was calculated successfully.
        let origin: VivoKineticOrigin = estimate.containsAssumptions ? .assumed : .calculated
        let parameter = VivoKineticParameter(value: estimate.ratePerSecond, unit: .perSecond,
            origin: origin, uncertainty: .unknown,
            evidence: .init(source: "NumiVivo conditional transition-state derivation",
                            locator: "request + estimate; full kinetic uncertainty is unresolved",
                            sourceFingerprint: fingerprint.hex))
        return .init(request: request, estimate: estimate, parameter: parameter,
                     evidenceFingerprint: fingerprint, evidenceData: data)
    }

    /// Explicitly replace ONLY the bound-complex conversion rate. Association,
    /// dissociation, target turnover and their evidence remain unchanged.
    public func applyingInactivation(to model: VivoCovalentKineticPack) throws -> VivoCovalentKineticPack {
        try model.validate()
        guard model.context == request.barrier.context,
              evidenceFingerprint == (try VivoCanonicalJSON.fingerprint(evidenceData)) else {
            throw VivoKineticsError.invalid("derived-rate context or evidence identity mismatch")
        }
        let result = VivoCovalentKineticPack(identifier: model.identifier, context: model.context,
            association: model.association, dissociation: model.dissociation, inactivation: parameter,
            targetTurnover: model.targetTurnover, baselineTarget: model.baselineTarget,
            competitor: model.competitor, maximumUnboundDrugM: model.maximumUnboundDrugM)
        try result.validate()
        return result
    }
}
