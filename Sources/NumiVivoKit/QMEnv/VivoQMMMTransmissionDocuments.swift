import Foundation

/// Serializable analysis envelope for already-computed paired dividing-surface
/// shooting evidence. The executable Born–Oppenheimer provider remains a runtime
/// resource and is deliberately not serialized into workflow artifacts.
public struct VivoQMMMDynamicalTransmissionAnalysisRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-dynamical-transmission-analysis/v2"
    public var schema:String
    public var transmissionRequest:VivoQMMMDynamicalTransmissionRequest
    public var retainedSystemFingerprint:VivoFingerprint
    public var providerFingerprint:VivoFingerprint
    public var pairs:[VivoQMMMTransmissionPair]

    public init(transmissionRequest:VivoQMMMDynamicalTransmissionRequest,
                retainedSystemFingerprint:VivoFingerprint,providerFingerprint:VivoFingerprint,
                pairs:[VivoQMMMTransmissionPair]) {
        schema=Self.schema;self.transmissionRequest=transmissionRequest
        self.retainedSystemFingerprint=retainedSystemFingerprint;self.providerFingerprint=providerFingerprint
        self.pairs=pairs
    }

    public func calculate() throws -> VivoQMMMDynamicalTransmissionResult {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("QM/MM transmission-analysis schema") }
        return try VivoQMMMDynamicalTransmission.analyze(request:transmissionRequest,
            systemFingerprint:retainedSystemFingerprint,providerFingerprint:providerFingerprint,pairs:pairs)
    }
}

/// Deterministic artifact-level application of a validated computed transmission
/// coefficient to the exact PMF rate request it was generated for.
public struct VivoQMMMComputedTransmissionApplicationRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/qmmm-computed-transmission-application/v2"
    public var schema:String
    public var analysisRequest:VivoQMMMDynamicalTransmissionAnalysisRequest
    public var result:VivoQMMMDynamicalTransmissionResult
    public var rateRequest:VivoQMMMFreeEnergyRateRequest

    public init(analysisRequest:VivoQMMMDynamicalTransmissionAnalysisRequest,
                result:VivoQMMMDynamicalTransmissionResult,rateRequest:VivoQMMMFreeEnergyRateRequest) {
        schema=Self.schema;self.analysisRequest=analysisRequest;self.result=result;self.rateRequest=rateRequest
    }

    public func calculate() throws -> VivoQMMMFreeEnergyRateRequest {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("QM/MM computed-transmission application schema") }
        let rebuilt=try analysisRequest.calculate()
        guard rebuilt==result else { throw VivoChemistryError.invalid("computed transmission result does not reconstruct from paired trajectory evidence") }
        return try VivoQMMMDynamicalTransmission.applying(result,
            transmissionRequest:analysisRequest.transmissionRequest,to:rateRequest)
    }
}
