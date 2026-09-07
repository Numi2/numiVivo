import Foundation

/// One independently converged QM/MM variant compared against a converged
/// replicated baseline. The dimension-specific invariants are delegated to the
/// existing chemical-qualification authority; this type only exposes one variant
/// at a time so an adaptive campaign can rank electronic-model, QM-region,
/// sampling, coordinate, transmission or finite-size work independently.
public struct VivoQMMMVariantSensitivityRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-variant-sensitivity/v1"
    public var schema: String
    public var baselineRequest: VivoQMMMReplicatedFreeEnergyRateRequest
    public var baselineResult: VivoQMMMReplicatedFreeEnergyRateResult
    public var variant: VivoQMMMQualificationVariant
    /// Used only to reconstruct the existing dimension-specific qualification
    /// contract. Adaptive acceptance remains a separate criterion.
    public var reconstructionToleranceAbsoluteLogRate: Double
    /// Conservative statistical guard. No independence assumption is made: the
    /// two standard deviations are added before multiplying by this factor.
    public var standardDeviationMultiplier: Double
    public init(baselineRequest: VivoQMMMReplicatedFreeEnergyRateRequest,
                baselineResult: VivoQMMMReplicatedFreeEnergyRateResult,
                variant: VivoQMMMQualificationVariant,
                reconstructionToleranceAbsoluteLogRate: Double = 1e12,
                standardDeviationMultiplier: Double = 2) {
        schema = Self.schema; self.baselineRequest = baselineRequest; self.baselineResult = baselineResult
        self.variant = variant; self.reconstructionToleranceAbsoluteLogRate = reconstructionToleranceAbsoluteLogRate
        self.standardDeviationMultiplier = standardDeviationMultiplier
    }
    public func validate() throws {
        guard schema == Self.schema,
              reconstructionToleranceAbsoluteLogRate.isFinite, reconstructionToleranceAbsoluteLogRate > 0,
              standardDeviationMultiplier.isFinite, standardDeviationMultiplier >= 0, standardDeviationMultiplier <= 10 else {
            throw VivoChemistryError.invalid("QM/MM adaptive variant sensitivity configuration")
        }
        try VivoQMMMReplicatedFreeEnergyRate.validate(baselineResult,request: baselineRequest)
        guard baselineResult.converged else { throw VivoChemistryError.invalid("QM/MM adaptive baseline is not converged") }
        try variant.validate()
    }
}

public struct VivoQMMMVariantSensitivityResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-variant-sensitivity-result/v1"
    public let schema: String
    public let dimension: VivoQMMMQualificationDimension
    public let baselineMeanLogRatePerSecond: Double
    public let variantMeanLogRatePerSecond: Double
    public let absoluteLogRateShift: Double
    public let statisticallyGuardedAbsoluteLogRateShift: Double?
    public let baselineConditionalLogRateStandardDeviation: Double?
    public let variantConditionalLogRateStandardDeviation: Double?
    public let qualificationEvidenceFingerprint: VivoFingerprint
    public let interpretation: String
}

public enum VivoQMMMVariantSensitivity {
    public static let interpretation = "One dimension-specific, independently replicated QM/MM rate variant reconstructed through the existing chemical-qualification invariants. The reported guarded shift is |delta ln k| + multiplier*(sigma_baseline + sigma_variant), a conservative standard-deviation inequality rather than a confidence interval, model-error bound or assumption of independent errors. The calculation does not transfer transmission evidence or equilibrium samples between different Hamiltonians."

    public static func calculate(_ request: VivoQMMMVariantSensitivityRequest) throws -> VivoQMMMVariantSensitivityResult {
        try request.validate()
        let criterion = VivoQMMMSensitivityCriterion(dimension: request.variant.dimension,requiredVariants: 1,
            maximumAbsoluteLogRateShift: request.reconstructionToleranceAbsoluteLogRate)
        let qualificationRequest = VivoQMMMChemicalQualificationRequest(identifier: "adaptive-single-variant:"+request.variant.identifier,
            baselineRequest: request.baselineRequest,baselineResult: request.baselineResult,variants: [request.variant],criteria: [criterion],
            externalTargets: [],maximumDirectValidationAbsoluteLogError: 1e12)
        let qualification = try VivoQMMMChemicalQualification.calculate(qualificationRequest)
        guard qualification.sensitivity.count == 1,
              qualification.sensitivity[0].identifier == request.variant.identifier,
              qualification.sensitivity[0].dimension == request.variant.dimension else {
            throw VivoChemistryError.invalid("QM/MM adaptive variant qualification did not reconstruct uniquely")
        }
        let shift = abs(request.variant.result.meanLogRatePerSecond-request.baselineResult.meanLogRatePerSecond)
        let a = request.baselineResult.combinedConditionalLogRateStandardDeviation
        let b = request.variant.result.combinedConditionalLogRateStandardDeviation
        let guarded: Double? = (a != nil && b != nil) ? shift+request.standardDeviationMultiplier*(a!+b!) : nil
        guard shift.isFinite, guarded == nil || guarded!.isFinite else {
            throw VivoChemistryError.convergence("QM/MM adaptive variant metric overflow")
        }
        return .init(schema: VivoQMMMVariantSensitivityResult.schema,dimension: request.variant.dimension,
            baselineMeanLogRatePerSecond: request.baselineResult.meanLogRatePerSecond,
            variantMeanLogRatePerSecond: request.variant.result.meanLogRatePerSecond,absoluteLogRateShift: shift,
            statisticallyGuardedAbsoluteLogRateShift: guarded,
            baselineConditionalLogRateStandardDeviation: a,variantConditionalLogRateStandardDeviation: b,
            qualificationEvidenceFingerprint: qualification.evidenceFingerprint,interpretation: interpretation)
    }

    public static func validate(_ result: VivoQMMMVariantSensitivityResult,
                                request: VivoQMMMVariantSensitivityRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoChemistryError.invalid("QM/MM adaptive variant sensitivity does not reconstruct")
        }
    }
}
