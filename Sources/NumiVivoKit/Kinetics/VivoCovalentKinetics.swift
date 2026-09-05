import Foundation

public enum VivoKineticsError: Error, LocalizedError, Sendable {
    case invalid(String)
    case unsupported(String)
    case numerical(String)
    case capacity(String)
    public var errorDescription: String? {
        switch self {
        case .invalid(let message): return "Invalid kinetic model: \(message)"
        case .unsupported(let message): return "Unsupported kinetic model: \(message)"
        case .numerical(let message): return "Kinetic numerical failure: \(message)"
        case .capacity(let message): return "Kinetic capacity exceeded: \(message)"
        }
    }
}

public enum VivoKineticOrigin: String, Codable, Sendable {
    case measured, fitted, calculated, assumed
}

/// Canonical units shared with the ProgramPack unit registry. Values are never
/// interpreted as nanomolar, minutes or hours by convention or magnitude.
public enum VivoKineticUnit: String, Codable, Sendable {
    case perSecond = "1/s"
    case perMolarSecond = "M^-1 s^-1"
    case molar = "M"
}

public struct VivoKineticEvidence: Codable, Equatable, Sendable {
    public let source: String
    public let locator: String
    public let sourceFingerprint: String?
    public init(source: String, locator: String, sourceFingerprint: String? = nil) {
        self.source = source; self.locator = locator; self.sourceFingerprint = sourceFingerprint
    }
    public func validate(origin: VivoKineticOrigin) throws {
        guard !source.isEmpty, !locator.isEmpty, source.utf8.count <= 4096, locator.utf8.count <= 4096 else {
            throw VivoKineticsError.invalid("evidence requires a bounded source and exact locator")
        }
        if let sourceFingerprint {
            guard Self.isSHA256(sourceFingerprint) else { throw VivoKineticsError.invalid("evidence digest is not SHA-256") }
        } else if origin != .assumed {
            throw VivoKineticsError.invalid("measured/fitted/calculated values require an immutable evidence snapshot")
        }
    }
    public static func isSHA256(_ value: String) -> Bool {
        value.utf8.count == 64 && value.utf8.allSatisfy { (48...57).contains($0) || (97...102).contains($0) }
    }
}

/// Bounded ranges are sensitivity assumptions, not confidence intervals. A
/// log-normal spread is a declared marginal distribution, not a fitted posterior.
public enum VivoKineticUncertainty: Codable, Equatable, Sendable {
    case unknown
    case bounded(lower: Double, upper: Double)
    case logNormal(logStandardDeviation: Double)

    public func validate(value: Double) throws {
        switch self {
        case .unknown: break
        case .bounded(let lower, let upper):
            guard lower.isFinite, upper.isFinite, lower >= 0, lower <= value, value <= upper else {
                throw VivoKineticsError.invalid("parameter lies outside its declared uncertainty bounds")
            }
        case .logNormal(let deviation):
            guard value > 0, deviation.isFinite, deviation > 0 else {
                throw VivoKineticsError.invalid("log-normal uncertainty requires a positive median and spread")
            }
        }
    }
}

public struct VivoKineticParameter: Codable, Equatable, Sendable {
    public let value: Double
    public let unit: VivoKineticUnit
    public let origin: VivoKineticOrigin
    public let uncertainty: VivoKineticUncertainty
    public let evidence: VivoKineticEvidence
    public init(value: Double, unit: VivoKineticUnit, origin: VivoKineticOrigin,
                uncertainty: VivoKineticUncertainty = .unknown, evidence: VivoKineticEvidence) {
        self.value = value; self.unit = unit; self.origin = origin
        self.uncertainty = uncertainty; self.evidence = evidence
    }
    public func validate(unit expected: VivoKineticUnit, label: String, positive: Bool = false) throws {
        guard unit == expected, value.isFinite, positive ? value > 0 : value >= 0 else {
            throw VivoKineticsError.invalid("\(label) has invalid units, sign or finite representation")
        }
        try uncertainty.validate(value: value)
        try evidence.validate(origin: origin)
    }
}

/// Rates are conditional on this exact chemical and biological context. Version
/// one intentionally does not extrapolate rates across temperature, pH or hosts.
public struct VivoKineticContext: Codable, Equatable, Sendable {
    public let compound: String
    public let target: String
    public let targetVariant: String
    public let site: String
    public let chemicalState: String
    public let hostContext: String
    public let temperatureK: Double
    public let pH: Double
    public let ionicStrengthM: Double
    public init(compound: String, target: String, targetVariant: String, site: String,
                chemicalState: String, hostContext: String, temperatureK: Double,
                pH: Double, ionicStrengthM: Double) {
        self.compound = compound; self.target = target; self.targetVariant = targetVariant
        self.site = site; self.chemicalState = chemicalState; self.hostContext = hostContext
        self.temperatureK = temperatureK; self.pH = pH; self.ionicStrengthM = ionicStrengthM
    }
    public func validate() throws {
        guard [compound, target, targetVariant, site, chemicalState, hostContext].allSatisfy({
            !$0.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty && $0.utf8.count <= 1024
        }), temperatureK.isFinite, temperatureK > 0, pH.isFinite,
        ionicStrengthM.isFinite, ionicStrengthM >= 0 else {
            throw VivoKineticsError.invalid("chemical/target/context identity or environmental condition is missing")
        }
    }
}

public struct VivoKineticCompetitor: Codable, Equatable, Sendable {
    public let identifier: String
    public let association: VivoKineticParameter
    public let dissociation: VivoKineticParameter
    public let unboundConcentration: VivoKineticParameter
    public init(identifier: String, association: VivoKineticParameter,
                dissociation: VivoKineticParameter, unboundConcentration: VivoKineticParameter) {
        self.identifier = identifier; self.association = association
        self.dissociation = dissociation; self.unboundConcentration = unboundConcentration
    }
    public func validate() throws {
        guard !identifier.isEmpty, identifier.utf8.count <= 1024 else { throw VivoKineticsError.invalid("competitor identity") }
        try association.validate(unit: .perMolarSecond, label: "competitor association")
        try dissociation.validate(unit: .perSecond, label: "competitor dissociation")
        try unboundConcentration.validate(unit: .molar, label: "competitor unbound concentration")
    }
}

public enum VivoExposureTreatment: String, Codable, Sendable {
    /// Prescribed free drug is not depleted by target binding. Its source must
    /// justify the reservoir approximation; this is NOT finite drug mass balance.
    case externallyMaintainedUnbound
}
public enum VivoTargetTurnoverModel: String, Codable, Sendable {
    case equalStateLossAndConstantSynthesis
}

public struct VivoCovalentKineticPack: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let identifier: String
    public let context: VivoKineticContext
    public let association: VivoKineticParameter
    public let dissociation: VivoKineticParameter
    public let inactivation: VivoKineticParameter
    public let targetTurnover: VivoKineticParameter
    public let baselineTarget: VivoKineticParameter
    public let competitor: VivoKineticCompetitor?
    public let maximumUnboundDrugM: Double
    public let exposureTreatment: VivoExposureTreatment
    public let turnoverModel: VivoTargetTurnoverModel

    public init(identifier: String, context: VivoKineticContext, association: VivoKineticParameter,
                dissociation: VivoKineticParameter, inactivation: VivoKineticParameter,
                targetTurnover: VivoKineticParameter, baselineTarget: VivoKineticParameter,
                competitor: VivoKineticCompetitor? = nil, maximumUnboundDrugM: Double) {
        schemaVersion = 1; self.identifier = identifier; self.context = context
        self.association = association; self.dissociation = dissociation; self.inactivation = inactivation
        self.targetTurnover = targetTurnover; self.baselineTarget = baselineTarget; self.competitor = competitor
        self.maximumUnboundDrugM = maximumUnboundDrugM
        exposureTreatment = .externallyMaintainedUnbound; turnoverModel = .equalStateLossAndConstantSynthesis
    }
    public func validate() throws {
        guard schemaVersion == 1, !identifier.isEmpty, identifier.utf8.count <= 1024,
              maximumUnboundDrugM.isFinite, maximumUnboundDrugM >= 0 else {
            throw VivoKineticsError.invalid("schema, model identity or exposure domain")
        }
        try context.validate()
        try association.validate(unit: .perMolarSecond, label: "association")
        try dissociation.validate(unit: .perSecond, label: "dissociation")
        try inactivation.validate(unit: .perSecond, label: "inactivation")
        try targetTurnover.validate(unit: .perSecond, label: "target turnover")
        try baselineTarget.validate(unit: .molar, label: "target abundance", positive: true)
        try competitor?.validate()
        let maxAssociation = association.value * maximumUnboundDrugM
        let competitorAssociation = (competitor?.association.value ?? 0) * (competitor?.unboundConcentration.value ?? 0)
        guard [maxAssociation, competitorAssociation, maxAssociation + competitorAssociation + targetTurnover.value,
               dissociation.value + inactivation.value + targetTurnover.value,
               (competitor?.dissociation.value ?? 0) + targetTurnover.value,
               targetTurnover.value * baselineTarget.value].allSatisfy(\.isFinite) else {
            throw VivoKineticsError.invalid("derived propensities overflow FP64")
        }
    }
    public var parameters: [VivoKineticParameter] {
        [association, dissociation, inactivation, targetTurnover, baselineTarget]
            + (competitor.map { [$0.association, $0.dissociation, $0.unboundConcentration] } ?? [])
    }
    public var containsAssumptions: Bool { parameters.contains { $0.origin == .assumed } }
    public var hasUnknownParameterUncertainty: Bool { parameters.contains { $0.uncertainty == .unknown } }
}

public struct VivoExposureKnot: Codable, Equatable, Sendable {
    public let timeSeconds: Double
    public let unboundDrugM: Double
    public init(timeSeconds: Double, unboundDrugM: Double) {
        self.timeSeconds = timeSeconds; self.unboundDrugM = unboundDrugM
    }
}

public enum VivoExposureInterpolation: String, Codable, Sendable { case rightContinuousPiecewiseConstant }

public struct VivoUnboundExposureTrace: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let context: VivoKineticContext
    public let interpolation: VivoExposureInterpolation
    public let knots: [VivoExposureKnot]
    public let origin: VivoKineticOrigin
    public let evidence: VivoKineticEvidence
    public init(context: VivoKineticContext, knots: [VivoExposureKnot],
                origin: VivoKineticOrigin, evidence: VivoKineticEvidence) {
        schemaVersion = 1; self.context = context; self.knots = knots
        interpolation = .rightContinuousPiecewiseConstant; self.origin = origin; self.evidence = evidence
    }
    public func validate(for model: VivoCovalentKineticPack) throws {
        try model.validate(); try evidence.validate(origin: origin)
        guard schemaVersion == 1, context == model.context, (2...262_144).contains(knots.count),
              knots.first?.timeSeconds == 0 else {
            throw VivoKineticsError.invalid("exposure context, schema, capacity or initial time")
        }
        var previous = -Double.infinity
        for knot in knots {
            guard knot.timeSeconds.isFinite, knot.timeSeconds > previous,
                  knot.unboundDrugM.isFinite, knot.unboundDrugM >= 0,
                  knot.unboundDrugM <= model.maximumUnboundDrugM else {
                throw VivoKineticsError.invalid("exposure times must increase; free concentration must remain in domain")
            }
            previous = knot.timeSeconds
        }
    }
}

public struct VivoTargetFractions: Codable, Equatable, Sendable {
    public let free: Double
    public let reversible: Double
    public let covalent: Double
    public let competitor: Double
    public init(free: Double = 1, reversible: Double = 0, covalent: Double = 0, competitor: Double = 0) {
        self.free = free; self.reversible = reversible; self.covalent = covalent; self.competitor = competitor
    }
    public var values: [Double] { [free, reversible, covalent, competitor] }
    public var total: Double { free + reversible + covalent + competitor }
    public var drugOccupancy: Double { (reversible + covalent) / total }
    public func validate(hasCompetitor: Bool, tolerance: Double = 1e-10) throws {
        guard tolerance.isFinite, tolerance > 0, tolerance <= 1e-3,
              values.allSatisfy({ $0.isFinite && $0 >= 0 }), abs(total - 1) <= tolerance,
              hasCompetitor || competitor == 0 else {
            throw VivoKineticsError.invalid("initial target fractions must be nonnegative and sum to one")
        }
    }
}

public struct VivoTargetEngagementExperiment: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let kinetics: VivoCovalentKineticPack
    public let exposure: VivoUnboundExposureTrace
    public let initial: VivoTargetFractions
    public let sampleTimesSeconds: [Double]
    public init(kinetics: VivoCovalentKineticPack, exposure: VivoUnboundExposureTrace,
                initial: VivoTargetFractions = .init(), sampleTimesSeconds: [Double]) {
        schemaVersion = 1; self.kinetics = kinetics; self.exposure = exposure
        self.initial = initial; self.sampleTimesSeconds = sampleTimesSeconds
    }
    public func validate() throws {
        try exposure.validate(for: kinetics); try initial.validate(hasCompetitor: kinetics.competitor != nil)
        guard schemaVersion == 1, (1...262_144).contains(sampleTimesSeconds.count),
              let finalTime = exposure.knots.last?.timeSeconds else {
            throw VivoKineticsError.invalid("experiment schema or sample capacity")
        }
        var previous = -Double.infinity
        for time in sampleTimesSeconds {
            guard time.isFinite, time >= 0, time > previous, time <= finalTime else {
                throw VivoKineticsError.invalid("observation times must increase and lie inside exposure support")
            }
            previous = time
        }
    }
}
