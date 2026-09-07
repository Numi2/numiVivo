import Foundation

public struct VivoQMMMChemicalExchangeValidationTarget: Codable, Sendable, Equatable {
    public var identifier: String
    public var timeSeconds: Double
    public var survivalProbability: Double
    public var standardDeviationProbability: Double?
    public var context: VivoKineticContext
    public var source: VivoKineticEvidence

    public init(identifier: String, timeSeconds: Double, survivalProbability: Double,
                standardDeviationProbability: Double? = nil, context: VivoKineticContext,
                source: VivoKineticEvidence) {
        self.identifier = identifier
        self.timeSeconds = timeSeconds
        self.survivalProbability = survivalProbability
        self.standardDeviationProbability = standardDeviationProbability
        self.context = context
        self.source = source
    }

    public func validate() throws {
        try context.validate()
        try source.validate(origin: .measured)
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              timeSeconds.isFinite, timeSeconds >= 0,
              survivalProbability.isFinite, (0...1).contains(survivalProbability) else {
            throw VivoKineticsError.invalid("chemical-exchange validation target")
        }
        if let standardDeviationProbability {
            guard standardDeviationProbability.isFinite,
                  standardDeviationProbability > 0,
                  standardDeviationProbability <= 1 else {
                throw VivoKineticsError.invalid("chemical-exchange validation uncertainty")
            }
        }
    }
}

public struct VivoQMMMChemicalExchangeValidationRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-exchange-validation/v1"
    public var schema: String
    public var identifier: String
    public var networkRequest: VivoQMMMChemicalExchangeNetworkRequest
    public var networkResult: VivoQMMMChemicalExchangeNetworkResult
    public var targets: [VivoQMMMChemicalExchangeValidationTarget]
    public var maximumAbsoluteProbabilityError: Double
    public var maximumStandardizedResidual: Double

    public init(identifier: String,
                networkRequest: VivoQMMMChemicalExchangeNetworkRequest,
                networkResult: VivoQMMMChemicalExchangeNetworkResult,
                targets: [VivoQMMMChemicalExchangeValidationTarget],
                maximumAbsoluteProbabilityError: Double = 0.10,
                maximumStandardizedResidual: Double = 3.0) {
        self.schema = Self.schema
        self.identifier = identifier
        self.networkRequest = networkRequest
        self.networkResult = networkResult
        self.targets = targets
        self.maximumAbsoluteProbabilityError = maximumAbsoluteProbabilityError
        self.maximumStandardizedResidual = maximumStandardizedResidual
    }
}

public struct VivoQMMMChemicalExchangeValidationPoint: Codable, Sendable, Equatable {
    public let identifier: String
    public let timeSeconds: Double
    public let predictedSurvivalProbability: Double?
    public let measuredSurvivalProbability: Double
    public let absoluteProbabilityError: Double?
    public let standardizedResidual: Double?
    public let comparable: Bool
    public let passed: Bool?
    public let reason: String
}

public struct VivoQMMMChemicalExchangeValidationResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-chemical-exchange-validation-result/v1"
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let points: [VivoQMMMChemicalExchangeValidationPoint]
    public let converged: Bool
    public let issues: [String]
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Direct validation for transient, potentially non-single-exponential chemical
/// state kinetics. It compares condition-matched measured survival fractions only
/// at observation times explicitly present in the immutable exchange-network run.
public enum VivoQMMMChemicalExchangeValidation {
    public static let interpretation = "Condition-matched time-course validation of an explicit chemical-state exchange network against measured survival fractions. This validates the declared transient kinetic model at observed times; it does not convert a multi-exponential curve into a single first-order rate."

    private static func sameEnvironment(_ a: VivoKineticContext, _ b: VivoKineticContext) -> Bool {
        a.compound == b.compound && a.target == b.target && a.targetVariant == b.targetVariant &&
        a.site == b.site && a.hostContext == b.hostContext && a.temperatureK == b.temperatureK &&
        a.pH == b.pH && a.ionicStrengthM == b.ionicStrengthM
    }

    public static func calculate(_ request: VivoQMMMChemicalExchangeValidationRequest) throws -> VivoQMMMChemicalExchangeValidationResult {
        try VivoQMMMChemicalExchangeNetwork.validate(request.networkResult, request: request.networkRequest)
        guard request.schema == VivoQMMMChemicalExchangeValidationRequest.schema,
              !request.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.identifier.utf8.count <= 512,
              !request.targets.isEmpty, request.targets.count <= 16_384,
              Set(request.targets.map(\.identifier)).count == request.targets.count,
              request.maximumAbsoluteProbabilityError.isFinite,
              request.maximumAbsoluteProbabilityError > 0,
              request.maximumAbsoluteProbabilityError <= 1,
              request.maximumStandardizedResidual.isFinite,
              request.maximumStandardizedResidual > 0 else {
            throw VivoKineticsError.invalid("chemical-exchange validation request")
        }
        for target in request.targets { try target.validate() }
        guard let reference = request.networkRequest.states.first?.pathways.first?.request.replicas.first?.context else {
            throw VivoKineticsError.invalid("chemical-exchange validation has no simulated context")
        }
        let observations = Dictionary(uniqueKeysWithValues: request.networkResult.observations.map { ($0.timeSeconds.bitPattern, $0) })
        var points: [VivoQMMMChemicalExchangeValidationPoint] = []
        var issues: [String] = []
        points.reserveCapacity(request.targets.count)

        for target in request.targets {
            guard sameEnvironment(reference, target.context) else {
                points.append(.init(identifier: target.identifier, timeSeconds: target.timeSeconds,
                    predictedSurvivalProbability: nil, measuredSurvivalProbability: target.survivalProbability,
                    absoluteProbabilityError: nil, standardizedResidual: nil, comparable: false, passed: nil,
                    reason: "experimental kinetic context differs from the simulated chemical-state ensemble"))
                continue
            }
            guard let observation = observations[target.timeSeconds.bitPattern] else {
                points.append(.init(identifier: target.identifier, timeSeconds: target.timeSeconds,
                    predictedSurvivalProbability: nil, measuredSurvivalProbability: target.survivalProbability,
                    absoluteProbabilityError: nil, standardizedResidual: nil, comparable: false, passed: nil,
                    reason: "target time is not an explicit observation time in the immutable exchange-network result"))
                continue
            }
            let error = abs(observation.survivalProbability - target.survivalProbability)
            let standardized = target.standardDeviationProbability.map { error / $0 }
            let passed = error <= request.maximumAbsoluteProbabilityError &&
                (standardized == nil || standardized! <= request.maximumStandardizedResidual)
            points.append(.init(identifier: target.identifier, timeSeconds: target.timeSeconds,
                predictedSurvivalProbability: observation.survivalProbability,
                measuredSurvivalProbability: target.survivalProbability, absoluteProbabilityError: error,
                standardizedResidual: standardized, comparable: true, passed: passed,
                reason: "condition-matched measured survival fraction at an explicit simulated observation time"))
            if !passed { issues.append("external transient validation \(target.identifier) exceeds its predeclared probability-error gate") }
        }

        if !points.contains(where: { $0.comparable }) {
            issues.append("no condition-matched transient external validation target is directly comparable")
        }
        let requestID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request))
        struct Evidence: Codable {
            let schema: String
            let request: VivoQMMMChemicalExchangeValidationRequest
            let points: [VivoQMMMChemicalExchangeValidationPoint]
            let issues: [String]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/qmmm-chemical-exchange-validation-evidence/v1",
            request: request, points: points, issues: issues)))
        return .init(schema: VivoQMMMChemicalExchangeValidationResult.schema,
                     requestFingerprint: requestID, points: points,
                     converged: issues.isEmpty, issues: issues,
                     interpretation: interpretation, evidenceFingerprint: evidenceID)
    }

    public static func validate(_ result: VivoQMMMChemicalExchangeValidationResult,
                                request: VivoQMMMChemicalExchangeValidationRequest) throws {
        guard result == (try calculate(request)) else {
            throw VivoKineticsError.invalid("chemical-exchange external validation does not reconstruct")
        }
    }
}
