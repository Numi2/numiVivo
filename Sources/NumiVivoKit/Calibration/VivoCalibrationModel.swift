import Foundation

public enum VivoCalibrationTransform: String, Codable, Sendable, CaseIterable {
    case linear
    case logarithmic
    case logit
}

public enum VivoCalibrationPrior: Codable, Sendable, Equatable {
    case uniform
    case normal(mean: Double, standardDeviation: Double)
    case logNormal(logMean: Double, logStandardDeviation: Double)
    case triangular(mode: Double)

    private enum CodingKeys: String, CodingKey {
        case kind, mean, standardDeviation, logMean, logStandardDeviation, mode
    }

    private enum Kind: String, Codable {
        case uniform, normal, logNormal, triangular
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .uniform:
            self = .uniform
        case .normal:
            self = .normal(
                mean: try container.decode(Double.self, forKey: .mean),
                standardDeviation: try container.decode(Double.self, forKey: .standardDeviation)
            )
        case .logNormal:
            self = .logNormal(
                logMean: try container.decode(Double.self, forKey: .logMean),
                logStandardDeviation: try container.decode(Double.self, forKey: .logStandardDeviation)
            )
        case .triangular:
            self = .triangular(mode: try container.decode(Double.self, forKey: .mode))
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .uniform:
            try container.encode(Kind.uniform, forKey: .kind)
        case .normal(let mean, let standardDeviation):
            try container.encode(Kind.normal, forKey: .kind)
            try container.encode(mean, forKey: .mean)
            try container.encode(standardDeviation, forKey: .standardDeviation)
        case .logNormal(let logMean, let logStandardDeviation):
            try container.encode(Kind.logNormal, forKey: .kind)
            try container.encode(logMean, forKey: .logMean)
            try container.encode(logStandardDeviation, forKey: .logStandardDeviation)
        case .triangular(let mode):
            try container.encode(Kind.triangular, forKey: .kind)
            try container.encode(mode, forKey: .mode)
        }
    }

    public func validate(lower: Double, upper: Double, label: String) throws {
        switch self {
        case .uniform:
            return
        case .normal(let mean, let deviation):
            guard mean.isFinite, deviation.isFinite, deviation > 0 else {
                throw VivoArtifactValidationError.invalid("\(label) normal prior requires finite mean and positive deviation")
            }
        case .logNormal(let mean, let deviation):
            guard mean.isFinite, deviation.isFinite, deviation > 0, lower > 0 else {
                throw VivoArtifactValidationError.invalid("\(label) log-normal prior requires positive bounds and deviation")
            }
        case .triangular(let mode):
            guard mode.isFinite, (lower...upper).contains(mode) else {
                throw VivoArtifactValidationError.invalid("\(label) triangular mode must lie inside parameter bounds")
            }
        }
    }
}

public struct VivoCalibrationParameter: Codable, Sendable, Equatable {
    public var identifier: String
    public var lowerBound: VivoQuantity
    public var upperBound: VivoQuantity
    public var initialValue: VivoQuantity?
    public var transform: VivoCalibrationTransform
    public var prior: VivoCalibrationPrior
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        lowerBound: VivoQuantity,
        upperBound: VivoQuantity,
        initialValue: VivoQuantity? = nil,
        transform: VivoCalibrationTransform = .linear,
        prior: VivoCalibrationPrior = .uniform,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.initialValue = initialValue
        self.transform = transform
        self.prior = prior
        self.evidence = evidence
    }
}

public enum VivoObservationNoise: Codable, Sendable, Equatable {
    case gaussian(standardDeviation: Double)
    case logNormal(logStandardDeviation: Double)
    case studentT(scale: Double, degreesOfFreedom: Double)
    case interval(lower: Double, upper: Double, outsideScale: Double)

    private enum CodingKeys: String, CodingKey {
        case kind, standardDeviation, logStandardDeviation, scale, degreesOfFreedom, lower, upper, outsideScale
    }

    private enum Kind: String, Codable {
        case gaussian, logNormal, studentT, interval
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .gaussian:
            self = .gaussian(standardDeviation: try container.decode(Double.self, forKey: .standardDeviation))
        case .logNormal:
            self = .logNormal(logStandardDeviation: try container.decode(Double.self, forKey: .logStandardDeviation))
        case .studentT:
            self = .studentT(
                scale: try container.decode(Double.self, forKey: .scale),
                degreesOfFreedom: try container.decode(Double.self, forKey: .degreesOfFreedom)
            )
        case .interval:
            self = .interval(
                lower: try container.decode(Double.self, forKey: .lower),
                upper: try container.decode(Double.self, forKey: .upper),
                outsideScale: try container.decode(Double.self, forKey: .outsideScale)
            )
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .gaussian(let deviation):
            try container.encode(Kind.gaussian, forKey: .kind)
            try container.encode(deviation, forKey: .standardDeviation)
        case .logNormal(let deviation):
            try container.encode(Kind.logNormal, forKey: .kind)
            try container.encode(deviation, forKey: .logStandardDeviation)
        case .studentT(let scale, let degreesOfFreedom):
            try container.encode(Kind.studentT, forKey: .kind)
            try container.encode(scale, forKey: .scale)
            try container.encode(degreesOfFreedom, forKey: .degreesOfFreedom)
        case .interval(let lower, let upper, let outsideScale):
            try container.encode(Kind.interval, forKey: .kind)
            try container.encode(lower, forKey: .lower)
            try container.encode(upper, forKey: .upper)
            try container.encode(outsideScale, forKey: .outsideScale)
        }
    }

    public func validate(label: String) throws {
        switch self {
        case .gaussian(let deviation), .logNormal(let deviation):
            guard deviation.isFinite, deviation > 0 else {
                throw VivoArtifactValidationError.invalid("\(label) deviation must be finite and positive")
            }
        case .studentT(let scale, let degreesOfFreedom):
            guard scale.isFinite, scale > 0,
                  degreesOfFreedom.isFinite, degreesOfFreedom > 0 else {
                throw VivoArtifactValidationError.invalid("\(label) Student-t scale and degrees of freedom must be positive")
            }
        case .interval(let lower, let upper, let outsideScale):
            guard lower.isFinite, upper.isFinite, lower <= upper,
                  outsideScale.isFinite, outsideScale > 0 else {
                throw VivoArtifactValidationError.invalid("\(label) interval noise contract is invalid")
            }
        }
    }
}

public struct VivoCalibrationObservation: Codable, Sendable, Equatable {
    public var identifier: String
    public var measurementIdentifier: String
    public var timeSeconds: Double
    public var maximumTimeDifferenceSeconds: Double
    public var values: [Double]
    public var unit: String
    public var noise: VivoObservationNoise
    public var weight: Double
    public var replicateIndices: [UInt32]?
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        measurementIdentifier: String,
        timeSeconds: Double,
        maximumTimeDifferenceSeconds: Double,
        values: [Double],
        unit: String,
        noise: VivoObservationNoise,
        weight: Double = 1,
        replicateIndices: [UInt32]? = nil,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.measurementIdentifier = measurementIdentifier
        self.timeSeconds = timeSeconds
        self.maximumTimeDifferenceSeconds = maximumTimeDifferenceSeconds
        self.values = values
        self.unit = unit
        self.noise = noise
        self.weight = weight
        self.replicateIndices = replicateIndices
        self.evidence = evidence
    }
}

public struct VivoCalibrationStrategy: Codable, Sendable, Equatable {
    public var populationSize: Int
    public var eliteFraction: Double
    public var smoothing: Double
    public var explorationFraction: Double
    public var minimumNormalizedDeviation: Double
    public var maximumEvaluations: Int
    public var maximumGenerations: Int
    public var stagnationGenerations: Int
    public var objectiveTolerance: Double
    public var seed: UInt64

    public init(
        populationSize: Int = 128,
        eliteFraction: Double = 0.15,
        smoothing: Double = 0.35,
        explorationFraction: Double = 0.10,
        minimumNormalizedDeviation: Double = 0.005,
        maximumEvaluations: Int = 20_000,
        maximumGenerations: Int = 200,
        stagnationGenerations: Int = 20,
        objectiveTolerance: Double = 1e-6,
        seed: UInt64 = 0x4e756d695669766f
    ) {
        self.populationSize = populationSize
        self.eliteFraction = eliteFraction
        self.smoothing = smoothing
        self.explorationFraction = explorationFraction
        self.minimumNormalizedDeviation = minimumNormalizedDeviation
        self.maximumEvaluations = maximumEvaluations
        self.maximumGenerations = maximumGenerations
        self.stagnationGenerations = stagnationGenerations
        self.objectiveTolerance = objectiveTolerance
        self.seed = seed
    }

    public func validate(parameterCount: Int) throws {
        let doubled = parameterCount.multipliedReportingOverflow(by: 2)
        guard parameterCount > 0, !doubled.overflow,
              populationSize >= max(doubled.partialValue, 8),
              eliteFraction.isFinite, eliteFraction > 0, eliteFraction <= 0.5,
              smoothing.isFinite, smoothing > 0, smoothing <= 1,
              explorationFraction.isFinite, explorationFraction >= 0, explorationFraction <= 1,
              minimumNormalizedDeviation.isFinite, minimumNormalizedDeviation > 0, minimumNormalizedDeviation < 0.5,
              maximumEvaluations >= populationSize,
              maximumGenerations > 0,
              stagnationGenerations > 0,
              objectiveTolerance.isFinite, objectiveTolerance >= 0 else {
            throw VivoArtifactValidationError.invalid("calibration strategy is inconsistent with the parameter count")
        }
    }
}

public struct VivoCalibrationPack: Codable, Sendable, Equatable {
    public var programFingerprint: VivoFingerprint
    public var hostContextFingerprint: VivoFingerprint
    public var experimentFingerprint: VivoFingerprint
    public var parameters: [VivoCalibrationParameter]
    public var observations: [VivoCalibrationObservation]
    public var strategy: VivoCalibrationStrategy
    public var failurePenalty: Double
    public var labels: [String: String]

    public init(
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint,
        experimentFingerprint: VivoFingerprint,
        parameters: [VivoCalibrationParameter],
        observations: [VivoCalibrationObservation],
        strategy: VivoCalibrationStrategy = .init(),
        failurePenalty: Double = 1e30,
        labels: [String: String] = [:]
    ) {
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.parameters = parameters
        self.observations = observations
        self.strategy = strategy
        self.failurePenalty = failurePenalty
        self.labels = labels
    }
}

public struct PreparedVivoCalibrationParameter: Codable, Sendable, Equatable {
    public let identifier: String
    public let parameterIndex: UInt32
    public let unit: String
    public let lowerBound: Double
    public let upperBound: Double
    public let initialValue: Double
    public let transform: VivoCalibrationTransform
    public let prior: VivoCalibrationPrior

    public func value(normalized: Double) -> Double {
        let coordinate = min(max(normalized, 0), 1)
        switch transform {
        case .linear:
            return lowerBound + coordinate * (upperBound - lowerBound)
        case .logarithmic:
            return exp(log(lowerBound) + coordinate * (log(upperBound) - log(lowerBound)))
        case .logit:
            let lowerLogit = log(lowerBound / (1 - lowerBound))
            let upperLogit = log(upperBound / (1 - upperBound))
            let logit = lowerLogit + coordinate * (upperLogit - lowerLogit)
            return 1 / (1 + exp(-logit))
        }
    }

    public func normalized(value: Double) -> Double {
        switch transform {
        case .linear:
            return (value - lowerBound) / (upperBound - lowerBound)
        case .logarithmic:
            return (log(value) - log(lowerBound)) / (log(upperBound) - log(lowerBound))
        case .logit:
            let valueLogit = log(value / (1 - value))
            let lowerLogit = log(lowerBound / (1 - lowerBound))
            let upperLogit = log(upperBound / (1 - upperBound))
            return (valueLogit - lowerLogit) / (upperLogit - lowerLogit)
        }
    }
}

public struct PreparedVivoCalibration: Codable, Sendable, Equatable {
    public let fingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint
    public let experimentFingerprint: VivoFingerprint
    public let parameters: [PreparedVivoCalibrationParameter]
    public let observations: [VivoCalibrationObservation]
    public let strategy: VivoCalibrationStrategy
    public let failurePenalty: Double
    public let labels: [String: String]
}

public struct VivoCalibrationCompiler: Sendable {
    private let units: VivoUnitSystem

    public init(units: VivoUnitSystem = .standard) {
        self.units = units
    }

    public func compile(
        _ calibration: VivoCalibrationPack,
        programPack: VivoProgramPack,
        hostContext: PreparedVivoHostContext,
        experiment: PreparedVivoExperiment
    ) throws -> PreparedVivoCalibration {
        guard calibration.programFingerprint == programPack.header.contentFingerprint,
              calibration.hostContextFingerprint == hostContext.contextFingerprint,
              calibration.experimentFingerprint == experiment.fingerprint else {
            throw VivoArtifactValidationError.incompatible("calibration artifact fingerprints do not match program, host context, and experiment")
        }
        guard !calibration.parameters.isEmpty, !calibration.observations.isEmpty,
              calibration.failurePenalty.isFinite, calibration.failurePenalty > 0 else {
            throw VivoArtifactValidationError.invalid("calibration requires parameters, observations, and a positive finite failure penalty")
        }
        try calibration.strategy.validate(parameterCount: calibration.parameters.count)

        let metadata = try programPack.parameterMetadata()
        guard metadata.count <= Int(UInt32.max), Set(metadata.map(\.identifier)).count == metadata.count else {
            throw VivoArtifactValidationError.invalid("parameter metadata has duplicate identifiers or exceeds UInt32 indexing")
        }
        // ParameterMetadata intentionally has no independent index authority.
        // Its array order is the compiled parameter-table order.
        let byIdentifier = Dictionary(uniqueKeysWithValues: metadata.enumerated().map {
            ($0.element.identifier, (index: UInt32($0.offset), metadata: $0.element))
        })
        var identifiers = Set<String>()
        var prepared: [PreparedVivoCalibrationParameter] = []
        prepared.reserveCapacity(calibration.parameters.count)
        for parameter in calibration.parameters {
            guard !parameter.identifier.isEmpty, identifiers.insert(parameter.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("calibration parameter identifiers must be non-empty and unique")
            }
            guard let selected = byIdentifier[parameter.identifier] else {
                throw VivoArtifactValidationError.unresolved("calibration parameter \(parameter.identifier) is absent from ProgramPack")
            }
            let target = selected.metadata
            try parameter.lowerBound.validate(label: "calibration.parameters.\(parameter.identifier).lowerBound")
            try parameter.upperBound.validate(label: "calibration.parameters.\(parameter.identifier).upperBound")
            try parameter.initialValue?.validate(label: "calibration.parameters.\(parameter.identifier).initialValue")
            guard units.areCompatible(parameter.lowerBound.unit, target.unit),
                  units.areCompatible(parameter.upperBound.unit, target.unit),
                  parameter.initialValue.map({ units.areCompatible($0.unit, target.unit) }) ?? true else {
                throw VivoArtifactValidationError.incompatible("calibration parameter \(parameter.identifier) units are incompatible with \(target.unit)")
            }
            let lower = try units.convert(parameter.lowerBound, to: target.unit)
            let upper = try units.convert(parameter.upperBound, to: target.unit)
            let initial = try parameter.initialValue.map { try units.convert($0, to: target.unit) } ?? target.value
            guard lower.isFinite, upper.isFinite, lower < upper, (lower...upper).contains(initial),
                  lower >= target.minimum, upper <= target.maximum else {
                throw VivoArtifactValidationError.invalid("calibration bounds for \(parameter.identifier) are outside ProgramPack bounds or exclude the initial value")
            }
            if parameter.transform == .logarithmic, lower <= 0 {
                throw VivoArtifactValidationError.invalid("logarithmic calibration parameter \(parameter.identifier) requires positive bounds")
            }
            if parameter.transform == .logit,
               !(lower > 0 && upper < 1) {
                throw VivoArtifactValidationError.invalid("logit calibration parameter \(parameter.identifier) requires bounds strictly inside (0,1)")
            }
            try parameter.prior.validate(lower: lower, upper: upper, label: "calibration.parameters.\(parameter.identifier).prior")
            prepared.append(.init(
                identifier: parameter.identifier,
                parameterIndex: selected.index,
                unit: target.unit,
                lowerBound: lower,
                upperBound: upper,
                initialValue: initial,
                transform: parameter.transform,
                prior: parameter.prior
            ))
        }

        var observationIdentifiers = Set<String>()
        let measurementIdentifiers = Set(experiment.measurements.map { $0.identifier })
        for observation in calibration.observations {
            guard !observation.identifier.isEmpty,
                  observationIdentifiers.insert(observation.identifier).inserted,
                  measurementIdentifiers.contains(observation.measurementIdentifier),
                  observation.timeSeconds.isFinite, observation.timeSeconds >= 0,
                  observation.timeSeconds <= experiment.durationSeconds,
                  observation.maximumTimeDifferenceSeconds.isFinite,
                  observation.maximumTimeDifferenceSeconds >= 0,
                  !observation.values.isEmpty,
                  observation.values.allSatisfy(\.isFinite),
                  !observation.unit.isEmpty,
                  observation.weight.isFinite, observation.weight > 0 else {
                throw VivoArtifactValidationError.invalid("calibration observation \(observation.identifier) is invalid or unresolved")
            }
            if let replicateIndices = observation.replicateIndices {
                guard !replicateIndices.isEmpty,
                      Set(replicateIndices).count == replicateIndices.count,
                      replicateIndices.allSatisfy({ $0 < experiment.replicateCount }) else {
                    throw VivoArtifactValidationError.invalid("calibration observation \(observation.identifier) has invalid replicate selection")
                }
            }
            try observation.noise.validate(label: "calibration.observations.\(observation.identifier).noise")
        }

        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(calibration))
        return PreparedVivoCalibration(
            fingerprint: fingerprint,
            programFingerprint: calibration.programFingerprint,
            hostContextFingerprint: calibration.hostContextFingerprint,
            experimentFingerprint: calibration.experimentFingerprint,
            parameters: prepared.sorted { $0.parameterIndex < $1.parameterIndex },
            observations: calibration.observations.sorted { $0.identifier < $1.identifier },
            strategy: calibration.strategy,
            failurePenalty: calibration.failurePenalty,
            labels: calibration.labels
        )
    }
}
