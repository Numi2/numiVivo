import Foundation

public struct VivoCalibrationCandidate: Codable, Sendable, Equatable, Hashable {
    public let identifier: VivoFingerprint
    public let generation: Int
    public let ordinal: Int
    public let normalizedCoordinates: [Double]
    public let parameterValues: [Double]
    public let parameterIndices: [UInt32]

    public init(
        generation: Int,
        ordinal: Int,
        normalizedCoordinates: [Double],
        parameters: [PreparedVivoCalibrationParameter]
    ) throws {
        guard generation >= 0,
              ordinal >= 0,
              normalizedCoordinates.count == parameters.count,
              normalizedCoordinates.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }) else {
            throw VivoArtifactValidationError.invalid("calibration candidate coordinates are invalid")
        }
        self.generation = generation
        self.ordinal = ordinal
        self.normalizedCoordinates = normalizedCoordinates
        self.parameterValues = zip(parameters, normalizedCoordinates).map { $0.value(normalized: $1) }
        self.parameterIndices = parameters.map(\.parameterIndex)
        struct Identity: Codable {
            let generation: Int
            let ordinal: Int
            let coordinates: [Double]
            let parameterIndices: [UInt32]
        }
        self.identifier = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            generation: generation,
            ordinal: ordinal,
            coordinates: normalizedCoordinates,
            parameterIndices: parameterIndices
        )))
    }

    public var indexedValues: [UInt32: Double] {
        Dictionary(uniqueKeysWithValues: zip(parameterIndices, parameterValues))
    }
}

public struct VivoCalibrationTerm: Codable, Sendable, Equatable, Hashable {
    public let identifier: String
    public let contribution: Double
    public let matchedTimeSeconds: Double?
    public let predictionCount: Int
    public let observationCount: Int
    public let note: String
}

public enum VivoCalibrationEvaluationStatus: String, Codable, Sendable {
    case succeeded
    case runtimeRejected
    case missingObservation
    case invalidPrediction
    case failed
}

public struct VivoCalibrationEvaluation: Codable, Sendable, Equatable {
    public let candidate: VivoCalibrationCandidate
    public let status: VivoCalibrationEvaluationStatus
    public let objective: Double
    public let likelihoodObjective: Double
    public let priorObjective: Double
    public let terms: [VivoCalibrationTerm]
    public let resultFingerprint: VivoFingerprint?
    public let message: String

    public init(
        candidate: VivoCalibrationCandidate,
        status: VivoCalibrationEvaluationStatus,
        objective: Double,
        likelihoodObjective: Double,
        priorObjective: Double,
        terms: [VivoCalibrationTerm],
        resultFingerprint: VivoFingerprint?,
        message: String = ""
    ) {
        self.candidate = candidate
        self.status = status
        self.objective = objective
        self.likelihoodObjective = likelihoodObjective
        self.priorObjective = priorObjective
        self.terms = terms
        self.resultFingerprint = resultFingerprint
        self.message = message
    }

    public var isUsable: Bool {
        status == .succeeded && objective.isFinite
    }
}

public struct VivoLikelihoodEvaluator: Sendable {
    private let units: VivoUnitSystem

    public init(units: VivoUnitSystem = .standard) {
        self.units = units
    }

    public func evaluate(
        candidate: VivoCalibrationCandidate,
        result: VivoResultPack,
        calibration: PreparedVivoCalibration
    ) throws -> VivoCalibrationEvaluation {
        guard result.programFingerprint == calibration.programFingerprint,
              result.hostContextFingerprint == calibration.hostContextFingerprint,
              result.experimentFingerprint == calibration.experimentFingerprint else {
            throw VivoArtifactValidationError.incompatible("calibration result fingerprints do not match the prepared calibration")
        }
        guard candidate.parameterIndices == calibration.parameters.map(\.parameterIndex),
              candidate.parameterValues.count == calibration.parameters.count else {
            throw VivoArtifactValidationError.incompatible("candidate parameter layout differs from prepared calibration")
        }
        try result.validate()

        var terms: [VivoCalibrationTerm] = []
        var likelihood = 0.0
        for observation in calibration.observations {
            let match = try prediction(for: observation, in: result)
            guard let match else {
                return penalty(
                    candidate: candidate,
                    calibration: calibration,
                    status: .missingObservation,
                    message: "No sample matched observation \(observation.identifier).",
                    terms: terms
                )
            }
            guard let aligned = align(predictions: match.values, observations: observation.values) else {
                return penalty(
                    candidate: candidate,
                    calibration: calibration,
                    status: .invalidPrediction,
                    message: "Prediction shape does not match observation \(observation.identifier).",
                    terms: terms
                )
            }
            guard let contribution = negativeLogLikelihood(
                predictions: aligned.predictions,
                observations: aligned.observations,
                noise: observation.noise
            ), contribution.isFinite else {
                return penalty(
                    candidate: candidate,
                    calibration: calibration,
                    status: .invalidPrediction,
                    message: "Observation model rejected prediction for \(observation.identifier).",
                    terms: terms
                )
            }
            let weighted = contribution * observation.weight
            likelihood += weighted
            terms.append(.init(
                identifier: observation.identifier,
                contribution: weighted,
                matchedTimeSeconds: match.time,
                predictionCount: aligned.predictions.count,
                observationCount: aligned.observations.count,
                note: match.note
            ))
        }

        let prior = priorObjective(candidate: candidate, parameters: calibration.parameters)
        let objective = likelihood + prior
        guard objective.isFinite else {
            return penalty(
                candidate: candidate,
                calibration: calibration,
                status: .invalidPrediction,
                message: "Combined likelihood and prior objective is non-finite.",
                terms: terms
            )
        }
        return VivoCalibrationEvaluation(
            candidate: candidate,
            status: .succeeded,
            objective: objective,
            likelihoodObjective: likelihood,
            priorObjective: prior,
            terms: terms,
            resultFingerprint: try result.fingerprint()
        )
    }

    public func penalty(
        candidate: VivoCalibrationCandidate,
        calibration: PreparedVivoCalibration,
        status: VivoCalibrationEvaluationStatus,
        message: String,
        terms: [VivoCalibrationTerm] = []
    ) -> VivoCalibrationEvaluation {
        VivoCalibrationEvaluation(
            candidate: candidate,
            status: status,
            objective: calibration.failurePenalty,
            likelihoodObjective: calibration.failurePenalty,
            priorObjective: 0,
            terms: terms,
            resultFingerprint: nil,
            message: message
        )
    }

    private struct PredictionMatch {
        let values: [Double]
        let time: Double
        let note: String
    }

    private func prediction(
        for observation: VivoCalibrationObservation,
        in result: VivoResultPack
    ) throws -> PredictionMatch? {
        let replicateFilter = observation.replicateIndices.map(Set.init)
        let candidates = result.measurements.filter { sample in
            sample.measurementIdentifier == observation.measurementIdentifier &&
            (replicateFilter?.contains(sample.replicateIndex) ?? true)
        }
        guard let nearestDistance = candidates.map({ abs($0.timeSeconds - observation.timeSeconds) }).min(),
              nearestDistance <= observation.maximumTimeDifferenceSeconds else {
            return nil
        }
        let tolerance = max(1e-12, nearestDistance * 1e-9)
        let nearest = candidates.filter {
            abs(abs($0.timeSeconds - observation.timeSeconds) - nearestDistance) <= tolerance
        }
        guard !nearest.isEmpty else { return nil }

        var converted: [Double] = []
        for sample in nearest.sorted(by: {
            if $0.replicateIndex != $1.replicateIndex { return $0.replicateIndex < $1.replicateIndex }
            return $0.stepIndex < $1.stepIndex
        }) {
            guard units.areCompatible(sample.unit, observation.unit) else {
                throw VivoArtifactValidationError.incompatible(
                    "measurement \(observation.measurementIdentifier) unit \(sample.unit) is incompatible with observation unit \(observation.unit)"
                )
            }
            for value in sample.values {
                converted.append(try units.convert(
                    VivoQuantity(value: value, unit: sample.unit),
                    to: observation.unit
                ))
            }
        }
        let time = nearest.map(\.timeSeconds).reduce(0, +) / Double(nearest.count)
        return PredictionMatch(
            values: converted,
            time: time,
            note: "aggregated \(nearest.count) nearest sample record(s)"
        )
    }

    private func align(
        predictions: [Double],
        observations: [Double]
    ) -> (predictions: [Double], observations: [Double])? {
        guard !predictions.isEmpty, !observations.isEmpty else { return nil }
        if predictions.count == observations.count {
            return (predictions, observations)
        }
        if observations.count == 1 {
            return ([predictions.reduce(0, +) / Double(predictions.count)], observations)
        }
        if predictions.count == 1 {
            return (predictions, [observations.reduce(0, +) / Double(observations.count)])
        }
        return nil
    }

    private func negativeLogLikelihood(
        predictions: [Double],
        observations: [Double],
        noise: VivoObservationNoise
    ) -> Double? {
        guard predictions.count == observations.count,
              predictions.allSatisfy(\.isFinite),
              observations.allSatisfy(\.isFinite) else { return nil }
        let logTwoPi = log(2 * Double.pi)
        switch noise {
        case .gaussian(let standardDeviation):
            let variance = standardDeviation * standardDeviation
            return zip(predictions, observations).reduce(0) { total, pair in
                let residual = pair.0 - pair.1
                return total + 0.5 * (residual * residual / variance + logTwoPi + log(variance))
            }
        case .logNormal(let logStandardDeviation):
            let variance = logStandardDeviation * logStandardDeviation
            guard predictions.allSatisfy({ $0 > 0 }), observations.allSatisfy({ $0 > 0 }) else { return nil }
            return zip(predictions, observations).reduce(0) { total, pair in
                let residual = log(pair.0) - log(pair.1)
                return total + 0.5 * (residual * residual / variance + logTwoPi + log(variance)) + log(pair.0)
            }
        case .studentT(let scale, let degreesOfFreedom):
            let normalizer = lgamma((degreesOfFreedom + 1) / 2) - lgamma(degreesOfFreedom / 2) -
                0.5 * log(degreesOfFreedom * Double.pi) - log(scale)
            return zip(predictions, observations).reduce(0) { total, pair in
                let residual = (pair.0 - pair.1) / scale
                let logDensity = normalizer - ((degreesOfFreedom + 1) / 2) * log1p(residual * residual / degreesOfFreedom)
                return total - logDensity
            }
        case .interval(let lower, let upper, let outsideScale):
            let scaleSquared = outsideScale * outsideScale
            return predictions.reduce(0) { total, value in
                let distance: Double
                if value < lower { distance = lower - value }
                else if value > upper { distance = value - upper }
                else { distance = 0 }
                return total + 0.5 * distance * distance / scaleSquared
            }
        }
    }

    private func priorObjective(
        candidate: VivoCalibrationCandidate,
        parameters: [PreparedVivoCalibrationParameter]
    ) -> Double {
        zip(parameters, candidate.parameterValues).reduce(0) { total, pair in
            total + priorNegativeLogDensity(value: pair.1, parameter: pair.0)
        }
    }

    private func priorNegativeLogDensity(
        value: Double,
        parameter: PreparedVivoCalibrationParameter
    ) -> Double {
        guard value >= parameter.lowerBound, value <= parameter.upperBound else {
            return Double.greatestFiniteMagnitude
        }
        switch parameter.prior {
        case .uniform:
            return log(parameter.upperBound - parameter.lowerBound)
        case .normal(let mean, let deviation):
            let residual = (value - mean) / deviation
            return 0.5 * residual * residual + log(deviation) + 0.5 * log(2 * Double.pi)
        case .logNormal(let logMean, let deviation):
            guard value > 0 else { return Double.greatestFiniteMagnitude }
            let residual = (log(value) - logMean) / deviation
            return 0.5 * residual * residual + log(deviation) + 0.5 * log(2 * Double.pi) + log(value)
        case .triangular(let mode):
            let width = parameter.upperBound - parameter.lowerBound
            let density: Double
            if value == mode {
                density = 2 / width
            } else if value < mode {
                let denominator = width * (mode - parameter.lowerBound)
                density = denominator > 0 ? 2 * (value - parameter.lowerBound) / denominator : 0
            } else {
                let denominator = width * (parameter.upperBound - mode)
                density = denominator > 0 ? 2 * (parameter.upperBound - value) / denominator : 0
            }
            return density > 0 ? -log(density) : Double.greatestFiniteMagnitude
        }
    }
}
