import Foundation

public struct VivoCalibrationBatch: Sendable, Equatable {
    public let generation: Int
    public let candidates: [VivoCalibrationCandidate]

    public init(generation: Int, candidates: [VivoCalibrationCandidate]) {
        self.generation = generation
        self.candidates = candidates
    }

    public func parameterMajorFP32(
        baseParameters: [Double],
        parameterCount: Int
    ) throws -> [Float] {
        guard parameterCount > 0,
              baseParameters.count == parameterCount,
              !candidates.isEmpty,
              candidates.count <= Int(UInt32.max),
              candidates.allSatisfy({ $0.parameterIndices.count == $0.parameterValues.count }) else {
            throw VivoArtifactValidationError.invalid("calibration batch parameter layout is invalid")
        }
        let elementCount = parameterCount.multipliedReportingOverflow(by: candidates.count)
        guard !elementCount.overflow else {
            throw VivoArtifactValidationError.invalid("calibration parameter matrix size overflow")
        }
        var matrix = [Float](repeating: 0, count: elementCount.partialValue)
        for parameter in 0..<parameterCount {
            let base = baseParameters[parameter]
            guard base.isFinite, abs(base) <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoArtifactValidationError.invalid("base calibration parameter \(parameter) is not FP32 representable")
            }
            for environment in candidates.indices {
                matrix[parameter * candidates.count + environment] = Float(base)
            }
        }
        for (environment, candidate) in candidates.enumerated() {
            for (index, parameterIndex) in candidate.parameterIndices.enumerated() {
                guard parameterIndex < UInt32(parameterCount) else {
                    throw VivoArtifactValidationError.invalid("candidate parameter index exceeds base parameter count")
                }
                let value = candidate.parameterValues[index]
                guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else {
                    throw VivoArtifactValidationError.invalid("candidate value is not FP32 representable")
                }
                matrix[Int(parameterIndex) * candidates.count + environment] = Float(value)
            }
        }
        return matrix
    }
}

public struct VivoCalibrationGenerationSummary: Codable, Sendable, Equatable {
    public let generation: Int
    public let evaluatedCount: Int
    public let usableCount: Int
    public let bestObjective: Double
    public let medianObjective: Double
    public let eliteObjective: Double
    public let normalizedMean: [Double]
    public let normalizedDeviation: [Double]
}

public struct VivoCalibrationCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/calibration-checkpoint/v1"

    public let schema: String
    public let calibrationFingerprint: VivoFingerprint
    public let nextGeneration: Int
    public let evaluationCount: Int
    public let randomState: UInt64
    public let normalizedMean: [Double]
    public let normalizedDeviation: [Double]
    public let best: VivoCalibrationEvaluation?
    public let history: [VivoCalibrationGenerationSummary]
    public let elite: [VivoCalibrationEvaluation]

    public init(
        calibrationFingerprint: VivoFingerprint,
        nextGeneration: Int,
        evaluationCount: Int,
        randomState: UInt64,
        normalizedMean: [Double],
        normalizedDeviation: [Double],
        best: VivoCalibrationEvaluation?,
        history: [VivoCalibrationGenerationSummary],
        elite: [VivoCalibrationEvaluation]
    ) {
        self.schema = Self.schema
        self.calibrationFingerprint = calibrationFingerprint
        self.nextGeneration = nextGeneration
        self.evaluationCount = evaluationCount
        self.randomState = randomState
        self.normalizedMean = normalizedMean
        self.normalizedDeviation = normalizedDeviation
        self.best = best
        self.history = history
        self.elite = elite
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoCalibrationResult: Codable, Sendable, Equatable {
    public let calibrationFingerprint: VivoFingerprint
    public let startedAt: Date
    public let finishedAt: Date
    public let evaluationCount: Int
    public let generationCount: Int
    public let terminationReason: String
    public let best: VivoCalibrationEvaluation
    public let ensemble: [VivoCalibrationEvaluation]
    public let history: [VivoCalibrationGenerationSummary]
    public let limitations: [String]

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    public func envelope(
        metadata: VivoArtifactMetadata,
        createdBy: String,
        softwareVersion: String,
        sourceArtifacts: [VivoFingerprint]
    ) -> VivoArtifactEnvelope<VivoCalibrationResult> {
        VivoArtifactEnvelope(
            kind: .calibration,
            metadata: metadata,
            provenance: .init(
                createdBy: createdBy,
                softwareVersion: softwareVersion,
                sourceArtifacts: sourceArtifacts
            ),
            payload: self
        )
    }
}

public enum VivoCalibrationOptimizerError: Error, Sendable, CustomStringConvertible {
    case invalidCheckpoint(String)
    case evaluatorContract(String)
    case noUsableEvaluation

    public var description: String {
        switch self {
        case .invalidCheckpoint(let message): "Invalid calibration checkpoint: \(message)"
        case .evaluatorContract(let message): "Calibration evaluator contract violation: \(message)"
        case .noUsableEvaluation: "Calibration terminated without a usable evaluation."
        }
    }
}

public typealias VivoCalibrationBatchEvaluator = @Sendable (VivoCalibrationBatch) async throws -> [VivoCalibrationEvaluation]
public typealias VivoCalibrationProgressSink = @Sendable (VivoCalibrationCheckpoint) async throws -> Void

public actor VivoAdaptiveEnsembleOptimizer {
    public let calibration: PreparedVivoCalibration

    private var random: VivoSplitMix64
    private var nextGeneration: Int
    private var evaluationCount: Int
    private var normalizedMean: [Double]
    private var normalizedDeviation: [Double]
    private var best: VivoCalibrationEvaluation?
    private var history: [VivoCalibrationGenerationSummary]
    private var eliteArchive: [VivoCalibrationEvaluation]

    public init(calibration: PreparedVivoCalibration) {
        self.calibration = calibration
        self.random = VivoSplitMix64(state: calibration.strategy.seed)
        self.nextGeneration = 0
        self.evaluationCount = 0
        self.normalizedMean = calibration.parameters.map {
            min(max($0.normalized(value: $0.initialValue), 0), 1)
        }
        self.normalizedDeviation = [Double](repeating: 0.25, count: calibration.parameters.count)
        self.best = nil
        self.history = []
        self.eliteArchive = []
    }

    public init(
        calibration: PreparedVivoCalibration,
        checkpoint: VivoCalibrationCheckpoint
    ) throws {
        guard checkpoint.schema == VivoCalibrationCheckpoint.schema,
              checkpoint.calibrationFingerprint == calibration.fingerprint,
              checkpoint.nextGeneration >= 0,
              checkpoint.evaluationCount >= 0,
              checkpoint.normalizedMean.count == calibration.parameters.count,
              checkpoint.normalizedDeviation.count == calibration.parameters.count,
              checkpoint.normalizedMean.allSatisfy({ $0.isFinite && $0 >= 0 && $0 <= 1 }),
              checkpoint.normalizedDeviation.allSatisfy({ $0.isFinite && $0 > 0 && $0 <= 1 }) else {
            throw VivoCalibrationOptimizerError.invalidCheckpoint("identity, dimensions, or optimizer state is invalid")
        }
        self.calibration = calibration
        self.random = VivoSplitMix64(state: checkpoint.randomState)
        self.nextGeneration = checkpoint.nextGeneration
        self.evaluationCount = checkpoint.evaluationCount
        self.normalizedMean = checkpoint.normalizedMean
        self.normalizedDeviation = checkpoint.normalizedDeviation
        self.best = checkpoint.best
        self.history = checkpoint.history
        self.eliteArchive = checkpoint.elite
    }

    public func run(
        evaluate: VivoCalibrationBatchEvaluator,
        progress: VivoCalibrationProgressSink? = nil
    ) async throws -> VivoCalibrationResult {
        let startedAt = Date()
        let strategy = calibration.strategy
        var stagnantGenerations = 0
        var previousBest = best?.objective ?? Double.infinity
        var terminationReason = "maximum generations reached"

        while nextGeneration < strategy.maximumGenerations,
              evaluationCount < strategy.maximumEvaluations {
            let remaining = strategy.maximumEvaluations - evaluationCount
            let count = min(strategy.populationSize, remaining)
            let candidates = try makeCandidates(generation: nextGeneration, count: count)
            let batch = VivoCalibrationBatch(generation: nextGeneration, candidates: candidates)
            let returned = try await evaluate(batch)
            let evaluations = try validate(returned: returned, for: candidates)
            evaluationCount += evaluations.count

            let ranked = evaluations.sorted(by: evaluationOrdering)
            let usable = ranked.filter(\.isUsable)
            if let generationBest = usable.first,
               best == nil || generationBest.objective < best!.objective {
                best = generationBest
            }
            let adaptationPool = usable.isEmpty ? ranked : usable
            let eliteCount = min(
                adaptationPool.count,
                max(1, Int(ceil(Double(count) * strategy.eliteFraction)))
            )
            let elites = Array(adaptationPool.prefix(eliteCount))
            eliteArchive = mergeArchive(eliteArchive, with: elites, limit: strategy.populationSize)
            adapt(to: elites, smoothing: strategy.smoothing, floor: strategy.minimumNormalizedDeviation)

            let objectives = adaptationPool.map(\.objective).filter(\.isFinite).sorted()
            let median = objectives.isEmpty ? calibration.failurePenalty : objectives[objectives.count / 2]
            let eliteObjective = elites.isEmpty
                ? calibration.failurePenalty
                : elites.map(\.objective).reduce(0, +) / Double(elites.count)
            history.append(.init(
                generation: nextGeneration,
                evaluatedCount: evaluations.count,
                usableCount: usable.count,
                bestObjective: best?.objective ?? calibration.failurePenalty,
                medianObjective: median,
                eliteObjective: eliteObjective,
                normalizedMean: normalizedMean,
                normalizedDeviation: normalizedDeviation
            ))
            nextGeneration += 1

            let currentBest = best?.objective ?? previousBest
            if previousBest.isFinite,
               currentBest.isFinite,
               previousBest - currentBest <= strategy.objectiveTolerance {
                stagnantGenerations += 1
            } else {
                stagnantGenerations = 0
            }
            previousBest = currentBest

            if let progress {
                try await progress(checkpoint())
            }
            if stagnantGenerations >= strategy.stagnationGenerations {
                terminationReason = "objective stagnation threshold reached"
                break
            }
            if evaluationCount >= strategy.maximumEvaluations {
                terminationReason = "maximum evaluation budget reached"
                break
            }
        }

        guard let best, best.isUsable else {
            throw VivoCalibrationOptimizerError.noUsableEvaluation
        }
        return VivoCalibrationResult(
            calibrationFingerprint: calibration.fingerprint,
            startedAt: startedAt,
            finishedAt: Date(),
            evaluationCount: evaluationCount,
            generationCount: history.count,
            terminationReason: terminationReason,
            best: best,
            ensemble: eliteArchive.filter(\.isUsable).sorted(by: evaluationOrdering),
            history: history,
            limitations: [
                "The adaptive ensemble approximates a bounded objective landscape; it is not a proof of global optimality or parameter identifiability.",
                "Posterior interpretation requires a likelihood, prior, sampling procedure, and convergence analysis appropriate to the scientific question.",
                "Calibration agreement does not validate unobserved mechanisms or transfer to an undeclared host context."
            ]
        )
    }

    public func checkpoint() -> VivoCalibrationCheckpoint {
        VivoCalibrationCheckpoint(
            calibrationFingerprint: calibration.fingerprint,
            nextGeneration: nextGeneration,
            evaluationCount: evaluationCount,
            randomState: random.state,
            normalizedMean: normalizedMean,
            normalizedDeviation: normalizedDeviation,
            best: best,
            history: history,
            elite: eliteArchive
        )
    }

    private func makeCandidates(generation: Int, count: Int) throws -> [VivoCalibrationCandidate] {
        guard count > 0 else { return [] }
        if generation == 0 {
            return try latinHypercube(generation: generation, count: count)
        }

        let explorationCount = min(count, Int(ceil(Double(count) * calibration.strategy.explorationFraction)))
        var candidates: [VivoCalibrationCandidate] = []
        candidates.reserveCapacity(count)
        for ordinal in 0..<count {
            let coordinates: [Double]
            if ordinal < explorationCount {
                coordinates = calibration.parameters.map { _ in random.unitInterval() }
            } else {
                coordinates = normalizedMean.indices.map { index in
                    let proposal = normalizedMean[index] + random.normal() * normalizedDeviation[index]
                    return reflectUnit(proposal)
                }
            }
            candidates.append(try VivoCalibrationCandidate(
                generation: generation,
                ordinal: ordinal,
                normalizedCoordinates: coordinates,
                parameters: calibration.parameters
            ))
        }
        return candidates
    }

    private func latinHypercube(generation: Int, count: Int) throws -> [VivoCalibrationCandidate] {
        let dimensionCount = calibration.parameters.count
        var permutations: [[Int]] = []
        permutations.reserveCapacity(dimensionCount)
        for _ in 0..<dimensionCount {
            var values = Array(0..<count)
            random.shuffle(&values)
            permutations.append(values)
        }
        var candidates: [VivoCalibrationCandidate] = []
        candidates.reserveCapacity(count)
        for ordinal in 0..<count {
            var coordinates = [Double](repeating: 0, count: dimensionCount)
            for dimension in 0..<dimensionCount {
                coordinates[dimension] = (Double(permutations[dimension][ordinal]) + random.unitInterval()) / Double(count)
            }
            if ordinal == 0 {
                coordinates = calibration.parameters.map {
                    min(max($0.normalized(value: $0.initialValue), 0), 1)
                }
            }
            candidates.append(try VivoCalibrationCandidate(
                generation: generation,
                ordinal: ordinal,
                normalizedCoordinates: coordinates,
                parameters: calibration.parameters
            ))
        }
        return candidates
    }

    private func validate(
        returned: [VivoCalibrationEvaluation],
        for candidates: [VivoCalibrationCandidate]
    ) throws -> [VivoCalibrationEvaluation] {
        guard returned.count == candidates.count else {
            throw VivoCalibrationOptimizerError.evaluatorContract(
                "expected \(candidates.count) evaluations, received \(returned.count)"
            )
        }
        let expected = Set(candidates.map(\.identifier))
        let actual = Set(returned.map { $0.candidate.identifier })
        guard actual.count == returned.count, actual == expected else {
            throw VivoCalibrationOptimizerError.evaluatorContract(
                "candidate identities are missing, duplicated, or substituted"
            )
        }
        return returned.map { evaluation in
            guard evaluation.objective.isFinite else {
                return VivoCalibrationEvaluation(
                    candidate: evaluation.candidate,
                    status: .invalidPrediction,
                    objective: calibration.failurePenalty,
                    likelihoodObjective: calibration.failurePenalty,
                    priorObjective: 0,
                    terms: evaluation.terms,
                    resultFingerprint: evaluation.resultFingerprint,
                    message: evaluation.message.isEmpty ? "evaluator returned a non-finite objective" : evaluation.message
                )
            }
            return evaluation
        }
    }

    private func adapt(
        to elite: [VivoCalibrationEvaluation],
        smoothing: Double,
        floor: Double
    ) {
        guard !elite.isEmpty else {
            normalizedDeviation = normalizedDeviation.map { min(max($0 * 1.25, floor), 0.5) }
            return
        }
        for dimension in normalizedMean.indices {
            let values = elite.map { $0.candidate.normalizedCoordinates[dimension] }
            let mean = values.reduce(0, +) / Double(values.count)
            let variance = values.reduce(0) { partial, value in
                let difference = value - mean
                return partial + difference * difference
            } / Double(max(values.count - 1, 1))
            let deviation = max(sqrt(max(variance, 0)), floor)
            normalizedMean[dimension] = (1 - smoothing) * normalizedMean[dimension] + smoothing * mean
            normalizedDeviation[dimension] = min(
                max((1 - smoothing) * normalizedDeviation[dimension] + smoothing * deviation, floor),
                0.5
            )
        }
    }

    private func mergeArchive(
        _ current: [VivoCalibrationEvaluation],
        with additions: [VivoCalibrationEvaluation],
        limit: Int
    ) -> [VivoCalibrationEvaluation] {
        var byIdentifier = Dictionary(uniqueKeysWithValues: current.map { ($0.candidate.identifier, $0) })
        for value in additions {
            if let existing = byIdentifier[value.candidate.identifier] {
                if evaluationOrdering(value, existing) { byIdentifier[value.candidate.identifier] = value }
            } else {
                byIdentifier[value.candidate.identifier] = value
            }
        }
        return Array(byIdentifier.values.sorted(by: evaluationOrdering).prefix(limit))
    }

    private func evaluationOrdering(
        _ lhs: VivoCalibrationEvaluation,
        _ rhs: VivoCalibrationEvaluation
    ) -> Bool {
        if lhs.objective != rhs.objective { return lhs.objective < rhs.objective }
        return lhs.candidate.identifier.hex < rhs.candidate.identifier.hex
    }

    private func reflectUnit(_ value: Double) -> Double {
        guard value.isFinite else { return 0.5 }
        var coordinate = value
        while coordinate < 0 || coordinate > 1 {
            if coordinate < 0 { coordinate = -coordinate }
            if coordinate > 1 { coordinate = 2 - coordinate }
        }
        return min(max(coordinate, 0), 1)
    }
}

public struct VivoSplitMix64: Sendable, Codable, Equatable {
    public private(set) var state: UInt64

    public init(state: UInt64) {
        self.state = state
    }

    public mutating func next() -> UInt64 {
        state &+= 0x9E3779B97F4A7C15
        var value = state
        value = (value ^ (value >> 30)) &* 0xBF58476D1CE4E5B9
        value = (value ^ (value >> 27)) &* 0x94D049BB133111EB
        return value ^ (value >> 31)
    }

    public mutating func unitInterval() -> Double {
        let mantissa = next() >> 11
        return Double(mantissa) * 0x1.0p-53
    }

    public mutating func normal() -> Double {
        let u1 = max(unitInterval(), Double.leastNonzeroMagnitude)
        let u2 = unitInterval()
        return sqrt(-2 * log(u1)) * cos(2 * Double.pi * u2)
    }

    public mutating func shuffle<T>(_ values: inout [T]) {
        guard values.count > 1 else { return }
        for index in stride(from: values.count - 1, through: 1, by: -1) {
            let selected = Int(next() % UInt64(index + 1))
            if selected != index { values.swapAt(index, selected) }
        }
    }
}
