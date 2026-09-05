import Foundation

public struct VivoPosteriorPredictivePolicy: Codable, Equatable, Sendable {
    public let intervalProbability: Double
    public let maximumCaseParticleRuns: Int
    public let maximumBufferedValues: Int
    public init(intervalProbability: Double = 0.95, maximumCaseParticleRuns: Int = 65_536,
                maximumBufferedValues: Int = 2_097_152) {
        self.intervalProbability = intervalProbability; self.maximumCaseParticleRuns = maximumCaseParticleRuns
        self.maximumBufferedValues = maximumBufferedValues
    }
    public func validate() throws {
        guard intervalProbability.isFinite, (0.5...0.999).contains(intervalProbability),
              (1...1_048_576).contains(maximumCaseParticleRuns), (1...16_777_216).contains(maximumBufferedValues) else {
            throw VivoPosteriorError.invalid("predictive policy")
        }
    }
}

public struct VivoPosteriorParameterSummary: Codable, Equatable, Sendable {
    public let identifier: String
    public let unit: String
    public let mean: Double
    public let standardDeviation: Double
    public let lower: Double
    public let median: Double
    public let upper: Double
}

public struct VivoPosteriorPredictivePoint: Codable, Equatable, Sendable {
    public let observationIdentifier: String
    public let timeSeconds: Double
    public let observable: VivoOccupancyObservable
    public let mean: Double
    public let latentStandardDeviation: Double
    public let credibleLower: Double
    public let credibleMedian: Double
    public let credibleUpper: Double
    /// Adds the explicitly supplied Gaussian assay SD to the latent-prediction
    /// mixture. Nil means measurement uncertainty is unknown, not zero.
    public let predictiveLower: Double?
    public let predictiveUpper: Double?
    public let observation: Double
    public let observationOrigin: VivoKineticOrigin
    public let measurementSD: Double?
    public let logPredictiveDensity: Double?
    public let predictiveIntervalContainsObservation: Bool?
}

public struct VivoPosteriorPredictiveFailure: Codable, Equatable, Sendable {
    public let particleIndex: Int
    public let message: String
}

public struct VivoPosteriorPredictiveCase: Codable, Equatable, Sendable {
    public let identifier: String
    public let partition: VivoStudyPartition
    public let leakageGroup: String
    public let evaluatedParticles: Int
    public let points: [VivoPosteriorPredictivePoint]
    public let failures: [VivoPosteriorPredictiveFailure]
}

public struct VivoPosteriorPredictiveReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let posteriorFingerprint: VivoFingerprint
    public let policy: VivoPosteriorPredictivePolicy
    public let parameterSummaries: [VivoPosteriorParameterSummary]
    /// Correlation of prior coordinates in parameter order. Nil means undefined
    /// (a collapsed coordinate), never an invented zero correlation. Every
    /// forward evaluation uses an intact joint particle, not marginal draws.
    public let priorCoordinateCorrelationRowMajor: [Double?]
    public let cases: [VivoPosteriorPredictiveCase]
    public let limitations: [String]
}

public enum VivoTargetPosteriorPredictor {
    private struct ParticleOutput: Sendable {
        let index: Int
        let values: [Double]?
        let failure: String?
    }

    public static func predict(_ record: VivoTargetPosteriorRecord,
                               policy: VivoPosteriorPredictivePolicy = .init()) async throws -> VivoPosteriorPredictiveReport {
        try record.validate(requireComplete: true); try policy.validate()
        let prepared = try VivoPreparedTargetPosterior(record.problem)
        guard let checkpoint = record.posterior.checkpoint else { throw VivoPosteriorError.invalid("missing posterior population") }
        let particles = checkpoint.particles, n = particles.count
        let work = n.multipliedReportingOverflow(by: record.problem.study.cases.count)
        guard !work.overflow, work.partialValue <= policy.maximumCaseParticleRuns else {
            throw VivoPosteriorError.budget("posterior prediction case-particle budget")
        }
        let physical = try particles.map { try prepared.physicalValues($0) }
        let tail = (1 - policy.intervalProbability) * 0.5
        let summaries = try prepared.plan.parameters.indices.map { column -> VivoPosteriorParameterSummary in
            let values = physical.map { $0[column] }.sorted()
            let (mean, sd) = try moments(values)
            return .init(identifier: prepared.plan.parameters[column].identifier, unit: prepared.plan.parameters[column].unit,
                mean: mean, standardDeviation: sd, lower: VivoPosteriorNumerics.quantile(sorted: values, probability: tail),
                median: VivoPosteriorNumerics.quantile(sorted: values, probability: 0.5),
                upper: VivoPosteriorNumerics.quantile(sorted: values, probability: 1 - tail))
        }
        var cases: [VivoPosteriorPredictiveCase] = []
        for item in record.problem.study.cases {
            try Task.checkCancellation()
            let count = n.multipliedReportingOverflow(by: item.observations.count)
            guard !count.overflow, count.partialValue <= policy.maximumBufferedValues else {
                throw VivoPosteriorError.budget("posterior prediction buffering capacity")
            }
            let outputs = try await evaluate(item, values: physical, prepared: prepared)
            var failures = outputs.compactMap { output -> VivoPosteriorPredictiveFailure? in
                output.failure.map { .init(particleIndex: output.index, message: $0) }
            }
            var points: [VivoPosteriorPredictivePoint] = []
            if failures.isEmpty {
                do {
                    for (column, observation) in item.observations.enumerated() {
                        let values = try outputs.map { output -> Double in
                            guard let row = output.values, row.count == item.observations.count else {
                                throw VivoPosteriorError.numerical("missing predictive particle output")
                            }
                            return row[column]
                        }.sorted()
                        points.append(try summarize(values, observation: observation, tail: tail))
                    }
                } catch {
                    failures = [.init(particleIndex: -1, message: "Summary failure: " + String(error.localizedDescription.prefix(4096)))]
                    points = []
                }
            }
            cases.append(.init(identifier: item.identifier, partition: item.partition, leakageGroup: item.leakageGroup,
                evaluatedParticles: outputs.count, points: points, failures: failures))
        }
        return try .init(schemaVersion: 1, posteriorFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(record)),
            policy: policy, parameterSummaries: summaries,
            priorCoordinateCorrelationRowMajor: correlations(particles), cases: cases, limitations: [
                "Joint finite-particle posterior propagated without independently sampling parameter marginals.",
                "Credible intervals describe latent occupancy conditional on model, priors, exposure representation and training data.",
                "Predictive intervals additionally include declared independent Gaussian measurement noise; unmodelled biological discrepancy is not covered.",
                "Intervals are pointwise, not simultaneous across times, targets or tissues. No clinical efficacy/safety inference is made.",
                "Coverage on calibration cases is in-sample. Validation/test labels are retained and synthetic observations are not experimental evidence.",
                "Any failed particle suppresses that case's intervals; particle filtering would change the posterior distribution.",
                "All final SMC particles have equal weights, but may remain genealogically or statistically dependent. Repeated-run diagnostics remain necessary."
            ])
    }

    private static func evaluate(_ item: VivoTargetEngagementStudyCase, values: [[Double]],
                                 prepared: VivoPreparedTargetPosterior) async throws -> [ParticleOutput] {
        try await withThrowingTaskGroup(of: ParticleOutput.self) { group in
            func add(_ index: Int) {
                group.addTask {
                    do {
                        try Task.checkCancellation()
                        let experiment = try prepared.experiment(for: item, values: values[index])
                        let result = try VivoTargetEngagementReference.run(experiment, numerics: prepared.problem.numerics)
                        let samples = Dictionary(uniqueKeysWithValues: result.samples.map { ($0.timeSeconds, $0) })
                        let predictions = try item.observations.map { observation -> Double in
                            guard let sample = samples[observation.timeSeconds] else {
                                throw VivoPosteriorError.numerical("missing predictive time")
                            }
                            let prediction = observation.observable.value(sample)
                            guard prediction.isFinite else { throw VivoPosteriorError.numerical("nonfinite predictive observable") }
                            return prediction
                        }
                        return ParticleOutput(index: index, values: predictions, failure: nil)
                    } catch is CancellationError { throw CancellationError() }
                    catch { return ParticleOutput(index: index, values: nil, failure: String(error.localizedDescription.prefix(4096))) }
                }
            }
            var next = 0, results: [ParticleOutput] = []
            for _ in 0..<min(prepared.problem.parallelEvaluations, values.count) { add(next); next += 1 }
            while let result = try await group.next() {
                results.append(result)
                if next < values.count { add(next); next += 1 }
            }
            return results.sorted { $0.index < $1.index }
        }
    }

    private static func moments(_ values: [Double]) throws -> (Double, Double) {
        guard let anchor = values.first, values.allSatisfy(\.isFinite) else { throw VivoPosteriorError.invalid("empty/nonfinite posterior summary") }
        let n = Double(values.count)
        let mean = anchor + values.reduce(0) { $0 + ($1 - anchor) / n }
        let scale = values.map { abs($0 - mean) }.max() ?? 0
        let variance = scale > 0 ? values.reduce(0) { partial, value in
            let normalized = (value - mean) / scale
            return partial + normalized * normalized / n
        } : 0
        let sd = scale * sqrt(variance)
        guard mean.isFinite, sd.isFinite else { throw VivoPosteriorError.numerical("posterior summary overflow") }
        return (mean, sd)
    }

    private static func summarize(_ sorted: [Double], observation: VivoOccupancyObservation,
                                  tail: Double) throws -> VivoPosteriorPredictivePoint {
        let n = Double(sorted.count), (mean, sdLatent) = try moments(sorted)
        var lower: Double?, upper: Double?, density: Double?, covered: Bool?
        if let sd = observation.standardDeviation {
            lower = try VivoPosteriorNumerics.gaussianMixtureQuantile(means: sorted, sd: sd, probability: tail)
            upper = try VivoPosteriorNumerics.gaussianMixtureQuantile(means: sorted, sd: sd, probability: 1 - tail)
            let logs = sorted.map { prediction -> Double in
                let z = (observation.value - prediction) / sd
                return -0.5 * z * z - log(sd) - 0.5 * log(2 * Double.pi)
            }
            guard let maximum = logs.max(), maximum.isFinite else {
                throw VivoPosteriorError.numerical("predictive density outside finite numerical support")
            }
            let sum = logs.reduce(0) { $0 + exp($1 - maximum) / n }
            density = maximum + log(sum)
            guard density?.isFinite == true else { throw VivoPosteriorError.numerical("predictive density overflow") }
            covered = observation.value >= lower! && observation.value <= upper!
        }
        return .init(observationIdentifier: observation.identifier, timeSeconds: observation.timeSeconds,
            observable: observation.observable, mean: mean, latentStandardDeviation: sdLatent,
            credibleLower: VivoPosteriorNumerics.quantile(sorted: sorted, probability: tail),
            credibleMedian: VivoPosteriorNumerics.quantile(sorted: sorted, probability: 0.5),
            credibleUpper: VivoPosteriorNumerics.quantile(sorted: sorted, probability: 1 - tail),
            predictiveLower: lower, predictiveUpper: upper, observation: observation.value,
            observationOrigin: observation.origin, measurementSD: observation.standardDeviation,
            logPredictiveDensity: density, predictiveIntervalContainsObservation: covered)
    }

    private static func correlations(_ particles: [VivoPosteriorParticle]) throws -> [Double?] {
        let d = particles[0].coordinates.count, n = Double(particles.count)
        let columns = (0..<d).map { j in particles.map { $0.coordinates[j] } }
        let statistics = try columns.map(moments)
        var result = [Double?](repeating: nil, count: d * d)
        for i in 0..<d { for j in 0..<d {
            let (mi, si) = statistics[i], (mj, sj) = statistics[j]
            if si == 0 || sj == 0 { continue }
            let value = zip(columns[i], columns[j]).reduce(0) {
                $0 + (($1.0 - mi) / si) * (($1.1 - mj) / sj) / n
            }
            guard value.isFinite else { throw VivoPosteriorError.numerical("posterior correlation overflow") }
            result[i * d + j] = value
        } }
        return result
    }
}
