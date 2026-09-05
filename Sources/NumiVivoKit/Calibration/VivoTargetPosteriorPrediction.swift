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
    public let predictiveLower: Double?
    public let predictiveUpper: Double?
    public let observation: Double
    public let observationOrigin: VivoKineticOrigin
    /// Reported measurement SD only; inferred residual noise is in the joint
    /// posterior. A null reported SD is never relabelled as a laboratory estimate.
    public let measurementSD: Double?
    public let logPredictiveDensity: Double?
    public let predictiveIntervalContainsObservation: Bool?
    public let observationSupport: VivoAssaySupport?
    public let logObservationProbability: Double?
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
    /// Mixture of the complete case likelihood, including declared within-case
    /// correlation. Not the sum of marginal pointwise log predictive densities.
    public let jointLogObservationLikelihood: Double?
}

public struct VivoPosteriorPredictiveReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let posteriorFingerprint: VivoFingerprint
    public let policy: VivoPosteriorPredictivePolicy
    public let parameterSummaries: [VivoPosteriorParameterSummary]
    public let priorCoordinateCorrelationRowMajor: [Double?]
    public let cases: [VivoPosteriorPredictiveCase]
    public let limitations: [String]
}

public enum VivoTargetPosteriorPredictor {
    private struct ParticleOutput: Sendable {
        let index: Int
        let latent: [Double]?
        let noise: VivoGaussianAssayNoise?
        let standardDeviations: [Double?]?
        let jointLogLikelihood: Double?
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
            let values = physical.map { $0[column] }.sorted(), stats = try moments(values)
            return .init(identifier: prepared.plan.parameters[column].identifier, unit: prepared.plan.parameters[column].unit,
                mean: stats.0, standardDeviation: stats.1, lower: VivoPosteriorNumerics.quantile(sorted: values, probability: tail),
                median: VivoPosteriorNumerics.quantile(sorted: values, probability: 0.5),
                upper: VivoPosteriorNumerics.quantile(sorted: values, probability: 1 - tail))
        }
        var cases: [VivoPosteriorPredictiveCase] = []
        for item in record.problem.study.cases {
            try Task.checkCancellation()
            let count = n.multipliedReportingOverflow(by: item.observations.count)
            guard !count.overflow, count.partialValue <= policy.maximumBufferedValues / 4 else {
                throw VivoPosteriorError.budget("posterior prediction buffering capacity")
            }
            let outputs = try await evaluate(item, values: physical, prepared: prepared)
            var failures = outputs.compactMap { output in
                output.failure.map { VivoPosteriorPredictiveFailure(particleIndex: output.index, message: $0) }
            }
            var points: [VivoPosteriorPredictivePoint] = [], joint: Double?
            if failures.isEmpty {
                do {
                    for (column, observation) in item.observations.enumerated() {
                        let support = record.problem.assaySupport(caseIdentifier: item.identifier, observationIdentifier: observation.identifier)
                        points.append(try summarize(outputs, column: column, observation: observation, support: support, tail: tail))
                    }
                    let logs = outputs.compactMap(\.jointLogLikelihood)
                    if logs.count == n { joint = try logMeanExp(logs) }
                } catch {
                    failures = [.init(particleIndex: -1, message: String(error.localizedDescription.prefix(4096)))]
                    points = []; joint = nil
                }
            }
            cases.append(.init(identifier: item.identifier, partition: item.partition, leakageGroup: item.leakageGroup,
                evaluatedParticles: outputs.count, points: points, failures: failures, jointLogObservationLikelihood: joint))
        }
        return try .init(schemaVersion: 2, posteriorFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(record)),
            policy: policy, parameterSummaries: summaries, priorCoordinateCorrelationRowMajor: correlations(particles), cases: cases,
            limitations: [
                "Joint particles preserve associations between kinetic parameters, exposure, assay noise and bias.",
                "Latent intervals describe occupancy; predictive intervals additionally include the declared Gaussian assay model and its parameter uncertainty.",
                "Intervals are pointwise, not simultaneous. Joint case scoring includes declared serial correlation; marginal scores must not be added as a joint score.",
                "Censored observations yield event probabilities, not density values or coverage tests on substituted measurements.",
                "Unresolved measurement uncertainty leaves predictive fields null. Reported SD remains separate from inferred residual SD.",
                "Any failed particle suppresses the affected case rather than conditioning predictions on successful simulations.",
                "Model discrepancy is included only to the extent explicitly represented. No clinical efficacy, patient safety or universal uncertainty-coverage claim is made."
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
                        let noise = try prepared.problem.assayNoise(caseIdentifier: item.identifier, values: values[index])
                        let predictions = try item.observations.map { observation -> Double in
                            guard let sample = samples[observation.timeSeconds] else { throw VivoPosteriorError.numerical("missing predictive time") }
                            let value = observation.observable.value(sample)
                            guard value.isFinite else { throw VivoPosteriorError.numerical("nonfinite predictive observable") }
                            return value
                        }
                        let scales = try item.observations.map { observation -> Double? in
                            if observation.standardDeviation == nil && noise.additionalStandardDeviation == 0 { return nil }
                            return try noise.standardDeviation(reported: observation.standardDeviation)
                        }
                        var joint: Double?
                        if scales.allSatisfy({ $0 != nil }) {
                            let data = item.observations.map { observation in
                                VivoGaussianAssayDatum(timeSeconds: observation.timeSeconds, value: observation.value,
                                    reportedStandardDeviation: observation.standardDeviation,
                                    support: prepared.problem.assaySupport(caseIdentifier: item.identifier, observationIdentifier: observation.identifier))
                            }
                            joint = try VivoGaussianAssay.logLikelihood(predictions: predictions, observations: data, noise: noise)
                        }
                        return ParticleOutput(index: index, latent: predictions, noise: noise, standardDeviations: scales, jointLogLikelihood: joint, failure: nil)
                    } catch is CancellationError { throw CancellationError() }
                    catch { return ParticleOutput(index: index, latent: nil, noise: nil, standardDeviations: nil, jointLogLikelihood: nil, failure: String(error.localizedDescription.prefix(4096))) }
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

    private static func summarize(_ outputs: [ParticleOutput], column: Int, observation: VivoOccupancyObservation,
                                  support: VivoAssaySupport, tail: Double) throws -> VivoPosteriorPredictivePoint {
        var latent: [Double] = [], means: [Double] = [], scales: [Double] = []
        for output in outputs {
            guard let values = output.latent, column < values.count, let noise = output.noise,
                  let sd = output.standardDeviations, column < sd.count else { throw VivoPosteriorError.numerical("predictive particle shape") }
            latent.append(values[column]); means.append(values[column] + noise.bias)
            if let value = sd[column] { scales.append(value) }
        }
        latent.sort(); let stats = try moments(latent)
        var lower: Double?, upper: Double?, density: Double?, probability: Double?, covered: Bool?
        if scales.count == outputs.count {
            lower = try VivoGaussianAssay.mixtureQuantile(means: means, standardDeviations: scales, probability: tail)
            upper = try VivoGaussianAssay.mixtureQuantile(means: means, standardDeviations: scales, probability: 1 - tail)
            let logs = try means.indices.map { try VivoGaussianAssay.logContribution(mean: means[$0], sd: scales[$0], value: observation.value, support: support) }
            let score = try logMeanExp(logs)
            if support == .exact { density = score; covered = observation.value >= lower! && observation.value <= upper! }
            else { probability = score }
        }
        return .init(observationIdentifier: observation.identifier, timeSeconds: observation.timeSeconds, observable: observation.observable,
            mean: stats.0, latentStandardDeviation: stats.1, credibleLower: VivoPosteriorNumerics.quantile(sorted: latent, probability: tail),
            credibleMedian: VivoPosteriorNumerics.quantile(sorted: latent, probability: 0.5), credibleUpper: VivoPosteriorNumerics.quantile(sorted: latent, probability: 1 - tail),
            predictiveLower: lower, predictiveUpper: upper, observation: observation.value, observationOrigin: observation.origin,
            measurementSD: observation.standardDeviation, logPredictiveDensity: density, predictiveIntervalContainsObservation: covered,
            observationSupport: support, logObservationProbability: probability)
    }
    private static func logMeanExp(_ values: [Double]) throws -> Double {
        guard let maximum = values.max(), values.allSatisfy(\.isFinite) else { throw VivoPosteriorError.numerical("invalid predictive likelihood mixture") }
        let result = maximum + log(values.reduce(0) { $0 + exp($1 - maximum) / Double(values.count) })
        guard result.isFinite else { throw VivoPosteriorError.numerical("predictive mixture overflow") }
        return result
    }
    private static func moments(_ values: [Double]) throws -> (Double, Double) {
        guard let anchor = values.first, values.allSatisfy(\.isFinite) else { throw VivoPosteriorError.invalid("empty/nonfinite posterior summary") }
        let n = Double(values.count), mean = anchor + values.reduce(0) { $0 + ($1 - anchor) / n }
        let scale = values.map { abs($0 - mean) }.max() ?? 0
        let variance = scale > 0 ? values.reduce(0) { sum, value in
            let z = (value - mean) / scale; return sum + z * z / n
        } : 0
        let sd = scale * sqrt(variance)
        guard mean.isFinite, sd.isFinite else { throw VivoPosteriorError.numerical("posterior summary overflow") }
        return (mean, sd)
    }
    private static func correlations(_ particles: [VivoPosteriorParticle]) throws -> [Double?] {
        let d = particles[0].coordinates.count, n = Double(particles.count)
        let columns = (0..<d).map { j in particles.map { $0.coordinates[j] } }, statistics = try columns.map(moments)
        var result = [Double?](repeating: nil, count: d * d)
        for i in 0..<d { for j in 0..<d {
            let (mi, si) = statistics[i], (mj, sj) = statistics[j]
            if si == 0 || sj == 0 { continue }
            let value = zip(columns[i], columns[j]).reduce(0) { $0 + (($1.0 - mi) / si) * (($1.1 - mj) / sj) / n }
            guard value.isFinite else { throw VivoPosteriorError.numerical("posterior correlation overflow") }
            result[i * d + j] = value
        } }
        return result
    }
}
