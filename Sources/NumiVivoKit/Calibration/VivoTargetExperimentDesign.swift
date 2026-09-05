import Foundation

/// A proposed measurement, with no observed outcome. Template selection explicitly
/// identifies which fitted context, shared parameters and assay model apply.
public struct VivoTargetDesignCandidate: Codable, Equatable, Sendable {
    public let identifier: String
    public let templateCaseIdentifier: String
    public let exposureKnots: [VivoExposureKnot]
    public let sampleTimeSeconds: Double
    public let observable: VivoOccupancyObservable
    public let reportedStandardDeviation: Double?
    public let relativeCost: Double
    public let evidence: VivoKineticEvidence
    public init(identifier: String, templateCaseIdentifier: String, exposureKnots: [VivoExposureKnot],
                sampleTimeSeconds: Double, observable: VivoOccupancyObservable,
                reportedStandardDeviation: Double?, relativeCost: Double = 1, evidence: VivoKineticEvidence) {
        self.identifier = identifier; self.templateCaseIdentifier = templateCaseIdentifier; self.exposureKnots = exposureKnots
        self.sampleTimeSeconds = sampleTimeSeconds; self.observable = observable
        self.reportedStandardDeviation = reportedStandardDeviation; self.relativeCost = relativeCost; self.evidence = evidence
    }
}

public struct VivoTargetDesignPolicy: Codable, Equatable, Sendable {
    public let simulatedMeasurements: Int
    public let maximumForwardEvaluations: Int
    public let maximumDensityEvaluations: Int
    public let seed: UInt64
    public init(simulatedMeasurements: Int = 512, maximumForwardEvaluations: Int = 16384,
                maximumDensityEvaluations: Int = 16_777_216, seed: UInt64 = 0x4e5644455349474e) {
        self.simulatedMeasurements = simulatedMeasurements; self.maximumForwardEvaluations = maximumForwardEvaluations
        self.maximumDensityEvaluations = maximumDensityEvaluations; self.seed = seed
    }
}

public struct VivoTargetDesignScore: Codable, Equatable, Sendable {
    public let identifier: String
    public let expectedInformationNats: Double
    public let conditionalMonteCarloStandardError: Double
    public let informationPerRelativeCost: Double
    public let firstHalfInformationNats: Double
    public let secondHalfInformationNats: Double
    public let meanPredictedOccupancy: Double
    public let minimumMeasurementSD: Double
    public let maximumMeasurementSD: Double
}

public struct VivoTargetDesignReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let posteriorFingerprint: VivoFingerprint
    public let candidates: [VivoTargetDesignCandidate]
    public let policy: VivoTargetDesignPolicy
    public let ranking: [VivoTargetDesignScore]
    public let limitations: [String]
}

/// Expected information gain for one future Gaussian measurement, conditional
/// on the empirical JOINT posterior. Scores use E log p(y|theta)/p(y) and retain
/// finite Monte Carlo error, including a potentially negative sample estimate.
/// This does not select a treatment or automate any laboratory experiment.
public enum VivoTargetExperimentDesigner {
    public static func rank(record: VivoTargetPosteriorRecord, candidates: [VivoTargetDesignCandidate],
                            policy: VivoTargetDesignPolicy = .init()) throws -> VivoTargetDesignReport {
        try record.validate(requireComplete: true)
        let prepared = try VivoPreparedTargetPosterior(record.problem)
        guard let checkpoint = record.posterior.checkpoint,
              !candidates.isEmpty, candidates.count <= 256,
              Set(candidates.map(\.identifier)).count == candidates.count,
              (64...65536).contains(policy.simulatedMeasurements),
              (1...1_048_576).contains(policy.maximumForwardEvaluations),
              (1...268_435_456).contains(policy.maximumDensityEvaluations) else {
            throw VivoPosteriorError.invalid("experiment-design policy or candidates")
        }
        let n = checkpoint.particles.count
        let forwardCount = n.multipliedReportingOverflow(by: candidates.count)
        let densityCount = forwardCount.partialValue.multipliedReportingOverflow(by: policy.simulatedMeasurements)
        guard !forwardCount.overflow, !densityCount.overflow, forwardCount.partialValue <= policy.maximumForwardEvaluations,
              densityCount.partialValue <= policy.maximumDensityEvaluations else {
            throw VivoPosteriorError.budget("experiment-design forward or mixture work")
        }
        let values = try checkpoint.particles.map { try prepared.physicalValues($0) }
        var scores: [VivoTargetDesignScore] = [], totalKnots = 0
        for candidate in candidates {
            try Task.checkCancellation(); try candidate.evidence.validate(origin: .assumed)
            guard !candidate.identifier.isEmpty, candidate.identifier.utf8.count <= 1024,
                  candidate.sampleTimeSeconds.isFinite, candidate.sampleTimeSeconds >= 0,
                  candidate.relativeCost.isFinite, candidate.relativeCost > 0,
                  let template = record.problem.study.cases.first(where: { $0.identifier == candidate.templateCaseIdentifier }) else {
                throw VivoPosteriorError.invalid("design candidate context, cost or time")
            }
            totalKnots += candidate.exposureKnots.count
            guard totalKnots <= 262_144 else { throw VivoPosteriorError.budget("design exposure table capacity") }
            try VivoUnboundExposureTrace(context: template.experiment.kinetics.context, knots: candidate.exposureKnots,
                origin: .assumed, evidence: candidate.evidence).validate(for: template.experiment.kinetics)
            var means: [Double] = [], scales: [Double] = [], latentMean = 0.0
            for value in values {
                try Task.checkCancellation()
                let original = try prepared.experiment(for: template, values: value)
                var exposureScale = 1.0
                for (binding, x) in zip(record.problem.bindings, value)
                    where binding.field == .exposureScale && binding.caseIdentifiers.contains(template.identifier) { exposureScale = x }
                let knots = candidate.exposureKnots.map {
                    VivoExposureKnot(timeSeconds: $0.timeSeconds, unboundDrugM: $0.unboundDrugM * exposureScale)
                }
                let experiment = VivoTargetEngagementExperiment(kinetics: original.kinetics,
                    exposure: .init(context: original.kinetics.context, knots: knots, origin: .assumed, evidence: candidate.evidence),
                    initial: original.initial, sampleTimesSeconds: [candidate.sampleTimeSeconds])
                let output = try VivoTargetEngagementReference.run(experiment, numerics: record.problem.numerics)
                guard let sample = output.samples.first else { throw VivoPosteriorError.numerical("missing design prediction") }
                let noise = try record.problem.assayNoise(caseIdentifier: template.identifier, values: value)
                let mean = candidate.observable.value(sample)
                means.append(mean + noise.bias); latentMean += mean / Double(n)
                scales.append(try noise.standardDeviation(reported: candidate.reportedStandardDeviation))
            }
            let identity = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(candidate))
            var candidateSeed = policy.seed
            for byte in identity.bytes { candidateSeed = (candidateSeed ^ UInt64(byte)) &* 0x100000001b3 }
            var random = VivoSplitMix64(state: candidateSeed)
            let estimate = try information(means: means, scales: scales, draws: policy.simulatedMeasurements, random: &random)
            let score = estimate.mean / candidate.relativeCost
            guard score.isFinite else { throw VivoPosteriorError.numerical("design cost normalization overflow") }
            scores.append(.init(identifier: candidate.identifier, expectedInformationNats: estimate.mean,
                conditionalMonteCarloStandardError: estimate.se, informationPerRelativeCost: score,
                firstHalfInformationNats: estimate.first, secondHalfInformationNats: estimate.second,
                meanPredictedOccupancy: latentMean, minimumMeasurementSD: scales.min()!, maximumMeasurementSD: scales.max()!))
        }
        scores.sort { $0.informationPerRelativeCost == $1.informationPerRelativeCost ? $0.identifier < $1.identifier : $0.informationPerRelativeCost > $1.informationPerRelativeCost }
        return try .init(schemaVersion: 1, posteriorFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(record)),
            candidates: candidates, policy: policy, ranking: scores, limitations: [
                "Ranks single measurements for information about the full fitted joint parameter vector, including assay nuisance parameters.",
                "Monte Carlo error is conditional on the finite posterior population; it excludes posterior approximation and missing biological mechanisms.",
                "Assay noise and fitted bias transfer only through the explicitly selected template case. No undisclosed context transfer is inferred.",
                "No candidate outcome is used. Candidates are not treatment recommendations, clinical benefits or automatically scheduled experiments.",
                "All posterior particles must succeed; the designer aborts rather than rank a survivor-conditioned candidate.",
                "Serial measurement correlation is irrelevant for one future scalar; joint multitime experimental design is not implemented.",
                "Finite estimates may be negative. Close rankings relative to their Monte Carlo error need more computation, not a claim of optimality."
            ])
    }

    internal static func information(means: [Double], scales: [Double], draws: Int,
                                     random: inout VivoSplitMix64) throws -> (mean: Double, se: Double, first: Double, second: Double) {
        guard !means.isEmpty, means.count == scales.count, draws >= 2,
              means.allSatisfy(\.isFinite), scales.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoPosteriorError.invalid("design information mixture")
        }
        var mean = 0.0, m2 = 0.0, first = 0.0, second = 0.0
        let half = draws / 2, count = UInt64(means.count), threshold = (0 &- count) % count
        for draw in 0..<draws {
            try Task.checkCancellation()
            var word = random.next()
            while word < threshold { word = random.next() }
            let component = Int(word % count), z = random.normal(), y = means[component] + scales[component] * z
            let logs = try means.indices.map { try VivoGaussianAssay.logContribution(mean: means[$0], sd: scales[$0], value: y) }
            let maximum = logs.max()!
            let marginal = maximum + log(logs.reduce(0) { $0 + exp($1 - maximum) / Double(means.count) })
            let value = logs[component] - marginal
            guard value.isFinite else { throw VivoPosteriorError.numerical("design information overflow") }
            let delta = value - mean; mean += delta / Double(draw + 1); m2 += delta * (value - mean)
            if draw < half { first += value / Double(half) } else { second += value / Double(draws - half) }
        }
        return (mean, sqrt(max(0, m2) / Double(draws - 1) / Double(draws)), first, second)
    }
}
