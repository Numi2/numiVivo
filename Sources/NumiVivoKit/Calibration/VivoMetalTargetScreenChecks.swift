import Foundation

public struct VivoMetalScreenCheckPolicy: Codable, Equatable, Sendable {
    public let candidateCount: Int
    public let seed: UInt64
    public init(candidateCount: Int = 8, seed: UInt64 = 0x4e5653435245454e) {
        self.candidateCount = candidateCount; self.seed = seed
    }
}
public struct VivoMetalScreenComparisonPoint: Codable, Equatable, Sendable {
    public let candidate: VivoPosteriorCandidate
    public let metalScreen: Double
    public let diagonalFP64Screen: Double
    public let authoritativeFP64Likelihood: Double
}
public struct VivoMetalScreenCheckReport: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let screenDescription: VivoMetalTargetScreenDescription
    public let policy: VivoMetalScreenCheckPolicy
    public let points: [VivoMetalScreenComparisonPoint]
    public let repeatedBatchIdentical: Bool
    public let permutationInvariant: Bool
    public let partitionInvariant: Bool
    public let ordinalInvariant: Bool
    public let maximumDifferenceFromDiagonalFP64: Double
    public let maximumDifferenceFromAuthority: Double
    public let limitations: [String]
    public var deterministicProbeChecksPassed: Bool {
        repeatedBatchIdentical && permutationInvariant && partitionInvariant && ordinalInvariant
    }
}

/// Probe checks, not a global numerical certificate or performance benchmark.
/// No failed candidate is silently removed. These checks execute actual GPU work
/// on the user's selected device and must never be labelled passed before running.
public enum VivoMetalTargetScreenChecks {
    public static func run(problem: VivoTargetPosteriorProblem, screen: VivoMetalTargetLikelihoodScreen,
                           policy: VivoMetalScreenCheckPolicy = .init()) async throws -> VivoMetalScreenCheckReport {
        let prepared = try VivoPreparedTargetPosterior(problem)
        guard (2...256).contains(policy.candidateCount), policy.candidateCount <= problem.sampler.particleCount,
              screen.description.authoritativeLikelihoodFingerprint == prepared.plan.likelihoodFingerprint else {
            throw VivoPosteriorError.invalid("screen check problem identity or probe capacity")
        }
        var random = VivoSplitMix64(state: policy.seed)
        let candidates = try (0..<policy.candidateCount).map { index in
            try VivoPosteriorCandidate(ordinal: UInt64(index), coordinates: prepared.plan.parameters.indices.map { _ in
                index == 0 ? 0.5 : VivoPosteriorNumerics.openUnit(&random)
            }, parameters: prepared.plan.parameters)
        }
        func values(_ evaluations: [VivoPosteriorEvaluation]) throws -> [UInt64: Double] {
            guard Set(evaluations.map { $0.candidate.ordinal }).count == evaluations.count else {
                throw VivoPosteriorError.invalid("duplicate probe identity")
            }
            return try Dictionary(uniqueKeysWithValues: evaluations.map { value in
                guard value.failure == nil, let score = value.logLikelihood, score.isFinite else {
                    throw VivoPosteriorError.numerical("screen check evaluation failed")
                }
                return (value.candidate.ordinal, score)
            })
        }
        let baseline = try values(await screen.evaluate(candidates))
        let repeated = try values(await screen.evaluate(candidates))
        let reversed = try values(await screen.evaluate(Array(candidates.reversed())))
        var partitioned: [UInt64: Double] = [:]
        for candidate in candidates {
            let single = try values(await screen.evaluate([candidate]))
            partitioned[candidate.ordinal] = single[candidate.ordinal]
        }
        let renamed = try candidates.map { try VivoPosteriorCandidate(ordinal: $0.ordinal + 1_000_000,
            coordinates: $0.coordinates, parameters: prepared.plan.parameters) }
        let renamedValues = try values(await screen.evaluate(renamed))
        var points: [VivoMetalScreenComparisonPoint] = []
        var ordinalInvariant = true
        for candidate in candidates {
            guard let gpu = baseline[candidate.ordinal], let renamedValue = renamedValues[candidate.ordinal + 1_000_000] else {
                throw VivoPosteriorError.invalid("missing probe")
            }
            ordinalInvariant = ordinalInvariant && gpu.bitPattern == renamedValue.bitPattern
            let exact = try prepared.logLikelihood(values: candidate.values)
            let diagonal = try diagonalReference(prepared: prepared, candidate: candidate)
            points.append(.init(candidate: candidate, metalScreen: gpu,
                diagonalFP64Screen: diagonal, authoritativeFP64Likelihood: exact))
        }
        func same(_ a: [UInt64: Double], _ b: [UInt64: Double]) -> Bool {
            a.count == b.count && a.allSatisfy { key, value in b[key]?.bitPattern == value.bitPattern }
        }
        return .init(schemaVersion: 1, screenDescription: screen.description, policy: policy, points: points,
            repeatedBatchIdentical: same(baseline, repeated), permutationInvariant: same(baseline, reversed),
            partitionInvariant: same(baseline, partitioned), ordinalInvariant: ordinalInvariant,
            maximumDifferenceFromDiagonalFP64: points.map { abs($0.metalScreen - $0.diagonalFP64Screen) }.max() ?? 0,
            maximumDifferenceFromAuthority: points.map { abs($0.metalScreen - $0.authoritativeFP64Likelihood) }.max() ?? 0,
            limitations: [
                "Finite deterministic probes only; not a proof across all prior coordinates or a speed measurement.",
                "Diagonal FP64 comparison includes F1 integration and FP32 scoring differences; no universal error bound is inferred.",
                "The screen intentionally omits within-case assay correlation. That correlation remains in authoritative acceptance correction.",
                "Different batches, partial-capacity padding and candidate ordinals must not change a candidate's score.",
                "All input evidence classes are unchanged. Numerical comparison is not biological validation.",
                "Probe evaluations are setup work, separate from the SMC stage screening-evaluation ledger."
            ])
    }

    private static func diagonalReference(prepared: VivoPreparedTargetPosterior,
                                          candidate: VivoPosteriorCandidate) throws -> Double {
        var sum = 0.0
        for item in prepared.calibrationCases {
            let experiment = try prepared.experiment(for: item, values: candidate.values)
            let result = try VivoTargetEngagementReference.run(experiment, numerics: prepared.problem.numerics)
            let samples = Dictionary(uniqueKeysWithValues: result.samples.map { ($0.timeSeconds, $0) })
            let predictions = try item.observations.map { observation -> Double in
                guard let sample = samples[observation.timeSeconds] else { throw VivoPosteriorError.invalid("missing reference probe time") }
                return observation.observable.value(sample)
            }
            let actualNoise = try prepared.problem.assayNoise(caseIdentifier: item.identifier, values: candidate.values)
            let noise = VivoGaussianAssayNoise(scale: actualNoise.scale,
                additionalStandardDeviation: actualNoise.additionalStandardDeviation, bias: actualNoise.bias)
            let observations = item.observations.map {
                VivoGaussianAssayDatum(timeSeconds: $0.timeSeconds, value: $0.value,
                    reportedStandardDeviation: $0.standardDeviation, support: .exact)
            }
            sum += try VivoGaussianAssay.logLikelihood(predictions: predictions, observations: observations, noise: noise)
        }
        guard sum.isFinite else { throw VivoPosteriorError.numerical("reference probe likelihood overflow") }
        return sum
    }
}
