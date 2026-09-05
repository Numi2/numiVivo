import Foundation

public enum VivoPosteriorError: Error, LocalizedError, Sendable {
    case invalid(String), numerical(String), budget(String), busy
    public var errorDescription: String? {
        switch self {
        case .invalid(let message): return "Invalid posterior request: \(message)"
        case .numerical(let message): return "Posterior numerical failure: \(message)"
        case .budget(let message): return "Posterior work limit: \(message)"
        case .busy: return "Posterior sampler already has an operation in flight"
        }
    }
}

/// Probability distributions, not merely optimization coordinate choices.
/// Mutation targets likelihood(x(u))^beta in uniform prior coordinates u.
public enum VivoBoundedPosteriorPrior: String, Codable, Sendable {
    case uniformPhysical, uniformLogPhysical
}

public struct VivoPosteriorParameter: Codable, Equatable, Sendable {
    public let identifier: String
    public let unit: String
    public let lower: Double
    public let upper: Double
    public let prior: VivoBoundedPosteriorPrior
    public init(identifier: String, unit: String, lower: Double, upper: Double,
                prior: VivoBoundedPosteriorPrior = .uniformLogPhysical) {
        self.identifier = identifier; self.unit = unit; self.lower = lower
        self.upper = upper; self.prior = prior
    }
    public func validate() throws {
        guard !identifier.isEmpty, identifier.utf8.count <= 1024, !unit.isEmpty, unit.utf8.count <= 128,
              lower.isFinite, upper.isFinite, lower < upper, (upper - lower).isFinite,
              abs(lower) <= 1e150, abs(upper) <= 1e150,
              prior != .uniformLogPhysical || lower > 0 else {
            throw VivoPosteriorError.invalid("parameter identity, prior support or numerical range")
        }
    }
    public func value(at coordinate: Double) throws -> Double {
        guard coordinate.isFinite, coordinate >= 0, coordinate <= 1 else {
            throw VivoPosteriorError.invalid("prior coordinate outside [0,1]")
        }
        if coordinate == 0 { return lower }
        if coordinate == 1 { return upper }
        let value: Double
        switch prior {
        case .uniformPhysical: value = lower + coordinate * (upper - lower)
        case .uniformLogPhysical: value = exp(log(lower) + coordinate * (log(upper) - log(lower)))
        }
        guard value.isFinite, value >= lower, value <= upper else {
            throw VivoPosteriorError.numerical("prior transform overflow or loss of support")
        }
        return value
    }
}

public struct VivoSMCConfiguration: Codable, Equatable, Sendable {
    public let particleCount: Int
    public let maximumStages: Int
    public let mutationSweeps: Int
    public let targetESSFraction: Double
    public let proposalScale: Double
    public let covarianceRegularization: Double
    public let independentPriorProposalProbability: Double
    public let minimumTemperatureIncrement: Double
    public let maximumLikelihoodEvaluations: Int
    public let seed: UInt64
    public init(particleCount: Int = 256, maximumStages: Int = 128, mutationSweeps: Int = 6,
                targetESSFraction: Double = 0.8, proposalScale: Double = 1,
                covarianceRegularization: Double = 1e-6,
                independentPriorProposalProbability: Double = 0.1,
                minimumTemperatureIncrement: Double = 1e-12,
                maximumLikelihoodEvaluations: Int = 500_000, seed: UInt64 = 0x4e56504f53544552) {
        self.particleCount = particleCount; self.maximumStages = maximumStages; self.mutationSweeps = mutationSweeps
        self.targetESSFraction = targetESSFraction; self.proposalScale = proposalScale
        self.covarianceRegularization = covarianceRegularization
        self.independentPriorProposalProbability = independentPriorProposalProbability
        self.minimumTemperatureIncrement = minimumTemperatureIncrement
        self.maximumLikelihoodEvaluations = maximumLikelihoodEvaluations; self.seed = seed
    }
    public func validate(dimension: Int) throws {
        guard (1...32).contains(dimension), (32...8192).contains(particleCount), particleCount >= 4 * dimension,
              (1...1024).contains(maximumStages), (1...64).contains(mutationSweeps),
              targetESSFraction.isFinite, targetESSFraction >= 0.1, targetESSFraction < 1,
              proposalScale.isFinite, proposalScale > 0, proposalScale <= 10,
              covarianceRegularization.isFinite, covarianceRegularization >= 1e-12, covarianceRegularization <= 0.1,
              independentPriorProposalProbability.isFinite,
              (0...1).contains(independentPriorProposalProbability),
              minimumTemperatureIncrement.isFinite, minimumTemperatureIncrement >= 1e-15,
              minimumTemperatureIncrement <= 1e-3,
              maximumLikelihoodEvaluations >= particleCount, maximumLikelihoodEvaluations <= 10_000_000 else {
            throw VivoPosteriorError.invalid("SMC dimensions, policy or resource limits")
        }
    }
}

public struct VivoPosteriorScreeningPolicy: Codable, Equatable, Sendable {
    /// Binds a deterministic approximation, its code, precision, settings and
    /// device context. This is not the authoritative likelihood identity.
    public let fingerprint: VivoFingerprint
    public let maximumEvaluations: Int
    public init(fingerprint: VivoFingerprint, maximumEvaluations: Int = 2_000_000) {
        self.fingerprint = fingerprint; self.maximumEvaluations = maximumEvaluations
    }
    public func validate() throws {
        guard (1...20_000_000).contains(maximumEvaluations) else {
            throw VivoPosteriorError.invalid("screening evaluation budget")
        }
    }
}

public struct VivoPosteriorPlan: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    /// Training inputs, observation model and authoritative forward numerics.
    public let likelihoodFingerprint: VivoFingerprint
    public let parameters: [VivoPosteriorParameter]
    public let configuration: VivoSMCConfiguration
    /// Omitted in legacy plans: their canonical bytes and hashes remain intact.
    public let screening: VivoPosteriorScreeningPolicy?
    public init(likelihoodFingerprint: VivoFingerprint, parameters: [VivoPosteriorParameter],
                configuration: VivoSMCConfiguration = .init(), screening: VivoPosteriorScreeningPolicy? = nil) {
        schemaVersion = 1; self.likelihoodFingerprint = likelihoodFingerprint
        self.parameters = parameters; self.configuration = configuration; self.screening = screening
    }
    public func withScreening(_ value: VivoPosteriorScreeningPolicy?) -> Self {
        .init(likelihoodFingerprint: likelihoodFingerprint, parameters: parameters,
              configuration: configuration, screening: value)
    }
    public func validate() throws {
        guard schemaVersion == 1, Set(parameters.map(\.identifier)).count == parameters.count else {
            throw VivoPosteriorError.invalid("posterior schema or duplicate parameter")
        }
        try configuration.validate(dimension: parameters.count); try screening?.validate()
        for parameter in parameters { try parameter.validate() }
    }
    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoPosteriorCandidate: Codable, Equatable, Sendable {
    public let ordinal: UInt64
    public let coordinates: [Double]
    public let values: [Double]
    public init(ordinal: UInt64, coordinates: [Double], parameters: [VivoPosteriorParameter]) throws {
        guard coordinates.count == parameters.count else { throw VivoPosteriorError.invalid("candidate dimension") }
        self.ordinal = ordinal; self.coordinates = coordinates
        self.values = try zip(parameters, coordinates).map { try $0.value(at: $1) }
    }
}

public struct VivoPosteriorEvaluation: Codable, Equatable, Sendable {
    public let candidate: VivoPosteriorCandidate
    /// Nil is failure, never zero likelihood. Both exact and screening evaluators
    /// must provide finite deterministic values; no stochastic estimate is assumed.
    public let logLikelihood: Double?
    public let failure: String?
    public init(candidate: VivoPosteriorCandidate, logLikelihood: Double) {
        self.candidate = candidate; self.logLikelihood = logLikelihood; failure = nil
    }
    public init(candidate: VivoPosteriorCandidate, failure: String) {
        self.candidate = candidate; logLikelihood = nil; self.failure = failure
    }
}
public typealias VivoPosteriorBatchEvaluator = @Sendable ([VivoPosteriorCandidate]) async throws -> [VivoPosteriorEvaluation]
public typealias VivoPosteriorCheckpointSink = @Sendable (VivoPosteriorCheckpoint) async throws -> Void

public struct VivoPosteriorScreen: Sendable {
    public let policy: VivoPosteriorScreeningPolicy
    /// Must depend only on coordinates/values and immutable context, not batch
    /// companions, ordinal, call order or previous evaluations. No adaptive refits.
    public let evaluate: VivoPosteriorBatchEvaluator
    public init(policy: VivoPosteriorScreeningPolicy, evaluate: @escaping VivoPosteriorBatchEvaluator) {
        self.policy = policy; self.evaluate = evaluate
    }
}

public struct VivoPosteriorParticle: Codable, Equatable, Sendable {
    public let coordinates: [Double]
    public let logLikelihood: Double
    public let initialAncestor: Int
}

public struct VivoSMCStage: Codable, Equatable, Sendable {
    public let index: Int
    public let betaBefore: Double
    public let betaAfter: Double
    public let preResamplingESS: Double
    public let likelihoodEvaluations: Int
    public let proposedMoves: Int
    public let acceptedMoves: Int
    public let outOfPriorMoves: Int
    public let distinctInitialAncestors: Int
    public let screenedOutMoves: Int?
    public let screeningEvaluations: Int?
    public init(index: Int, betaBefore: Double, betaAfter: Double, preResamplingESS: Double,
                likelihoodEvaluations: Int, proposedMoves: Int, acceptedMoves: Int, outOfPriorMoves: Int,
                distinctInitialAncestors: Int, screenedOutMoves: Int? = nil, screeningEvaluations: Int? = nil) {
        self.index = index; self.betaBefore = betaBefore; self.betaAfter = betaAfter
        self.preResamplingESS = preResamplingESS; self.likelihoodEvaluations = likelihoodEvaluations
        self.proposedMoves = proposedMoves; self.acceptedMoves = acceptedMoves; self.outOfPriorMoves = outOfPriorMoves
        self.distinctInitialAncestors = distinctInitialAncestors; self.screenedOutMoves = screenedOutMoves
        self.screeningEvaluations = screeningEvaluations
    }
}

/// Fully committed tempering boundary. Screening values are recomputed for the
/// resampled population at the next stage; no partially filled GPU cache is state.
public struct VivoPosteriorCheckpoint: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let planFingerprint: VivoFingerprint
    public let beta: Double
    public let particles: [VivoPosteriorParticle]
    public let randomState: UInt64
    public let nextCandidateOrdinal: UInt64
    public let likelihoodEvaluationCount: Int
    public let stages: [VivoSMCStage]
    public func validate(for plan: VivoPosteriorPlan) throws {
        try plan.validate()
        let c = plan.configuration, d = plan.parameters.count
        guard schemaVersion == 1, planFingerprint == (try plan.fingerprint()),
              beta.isFinite, (0...1).contains(beta), particles.count == c.particleCount,
              stages.count <= c.maximumStages, likelihoodEvaluationCount >= c.particleCount,
              likelihoodEvaluationCount <= c.maximumLikelihoodEvaluations,
              nextCandidateOrdinal == UInt64(likelihoodEvaluationCount) else {
            throw VivoPosteriorError.invalid("checkpoint identity, population, clock or evaluation count")
        }
        var previous = 0.0, evaluations = c.particleCount, screens = 0
        for (index, stage) in stages.enumerated() {
            guard stage.index == index, stage.betaBefore == previous, stage.betaAfter.isFinite,
                  stage.betaAfter > previous, stage.betaAfter <= 1,
                  stage.preResamplingESS.isFinite, stage.preResamplingESS >= 1 - 1e-10,
                  stage.preResamplingESS <= Double(c.particleCount) + 1e-8,
                  stage.proposedMoves == c.particleCount * c.mutationSweeps,
                  stage.outOfPriorMoves >= 0, stage.outOfPriorMoves <= stage.proposedMoves,
                  stage.acceptedMoves >= 0, stage.acceptedMoves <= stage.likelihoodEvaluations,
                  (1...c.particleCount).contains(stage.distinctInitialAncestors) else {
                throw VivoPosteriorError.invalid("checkpoint stage ledger")
            }
            let valid = stage.proposedMoves - stage.outOfPriorMoves
            if let policy = plan.screening {
                guard let rejected = stage.screenedOutMoves, let count = stage.screeningEvaluations,
                      rejected >= 0, rejected <= valid, count == c.particleCount + valid,
                      stage.likelihoodEvaluations == valid - rejected else {
                    throw VivoPosteriorError.invalid("delayed-acceptance evaluation ledger")
                }
                screens += count
                guard screens <= policy.maximumEvaluations else { throw VivoPosteriorError.invalid("screening budget exceeded") }
            } else {
                guard stage.screenedOutMoves == nil, stage.screeningEvaluations == nil,
                      stage.likelihoodEvaluations == valid else { throw VivoPosteriorError.invalid("unexpected screening ledger") }
            }
            evaluations += stage.likelihoodEvaluations; previous = stage.betaAfter
        }
        guard previous == beta, evaluations == likelihoodEvaluationCount,
              beta != 1 || !stages.isEmpty else { throw VivoPosteriorError.invalid("checkpoint terminal stage") }
        for p in particles {
            guard p.coordinates.count == d, p.coordinates.allSatisfy({ $0.isFinite && (0...1).contains($0) }),
                  p.logLikelihood.isFinite, abs(p.logLikelihood) <= 1e250,
                  (0..<c.particleCount).contains(p.initialAncestor) else {
                throw VivoPosteriorError.invalid("checkpoint particle support")
            }
        }
        if let last = stages.last {
            guard Set(particles.map(\.initialAncestor)).count == last.distinctInitialAncestors else {
                throw VivoPosteriorError.invalid("checkpoint ancestry ledger")
            }
        }
    }
}

public struct VivoPosteriorFailure: Codable, Equatable, Sendable {
    public let message: String
    public let attemptedBeta: Double?
    public let failedEvaluations: [VivoPosteriorEvaluation]
}

public struct VivoPosteriorRun: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let plan: VivoPosteriorPlan
    public let checkpoint: VivoPosteriorCheckpoint?
    public let failure: VivoPosteriorFailure?
    public let limitations: [String]
    public var completed: Bool { schemaVersion == 1 && failure == nil && checkpoint?.beta == 1 }
    public func validate(requireComplete: Bool = false) throws {
        guard schemaVersion == 1 else { throw VivoPosteriorError.invalid("posterior run schema") }
        try plan.validate(); try checkpoint?.validate(for: plan)
        if requireComplete, !completed { throw VivoPosteriorError.invalid("tempered or failed population is not a posterior result") }
        guard checkpoint != nil || failure != nil else { throw VivoPosteriorError.invalid("empty run") }
    }
}
