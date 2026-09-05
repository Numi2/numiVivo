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

/// These are probability distributions, not optimization coordinate choices.
/// Independent uniform prior coordinates u in (0,1) map either to uniform x or
/// log-uniform x. Mutation targets likelihood(x(u))^beta in u coordinates, so no
/// Jacobian is omitted: the physical prior is defined by this pushforward.
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

public struct VivoPosteriorPlan: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    /// Must bind the actual training inputs, observation model, forward numerical
    /// policy and backend version. Held-out observations do not belong here.
    public let likelihoodFingerprint: VivoFingerprint
    public let parameters: [VivoPosteriorParameter]
    public let configuration: VivoSMCConfiguration
    public init(likelihoodFingerprint: VivoFingerprint, parameters: [VivoPosteriorParameter],
                configuration: VivoSMCConfiguration = .init()) {
        schemaVersion = 1; self.likelihoodFingerprint = likelihoodFingerprint
        self.parameters = parameters; self.configuration = configuration
    }
    public func validate() throws {
        guard schemaVersion == 1, Set(parameters.map(\.identifier)).count == parameters.count else {
            throw VivoPosteriorError.invalid("posterior schema or duplicate parameter")
        }
        try configuration.validate(dimension: parameters.count)
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
    /// Nil is failure, never zero likelihood. This sampler's current contract is
    /// finite deterministic log likelihood; stochastic estimates need pseudo-marginal semantics.
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
}

/// A fully committed tempering boundary, never a partially mutated population.
/// All particles have equal weight because each completed stage includes
/// resampling followed by fixed-kernel invariant mutation.
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
        var previous = 0.0, evaluations = c.particleCount
        for (index, stage) in stages.enumerated() {
            guard stage.index == index, stage.betaBefore == previous, stage.betaAfter.isFinite,
                  stage.betaAfter > previous, stage.betaAfter <= 1,
                  stage.preResamplingESS.isFinite, stage.preResamplingESS >= 1 - 1e-10,
                  stage.preResamplingESS <= Double(c.particleCount) + 1e-8,
                  stage.proposedMoves == c.particleCount * c.mutationSweeps,
                  stage.outOfPriorMoves >= 0, stage.outOfPriorMoves <= stage.proposedMoves,
                  stage.likelihoodEvaluations == stage.proposedMoves - stage.outOfPriorMoves,
                  stage.acceptedMoves >= 0, stage.acceptedMoves <= stage.likelihoodEvaluations,
                  (1...c.particleCount).contains(stage.distinctInitialAncestors) else {
                throw VivoPosteriorError.invalid("checkpoint stage ledger")
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
