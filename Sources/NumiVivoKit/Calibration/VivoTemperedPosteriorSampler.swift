import Foundation

/// Adaptive tempered sequential Monte Carlo, not an elite optimizer. Targets
/// p_beta(u) proportional to L(x(u))^beta for uniform prior coordinates. Each
/// stage weights, resamples and applies invariant Metropolis kernels. The full
/// proposal covariance is frozen before that stage's mutations; proposals outside
/// prior support are rejected, never reflected with an invalid correlated kernel.
/// Method family: Del Moral, Doucet & Jasra (2006), doi:10.1111/j.1467-9868.2006.00553.x.
public actor VivoTemperedPosteriorSampler {
    public nonisolated let plan: VivoPosteriorPlan
    private var committed: VivoPosteriorCheckpoint?
    private var inFlight = false

    private struct EvaluationFailure: Error, Sendable {
        let message: String
        let evaluations: [VivoPosteriorEvaluation]
    }

    public init(plan: VivoPosteriorPlan, checkpoint: VivoPosteriorCheckpoint? = nil) throws {
        try plan.validate(); try checkpoint?.validate(for: plan)
        self.plan = plan; committed = checkpoint
    }

    /// May be read during a run: this is always the last COMPLETED stage, not the
    /// scratch RNG state or a partial population. The returned value is immutable.
    public func checkpoint() -> VivoPosteriorCheckpoint? { committed }

    public func run(evaluate: VivoPosteriorBatchEvaluator,
                    progress: VivoPosteriorCheckpointSink? = nil) async throws -> VivoPosteriorRun {
        guard !inFlight else { throw VivoPosteriorError.busy }
        inFlight = true
        defer { inFlight = false }
        var attemptedBeta: Double?
        do {
            try Task.checkCancellation()
            if committed == nil {
                let initial = try await initialize(evaluate: evaluate)
                try Task.checkCancellation()
                committed = initial
                if let progress { try await progress(initial) }
            }
            while let current = committed, current.beta < 1 {
                try Task.checkCancellation()
                guard current.stages.count < plan.configuration.maximumStages else {
                    throw VivoPosteriorError.budget("maximum tempering stages reached before beta=1")
                }
                let nextBeta = try VivoPosteriorNumerics.nextTemperature(
                    logLikelihoods: current.particles.map(\.logLikelihood), beta: current.beta,
                    configuration: plan.configuration)
                attemptedBeta = nextBeta
                let next = try await advance(current, beta: nextBeta, evaluate: evaluate)
                try Task.checkCancellation()
                try next.validate(for: plan)
                committed = next
                if let progress { try await progress(next) }
            }
            return result(failure: nil)
        } catch is CancellationError {
            // Caller can persist checkpoint(); no partial stage is published.
            throw CancellationError()
        } catch let failure as EvaluationFailure {
            return result(failure: .init(message: failure.message, attemptedBeta: attemptedBeta,
                                        failedEvaluations: failure.evaluations))
        } catch {
            return result(failure: .init(message: String(error.localizedDescription.prefix(8192)),
                                        attemptedBeta: attemptedBeta, failedEvaluations: []))
        }
    }

    private func initialize(evaluate: VivoPosteriorBatchEvaluator) async throws -> VivoPosteriorCheckpoint {
        var random = VivoSplitMix64(state: plan.configuration.seed)
        let n = plan.configuration.particleCount, d = plan.parameters.count
        let candidates = try (0..<n).map { index in
            try VivoPosteriorCandidate(ordinal: UInt64(index),
                coordinates: (0..<d).map { _ in VivoPosteriorNumerics.openUnit(&random) }, parameters: plan.parameters)
        }
        let likelihoods = try await checkedEvaluation(candidates, evaluate: evaluate)
        let particles = candidates.indices.map {
            VivoPosteriorParticle(coordinates: candidates[$0].coordinates, logLikelihood: likelihoods[$0], initialAncestor: $0)
        }
        return try .init(schemaVersion: 1, planFingerprint: plan.fingerprint(), beta: 0,
                         particles: particles, randomState: random.state, nextCandidateOrdinal: UInt64(n),
                         likelihoodEvaluationCount: n, stages: [])
    }

    private func advance(_ current: VivoPosteriorCheckpoint, beta: Double,
                         evaluate: VivoPosteriorBatchEvaluator) async throws -> VivoPosteriorCheckpoint {
        let c = plan.configuration, n = c.particleCount, d = plan.parameters.count
        var random = VivoSplitMix64(state: current.randomState)
        let (weights, ess) = try VivoPosteriorNumerics.normalizedWeights(
            logLikelihoods: current.particles.map(\.logLikelihood), delta: beta - current.beta)
        let factor = try VivoPosteriorNumerics.proposalCholesky(particles: current.particles,
                                                               weights: weights, configuration: c)
        var particles = VivoPosteriorNumerics.resample(current.particles, weights: weights, random: &random)
        var ordinal = current.nextCandidateOrdinal, evaluations = 0, accepted = 0, outside = 0
        for _ in 0..<c.mutationSweeps {
            try Task.checkCancellation()
            var candidates: [VivoPosteriorCandidate] = [], destinations: [Int] = [], logUniforms: [Double] = []
            candidates.reserveCapacity(n); destinations.reserveCapacity(n); logUniforms.reserveCapacity(n)
            for index in 0..<n {
                let coordinates: [Double]
                if VivoPosteriorNumerics.openUnit(&random) < c.independentPriorProposalProbability {
                    coordinates = (0..<d).map { _ in VivoPosteriorNumerics.openUnit(&random) }
                } else {
                    let normal = (0..<d).map { _ in random.normal() }
                    coordinates = (0..<d).map { row in
                        var value = particles[index].coordinates[row]
                        for column in 0...row { value += factor[row * d + column] * normal[column] }
                        return value
                    }
                }
                guard coordinates.allSatisfy({ $0.isFinite && (0...1).contains($0) }) else { outside += 1; continue }
                guard ordinal < UInt64.max else { throw VivoPosteriorError.budget("candidate ordinal overflow") }
                candidates.append(try .init(ordinal: ordinal, coordinates: coordinates, parameters: plan.parameters))
                destinations.append(index); logUniforms.append(log(VivoPosteriorNumerics.openUnit(&random)))
                ordinal += 1
            }
            guard current.likelihoodEvaluationCount + evaluations + candidates.count <= c.maximumLikelihoodEvaluations else {
                throw VivoPosteriorError.budget("likelihood budget would be exceeded; last completed stage retained")
            }
            if candidates.isEmpty { continue }
            let likelihoods = try await checkedEvaluation(candidates, evaluate: evaluate)
            evaluations += candidates.count
            for index in candidates.indices {
                let destination = destinations[index]
                let logRatio = beta * (likelihoods[index] - particles[destination].logLikelihood)
                if logUniforms[index] < min(0, logRatio) {
                    particles[destination] = .init(coordinates: candidates[index].coordinates,
                        logLikelihood: likelihoods[index], initialAncestor: particles[destination].initialAncestor)
                    accepted += 1
                }
            }
        }
        let stage = VivoSMCStage(index: current.stages.count, betaBefore: current.beta, betaAfter: beta,
            preResamplingESS: ess, likelihoodEvaluations: evaluations, proposedMoves: n * c.mutationSweeps,
            acceptedMoves: accepted, outOfPriorMoves: outside,
            distinctInitialAncestors: Set(particles.map(\.initialAncestor)).count)
        return .init(schemaVersion: 1, planFingerprint: current.planFingerprint, beta: beta, particles: particles,
            randomState: random.state, nextCandidateOrdinal: ordinal,
            likelihoodEvaluationCount: current.likelihoodEvaluationCount + evaluations, stages: current.stages + [stage])
    }

    private func checkedEvaluation(_ candidates: [VivoPosteriorCandidate],
                                   evaluate: VivoPosteriorBatchEvaluator) async throws -> [Double] {
        let returned: [VivoPosteriorEvaluation]
        do { returned = try await evaluate(candidates) }
        catch is CancellationError { throw CancellationError() }
        catch {
            let message = String(error.localizedDescription.prefix(4096))
            throw EvaluationFailure(message: "Batch evaluator failed; no support truncation or redraw was applied.",
                                    evaluations: candidates.map { .init(candidate: $0, failure: message) })
        }
        guard returned.count == candidates.count,
              Set(returned.map { $0.candidate.ordinal }).count == candidates.count else {
            throw EvaluationFailure(message: "Evaluator omitted or duplicated candidate identities.", evaluations: [])
        }
        let byOrdinal = Dictionary(uniqueKeysWithValues: returned.map { ($0.candidate.ordinal, $0) })
        var ordered: [Double] = [], failures: [VivoPosteriorEvaluation] = []
        for candidate in candidates {
            guard let result = byOrdinal[candidate.ordinal], result.candidate == candidate else {
                failures.append(.init(candidate: candidate, failure: "candidate values/identity were substituted")); continue
            }
            guard result.failure == nil, let logLikelihood = result.logLikelihood,
                  logLikelihood.isFinite, abs(logLikelihood) <= 1e250 else {
                failures.append(.init(candidate: candidate, failure: String((result.failure ?? "nonfinite/unsupported log likelihood").prefix(4096))))
                continue
            }
            ordered.append(logLikelihood)
        }
        guard failures.isEmpty else {
            throw EvaluationFailure(message: "Likelihood failure stops inference; failed parameter regions were not discarded.",
                                    evaluations: failures)
        }
        return ordered
    }

    private func result(failure: VivoPosteriorFailure?) -> VivoPosteriorRun {
        .init(schemaVersion: 1, plan: plan, checkpoint: committed, failure: failure, limitations: [
            "Finite-particle approximation conditional on declared independent bounded priors and a deterministic likelihood.",
            "Uniform physical and uniform log-physical priors are different distributions; mutation uses their uniform prior coordinates.",
            "All completed-stage particles have equal weights after resampling and mutation; particle count is not an independent-sample ESS.",
            "Pre-resampling ESS measures weight concentration, not posterior convergence or interval coverage.",
            "A completed beta=1 schedule is not proof of adequate exploration; use repeated seeds, particle/sweep sensitivity and independent validation.",
            "Numerical, batch-contract and resource failures stop the run. An incomplete tempered population must not be published as a posterior.",
            "Checkpoints preserve completed stages, RNG state and candidate order. External progress writes are not process-crash-atomic with evaluator side effects.",
            "No marginal-likelihood estimate, clinical validity or global identifiability claim is produced."
        ])
    }
}
