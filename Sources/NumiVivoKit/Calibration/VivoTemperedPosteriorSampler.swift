import Foundation

/// Adaptive tempered SMC in uniform prior coordinates. An optional deterministic
/// screen changes proposal efficiency, not the authoritative target: stage-one
/// acceptance uses beta*(s'-s); stage two uses beta*((l'-l)-(s'-s)). Weights,
/// temperature selection and stored particle likelihoods always use l.
/// Del Moral, Doucet & Jasra (2006), doi:10.1111/j.1467-9868.2006.00553.x;
/// Christen & Fox (2005), JCGS 14(4):795-810, doi:10.1198/106186005X76983.
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
    public func checkpoint() -> VivoPosteriorCheckpoint? { committed }

    public func run(evaluate: VivoPosteriorBatchEvaluator,
                    progress: VivoPosteriorCheckpointSink? = nil,
                    screen: VivoPosteriorScreen? = nil) async throws -> VivoPosteriorRun {
        guard !inFlight else { throw VivoPosteriorError.busy }
        guard screen?.policy == plan.screening else {
            throw VivoPosteriorError.invalid("screen implementation/budget differs from the checkpoint-bound plan")
        }
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
                let next = try await advance(current, beta: nextBeta, evaluate: evaluate, screen: screen)
                try Task.checkCancellation()
                try next.validate(for: plan)
                committed = next
                if let progress { try await progress(next) }
            }
            return result(failure: nil)
        } catch is CancellationError {
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
        let likelihoods = try await checkedEvaluation(candidates, evaluate: evaluate, label: "authoritative")
        let particles = candidates.indices.map {
            VivoPosteriorParticle(coordinates: candidates[$0].coordinates, logLikelihood: likelihoods[$0], initialAncestor: $0)
        }
        return try .init(schemaVersion: 1, planFingerprint: plan.fingerprint(), beta: 0,
                         particles: particles, randomState: random.state, nextCandidateOrdinal: UInt64(n),
                         likelihoodEvaluationCount: n, stages: [])
    }

    private func advance(_ current: VivoPosteriorCheckpoint, beta: Double,
                         evaluate: VivoPosteriorBatchEvaluator, screen: VivoPosteriorScreen?) async throws -> VivoPosteriorCheckpoint {
        let c = plan.configuration, n = c.particleCount, d = plan.parameters.count
        var random = VivoSplitMix64(state: current.randomState)
        let (weights, ess) = try VivoPosteriorNumerics.normalizedWeights(
            logLikelihoods: current.particles.map(\.logLikelihood), delta: beta - current.beta)
        let factor = try VivoPosteriorNumerics.proposalCholesky(particles: current.particles,
                                                               weights: weights, configuration: c)
        var particles = VivoPosteriorNumerics.resample(current.particles, weights: weights, random: &random)
        var ordinal = current.nextCandidateOrdinal, evaluations = 0, accepted = 0, outside = 0
        var screenValues = [Double](repeating: 0, count: n), screens = 0, screenedOut = 0
        let pastScreens = current.stages.reduce(0) { $0 + ($1.screeningEvaluations ?? 0) }
        if let screen {
            guard pastScreens + n <= screen.policy.maximumEvaluations else {
                throw VivoPosteriorError.budget("screening budget before stage initialization")
            }
            let initial = try particles.enumerated().map {
                try VivoPosteriorCandidate(ordinal: UInt64(pastScreens + $0.offset),
                    coordinates: $0.element.coordinates, parameters: plan.parameters)
            }
            screenValues = try await checkedEvaluation(initial, evaluate: screen.evaluate, label: "screening")
            screens = n
        }
        for _ in 0..<c.mutationSweeps {
            try Task.checkCancellation()
            var proposals: [VivoPosteriorCandidate] = [], destinations: [Int] = []
            var firstUniforms: [Double] = [], secondUniforms: [Double] = []
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
                // Screening ordinals have their own monotonically increasing
                // namespace. They are never cache keys without input values.
                let proposalOrdinal = screen == nil ? ordinal + UInt64(proposals.count)
                    : UInt64(pastScreens + screens + proposals.count)
                proposals.append(try .init(ordinal: proposalOrdinal, coordinates: coordinates, parameters: plan.parameters))
                destinations.append(index)
                firstUniforms.append(log(VivoPosteriorNumerics.openUnit(&random)))
                if screen != nil { secondUniforms.append(log(VivoPosteriorNumerics.openUnit(&random))) }
            }
            if proposals.isEmpty { continue }
            var proposedScreens = [Double](repeating: 0, count: proposals.count)
            var eligible = Array(proposals.indices)
            if let screen {
                guard pastScreens + screens + proposals.count <= screen.policy.maximumEvaluations else {
                    throw VivoPosteriorError.budget("screening proposal budget; last completed stage retained")
                }
                proposedScreens = try await checkedEvaluation(proposals, evaluate: screen.evaluate, label: "screening")
                screens += proposals.count
                eligible = proposals.indices.filter { i in
                    firstUniforms[i] < min(0, beta * (proposedScreens[i] - screenValues[destinations[i]]))
                }
                screenedOut += proposals.count - eligible.count
            }
            guard current.likelihoodEvaluationCount + evaluations + eligible.count <= c.maximumLikelihoodEvaluations else {
                throw VivoPosteriorError.budget("authoritative likelihood budget; last completed stage retained")
            }
            if eligible.isEmpty { continue }
            let exact = try eligible.enumerated().map { offset, source in
                try VivoPosteriorCandidate(ordinal: ordinal + UInt64(offset),
                    coordinates: proposals[source].coordinates, parameters: plan.parameters)
            }
            let likelihoods = try await checkedEvaluation(exact, evaluate: evaluate, label: "authoritative")
            ordinal += UInt64(exact.count); evaluations += exact.count
            for (index, source) in eligible.enumerated() {
                let destination = destinations[source]
                let exactDifference = likelihoods[index] - particles[destination].logLikelihood
                let correction = screen == nil ? exactDifference
                    : exactDifference - (proposedScreens[source] - screenValues[destination])
                let logUniform = screen == nil ? firstUniforms[source] : secondUniforms[source]
                if logUniform < min(0, beta * correction) {
                    particles[destination] = .init(coordinates: exact[index].coordinates,
                        logLikelihood: likelihoods[index], initialAncestor: particles[destination].initialAncestor)
                    screenValues[destination] = proposedScreens[source]
                    accepted += 1
                }
            }
        }
        let stage = VivoSMCStage(index: current.stages.count, betaBefore: current.beta, betaAfter: beta,
            preResamplingESS: ess, likelihoodEvaluations: evaluations, proposedMoves: n * c.mutationSweeps,
            acceptedMoves: accepted, outOfPriorMoves: outside,
            distinctInitialAncestors: Set(particles.map(\.initialAncestor)).count,
            screenedOutMoves: screen == nil ? nil : screenedOut, screeningEvaluations: screen == nil ? nil : screens)
        return .init(schemaVersion: 1, planFingerprint: current.planFingerprint, beta: beta, particles: particles,
            randomState: random.state, nextCandidateOrdinal: ordinal,
            likelihoodEvaluationCount: current.likelihoodEvaluationCount + evaluations, stages: current.stages + [stage])
    }

    private func checkedEvaluation(_ candidates: [VivoPosteriorCandidate], evaluate: VivoPosteriorBatchEvaluator,
                                   label: String) async throws -> [Double] {
        let returned: [VivoPosteriorEvaluation]
        do { returned = try await evaluate(candidates) }
        catch is CancellationError { throw CancellationError() }
        catch {
            let message = String(error.localizedDescription.prefix(4096))
            throw EvaluationFailure(message: "\(label) batch failed; no support truncation or redraw applied.",
                                    evaluations: candidates.map { .init(candidate: $0, failure: message) })
        }
        guard returned.count == candidates.count,
              Set(returned.map { $0.candidate.ordinal }).count == candidates.count else {
            throw EvaluationFailure(message: "\(label) evaluator omitted or duplicated identities.", evaluations: [])
        }
        let byOrdinal = Dictionary(uniqueKeysWithValues: returned.map { ($0.candidate.ordinal, $0) })
        var ordered: [Double] = [], failures: [VivoPosteriorEvaluation] = []
        for candidate in candidates {
            guard let result = byOrdinal[candidate.ordinal], result.candidate == candidate else {
                failures.append(.init(candidate: candidate, failure: "candidate values/identity substituted")); continue
            }
            guard result.failure == nil, let logLikelihood = result.logLikelihood,
                  logLikelihood.isFinite, abs(logLikelihood) <= 1e250 else {
                failures.append(.init(candidate: candidate, failure: String((result.failure ?? "nonfinite/unsupported log likelihood").prefix(4096))))
                continue
            }
            ordered.append(logLikelihood)
        }
        guard failures.isEmpty else {
            throw EvaluationFailure(message: "\(label) failure stops inference; failed parameter regions were not discarded.", evaluations: failures)
        }
        return ordered
    }

    private func result(failure: VivoPosteriorFailure?) -> VivoPosteriorRun {
        .init(schemaVersion: 1, plan: plan, checkpoint: committed, failure: failure, limitations: [
            "Finite-particle approximation conditional on declared independent bounded priors and a deterministic authoritative likelihood.",
            "Uniform physical and uniform log-physical priors are different distributions; mutation uses uniform prior coordinates.",
            "Completed-stage particles have equal weights; particle count and incremental-weight ESS are not independent-sample ESS or convergence certificates.",
            "When screening is present, weights and stored likelihoods remain authoritative. Two independent accept/reject draws correct the deterministic approximation.",
            "A screen must be an immutable function of a single candidate and context, independent of batch companions; identity alone cannot prove that property.",
            "Screening and authoritative failures stop the run. No fallback, support filtering or approximation-only acceptance is performed.",
            "Screen values are recomputed at completed-stage restart. Screening identity, evaluation budget and RNG continuation are bound by the plan.",
            "Completed beta=1 does not prove adequate exploration. Use repeated seeds, particle/sweep sensitivity and independent validation.",
            "Checkpoints preserve completed work, not evaluator side effects or cross-filesystem crash atomicity.",
            "No marginal-likelihood estimate, clinical validity, measured acceleration or global identifiability claim is produced."
        ])
    }
}
