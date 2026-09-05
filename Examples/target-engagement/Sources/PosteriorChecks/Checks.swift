import Foundation
import NumiVivoKit

private struct CheckFailure: Error, LocalizedError {
    let message: String
    var errorDescription: String? { message }
}
private struct PauseAfterStage: Error, Sendable {}

/// Regression executable for the user's Apple environment. Its presence does
/// not establish execution, statistical calibration, or biological validation.
@main struct PosteriorChecks {
    static func main() async throws {
        var checks = 0
        func check(_ condition: Bool, _ message: String) throws {
            guard condition else { throw CheckFailure(message: message) }
            checks += 1
        }
        let evidence = try VivoCanonicalJSON.fingerprint(Data("synthetic correlated Gaussian likelihood, version 1".utf8))
        let parameters = [
            VivoPosteriorParameter(identifier: "x", unit: "1", lower: 0, upper: 1, prior: .uniformPhysical),
            VivoPosteriorParameter(identifier: "y", unit: "1", lower: 0, upper: 1, prior: .uniformPhysical)
        ]
        let plan = VivoPosteriorPlan(likelihoodFingerprint: evidence, parameters: parameters,
            configuration: .init(particleCount: 256, mutationSweeps: 8, seed: 1234567))
        let linear = try VivoPosteriorParameter(identifier: "linear", unit: "1", lower: 1, upper: 100,
            prior: .uniformPhysical).value(at: 0.5)
        let logarithmic = try VivoPosteriorParameter(identifier: "log", unit: "1", lower: 1, upper: 100,
            prior: .uniformLogPhysical).value(at: 0.5)
        try check(abs(linear - 50.5) < 1e-12 && abs(logarithmic - 10) < 1e-12, "prior distributions differ explicitly")

        let evaluate: VivoPosteriorBatchEvaluator = { candidates in
            // Deliberately reverse results: asynchronous completion order must
            // not select a different random trajectory or reassign likelihoods.
            candidates.reversed().map { candidate in
                let x = (candidate.values[0] - 0.35) / 0.08
                let y = (candidate.values[1] - 0.65) / 0.05
                let rho = 0.75
                return .init(candidate: candidate, logLikelihood: -0.5 * (x*x - 2*rho*x*y + y*y) / (1-rho*rho))
            }
        }
        let sampler = try VivoTemperedPosteriorSampler(plan: plan)
        let complete = try await sampler.run(evaluate: evaluate)
        try complete.validate(requireComplete: true)
        guard let final = complete.checkpoint else { throw CheckFailure(message: "missing final checkpoint") }
        try check(final.beta == 1, "tempering completed")
        try check(final.stages.count > 1, "nontrivial tempering schedule")
        let n = Double(final.particles.count)
        let mx = final.particles.reduce(0) { $0 + $1.coordinates[0] / n }
        let my = final.particles.reduce(0) { $0 + $1.coordinates[1] / n }
        let vx = final.particles.reduce(0) { $0 + pow($1.coordinates[0] - mx, 2) / n }
        let vy = final.particles.reduce(0) { $0 + pow($1.coordinates[1] - my, 2) / n }
        let covariance = final.particles.reduce(0) { $0 + ($1.coordinates[0] - mx) * ($1.coordinates[1] - my) / n }
        try check(abs(mx - 0.35) < 0.04 && abs(my - 0.65) < 0.04, "correlated Gaussian posterior means")
        try check(covariance / sqrt(vx * vy) > 0.4, "joint posterior correlation retained")
        try check(final.stages.allSatisfy { $0.acceptedMoves <= $0.likelihoodEvaluations }, "mutation ledger")

        let interrupted = try VivoTemperedPosteriorSampler(plan: plan)
        let partial = try await interrupted.run(evaluate: evaluate, progress: { checkpoint in
            if checkpoint.stages.count == 1 { throw PauseAfterStage() }
        })
        try check(!partial.completed && partial.failure != nil, "interrupted output is not a posterior")
        guard let boundary = partial.checkpoint else { throw CheckFailure(message: "no completed stage retained") }
        let saved = try VivoCanonicalJSON.decode(VivoPosteriorCheckpoint.self, from: VivoCanonicalJSON.encode(boundary))
        let restored = try VivoTemperedPosteriorSampler(plan: plan, checkpoint: saved)
        let resumed = try await restored.run(evaluate: evaluate)
        try check(resumed.checkpoint == complete.checkpoint && resumed.completed, "stage-boundary resume reproduces particle/RNG history")
        do {
            let changed = VivoPosteriorPlan(likelihoodFingerprint: evidence,
                parameters: [parameters[0], .init(identifier: "y", unit: "1", lower: 0, upper: 0.9, prior: .uniformPhysical)],
                configuration: plan.configuration)
            _ = try VivoTemperedPosteriorSampler(plan: changed, checkpoint: saved)
            throw CheckFailure(message: "foreign prior accepted at resume")
        } catch VivoPosteriorError.invalid { checks += 1 }
        do {
            try partial.validate(requireComplete: true)
            throw CheckFailure(message: "partial posterior accepted for prediction")
        } catch VivoPosteriorError.invalid { checks += 1 }

        let failing = try VivoTemperedPosteriorSampler(plan: plan)
        let failed = try await failing.run(evaluate: { candidates in
            candidates.map {
                $0.ordinal == 0 ? .init(candidate: $0, failure: "synthetic numerical failure")
                    : .init(candidate: $0, logLikelihood: 0)
            }
        })
        try check(!failed.completed && failed.checkpoint == nil, "initial evaluation failure not discarded")
        try check(failed.failure?.failedEvaluations.count == 1, "failed parameter candidate remains visible")
        let missing = try VivoTemperedPosteriorSampler(plan: plan)
        let omitted = try await missing.run(evaluate: { candidates in
            candidates.dropFirst().map { .init(candidate: $0, logLikelihood: 0) }
        })
        try check(omitted.failure != nil && !omitted.completed, "omitted likelihood cannot pass")

        let root = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
            .deletingLastPathComponent().deletingLastPathComponent()
        let problem = try VivoKineticsDocumentIO.read(VivoTargetPosteriorProblem.self,
            from: root.appendingPathComponent("synthetic-inference.json"))
        let prepared = try VivoPreparedTargetPosterior(problem)
        let changedHeldout = try changingObservation(problem, caseIndex: 2, value: 0.99)
        let heldoutPlan = try VivoPreparedTargetPosterior(changedHeldout).plan
        try check(heldoutPlan == prepared.plan, "held-out values do not enter likelihood identity")
        try check(try prepared.logLikelihood(values: [100000, 0.2]) == VivoPreparedTargetPosterior(changedHeldout)
            .logLikelihood(values: [100000, 0.2]), "held-out values do not enter likelihood")
        let changedTraining = try changingObservation(problem, caseIndex: 0, value: 0.99)
        let changedPlan = try VivoPreparedTargetPosterior(changedTraining).plan
        try check(changedPlan.likelihoodFingerprint != prepared.plan.likelihoodFingerprint, "training changes invalidate resume identity")

        let fitted = try await VivoTargetPosteriorFitter.run(problem)
        try fitted.validate(requireComplete: true)
        checks += 1
        let predictive = try await VivoTargetPosteriorPredictor.predict(fitted)
        try check(predictive.cases.count == 4 && predictive.cases.allSatisfy { $0.failures.isEmpty }, "all prediction cases retained")
        try check(predictive.parameterSummaries.count == 2 && predictive.priorCoordinateCorrelationRowMajor.count == 4,
                  "joint parameter summaries")
        try check(predictive.cases.filter { $0.partition != .calibration }.count == 2, "prediction keeps held-out labels")
        for item in predictive.cases {
            for point in item.points {
                try check(point.credibleLower <= point.credibleMedian && point.credibleMedian <= point.credibleUpper,
                          "ordered latent credible interval")
                guard let lower = point.predictiveLower, let upper = point.predictiveUpper else {
                    throw CheckFailure(message: "explicit assay noise omitted")
                }
                try check(lower < upper && point.logPredictiveDensity?.isFinite == true,
                          "finite measurement predictive distribution")
                try check(point.observationOrigin == .assumed, "synthetic evidence was not relabelled measured")
            }
        }
        let sensitivity = try VivoTargetPosteriorSensitivity.evaluate(fitted)
        try check(sensitivity.derivativeChecksPassed && sensitivity.numericalRank == 2, "synthetic local two-rate information")
        print("\(checks) posterior regression checks passed in this invocation. These are numerical/software checks, not biological validation.")
    }

    private static func changingObservation(_ problem: VivoTargetPosteriorProblem, caseIndex: Int,
                                            value: Double) throws -> VivoTargetPosteriorProblem {
        var json = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(problem)) as! [String: Any]
        var study = json["study"] as! [String: Any]
        var cases = study["cases"] as! [[String: Any]]
        var observations = cases[caseIndex]["observations"] as! [[String: Any]]
        observations[0]["value"] = value
        cases[caseIndex]["observations"] = observations
        study["cases"] = cases; json["study"] = study
        return try VivoCanonicalJSON.decode(VivoTargetPosteriorProblem.self,
                                             from: JSONSerialization.data(withJSONObject: json))
    }
}
