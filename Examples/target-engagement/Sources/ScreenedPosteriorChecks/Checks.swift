import Foundation
import NumiVivoKit

private struct CheckFailure: Error, LocalizedError {
    let message: String
    var errorDescription: String? { message }
}
private struct Pause: Error, Sendable {}

/// Source for user-side execution. No result or passed-check count is implied
/// merely by committing this executable. --metal requires actual Apple hardware.
@main struct ScreenedPosteriorChecks {
    static func target(_ x: Double) -> Double { -0.5 * pow((x - 0.35) / 0.08, 2) }
    static func approximation(_ x: Double) -> Double { -0.5 * pow((x - 0.50) / 0.22, 2) }

    static func main() async throws {
        var checks = 0
        func check(_ condition: Bool, _ message: String) throws {
            guard condition else { throw CheckFailure(message: message) }
            checks += 1
        }
        let targetID = try VivoCanonicalJSON.fingerprint(Data("synthetic bounded Gaussian authority v1".utf8))
        let screenID = try VivoCanonicalJSON.fingerprint(Data("synthetic displaced wider Gaussian screen v1".utf8))
        let base = VivoPosteriorPlan(likelihoodFingerprint: targetID,
            parameters: [.init(identifier: "x", unit: "1", lower: 0, upper: 1, prior: .uniformPhysical)],
            configuration: .init(particleCount: 512, mutationSweeps: 8, seed: 41071))
        let policy = VivoPosteriorScreeningPolicy(fingerprint: screenID)
        let plan = base.withScreening(policy)
        let exact: VivoPosteriorBatchEvaluator = { values in
            values.map { .init(candidate: $0, logLikelihood: target($0.values[0])) }
        }
        let approximate: VivoPosteriorBatchEvaluator = { values in
            values.reversed().map { .init(candidate: $0, logLikelihood: approximation($0.values[0])) }
        }
        let screen = VivoPosteriorScreen(policy: policy, evaluate: approximate)
        let sampler = try VivoTemperedPosteriorSampler(plan: plan)
        let completed = try await sampler.run(evaluate: exact, screen: screen)
        try completed.validate(requireComplete: true)
        guard let final = completed.checkpoint else { throw CheckFailure(message: "missing completed population") }
        try check(final.stages.count > 1, "nontrivial tempering")
        try check(final.stages.reduce(0, { $0 + ($1.screenedOutMoves ?? 0) }) > 0, "screen bypasses some authoritative evaluations")
        try check(final.stages.allSatisfy {
            $0.likelihoodEvaluations + ($0.screenedOutMoves ?? 0) + $0.outOfPriorMoves == $0.proposedMoves
        }, "separate screening and authoritative ledgers")
        try check(final.particles.allSatisfy { $0.logLikelihood == target($0.coordinates[0]) }, "stored likelihoods remain authoritative")
        let mean = final.particles.reduce(0) { $0 + $1.coordinates[0] / Double(final.particles.count) }
        try check(abs(mean - 0.35) < 0.04, "imperfect screen does not replace target distribution")
        try check(plan.likelihoodFingerprint == base.likelihoodFingerprint, "authoritative likelihood identity retained")
        try check(try plan.fingerprint() != base.fingerprint(), "screen execution changes checkpoint identity")

        let pausedSampler = try VivoTemperedPosteriorSampler(plan: plan)
        let paused = try await pausedSampler.run(evaluate: exact, progress: { checkpoint in
            if checkpoint.stages.count == 1 { throw Pause() }
        }, screen: screen)
        try check(!paused.completed && paused.failure != nil, "paused run is not presented as complete")
        guard let boundary = paused.checkpoint else { throw CheckFailure(message: "no completed boundary retained") }
        let saved = try VivoCanonicalJSON.decode(VivoPosteriorCheckpoint.self, from: VivoCanonicalJSON.encode(boundary))
        let resumedSampler = try VivoTemperedPosteriorSampler(plan: plan, checkpoint: saved)
        let resumed = try await resumedSampler.run(evaluate: exact, screen: screen)
        try check(resumed.checkpoint == completed.checkpoint, "screened checkpoint continuation is deterministic")
        do {
            _ = try VivoTemperedPosteriorSampler(plan: base, checkpoint: saved)
            throw CheckFailure(message: "screened checkpoint resumed without screen identity")
        } catch VivoPosteriorError.invalid { checks += 1 }
        do {
            _ = try await sampler.run(evaluate: exact)
            throw CheckFailure(message: "missing bound screen accepted")
        } catch VivoPosteriorError.invalid { checks += 1 }

        let exactScreenPolicy = VivoPosteriorScreeningPolicy(fingerprint: targetID)
        let exactScreenSampler = try VivoTemperedPosteriorSampler(plan: base.withScreening(exactScreenPolicy))
        let perfect = try await exactScreenSampler.run(evaluate: exact,
            screen: .init(policy: exactScreenPolicy, evaluate: exact))
        try perfect.validate(requireComplete: true)
        try check(perfect.checkpoint!.stages.allSatisfy { $0.acceptedMoves == $0.likelihoodEvaluations },
                  "perfect screen requires no second-stage rejection")

        let constantPolicy = VivoPosteriorScreeningPolicy(fingerprint: try VivoCanonicalJSON.fingerprint(Data("constant-zero-screen".utf8)))
        let constantSampler = try VivoTemperedPosteriorSampler(plan: base.withScreening(constantPolicy))
        let constant = try await constantSampler.run(evaluate: exact, screen: .init(policy: constantPolicy, evaluate: { values in
            values.map { .init(candidate: $0, logLikelihood: 0) }
        }))
        try constant.validate(requireComplete: true)
        try check(constant.checkpoint!.stages.allSatisfy { $0.screenedOutMoves == 0 }, "constant screen retains every valid proposal")

        let failedSampler = try VivoTemperedPosteriorSampler(plan: plan)
        let failed = try await failedSampler.run(evaluate: exact, screen: .init(policy: policy, evaluate: { values in
            values.enumerated().map { index, value in
                index == 0 ? .init(candidate: value, failure: "synthetic screen failure")
                    : .init(candidate: value, logLikelihood: 0)
            }
        }))
        try check(failed.failure?.failedEvaluations.count == 1 && failed.checkpoint?.beta == 0,
                  "screen failure retains exact initial population without filtering")
        let shortPolicy = VivoPosteriorScreeningPolicy(fingerprint: screenID, maximumEvaluations: 1)
        let limited = try VivoTemperedPosteriorSampler(plan: base.withScreening(shortPolicy))
        let exhausted = try await limited.run(evaluate: exact, screen: .init(policy: shortPolicy, evaluate: approximate))
        try check(!exhausted.completed && exhausted.failure != nil && exhausted.checkpoint?.beta == 0,
                  "screen budget cannot force completion")
        let json = try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(base)) as! [String: Any]
        try check(json["screening"] == nil, "legacy plan omits optional screening field")
        let legacySampler = try VivoTemperedPosteriorSampler(plan: base)
        let legacy = try await legacySampler.run(evaluate: exact)
        try legacy.validate(requireComplete: true)
        try check(legacy.checkpoint!.stages.allSatisfy { $0.screeningEvaluations == nil && $0.screenedOutMoves == nil },
                  "legacy stage shape retained")

        if CommandLine.arguments.contains("--metal") {
            let root = URL(fileURLWithPath: #filePath).deletingLastPathComponent()
                .deletingLastPathComponent().deletingLastPathComponent()
            let problem = try VivoKineticsDocumentIO.read(VivoTargetPosteriorProblem.self,
                from: root.appendingPathComponent("synthetic-inference.json"))
            let gpu = try await VivoMetalTargetLikelihoodScreen.make(problem: problem)
            let probes = try await VivoMetalTargetScreenChecks.run(problem: problem, screen: gpu)
            try check(probes.deterministicProbeChecksPassed, "actual Metal repeat/permutation/partition/ordinal independence")
            let record = try await VivoTargetPosteriorFitter.run(problem, screen: gpu.screen())
            try record.validate(requireComplete: true)
            try check(record.posterior.plan.screening?.fingerprint == gpu.screeningPolicy.fingerprint,
                      "GPU description remains bound in result")
            let prediction = try await VivoTargetPosteriorPredictor.predict(record)
            try check(prediction.cases.allSatisfy { $0.failures.isEmpty }, "screened record uses existing FP64 predictive workflow")
        }
        print("\(checks) screened inference checks passed in this invocation; Metal executed: \(CommandLine.arguments.contains("--metal")).")
    }
}
