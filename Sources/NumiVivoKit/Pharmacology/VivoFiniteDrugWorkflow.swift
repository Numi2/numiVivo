import Foundation

public struct VivoFiniteDrugExperiment: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let model: VivoFiniteDrugModel
    public let initial: VivoFiniteDrugState
    public let durationSeconds: Double
    public let policy: VivoFiniteDrugPolicy
    public init(model: VivoFiniteDrugModel, initial: VivoFiniteDrugState, durationSeconds: Double,
                policy: VivoFiniteDrugPolicy = .init()) {
        schemaVersion = 1; self.model = model; self.initial = initial
        self.durationSeconds = durationSeconds; self.policy = policy
    }
}

public struct VivoFiniteDrugRunRecord: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let numericalMethod: String
    public let experimentFingerprint: VivoFingerprint
    public let experiment: VivoFiniteDrugExperiment
    public let result: VivoFiniteDrugStep
    public let limitations: [String]
    public static func run(_ experiment: VivoFiniteDrugExperiment) throws -> Self {
        guard experiment.schemaVersion == 1 else { throw VivoKineticsError.invalid("finite-drug experiment schema") }
        let step = try VivoFiniteDrugReactions.advance(experiment.initial, model: experiment.model,
            durationSeconds: experiment.durationSeconds, policy: experiment.policy)
        return try .init(schemaVersion: 1, numericalMethod: "numivivo.finite-drug.positive-split-fp64.v1",
            experimentFingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(experiment)),
            experiment: experiment, result: step, limitations: [
                "Fixed-volume local chemical operator: concentration in M, elapsed time in seconds. No dosing or circulation model is run.",
                "Drug is consumed by binding. Removed bound target transfers its drug to one declared inactive metabolite; clearance transfers drug equivalents to a ledger.",
                "Target synthesis and removal are explicitly accounted; this model uses one shared target turnover rate and constant synthesis.",
                "Step-doubling controls local splitting error, not a rigorous global solution error. Material balance is checked independently.",
                "Caller owns coupled simulation state and commit/rollback; applying this operator alongside another owner of the same reactions would double count.",
                "No active-metabolite binding, spatial transport, clinical efficacy, toxicity or patient-specific prediction is inferred."
            ])
    }
}
