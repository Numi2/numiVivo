import Foundation

/// Explicit, caller-declared sample-support policy. Passing it permits an experiment,
/// not an accuracy claim. No universal minimum-positive threshold is implied.
public enum VivoBinderTrainingSupport {
    public struct Policy: Codable, Sendable {
        public let schemaVersion: Int
        public let identifier: String
        public let rationale: String
        public let minimumPositiveGroups: Int
        public let minimumNegativeGroups: Int
        public let minimumTrainingTargets: Int
        public let minimumTargetsWithBothClasses: Int
        public let maximumFeatures: Int
        public init(identifier: String, rationale: String, minimumPositiveGroups: Int,
                    minimumNegativeGroups: Int, minimumTrainingTargets: Int,
                    minimumTargetsWithBothClasses: Int, maximumFeatures: Int) {
            schemaVersion = 1; self.identifier = identifier; self.rationale = rationale
            self.minimumPositiveGroups = minimumPositiveGroups; self.minimumNegativeGroups = minimumNegativeGroups
            self.minimumTrainingTargets = minimumTrainingTargets
            self.minimumTargetsWithBothClasses = minimumTargetsWithBothClasses; self.maximumFeatures = maximumFeatures
        }
    }
    public struct Counts: Codable, Sendable {
        public let candidates: Int
        public let positiveRows: Int
        public let negativeRows: Int
        /// A mixed-outcome group is neither an exclusively positive nor negative group.
        public let positiveGroups: Int
        public let negativeGroups: Int
        public let mixedGroups: Int
    }
    public struct Report: Codable, Sendable {
        public let schemaVersion: Int
        public let sourceSHA256: String
        public let planSHA256: String
        public let policy: Policy
        public let groupingMethod: String
        public let trainingIDs: [String]
        public let trainingExclusions: [String: String]
        public let overall: Counts
        public let byTarget: [String: Counts]
        public let featureCount: Int
        public let eligibleForExperimentalFit: Bool
        public let unmetRequirements: [String]
        public let limitations: [String]
    }
    public struct Evaluation: Codable, Sendable {
        public let schemaVersion: Int
        public let disposition: String
        public let support: Report
        /// Absent when support is insufficient; no unsupported fit is run or published.
        public let fittedReport: VivoBinderBenchmark.Report?
    }

    public static func decodePolicy(_ data: Data) throws -> Policy {
        guard let object = try JSONSerialization.jsonObject(with: data) as? [String: Any],
              Set(object.keys) == ["schemaVersion", "identifier", "rationale", "minimumPositiveGroups",
                  "minimumNegativeGroups", "minimumTrainingTargets", "minimumTargetsWithBothClasses", "maximumFeatures"] else {
            throw invalid("unknown or missing training-support policy fields")
        }
        return try VivoCanonicalJSON.decode(Policy.self, from: data)
    }

    public static func assess(_ dataset: VivoBinderBenchmark.Dataset, plan: VivoBinderBenchmark.Plan,
                              policy: Policy) throws -> Report {
        guard policy.schemaVersion == 1, !policy.identifier.isEmpty, policy.identifier.utf8.count <= 1024,
              !policy.rationale.isEmpty, policy.rationale.utf8.count <= 8192,
              (1...100_000).contains(policy.minimumPositiveGroups),
              (1...100_000).contains(policy.minimumNegativeGroups),
              (1...10_000).contains(policy.minimumTrainingTargets),
              (1...10_000).contains(policy.minimumTargetsWithBothClasses),
              (1...32).contains(policy.maximumFeatures) else { throw invalid("invalid explicit training-support policy") }
        let cohort = try VivoBinderBenchmark.partition(dataset, plan: plan)
        let training = cohort.training, overall = counts(training)
        let byTarget = Dictionary(grouping: training, by: \.target).mapValues(counts)
        let both = byTarget.values.filter { $0.positiveGroups > 0 && $0.negativeGroups > 0 }.count
        let declaredTrainingIDs = Set(dataset.records.filter { plan.trainingTargets.contains($0.target) }.map(\.id))
        var unmet: [String] = []
        if training.count < 4 { unmet.append("fewer than four eligible training candidates") }
        if overall.positiveGroups < policy.minimumPositiveGroups { unmet.append("insufficient exclusively positive training groups") }
        if overall.negativeGroups < policy.minimumNegativeGroups { unmet.append("insufficient exclusively negative training groups") }
        if byTarget.count < policy.minimumTrainingTargets { unmet.append("insufficient represented training targets") }
        if both < policy.minimumTargetsWithBothClasses { unmet.append("insufficient training targets containing both group classes") }
        if plan.modelFeatures.count > policy.maximumFeatures { unmet.append("feature count exceeds declared support policy") }
        return Report(schemaVersion: 1, sourceSHA256: dataset.sourceSHA256,
            planSHA256: try VivoBinderRanking.fingerprint(plan), policy: policy, groupingMethod: dataset.groupingMethod,
            trainingIDs: training.map(\.id), trainingExclusions: cohort.excluded.filter { declaredTrainingIDs.contains($0.key) },
            overall: overall, byTarget: byTarget, featureCount: plan.modelFeatures.count,
            eligibleForExperimentalFit: unmet.isEmpty, unmetRequirements: unmet,
            limitations: ["Thresholds are a declared development policy, not calibrated sample-size or accuracy guarantees.",
                "Groups are those supplied by the existing dataset; exact-sequence grouping does not establish homology or experimental independence.",
                "Mixed-outcome groups remain visible and do not count as exclusively positive or negative support.",
                "Only training counts and training exclusions determine admission; test outcomes do not grant support.",
                "This check does not establish coefficient stability, uncertainty calibration, causal validity or generalization."])
    }

    public static func evaluate(_ dataset: VivoBinderBenchmark.Dataset, plan: VivoBinderBenchmark.Plan,
                                policy: Policy) throws -> Evaluation {
        let support = try assess(dataset, plan: plan, policy: policy)
        guard support.eligibleForExperimentalFit else {
            return Evaluation(schemaVersion: 1, disposition: "retainFixedRanking", support: support, fittedReport: nil)
        }
        // Numerical failures remain errors. A passed sample policy must never hide a failed fit.
        return Evaluation(schemaVersion: 1, disposition: "experimentalFitOnly", support: support,
                          fittedReport: try VivoBinderBenchmark.evaluate(dataset, plan: plan))
    }

    private static func counts(_ rows: [VivoBinderBenchmark.Record]) -> Counts {
        let groups = Dictionary(grouping: rows, by: \.leakageGroup)
        var positive = 0, negative = 0, mixed = 0
        for members in groups.values {
            let n = members.filter { $0.outcome == .binder }.count
            if n == members.count { positive += 1 }
            else if n == 0 { negative += 1 }
            else { mixed += 1 }
        }
        let positives = rows.filter { $0.outcome == .binder }.count
        return Counts(candidates: rows.count, positiveRows: positives, negativeRows: rows.count - positives,
                      positiveGroups: positive, negativeGroups: negative, mixedGroups: mixed)
    }
    private static func invalid(_ text: String) -> VivoBinderBenchmark.Failure { .invalid(text) }
}
