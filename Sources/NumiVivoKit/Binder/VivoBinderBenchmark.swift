import Foundation

/// Retrospective candidate ranking, not affinity, simulation or clinical evidence.
public enum VivoBinderBenchmark {
    public enum Failure: Error, CustomStringConvertible {
        case invalid(String)
        public var description: String { switch self { case .invalid(let reason): return reason } }
    }
    public enum Outcome: String, Codable, Sendable {
        case binder, nonBinder, inconclusive, notTested, expressionFailure
        public var binary: Double? {
            switch self { case .binder: return 1; case .nonBinder: return 0; default: return nil }
        }
    }
    public struct Record: Codable, Sendable {
        public let id: String
        public let target: String
        public let leakageGroup: String
        public let outcome: Outcome
        public let rawOutcome: String
        /// Numeric predictors only. Missing values must be omitted, never replaced with zero.
        public let features: [String: Double]
        public let sourceFields: [String: String]
        public init(id: String, target: String, leakageGroup: String, outcome: Outcome,
                    rawOutcome: String, features: [String: Double], sourceFields: [String: String] = [:]) {
            self.id = id; self.target = target; self.leakageGroup = leakageGroup
            self.outcome = outcome; self.rawOutcome = rawOutcome
            self.features = features; self.sourceFields = sourceFields
        }
    }
    public struct Dataset: Codable, Sendable {
        public let schemaVersion: Int
        public let sourceSHA256: String
        public let assay: String
        public let groupingMethod: String
        public let records: [Record]
        public init(sourceSHA256: String, assay: String, groupingMethod: String, records: [Record]) {
            self.schemaVersion = 1; self.sourceSHA256 = sourceSHA256; self.assay = assay
            self.groupingMethod = groupingMethod; self.records = records
        }
    }
    public struct Plan: Codable, Sendable {
        public let schemaVersion: Int
        public let sourceSHA256: String
        public let trainingTargets: [String]
        public let testTargets: [String]
        public let baselineFeature: String
        public let modelFeatures: [String]
        public let topK: Int
        /// Objective: mean binary cross-entropy + penalty/2 * squared coefficient norm.
        /// The intercept is not penalized. All scaling is learned from training rows only.
        public let ridgePenalty: Double
        public init(sourceSHA256: String, trainingTargets: [String], testTargets: [String],
                    baselineFeature: String, modelFeatures: [String], topK: Int = 10,
                    ridgePenalty: Double = 0.1) {
            self.schemaVersion = 1; self.sourceSHA256 = sourceSHA256
            self.trainingTargets = trainingTargets; self.testTargets = testTargets
            self.baselineFeature = baselineFeature; self.modelFeatures = modelFeatures
            self.topK = topK; self.ridgePenalty = ridgePenalty
        }
    }
    public struct Metrics: Codable, Sendable {
        public let count: Int
        public let positives: Int
        public let effectiveK: Int
        public let expectedHitsAtK: Double
        public let precisionAtK: Double
        public let auroc: Double?
        /// Threshold-block average precision; equal scores are evaluated together.
        public let averagePrecision: Double?
        public let brier: Double?
    }
    public struct Model: Codable, Sendable {
        public let featureNames: [String]
        public let means: [Double]
        public let scales: [Double]
        public let coefficients: [Double]
        public let intercept: Double
        public let iterations: Int
        public let gradientInfinityNorm: Double
    }
    public struct TargetResult: Codable, Sendable {
        public let target: String
        public let candidateIDs: [String]
        public let labels: [Double]
        public let baselineScores: [Double]
        public let modelProbabilities: [Double]
        public let baseline: Metrics
        public let learned: Metrics
        public let trainingPrevalence: Metrics
    }
    public struct Report: Codable, Sendable {
        public let schemaVersion: Int
        public let evidenceStatus: String
        public let sourceSHA256: String
        public let assay: String
        public let groupingMethod: String
        public let plan: Plan
        public let trainingIDs: [String]
        public let excluded: [String: String]
        public let model: Model
        public let targets: [TargetResult]
        public let limitations: [String]
    }

    public static func validate(_ dataset: Dataset) throws {
        try require(dataset.schemaVersion == 1, "unsupported dataset schema")
        try require(isSHA256(dataset.sourceSHA256), "sourceSHA256 must be lowercase SHA-256")
        try require(!dataset.assay.isEmpty && !dataset.groupingMethod.isEmpty, "missing assay or grouping method")
        try require(!dataset.records.isEmpty && dataset.records.count <= 100_000, "record count outside 1...100000")
        var ids = Set<String>()
        for row in dataset.records {
            try require(!row.id.isEmpty && ids.insert(row.id).inserted, "empty or duplicate candidate ID: \(row.id)")
            try require(!row.target.isEmpty && !row.leakageGroup.isEmpty && !row.rawOutcome.isEmpty, "missing identity/outcome: \(row.id)")
            try require(row.features.count <= 128, "too many features: \(row.id)")
            for (name, value) in row.features {
                try require(!name.isEmpty && value.isFinite && abs(value) <= 1e12, "invalid feature \(name): \(row.id)")
            }
        }
    }

    public static func evaluate(_ dataset: Dataset, plan: Plan) throws -> Report {
        try validate(dataset)
        try require(plan.schemaVersion == 1 && plan.sourceSHA256 == dataset.sourceSHA256, "plan/source mismatch")
        let trainTargets = Set(plan.trainingTargets), testTargets = Set(plan.testTargets)
        try require(!trainTargets.isEmpty && !testTargets.isEmpty, "explicit nonempty training and test targets required")
        try require(trainTargets.count == plan.trainingTargets.count && testTargets.count == plan.testTargets.count,
                    "duplicate target in plan")
        try require(trainTargets.isDisjoint(with: testTargets), "training/test target overlap")
        let known = Set(dataset.records.map(\.target))
        try require(trainTargets.union(testTargets).isSubset(of: known), "unknown target in plan")
        try require((1...100_000).contains(plan.topK), "topK outside 1...100000")
        try require(!plan.baselineFeature.isEmpty && (1...32).contains(plan.modelFeatures.count)
                    && Set(plan.modelFeatures).count == plan.modelFeatures.count
                    && plan.modelFeatures.allSatisfy { !$0.isEmpty }, "invalid model features")
        try require(plan.ridgePenalty.isFinite && plan.ridgePenalty >= 0.001 && plan.ridgePenalty <= 100,
                    "ridge penalty outside 0.001...100")
        let required = Set(plan.modelFeatures + [plan.baselineFeature])
        // Purge training groups matching ANY held-out candidate, before looking at outcomes.
        let testGroups = Set(dataset.records.filter { testTargets.contains($0.target) }.map(\.leakageGroup))
        var excluded: [String: String] = [:], train: [Record] = [], test: [Record] = []
        for row in dataset.records.sorted(by: { $0.id < $1.id }) {
            if !trainTargets.contains(row.target) && !testTargets.contains(row.target) {
                excluded[row.id] = "target outside plan"; continue
            }
            if trainTargets.contains(row.target) && testGroups.contains(row.leakageGroup) {
                excluded[row.id] = "training group overlaps test"; continue
            }
            guard row.outcome.binary != nil else { excluded[row.id] = "outcome:\(row.outcome.rawValue)"; continue }
            let missing = required.subtracting(row.features.keys).sorted()
            guard missing.isEmpty else { excluded[row.id] = "missing features:\(missing.joined(separator: ","))"; continue }
            if trainTargets.contains(row.target) { train.append(row) } else { test.append(row) }
        }
        try require(train.count >= 4, "fewer than four eligible training candidates")
        let positives = train.reduce(0.0) { $0 + ($1.outcome.binary ?? 0) }
        try require(positives > 0 && positives < Double(train.count), "training requires both outcome classes")
        for target in testTargets { try require(test.contains { $0.target == target }, "no matched eligible candidates for \(target)") }
        let model = try fit(train, features: plan.modelFeatures, penalty: plan.ridgePenalty)
        var results: [TargetResult] = []
        for target in testTargets.sorted() {
            let rows = test.filter { $0.target == target }
            let labels = rows.map { $0.outcome.binary! }
            let baseline = rows.map { $0.features[plan.baselineFeature]! }
            let predicted = rows.map { row in
                sigmoid(zip(model.featureNames.indices, model.coefficients).reduce(model.intercept) { sum, item in
                    let j = item.0
                    return sum + item.1 * ((row.features[model.featureNames[j]]! - model.means[j]) / model.scales[j])
                })
            }
            results.append(TargetResult(target: target, candidateIDs: rows.map(\.id), labels: labels,
                baselineScores: baseline, modelProbabilities: predicted,
                baseline: try metrics(labels: labels, scores: baseline, topK: plan.topK),
                learned: try metrics(labels: labels, scores: predicted, topK: plan.topK, probabilities: true),
                trainingPrevalence: try metrics(labels: labels,
                    scores: Array(repeating: positives / Double(train.count), count: rows.count),
                    topK: plan.topK, probabilities: true)))
        }
        return Report(schemaVersion: 1, evidenceStatus: "retrospective-development-only",
            sourceSHA256: dataset.sourceSHA256, assay: dataset.assay, groupingMethod: dataset.groupingMethod,
            plan: plan, trainingIDs: train.map(\.id), excluded: excluded, model: model, targets: results,
            limitations: ["No prospective, affinity, causal, physics or clinical qualification.",
                "Both methods use the same complete-case candidates; inspect coverage exclusions.",
                "Grouping is caller-declared; exact-sequence groups do not control near homology.",
                "No independence claim or confidence intervals across related designs/targets.",
                "The plan must be frozen before test-outcome inspection; code cannot attest that history.",
                "Raw baseline scores are rankings, not calibrated binding probabilities."])
    }

    /// Tie-invariant metrics. Single-class AUROC is unavailable, not zero or one.
    public static func metrics(labels: [Double], scores: [Double], topK: Int,
                               probabilities: Bool = false) throws -> Metrics {
        try require(!labels.isEmpty && labels.count == scores.count && labels.count <= 100_000 && topK > 0,
                    "invalid metric dimensions")
        try require(labels.allSatisfy { $0 == 0 || $0 == 1 } && scores.allSatisfy(\.isFinite), "invalid metric values")
        if probabilities { try require(scores.allSatisfy { (0...1).contains($0) }, "probability outside [0,1]") }
        let order = scores.indices.sorted { scores[$0] > scores[$1] }
        let positive = Int(labels.reduce(0, +)), negative = labels.count - positive, k = min(topK, labels.count)
        var seen = 0, seenPositive = 0.0, hits = 0.0, aucNumerator = 0.0, ap = 0.0
        while seen < order.count {
            var end = seen + 1
            while end < order.count && scores[order[end]] == scores[order[seen]] { end += 1 }
            let groupPositive = order[seen..<end].reduce(0.0) { $0 + labels[$1] }
            let size = end - seen, groupNegative = Double(size) - groupPositive
            let negativeBelow = Double(negative) - (Double(seen) - seenPositive) - groupNegative
            aucNumerator += groupPositive * (negativeBelow + 0.5 * groupNegative)
            if seen < k { hits += Double(min(end, k) - seen) * groupPositive / Double(size) }
            seenPositive += groupPositive
            if positive > 0 { ap += groupPositive / Double(positive) * seenPositive / Double(end) }
            seen = end
        }
        let brier = probabilities ? zip(labels, scores).reduce(0.0) { $0 + pow($1.0 - $1.1, 2) } / Double(labels.count) : nil
        return Metrics(count: labels.count, positives: positive, effectiveK: k, expectedHitsAtK: hits,
            precisionAtK: hits / Double(k), auroc: positive > 0 && negative > 0 ? aucNumerator / (Double(positive) * Double(negative)) : nil,
            averagePrecision: positive > 0 ? ap : nil, brier: brier)
    }

    private static func fit(_ rows: [Record], features: [String], penalty: Double) throws -> Model {
        let n = Double(rows.count), p = features.count
        let means = features.map { name in rows.reduce(0.0) { $0 + $1.features[name]! } / n }
        let scales = features.indices.map { j -> Double in
            let variance = rows.reduce(0.0) { $0 + pow($1.features[features[j]]! - means[j], 2) } / n
            return variance > 1e-20 ? sqrt(variance) : 1
        }
        let x = rows.map { row in features.indices.map { (row.features[features[$0]]! - means[$0]) / scales[$0] } }
        let y = rows.map { $0.outcome.binary! }, prevalence = rows.reduce(0.0) { $0 + $1.outcome.binary! } / n
        var b = log(prevalence / (1 - prevalence)), w = Array(repeating: 0.0, count: p)
        // Trace bound on the standardized design's spectral norm gives a safe fixed step.
        let step = 1 / (0.25 * Double(p + 1) + penalty)
        var norm = Double.infinity
        for iteration in 0..<20_000 {
            var gb = 0.0, gw = w.map { penalty * $0 }
            for i in rows.indices {
                let error = sigmoid(zip(w, x[i]).reduce(b) { $0 + $1.0 * $1.1 }) - y[i]
                gb += error / n
                for j in 0..<p { gw[j] += error * x[i][j] / n }
            }
            norm = max(abs(gb), gw.map(abs).max() ?? 0)
            try require(norm.isFinite && b.isFinite && w.allSatisfy(\.isFinite), "nonfinite logistic fit")
            if norm <= 1e-7 {
                return Model(featureNames: features, means: means, scales: scales, coefficients: w,
                             intercept: b, iterations: iteration, gradientInfinityNorm: norm)
            }
            b -= step * gb
            for j in 0..<p { w[j] -= step * gw[j] }
        }
        throw Failure.invalid("logistic fit did not converge; gradient=\(norm)")
    }
    private static func sigmoid(_ x: Double) -> Double {
        if x >= 0 { return 1 / (1 + exp(-x)) }
        let e = exp(x); return e / (1 + e)
    }
    private static func isSHA256(_ value: String) -> Bool {
        value.utf8.count == 64 && value.utf8.allSatisfy { (48...57).contains($0) || (97...102).contains($0) }
    }
    private static func require(_ value: Bool, _ message: String) throws {
        if !value { throw Failure.invalid(message) }
    }
}
