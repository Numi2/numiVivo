import Foundation

/// Outcome-blind fixed ranking. Experimental labels and source fields are not inputs.
/// Scores are rankings, not probabilities; target normalization is pool-dependent.
public enum VivoBinderRanking {
    public struct Candidate: Codable, Sendable {
        public let id: String
        public let target: String
        public let scores: [String: Double]
        public init(id: String, target: String, scores: [String: Double]) {
            self.id = id; self.target = target; self.scores = scores
        }
    }
    public struct Query: Codable, Sendable {
        public let schemaVersion: Int
        public let sourceSHA256: String
        public let candidates: [Candidate]
        public init(sourceSHA256: String, candidates: [Candidate]) {
            schemaVersion = 1; self.sourceSHA256 = sourceSHA256; self.candidates = candidates
        }
    }
    public enum Method: String, Codable, Sendable {
        case fixedScore
        /// Weighted mean of within-target z scores; population variance, complete cases.
        case targetStandardizedMean
    }
    public struct Component: Codable, Sendable {
        public let feature: String
        public let weight: Double
        public init(feature: String, weight: Double = 1) { self.feature = feature; self.weight = weight }
    }
    public struct Plan: Codable, Sendable {
        public let schemaVersion: Int
        public let method: Method
        public let components: [Component]
        public let topK: Int
        public init(method: Method, components: [Component], topK: Int = 10) {
            schemaVersion = 1; self.method = method; self.components = components; self.topK = topK
        }
    }
    public struct Normalization: Codable, Sendable {
        public let feature: String
        public let mean: Double
        public let populationSD: Double
        public let constant: Bool
    }
    public struct RankedCandidate: Codable, Sendable {
        public let id: String
        public let score: Double
        /// Expected inclusion at a tied cutoff, not an arbitrary ID-based winner.
        public let selectionWeight: Double
    }
    public struct TargetRanking: Codable, Sendable {
        public let target: String
        public let inputCount: Int
        public let effectiveK: Int
        public let normalization: [Normalization]
        public let candidates: [RankedCandidate]
        public let excluded: [String: [String]]
    }
    public struct Report: Codable, Sendable {
        public let schemaVersion: Int
        public let methodVersion: String
        public let querySHA256: String
        public let planSHA256: String
        public let targets: [TargetRanking]
        public let limitations: [String]
    }
    public struct TargetAssessment: Codable, Sendable {
        public let target: String
        public let rankedCount: Int
        public let effectiveK: Int
        public let binaryOutcomeCount: Int
        public let selectedWeightByOutcome: [String: Double]
        public let measuredBinderHits: Double
        public let selectedNonbinaryWeight: Double
        public let precisionAmongBinarySelected: Double?
        /// Discrimination over ranked candidates with binary outcomes; no top-K refill.
        public let binaryAUROC: Double?
        public let binaryAveragePrecision: Double?
    }
    public struct Assessment: Codable, Sendable {
        public let schemaVersion: Int
        public let evidenceStatus: String
        public let assay: String
        public let sourceSHA256: String
        public let rankingSHA256: String
        public let targets: [TargetAssessment]
        public let limitations: [String]
    }

    /// Strips ALL experimental fields. All source candidates remain, including untested ones.
    public static func query(from dataset: VivoBinderBenchmark.Dataset) throws -> Query {
        try VivoBinderBenchmark.validate(dataset)
        return Query(sourceSHA256: dataset.sourceSHA256, candidates: dataset.records.map {
            Candidate(id: $0.id, target: $0.target,
                      scores: $0.features.filter { VivoBinderAnthropicImport.allowedFeatures.contains($0.key) })
        }.sorted { $0.id < $1.id })
    }

    /// Rejects outcome-bearing/unknown JSON fields instead of silently discarding them.
    public static func decodeQuery(_ data: Data) throws -> Query {
        let value = try JSONSerialization.jsonObject(with: data)
        guard let object = value as? [String: Any], Set(object.keys) == ["schemaVersion", "sourceSHA256", "candidates"],
              let candidates = object["candidates"] as? [[String: Any]],
              candidates.allSatisfy({ Set($0.keys) == ["id", "target", "scores"] }) else {
            throw invalid("ranking query contains unknown or outcome-bearing fields")
        }
        let result = try VivoCanonicalJSON.decode(Query.self, from: data)
        try validate(result)
        return result
    }

    public static func decodePlan(_ data: Data) throws -> Plan {
        let value = try JSONSerialization.jsonObject(with: data)
        guard let object = value as? [String: Any],
              Set(object.keys) == ["schemaVersion", "method", "components", "topK"],
              let components = object["components"] as? [[String: Any]],
              components.allSatisfy({ Set($0.keys) == ["feature", "weight"] }) else {
            throw invalid("ranking plan contains unknown fields")
        }
        return try VivoCanonicalJSON.decode(Plan.self, from: data)
    }

    public static func rank(_ query: Query, plan: Plan) throws -> Report {
        try validate(query)
        try require(plan.schemaVersion == 1 && (1...100_000).contains(plan.topK), "invalid ranking plan schema/K")
        try require((1...32).contains(plan.components.count)
                    && Set(plan.components.map(\.feature)).count == plan.components.count,
                    "duplicate/empty ranking components")
        for component in plan.components {
            try require(VivoBinderAnthropicImport.allowedFeatures.contains(component.feature)
                        && component.weight.isFinite && component.weight > 0 && component.weight <= 1_000,
                        "ranking requires declared in-silico features and positive finite weights")
        }
        if plan.method == .fixedScore {
            try require(plan.components.count == 1 && plan.components[0].weight == 1,
                        "fixedScore requires one unit-weight component")
        }
        let required = Set(plan.components.map(\.feature))
        let byTarget = Dictionary(grouping: query.candidates, by: \.target)
        var targets: [TargetRanking] = []
        for target in byTarget.keys.sorted() {
            let pool = byTarget[target]!.sorted { $0.id < $1.id }
            var rows: [Candidate] = [], excluded: [String: [String]] = [:]
            for row in pool {
                let missing = required.subtracting(row.scores.keys).sorted()
                if missing.isEmpty { rows.append(row) } else { excluded[row.id] = missing }
            }
            // A target with no complete candidates is recorded, not silently removed.
            var norms: [Normalization] = []
            if !rows.isEmpty && plan.method == .targetStandardizedMean {
                for component in plan.components {
                    let values = rows.map { $0.scores[component.feature]! }
                    let anchor = values[0]
                    let mean = anchor + values.reduce(0.0) { $0 + ($1 - anchor) } / Double(values.count)
                    let variance = values.reduce(0.0) { $0 + ($1 - mean) * ($1 - mean) } / Double(values.count)
                    let sd = sqrt(variance)
                    try require(mean.isFinite && sd.isFinite, "nonfinite ranking normalization")
                    norms.append(.init(feature: component.feature, mean: mean, populationSD: sd, constant: sd == 0))
                }
            }
            let weightSum = plan.components.reduce(0.0) { $0 + $1.weight }
            var scores: [(String, Double)] = []
            for row in rows {
                let score: Double
                if plan.method == .fixedScore { score = row.scores[plan.components[0].feature]! }
                else {
                    score = plan.components.indices.reduce(0.0) { sum, j in
                        let value = norms[j].constant ? 0 : (row.scores[plan.components[j].feature]! - norms[j].mean) / norms[j].populationSD
                        return sum + plan.components[j].weight / weightSum * value
                    }
                }
                try require(score.isFinite, "nonfinite ranking score")
                scores.append((row.id, score))
            }
            scores.sort { $0.1 == $1.1 ? $0.0 < $1.0 : $0.1 > $1.1 }
            let k = min(plan.topK, scores.count)
            var ranked: [RankedCandidate] = [], start = 0
            while start < scores.count {
                var end = start + 1
                while end < scores.count && scores[end].1 == scores[start].1 { end += 1 }
                let selected = max(0, min(end, k) - start)
                let inclusion = Double(selected) / Double(end - start)
                for i in start..<end { ranked.append(.init(id: scores[i].0, score: scores[i].1, selectionWeight: inclusion)) }
                start = end
            }
            targets.append(.init(target: target, inputCount: pool.count, effectiveK: k,
                normalization: norms, candidates: ranked, excluded: excluded))
        }
        return Report(schemaVersion: 1, methodVersion: "numivivo-fixed-ranking-v1",
            querySHA256: try fingerprint(query), planSHA256: try fingerprint(plan), targets: targets,
            limitations: ["No experimental outcomes, fitting or binding probabilities enter this operation.",
                "Within-target normalization uses all complete-case query candidates, including untested designs; rankings depend on the query pool.",
                "Missing component values exclude a candidate explicitly; no zero imputation or unreported predictor substitution.",
                "Constant standardized components contribute zero. Cutoff ties receive equal fractional inclusion.",
                "This is scoring of supplied predictor outputs, not model-weight inference or biological qualification."])
    }

    /// Only this separate operation sees outcomes. It verifies exact query reconstruction
    /// and recomputes ranking before joining labels; missing outcomes never refill top K.
    public static func assess(_ dataset: VivoBinderBenchmark.Dataset, query: Query,
                              plan: Plan, ranking: Report) throws -> Assessment {
        try require(try fingerprint(self.query(from: dataset)) == fingerprint(query), "ranking query does not reconstruct from experimental import")
        try require(try fingerprint(rank(query, plan: plan)) == fingerprint(ranking), "ranking report does not replay")
        let lookup = Dictionary(uniqueKeysWithValues: dataset.records.map { ($0.id, $0) })
        var results: [TargetAssessment] = []
        for target in ranking.targets {
            var weights: [String: Double] = [:], labels: [Double] = [], scores: [Double] = []
            for candidate in target.candidates {
                let row = lookup[candidate.id]!
                weights[row.outcome.rawValue, default: 0] += candidate.selectionWeight
                if let label = row.outcome.binary { labels.append(label); scores.append(candidate.score) }
            }
            let hits = weights[VivoBinderBenchmark.Outcome.binder.rawValue, default: 0]
            let binarySelected = hits + weights[VivoBinderBenchmark.Outcome.nonBinder.rawValue, default: 0]
            let nonbinary = weights.filter { !["binder", "nonBinder"].contains($0.key) }.values.reduce(0, +)
            let metrics = labels.isEmpty ? nil : try VivoBinderBenchmark.metrics(labels: labels, scores: scores, topK: 1)
            results.append(.init(target: target.target, rankedCount: target.candidates.count, effectiveK: target.effectiveK,
                binaryOutcomeCount: labels.count, selectedWeightByOutcome: weights, measuredBinderHits: hits,
                selectedNonbinaryWeight: nonbinary,
                precisionAmongBinarySelected: binarySelected > 0 ? hits / binarySelected : nil,
                binaryAUROC: metrics?.auroc, binaryAveragePrecision: metrics?.averagePrecision))
        }
        return Assessment(schemaVersion: 1, evidenceStatus: "retrospective-fixed-ranking; not prospective validation",
            assay: dataset.assay, sourceSHA256: dataset.sourceSHA256, rankingSHA256: try fingerprint(ranking), targets: results,
            limitations: ["Selection is fixed before joining outcomes. Unknown/untested/expression-failure outcomes are not non-binders and do not refill top K.",
                "Measured binder hits are known observed selection weight, not total biological success when selected outcomes are unavailable.",
                "Binary discrimination and conditional precision exclude nonbinary outcomes; inspect their coverage separately.",
                "Separate laboratory endpoints measure overlapping designs, not independent biological replicates.",
                "Integrity and label isolation do not prove historical pre-registration, independence or new-target accuracy."])
    }

    private static func validate(_ query: Query) throws {
        try require(query.schemaVersion == 1 && query.sourceSHA256.utf8.count == 64
                    && query.sourceSHA256.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }), "invalid query schema/source digest")
        try require(!query.candidates.isEmpty && query.candidates.count <= 100_000, "query capacity outside 1...100000")
        var ids = Set<String>()
        for row in query.candidates {
            try require(!row.id.isEmpty && row.id.utf8.count <= 1024 && !row.target.isEmpty && row.target.utf8.count <= 1024
                        && ids.insert(row.id).inserted, "invalid or duplicate ranking candidate")
            try require(row.scores.count <= 128, "too many query scores")
            for (name, value) in row.scores {
                try require(VivoBinderAnthropicImport.allowedFeatures.contains(name) && value.isFinite && (0...1).contains(value),
                            "query score is undeclared, nonfinite or outside [0,1]")
            }
        }
    }
    static func fingerprint<T: Encodable>(_ value: T) throws -> String {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(value)).hex
    }
    private static func require(_ condition: Bool, _ reason: String) throws {
        if !condition { throw invalid(reason) }
    }
    private static func invalid(_ reason: String) -> VivoBinderBenchmark.Failure { .invalid(reason) }
}
