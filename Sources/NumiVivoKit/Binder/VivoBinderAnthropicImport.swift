import Foundation

/// Admits the published design_summary.csv without using wet-lab fields as predictors.
/// Source-byte hashing and immutable storage are owned by the calling I/O layer.
public enum VivoBinderAnthropicImport {
    public enum Assay: String, Codable, Sendable { case adaptyv, twist }
    public struct Configuration: Codable, Sendable {
        public let schemaVersion: Int
        public let assay: Assay
        public let targets: [String]
        public let features: [String]
        public init(assay: Assay, targets: [String], features: [String]) {
            schemaVersion = 1; self.assay = assay; self.targets = targets; self.features = features
        }
    }
    public struct Result: Codable, Sendable {
        public let dataset: VivoBinderBenchmark.Dataset
        public let sourceRowByID: [String: Int]
        public let excludedTargetCounts: [String: Int]
        public let unknownOutcomeCounts: [String: Int]
    }
    public static let repository = "Anthropic/claude-protein-binder-design"
    public static let revision = "9e1b81696da46835e9e9cde9a3da976e0abc92ab"
    public static let sourcePath = "data/tables/design_summary.csv"
    public static let sourceGitBlobSHA1 = "d1573ba03e8322c70ccb3a40e46418e86b40e2dc"
    public static let allowedFeatures: Set<String> = Set(
        ["ef2fast", "ef2full", "ptxv2", "odde", "afm3", "boltz2", "chai1", "of3", "rf3", "af3of3"]
            .flatMap { ["ipsae_min_\($0)", "sc_dockq_\($0)"] })

    public static func parse(_ data: Data, sourceSHA256: String,
                             configuration: Configuration) throws -> Result {
        typealias B = VivoBinderBenchmark
        let c = configuration
        guard c.schemaVersion == 1, !c.targets.isEmpty, Set(c.targets).count == c.targets.count,
              !c.features.isEmpty, Set(c.features).count == c.features.count,
              Set(c.features).isSubset(of: allowedFeatures) else {
            throw B.Failure.invalid("invalid import config; targets explicit and predictors restricted to in-silico scores")
        }
        let table = try csv(data)
        guard let header = table.first, Set(header).count == header.count,
              !header.contains("") else { throw B.Failure.invalid("empty or duplicate CSV header") }
        let binding = "\(c.assay.rawValue)_binding"
        let needed = Set(["uuid", "target", "sequence", binding] + c.features)
        guard needed.isSubset(of: Set(header)) else { throw B.Failure.invalid("missing required CSV columns") }
        var rows: [B.Record] = [], rowNumbers: [String: Int] = [:]
        var excluded: [String: Int] = [:], unknown: [String: Int] = [:], found = Set<String>()
        let targets = Set(c.targets), alphabet = Set("ACDEFGHIKLMNPQRSTVWY")
        for (offset, values) in table.dropFirst().enumerated() {
            guard values.count == header.count else { throw B.Failure.invalid("CSV row width at logical record \(offset + 2)") }
            let fields = Dictionary(uniqueKeysWithValues: zip(header, values))
            let target = fields["target"]!
            guard targets.contains(target) else { excluded[target, default: 0] += 1; continue }
            found.insert(target)
            let sequence = fields["sequence"]!, id = fields["uuid"]!, raw = fields[binding]!
            guard !sequence.isEmpty, sequence.count <= 5000, sequence.allSatisfy({ alphabet.contains($0) }) else {
                throw B.Failure.invalid("unsupported/missing amino-acid sequence: \(id)")
            }
            let outcome: B.Outcome
            switch raw {
            case "binder": outcome = .binder
            case "non_binder": outcome = .nonBinder
            case "not_tested": outcome = .notTested
            case "expression_failure": outcome = .expressionFailure
            default: outcome = .inconclusive; unknown[raw, default: 0] += 1
            }
            var features: [String: Double] = [:]
            for name in c.features {
                let text = fields[name]!
                if text.isEmpty { continue }
                guard let number = Double(text), number.isFinite, (0...1).contains(number) else {
                    throw B.Failure.invalid("invalid normalized in-silico score \(name): \(id)")
                }
                features[name] = number
            }
            rows.append(B.Record(id: id, target: target, leakageGroup: "exact-sequence:\(sequence)",
                outcome: outcome, rawOutcome: raw.isEmpty ? "<missing>" : raw,
                features: features, sourceFields: fields))
            rowNumbers[id] = offset + 2
        }
        guard found == targets else { throw B.Failure.invalid("requested targets absent from source") }
        let dataset = B.Dataset(sourceSHA256: sourceSHA256, assay: "\(c.assay.rawValue):reported-binding-class",
            groupingMethod: "exact-sequence-only; near-homology not controlled", records: rows.sorted { $0.id < $1.id })
        try B.validate(dataset)
        return Result(dataset: dataset, sourceRowByID: rowNumbers, excludedTargetCounts: excluded, unknownOutcomeCounts: unknown)
    }

    /// Strict bounded UTF-8 RFC-4180-style reader (LF/CRLF, quoted newlines, doubled quotes).
    static func csv(_ data: Data) throws -> [[String]] {
        typealias F = VivoBinderBenchmark.Failure
        guard !data.isEmpty, data.count <= 64 * 1024 * 1024, String(data: data, encoding: .utf8) != nil else {
            throw F.invalid("CSV must be UTF-8 and at most 64 MiB")
        }
        let bytes = Array(data), start = data.starts(with: [0xef, 0xbb, 0xbf]) ? 3 : 0
        var i = start, quoted = false, closed = false, field: [UInt8] = [], row: [String] = [], rows: [[String]] = []
        func finishField() throws {
            guard field.count <= 65_536, row.count < 128 else { throw F.invalid("CSV field/column capacity") }
            row.append(String(decoding: field, as: UTF8.self)); field.removeAll(keepingCapacity: true); closed = false
        }
        func finishRow() throws {
            try finishField()
            guard rows.count <= 100_000 else { throw F.invalid("CSV record capacity") }
            rows.append(row); row.removeAll(keepingCapacity: true)
        }
        while i < bytes.count {
            let b = bytes[i]
            if quoted {
                if b == 34 {
                    if i + 1 < bytes.count && bytes[i + 1] == 34 { field.append(34); i += 1 }
                    else { quoted = false; closed = true }
                } else { field.append(b) }
            } else if b == 44 { try finishField() }
            else if b == 10 || b == 13 {
                if b == 13 {
                    guard i + 1 < bytes.count && bytes[i + 1] == 10 else { throw F.invalid("bare CR outside quoted CSV field") }
                    i += 1
                }
                try finishRow()
            } else if closed { throw F.invalid("characters after closing CSV quote") }
            else if b == 34 {
                guard field.isEmpty else { throw F.invalid("quote inside unquoted CSV field") }
                quoted = true
            } else { field.append(b) }
            guard field.count <= 65_536 else { throw F.invalid("CSV field capacity") }
            i += 1
        }
        guard !quoted else { throw F.invalid("unterminated quoted CSV field") }
        if !field.isEmpty || !row.isEmpty || closed { try finishRow() }
        return rows
    }
}
