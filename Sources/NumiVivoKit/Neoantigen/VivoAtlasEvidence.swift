import Foundation

public enum VivoGenomicEvidenceError: Error, Sendable, CustomStringConvertible {
    case invalid(String)
    public var description: String { switch self { case .invalid(let message): message } }
}

/// Shared between external captures and the native importer. This is a declared
/// case binding, not genetic sample authentication or proof of sequencing QC.
public struct VivoGenomicCaseBinding: Codable, Sendable, Equatable {
    public var parentReportSHA256: String
    public var caseFingerprint: String
    public var dataClass: String
    public var assembly: String
    public var referenceSHA256: String
    public var hlaAlleles: [String]
    public var tumorRNASampleID: String
    public init(parentReportSHA256: String, caseFingerprint: String, dataClass: String,
                assembly: String, referenceSHA256: String, hlaAlleles: [String], tumorRNASampleID: String) {
        self.parentReportSHA256 = parentReportSHA256; self.caseFingerprint = caseFingerprint
        self.dataClass = dataClass; self.assembly = assembly; self.referenceSHA256 = referenceSHA256
        self.hlaAlleles = hlaAlleles; self.tumorRNASampleID = tumorRNASampleID
    }
}

/// Values are preserved as JSON, not coerced to a common measurement or unit.
public indirect enum VivoGenomicJSON: Codable, Sendable, Equatable {
    case null, bool(Bool), number(Double), string(String)
    case array([VivoGenomicJSON]), object([String: VivoGenomicJSON])
    public init(from decoder: Decoder) throws {
        let c = try decoder.singleValueContainer()
        if c.decodeNil() { self = .null }
        else if let v = try? c.decode(Bool.self) { self = .bool(v) }
        else if let v = try? c.decode(Double.self) { self = .number(v) }
        else if let v = try? c.decode(String.self) { self = .string(v) }
        else if let v = try? c.decode([VivoGenomicJSON].self) { self = .array(v) }
        else { self = .object(try c.decode([String: VivoGenomicJSON].self)) }
    }
    public func encode(to encoder: Encoder) throws {
        var c = encoder.singleValueContainer()
        switch self {
        case .null: try c.encodeNil()
        case .bool(let v): try c.encode(v)
        case .number(let v): try c.encode(v)
        case .string(let v): try c.encode(v)
        case .array(let v): try c.encode(v)
        case .object(let v): try c.encode(v)
        }
    }
}

public struct VivoAtlasVariantInput: Codable, Sendable, Equatable {
    public var candidateID: String
    public var chromosome: String
    public var start0: Int
    public var stop0: Int
    public var reference: String
    public var alternate: String
    public init(candidateID: String, chromosome: String, start0: Int, stop0: Int, reference: String, alternate: String) {
        self.candidateID = candidateID; self.chromosome = chromosome; self.start0 = start0; self.stop0 = stop0
        self.reference = reference; self.alternate = alternate
    }
}
public struct VivoAtlasQuery: Codable, Sendable, Equatable {
    public var key: String
    public var chromosome: String
    public var position1: Int
    public var reference: String
    public var alternate: String
    public var candidateIDs: [String]
}
public struct VivoAtlasRequest: Codable, Sendable, Equatable {
    public var schema: String
    public var binding: VivoGenomicCaseBinding
    public var requestedScorers: [String]
    public var ontologyTerms: [String]
    public var variants: [VivoAtlasQuery]
    public var excludedCandidates: [String: String]
}
public struct VivoAtlasObservation: Codable, Sendable, Equatable {
    public var variantKey: String
    public var metadata: [String: VivoGenomicJSON]
}
public struct VivoAtlasTable: Codable, Sendable, Equatable {
    public var scorer: String
    public var isSigned: Bool
    public var observations: [VivoAtlasObservation]
    public var tracks: [[String: VivoGenomicJSON]]
    public var rawScores: [[Double?]]
    public var quantiles: [[Double?]]?
}
public struct VivoAtlasQueryOutcome: Codable, Sendable, Equatable {
    public enum Status: String, Codable, Sendable { case available, noData, failed }
    public var variantKey: String
    public var status: Status
    public var diagnostic: String?
}
public struct VivoAtlasCapture: Codable, Sendable, Equatable {
    public var schema: String
    /// SHA-256 of exact input request bytes, not a cross-language JSON re-encoding.
    public var requestSHA256: String
    public var mode: String
    public var provider: String
    public var clientCommit: String
    public var clientVersion: String
    public var adapterSHA256: String
    public var retrievedAt: String
    /// Atlas currently exposes no immutable model/dataset version selector.
    /// Capture bytes are pinned; an invented model version must not be substituted.
    public var serviceVersion: String?
    public var serviceVersionStatus: String
    public var termsURI: String
    public var permittedUse: String
    public var trainingPermitted: Bool
    public var referenceCheck: String
    public var referenceSHA256: String
    public var outcomes: [VivoAtlasQueryOutcome]
    public var tables: [VivoAtlasTable]
}

public enum VivoAtlasEvidence {
    public static let clientCommit = "aa6fc8f6faadcb8c910fa2b85b57386fbd5c7b5d"
    public static let maximumDocumentBytes = 32 * 1024 * 1024
    public static let maximumCells = 250_000
    public static let maximumVariants = 128
    public static let limitations = [
        "Atlas molecular-impact predictions are not evidence of antigen presentation, immune recognition or clinical benefit.",
        "Signed quantiles retain their direction; a quantile is not a probability of vaccine effectiveness.",
        "AVI and its AlphaGenome/AlphaMissense components are dependent evidence, not independent confirmations.",
        "No tissue match, service response, or score is substituted for missing measured tumor RNA.",
        "The live Atlas service has no model-version selector in the pinned client; archived responses, not future re-queries, are reproducible.",
        "Research output usage restrictions remain attached to the capture. No training-data export is implemented."
    ]
    static func require(_ condition: Bool, _ message: String) throws {
        guard condition else { throw VivoGenomicEvidenceError.invalid(message) }
    }
    static func matches(_ text: String, _ expression: String) -> Bool {
        text.range(of: expression, options: .regularExpression) != nil
    }
    public static func digest(_ bytes: Data) throws -> String { try VivoCanonicalJSON.fingerprint(bytes).hex }
    static func isDigest(_ text: String) -> Bool { matches(text, "^[a-f0-9]{64}$") }
    static func validToken(_ text: String) -> Bool { matches(text, "^[A-Za-z0-9][A-Za-z0-9_.:-]{0,127}$") }
    public static func chromosome(_ text: String) -> String? {
        let value = text.hasPrefix("chr") ? String(text.dropFirst(3)) : text
        guard matches(value, "^([1-9]|1[0-9]|2[0-2]|X|Y)$") else { return nil }
        return "chr" + value
    }
    public static func validate(_ binding: VivoGenomicCaseBinding) throws {
        try require(isDigest(binding.parentReportSHA256) && isDigest(binding.caseFingerprint) && isDigest(binding.referenceSHA256), "Invalid case digest binding.")
        try require(["synthetic", "publicReference"].contains(binding.dataClass), "Only synthetic and public-reference research cases are supported.")
        try require(["GRCh37", "GRCh38"].contains(binding.assembly) && validToken(binding.tumorRNASampleID), "Invalid assembly or RNA sample binding.")
        try require(!binding.hlaAlleles.isEmpty && binding.hlaAlleles.count <= 6 && Set(binding.hlaAlleles).count == binding.hlaAlleles.count && binding.hlaAlleles.allSatisfy { matches($0, "^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}$") }, "This profile requires explicit class-I HLA alleles.")
    }
    /// Only exact single-base substitutions are converted. No liftover, indel
    /// normalization, reverse complementation or contig guessing is performed.
    public static func plan(binding: VivoGenomicCaseBinding, inputs: [VivoAtlasVariantInput],
                            requestedScorers: [String], ontologyTerms: [String]) throws -> VivoAtlasRequest {
        try validate(binding)
        try require(inputs.count <= 10_000 && Set(inputs.map(\.candidateID)).count == inputs.count, "Candidate identifiers must be unique and bounded.")
        try require(!requestedScorers.isEmpty && requestedScorers.count <= 16 && Set(requestedScorers).count == requestedScorers.count && requestedScorers.allSatisfy { validToken($0) }, "Specify exact, unique Atlas scorer names from a metadata snapshot.")
        try require(ontologyTerms.count <= 32 && Set(ontologyTerms).count == ontologyTerms.count && ontologyTerms.allSatisfy { matches($0, "^(UBERON|CL|CLO|EFO):[0-9]+$") }, "Unsupported or duplicate ontology selector.")
        var queries: [String: VivoAtlasQuery] = [:], excluded: [String: String] = [:]
        for input in inputs {
            try require(isDigest(input.candidateID), "Invalid candidate identifier.")
            guard binding.assembly == "GRCh38" else { excluded[input.candidateID] = "unsupportedAssembly: GRCh38 is required; no automatic liftover"; continue }
            guard let chrom = chromosome(input.chromosome), input.start0 >= 0, input.start0 < 300_000_000,
                  input.stop0 == input.start0 + 1, matches(input.reference, "^[ACGT]$"),
                  matches(input.alternate, "^[ACGT]$"), input.reference != input.alternate else {
                excluded[input.candidateID] = "unsupportedVariant: exact SNV on chr1–22/X/Y required"; continue
            }
            let position = input.start0 + 1
            let key = "GRCh38|\(chrom)|\(position)|\(input.reference)|\(input.alternate)"
            if queries[key] == nil { queries[key] = .init(key: key, chromosome: chrom, position1: position, reference: input.reference, alternate: input.alternate, candidateIDs: []) }
            queries[key]!.candidateIDs.append(input.candidateID)
        }
        try require(queries.count <= maximumVariants, "This bounded request supports at most 128 distinct SNVs.")
        let variants = queries.keys.sorted().map { key in
            var q = queries[key]!; q.candidateIDs.sort(); return q
        }
        return .init(schema: "numivivo.org/atlas-request/v1", binding: binding,
            requestedScorers: requestedScorers.sorted(), ontologyTerms: ontologyTerms.sorted(), variants: variants, excludedCandidates: excluded)
    }
    /// Replanning checks every key, coordinate and duplicate candidate mapping.
    public static func validate(_ request: VivoAtlasRequest) throws {
        try require(request.schema == "numivivo.org/atlas-request/v1", "Unsupported Atlas request schema.")
        try require(request.variants.count <= maximumVariants && request.excludedCandidates.count <= 10_000, "Atlas request exceeds limits.")
        var inputs: [VivoAtlasVariantInput] = [], seen = Set<String>()
        for q in request.variants {
            try require(q.candidateIDs.count > 0 && q.candidateIDs.count <= 10_000 && q.position1 > 0 && q.position1 <= 300_000_000, "Invalid query size or one-based position.")
            for id in q.candidateIDs {
                try require(seen.insert(id).inserted && isDigest(id), "Duplicate/invalid candidate binding.")
                inputs.append(.init(candidateID: id, chromosome: q.chromosome, start0: q.position1 - 1, stop0: q.position1, reference: q.reference, alternate: q.alternate))
            }
        }
        try require(request.excludedCandidates.allSatisfy { isDigest($0.key) && !seen.contains($0.key) && !$0.value.isEmpty && $0.value.utf8.count <= 512 }, "Invalid unsupported-candidate record.")
        let generated = try plan(binding: request.binding, inputs: inputs, requestedScorers: request.requestedScorers, ontologyTerms: request.ontologyTerms)
        try require(generated.variants == request.variants && generated.excludedCandidates.isEmpty && generated.requestedScorers == request.requestedScorers && generated.ontologyTerms == request.ontologyTerms, "Request keys, alleles, ordering or selectors do not reconstruct.")
    }
    static func validateMetadata(_ metadata: [String: VivoGenomicJSON]) throws {
        func check(_ value: VivoGenomicJSON, depth: Int) throws {
            try require(depth <= 8, "Metadata nesting exceeds limits.")
            switch value {
            case .string(let s): try require(s.utf8.count <= 4096 && !s.contains("\u{0}"), "Oversized metadata value.")
            case .number(let n): try require(n.isFinite, "Nonfinite metadata number.")
            case .array(let a): try require(a.count <= 128, "Oversized metadata array."); for v in a { try check(v, depth: depth + 1) }
            case .object(let o): try require(o.count <= 128, "Oversized metadata object."); for (k, v) in o { try require(k.utf8.count <= 256, "Oversized metadata key."); try check(v, depth: depth + 1) }
            default: break
            }
        }
        try check(.object(metadata), depth: 0)
    }
    public static func validate(capture: VivoAtlasCapture, request: VivoAtlasRequest, requestData: Data) throws {
        try require(requestData.count <= 1024 * 1024, "Atlas request exceeds 1 MiB.")
        try require(try VivoCanonicalJSON.decode(VivoAtlasRequest.self, from: requestData) == request, "Request bytes do not encode the supplied request.")
        try validate(request)
        try require(capture.schema == "numivivo.org/atlas-capture/v1" && capture.requestSHA256 == (try digest(requestData)), "Atlas capture belongs to another request.")
        try require(capture.provider == "google-deepmind/alphagenome-atlas" && capture.clientCommit == clientCommit && validToken(capture.clientVersion) && isDigest(capture.adapterSHA256), "Unqualified Atlas client/provider identity.")
        try require(ISO8601DateFormatter().date(from: capture.retrievedAt) != nil, "Invalid retrieval timestamp.")
        try require(capture.serviceVersion == nil && capture.serviceVersionStatus == "not-exposed-by-pinned-api", "Do not invent a version for the unversioned Atlas endpoint.")
        try require(capture.permittedUse == "noncommercialResearch" && !capture.trainingPermitted && capture.termsURI == "https://deepmind.google.com/science/alphagenome/terms", "Unsupported output-usage contract.")
        try require(capture.referenceSHA256 == request.binding.referenceSHA256, "Reference identity changed during retrieval.")
        try require((capture.mode == "live" && request.binding.dataClass == "publicReference" && capture.referenceCheck == "local-fasta-sha256-and-ref-bases") || (capture.mode == "fixture" && request.binding.dataClass == "synthetic" && capture.referenceCheck == "synthetic-not-biological"), "Live/fixture or reference-check provenance mismatch.")
        let keys = Set(request.variants.map(\.key)), outcomes = capture.outcomes
        try require(outcomes.count == keys.count && Set(outcomes.map(\.variantKey)) == keys && Set(outcomes.map(\.variantKey)).count == outcomes.count, "Each requested variant needs exactly one explicit outcome.")
        try require(capture.tables.count <= maximumVariants * request.requestedScorers.count, "Too many scorer tables.")
        let statuses = Dictionary(uniqueKeysWithValues: outcomes.map { ($0.variantKey, $0.status) })
        var observedKeys = Set<String>(), cellKeys = Set<String>(), cells = 0
        for outcome in outcomes {
            try require((outcome.diagnostic?.utf8.count ?? 0) <= 512, "Oversized query diagnostic.")
            if outcome.status == .failed { try require(!(outcome.diagnostic ?? "").isEmpty, "A failed query must not be represented as an empty scientific result.") }
        }
        for table in capture.tables {
            try require(request.requestedScorers.contains(table.scorer) && !table.observations.isEmpty && !table.tracks.isEmpty && table.observations.count <= 10_000 && table.tracks.count <= 10_000, "Unknown/empty/oversized scorer table.")
            let product = table.observations.count.multipliedReportingOverflow(by: table.tracks.count)
            try require(!product.overflow && product.partialValue <= maximumCells - cells, "Capture exceeds the 250,000-cell budget.")
            cells += product.partialValue
            try require(table.rawScores.count == table.observations.count && (table.quantiles == nil || table.quantiles!.count == table.observations.count), "Score matrix dimensions disagree.")
            for metadata in table.tracks { try validateMetadata(metadata) }
            let trackDigests = try table.tracks.map { try digest(VivoCanonicalJSON.encode($0)) }
            try require(Set(trackDigests).count == trackDigests.count, "Duplicate track identities.")
            for (i, observation) in table.observations.enumerated() {
                try require(statuses[observation.variantKey] == .available, "Foreign or unavailable variant has scores.")
                try validateMetadata(observation.metadata)
                observedKeys.insert(observation.variantKey)
                let observationID = try digest(VivoCanonicalJSON.encode(observation.metadata))
                try require(table.rawScores[i].count == table.tracks.count && (table.quantiles == nil || table.quantiles![i].count == table.tracks.count), "Score row dimensions disagree.")
                for j in table.tracks.indices {
                    let identity = [observation.variantKey, table.scorer, observationID, trackDigests[j]].joined(separator: "|")
                    try require(cellKeys.insert(identity).inserted, "Duplicate variant/scorer/observation/track cell.")
                    if let raw = table.rawScores[i][j] { try require(raw.isFinite && (table.isSigned || raw >= 0), "Invalid raw score for scorer sign metadata.") }
                    if let quantile = table.quantiles?[i][j] {
                        try require(quantile.isFinite && quantile <= 1 && quantile >= (table.isSigned ? -1 : 0), "Quantile outside its signed/unsigned range.")
                        try require(table.rawScores[i][j] != nil, "A calibrated score cannot conceal an absent raw score.")
                    }
                }
            }
        }
        try require(observedKeys == Set(outcomes.filter { $0.status == .available }.map(\.variantKey)), "Available status requires an actual score table; missing is not zero.")
    }
}
