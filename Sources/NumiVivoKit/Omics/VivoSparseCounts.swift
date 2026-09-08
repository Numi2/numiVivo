import Foundation

public enum VivoOmicsError: Error, LocalizedError, Sendable {
    case invalid(String)
    case limit(String)
    public var errorDescription: String? {
        switch self {
        case .invalid(let value): return "omics: \(value)"
        case .limit(let value): return "omics resource limit: \(value)"
        }
    }
}

struct VivoOmicsJSONKey: CodingKey {
    let stringValue: String
    let intValue: Int? = nil
    init?(stringValue: String) { self.stringValue = stringValue }
    init?(intValue: Int) { return nil }
}
func vivoOmicsRejectUnknownKeys(_ decoder: Decoder, allowed: Set<String>) throws {
    let fields = try decoder.container(keyedBy: VivoOmicsJSONKey.self)
    let unknown = Set(fields.allKeys.map(\.stringValue)).subtracting(allowed)
    guard unknown.isEmpty else {
        throw DecodingError.dataCorrupted(.init(codingPath: decoder.codingPath,
            debugDescription: "Unknown omics fields: " + unknown.sorted().joined(separator: ", ")))
    }
}

/// Limits are checked before allocating from untrusted dimensions. This profile
/// is bounded in-memory processing, not an out-of-core or atlas-scale claim.
public struct VivoOmicsLimits: Codable, Sendable, Equatable {
    public var maximumCells: Int = 100_000
    public var maximumFeatures: Int = 100_000
    public var maximumNonzeros: Int = 2_000_000
    public var maximumInputBytes: Int = 64 * 1_024 * 1_024
    public var maximumLineBytes: Int = 16_384
    private enum CodingKeys: String, CodingKey { case maximumCells, maximumFeatures, maximumNonzeros, maximumInputBytes, maximumLineBytes }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["maximumCells", "maximumFeatures", "maximumNonzeros", "maximumInputBytes", "maximumLineBytes"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        maximumCells = try values.decodeIfPresent(Int.self, forKey: .maximumCells) ?? 100_000
        maximumFeatures = try values.decodeIfPresent(Int.self, forKey: .maximumFeatures) ?? 100_000
        maximumNonzeros = try values.decodeIfPresent(Int.self, forKey: .maximumNonzeros) ?? 2_000_000
        maximumInputBytes = try values.decodeIfPresent(Int.self, forKey: .maximumInputBytes) ?? 64 * 1_024 * 1_024
        maximumLineBytes = try values.decodeIfPresent(Int.self, forKey: .maximumLineBytes) ?? 16_384
    }
    public init() {}
    public func validate() throws {
        guard maximumCells > 0, maximumCells < Int.max,
              maximumFeatures > 0, maximumNonzeros >= 0,
              maximumInputBytes > 0, maximumInputBytes < Int.max,
              maximumLineBytes > 0 else { throw VivoOmicsError.invalid("nonpositive or overflowing limits") }
    }
}

func vivoOmicsID(_ value: String) -> Bool {
    !value.isEmpty && value.utf8.count <= 1_024 &&
    value == value.trimmingCharacters(in: .whitespacesAndNewlines) &&
    !value.unicodeScalars.contains(where: { CharacterSet.controlCharacters.contains($0) })
}
func vivoOmicsSum(_ a: UInt64, _ b: UInt64) throws -> UInt64 {
    let result = a.addingReportingOverflow(b)
    guard !result.overflow else { throw VivoOmicsError.invalid("raw count accumulation exceeds UInt64") }
    return result.partialValue
}

/// Measurement units describe assay counts, not a claim of absolute molecules.
public enum VivoOmicsCountUnit: String, Codable, Sendable { case umiCount, readCount }
public enum VivoOmicsEvidence: String, Codable, Sendable { case measured, synthetic, simulated }

public struct VivoOmicsSample: Codable, Sendable, Equatable {
    public let id: String
    public let biologicalReplicateID: String
    public let donorID: String?
    public let condition: String
    public let batchID: String
    public let organism: String
    private enum CodingKeys: String, CodingKey { case id, biologicalReplicateID, donorID, condition, batchID, organism }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["id", "biologicalReplicateID", "donorID", "condition", "batchID", "organism"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        id = try values.decode(String.self, forKey: .id)
        biologicalReplicateID = try values.decode(String.self, forKey: .biologicalReplicateID)
        donorID = try values.decodeIfPresent(String.self, forKey: .donorID)
        condition = try values.decode(String.self, forKey: .condition)
        batchID = try values.decode(String.self, forKey: .batchID)
        organism = try values.decode(String.self, forKey: .organism)
    }
    public init(id: String, biologicalReplicateID: String, donorID: String? = nil,
                condition: String, batchID: String, organism: String) {
        self.id = id; self.biologicalReplicateID = biologicalReplicateID; self.donorID = donorID
        self.condition = condition; self.batchID = batchID; self.organism = organism
    }
    public func validate() throws {
        guard [id, biologicalReplicateID, condition, batchID, organism].allSatisfy(vivoOmicsID),
              donorID.map(vivoOmicsID) ?? true else { throw VivoOmicsError.invalid("sample identity or metadata") }
    }
}
public struct VivoOmicsFeature: Codable, Sendable, Equatable {
    public let id: String
    public let name: String
    /// Explicit annotation; no species-dependent gene-name heuristic.
    public let mitochondrial: Bool
    public init(id: String, name: String, mitochondrial: Bool = false) {
        self.id = id; self.name = name; self.mitochondrial = mitochondrial
    }
}
public struct VivoOmicsCell: Codable, Sendable, Equatable, Hashable {
    public let barcode: String
    public let sampleID: String
    /// Supplied annotation, never an inferred cell type in this implementation.
    public let group: String?
    public init(barcode: String, sampleID: String, group: String? = nil) {
        self.barcode = barcode; self.sampleID = sampleID; self.group = group
    }
}

/// Canonical CSR: cells are rows, features are columns, sorted unique positive
/// entries per row. Raw counts never pass through a floating-point representation.
public struct VivoSparseCounts: Codable, Sendable, Equatable {
    public let cellCount: Int
    public let featureCount: Int
    public let rowOffsets: [Int]
    public let featureIndices: [Int]
    public let counts: [UInt64]
    public init(cellCount: Int, featureCount: Int, rowOffsets: [Int], featureIndices: [Int], counts: [UInt64]) {
        self.cellCount = cellCount; self.featureCount = featureCount
        self.rowOffsets = rowOffsets; self.featureIndices = featureIndices; self.counts = counts
    }
    public func validate(limits: VivoOmicsLimits = .init()) throws {
        try limits.validate()
        guard cellCount >= 0, cellCount <= limits.maximumCells, featureCount > 0,
              featureCount <= limits.maximumFeatures, counts.count <= limits.maximumNonzeros else {
            throw VivoOmicsError.limit("sparse matrix shape or nonzero count")
        }
        guard rowOffsets.count == cellCount + 1, rowOffsets.first == 0,
              rowOffsets.last == counts.count, featureIndices.count == counts.count else {
            throw VivoOmicsError.invalid("CSR array lengths or endpoints")
        }
        for row in 0..<cellCount {
            let start = rowOffsets[row], end = rowOffsets[row + 1]
            guard start >= 0, end >= start, end <= counts.count else { throw VivoOmicsError.invalid("CSR row offsets") }
            var previous = -1, total: UInt64 = 0
            for k in start..<end {
                let column = featureIndices[k]
                guard column > previous, column < featureCount, counts[k] > 0 else {
                    throw VivoOmicsError.invalid("CSR indices must be sorted, unique and in range; counts positive")
                }
                total = try vivoOmicsSum(total, counts[k]); previous = column
            }
        }
    }
}

/// This is a count-assay payload for the existing artifact store, not another
/// artifact/provenance database. Wrappers retain original sources by fingerprint.
public struct VivoSingleCellDataset: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let countUnit: VivoOmicsCountUnit
    public let samples: [VivoOmicsSample]
    public let features: [VivoOmicsFeature]
    public let cells: [VivoOmicsCell]
    public let matrix: VivoSparseCounts
    public init(id: String, evidence: VivoOmicsEvidence, sourceDescription: String,
                countUnit: VivoOmicsCountUnit, samples: [VivoOmicsSample], features: [VivoOmicsFeature],
                cells: [VivoOmicsCell], matrix: VivoSparseCounts) {
        self.schemaVersion = 1; self.id = id; self.evidence = evidence; self.sourceDescription = sourceDescription
        self.countUnit = countUnit; self.samples = samples; self.features = features; self.cells = cells; self.matrix = matrix
    }
    public func validate(limits: VivoOmicsLimits = .init()) throws {
        try matrix.validate(limits: limits)
        guard schemaVersion == 1, vivoOmicsID(id), !sourceDescription.isEmpty,
              sourceDescription.utf8.count <= 16_384, !samples.isEmpty,
              samples.count <= limits.maximumCells, cells.count == matrix.cellCount,
              features.count == matrix.featureCount else { throw VivoOmicsError.invalid("dataset schema or shape") }
        var sampleByID: [String: VivoOmicsSample] = [:]
        var replicateByID: [String: VivoOmicsSample] = [:]
        for sample in samples {
            try sample.validate()
            guard sampleByID.updateValue(sample, forKey: sample.id) == nil else { throw VivoOmicsError.invalid("duplicate sample ID") }
            if let old = replicateByID[sample.biologicalReplicateID] {
                guard old.donorID == sample.donorID, old.organism == sample.organism else {
                    throw VivoOmicsError.invalid("replicate maps to inconsistent donor or organism")
                }
            }
            replicateByID[sample.biologicalReplicateID] = sample
        }
        var featureIDs = Set<String>()
        for feature in features {
            guard vivoOmicsID(feature.id), vivoOmicsID(feature.name), featureIDs.insert(feature.id).inserted else {
                throw VivoOmicsError.invalid("invalid or duplicate feature ID; duplicate names are permitted")
            }
        }
        struct Identity: Hashable { let sample: String; let barcode: String }
        var identities = Set<Identity>()
        for cell in cells {
            guard vivoOmicsID(cell.barcode), sampleByID[cell.sampleID] != nil,
                  cell.group.map(vivoOmicsID) ?? true,
                  identities.insert(.init(sample: cell.sampleID, barcode: cell.barcode)).inserted else {
                throw VivoOmicsError.invalid("cell identity, sample reference or group")
            }
        }
    }
}
