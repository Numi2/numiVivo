import Foundation

public enum VivoVisiumPositionsFormat: String, Codable, Sendable {
    case legacyHeaderlessCSV, csvWithHeader
    var filename: String {
        self == .legacyHeaderlessCSV ? "tissue_positions_list.csv" : "tissue_positions.csv"
    }
}

/// Standard filtered Visium HDF5 plus an explicitly selected Space Ranger CSV.
/// Image files and scale factors are not required for full-resolution pixels.
public struct VivoVisiumPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let counts: VivoTenXMultiAssayPlan
    public let positionsFormat: VivoVisiumPositionsFormat
    public let frame: VivoSpatialFrame
    public init(schemaVersion: Int, counts: VivoTenXMultiAssayPlan, positionsFormat: VivoVisiumPositionsFormat, frame: VivoSpatialFrame) {
        self.schemaVersion = schemaVersion; self.counts = counts
        self.positionsFormat = positionsFormat; self.frame = frame
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, counts, positionsFormat, frame }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "counts", "positionsFormat", "frame"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decode(Int.self, forKey: .schemaVersion)
        counts = try v.decode(VivoTenXMultiAssayPlan.self, forKey: .counts)
        positionsFormat = try v.decode(VivoVisiumPositionsFormat.self, forKey: .positionsFormat)
        frame = try v.decode(VivoSpatialFrame.self, forKey: .frame)
    }
}

public enum VivoMultiAssayVisium {
    static let maximumPositionBytes = 16_777_216
    struct Position: Equatable {
        let inTissue: Bool
        let arrayRow: UInt64
        let arrayColumn: UInt64
        let x: Double
        let y: Double
    }
    public static func validate(_ plan: VivoVisiumPlan) throws {
        try VivoMultiAssayTenX.validate(plan.counts)
        guard plan.schemaVersion == 1, plan.counts.assays.count == 1,
              plan.counts.assays[0].kind == .rna, plan.counts.assays[0].countUnit == .umiCount,
              vivoOmicsID(plan.frame.id), plan.frame.unit == .pixel, plan.frame.axes == ["x", "y"],
              !plan.frame.sourceDescription.isEmpty, plan.frame.sourceDescription.utf8.count <= 16_384 else {
            throw VivoOmicsError.invalid("Visium requires one RNA UMI assay and an x,y full-resolution pixel frame")
        }
    }
    /// Parse only the publisher's six-column numeric format, never arbitrary CSV
    /// or inferred headers. Retain outside-tissue rows in the source snapshot.
    static func positions(_ bytes: Data, format: VivoVisiumPositionsFormat) throws -> [String: Position] {
        guard bytes.count <= maximumPositionBytes, let text = String(data: bytes, encoding: .utf8) else {
            throw VivoOmicsError.limit("Visium position bytes or UTF-8")
        }
        var lines = text.components(separatedBy: "\n")
        if lines.last == "" { lines.removeLast() }
        guard !lines.isEmpty, lines.count <= VivoMultiAssayDataset.limits.maximumCells + 1 else {
            throw VivoOmicsError.limit("Visium position rows")
        }
        func fields(_ line: String) -> [Substring] {
            let clean = line.hasSuffix("\r") ? line.dropLast() : line[...]
            return clean.split(separator: ",", omittingEmptySubsequences: false)
        }
        if format == .csvWithHeader {
            guard fields(lines.removeFirst()).map(String.init) == ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"] else {
                throw VivoOmicsError.invalid("Visium position header")
            }
        }
        guard !lines.isEmpty, lines.count <= VivoMultiAssayDataset.limits.maximumCells else { throw VivoOmicsError.limit("Visium position rows") }
        var result: [String: Position] = [:]
        for line in lines {
            try Task.checkCancellation()
            let v = fields(line)
            guard v.count == 6, vivoOmicsID(String(v[0])), result[String(v[0])] == nil,
                  v[1] == "0" || v[1] == "1" else { throw VivoOmicsError.invalid("Visium position fields or duplicate barcode") }
            func number(_ i: Int) throws -> UInt64 {
                guard !v[i].isEmpty, v[i].utf8.allSatisfy({ $0 >= 48 && $0 <= 57 }),
                      let value = UInt64(v[i]), value <= 9_007_199_254_740_992 else {
                    throw VivoOmicsError.invalid("Visium coordinates must be nonnegative exact integers")
                }
                return value
            }
            result[String(v[0])] = try .init(inTissue: v[1] == "1", arrayRow: number(2), arrayColumn: number(3), x: Double(number(5)), y: Double(number(4)))
        }
        return result
    }
    static func attach(_ positions: [String: Position], to data: VivoMultiAssayDataset, frame: VivoSpatialFrame) throws -> VivoMultiAssayDataset {
        let expected = Set(positions.filter { $0.value.inTissue }.keys)
        guard expected == Set(data.observations.map { $0.identity.barcode }) else {
            throw VivoOmicsError.invalid("Visium filtered matrix barcodes must equal all in-tissue positions")
        }
        let observations: [VivoAssayObservation] = data.observations.map {
            let p = positions[$0.identity.barcode]!
            return .init(identity: $0.identity, kind: .spot, position: .init(frameID: frame.id, coordinates: [p.x, p.y]))
        }
        let result = VivoMultiAssayDataset(schemaVersion: data.schemaVersion, id: data.id, evidence: data.evidence,
            sourceDescription: data.sourceDescription, samples: data.samples, observations: observations,
            spatialFrames: [frame], assays: data.assays)
        try result.validate(); return result
    }
    /// Both inputs must be immutable snapshots. Counts use the existing native
    /// sparse 10x owner; coordinate rows join by barcode, never by source order.
    public static func readSnapshot(counts: URL, positions: URL, plan: VivoVisiumPlan) throws -> VivoMultiAssayDataset {
        try validate(plan)
        let bytes = try VivoSingleCellCampaignIO.readDocument(positions, maximumBytes: maximumPositionBytes)
        let table = try self.positions(bytes, format: plan.positionsFormat)
        return try attach(table, to: VivoMultiAssayTenX.readSnapshot(counts, plan: plan.counts), frame: plan.frame)
    }
}
