import Foundation

public struct VivoSingleCellImport: Codable, Sendable, Equatable {
    public let datasetID: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let countUnit: VivoOmicsCountUnit
    public let sample: VivoOmicsSample
    public let mitochondrialFeatureIDs: [String]
    public let cellGroups: [String: String]
    private enum CodingKeys: String, CodingKey { case datasetID, evidence, sourceDescription, countUnit, sample, mitochondrialFeatureIDs, cellGroups }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["datasetID", "evidence", "sourceDescription", "countUnit", "sample", "mitochondrialFeatureIDs", "cellGroups"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        datasetID = try values.decode(String.self, forKey: .datasetID)
        evidence = try values.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try values.decode(String.self, forKey: .sourceDescription)
        countUnit = try values.decode(VivoOmicsCountUnit.self, forKey: .countUnit)
        sample = try values.decode(VivoOmicsSample.self, forKey: .sample)
        mitochondrialFeatureIDs = try values.decodeIfPresent([String].self, forKey: .mitochondrialFeatureIDs) ?? []
        cellGroups = try values.decodeIfPresent([String: String].self, forKey: .cellGroups) ?? [:]
    }
    public init(datasetID: String, evidence: VivoOmicsEvidence, sourceDescription: String,
                countUnit: VivoOmicsCountUnit, sample: VivoOmicsSample,
                mitochondrialFeatureIDs: [String] = [], cellGroups: [String: String] = [:]) {
        self.datasetID = datasetID; self.evidence = evidence; self.sourceDescription = sourceDescription
        self.countUnit = countUnit; self.sample = sample
        self.mitochondrialFeatureIDs = mitochondrialFeatureIDs; self.cellGroups = cellGroups
    }
}

/// Strict uncompressed 10x MEX profile: integer/general coordinate Matrix Market,
/// three-column Gene Expression features, one barcode per line. Mixed modalities,
/// floating-point counts, duplicate coordinates and silently dropped rows are rejected.
public enum VivoMatrixMarketCounts {
    private static func lines(_ data: Data, maximumLines: Int, limits: VivoOmicsLimits) throws -> [Substring] {
        guard data.count <= limits.maximumInputBytes else { throw VivoOmicsError.limit("text bytes") }
        guard let text = String(data: data, encoding: .utf8), !text.contains("\0") else {
            throw VivoOmicsError.invalid("input must be uncompressed UTF-8 without NUL")
        }
        var rows = text.split(maxSplits: maximumLines, omittingEmptySubsequences: false, whereSeparator: { $0 == "\n" || $0 == "\r\n" })
        if rows.last?.isEmpty == true { rows.removeLast() }
        guard rows.count <= maximumLines else { throw VivoOmicsError.limit("text record count") }
        return try rows.map { row in
            let line = row.last == "\r" ? row.dropLast() : row
            guard line.utf8.count <= limits.maximumLineBytes else { throw VivoOmicsError.limit("text line bytes") }
            return line
        }
    }
    public static func decode(matrix: Data, features: Data, barcodes: Data,
                              metadata: VivoSingleCellImport, limits: VivoOmicsLimits = .init()) throws -> VivoSingleCellDataset {
        try limits.validate(); try metadata.sample.validate()
        var remaining = limits.maximumInputBytes
        for data in [matrix, features, barcodes] {
            guard data.count <= remaining else { throw VivoOmicsError.limit("aggregate source bytes") }
            remaining -= data.count
        }
        guard metadata.mitochondrialFeatureIDs.count <= limits.maximumFeatures,
              metadata.cellGroups.count <= limits.maximumCells else { throw VivoOmicsError.limit("annotation count") }
        let mitochondrial = Set(metadata.mitochondrialFeatureIDs)
        guard mitochondrial.count == metadata.mitochondrialFeatureIDs.count else { throw VivoOmicsError.invalid("duplicate mitochondrial annotation") }
        let featureLines = try lines(features, maximumLines: limits.maximumFeatures, limits: limits)
        guard !featureLines.isEmpty, featureLines.count <= limits.maximumFeatures else { throw VivoOmicsError.limit("features") }
        let decodedFeatures: [VivoOmicsFeature] = try featureLines.map { line in
            let fields = line.split(separator: "\t", omittingEmptySubsequences: false)
            guard fields.count == 3, fields[2] == "Gene Expression" else {
                throw VivoOmicsError.invalid("features.tsv requires ID, name and Gene Expression; no automatic modality filtering")
            }
            return .init(id: String(fields[0]), name: String(fields[1]), mitochondrial: mitochondrial.contains(String(fields[0])))
        }
        guard mitochondrial.isSubset(of: Set(decodedFeatures.map(\.id))) else { throw VivoOmicsError.invalid("unknown mitochondrial feature ID") }
        let barcodeLines = try lines(barcodes, maximumLines: limits.maximumCells, limits: limits)
        guard barcodeLines.count <= limits.maximumCells else { throw VivoOmicsError.limit("barcodes") }
        let cells = barcodeLines.map { VivoOmicsCell(barcode: String($0), sampleID: metadata.sample.id, group: metadata.cellGroups[String($0)]) }
        guard Set(metadata.cellGroups.keys).isSubset(of: Set(cells.map(\.barcode))) else { throw VivoOmicsError.invalid("annotation names an unknown barcode") }
        let (recordLimit, overflow) = limits.maximumNonzeros.addingReportingOverflow(4096)
        let matrixLines = try lines(matrix, maximumLines: overflow ? Int.max : recordLimit, limits: limits)
        guard matrixLines.first == "%%MatrixMarket matrix coordinate integer general" else {
            throw VivoOmicsError.invalid("requires Matrix Market coordinate integer general")
        }
        var dimensions: [Int]?, entries: [(row: Int, column: Int, count: UInt64)] = []
        for line in matrixLines.dropFirst() {
            if line.hasPrefix("%") { continue }
            let fields = line.split(whereSeparator: { $0 == " " || $0 == "\t" })
            guard fields.count == 3 else { throw VivoOmicsError.invalid("matrix record must contain exactly three fields") }
            if dimensions == nil {
                guard let rows = Int(fields[0]), let columns = Int(fields[1]), let nnz = Int(fields[2]),
                      rows == decodedFeatures.count, columns == cells.count, nnz >= 0,
                      nnz <= limits.maximumNonzeros else { throw VivoOmicsError.invalid("matrix dimensions/nonzeros disagree with metadata or limits") }
                dimensions = [rows, columns, nnz]
                // Do not allocate from a declared count until the number of source
                // records can support it. Sparse arrays grow with accepted input.
            } else {
                guard let shape = dimensions, entries.count < shape[2],
                      let feature = Int(fields[0]), let cell = Int(fields[1]),
                      feature >= 1, feature <= shape[0], cell >= 1, cell <= shape[1],
                      let count = UInt64(fields[2]), count > 0 else {
                    throw VivoOmicsError.invalid("invalid coordinate, count, or excess matrix records")
                }
                entries.append((cell - 1, feature - 1, count))
            }
        }
        guard let shape = dimensions, entries.count == shape[2] else { throw VivoOmicsError.invalid("truncated matrix or missing dimensions") }
        entries.sort { $0.row == $1.row ? $0.column < $1.column : $0.row < $1.row }
        var offsets = [Int](repeating: 0, count: cells.count + 1)
        var columns: [Int] = [], values: [UInt64] = []
        var previousRow = -1, previousColumn = -1
        for entry in entries {
            guard entry.row != previousRow || entry.column != previousColumn else { throw VivoOmicsError.invalid("duplicate matrix coordinate") }
            offsets[entry.row + 1] += 1; columns.append(entry.column); values.append(entry.count)
            previousRow = entry.row; previousColumn = entry.column
        }
        for index in 1..<offsets.count { offsets[index] += offsets[index - 1] }
        let result = VivoSingleCellDataset(id: metadata.datasetID, evidence: metadata.evidence,
            sourceDescription: metadata.sourceDescription, countUnit: metadata.countUnit, samples: [metadata.sample],
            features: decodedFeatures, cells: cells,
            matrix: .init(cellCount: cells.count, featureCount: decodedFeatures.count, rowOffsets: offsets, featureIndices: columns, counts: values))
        try result.validate(limits: limits); return result
    }
}
