import Foundation

/// An explicit assay/design mapping: AnnData does not standardize count layers,
/// donor columns, organisms, or biological replication. Never infer these.
public struct VivoH5ADImportPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let countUnit: VivoOmicsCountUnit
    public let matrixPath: String
    public let samples: [VivoOmicsSample]
    public let sampleColumn: String
    public let barcodeColumn: String?
    public let groupColumn: String?
    public let featureIDColumn: String?
    public let featureNameColumn: String?
    public let mitochondrialFeatureIDs: [String]
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, evidence, sourceDescription, countUnit, matrixPath, samples, sampleColumn, barcodeColumn, groupColumn, featureIDColumn, featureNameColumn, mitochondrialFeatureIDs }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "evidence", "sourceDescription", "countUnit", "matrixPath", "samples", "sampleColumn", "barcodeColumn", "groupColumn", "featureIDColumn", "featureNameColumn", "mitochondrialFeatureIDs"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        id = try values.decode(String.self, forKey: .id)
        evidence = try values.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try values.decode(String.self, forKey: .sourceDescription)
        countUnit = try values.decode(VivoOmicsCountUnit.self, forKey: .countUnit)
        matrixPath = try values.decode(String.self, forKey: .matrixPath)
        samples = try values.decode([VivoOmicsSample].self, forKey: .samples)
        sampleColumn = try values.decode(String.self, forKey: .sampleColumn)
        barcodeColumn = try values.decodeIfPresent(String.self, forKey: .barcodeColumn)
        groupColumn = try values.decodeIfPresent(String.self, forKey: .groupColumn)
        featureIDColumn = try values.decodeIfPresent(String.self, forKey: .featureIDColumn)
        featureNameColumn = try values.decodeIfPresent(String.self, forKey: .featureNameColumn)
        mitochondrialFeatureIDs = try values.decodeIfPresent([String].self, forKey: .mitochondrialFeatureIDs) ?? []
    }
    public init(id: String, evidence: VivoOmicsEvidence, sourceDescription: String, countUnit: VivoOmicsCountUnit,
                matrixPath: String, samples: [VivoOmicsSample], sampleColumn: String, barcodeColumn: String? = nil,
                groupColumn: String? = nil, featureIDColumn: String? = nil, featureNameColumn: String? = nil, mitochondrialFeatureIDs: [String] = []) {
        self.schemaVersion = 1; self.id = id; self.evidence = evidence; self.sourceDescription = sourceDescription; self.countUnit = countUnit
        self.matrixPath = matrixPath; self.samples = samples; self.sampleColumn = sampleColumn; self.barcodeColumn = barcodeColumn
        self.groupColumn = groupColumn; self.featureIDColumn = featureIDColumn; self.featureNameColumn = featureNameColumn; self.mitochondrialFeatureIDs = mitochondrialFeatureIDs
    }
}

/// The native count projection is deliberately distinct from the AnnData source.
/// Keeping the entire original file retains categories, nullable columns, raw,
/// layers, embeddings, graphs and extension encodings without coercion or loss.
/// Re-export of this source never pretends to contain downstream count edits.
public struct VivoH5ADDocument: Sendable {
    public let source: Data
    public let hdf5Version: String
    public let dataset: VivoSingleCellDataset
    public let plan: VivoH5ADImportPlan
    public func exportOriginal(to destination: URL) throws {
        try VivoSingleCellH5AD.publish(source, to: destination)
    }
}

public enum VivoSingleCellH5AD {
    static func publish(_ bytes: Data, to destination: URL) throws {
        guard destination.isFileURL else { throw VivoOmicsError.invalid("local H5AD output required") }
        let files = try VivoRootedFileStore(rootURL: destination.deletingLastPathComponent(), createIfNeeded: false)
        guard try files.writeFile(bytes, relative: destination.lastPathComponent, immutable: true) else {
            throw VivoOmicsError.invalid("H5AD output already exists")
        }
    }
    public static func read(_ url: URL, plan: VivoH5ADImportPlan, limits: VivoOmicsLimits = .init()) throws -> VivoH5ADDocument {
        try limits.validate()
        guard plan.schemaVersion == 1 else { throw VivoOmicsError.invalid("unsupported H5AD mapping schema") }
        let bytes = try VivoSingleCellCampaignIO.readDocument(url, maximumBytes: limits.maximumInputBytes)
        // Snapshot before HDF5 reads so the retained source and projection always
        // refer to the same file, including during concurrent external changes.
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false,
                                                attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: directory) }
        let snapshot = directory.appendingPathComponent("source.h5ad")
        try bytes.write(to: snapshot, options: .withoutOverwriting)
        var metadata: VivoSingleCellCountMetadata?
        var rows: [[Int: UInt64]] = []
        let version = try scanSnapshot(snapshot,plan: plan,limits: limits,onMetadata: {
            metadata = $0; rows = Array(repeating: [:],count: $0.cells.count)
        },onEntry: { row,feature,count in rows[row][feature] = count })
        guard let metadata else { throw VivoOmicsError.invalid("missing count metadata") }
        var offsets = [0], columns: [Int] = [], counts: [UInt64] = []
        for row in rows {
            for feature in row.keys.sorted() { columns.append(feature); counts.append(row[feature]!) }
            offsets.append(counts.count)
        }
        let result = VivoSingleCellDataset(id: metadata.id,evidence: metadata.evidence,sourceDescription: metadata.sourceDescription,
            countUnit: metadata.countUnit,samples: metadata.samples,features: metadata.features,cells: metadata.cells,
            matrix: .init(cellCount: metadata.cells.count,featureCount: metadata.features.count,rowOffsets: offsets,featureIndices: columns,counts: counts))
        try result.validate(limits: limits)
        return .init(source: bytes,hdf5Version: version,dataset: result,plan: plan)
    }
    /// The caller owns an immutable source snapshot. Sparse arrays are read in
    /// bounded slices and duplicate coordinates merged per major segment.
    static func scanSnapshot(_ url: URL,plan: VivoH5ADImportPlan,limits: VivoOmicsLimits,
        onMetadata: (VivoSingleCellCountMetadata) throws -> Void,
        onEntry: (Int,Int,UInt64) throws -> Void) throws -> String {
        try limits.validate()
        guard plan.schemaVersion == 1 else { throw VivoOmicsError.invalid("unsupported H5AD mapping schema") }
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(url.path)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "anndata", try h.text(file, "encoding-version") == "0.1.0" else {
                throw VivoOmicsError.invalid("unsupported AnnData root encoding")
            }
            let reader = VivoH5ADFrameReader(h: h, file: file)
            func component(_ name: String) throws -> String { try reader.component(name) }
            func column(_ frame: String, _ name: String, maximum: Int) throws -> [String?] { try reader.column(frame, name, maximum: maximum) }
            func required(_ values: [String?]) throws -> [String] { try reader.required(values) }
            func index(_ frame: String, maximum: Int) throws -> [String] { try reader.index(frame, maximum: maximum) }
            let featureFrame = plan.matrixPath == "raw/X" ? "raw/var" : "var"
            if plan.matrixPath == "raw/X" {
                let raw = try h.object(file, "raw"); defer { h.close(raw, "H5Oclose") }
                guard try h.text(raw, "encoding-type") == "raw", try h.text(raw, "encoding-version") == "0.1.0" else {
                    throw VivoOmicsError.invalid("unsupported AnnData raw encoding")
                }
            }
            let observationCount = try reader.indexLength("obs", maximum: limits.maximumCells)
            let featureCount = try reader.indexLength(featureFrame, maximum: limits.maximumFeatures)
            let features = try plan.featureIDColumn.map { try required(column(featureFrame, $0, maximum: limits.maximumFeatures)) } ?? index(featureFrame, maximum: limits.maximumFeatures)
            let sampleIDs = try required(column("obs", plan.sampleColumn, maximum: limits.maximumCells))
            let barcodes = try plan.barcodeColumn.map { try required(column("obs", $0, maximum: limits.maximumCells)) } ?? index("obs", maximum: limits.maximumCells)
            let groups = try plan.groupColumn.map { try column("obs", $0, maximum: limits.maximumCells) }
            let names = try plan.featureNameColumn.map { try required(column(featureFrame, $0, maximum: limits.maximumFeatures)) } ?? index(featureFrame, maximum: limits.maximumFeatures)
            guard sampleIDs.count == observationCount, barcodes.count == observationCount, groups == nil || groups!.count == observationCount,
                  names.count == features.count, features.count == featureCount else { throw VivoOmicsError.invalid("AnnData annotation dimensions disagree") }
            guard plan.matrixPath == "X" || plan.matrixPath == "raw/X" || (plan.matrixPath.hasPrefix("layers/") && plan.matrixPath.split(separator: "/", omittingEmptySubsequences: false).count == 2) else {
                throw VivoOmicsError.invalid("select X, raw/X or layers/<name> explicitly")
            }
            if plan.matrixPath.hasPrefix("layers/") { _ = try component(String(plan.matrixPath.dropFirst(7))) }
            let mitochondrial = Set(plan.mitochondrialFeatureIDs)
            guard mitochondrial.isSubset(of: Set(features)) else { throw VivoOmicsError.invalid("unknown mitochondrial feature ID") }
            let metadata = VivoSingleCellCountMetadata(id: plan.id,evidence: plan.evidence,sourceDescription: plan.sourceDescription,
                countUnit: plan.countUnit,samples: plan.samples,
                features: features.enumerated().map { .init(id: $0.element,name: names[$0.offset],mitochondrial: mitochondrial.contains($0.element)) },
                cells: (0..<observationCount).map { .init(barcode: barcodes[$0],sampleID: sampleIDs[$0],group: groups?[$0]) })
            try metadata.validate(limits: limits)
            try onMetadata(metadata)
            try VivoH5ADCountReader.scan(h, file: file, path: plan.matrixPath, rows: observationCount, features: features.count, limits: limits, onEntry: onEntry)
            return try h.version()
        }
    }
}

extension VivoSingleCellH5AD {
    /// Creates a new AnnData count object from the native dataset. Its metadata
    /// describes this dataset, not any separately retained source AnnData object.
    public static func write(_ dataset: VivoSingleCellDataset, to destination: URL) throws {
        try dataset.validate()
        guard destination.isFileURL else { throw VivoOmicsError.invalid("local H5AD destination required") }
        let directory = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5ad-write-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: directory) }
        let staging = directory.appendingPathComponent("counts.h5ad")
        try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(staging.path, create: true)
            defer { h.close(file, "H5Fclose") }
            try h.encoding(file, "anndata", "0.1.0")
            func strings(_ group: Int64, _ name: String, _ values: [String]) throws {
                try h.writeStrings(group, name, values)
                let d = try h.dataset(group, name); defer { h.close(d, "H5Dclose") }
                try h.encoding(d, "string-array", "0.2.0")
            }
            func frame(_ name: String, index: [String], columns: [(String, [String])]) throws {
                let g = try h.group(file, name); defer { h.close(g, "H5Gclose") }
                try h.encoding(g, "dataframe", "0.2.0")
                try h.writeStrings(g, "_index", ["_index"], attribute: true, scalar: true)
                try h.writeStrings(g, "column-order", columns.map(\.0), attribute: true)
                try strings(g, "_index", index)
                for (name, values) in columns { try strings(g, name, values) }
            }
            let sampleByID = Dictionary(uniqueKeysWithValues: dataset.samples.map { ($0.id, $0) })
            let cells = dataset.cells
            // Unique generated observation index; original barcode and sample
            // identity stay separate and exact even across reused barcodes.
            var obs: [(String, [String])] = [
                ("barcode", cells.map(\.barcode)), ("sample", cells.map(\.sampleID)),
                ("biological_replicate", cells.map { sampleByID[$0.sampleID]!.biologicalReplicateID }),
                ("condition", cells.map { sampleByID[$0.sampleID]!.condition }),
                ("batch", cells.map { sampleByID[$0.sampleID]!.batchID }),
                ("organism", cells.map { sampleByID[$0.sampleID]!.organism })]
            // Optional values use nullable string-array encoding with a mask.
            try frame("obs", index: cells.indices.map { "cell-\($0)" }, columns: obs)
            func nullable(_ name: String, _ values: [String?]) throws {
                let g = try h.group(file, "obs/" + name); defer { h.close(g, "H5Gclose") }
                try h.encoding(g, "nullable-string-array", "0.1.0")
                try strings(g, "values", values.map { $0 ?? "" })
                // AnnData expects a boolean mask, represented by an HDF5 enum.
                let enumCreate: @convention(c) (Int64) -> Int64 = try h.symbol("H5Tenum_create")
                let t = try h.id(enumCreate(h.native("NATIVE_UCHAR")), "create boolean type"); defer { h.close(t, "H5Tclose") }
                let insert: @convention(c) (Int64, UnsafePointer<CChar>, UnsafeRawPointer) -> Int32 = try h.symbol("H5Tenum_insert")
                var zero: UInt8 = 0, one: UInt8 = 1
                try h.check(insert(t, "FALSE", &zero), "false enum"); try h.check(insert(t, "TRUE", &one), "true enum")
                let mask: [UInt8] = values.map { $0 == nil ? 1 : 0 }
                try mask.withUnsafeBytes { try h.write(g, "mask", type: t, dimensions: [UInt64(values.count)], attribute: false, buffer: $0.baseAddress) }
            }
            try nullable("donor", cells.map { sampleByID[$0.sampleID]!.donorID })
            try nullable("group", cells.map(\.group))
            // Replace column-order to include nullable columns.
            obs += [("donor", []), ("group", [])]
            let og = try h.object(file, "obs"); defer { h.close(og, "H5Oclose") }
            let delete: @convention(c) (Int64, UnsafePointer<CChar>) -> Int32 = try h.symbol("H5Adelete")
            try h.check(delete(og, "column-order"), "replace column order")
            try h.writeStrings(og, "column-order", obs.map(\.0), attribute: true)
            try frame("var", index: dataset.features.map(\.id), columns: [
                ("name", dataset.features.map(\.name)), ("mitochondrial", dataset.features.map { String($0.mitochondrial) })])
            let x = try h.group(file, "X"); defer { h.close(x, "H5Gclose") }
            try h.encoding(x, "csr_matrix", "0.1.0")
            try h.writeIntegers(x, "shape", [UInt64(cells.count), UInt64(dataset.features.count)], attribute: true)
            try h.writeIndices(x, "indptr", dataset.matrix.rowOffsets)
            try h.writeIndices(x, "indices", dataset.matrix.featureIndices)
            try h.writeIntegers(x, "data", dataset.matrix.counts)
            for name in ["layers", "obsm", "obsp", "varm", "varp", "uns"] {
                let g = try h.group(file, name); defer { h.close(g, "H5Gclose") }
                try h.encoding(g, "dict", "0.1.0")
                if name == "uns" {
                    // Complete metadata remains recoverable, including samples
                    // with zero cells and missing-vs-empty optional annotations.
                    struct Metadata: Encodable {
                        let schema = "numivivo.org/h5ad-count-metadata/v1"
                        let id: String; let evidence: VivoOmicsEvidence; let sourceDescription: String
                        let countUnit: VivoOmicsCountUnit; let samples: [VivoOmicsSample]
                    }
                    let metadata = Metadata(id: dataset.id, evidence: dataset.evidence, sourceDescription: dataset.sourceDescription,
                                            countUnit: dataset.countUnit, samples: dataset.samples)
                    let json = try String(decoding: VivoCanonicalJSON.encode(metadata), as: UTF8.self)
                    try h.writeStrings(g, "numivivo", [json], scalar: true)
                    let d = try h.dataset(g, "numivivo"); defer { h.close(d, "H5Dclose") }
                    try h.encoding(d, "string", "0.2.0")
                }
            }
            let flush: @convention(c) (Int64, Int32) -> Int32 = try h.symbol("H5Fflush")
            try h.check(flush(file, 1), "flush output")
        }
        try Task.checkCancellation()
        let bytes = try VivoSingleCellCampaignIO.readDocument(staging, maximumBytes: 128 * 1_024 * 1_024)
        try publish(bytes, to: destination)
    }
}
