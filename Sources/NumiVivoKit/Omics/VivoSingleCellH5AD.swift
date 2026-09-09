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
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(snapshot.path)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "anndata", try h.text(file, "encoding-version") == "0.1.0" else {
                throw VivoOmicsError.invalid("unsupported AnnData root encoding")
            }
            func component(_ value: String) throws -> String {
                guard !value.isEmpty, !value.contains("/"), value != ".", value != "..", !value.contains("\0") else {
                    throw VivoOmicsError.invalid("invalid AnnData component")
                }; return value
            }
            func column(_ frame: String, _ name: String, maximum: Int) throws -> [String?] {
                let path = frame + "/" + (try component(name)), object = try h.object(file, path)
                defer { h.close(object, "H5Oclose") }
                let encoding = try h.text(object, "encoding-type"), version = try h.text(object, "encoding-version")
                if encoding == "categorical", version == "0.2.0" {
                    let categories = try h.dataset(file, path + "/categories"), codes = try h.dataset(file, path + "/codes")
                    defer { h.close(categories, "H5Dclose"); h.close(codes, "H5Dclose") }
                    guard try h.shape(categories, attribute: false).count == 1, try h.shape(codes, attribute: false).count == 1 else {
                        throw VivoOmicsError.invalid("categorical arrays must be vectors")
                    }
                    let labels = try h.strings(categories, maximum: maximum)
                    // Missing categorical values stay missing; required design identities reject them below.
                    return try h.categoricalCodes(codes, maximum: maximum).map {
                        if $0 == -1 { return nil }
                        guard $0 >= 0, $0 < labels.count else { throw VivoOmicsError.invalid("categorical code out of range") }
                        return labels[Int($0)]
                    }
                }
                if encoding == "nullable-string-array", version == "0.1.0" {
                    let values = try h.dataset(file, path + "/values"), mask = try h.dataset(file, path + "/mask")
                    defer { h.close(values, "H5Dclose"); h.close(mask, "H5Dclose") }
                    guard try h.shape(values, attribute: false).count == 1, try h.shape(mask, attribute: false).count == 1 else {
                        throw VivoOmicsError.invalid("nullable arrays must be vectors")
                    }
                    let strings = try h.strings(values, maximum: maximum), missing = try h.mask(mask, maximum: maximum)
                    guard strings.count == missing.count else { throw VivoOmicsError.invalid("nullable mask shape differs") }
                    return strings.indices.map { missing[$0] ? nil : strings[$0] }
                }
                guard encoding == "string-array", version == "0.2.0" else { throw VivoOmicsError.invalid("selected identity column must be string-array or categorical") }
                let d = try h.dataset(file, path); defer { h.close(d, "H5Dclose") }
                guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("identity column must be a vector") }
                return try h.strings(d, maximum: maximum)
            }
            func required(_ values: [String?]) throws -> [String] {
                try values.map { guard let value = $0 else { throw VivoOmicsError.invalid("missing selected design identity") }; return value }
            }
            func index(_ frame: String, maximum: Int) throws -> [String] {
                let group = try h.object(file, frame); defer { h.close(group, "H5Oclose") }
                guard try h.text(group, "encoding-type") == "dataframe", try h.text(group, "encoding-version") == "0.2.0" else {
                    throw VivoOmicsError.invalid("unsupported AnnData dataframe")
                }
                return try required(column(frame, h.text(group, "_index"), maximum: maximum))
            }
            let featureFrame = plan.matrixPath == "raw/X" ? "raw/var" : "var"
            if plan.matrixPath == "raw/X" {
                let raw = try h.object(file, "raw"); defer { h.close(raw, "H5Oclose") }
                guard try h.text(raw, "encoding-type") == "raw", try h.text(raw, "encoding-version") == "0.1.0" else {
                    throw VivoOmicsError.invalid("unsupported AnnData raw encoding")
                }
            }
            let obs = try index("obs", maximum: limits.maximumCells), featureIndex = try index(featureFrame, maximum: limits.maximumFeatures)
            let features = try plan.featureIDColumn.map { try required(column(featureFrame, $0, maximum: limits.maximumFeatures)) } ?? featureIndex
            let sampleIDs = try required(column("obs", plan.sampleColumn, maximum: limits.maximumCells))
            let barcodes = try plan.barcodeColumn.map { try required(column("obs", $0, maximum: limits.maximumCells)) } ?? obs
            let groups = try plan.groupColumn.map { try column("obs", $0, maximum: limits.maximumCells) }
            let names = try plan.featureNameColumn.map { try required(column(featureFrame, $0, maximum: limits.maximumFeatures)) } ?? featureIndex
            guard sampleIDs.count == obs.count, barcodes.count == obs.count, groups == nil || groups!.count == obs.count,
                  names.count == features.count, features.count == featureIndex.count else { throw VivoOmicsError.invalid("AnnData annotation dimensions disagree") }
            guard plan.matrixPath == "X" || plan.matrixPath == "raw/X" || (plan.matrixPath.hasPrefix("layers/") && plan.matrixPath.split(separator: "/", omittingEmptySubsequences: false).count == 2) else {
                throw VivoOmicsError.invalid("select X, raw/X or layers/<name> explicitly")
            }
            if plan.matrixPath.hasPrefix("layers/") { _ = try component(String(plan.matrixPath.dropFirst(7))) }
            let matrix = try h.object(file, plan.matrixPath); defer { h.close(matrix, "H5Oclose") }
            let encoding = try h.text(matrix, "encoding-type")
            var rows = [[Int: UInt64]](repeating: [:], count: obs.count)
            if encoding == "array" {
                guard try h.text(matrix, "encoding-version") == "0.2.0" else { throw VivoOmicsError.invalid("unsupported dense encoding") }
                let d = try h.dataset(file, plan.matrixPath); defer { h.close(d, "H5Dclose") }
                guard try h.shape(d, attribute: false) == [UInt64(obs.count), UInt64(features.count)] else {
                    throw VivoOmicsError.invalid("dense matrix shape differs from obs/var")
                }
                var retained = 0
                for row in obs.indices {
                    try Task.checkCancellation()
                    let values = try h.integers(d, maximum: limits.maximumFeatures, allowFloat: true, row: row)
                    for (feature, count) in values.enumerated() where count > 0 {
                        guard retained < limits.maximumNonzeros else { throw VivoOmicsError.limit("dense input exceeds sparse nonzero allowance") }
                        rows[row][feature] = count; retained += 1
                    }
                }
            } else {
                guard ["csr_matrix", "csc_matrix"].contains(encoding), try h.text(matrix, "encoding-version") == "0.1.0" else {
                    throw VivoOmicsError.invalid("count import requires CSR, CSC or dense numeric array")
                }
                let shapeAttribute = try h.attribute(matrix, "shape"); defer { h.close(shapeAttribute, "H5Aclose") }
                guard try h.shape(shapeAttribute, attribute: true) == [2] else { throw VivoOmicsError.invalid("sparse shape must have two dimensions") }
                let shape = try h.integers(shapeAttribute, attribute: true, maximum: 2)
                guard shape == [UInt64(obs.count), UInt64(features.count)] else { throw VivoOmicsError.invalid("matrix shape differs from obs/var") }
                func integers(_ name: String, maximum: Int, counts: Bool = false) throws -> [UInt64] {
                    let d = try h.dataset(file, plan.matrixPath + "/" + name); defer { h.close(d, "H5Dclose") }
                    guard try h.shape(d, attribute: false).count == 1 else { throw VivoOmicsError.invalid("sparse arrays must be vectors") }
                    return try h.integers(d, maximum: maximum, allowFloat: counts)
                }
                let major = encoding == "csr_matrix" ? obs.count : features.count
                let minor = encoding == "csr_matrix" ? features.count : obs.count
                let offsets = try integers("indptr", maximum: major + 1)
                let indices = try integers("indices", maximum: limits.maximumNonzeros)
                let counts = try integers("data", maximum: limits.maximumNonzeros, counts: true)
                guard offsets.count == major + 1, offsets.first == 0, offsets.last == UInt64(counts.count), indices.count == counts.count else {
                    throw VivoOmicsError.invalid("malformed sparse array lengths")
                }
                // Canonicalize unsorted indices, duplicate entries and explicit zeros
                // using O(cells + nnz) storage; never allocate cells x genes.
                for i in 0..<major {
                    try Task.checkCancellation()
                    guard offsets[i] <= offsets[i + 1], offsets[i + 1] <= counts.count else { throw VivoOmicsError.invalid("malformed sparse offsets") }
                    for k in Int(offsets[i])..<Int(offsets[i + 1]) {
                        guard indices[k] < minor else { throw VivoOmicsError.invalid("sparse index out of range") }
                        let row = encoding == "csr_matrix" ? i : Int(indices[k]), feature = encoding == "csr_matrix" ? Int(indices[k]) : i
                        if counts[k] > 0 { rows[row][feature] = try vivoOmicsSum(rows[row][feature] ?? 0, counts[k]) }
                    }
                }
            }
            var rowOffsets = [0], featureIndices: [Int] = [], values: [UInt64] = []
            for row in rows {
                for feature in row.keys.sorted() { featureIndices.append(feature); values.append(row[feature]!) }
                rowOffsets.append(values.count)
            }
            let mitochondrial = Set(plan.mitochondrialFeatureIDs)
            guard mitochondrial.isSubset(of: Set(features)) else { throw VivoOmicsError.invalid("unknown mitochondrial feature ID") }
            let result = VivoSingleCellDataset(id: plan.id, evidence: plan.evidence, sourceDescription: plan.sourceDescription,
                countUnit: plan.countUnit, samples: plan.samples,
                features: features.enumerated().map { .init(id: $0.element, name: names[$0.offset], mitochondrial: mitochondrial.contains($0.element)) },
                cells: obs.indices.map { .init(barcode: barcodes[$0], sampleID: sampleIDs[$0], group: groups?[$0]) },
                matrix: .init(cellCount: obs.count, featureCount: features.count, rowOffsets: rowOffsets, featureIndices: featureIndices, counts: values))
            try result.validate(limits: limits)
            return try .init(source: bytes, hdf5Version: h.version(), dataset: result, plan: plan)
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
