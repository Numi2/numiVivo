import Foundation

/// Native MuData export. Each modality is its own AnnData with exact integer
/// counts; one-based maps preserve missing measurements and independent row order.
public enum VivoMultiAssayH5MU {
    public static func writeSnapshot(_ data: VivoMultiAssayDataset, to url: URL) throws {
        try data.validate()
        for a in data.assays { try vivoH5ADComponent(a.id) }
        try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.projectionFile(url.path)
            defer { h.close(file, "H5Fclose") }
            try h.encoding(file, "MuData", "0.1.0")
            func group(_ parent: Int64, _ name: String, encoding: String? = "dict") throws -> Int64 {
                let g = try h.projectionGroup(parent, name)
                var complete = false
                defer { if !complete { h.close(g, "H5Gclose") } }
                if let encoding { try h.encoding(g, encoding, encoding == "dict" ? "0.1.0" : "0.2.0") }
                complete = true; return g
            }
            func raw(_ parent: Int64, _ name: String, type: Int64, shape: [UInt64], buffer: UnsafeRawPointer?) throws {
                let d = try h.projectionDataset(parent, name, type: type, shape: shape); defer { h.close(d, "H5Dclose") }
                let write: @convention(c) (Int64, Int64, Int64, Int64, Int64, UnsafeRawPointer?) -> Int32 = try h.symbol("H5Dwrite")
                try h.check(write(d, type, 0, 0, 0, buffer), "multi-assay write")
                try h.encoding(d, "array", "0.2.0")
            }
            func integers(_ parent: Int64, _ name: String, _ values: [UInt64]) throws {
                try values.withUnsafeBytes { try raw(parent, name, type: h.native("NATIVE_ULLONG"), shape: [UInt64(values.count)], buffer: $0.baseAddress) }
            }
            func strings(_ parent: Int64, _ name: String, _ values: [String], scalar: Bool = false) throws {
                // A timestamp-free dataset with a variable-length UTF-8 type.
                let copy: @convention(c) (Int64) -> Int64 = try h.symbol("H5Tcopy")
                let type = try h.id(copy(h.native("C_S1")), "multi-assay string type"); defer { h.close(type, "H5Tclose") }
                let size: @convention(c) (Int64, Int) -> Int32 = try h.symbol("H5Tset_size")
                let charset: @convention(c) (Int64, Int32) -> Int32 = try h.symbol("H5Tset_cset")
                try h.check(size(type, -1), "variable string"); try h.check(charset(type, 1), "UTF8 string")
                let d = try h.projectionDataset(parent, name, type: type, shape: scalar ? [] : [UInt64(values.count)])
                defer { h.close(d, "H5Dclose") }
                let pointers = values.map { strdup($0) }; defer { pointers.forEach { free($0) } }
                guard pointers.allSatisfy({ $0 != nil }) else { throw VivoOmicsError.limit("multi-assay strings") }
                let write: @convention(c) (Int64, Int64, Int64, Int64, Int64, UnsafeRawPointer?) -> Int32 = try h.symbol("H5Dwrite")
                try pointers.withUnsafeBytes { try h.check(write(d, type, 0, 0, 0, $0.baseAddress), "write strings") }
                try h.encoding(d, scalar ? "string" : "string-array", "0.2.0")
            }
            func frame(_ parent: Int64, _ name: String, index: [String], columns: [(String, [String])]) throws {
                let g = try group(parent, name, encoding: "dataframe"); defer { h.close(g, "H5Gclose") }
                try h.writeStrings(g, "_index", ["_index"], attribute: true, scalar: true)
                try h.writeStrings(g, "column-order", columns.map(\.0), attribute: true)
                try strings(g, "_index", index)
                for (name, values) in columns { try strings(g, name, values) }
            }
            let enumCreate: @convention(c) (Int64) -> Int64 = try h.symbol("H5Tenum_create")
            let boolean = try h.id(enumCreate(h.native("NATIVE_UCHAR")), "boolean type"); defer { h.close(boolean, "H5Tclose") }
            let insert: @convention(c) (Int64, UnsafePointer<CChar>, UnsafeRawPointer) -> Int32 = try h.symbol("H5Tenum_insert")
            var zero: UInt8 = 0, one: UInt8 = 1
            try h.check(insert(boolean, "FALSE", &zero), "boolean false"); try h.check(insert(boolean, "TRUE", &one), "boolean true")
            func mask(_ parent: Int64, _ name: String, _ values: [UInt64]) throws {
                let bytes: [UInt8] = values.map { $0 == 0 ? 0 : 1 }
                try bytes.withUnsafeBytes { try raw(parent, name, type: boolean, shape: [UInt64(values.count), 1], buffer: $0.baseAddress) }
            }
            let obsNames = data.observations.indices.map { "observation-\($0)" }
            func obs(_ parent: Int64, rows: [Int]) throws {
                try frame(parent, "obs", index: rows.map { obsNames[$0] }, columns: [
                    ("barcode", rows.map { data.observations[$0].identity.barcode }),
                    ("sample", rows.map { data.observations[$0].identity.sampleID }),
                    ("observation_kind", rows.map { data.observations[$0].kind.rawValue })])
            }
            try obs(file, rows: Array(data.observations.indices))
            let allFeatures = data.assays.flatMap(\.features)
            try frame(file, "var", index: allFeatures.map(\.id), columns: [("feature_space", data.assays.flatMap { a in a.features.map { _ in a.id } })])
            let mod = try group(file, "mod", encoding: nil); defer { h.close(mod, "H5Gclose") }
            try h.writeStrings(mod, "mod-order", data.assays.map(\.id), attribute: true)
            let obsmap = try group(file, "obsmap"), varmap = try group(file, "varmap")
            let obsm = try group(file, "obsm"), varm = try group(file, "varm")
            defer { for g in [obsmap, varmap, obsm, varm] { h.close(g, "H5Gclose") } }
            for name in ["obsp", "varp"] { h.close(try group(file, name), "H5Gclose") }
            var featureOffset = 0
            for a in data.assays {
                try Task.checkCancellation()
                let g = try group(mod, a.id, encoding: nil); defer { h.close(g, "H5Gclose") }
                try h.encoding(g, "anndata", "0.1.0")
                try obs(g, rows: a.observationIndices)
                try frame(g, "var", index: a.features.map(\.id), columns: [("name", a.features.map(\.name)),
                    ("feature_namespace", a.features.map { _ in a.featureNamespace }), ("assay_kind", a.features.map { _ in a.kind.rawValue })])
                let x = try group(g, "X", encoding: nil); defer { h.close(x, "H5Gclose") }
                try h.encoding(x, "csr_matrix", "0.1.0")
                try h.writeIntegers(x, "shape", [UInt64(a.matrix.cellCount), UInt64(a.matrix.featureCount)], attribute: true)
                // Sparse indices are signed as expected by scipy; counts stay UInt64.
                for (name, values) in [("indptr", a.matrix.rowOffsets), ("indices", a.matrix.featureIndices)] {
                    let signed = values.map(Int64.init)
                    try signed.withUnsafeBytes { try raw(x, name, type: h.native("NATIVE_LLONG"), shape: [UInt64(signed.count)], buffer: $0.baseAddress) }
                }
                try integers(x, "data", a.matrix.counts)
                for name in ["layers", "obsm", "obsp", "varm", "varp", "uns"] { h.close(try group(g, name), "H5Gclose") }
                var om = [UInt64](repeating: 0, count: data.observations.count)
                for (local, global) in a.observationIndices.enumerated() { om[global] = UInt64(local + 1) }
                var vm = [UInt64](repeating: 0, count: allFeatures.count)
                for i in a.features.indices { vm[featureOffset + i] = UInt64(i + 1) }
                featureOffset += a.features.count
                try integers(obsmap, a.id, om); try integers(varmap, a.id, vm)
                try mask(obsm, a.id, om); try mask(varm, a.id, vm)
            }
            struct AssayMetadata: Encodable {
                let id: String; let kind: VivoAssayKind; let featureNamespace: String; let countUnit: VivoAssayCountUnit
                let genomeAssembly: String?; let sourceDescription: String; let features: [VivoAssayFeature]; let observationIndices: [Int]
            }
            struct Metadata: Encodable {
                let schema = "numivivo.org/multi-assay/v1"
                let id: String; let evidence: VivoOmicsEvidence; let sourceDescription: String
                let samples: [VivoOmicsSample]; let observations: [VivoAssayObservation]; let spatialFrames: [VivoSpatialFrame]; let assays: [AssayMetadata]
            }
            let metadata = Metadata(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
                samples: data.samples, observations: data.observations, spatialFrames: data.spatialFrames,
                assays: data.assays.map { .init(id: $0.id, kind: $0.kind, featureNamespace: $0.featureNamespace, countUnit: $0.countUnit,
                    genomeAssembly: $0.genomeAssembly, sourceDescription: $0.sourceDescription, features: $0.features, observationIndices: $0.observationIndices) })
            let uns = try group(file, "uns"); defer { h.close(uns, "H5Gclose") }
            try strings(uns, "numivivo", [String(decoding: try VivoCanonicalJSON.encode(metadata), as: UTF8.self)], scalar: true)
            try h.projectionStorageLimit(file)
        }
    }
}

extension VivoMultiAssayH5MU {
    public static func write(_ data: VivoMultiAssayDataset, to destination: URL) throws {
        try data.validate()
        guard destination.isFileURL else { throw VivoOmicsError.invalid("H5MU requires local output") }
        let temp = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-h5mu-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temp, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temp) }
        let file = temp.appendingPathComponent("data.h5mu")
        try writeSnapshot(data, to: file)
        try VivoSingleCellH5AD.publish(VivoSingleCellCampaignIO.readDocument(file, maximumBytes: 1_073_741_824), to: destination)
    }
}
