import Foundation

/// Explicit mapping for one continuous-valued MuData modality.  The mapping
/// keeps the matrix and feature identity choices separate from count-assay
/// plans, so a fractional intensity can never be decoded as a molecule count.
public struct VivoH5MUQuantitativeAssayMapping: Codable, Sendable, Equatable {
    public let sourceName: String
    public let id: String
    public let kind: VivoQuantitativeAssayKind
    public let featureNamespace: String
    public let unit: String
    public let sourceDescription: String
    public let matrixPath: String
    public let featureIDColumn: String?
    public let featureNameColumn: String?
    public let intervalContigColumn: String?
    public let intervalStartColumn: String?
    public let intervalEndColumn: String?
    public let intervalPresentColumn: String?

    public init(sourceName: String, id: String, kind: VivoQuantitativeAssayKind,
                featureNamespace: String, unit: String, sourceDescription: String,
                matrixPath: String = "X", featureIDColumn: String? = nil,
                featureNameColumn: String? = nil, intervalContigColumn: String? = nil,
                intervalStartColumn: String? = nil, intervalEndColumn: String? = nil,
                intervalPresentColumn: String? = nil) {
        self.sourceName = sourceName; self.id = id; self.kind = kind
        self.featureNamespace = featureNamespace; self.unit = unit
        self.sourceDescription = sourceDescription; self.matrixPath = matrixPath
        self.featureIDColumn = featureIDColumn; self.featureNameColumn = featureNameColumn
        self.intervalContigColumn = intervalContigColumn; self.intervalStartColumn = intervalStartColumn
        self.intervalEndColumn = intervalEndColumn; self.intervalPresentColumn = intervalPresentColumn
    }

    private enum CodingKeys: String, CodingKey {
        case sourceName, id, kind, featureNamespace, unit, sourceDescription, matrixPath
        case featureIDColumn, featureNameColumn, intervalContigColumn, intervalStartColumn
        case intervalEndColumn, intervalPresentColumn
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "sourceName", "id", "kind", "featureNamespace", "unit", "sourceDescription", "matrixPath",
            "featureIDColumn", "featureNameColumn", "intervalContigColumn", "intervalStartColumn",
            "intervalEndColumn", "intervalPresentColumn"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        sourceName = try v.decode(String.self, forKey: .sourceName)
        id = try v.decode(String.self, forKey: .id)
        kind = try v.decode(VivoQuantitativeAssayKind.self, forKey: .kind)
        featureNamespace = try v.decode(String.self, forKey: .featureNamespace)
        unit = try v.decode(String.self, forKey: .unit)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        matrixPath = try v.decodeIfPresent(String.self, forKey: .matrixPath) ?? "X"
        featureIDColumn = try v.decodeIfPresent(String.self, forKey: .featureIDColumn)
        featureNameColumn = try v.decodeIfPresent(String.self, forKey: .featureNameColumn)
        intervalContigColumn = try v.decodeIfPresent(String.self, forKey: .intervalContigColumn)
        intervalStartColumn = try v.decodeIfPresent(String.self, forKey: .intervalStartColumn)
        intervalEndColumn = try v.decodeIfPresent(String.self, forKey: .intervalEndColumn)
        intervalPresentColumn = try v.decodeIfPresent(String.self, forKey: .intervalPresentColumn)
    }
}

public struct VivoH5MUQuantitativeSpatialMapping: Codable, Sendable, Equatable {
    public let path: String
    public let frame: VivoSpatialFrame
    public init(path: String, frame: VivoSpatialFrame) {
        self.path = path; self.frame = frame
    }
    private enum CodingKeys: String, CodingKey { case path, frame }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["path", "frame"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        path = try v.decode(String.self, forKey: .path)
        frame = try v.decode(VivoSpatialFrame.self, forKey: .frame)
    }
}

/// H5MU plans for continuous assays.  Every design column and modality is
/// selected explicitly; no `uns` metadata is treated as an authority for
/// sample, feature or unit identity.
public struct VivoH5MUQuantitativePlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let id: String
    public let evidence: VivoOmicsEvidence
    public let sourceDescription: String
    public let samples: [VivoOmicsSample]
    public let sampleColumn: String
    public let barcodeColumn: String?
    public let groupColumn: String?
    public let defaultObservationKind: VivoObservationKind
    public let observationKindColumn: String?
    public let spatial: [VivoH5MUQuantitativeSpatialMapping]
    public let assays: [VivoH5MUQuantitativeAssayMapping]

    public init(schemaVersion: Int = 1, id: String, evidence: VivoOmicsEvidence,
                sourceDescription: String, samples: [VivoOmicsSample], sampleColumn: String,
                barcodeColumn: String? = nil, groupColumn: String? = nil,
                defaultObservationKind: VivoObservationKind = .cell,
                observationKindColumn: String? = nil,
                spatial: [VivoH5MUQuantitativeSpatialMapping] = [],
                assays: [VivoH5MUQuantitativeAssayMapping]) {
        self.schemaVersion = schemaVersion; self.id = id; self.evidence = evidence
        self.sourceDescription = sourceDescription; self.samples = samples
        self.sampleColumn = sampleColumn; self.barcodeColumn = barcodeColumn
        self.groupColumn = groupColumn; self.defaultObservationKind = defaultObservationKind
        self.observationKindColumn = observationKindColumn; self.spatial = spatial; self.assays = assays
    }

    private enum CodingKeys: String, CodingKey {
        case schemaVersion, id, evidence, sourceDescription, samples, sampleColumn, barcodeColumn
        case groupColumn, defaultObservationKind, observationKindColumn, spatial, assays
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: [
            "schemaVersion", "id", "evidence", "sourceDescription", "samples", "sampleColumn",
            "barcodeColumn", "groupColumn", "defaultObservationKind", "observationKindColumn", "spatial", "assays"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decodeIfPresent(Int.self, forKey: .schemaVersion) ?? 1
        id = try v.decode(String.self, forKey: .id)
        evidence = try v.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        samples = try v.decode([VivoOmicsSample].self, forKey: .samples)
        sampleColumn = try v.decode(String.self, forKey: .sampleColumn)
        barcodeColumn = try v.decodeIfPresent(String.self, forKey: .barcodeColumn)
        groupColumn = try v.decodeIfPresent(String.self, forKey: .groupColumn)
        defaultObservationKind = try v.decodeIfPresent(VivoObservationKind.self, forKey: .defaultObservationKind) ?? .cell
        observationKindColumn = try v.decodeIfPresent(String.self, forKey: .observationKindColumn)
        spatial = try v.decodeIfPresent([VivoH5MUQuantitativeSpatialMapping].self, forKey: .spatial) ?? []
        assays = try v.decode([VivoH5MUQuantitativeAssayMapping].self, forKey: .assays)
    }
}

public enum VivoQuantitativeH5MU {
    public static func writeSnapshot(_ data: VivoQuantitativeAssayDataset, to url: URL) throws {
        try data.validate()
        guard url.isFileURL else { throw VivoOmicsError.invalid("quantitative H5MU requires local output") }
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
                try h.check(write(d, type, 0, 0, 0, buffer), "quantitative H5MU write")
                try h.encoding(d, "array", "0.2.0")
            }
            func strings(_ parent: Int64, _ name: String, _ values: [String], scalar: Bool = false) throws {
                try h.writeStrings(parent, name, values, scalar: scalar)
                let d = try h.dataset(parent, name); defer { h.close(d, "H5Dclose") }
                try h.encoding(d, scalar ? "string" : "string-array", "0.2.0")
            }
            func frame(_ parent: Int64, _ name: String, index: [String], columns: [(String, [String])]) throws {
                let g = try group(parent, name, encoding: "dataframe"); defer { h.close(g, "H5Gclose") }
                try h.writeStrings(g, "_index", ["_index"], attribute: true, scalar: true)
                try h.writeStrings(g, "column-order", columns.map(\.0), attribute: true)
                try strings(g, "_index", index)
                for (column, values) in columns { guard values.count == index.count else { throw VivoOmicsError.invalid("quantitative H5MU dataframe column length") }; try strings(g, column, values) }
            }
            let enumCreate: @convention(c) (Int64) -> Int64 = try h.symbol("H5Tenum_create")
            let boolean = try h.id(enumCreate(h.native("NATIVE_UCHAR")), "quantitative boolean type")
            defer { h.close(boolean, "H5Tclose") }
            let insert: @convention(c) (Int64, UnsafePointer<CChar>, UnsafeRawPointer) -> Int32 = try h.symbol("H5Tenum_insert")
            var zero: UInt8 = 0, one: UInt8 = 1
            try h.check(insert(boolean, "FALSE", &zero), "quantitative boolean false")
            try h.check(insert(boolean, "TRUE", &one), "quantitative boolean true")
            func nullable(_ parent: Int64, _ name: String, values: [String?]) throws {
                let g = try group(parent, name, encoding: nil); defer { h.close(g, "H5Gclose") }
                try h.encoding(g, "nullable-string-array", "0.1.0")
                try strings(g, "values", values.map { $0 ?? "" })
                let mask = values.map { $0 == nil ? UInt8(1) : UInt8(0) }
                try mask.withUnsafeBytes { try raw(g, "mask", type: boolean, shape: [UInt64(mask.count)], buffer: $0.baseAddress) }
            }
            func mask(_ parent: Int64, _ name: String, values: [UInt64]) throws {
                let bytes = values.map { $0 == 0 ? UInt8(0) : UInt8(1) }
                try bytes.withUnsafeBytes { try raw(parent, name, type: boolean, shape: [UInt64(bytes.count), 1], buffer: $0.baseAddress) }
            }
            func observationColumns(_ rows: [Int]) -> [(String, [String])] {
                [("barcode", rows.map { data.observations[$0].identity.barcode }),
                 ("sample", rows.map { data.observations[$0].identity.sampleID }),
                 ("observation_kind", rows.map { data.observations[$0].kind.rawValue })]
            }

            let observationNames = data.observations.indices.map { "observation-\($0)" }
            try frame(file, "obs", index: observationNames, columns: observationColumns(Array(data.observations.indices)))
            let rootObs = try h.object(file, "obs"); defer { h.close(rootObs, "H5Oclose") }
            try nullable(rootObs, "group", values: data.observations.map { $0.identity.group })
            let delete: @convention(c) (Int64, UnsafePointer<CChar>) -> Int32 = try h.symbol("H5Adelete")
            try h.check(delete(rootObs, "column-order"), "replace quantitative observation column order")
            try h.writeStrings(rootObs, "column-order", ["barcode", "sample", "observation_kind", "group"], attribute: true)

            let allFeatures = data.assays.flatMap(\.features)
            let globalFeatureSpaces = data.assays.flatMap { assay in assay.features.map { _ in assay.featureNamespace } }
            let globalAssayKinds = data.assays.flatMap { assay in assay.features.map { _ in assay.kind.rawValue } }
            let globalUnits = data.assays.flatMap { assay in assay.features.map { _ in assay.unit } }
            try frame(file, "var", index: allFeatures.map(\.id), columns: [
                ("feature_space", globalFeatureSpaces), ("assay_kind", globalAssayKinds), ("unit", globalUnits)])
            let mod = try group(file, "mod", encoding: nil); defer { h.close(mod, "H5Gclose") }
            try h.writeStrings(mod, "mod-order", data.assays.map(\.id), attribute: true)
            let obsmap = try group(file, "obsmap"), varmap = try group(file, "varmap")
            let obsm = try group(file, "obsm"), varm = try group(file, "varm")
            defer { for g in [obsmap, varmap, obsm, varm] { h.close(g, "H5Gclose") } }
            for name in ["obsp", "varp"] { h.close(try group(file, name), "H5Gclose") }

            struct SpatialArray: Encodable { let frameID: String; let path: String }
            var spatialArrays: [SpatialArray] = [], occupied = Set(data.assays.map(\.id))
            for frame in data.spatialFrames {
                var name = "spatial", suffix = 0
                while occupied.contains(name) { suffix += 1; name = "spatial_" + String(suffix) }
                occupied.insert(name)
                var coordinates = [Double](repeating: .nan, count: data.observations.count * frame.axes.count)
                for (row, observation) in data.observations.enumerated() where observation.position?.frameID == frame.id {
                    guard let position = observation.position, position.coordinates.count == frame.axes.count else { throw VivoOmicsError.invalid("quantitative spatial frame shape") }
                    for axis in frame.axes.indices { coordinates[row * frame.axes.count + axis] = position.coordinates[axis] }
                }
                try coordinates.withUnsafeBytes { try raw(obsm, name, type: h.native("NATIVE_DOUBLE"), shape: [UInt64(data.observations.count), UInt64(frame.axes.count)], buffer: $0.baseAddress) }
                spatialArrays.append(.init(frameID: frame.id, path: "obsm/" + name))
            }

            var featureOffset = 0
            for assay in data.assays {
                try Task.checkCancellation()
                let modality = try group(mod, assay.id, encoding: nil); defer { h.close(modality, "H5Gclose") }
                try h.encoding(modality, "anndata", "0.1.0")
                try frame(modality, "obs", index: assay.observationIndices.map { observationNames[$0] }, columns: observationColumns(assay.observationIndices))
                let intervalContig = assay.features.map { $0.interval?.contig ?? "" }
                let intervalStart = assay.features.map { $0.interval.map { String($0.start) } ?? "" }
                let intervalEnd = assay.features.map { $0.interval.map { String($0.end) } ?? "" }
                let intervalPresent = assay.features.map { $0.interval == nil ? "false" : "true" }
                try frame(modality, "var", index: assay.features.map(\.id), columns: [
                    ("name", assay.features.map(\.name)),
                    ("feature_namespace", assay.features.map { _ in assay.featureNamespace }),
                    ("assay_kind", assay.features.map { _ in assay.kind.rawValue }),
                    ("unit", assay.features.map { _ in assay.unit }),
                    ("interval_contig", intervalContig), ("interval_start", intervalStart),
                    ("interval_end", intervalEnd), ("interval_present", intervalPresent)])
                let localObs = try h.object(file, "mod/\(assay.id)/obs"); defer { h.close(localObs, "H5Oclose") }
                try nullable(localObs, "group", values: assay.observationIndices.map { data.observations[$0].identity.group })
                try h.check(delete(localObs, "column-order"), "replace quantitative local column order")
                try h.writeStrings(localObs, "column-order", ["barcode", "sample", "observation_kind", "group"], attribute: true)
                let x = try group(modality, "X", encoding: nil); defer { h.close(x, "H5Gclose") }
                try h.encoding(x, "csr_matrix", "0.1.0")
                try h.writeIntegers(x, "shape", [UInt64(assay.matrix.observationCount), UInt64(assay.matrix.featureCount)], attribute: true)
                try h.writeIndices(x, "indptr", assay.matrix.rowOffsets)
                try h.writeIndices(x, "indices", assay.matrix.featureIndices)
                try h.writeDoubles(x, "data", assay.matrix.values)
                for name in ["layers", "obsm", "obsp", "varm", "varp", "uns"] { h.close(try group(modality, name), "H5Gclose") }
                var om = [UInt64](repeating: 0, count: data.observations.count)
                for (local, global) in assay.observationIndices.enumerated() { om[global] = UInt64(local + 1) }
                var vm = [UInt64](repeating: 0, count: allFeatures.count)
                for i in assay.features.indices { vm[featureOffset + i] = UInt64(i + 1) }
                featureOffset += assay.features.count
                try h.writeIntegers(obsmap, assay.id, om); try h.writeIntegers(varmap, assay.id, vm)
                try mask(obsm, assay.id, values: om); try mask(varm, assay.id, values: vm)
            }

            struct AssayMetadata: Encodable {
                let id: String; let kind: VivoQuantitativeAssayKind; let featureNamespace: String
                let unit: String; let sourceDescription: String; let features: [VivoAssayFeature]
                let observationIndices: [Int]
            }
            struct Metadata: Encodable {
                let schema = "numivivo.org/quantitative-h5mu/v1"
                let id: String; let evidence: VivoOmicsEvidence; let sourceDescription: String
                let samples: [VivoOmicsSample]; let observations: [VivoAssayObservation]
                let spatialFrames: [VivoSpatialFrame]; let spatialArrays: [SpatialArray]?
                let assays: [AssayMetadata]
            }
            let metadata = Metadata(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
                samples: data.samples, observations: data.observations, spatialFrames: data.spatialFrames,
                spatialArrays: spatialArrays.isEmpty ? nil : spatialArrays,
                assays: data.assays.map { .init(id: $0.id, kind: $0.kind, featureNamespace: $0.featureNamespace,
                    unit: $0.unit, sourceDescription: $0.sourceDescription, features: $0.features,
                    observationIndices: $0.observationIndices) })
            let uns = try group(file, "uns"); defer { h.close(uns, "H5Gclose") }
            let json = try String(decoding: VivoCanonicalJSON.encode(metadata), as: UTF8.self)
            try strings(uns, "numivivo", [json], scalar: true)
            try h.projectionStorageLimit(file)
            let flush: @convention(c) (Int64, Int32) -> Int32 = try h.symbol("H5Fflush")
            try h.check(flush(file, 1), "flush quantitative H5MU")
        }
    }

    public static func write(_ data: VivoQuantitativeAssayDataset, to destination: URL) throws {
        guard destination.isFileURL, !FileManager.default.fileExists(atPath: destination.path) else {
            throw VivoOmicsError.invalid("quantitative H5MU destination exists or is not local")
        }
        try writeSnapshot(data, to: destination)
    }
}

/// Native continuous-valued H5MU import.  The scanner retains explicit zero
/// values, treats NaN in dense arrays as missing, and rejects infinities.
public enum VivoQuantitativeH5MUImport {
    public static func validate(_ plan: VivoH5MUQuantitativePlan) throws {
        guard plan.schemaVersion == 1, vivoOmicsID(plan.id), !plan.sourceDescription.isEmpty,
              plan.sourceDescription.utf8.count <= 16_384, !plan.samples.isEmpty,
              plan.samples.count <= VivoQuantitativeAssayDataset.limits.maximumCells,
              !plan.assays.isEmpty, plan.assays.count <= 16,
              Set(plan.assays.map(\.id)).count == plan.assays.count,
              Set(plan.assays.map(\.sourceName)).count == plan.assays.count else {
            throw VivoOmicsError.invalid("quantitative H5MU plan bounds or identities")
        }
        try vivoH5ADComponent(plan.sampleColumn)
        for column in [plan.barcodeColumn, plan.groupColumn, plan.observationKindColumn].compactMap({ $0 }) { try vivoH5ADComponent(column) }
        try vivoOmicsValidateIdentities(samples: plan.samples, cells: [])
        var spatialIDs = Set<String>(), spatialPaths = Set<String>()
        for spatial in plan.spatial {
            guard spatial.path.hasPrefix("obsm/"), spatial.path.split(separator: "/", omittingEmptySubsequences: false).count == 2,
                  spatialIDs.insert(spatial.frame.id).inserted, spatialPaths.insert(spatial.path).inserted else {
                throw VivoOmicsError.invalid("quantitative H5MU spatial mapping")
            }
            try vivoH5ADComponent(String(spatial.path.dropFirst(5)))
            let frame = spatial.frame
            guard vivoOmicsID(frame.id), (2...3).contains(frame.axes.count), Set(frame.axes).count == frame.axes.count,
                  frame.axes.allSatisfy(vivoOmicsID), !frame.sourceDescription.isEmpty,
                  frame.sourceDescription.utf8.count <= 16_384 else { throw VivoOmicsError.invalid("quantitative H5MU spatial frame") }
        }
        for assay in plan.assays {
            try vivoH5ADComponent(assay.sourceName); try vivoH5ADComponent(assay.id)
            guard vivoOmicsID(assay.featureNamespace), vivoOmicsID(assay.unit),
                  !assay.sourceDescription.isEmpty, assay.sourceDescription.utf8.count <= 16_384,
                  assay.matrixPath == "X" || (assay.matrixPath.hasPrefix("layers/") && assay.matrixPath.split(separator: "/", omittingEmptySubsequences: false).count == 2) else {
                throw VivoOmicsError.invalid("quantitative H5MU assay mapping")
            }
            if assay.matrixPath != "X" { try vivoH5ADComponent(String(assay.matrixPath.dropFirst(7))) }
            for column in [assay.featureIDColumn, assay.featureNameColumn, assay.intervalContigColumn,
                           assay.intervalStartColumn, assay.intervalEndColumn, assay.intervalPresentColumn].compactMap({ $0 }) {
                try vivoH5ADComponent(column)
            }
            let intervalColumns = [assay.intervalContigColumn, assay.intervalStartColumn, assay.intervalEndColumn].compactMap({ $0 })
            guard intervalColumns.count == 0 || intervalColumns.count == 3 else { throw VivoOmicsError.invalid("quantitative interval columns must be complete") }
            guard assay.intervalPresentColumn == nil || intervalColumns.count == 3 else { throw VivoOmicsError.invalid("quantitative interval presence column requires coordinates") }
        }
    }

    public static func readSnapshot(_ source: URL, plan: VivoH5MUQuantitativePlan) throws -> VivoQuantitativeAssayDataset {
        try validate(plan)
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(source.path)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "MuData", try h.text(file, "encoding-version") == "0.1.0" else {
                throw VivoOmicsError.invalid("unsupported quantitative MuData root encoding")
            }
            let mod = try h.object(file, "mod"); defer { h.close(mod, "H5Oclose") }
            guard Set(try h.projectionNames(mod)) == Set(plan.assays.map(\.sourceName)) else {
                throw VivoOmicsError.invalid("every quantitative H5MU modality must be mapped exactly once")
            }
            let limits = VivoQuantitativeAssayDataset.limits, reader = VivoH5ADFrameReader(h: h, file: file)
            let observationIDs = try reader.index("obs", maximum: limits.maximumCells)
            let variableIDs = try reader.index("var", maximum: limits.maximumFeatures)
            guard Set(observationIDs).count == observationIDs.count else { throw VivoOmicsError.invalid("duplicate quantitative global observation index") }
            func global(_ name: String) throws -> [String?] {
                let values = try reader.column("obs", name, maximum: limits.maximumCells)
                guard values.count == observationIDs.count else { throw VivoOmicsError.invalid("quantitative observation column length") }
                return values
            }
            let samples = try reader.required(global(plan.sampleColumn))
            let barcodes = try plan.barcodeColumn.map { try reader.required(global($0)) } ?? observationIDs
            let groups = try plan.groupColumn.map { try global($0) }
            let kinds = try plan.observationKindColumn.map { try reader.required(global($0)) }

            var positions = [VivoSpatialPosition?](repeating: nil, count: observationIDs.count)
            for spatial in plan.spatial {
                let dataset = try h.dataset(file, spatial.path); defer { h.close(dataset, "H5Dclose") }
                let dimensions = try h.shape(dataset, attribute: false)
                guard dimensions == [UInt64(observationIDs.count), UInt64(spatial.frame.axes.count)],
                      try h.text(dataset, "encoding-type") == "array", try h.text(dataset, "encoding-version") == "0.2.0" else {
                    throw VivoOmicsError.invalid("quantitative spatial array shape or encoding")
                }
                let count = observationIDs.count * spatial.frame.axes.count
                let values = try h.doubles(dataset, maximum: count, allowInteger: true)
                for row in observationIDs.indices {
                    let slice = Array(values[(row * spatial.frame.axes.count)..<((row + 1) * spatial.frame.axes.count)])
                    if slice.allSatisfy(\.isNaN) { continue }
                    guard slice.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("quantitative spatial position is partial or nonfinite") }
                    guard positions[row] == nil else { throw VivoOmicsError.invalid("observation has multiple quantitative spatial frames") }
                    positions[row] = .init(frameID: spatial.frame.id, coordinates: slice)
                }
            }
            let observations: [VivoAssayObservation] = try observationIDs.indices.map { index in
                let kind: VivoObservationKind
                if let kinds { guard let parsed = VivoObservationKind(rawValue: kinds[index]) else { throw VivoOmicsError.invalid("unknown quantitative observation kind") }; kind = parsed }
                else { kind = plan.defaultObservationKind }
                return .init(identity: .init(barcode: barcodes[index], sampleID: samples[index], group: groups?[index]), kind: kind, position: positions[index])
            }
            try vivoOmicsValidateIdentities(samples: plan.samples, cells: observations.map(\.identity))

            func map(_ parent: String, name: String, global: [String], local: [String]) throws -> [Int] {
                let object = try h.object(file, parent); defer { h.close(object, "H5Oclose") }
                guard try h.text(object, "encoding-type") == "dict", try h.text(object, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("quantitative H5MU map encoding") }
                let dataset = try h.dataset(file, parent + "/" + name); defer { h.close(dataset, "H5Dclose") }
                let shape = try h.shape(dataset, attribute: false)
                guard shape == [UInt64(global.count)] || shape == [UInt64(global.count), 1] else { throw VivoOmicsError.invalid("quantitative H5MU map dimensions") }
                let values = try h.integers(dataset, maximum: global.count)
                var result = [Int](repeating: -1, count: local.count)
                for index in values.indices where values[index] != 0 {
                    guard values[index] <= UInt64(local.count) else { throw VivoOmicsError.invalid("quantitative H5MU map index out of range") }
                    let localIndex = Int(values[index]) - 1
                    guard result[localIndex] == -1, global[index] == local[localIndex] else { throw VivoOmicsError.invalid("quantitative H5MU map duplicate or identity mismatch") }
                    result[localIndex] = index
                }
                guard result.allSatisfy({ $0 >= 0 }) else { throw VivoOmicsError.invalid("quantitative H5MU map omits local rows") }
                return result
            }
            func parseBool(_ value: String) throws -> Bool {
                switch value.lowercased() { case "true", "1": return true; case "false", "0": return false; default: throw VivoOmicsError.invalid("quantitative interval presence value") }
            }

            var assays: [VivoQuantitativeAssaySpace] = [], remainingFeatures = limits.maximumFeatures
            var remainingEntries = limits.maximumNonzeros, remainingDenseWork = 500_000_000
            var coveredVariables = Set<Int>()
            for assay in plan.assays {
                try Task.checkCancellation()
                let base = "mod/" + assay.sourceName, modality = try h.object(file, base); defer { h.close(modality, "H5Oclose") }
                guard try h.text(modality, "encoding-type") == "anndata", try h.text(modality, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("quantitative H5MU modality must be AnnData") }
                let localObservations = try reader.index(base + "/obs", maximum: limits.maximumCells)
                let localVariables = try reader.index(base + "/var", maximum: remainingFeatures)
                remainingFeatures -= localVariables.count
                let rows = try map("obsmap", name: assay.sourceName, global: observationIDs, local: localObservations)
                let variables = try map("varmap", name: assay.sourceName, global: variableIDs, local: localVariables)
                guard variables.allSatisfy({ coveredVariables.insert($0).inserted }) else { throw VivoOmicsError.invalid("quantitative H5MU feature maps overlap") }
                let localObsObject = try h.object(file, base + "/obs"); defer { h.close(localObsObject, "H5Oclose") }
                let localColumns = Set(try h.projectionNames(localObsObject))
                for (name, expected) in [(Optional(plan.sampleColumn), samples.map(Optional.init)), (plan.barcodeColumn, barcodes.map(Optional.init)),
                                          (plan.groupColumn, groups ?? []), (plan.observationKindColumn, observations.map { Optional($0.kind.rawValue) })] {
                    if let name, localColumns.contains(name) {
                        let values = try reader.column(base + "/obs", name, maximum: limits.maximumCells)
                        guard values.count == rows.count, values == rows.map({ expected[$0] }) else { throw VivoOmicsError.invalid("quantitative local/global observation metadata conflict") }
                    }
                }
                let localVarObject = try h.object(file, base + "/var"); defer { h.close(localVarObject, "H5Oclose") }
                let featureColumns = Set(try h.projectionNames(localVarObject))
                for (name, expected) in [("feature_namespace", assay.featureNamespace), ("assay_kind", assay.kind.rawValue), ("unit", assay.unit)] where featureColumns.contains(name) {
                    let values = try reader.required(reader.column(base + "/var", name, maximum: limits.maximumFeatures))
                    guard values.count == localVariables.count, values.allSatisfy({ $0 == expected }) else { throw VivoOmicsError.invalid("quantitative feature metadata conflict") }
                }
                let ids = try assay.featureIDColumn.map { try reader.required(reader.column(base + "/var", $0, maximum: limits.maximumFeatures)) } ?? localVariables
                let names = try assay.featureNameColumn.map { try reader.required(reader.column(base + "/var", $0, maximum: limits.maximumFeatures)) } ?? localVariables
                guard ids.count == localVariables.count, names.count == localVariables.count else { throw VivoOmicsError.invalid("quantitative feature column length") }
                var intervalContigs: [String?]?, intervalStarts: [String?]?, intervalEnds: [String?]?, intervalPresence: [Bool]?
                if let contigColumn = assay.intervalContigColumn, let startColumn = assay.intervalStartColumn, let endColumn = assay.intervalEndColumn {
                    let contigs = try reader.column(base + "/var", contigColumn, maximum: limits.maximumFeatures)
                    let starts = try reader.column(base + "/var", startColumn, maximum: limits.maximumFeatures)
                    let ends = try reader.column(base + "/var", endColumn, maximum: limits.maximumFeatures)
                    guard contigs.count == localVariables.count, starts.count == localVariables.count, ends.count == localVariables.count else { throw VivoOmicsError.invalid("quantitative interval column length") }
                    intervalContigs = contigs; intervalStarts = starts; intervalEnds = ends
                    if let presenceColumn = assay.intervalPresentColumn {
                        let values = try reader.required(reader.column(base + "/var", presenceColumn, maximum: limits.maximumFeatures))
                        guard values.count == localVariables.count else { throw VivoOmicsError.invalid("quantitative interval presence length") }
                        intervalPresence = try values.map(parseBool)
                    }
                }
                let features = try ids.indices.map { index in
                    guard let contigs = intervalContigs, let starts = intervalStarts, let ends = intervalEnds else {
                        return VivoAssayFeature(id: ids[index], name: names[index], interval: nil)
                    }
                    let present = intervalPresence?[index] ?? (contigs[index] != nil || starts[index] != nil || ends[index] != nil)
                    guard present else {
                        guard contigs[index] == nil || contigs[index] == "", starts[index] == nil || starts[index] == "", ends[index] == nil || ends[index] == "" else { throw VivoOmicsError.invalid("quantitative absent interval has values") }
                        return VivoAssayFeature(id: ids[index], name: names[index], interval: nil)
                    }
                    guard let contig = contigs[index], !contig.isEmpty, let startText = starts[index], let endText = ends[index],
                          let start = UInt64(startText), let end = UInt64(endText), start < end, vivoOmicsID(contig) else {
                        throw VivoOmicsError.invalid("quantitative interval value")
                    }
                    return VivoAssayFeature(id: ids[index], name: names[index], interval: .init(contig: contig, start: start, end: end))
                }
                let path = base + "/" + assay.matrixPath, matrix = try h.object(file, path); defer { h.close(matrix, "H5Oclose") }
                let encoding = try h.text(matrix, "encoding-type")
                var sourceEntries = 0
                if encoding == "array" {
                    guard try h.text(matrix, "encoding-version") == "0.2.0", localVariables.count > 0,
                          localObservations.count <= remainingDenseWork / localVariables.count,
                          localObservations.count * localVariables.count <= remainingEntries else {
                        throw VivoOmicsError.limit("quantitative H5MU dense bounds")
                    }
                    sourceEntries = localObservations.count * localVariables.count
                    remainingDenseWork -= sourceEntries
                } else {
                    guard encoding == "csr_matrix" || encoding == "csc_matrix", try h.text(matrix, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("quantitative matrix requires CSR, CSC or dense array") }
                    let shapeAttribute = try h.attribute(matrix, "shape"); defer { h.close(shapeAttribute, "H5Aclose") }
                    guard try h.shape(shapeAttribute, attribute: true) == [2], try h.integers(shapeAttribute, attribute: true, maximum: 2) == [UInt64(localObservations.count), UInt64(localVariables.count)] else { throw VivoOmicsError.invalid("quantitative sparse matrix shape") }
                    let data = try h.dataset(file, path + "/data"); defer { h.close(data, "H5Dclose") }
                    let dimensions = try h.shape(data, attribute: false)
                    guard dimensions.count == 1, dimensions[0] <= UInt64(remainingEntries) else { throw VivoOmicsError.limit("quantitative sparse source entries") }
                    sourceEntries = Int(dimensions[0])
                }
                var valuesByRow = [[Int: Double]](repeating: [:], count: rows.count), admitted = 0
                func admit(_ row: Int, _ feature: Int, _ value: Double) throws {
                    guard value.isFinite, row >= 0, row < rows.count, feature >= 0, feature < features.count else { throw VivoOmicsError.invalid("quantitative matrix value or index") }
                    guard admitted < remainingEntries else { throw VivoOmicsError.limit("quantitative H5MU aggregate entries") }
                    let previous = valuesByRow[row].updateValue(value, forKey: feature)
                    if let previous {
                        let sum = previous + value
                        guard sum.isFinite else { throw VivoOmicsError.invalid("quantitative duplicate accumulation is nonfinite") }
                        valuesByRow[row][feature] = sum
                    }
                    admitted += 1
                }
                if encoding == "array" {
                    let dataset = try h.dataset(file, path); defer { h.close(dataset, "H5Dclose") }
                    for row in 0..<localObservations.count {
                        try Task.checkCancellation()
                        let values = try h.doubles(dataset, maximum: limits.maximumFeatures, allowInteger: true, row: row)
                        for (feature, value) in values.enumerated() {
                            if value.isNaN { continue }
                            try admit(row, feature, value)
                        }
                    }
                } else {
                    let indexDataset = try h.dataset(file, path + "/indices"), dataDataset = try h.dataset(file, path + "/data")
                    defer { h.close(indexDataset, "H5Dclose"); h.close(dataDataset, "H5Dclose") }
                    let offsetDataset = try h.dataset(file, path + "/indptr"); defer { h.close(offsetDataset, "H5Dclose") }
                    let isCSR = encoding == "csr_matrix", major = isCSR ? localObservations.count : localVariables.count, minor = isCSR ? localVariables.count : localObservations.count
                    guard try h.shape(offsetDataset, attribute: false).count == 1,
                          try h.shape(indexDataset, attribute: false) == [UInt64(sourceEntries)],
                          try h.shape(dataDataset, attribute: false) == [UInt64(sourceEntries)] else {
                        throw VivoOmicsError.invalid("quantitative sparse arrays must be vectors with matching lengths")
                    }
                    let offsets = try h.integers(offsetDataset, maximum: major + 1)
                    guard offsets.count == major + 1, offsets.first == 0, offsets.last == UInt64(sourceEntries) else { throw VivoOmicsError.invalid("quantitative sparse offsets") }
                    for majorIndex in 0..<major {
                        try Task.checkCancellation()
                        guard offsets[majorIndex] <= offsets[majorIndex + 1], offsets[majorIndex + 1] <= UInt64(sourceEntries) else { throw VivoOmicsError.invalid("quantitative sparse offset range") }
                        var cursor = Int(offsets[majorIndex]), canonical: [Int: Double] = [:]
                        let end = Int(offsets[majorIndex + 1])
                        while cursor < end {
                            let next = min(end, cursor + 65_536)
                            let indices = try h.integers(indexDataset, maximum: 65_536, range: cursor..<next)
                            let values = try h.doubles(dataDataset, maximum: 65_536, allowInteger: true, range: cursor..<next)
                            guard indices.count == values.count else { throw VivoOmicsError.invalid("quantitative sparse slice shape") }
                            for i in indices.indices {
                                guard indices[i] < UInt64(minor), values[i].isFinite else { throw VivoOmicsError.invalid("quantitative sparse index or value") }
                                let key = Int(indices[i]), previous = canonical[key] ?? 0, sum = previous + values[i]
                                guard sum.isFinite else { throw VivoOmicsError.invalid("quantitative duplicate accumulation is nonfinite") }
                                canonical[key] = sum
                            }
                            cursor = next
                        }
                        for key in canonical.keys.sorted() {
                            if isCSR { try admit(majorIndex, key, canonical[key]!) }
                            else { try admit(key, majorIndex, canonical[key]!) }
                        }
                    }
                }
                remainingEntries -= max(sourceEntries, admitted)
                var offsets = [0], indices: [Int] = [], values: [Double] = []
                for row in valuesByRow {
                    for feature in row.keys.sorted() { indices.append(feature); values.append(row[feature]!) }
                    offsets.append(values.count)
                }
                assays.append(.init(id: assay.id, kind: assay.kind, featureNamespace: assay.featureNamespace,
                    unit: assay.unit, sourceDescription: assay.sourceDescription, features: features,
                    observationIndices: rows, matrix: .init(observationCount: rows.count, featureCount: features.count,
                        rowOffsets: offsets, featureIndices: indices, values: values)))
            }
            guard coveredVariables.count == variableIDs.count else { throw VivoOmicsError.invalid("unmapped quantitative H5MU features") }
            let result = VivoQuantitativeAssayDataset(id: plan.id, evidence: plan.evidence,
                sourceDescription: plan.sourceDescription, samples: plan.samples, observations: observations,
                spatialFrames: plan.spatial.map(\.frame), assays: assays)
            try result.validate(); return result
        }
    }
}
