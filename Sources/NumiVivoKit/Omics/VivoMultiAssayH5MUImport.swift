import Foundation

public struct VivoH5MUAssayMapping: Codable, Sendable, Equatable {
    public let sourceName: String
    public let id: String
    public let kind: VivoAssayKind
    public let featureNamespace: String
    public let countUnit: VivoAssayCountUnit
    public let genomeAssembly: String?
    public let matrixPath: String
    public let featureIDColumn: String?
    public let featureNameColumn: String?
    public let peakIDConvention: String?
    public init(sourceName: String, id: String, kind: VivoAssayKind, featureNamespace: String, countUnit: VivoAssayCountUnit, genomeAssembly: String?, matrixPath: String, featureIDColumn: String?, featureNameColumn: String?, peakIDConvention: String?) {
        self.sourceName = sourceName
        self.id = id
        self.kind = kind
        self.featureNamespace = featureNamespace
        self.countUnit = countUnit
        self.genomeAssembly = genomeAssembly
        self.matrixPath = matrixPath
        self.featureIDColumn = featureIDColumn
        self.featureNameColumn = featureNameColumn
        self.peakIDConvention = peakIDConvention
    }
    private enum CodingKeys: String, CodingKey { case sourceName, id, kind, featureNamespace, countUnit, genomeAssembly, matrixPath, featureIDColumn, featureNameColumn, peakIDConvention }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["sourceName", "id", "kind", "featureNamespace", "countUnit", "genomeAssembly", "matrixPath", "featureIDColumn", "featureNameColumn", "peakIDConvention"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        sourceName = try v.decode(String.self, forKey: .sourceName)
        id = try v.decode(String.self, forKey: .id)
        kind = try v.decode(VivoAssayKind.self, forKey: .kind)
        featureNamespace = try v.decode(String.self, forKey: .featureNamespace)
        countUnit = try v.decode(VivoAssayCountUnit.self, forKey: .countUnit)
        genomeAssembly = try v.decodeIfPresent(String.self, forKey: .genomeAssembly)
        matrixPath = try v.decode(String.self, forKey: .matrixPath)
        featureIDColumn = try v.decodeIfPresent(String.self, forKey: .featureIDColumn)
        featureNameColumn = try v.decodeIfPresent(String.self, forKey: .featureNameColumn)
        peakIDConvention = try v.decodeIfPresent(String.self, forKey: .peakIDConvention)
    }
}

public struct VivoH5MUSpatialMapping: Codable, Sendable, Equatable {
    public let path: String
    public let frame: VivoSpatialFrame
    public init(path: String, frame: VivoSpatialFrame) {
        self.path = path
        self.frame = frame
    }
    private enum CodingKeys: String, CodingKey { case path, frame }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["path", "frame"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        path = try v.decode(String.self, forKey: .path)
        frame = try v.decode(VivoSpatialFrame.self, forKey: .frame)
    }
}

public struct VivoH5MUMultiAssayPlan: Codable, Sendable, Equatable {
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
    public let spatial: VivoH5MUSpatialMapping?
    public let assays: [VivoH5MUAssayMapping]
    public init(schemaVersion: Int, id: String, evidence: VivoOmicsEvidence, sourceDescription: String, samples: [VivoOmicsSample], sampleColumn: String, barcodeColumn: String?, groupColumn: String?, defaultObservationKind: VivoObservationKind, observationKindColumn: String?, spatial: VivoH5MUSpatialMapping?, assays: [VivoH5MUAssayMapping]) {
        self.schemaVersion = schemaVersion
        self.id = id
        self.evidence = evidence
        self.sourceDescription = sourceDescription
        self.samples = samples
        self.sampleColumn = sampleColumn
        self.barcodeColumn = barcodeColumn
        self.groupColumn = groupColumn
        self.defaultObservationKind = defaultObservationKind
        self.observationKindColumn = observationKindColumn
        self.spatial = spatial
        self.assays = assays
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, id, evidence, sourceDescription, samples, sampleColumn, barcodeColumn, groupColumn, defaultObservationKind, observationKindColumn, spatial, assays }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "id", "evidence", "sourceDescription", "samples", "sampleColumn", "barcodeColumn", "groupColumn", "defaultObservationKind", "observationKindColumn", "spatial", "assays"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decode(Int.self, forKey: .schemaVersion)
        id = try v.decode(String.self, forKey: .id)
        evidence = try v.decode(VivoOmicsEvidence.self, forKey: .evidence)
        sourceDescription = try v.decode(String.self, forKey: .sourceDescription)
        samples = try v.decode([VivoOmicsSample].self, forKey: .samples)
        sampleColumn = try v.decode(String.self, forKey: .sampleColumn)
        barcodeColumn = try v.decodeIfPresent(String.self, forKey: .barcodeColumn)
        groupColumn = try v.decodeIfPresent(String.self, forKey: .groupColumn)
        defaultObservationKind = try v.decode(VivoObservationKind.self, forKey: .defaultObservationKind)
        observationKindColumn = try v.decodeIfPresent(String.self, forKey: .observationKindColumn)
        spatial = try v.decodeIfPresent(VivoH5MUSpatialMapping.self, forKey: .spatial)
        assays = try v.decode([VivoH5MUAssayMapping].self, forKey: .assays)
    }
}

public enum VivoMultiAssayH5MUImport {
    public static func validate(_ plan: VivoH5MUMultiAssayPlan) throws {
        guard plan.schemaVersion == 1, vivoOmicsID(plan.id), !plan.sourceDescription.isEmpty,
              plan.sourceDescription.utf8.count <= 16_384, !plan.assays.isEmpty, plan.assays.count <= 16,
              !plan.samples.isEmpty, plan.samples.count <= VivoMultiAssayDataset.limits.maximumCells,
              Set(plan.assays.map(\.id)).count == plan.assays.count,
              Set(plan.assays.map(\.sourceName)).count == plan.assays.count else { throw VivoOmicsError.invalid("H5MU plan bounds or identities") }
        try vivoH5ADComponent(plan.sampleColumn)
        for c in [plan.barcodeColumn, plan.groupColumn, plan.observationKindColumn].compactMap({ $0 }) { try vivoH5ADComponent(c) }
        try vivoOmicsValidateIdentities(samples: plan.samples, cells: [])
        for a in plan.assays {
            try vivoH5ADComponent(a.sourceName); try vivoH5ADComponent(a.id)
            guard vivoOmicsID(a.featureNamespace), a.genomeAssembly.map(vivoOmicsID) ?? true,
                  a.matrixPath == "X" || (a.matrixPath.hasPrefix("layers/") && a.matrixPath.split(separator: "/", omittingEmptySubsequences: false).count == 2) else { throw VivoOmicsError.invalid("H5MU assay mapping") }
            if a.matrixPath != "X" { try vivoH5ADComponent(String(a.matrixPath.dropFirst(7))) }
            for c in [a.featureIDColumn, a.featureNameColumn].compactMap({ $0 }) { try vivoH5ADComponent(c) }
            if a.kind == .chromatinAccessibility {
                guard a.genomeAssembly != nil, a.countUnit != .umiCount,
                      a.peakIDConvention == "contig:start-end:zero-based-half-open" else { throw VivoOmicsError.invalid("H5MU peak identity convention/assembly/unit required") }
            } else {
                guard a.peakIDConvention == nil, a.countUnit != .fragmentCount else { throw VivoOmicsError.invalid("H5MU non-peak mapping has peak convention/units") }
            }
        }
        if let spatial = plan.spatial {
            guard spatial.path.hasPrefix("obsm/"), spatial.path.split(separator: "/", omittingEmptySubsequences: false).count == 2 else { throw VivoOmicsError.invalid("H5MU spatial path must name obsm array") }
            try vivoH5ADComponent(String(spatial.path.dropFirst(5)))
            let f = spatial.frame
            guard vivoOmicsID(f.id), (2...3).contains(f.axes.count), Set(f.axes).count == f.axes.count,
                  f.axes.allSatisfy(vivoOmicsID), !f.sourceDescription.isEmpty, f.sourceDescription.utf8.count <= 16_384 else { throw VivoOmicsError.invalid("H5MU spatial frame") }
        }
    }
    public static func readSnapshot(_ source: URL, plan: VivoH5MUMultiAssayPlan) throws -> VivoMultiAssayDataset {
        try validate(plan)
        return try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(source.path)
            defer { h.close(file, "H5Fclose") }
            guard try h.text(file, "encoding-type") == "MuData", try h.text(file, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("unsupported MuData root encoding") }
            let attributes = try h.projectionNames(file, attributes: true)
            if attributes.contains("axis") {
                let axis = try h.attribute(file, "axis"); defer { h.close(axis, "H5Aclose") }
                guard try h.integers(axis, attribute: true, maximum: 1) == [0] else { throw VivoOmicsError.invalid("H5MU requires shared-observation axis 0") }
            }
            let mod = try h.object(file, "mod"); defer { h.close(mod, "H5Oclose") }
            guard Set(try h.projectionNames(mod)) == Set(plan.assays.map(\.sourceName)) else { throw VivoOmicsError.invalid("every H5MU modality must be mapped exactly once") }
            let limits = VivoMultiAssayDataset.limits, reader = VivoH5ADFrameReader(h: h, file: file)
            let obsIDs = try reader.index("obs", maximum: limits.maximumCells)
            let varIDs = try reader.index("var", maximum: limits.maximumFeatures)
            guard Set(obsIDs).count == obsIDs.count else { throw VivoOmicsError.invalid("duplicate global H5MU observation index") }
            func global(_ name: String) throws -> [String?] {
                let v = try reader.column("obs", name, maximum: limits.maximumCells)
                guard v.count == obsIDs.count else { throw VivoOmicsError.invalid("global observation column length") }; return v
            }
            let samples = try reader.required(global(plan.sampleColumn))
            let barcodes = try plan.barcodeColumn.map { try reader.required(global($0)) } ?? obsIDs
            let groups = try plan.groupColumn.map { try global($0) }
            let kinds = try plan.observationKindColumn.map { try reader.required(global($0)) }
            var positions = [VivoSpatialPosition?](repeating: nil, count: obsIDs.count)
            if let spatial = plan.spatial {
                let d = try h.dataset(file, spatial.path); defer { h.close(d, "H5Dclose") }
                let dims = try h.shape(d, attribute: false)
                guard dims == [UInt64(obsIDs.count), UInt64(spatial.frame.axes.count)],
                      try h.text(d, "encoding-type") == "array", try h.text(d, "encoding-version") == "0.2.0" else { throw VivoOmicsError.invalid("spatial array shape or encoding") }
                let n = obsIDs.count * spatial.frame.axes.count
                let t = try h.type(d, attribute: false); defer { h.close(t, "H5Tclose") }
                let kind: @convention(c) (Int64) -> Int32 = try h.symbol("H5Tget_class")
                var values: [Double]
                if kind(t) == 0 {
                    let v = try h.categoricalCodes(d, maximum: n)
                    guard v.allSatisfy({ $0 >= -9_007_199_254_740_992 && $0 <= 9_007_199_254_740_992 }) else { throw VivoOmicsError.invalid("spatial integers exceed exact Double range") }
                    values = v.map(Double.init)
                } else {
                    guard kind(t) == 1 else { throw VivoOmicsError.invalid("spatial array must be numeric") }
                    values = [Double](repeating: 0, count: n)
                    try values.withUnsafeMutableBytes { try h.read(d, type: h.native("NATIVE_DOUBLE"), attribute: false, into: $0.baseAddress) }
                }
                for row in obsIDs.indices {
                    let v = Array(values[(row * spatial.frame.axes.count)..<((row + 1) * spatial.frame.axes.count)])
                    if v.allSatisfy(\.isNaN) { continue }
                    guard v.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("partially missing or infinite spatial position") }
                    positions[row] = .init(frameID: spatial.frame.id, coordinates: v)
                }
            }
            let observations: [VivoAssayObservation] = try obsIDs.indices.map { i in
                let kind: VivoObservationKind
                if let kinds { guard let parsed = VivoObservationKind(rawValue: kinds[i]) else { throw VivoOmicsError.invalid("unknown observation kind") }; kind = parsed }
                else { kind = plan.defaultObservationKind }
                return .init(identity: .init(barcode: barcodes[i], sampleID: samples[i], group: groups?[i]), kind: kind, position: positions[i])
            }
            try vivoOmicsValidateIdentities(samples: plan.samples, cells: observations.map(\.identity))
            func map(_ parent: String, name: String, global: [String], local: [String]) throws -> [Int] {
                let g = try h.object(file, parent); defer { h.close(g, "H5Oclose") }
                guard try h.text(g, "encoding-type") == "dict", try h.text(g, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("H5MU map encoding") }
                let d = try h.dataset(file, parent + "/" + name); defer { h.close(d, "H5Dclose") }
                let shape = try h.shape(d, attribute: false)
                guard shape == [UInt64(global.count)] || shape == [UInt64(global.count), 1] else { throw VivoOmicsError.invalid("H5MU map dimensions") }
                let values = try h.integers(d, maximum: global.count)
                var result = [Int](repeating: -1, count: local.count)
                for i in values.indices where values[i] != 0 {
                    guard values[i] <= UInt64(local.count) else { throw VivoOmicsError.invalid("H5MU map index out of range") }
                    let j = Int(values[i]) - 1
                    guard result[j] == -1, global[i] == local[j] else { throw VivoOmicsError.invalid("H5MU map duplicate or identity mismatch") }
                    result[j] = i
                }
                guard result.allSatisfy({ $0 >= 0 }) else { throw VivoOmicsError.invalid("H5MU map omits local rows") }
                return result
            }
            var assays: [VivoAssaySpace] = [], remainingEntries = limits.maximumNonzeros
            var remainingFeatures = limits.maximumFeatures, remainingDenseWork = 500_000_000
            var coveredVariables = Set<Int>()
            for a in plan.assays {
                try Task.checkCancellation()
                let base = "mod/" + a.sourceName, g = try h.object(file, base); defer { h.close(g, "H5Oclose") }
                guard try h.text(g, "encoding-type") == "anndata", try h.text(g, "encoding-version") == "0.1.0" else { throw VivoOmicsError.invalid("H5MU modality must be AnnData") }
                let localObs = try reader.index(base + "/obs", maximum: limits.maximumCells)
                let localVar = try reader.index(base + "/var", maximum: remainingFeatures)
                remainingFeatures -= localVar.count
                let rows = try map("obsmap", name: a.sourceName, global: obsIDs, local: localObs)
                let variables = try map("varmap", name: a.sourceName, global: varIDs, local: localVar)
                guard variables.allSatisfy({ coveredVariables.insert($0).inserted }) else { throw VivoOmicsError.invalid("shared H5MU feature maps require a different axis model") }
                // Selected identities in local obs may not contradict the global frame.
                let localFrame = try h.object(file, base + "/obs"); defer { h.close(localFrame, "H5Oclose") }
                let localColumns = Set(try h.projectionNames(localFrame))
                for (name, globalValues) in [(Optional(plan.sampleColumn), samples.map(Optional.init)), (plan.barcodeColumn, barcodes.map(Optional.init)), (plan.groupColumn, groups ?? []), (plan.observationKindColumn, observations.map { Optional($0.kind.rawValue) })] {
                    if let name, localColumns.contains(name) {
                        let values = try reader.column(base + "/obs", name, maximum: limits.maximumCells)
                        guard values.count == rows.count, values == rows.map({ globalValues[$0] }) else { throw VivoOmicsError.invalid("local/global H5MU observation metadata conflict") }
                    }
                }
                let featureFrame = try h.object(file, base + "/var"); defer { h.close(featureFrame, "H5Oclose") }
                let featureColumns = Set(try h.projectionNames(featureFrame))
                for (name, expected) in [("feature_namespace", a.featureNamespace), ("assay_kind", a.kind.rawValue)] where featureColumns.contains(name) {
                    let values = try reader.required(reader.column(base + "/var", name, maximum: limits.maximumFeatures))
                    guard values.count == localVar.count, values.allSatisfy({ $0 == expected }) else { throw VivoOmicsError.invalid("H5MU feature-space metadata conflict") }
                }
                if featureColumns.contains("genome") {
                    let values = try reader.column(base + "/var", "genome", maximum: limits.maximumFeatures)
                    guard values.count == localVar.count, values.allSatisfy({ $0 == nil || $0 == "" || $0 == a.genomeAssembly }) else { throw VivoOmicsError.invalid("H5MU source genome conflicts with plan") }
                }
                let ids = try a.featureIDColumn.map { try reader.required(reader.column(base + "/var", $0, maximum: limits.maximumFeatures)) } ?? localVar
                let names = try a.featureNameColumn.map { try reader.required(reader.column(base + "/var", $0, maximum: limits.maximumFeatures)) } ?? localVar
                guard ids.count == localVar.count, names.count == localVar.count else { throw VivoOmicsError.invalid("H5MU feature column length") }
                let features: [VivoAssayFeature] = try ids.indices.map { i in
                    var interval: VivoGenomicInterval?
                    if a.kind == .chromatinAccessibility {
                        let p = ids[i].split(separator: ":", omittingEmptySubsequences: false)
                        guard p.count == 2 else { throw VivoOmicsError.invalid("H5MU peak ID format") }
                        let b = p[1].split(separator: "-", omittingEmptySubsequences: false)
                        guard b.count == 2, let start = UInt64(b[0]), let end = UInt64(b[1]), start < end else { throw VivoOmicsError.invalid("H5MU peak interval") }
                        interval = .init(contig: String(p[0]), start: start, end: end)
                    }
                    return .init(id: ids[i], name: names[i], interval: interval)
                }
                let path = base + "/" + a.matrixPath, matrix = try h.object(file, path); defer { h.close(matrix, "H5Oclose") }
                var sourceEntries = 0
                if try h.text(matrix, "encoding-type") == "array" {
                    guard localVar.count > 0, localObs.count <= remainingDenseWork / localVar.count else { throw VivoOmicsError.limit("H5MU dense scan work") }
                    remainingDenseWork -= localObs.count * localVar.count
                } else {
                    let d = try h.dataset(file, path + "/data"); defer { h.close(d, "H5Dclose") }
                    let shape = try h.shape(d, attribute: false)
                    guard shape.count == 1, shape[0] <= UInt64(remainingEntries) else { throw VivoOmicsError.limit("H5MU aggregate source entries") }
                    sourceEntries = Int(shape[0])
                }
                var counts = [[Int: UInt64]](repeating: [:], count: rows.count), admitted = 0
                var selectedLimits = limits; selectedLimits.maximumNonzeros = remainingEntries
                try VivoH5ADCountReader.scan(h, file: file, path: path, rows: rows.count, features: features.count, limits: selectedLimits) { row, feature, value in
                    guard admitted < remainingEntries else { throw VivoOmicsError.limit("H5MU aggregate count entries") }
                    counts[row][feature] = value; admitted += 1
                }
                remainingEntries -= max(sourceEntries, admitted)
                var offsets = [0], columns: [Int] = [], values: [UInt64] = []
                for row in counts {
                    for j in row.keys.sorted() { columns.append(j); values.append(row[j]!) }; offsets.append(values.count)
                }
                assays.append(.init(id: a.id, kind: a.kind, featureNamespace: a.featureNamespace, countUnit: a.countUnit, genomeAssembly: a.genomeAssembly,
                    sourceDescription: plan.sourceDescription, features: features, observationIndices: rows,
                    matrix: .init(cellCount: rows.count, featureCount: features.count, rowOffsets: offsets, featureIndices: columns, counts: values)))
            }
            guard coveredVariables.count == varIDs.count else { throw VivoOmicsError.invalid("unmapped global H5MU features") }
            let result = VivoMultiAssayDataset(schemaVersion: 1, id: plan.id, evidence: plan.evidence, sourceDescription: plan.sourceDescription,
                samples: plan.samples, observations: observations, spatialFrames: plan.spatial.map { [$0.frame] } ?? [], assays: assays)
            try result.validate(); return result
        }
    }
}
