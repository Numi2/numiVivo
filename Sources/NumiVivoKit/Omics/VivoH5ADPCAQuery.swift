import Foundation

public struct VivoH5ADPCAQueryPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public let featureNamespace: String
    public let maximumProjectionUpdates: Int
    public init(mapping: VivoH5ADImportPlan, featureNamespace: String, maximumProjectionUpdates: Int = 2_000_000_000) {
        schemaVersion = 1; self.mapping = mapping; self.featureNamespace = featureNamespace
        self.maximumProjectionUpdates = maximumProjectionUpdates
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, mapping, featureNamespace, maximumProjectionUpdates }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "mapping", "featureNamespace", "maximumProjectionUpdates"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        mapping = try c.decode(VivoH5ADImportPlan.self, forKey: .mapping)
        featureNamespace = try c.decode(String.self, forKey: .featureNamespace)
        maximumProjectionUpdates = try c.decodeIfPresent(Int.self, forKey: .maximumProjectionUpdates) ?? 2_000_000_000
    }
    public func validate() throws {
        guard schemaVersion == 1, vivoOmicsID(featureNamespace), mapping.groupColumn == nil,
              (1...20_000_000_000).contains(maximumProjectionUpdates) else {
            throw VivoOmicsError.invalid("PCA query schema, namespace, label mapping or work budget")
        }
    }
}
public struct VivoH5ADPCAQueryReport: Codable, Sendable, Equatable {
    public let method: String
    public let cells: Int
    public let components: Int
    public let selectedEntries: Int
    public let projectionUpdates: Int
    public let scratchBytes: Int
    public let sourcePasses: Int
    public let emptyLibraries: Int
    public let overlappingDonorIDs: [String]
    public let hdf5Version: String
    public let qualification: String
}
public struct VivoH5ADPCAQueryReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let matrixFormat: String
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let reference: VivoFingerprint
    public let metadata: VivoFingerprint
    public let quality: VivoFingerprint
    public let report: VivoFingerprint
    public let scores: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoH5ADPCAQuery {
    private static func staging(_ parent: URL) throws -> URL {
        let url = parent.appendingPathComponent(".numivivo-pca-query-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return url
    }
    private static func read<T: Decodable>(_ type: T.Type, root: URL, name: String, maximum: Int) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoH5ADCountStore.read(root, name, maximum: maximum))
    }
    private static func write<T: Encodable>(_ value: T, root: URL, name: String, maximum: Int) throws -> VivoFingerprint {
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximum else { throw VivoOmicsError.limit("PCA query artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    private static func snapshotReference(_ source: URL, to destination: URL) throws {
        try FileManager.default.createDirectory(at: destination, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        for (name, limit) in [("original.h5ad", 1_073_741_824), ("plan.json", 2_097_152), ("receipt.json", 65_536),
            ("metadata.json", 536_870_912), ("quality.json", 268_435_456), ("model.json", 67_108_864),
            ("scores.bin", 1_024_000_000), ("loadings.bin", 10_240_000)] {
            _ = try VivoOmicsFileSnapshot.fingerprint(source.appendingPathComponent(name), copyTo: destination.appendingPathComponent(name), maximumBytes: limit)
        }
    }
    public static func publish(source: URL, plan: VivoH5ADPCAQueryPlan, reference: URL, implementation: VivoFingerprint, to destination: URL) throws -> VivoH5ADPCAQueryReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let referenceCopy = temp.appendingPathComponent("reference")
        try snapshotReference(reference, to: referenceCopy)
        let referenceReceipt = try VivoH5ADPCA.verify(referenceCopy, implementation: implementation)
        let trainingPlan = try read(VivoH5ADPCAPlan.self, root: referenceCopy, name: "plan.json", maximum: 2_097_152)
        guard trainingPlan.featureNamespace == plan.featureNamespace, trainingPlan.mapping.countUnit == plan.mapping.countUnit else {
            throw VivoOmicsError.invalid("PCA query feature namespace or count unit mismatch; training namespace required")
        }
        let model = try read(VivoH5ADPCAModel.self, root: referenceCopy, name: "model.json", maximum: 67_108_864)
        let training = try read(VivoSingleCellCountMetadata.self, root: referenceCopy, name: "metadata.json", maximum: 536_870_912)
        let d = model.options.components, width = model.selectedFeatureIndices.count
        let records = try VivoWindowedCountRecords(referenceCopy.appendingPathComponent("loadings.bin"), entries: width * d)
        var loadings = Array(repeating: Array(repeating: 0.0, count: d), count: width)
        for i in 0..<records.count {
            let record = try records.record(i), value = Double(bitPattern: record.bits)
            guard record.row == i / d, record.feature == i % d, value.isFinite else { throw VivoOmicsError.invalid("PCA loading coordinates or values") }
            loadings[record.row][record.feature] = value
        }
        let projection = try VivoFrozenPCAProjection(centers: model.projectionCenters, loadings: loadings, target: model.normalizationTarget)
        let sourceHash = try VivoOmicsFileSnapshot.fingerprint(source, copyTo: temp.appendingPathComponent("original.h5ad"), maximumBytes: VivoH5ADPseudobulk.sourceLimits.maximumInputBytes)
        let snapshot = temp.appendingPathComponent("original.h5ad")
        let trainingIDs = Set(training.features.map(\.id)), organisms = Set(training.samples.map(\.organism))
        let trainingCells = Set(training.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) })
        let selected = Dictionary(uniqueKeysWithValues: model.selectedFeatureIndices.enumerated().map { (training.features[$0.element].id, $0.offset) })
        var qc: VivoSingleCellQualityAccumulator?, local: [Int] = [], selectedEntries = 0
        let version = try VivoSingleCellH5AD.scanSnapshot(snapshot, plan: plan.mapping, limits: VivoH5ADPseudobulk.sourceLimits, onMetadata: { metadata in
            guard organisms.count == 1, Set(metadata.samples.map(\.organism)) == organisms,
                  metadata.features.count == training.features.count, Set(metadata.features.map(\.id)) == trainingIDs else {
                throw VivoOmicsError.invalid("PCA query must match complete training gene universe and organism")
            }
            guard metadata.cells.allSatisfy({ !trainingCells.contains(.init(sampleID: $0.sampleID, barcode: $0.barcode)) }) else {
                throw VivoOmicsError.invalid("PCA query overlaps training cell identities")
            }
            qc = VivoSingleCellQualityAccumulator(metadata); local = metadata.features.map { selected[$0.id] ?? -1 }
        }, onEntry: { row, column, count in
            guard let qc else { throw VivoOmicsError.invalid("PCA query metadata missing") }
            try qc.add(row: row, feature: column, count: count)
            if local[column] >= 0 {
                guard selectedEntries < plan.maximumProjectionUpdates / d else { throw VivoOmicsError.limit("PCA query projection-update budget") }
                selectedEntries += 1
            }
        })
        guard let qc else { throw VivoOmicsError.invalid("PCA query metadata missing") }
        let quality = qc.finish(), metadata = qc.metadata
        let scratchURL = temp.appendingPathComponent(".scores-scratch")
        let scores = try VivoWindowedPCAScores(scratchURL, rows: metadata.cells.count, initial: projection.initial)
        var written = 0
        _ = try VivoSingleCellH5AD.scanSnapshot(snapshot, plan: plan.mapping, limits: VivoH5ADPseudobulk.sourceLimits, onMetadata: {
            guard $0 == metadata else { throw VivoOmicsError.invalid("PCA query metadata changed") }
        }, onEntry: { row, column, count in
            let index = local[column]; if index < 0 { return }
            guard written < selectedEntries else { throw VivoOmicsError.invalid("PCA query entry count changed") }
            try scores.withRow(row) { values in
                try projection.add(count: count, total: quality[row].totalCounts, selected: index, to: values)
            }; written += 1
        })
        guard written == selectedEntries else { throw VivoOmicsError.invalid("PCA query entry count changed") }
        let scoreHash = try scores.write(to: temp.appendingPathComponent("scores.bin"))
        try FileManager.default.removeItem(at: scratchURL)
        let report = VivoH5ADPCAQueryReport(method: "frozen-training-log-PCA-windowed-query-v1", cells: metadata.cells.count, components: d,
            selectedEntries: selectedEntries, projectionUpdates: selectedEntries * d, scratchBytes: scores.byteCount, sourcePasses: 2,
            emptyLibraries: quality.filter { $0.totalCounts == 0 }.count,
            overlappingDonorIDs: Set(training.samples.compactMap(\.donorID)).intersection(metadata.samples.compactMap(\.donorID)).sorted(),
            hdf5Version: version,
            qualification: "Training-only features, centers and loadings; exact full gene universe, namespace, organism and count units required. Query cells must be disjoint by sampleID/barcode. Empty libraries retain the mathematical centered-zero score and are flagged in QC. No labels, biological generalization or million-cell qualification. Query scores use 16 MiB mapping windows; metadata/QC and reference reconstruction remain resident.")
        let receipt = try VivoH5ADPCAQueryReceipt(schemaVersion: 1, matrixFormat: VivoH5ADPCA.matrixFormat, source: sourceHash,
            plan: write(plan, root: temp, name: "plan.json", maximum: 2_097_152), reference: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(referenceReceipt)),
            metadata: write(metadata, root: temp, name: "metadata.json", maximum: 536_870_912), quality: write(quality, root: temp, name: "quality.json", maximum: 268_435_456),
            report: write(report, root: temp, name: "report.json", maximum: 1_048_576), scores: scoreHash, implementation: implementation)
        _ = try write(receipt, root: temp, name: "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    public static func verify(_ root: URL, implementation: VivoFingerprint) throws -> VivoH5ADPCAQueryReceipt {
        let bytes = try VivoH5ADCountStore.read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoH5ADPCAQueryReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.matrixFormat == VivoH5ADPCA.matrixFormat, receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("PCA query receipt") }
        let plan = try read(VivoH5ADPCAQueryPlan.self, root: root, name: "plan.json", maximum: 2_097_152)
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(source: root.appendingPathComponent("original.h5ad"), plan: plan, reference: root.appendingPathComponent("reference"), implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("PCA query source reconstruction differs") }
        for (name, hash, maximum) in [("plan.json", receipt.plan, 2_097_152), ("metadata.json", receipt.metadata, 536_870_912),
            ("quality.json", receipt.quality, 268_435_456), ("report.json", receipt.report, 1_048_576), ("scores.bin", receipt.scores, 1_024_000_000)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: maximum) == hash else {
                throw VivoOmicsError.invalid("PCA query artifact fingerprint differs")
            }
        }
        return receipt
    }
}
