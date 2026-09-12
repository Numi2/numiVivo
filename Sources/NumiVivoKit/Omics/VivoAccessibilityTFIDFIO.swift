import Foundation

public struct VivoAccessibilityTFIDFReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let format: String
    public let assayID: String
    public let source: VivoFingerprint
    public let mapping: VivoFingerprint
    public let axes: VivoFingerprint
    public let statistics: VivoFingerprint
    public let values: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoAccessibilityTFIDFIO {
    struct VerifiedInput {
        let receipt: VivoAccessibilityTFIDFReceipt
        let receiptData: Data
        let axes: Axes
        let axesData: Data
        let statistics: VivoAccessibilityTFIDFReport
        let values: Data
    }
    public struct Verification: Codable, Sendable {
        public let receipt: VivoAccessibilityTFIDFReceipt
        public let verifier: VivoFingerprint
    }
    /// Reconstructs source counts and normalization, including across producer
    /// revisions. The recorded producer is provenance, not authentication.
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> Verification {
        let input = try loadVerified(directory)
        return .init(receipt: input.receipt, verifier: implementation)
    }
    static func loadVerified(_ directory: URL) throws -> VerifiedInput {
        func read(_ name: String, _ maximum: Int) throws -> Data {
            try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name), maximumBytes: maximum)
        }
        let receiptData = try read("receipt.json", 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoAccessibilityTFIDFReceipt.self, from: receiptData)
        guard receipt.schemaVersion == 1, receipt.method == "signac-method1-scale10000-count-idf-v1",
              receipt.format == "sparse-row-u32-feature-u32-value-f64-le/v1" else { throw VivoOmicsError.invalid("TF-IDF receipt schema or method") }
        let planData = try read("mapping.json", 131_072), axesData = try read("axes.json", 536_870_912)
        let statisticsData = try read("statistics.json", 536_870_912)
        let temp = FileManager.default.temporaryDirectory.appendingPathComponent(".numivivo-tfidf-verify-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temp, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temp) }
        // Reuse the native immutable clone/copy boundary. Mapped Data retains
        // its backing inode after the private snapshot directory is removed.
        let valuesURL = temp.appendingPathComponent("values.bin")
        guard try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent("values.bin"), copyTo: valuesURL,
            maximumBytes: 536_870_912) == receipt.values else { throw VivoOmicsError.invalid("TF-IDF value hash") }
        let values = try Data(contentsOf: valuesURL, options: .mappedIfSafe)
        guard try VivoCanonicalJSON.fingerprint(planData) == receipt.mapping,
              try VivoCanonicalJSON.fingerprint(axesData) == receipt.axes,
              try VivoCanonicalJSON.fingerprint(statisticsData) == receipt.statistics,
              try VivoCanonicalJSON.fingerprint(values) == receipt.values else { throw VivoOmicsError.invalid("TF-IDF artifact hash") }
        let plan = try VivoCanonicalJSON.decode(VivoTenXMultiAssayPlan.self, from: planData)
        let rebuilt = try build(source: directory.appendingPathComponent("original.h5"), plan: plan,
            assayID: receipt.assayID, implementation: receipt.implementation, temp: temp, retainValues: false)
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("TF-IDF differs from source reconstruction") }
        return .init(receipt: receipt, receiptData: receiptData,
            axes: try VivoCanonicalJSON.decode(Axes.self, from: axesData), axesData: axesData,
            statistics: try VivoCanonicalJSON.decode(VivoAccessibilityTFIDFReport.self, from: statisticsData), values: values)
    }
    struct Axes: Codable {
        let assayID: String
        let kind: VivoAssayKind
        let featureNamespace: String
        let countUnit: VivoAssayCountUnit
        let genomeAssembly: String?
        let features: [VivoAssayFeature]
        let sourceObservationIndices: [Int]
        let observations: [VivoAssayObservation]
        let samples: [VivoOmicsSample]
    }
    public static func publishTenX(source: URL, plan: VivoTenXMultiAssayPlan,
                                    assayID: String, implementation: VivoFingerprint,
                                    to destination: URL) throws -> VivoAccessibilityTFIDFReceipt {
        try VivoMultiAssayTenX.validate(plan)
        try VivoH5ADCountStore.requireNew(destination)
        guard plan.assays.contains(where: { $0.id == assayID && $0.kind == .chromatinAccessibility }) else {
            throw VivoOmicsError.invalid("TF-IDF requires an explicitly mapped accessibility assay")
        }
        let temp = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-tfidf-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temp, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temp) }
        let receipt = try build(source: source, plan: plan, assayID: assayID, implementation: implementation, temp: temp, retainValues: true)
        try VivoCanonicalJSON.encode(receipt).write(to: temp.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }

    private static func build(source: URL, plan: VivoTenXMultiAssayPlan, assayID: String,
                              implementation: VivoFingerprint, temp: URL, retainValues: Bool) throws -> VivoAccessibilityTFIDFReceipt {
        try VivoMultiAssayTenX.validate(plan)
        let snapshot = temp.appendingPathComponent("original.h5")
        let sourceHash = try VivoH5ADPseudobulk.fingerprint(source, copyTo: snapshot)
        let dataset = try VivoMultiAssayTenX.readSnapshot(snapshot, plan: plan)
        guard let assay = dataset.assays.first(where: { $0.id == assayID && $0.kind == .chromatinAccessibility }) else {
            throw VivoOmicsError.invalid("TF-IDF accessibility assay missing")
        }
        func write<T: Encodable>(_ object: T, _ name: String) throws -> VivoFingerprint {
            let data = try VivoCanonicalJSON.encode(object)
            guard data.count <= 536_870_912 else { throw VivoOmicsError.limit("TF-IDF metadata bytes") }
            try data.write(to: temp.appendingPathComponent(name), options: .withoutOverwriting)
            return try VivoCanonicalJSON.fingerprint(data)
        }
        let mapping = try write(plan, "mapping.json")
        let axes = try write(Axes(assayID: assay.id, kind: assay.kind, featureNamespace: assay.featureNamespace,
            countUnit: assay.countUnit, genomeAssembly: assay.genomeAssembly, features: assay.features,
            sourceObservationIndices: assay.observationIndices,
            observations: assay.observationIndices.map { dataset.observations[$0] }, samples: dataset.samples), "axes.json")
        let writer = try VivoCountRecordWriter(retainValues ? temp.appendingPathComponent("values.bin") : nil)
        let matrix = assay.matrix
        let report = try VivoAccessibilityTFIDF.normalize(cells: assay.observationIndices.count,
            features: assay.features.count, maximumRecords: matrix.counts.count,
            scan: { accept in
                for row in assay.observationIndices.indices {
                    try Task.checkCancellation()
                    for k in matrix.rowOffsets[row]..<matrix.rowOffsets[row + 1] {
                        try accept(row, matrix.featureIndices[k], matrix.counts[k])
                    }
                }
            }, emit: { row, feature, value in try writer.append(row: row, feature: feature, bits: value.bitPattern) })
        let receipt = VivoAccessibilityTFIDFReceipt(schemaVersion: 1, method: "signac-method1-scale10000-count-idf-v1",
            format: "sparse-row-u32-feature-u32-value-f64-le/v1", assayID: assayID, source: sourceHash,
            mapping: mapping, axes: axes, statistics: try write(report, "statistics.json"),
            values: try writer.finish(), implementation: implementation)
        return receipt
    }
}
