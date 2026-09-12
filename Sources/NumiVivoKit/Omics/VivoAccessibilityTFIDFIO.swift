import Foundation

public struct VivoAccessibilityTFIDFReceipt: Codable, Sendable {
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
    private struct Axes: Codable {
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
        let writer = try VivoCountRecordWriter(temp.appendingPathComponent("values.bin"))
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
        _ = try write(receipt, "receipt.json")
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }
}
