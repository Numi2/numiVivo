import Foundation

public struct VivoH5ADPCAPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let mapping: VivoH5ADImportPlan
    public var reduction: VivoH5ADReductionOptions
    public let featureNamespace: String?
    public init(mapping: VivoH5ADImportPlan, reduction: VivoH5ADReductionOptions = .init(), featureNamespace: String? = nil) {
        schemaVersion = 1; self.mapping = mapping; self.reduction = reduction; self.featureNamespace = featureNamespace
        if self.reduction.pca.retainProjectionCenters == nil { self.reduction.pca.retainProjectionCenters = true }
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, mapping, reduction, featureNamespace }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "mapping", "reduction", "featureNamespace"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        mapping = try c.decode(VivoH5ADImportPlan.self, forKey: .mapping)
        featureNamespace = try c.decodeIfPresent(String.self, forKey: .featureNamespace)
        reduction = try c.decodeIfPresent(VivoH5ADReductionOptions.self, forKey: .reduction) ?? .init()
        if reduction.pca.retainProjectionCenters == nil { reduction.pca.retainProjectionCenters = true }
    }
    public func validate() throws {
        guard schemaVersion == 1, reduction.pca.retainProjectionCenters == true,
              featureNamespace.map(vivoOmicsID) ?? true else {
            throw VivoOmicsError.invalid("standalone PCA schema or required projection centers")
        }
        try reduction.validate()
    }
}

/// Small model state; large cell scores and selected-feature loadings live in
/// separately fingerprinted row-major binary records. Cell/feature axes are
/// bound by metadata and selectedFeatureIndices, never inferred by a consumer.
public struct VivoH5ADPCAModel: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let options: VivoSingleCellReductionOptions
    public let normalizationTarget: Double
    public let cells: Int
    public let canonicalNonzeros: Int
    public let hdf5Version: String
    public let features: [VivoSingleCellVariableFeature]
    public let selectedFeatureIndices: [Int]
    public let explainedVariance: [Double]
    public let explainedVarianceRatio: [Double]
    public let relativeResiduals: [Double]
    public let maximumLoadingOrthogonalityError: Double
    public let basisSize: Int
    public let projectionCenters: [Double]
    public let storage: VivoH5ADReductionStorage
    public let qualification: String
}
public struct VivoH5ADPCAReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let matrixFormat: String
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let metadata: VivoFingerprint
    public let quality: VivoFingerprint
    public let model: VivoFingerprint
    public let scores: VivoFingerprint
    public let loadings: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoH5ADPCA {
    static let matrixFormat = "complete-row-major-u32-row-u32-component-f64-le/v1"
    private static func staging(_ parent: URL) throws -> URL {
        let url = parent.appendingPathComponent(".numivivo-pca-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        return url
    }
    /// Every coordinate, including zeros, is emitted in row-major order. The
    /// record primitive buffers 1 MiB, and never encodes the score matrix as JSON.
    static func writeMatrix(_ values: [[Double]], columns: Int, to url: URL) throws -> VivoFingerprint {
        guard columns > 0, columns <= 64, values.count <= 1_000_000 else { throw VivoOmicsError.limit("PCA matrix dimensions") }
        let writer = try VivoCountRecordWriter(url)
        for (row, vector) in values.enumerated() {
            guard vector.count == columns else { throw VivoOmicsError.invalid("PCA matrix row width") }
            for (column, value) in vector.enumerated() {
                guard value.isFinite else { throw VivoOmicsError.invalid("PCA matrix nonfinite value") }
                try writer.append(row: row, feature: column, bits: value.bitPattern)
            }
        }
        return try writer.finish()
    }
    private static func writeJSON<T: Encodable>(_ value: T, name: String, root: URL, maximum: Int) throws -> VivoFingerprint {
        let data = try VivoCanonicalJSON.encode(value)
        guard data.count <= maximum else { throw VivoOmicsError.limit("PCA artifact bytes") }
        try data.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(data)
    }
    public static func publish(source: URL, plan: VivoH5ADPCAPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoH5ADPCAReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let planHash = try writeJSON(plan, name: "plan.json", root: temp, maximum: 2_097_152)
        let snapshot = temp.appendingPathComponent("original.h5ad")
        let sourceHash = try VivoOmicsFileSnapshot.fingerprint(source, copyTo: snapshot, maximumBytes: VivoH5ADPseudobulk.sourceLimits.maximumInputBytes)
        var accumulator: VivoSingleCellQualityAccumulator?
        let version = try VivoSingleCellH5AD.scanSnapshot(snapshot, plan: plan.mapping, limits: VivoH5ADPseudobulk.sourceLimits, onMetadata: {
            accumulator = VivoSingleCellQualityAccumulator($0)
        }, onEntry: { row, feature, count in
            guard let accumulator else { throw VivoOmicsError.invalid("PCA source metadata missing") }
            try accumulator.add(row: row, feature: feature, count: count)
        })
        guard let accumulator else { throw VivoOmicsError.invalid("PCA source metadata missing") }
        let metadata = accumulator.metadata, quality = accumulator.finish()
        let (result, storage) = try VivoH5ADReduction.run(snapshot: snapshot, mapping: plan.mapping, metadata: metadata, quality: quality, options: plan.reduction)
        guard let centers = result.projectionCenters, centers.count == result.selectedFeatureIndices.count else { throw VivoOmicsError.invalid("PCA centers missing") }
        let model = VivoH5ADPCAModel(schemaVersion: 1, method: result.method, options: result.options,
            normalizationTarget: plan.reduction.normalizationTarget, cells: result.cells.count, canonicalNonzeros: accumulator.nonzeros,
            hdf5Version: version, features: result.features, selectedFeatureIndices: result.selectedFeatureIndices,
            explainedVariance: result.explainedVariance, explainedVarianceRatio: result.explainedVarianceRatio,
            relativeResiduals: result.relativeResiduals, maximumLoadingOrthogonalityError: result.maximumLoadingOrthogonalityError,
            basisSize: result.basisSize, projectionCenters: centers, storage: storage,
            qualification: "Shared native QC/HVG/PCA without condition aggregates. Binary scores/loadings; metadata, quality, feature statistics and fitted arrays remain resident. No million-cell, parallel, GPU or held-out biological qualification.")
        let metadataHash = try writeJSON(metadata, name: "metadata.json", root: temp, maximum: 536_870_912)
        let qualityHash = try writeJSON(quality, name: "quality.json", root: temp, maximum: 268_435_456)
        let modelHash = try writeJSON(model, name: "model.json", root: temp, maximum: 67_108_864)
        let scoresHash = try writeMatrix(result.scores, columns: result.options.components, to: temp.appendingPathComponent("scores.bin"))
        let loadingsHash = try writeMatrix(result.loadings, columns: result.options.components, to: temp.appendingPathComponent("loadings.bin"))
        let receipt = VivoH5ADPCAReceipt(schemaVersion: 1, matrixFormat: matrixFormat, source: sourceHash, plan: planHash,
            metadata: metadataHash, quality: qualityHash, model: modelHash, scores: scoresHash, loadings: loadingsHash, implementation: implementation)
        _ = try writeJSON(receipt, name: "receipt.json", root: temp, maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> VivoH5ADPCAReceipt {
        let bytes = try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("receipt.json"), maximumBytes: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoH5ADPCAReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.matrixFormat == matrixFormat, receipt.implementation == implementation,
              try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("PCA receipt") }
        let planBytes = try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent("plan.json"), maximumBytes: 2_097_152)
        let plan = try VivoCanonicalJSON.decode(VivoH5ADPCAPlan.self, from: planBytes)
        guard try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan, try VivoCanonicalJSON.encode(plan) == planBytes else { throw VivoOmicsError.invalid("PCA plan") }
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(source: directory.appendingPathComponent("original.h5ad"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("PCA source reconstruction differs") }
        for (name, hash, maximum) in [("metadata.json", receipt.metadata, 536_870_912), ("quality.json", receipt.quality, 268_435_456),
            ("model.json", receipt.model, 67_108_864), ("scores.bin", receipt.scores, 1_024_000_000), ("loadings.bin", receipt.loadings, 10_240_000)] {
            guard try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent(name), maximumBytes: maximum) == hash else {
                throw VivoOmicsError.invalid("PCA artifact fingerprint differs")
            }
        }
        return receipt
    }
}
