import Foundation

public struct VivoAccessibilityLSIPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let components: Int
    public let maximumBasis: Int
    public let relativeResidualTolerance: Double
    public let seed: UInt64
    public init(components: Int = 30, maximumBasis: Int = 256, relativeResidualTolerance: Double = 1e-5, seed: UInt64 = 7) {
        schemaVersion = 1; self.components = components; self.maximumBasis = maximumBasis
        self.relativeResidualTolerance = relativeResidualTolerance; self.seed = seed
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, components, maximumBasis, relativeResidualTolerance, seed }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "components", "maximumBasis", "relativeResidualTolerance", "seed"])
        let values = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try values.decode(Int.self, forKey: .schemaVersion)
        components = try values.decode(Int.self, forKey: .components)
        maximumBasis = try values.decode(Int.self, forKey: .maximumBasis)
        relativeResidualTolerance = try values.decode(Double.self, forKey: .relativeResidualTolerance)
        seed = try values.decode(UInt64.self, forKey: .seed)
    }
    public func validate() throws {
        guard schemaVersion == 1 else { throw VivoOmicsError.invalid("LSI plan schema") }
        var options = VivoSingleCellReductionOptions()
        options.components = components; options.maximumBasis = maximumBasis
        options.relativeResidualTolerance = relativeResidualTolerance
        try options.validate()
    }
}

public struct VivoAccessibilityLSIReceipt: Codable, Sendable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let axes: VivoFingerprint
    public let plan: VivoFingerprint
    public let model: VivoFingerprint
    public let leftSingularVectors: VivoFingerprint
    public let embeddings: VivoFingerprint
    public let loadings: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoAccessibilityLSIIO {
    private struct Model: Codable {
        let method: String
        let matrixFormat: String
        let cells: Int
        let features: Int
        let components: Int
        let singularValues: [Double]
        let relativeResiduals: [Double]
        let maximumLoadingOrthogonalityError: Double
        let retainedEnergyFraction: Double
        let zeroRows: [Int]
        let basisSize: Int
        let operatorScans: Int
        let logDepthPearsonCorrelations: [Double?]
    }
    public static func publish(tfidf: URL, plan: VivoAccessibilityLSIPlan,
                               implementation: VivoFingerprint, to destination: URL) throws -> VivoAccessibilityLSIReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let input = try VivoAccessibilityTFIDFIO.loadVerified(tfidf)
        guard input.values.count == input.statistics.records * 16 else { throw VivoOmicsError.invalid("LSI input records") }
        let temp = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-lsi-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temp, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temp) }
        let cells = input.axes.observations.map { VivoOmicsCellIdentity(sampleID: $0.identity.sampleID, barcode: $0.identity.barcode) }
        let result = try VivoAccessibilityLSI.fit(cells: cells, featureIDs: input.axes.features.map(\.id),
            components: plan.components, maximumBasis: plan.maximumBasis,
            relativeResidualTolerance: plan.relativeResidualTolerance, seed: plan.seed,
            scan: { accept in
                try Task.checkCancellation()
                try input.values.withUnsafeBytes { bytes in
                    for offset in stride(from: 0, to: bytes.count, by: 16) {
                        let coordinate = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset, as: UInt64.self))
                        let bits = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset + 8, as: UInt64.self))
                        try accept(Int(coordinate & 0xffff_ffff), Int(coordinate >> 32), Double(bitPattern: bits))
                    }
                }
            })
        let depths = input.statistics.cellTotals.map { log1p(Double($0)) }
        let depthMean = depths.reduce(0, +) / Double(depths.count)
        let depthEnergy = depths.reduce(0) { $0 + ($1 - depthMean) * ($1 - depthMean) }
        let correlations: [Double?] = (0..<plan.components).map { k in
            let values = result.leftSingularVectors.map { $0[k] }
            let mean = values.reduce(0, +) / Double(values.count)
            let energy = values.reduce(0) { $0 + ($1 - mean) * ($1 - mean) }
            guard energy > 0, depthEnergy > 0 else { return nil }
            let covariance = values.indices.reduce(0.0) { $0 + (values[$1] - mean) * (depths[$1] - depthMean) }
            return covariance / sqrt(energy * depthEnergy)
        }
        func write<T: Encodable>(_ value: T, _ name: String) throws -> VivoFingerprint {
            let bytes = try VivoCanonicalJSON.encode(value)
            guard bytes.count <= 536_870_912 else { throw VivoOmicsError.limit("LSI metadata bytes") }
            try bytes.write(to: temp.appendingPathComponent(name), options: .withoutOverwriting)
            return try VivoCanonicalJSON.fingerprint(bytes)
        }
        try input.receiptData.write(to: temp.appendingPathComponent("input-receipt.json"), options: .withoutOverwriting)
        try input.axesData.write(to: temp.appendingPathComponent("axes.json"), options: .withoutOverwriting)
        let model = Model(method: result.method, matrixFormat: "complete-row-u32-column-u32-f64-le/v1",
            cells: cells.count, features: result.featureIDs.count, components: plan.components,
            singularValues: result.singularValues, relativeResiduals: result.relativeResiduals,
            maximumLoadingOrthogonalityError: result.maximumLoadingOrthogonalityError,
            retainedEnergyFraction: result.retainedEnergyFraction, zeroRows: result.zeroRows,
            basisSize: result.basisSize, operatorScans: result.operatorScans, logDepthPearsonCorrelations: correlations)
        let receipt = VivoAccessibilityLSIReceipt(schemaVersion: 1,
            input: try VivoCanonicalJSON.fingerprint(input.receiptData), axes: input.receipt.axes,
            plan: try write(plan, "plan.json"), model: try write(model, "model.json"),
            leftSingularVectors: try VivoH5ADPCA.writeMatrix(result.leftSingularVectors, columns: plan.components, to: temp.appendingPathComponent("left-singular-vectors.bin")),
            embeddings: try VivoH5ADPCA.writeMatrix(result.standardizedEmbeddings, columns: plan.components, to: temp.appendingPathComponent("embeddings.bin")),
            loadings: try VivoH5ADPCA.writeMatrix(result.featureLoadings, columns: plan.components, to: temp.appendingPathComponent("loadings.bin")),
            implementation: implementation)
        _ = try write(receipt, "receipt.json")
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination)
        return receipt
    }
}
