import Foundation

public struct VivoPCANeighborPlan: Codable, Sendable, Equatable {
    public enum InputKind: String, Codable, Sendable { case fitted, query }
    public let schemaVersion: Int
    public let inputKind: InputKind
    public let neighbors: VivoSingleCellNeighborOptions
    public let execution: VivoPCANeighborExecution
    public init(inputKind: InputKind = .fitted, neighbors: VivoSingleCellNeighborOptions = .init(), execution: VivoPCANeighborExecution = .init()) {
        schemaVersion = 1; self.inputKind = inputKind; self.neighbors = neighbors; self.execution = execution
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, inputKind, neighbors, execution }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "inputKind", "neighbors", "execution"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        inputKind = try c.decodeIfPresent(InputKind.self, forKey: .inputKind) ?? .fitted
        neighbors = try c.decodeIfPresent(VivoSingleCellNeighborOptions.self, forKey: .neighbors) ?? .init()
        execution = try c.decodeIfPresent(VivoPCANeighborExecution.self, forKey: .execution) ?? .init()
    }
    public func validate() throws {
        try neighbors.validate(); try execution.validate()
        guard schemaVersion == 1, neighbors.representation != .integrated else { throw VivoOmicsError.invalid("PCA neighbor plan schema or representation") }
    }
}
public struct VivoPCANeighborExecutionReport: Codable, Sendable, Equatable {
    public let method: String
    public let cells: Int
    public let dimensions: Int
    public let directedDistanceEvaluations: Int
    public let scalarDistanceTerms: Int
    public let scoreRecordReads: Int
    public let scoreFileBytes: Int
    public let maximumMappedBytesPerWorker: Int
    public let execution: VivoPCANeighborExecution
    public let qualification: String
}
public struct VivoPCANeighborReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let plan: VivoFingerprint
    public let graph: VivoFingerprint
    public let executionReport: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoPCANeighborBundle {
    private static func staging(_ parent: URL) throws -> URL {
        let url = parent.appendingPathComponent(".numivivo-pca-neighbors-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: url, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return url
    }
    private static func read<T: Decodable>(_ type: T.Type, root: URL, name: String, maximum: Int) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoH5ADCountStore.read(root, name, maximum: maximum))
    }
    private static func write<T: Encodable>(_ value: T, root: URL, name: String, maximum: Int) throws -> VivoFingerprint {
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximum else { throw VivoOmicsError.limit("PCA neighbor artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    private static func snapshot(_ input: URL, kind: VivoPCANeighborPlan.InputKind, to output: URL) throws {
        if kind == .fitted { try VivoH5ADPCAQuery.snapshotReference(input, to: output); return }
        try FileManager.default.createDirectory(at: output, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        try VivoH5ADPCAQuery.snapshotReference(input.appendingPathComponent("reference"), to: output.appendingPathComponent("reference"))
        for (name, limit) in [("original.h5ad", 1_073_741_824), ("plan.json", 2_097_152), ("receipt.json", 65_536),
            ("metadata.json", 536_870_912), ("quality.json", 268_435_456), ("report.json", 1_048_576), ("scores.bin", 1_024_000_000)] {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: limit)
        }
    }
    public static func publish(input: URL, plan: VivoPCANeighborPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoPCANeighborReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let source = temp.appendingPathComponent("input")
        try snapshot(input, kind: plan.inputKind, to: source)
        let inputHash: VivoFingerprint, dimensions: Int
        switch plan.inputKind {
        case .fitted:
            let receipt = try VivoH5ADPCA.verify(source, implementation: implementation)
            inputHash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
            dimensions = try read(VivoH5ADPCAModel.self, root: source, name: "model.json", maximum: 67_108_864).options.components
        case .query:
            let receipt = try VivoH5ADPCAQuery.verify(source, implementation: implementation)
            inputHash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
            dimensions = try read(VivoH5ADPCAQueryReport.self, root: source, name: "report.json", maximum: 1_048_576).components
        }
        let metadata = try read(VivoSingleCellCountMetadata.self, root: source, name: "metadata.json", maximum: 536_870_912)
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }, n = cells.count
        let graph = try VivoWindowedPCANeighbors.run(source: source.appendingPathComponent("scores.bin"), cells: cells, dimensions: dimensions, options: plan.neighbors, execution: plan.execution)
        let report = VivoPCANeighborExecutionReport(method: "exact-row-owned-dispatch-windowed-PCA-knn-v1", cells: n, dimensions: dimensions,
            directedDistanceEvaluations: n*(n-1), scalarDistanceTerms: n*(n-1)*dimensions,
            scoreRecordReads: n*dimensions*((n+plan.execution.queryBlockRows-1)/plan.execution.queryBlockRows+1), scoreFileBytes: n*dimensions*16,
            maximumMappedBytesPerWorker: VivoWindowedCountRecords.windowBytes, execution: plan.execution,
            qualification: "Independent row heaps and score readers; deterministic ordered merge, shared fuzzy graph. Exact search evaluates both directions; graph distancePairs counts unique unordered pairs. Scores use bounded tiles and 16 MiB mapping windows per worker; input reconstruction, identities and final graph remain resident. No approximate search, million-cell, Metal, embedding or biological qualification.")
        let receipt = try VivoPCANeighborReceipt(schemaVersion: 1, input: inputHash,
            plan: write(plan, root: temp, name: "plan.json", maximum: 65_536), graph: write(graph, root: temp, name: "graph.json", maximum: 536_870_912),
            executionReport: write(report, root: temp, name: "execution.json", maximum: 65_536), implementation: implementation)
        _ = try write(receipt, root: temp, name: "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    public static func verify(_ root: URL, implementation: VivoFingerprint) throws -> VivoPCANeighborReceipt {
        let bytes = try VivoH5ADCountStore.read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoPCANeighborReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == bytes else {
            throw VivoOmicsError.invalid("PCA neighbor receipt")
        }
        let plan = try read(VivoPCANeighborPlan.self, root: root, name: "plan.json", maximum: 65_536)
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(input: root.appendingPathComponent("input"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("PCA neighbor source reconstruction differs") }
        for (name, hash, maximum) in [("plan.json", receipt.plan, 65_536), ("graph.json", receipt.graph, 536_870_912), ("execution.json", receipt.executionReport, 65_536)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: maximum) == hash else { throw VivoOmicsError.invalid("PCA neighbor artifact fingerprint differs") }
        }
        return receipt
    }
}
