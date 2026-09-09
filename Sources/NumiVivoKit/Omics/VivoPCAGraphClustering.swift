import Foundation

public struct VivoPCAGraphClusteringPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let clustering: VivoSingleCellClusteringOptions
    public let maximumEdgeVisits: Int
    public init(clustering: VivoSingleCellClusteringOptions = .init(), maximumEdgeVisits: Int = 1_000_000_000) {
        schemaVersion = 1; self.clustering = clustering; self.maximumEdgeVisits = maximumEdgeVisits
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, clustering, maximumEdgeVisits }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "clustering", "maximumEdgeVisits"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        clustering = try c.decodeIfPresent(VivoSingleCellClusteringOptions.self, forKey: .clustering) ?? .init()
        maximumEdgeVisits = try c.decodeIfPresent(Int.self, forKey: .maximumEdgeVisits) ?? 1_000_000_000
    }
    public func validate() throws {
        try clustering.validate()
        guard schemaVersion == 1, (1...20_000_000_000).contains(maximumEdgeVisits) else { throw VivoOmicsError.invalid("graph clustering plan or work budget") }
    }
}
public struct VivoPCAGraphClusteringExecution: Codable, Sendable, Equatable {
    public let method: String
    public let edgeVisits: Int
    public let rowReads: Int
    public let aggregatedLevels: Int
    public let maximumEdgeMapEntries: Int
    public let maximumMappedBytesPerReader: Int
    public let qualification: String
}
public struct VivoPCAGraphClusteringReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let plan: VivoFingerprint
    public let result: VivoFingerprint
    public let executionReport: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoPCAGraphClustering {
    private static func read<T: Decodable>(_ type: T.Type, _ root: URL, _ name: String, maximum: Int) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoH5ADCountStore.read(root, name, maximum: maximum))
    }
    private static func write<T: Encodable>(_ value: T, _ root: URL, _ name: String, maximum: Int) throws -> VivoFingerprint {
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximum else { throw VivoOmicsError.limit("graph clustering artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    private static func staging(_ parent: URL) throws -> URL {
        let root = parent.appendingPathComponent(".numivivo-graph-clustering-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return root
    }
    private static func snapshot(_ input: URL, to output: URL) throws {
        let plan = try read(VivoPCANeighborPlan.self, input, "plan.json", maximum: 65_536)
        try plan.validate()
        guard plan.storage == .binary else { throw VivoOmicsError.invalid("file-backed clustering requires a binary graph store") }
        try FileManager.default.createDirectory(at: output, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        try VivoPCANeighborBundle.snapshot(input.appendingPathComponent("input"), kind: plan.inputKind, to: output.appendingPathComponent("input"))
        for name in ["plan.json", "graph.json", "receipt.json", "execution.json"] {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: 65_536)
        }
        for name in VivoPCAGraphStore.files {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: 4_096_000_000)
        }
    }
    public static func publish(input: URL, plan: VivoPCAGraphClusteringPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoPCAGraphClusteringReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let source = temp.appendingPathComponent("input")
        try snapshot(input, to: source)
        let parent = try VivoPCANeighborBundle.verify(source, implementation: implementation)
        let graph = try read(VivoPCAGraphStoreReport.self, source, "graph.json", maximum: 65_536)
        let metadata = try read(VivoSingleCellCountMetadata.self, source.appendingPathComponent("input"), "metadata.json", maximum: 536_870_912)
        guard graph.cells == metadata.cells.count else { throw VivoOmicsError.invalid("clustering graph identity count") }
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        let work = VivoClusteringWork(maximumEdgeVisits: plan.maximumEdgeVisits)
        let reader = try VivoFileClusteringGraph(root: source, rows: graph.cells, entries: graph.connectivityEntries, work: work)
        let result = try VivoSingleCellClustering.run(original: reader, cells: cells, options: plan.clustering, scratch: temp)
        let execution = VivoPCAGraphClusteringExecution(method: "file-backed-shared-Louvain-stable-scatter-v1", edgeVisits: work.edgeVisits,
            rowReads: work.rowReads, aggregatedLevels: work.aggregatedLevels, maximumEdgeMapEntries: work.maximumEdgeMapEntries,
            maximumMappedBytesPerReader: VivoWindowedCountRecords.windowBytes,
            qualification: "Original and aggregated CSR edges use 16 MiB file windows. Stable disk scatter preserves resident source-row/column summation order; each aggregate row uses one edge map. edgeVisits includes repeated graph reads and temporary aggregation records. Cell identities, offsets, labels, degrees, totals, permutation, traversal queues and output JSON remain resident. No million-cell, Leiden, biological or Metal qualification.")
        let receipt = try VivoPCAGraphClusteringReceipt(schemaVersion: 1, input: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(parent)),
            plan: write(plan, temp, "plan.json", maximum: 65_536), result: write(result, temp, "result.json", maximum: 268_435_456),
            executionReport: write(execution, temp, "execution.json", maximum: 65_536), implementation: implementation)
        _ = try write(receipt, temp, "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    public static func verify(_ root: URL, implementation: VivoFingerprint) throws -> VivoPCAGraphClusteringReceipt {
        let bytes = try VivoH5ADCountStore.read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoPCAGraphClusteringReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("graph clustering receipt") }
        let plan = try read(VivoPCAGraphClusteringPlan.self, root, "plan.json", maximum: 65_536)
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(input: root.appendingPathComponent("input"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("graph clustering reconstruction differs") }
        for (name, hash, limit) in [("plan.json", receipt.plan, 65_536), ("result.json", receipt.result, 268_435_456), ("execution.json", receipt.executionReport, 65_536)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: limit) == hash else { throw VivoOmicsError.invalid("graph clustering fingerprint differs") }
        }
        return receipt
    }
}
