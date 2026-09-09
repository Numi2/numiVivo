import Foundation

public struct VivoPCAGraphEmbeddingPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let embedding: VivoSingleCellEmbeddingOptions
    public init(embedding: VivoSingleCellEmbeddingOptions = .init()) { schemaVersion = 1; self.embedding = embedding }
    private enum CodingKeys: String, CodingKey { case schemaVersion, embedding }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "embedding"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        embedding = try c.decodeIfPresent(VivoSingleCellEmbeddingOptions.self, forKey: .embedding) ?? .init()
    }
    public func validate() throws {
        try embedding.validate()
        guard schemaVersion == 1 else { throw VivoOmicsError.invalid("graph embedding schema") }
    }
}
public struct VivoPCAGraphEmbeddingExecution: Codable, Sendable, Equatable {
    public let method: String
    public let scheduleBytes: Int
    public let maximumScheduleMappedBytes: Int
    public let maximumInputMappedBytesPerReader: Int
    public let scoreTileRows: Int
    public let qualification: String
}
public struct VivoPCAGraphEmbeddingReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let plan: VivoFingerprint
    public let result: VivoFingerprint
    public let executionReport: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoPCAGraphEmbedding {
    private static func read<T: Decodable>(_ type: T.Type, _ root: URL, _ name: String, maximum: Int) throws -> T {
        try VivoCanonicalJSON.decode(type, from: VivoH5ADCountStore.read(root, name, maximum: maximum))
    }
    private static func write<T: Encodable>(_ value: T, _ root: URL, _ name: String, maximum: Int) throws -> VivoFingerprint {
        let bytes = try VivoCanonicalJSON.encode(value)
        guard bytes.count <= maximum else { throw VivoOmicsError.limit("graph embedding artifact bytes") }
        try bytes.write(to: root.appendingPathComponent(name), options: .withoutOverwriting)
        return try VivoCanonicalJSON.fingerprint(bytes)
    }
    private static func staging(_ parent: URL) throws -> URL {
        let root = parent.appendingPathComponent(".numivivo-graph-embedding-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700]); return root
    }
    public static func publish(input: URL, plan: VivoPCAGraphEmbeddingPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoPCAGraphEmbeddingReceipt {
        try plan.validate(); try VivoH5ADCountStore.requireNew(destination)
        let temp = try staging(destination.deletingLastPathComponent()); defer { try? FileManager.default.removeItem(at: temp) }
        let source = temp.appendingPathComponent("input")
        try VivoPCANeighborBundle.snapshotGraph(input, to: source)
        let parent = try VivoPCANeighborBundle.verify(source, implementation: implementation)
        let graph = try read(VivoPCAGraphStoreReport.self, source, "graph.json", maximum: 65_536)
        let metadata = try read(VivoSingleCellCountMetadata.self, source.appendingPathComponent("input"), "metadata.json", maximum: 536_870_912)
        guard graph.cells == metadata.cells.count else { throw VivoOmicsError.invalid("embedding graph identity count") }
        let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
        guard graph.dimensions >= plan.embedding.dimensions else { throw VivoOmicsError.invalid("embedding PCA dimensions") }
        let scoreReader = try VivoPCAScoreReader(source.appendingPathComponent("input/scores.bin"), rows: graph.cells, columns: graph.dimensions)
        var scores: [[Double]] = []; scores.reserveCapacity(graph.cells)
        for start in stride(from: 0, to: graph.cells, by: 2_048) {
            let end = min(graph.cells, start + 2_048), tile = try scoreReader.readRows(start..<end)
            for row in 0..<(end-start) { scores.append(Array(tile[(row*graph.dimensions)..<(row*graph.dimensions+plan.embedding.dimensions)])) }
        }
        let schedule = try VivoFileEmbeddingSchedule(edges: source.appendingPathComponent("edges.bin"), entries: graph.connectivityEntries,
            cells: graph.cells, options: plan.embedding, scratch: temp)
        let result = try VivoSingleCellEmbedding.run(cells: cells, scores: scores, schedule: schedule, options: plan.embedding)
        let execution = VivoPCAGraphEmbeddingExecution(method: "file-backed-shared-UMAP-SplitMix64-schedule-v1", scheduleBytes: schedule.fileBytes,
            maximumScheduleMappedBytes: schedule.windowBytes, maximumInputMappedBytesPerReader: VivoWindowedCountRecords.windowBytes, scoreTileRows: 2_048,
            qualification: "Edges are scanned through 16 MiB windows; mutable 32-byte edge schedules use one shared 16 MiB file window. Only the requested 2/3 PCA columns, cell identities, initial/final coordinates and result JSON remain resident. Sequential edge order, floating-point schedule increments and SplitMix64 sampling match the resident optimizer. Scratch is not a restart checkpoint. No million-cell, biological, convergence or Metal qualification.")
        try schedule.remove()
        let receipt = try VivoPCAGraphEmbeddingReceipt(schemaVersion: 1, input: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(parent)),
            plan: write(plan, temp, "plan.json", maximum: 65_536), result: write(result, temp, "result.json", maximum: 268_435_456),
            executionReport: write(execution, temp, "execution.json", maximum: 65_536), implementation: implementation)
        _ = try write(receipt, temp, "receipt.json", maximum: 65_536)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: temp, to: destination); return receipt
    }
    public static func verify(_ root: URL, implementation: VivoFingerprint) throws -> VivoPCAGraphEmbeddingReceipt {
        let bytes = try VivoH5ADCountStore.read(root, "receipt.json", maximum: 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoPCAGraphEmbeddingReceipt.self, from: bytes)
        guard receipt.schemaVersion == 1, receipt.implementation == implementation, try VivoCanonicalJSON.encode(receipt) == bytes else { throw VivoOmicsError.invalid("graph embedding receipt") }
        let plan = try read(VivoPCAGraphEmbeddingPlan.self, root, "plan.json", maximum: 65_536)
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let rebuilt = try publish(input: root.appendingPathComponent("input"), plan: plan, implementation: implementation, to: temp.appendingPathComponent("rebuilt"))
        guard rebuilt == receipt else { throw VivoOmicsError.invalid("graph embedding reconstruction differs") }
        for (name, hash, limit) in [("plan.json", receipt.plan, 65_536), ("result.json", receipt.result, 268_435_456), ("execution.json", receipt.executionReport, 65_536)] {
            guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: limit) == hash else { throw VivoOmicsError.invalid("graph embedding fingerprint differs") }
        }
        return receipt
    }
}
