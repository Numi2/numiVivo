import Foundation

public struct VivoPCANeighborPlan: Codable, Sendable, Equatable {
    public enum InputKind: String, Codable, Sendable { case fitted, query, integrated }
    public enum Storage: String, Codable, Sendable { case json, binary }
    public let schemaVersion: Int
    public let inputKind: InputKind
    public let neighbors: VivoSingleCellNeighborOptions
    public let execution: VivoPCANeighborExecution
    public let approximation: VivoHNSWOptions?
    public let storage: Storage?
    public init(inputKind: InputKind = .fitted, neighbors: VivoSingleCellNeighborOptions = .init(), execution: VivoPCANeighborExecution = .init(), approximation: VivoHNSWOptions? = nil, storage: Storage? = nil) {
        schemaVersion = 1; self.inputKind = inputKind; self.neighbors = neighbors; self.execution = execution; self.approximation = approximation; self.storage = storage
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, inputKind, neighbors, execution, approximation, storage }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "inputKind", "neighbors", "execution", "approximation", "storage"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        inputKind = try c.decodeIfPresent(InputKind.self, forKey: .inputKind) ?? .fitted
        approximation = try c.decodeIfPresent(VivoHNSWOptions.self, forKey: .approximation)
        storage = try c.decodeIfPresent(Storage.self, forKey: .storage)
        neighbors = try c.decodeIfPresent(VivoSingleCellNeighborOptions.self, forKey: .neighbors) ?? .init()
        execution = try c.decodeIfPresent(VivoPCANeighborExecution.self, forKey: .execution) ?? .init()
    }
    public func validate() throws {
        try neighbors.validate(); try execution.validate()
        try approximation?.validate(neighbors: neighbors.neighbors)
        guard approximation == nil || execution.workers == 1 else { throw VivoOmicsError.invalid("HNSW currently requires one serial worker") }
        guard schemaVersion == 1, (inputKind == .integrated) == (neighbors.representation == .integrated) else { throw VivoOmicsError.invalid("PCA neighbor plan schema or representation") }
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
    public var hnsw: VivoHNSWReport? = nil
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
    static func snapshot(_ input: URL, kind: VivoPCANeighborPlan.InputKind, to output: URL) throws {
        if kind == .integrated { try VivoPCAIntegration.snapshot(input, to: output); return }
        if kind == .fitted { try VivoH5ADPCAQuery.snapshotReference(input, to: output); return }
        try FileManager.default.createDirectory(at: output, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        try VivoH5ADPCAQuery.snapshotReference(input.appendingPathComponent("reference"), to: output.appendingPathComponent("reference"))
        for (name, limit) in [("original.h5ad", 1_073_741_824), ("plan.json", 2_097_152), ("receipt.json", 65_536),
            ("metadata.json", 536_870_912), ("quality.json", 268_435_456), ("report.json", 1_048_576), ("scores.bin", 1_024_000_000)] {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: limit)
        }
    }
    static func snapshotGraph(_ input: URL, to output: URL) throws {
        let plan = try read(VivoPCANeighborPlan.self, root: input, name: "plan.json", maximum: 65_536)
        try plan.validate()
        guard plan.storage == .binary else { throw VivoOmicsError.invalid("file-backed graph analysis requires a binary graph store") }
        try FileManager.default.createDirectory(at: output, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        try VivoPCANeighborBundle.snapshot(input.appendingPathComponent("input"), kind: plan.inputKind, to: output.appendingPathComponent("input"))
        for name in ["plan.json", "graph.json", "receipt.json", "execution.json"] {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: 65_536)
        }
        for name in VivoPCAGraphStore.files {
            _ = try VivoOmicsFileSnapshot.fingerprint(input.appendingPathComponent(name), copyTo: output.appendingPathComponent(name), maximumBytes: 4_096_000_000)
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
        case .integrated:
            let receipt = try VivoPCAIntegration.verify(source, implementation: implementation)
            inputHash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
            dimensions = try read(VivoPCAIntegrationReport.self, root: source, name: "report.json", maximum: 16_777_216).components
        case .query:
            let receipt = try VivoH5ADPCAQuery.verify(source, implementation: implementation)
            inputHash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(receipt))
            dimensions = try read(VivoH5ADPCAQueryReport.self, root: source, name: "report.json", maximum: 1_048_576).components
        }
        let metadata = try read(VivoSingleCellCountMetadata.self, root: source, name: "metadata.json", maximum: 536_870_912)
        let n = metadata.cells.count
        let graphHash: VivoFingerprint
        var report: VivoPCANeighborExecutionReport
        if plan.storage == .binary {
            let writer = try VivoCountRecordWriter(temp.appendingPathComponent("neighbors.bin"))
            let sink: (Int, [Int], [Double]) throws -> Void = { row, indices, distances in
                guard indices.count == plan.neighbors.neighbors, distances.count == indices.count else { throw VivoOmicsError.invalid("neighbor stream width") }
                for slot in indices.indices { try writer.append(row: row, feature: indices[slot], bits: distances[slot].bitPattern) }
            }
            let hnsw: VivoHNSWReport?
            if let approximation = plan.approximation {
                hnsw = try VivoHNSWNeighbors.stream(source: source.appendingPathComponent("scores.bin"), rows: n, dimensions: dimensions,
                    options: plan.neighbors, approximation: approximation, sink: sink)
            } else {
                try VivoWindowedPCANeighbors.stream(source: source.appendingPathComponent("scores.bin"), rows: n, dimensions: dimensions,
                    options: plan.neighbors, execution: plan.execution, sink: sink)
                hnsw = nil
            }
            let hash = try writer.finish()
            let work = hnsw.map { $0.constructionDistances + $0.queryDistances }
            let stored = try VivoPCAGraphStore.build(root: temp, rows: n, dimensions: dimensions, options: plan.neighbors,
                distancePairs: work ?? n*(n-1)/2, approximate: hnsw != nil, neighborsHash: hash)
            report = .init(method: hnsw == nil ? "exact-row-owned-dispatch-windowed-PCA-knn-v1" : "serial-HNSW-cached-FP64-PCA-knn-v1",
                cells: n, dimensions: dimensions, directedDistanceEvaluations: work ?? n*(n-1),
                scalarDistanceTerms: (work ?? n*(n-1))*dimensions,
                scoreRecordReads: hnsw.map { $0.scoreReadBytes/16 } ?? n*dimensions*((n+plan.execution.queryBlockRows-1)/plan.execution.queryBlockRows+1),
                scoreFileBytes: n*dimensions*16, maximumMappedBytesPerWorker: hnsw == nil ? VivoWindowedCountRecords.windowBytes : 0,
                execution: plan.execution, qualification: stored.qualification)
            report.hnsw = hnsw
            graphHash = try write(stored, root: temp, name: "graph.json", maximum: 65_536)
        } else {
            let cells = metadata.cells.map { VivoOmicsCellIdentity(sampleID: $0.sampleID, barcode: $0.barcode) }
            let graph: VivoSingleCellNeighborGraph
            if let approximation = plan.approximation {
                let result = try VivoHNSWNeighbors.run(source: source.appendingPathComponent("scores.bin"), cells: cells, dimensions: dimensions, options: plan.neighbors, approximation: approximation)
                graph = result.0
                let evaluations = result.1.constructionDistances + result.1.queryDistances
                report = .init(method: "serial-HNSW-cached-FP64-PCA-knn-v1", cells: n, dimensions: dimensions,
                    directedDistanceEvaluations: evaluations, scalarDistanceTerms: evaluations*dimensions,
                    scoreRecordReads: result.1.scoreReadBytes/16, scoreFileBytes: n*dimensions*16, maximumMappedBytesPerWorker: 0,
                    execution: plan.execution, qualification: result.1.qualification)
                report.hnsw = result.1
            } else {
                graph = try VivoWindowedPCANeighbors.run(source: source.appendingPathComponent("scores.bin"), cells: cells, dimensions: dimensions, options: plan.neighbors, execution: plan.execution)
                report = VivoPCANeighborExecutionReport(method: "exact-row-owned-dispatch-windowed-PCA-knn-v1", cells: n, dimensions: dimensions,
                directedDistanceEvaluations: n*(n-1), scalarDistanceTerms: n*(n-1)*dimensions,
                scoreRecordReads: n*dimensions*((n+plan.execution.queryBlockRows-1)/plan.execution.queryBlockRows+1), scoreFileBytes: n*dimensions*16,
                maximumMappedBytesPerWorker: VivoWindowedCountRecords.windowBytes, execution: plan.execution,
                qualification: "Independent row heaps and score readers; deterministic ordered merge, shared fuzzy graph. Exact search evaluates both directions; graph distancePairs counts unique unordered pairs. Scores use bounded tiles and 16 MiB mapping windows per worker; input reconstruction, identities and final graph remain resident. No approximate search, million-cell, Metal, embedding or biological qualification.")
            }
            graphHash = try write(graph, root: temp, name: "graph.json", maximum: 536_870_912)
        }
        let receipt = try VivoPCANeighborReceipt(schemaVersion: 1, input: inputHash,
            plan: write(plan, root: temp, name: "plan.json", maximum: 65_536), graph: graphHash,
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
        if plan.storage == .binary {
            let graph = try read(VivoPCAGraphStoreReport.self, root: root, name: "graph.json", maximum: 65_536)
            for (name, hash) in zip(VivoPCAGraphStore.files, [graph.neighbors, graph.bandwidths, graph.offsets, graph.edges]) {
                guard try VivoOmicsFileSnapshot.fingerprint(root.appendingPathComponent(name), maximumBytes: 4_096_000_000) == hash else {
                    throw VivoOmicsError.invalid("binary graph fingerprint differs")
                }
            }
        }
        return receipt
    }
}
