import Foundation
import NumiVivoCore

public struct VivoHNSWOptions: Codable, Sendable, Equatable {
    public var connections: Int = 16
    public var constructionWidth: Int = 200
    public var searchWidth: Int = 128
    public var seed: UInt64 = 7
    public var maximumDistanceEvaluations: Int = 500_000_000
    public var scoreCacheBytes: Int = 33_554_432
    public init() {}
    private enum CodingKeys: String, CodingKey { case connections, constructionWidth, searchWidth, seed, maximumDistanceEvaluations, scoreCacheBytes }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["connections", "constructionWidth", "searchWidth", "seed", "maximumDistanceEvaluations", "scoreCacheBytes"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        connections = try c.decodeIfPresent(Int.self, forKey: .connections) ?? 16
        constructionWidth = try c.decodeIfPresent(Int.self, forKey: .constructionWidth) ?? 200
        searchWidth = try c.decodeIfPresent(Int.self, forKey: .searchWidth) ?? 128
        seed = try c.decodeIfPresent(UInt64.self, forKey: .seed) ?? 7
        maximumDistanceEvaluations = try c.decodeIfPresent(Int.self, forKey: .maximumDistanceEvaluations) ?? 500_000_000
        scoreCacheBytes = try c.decodeIfPresent(Int.self, forKey: .scoreCacheBytes) ?? 33_554_432
    }
    public func validate(neighbors: Int) throws {
        guard (8...64).contains(connections), (connections...512).contains(constructionWidth),
              (neighbors...1024).contains(searchWidth), (1...2_000_000_000).contains(maximumDistanceEvaluations),
              (262_144...67_108_864).contains(scoreCacheBytes) else { throw VivoOmicsError.invalid("HNSW options or resource bounds") }
    }
}
public struct VivoHNSWReport: Codable, Sendable, Equatable {
    public let implementation: String
    public let options: VivoHNSWOptions
    public let constructionDistances: Int
    public let queryDistances: Int
    public let scoreReadBytes: Int
    public let scoreReadCalls: Int
    public let scoreCacheHits: Int
    public let peakCachedScoreBytes: Int
    public let indexStorageBytes: Int
    public let qualification: String
}

enum VivoHNSWNeighbors {
    private final class Output {
        let sink: (Int, [Int], [Double]) throws -> Void
        var error: Error?
        init(_ sink: @escaping (Int, [Int], [Double]) throws -> Void) { self.sink = sink }
    }
    static func stream(source: URL, rows n: Int, dimensions: Int, options: VivoSingleCellNeighborOptions,
                       approximation: VivoHNSWOptions, sink: @escaping (Int, [Int], [Double]) throws -> Void) throws -> VivoHNSWReport {
        try options.validate(); try approximation.validate(neighbors: options.neighbors)
        let k = options.neighbors
        guard n >= k, n <= VivoPCAStorageLimits.maximumRows, (1...64).contains(dimensions) else {
            throw VivoOmicsError.limit("HNSW axes or representation")
        }
        try Task.checkCancellation()
        var native = NVivoHNSWOptions()
        native.struct_size = UInt32(MemoryLayout<NVivoHNSWOptions>.size); native.abi_version = 1
        native.rows = UInt32(n); native.dimensions = UInt32(dimensions); native.neighbors = UInt32(k)
        native.connections = UInt32(approximation.connections); native.ef_construction = UInt32(approximation.constructionWidth)
        native.ef_search = UInt32(approximation.searchWidth); native.seed = approximation.seed
        native.maximum_distance_evaluations = UInt64(approximation.maximumDistanceEvaluations); native.score_cache_bytes = UInt64(approximation.scoreCacheBytes)
        var raw = NVivoHNSWReport()
        let output = Output(sink)
        // Retain the callback owner for the entire synchronous C invocation.
        let status = withExtendedLifetime(output) {
            source.path.withCString { path in
                nvivo_omics_hnsw_neighbors_stream(path, &native, { row, indices, distances, count, opaque in
                    guard let indices, let distances, let opaque else { return 1 }
                    let owner = Unmanaged<Output>.fromOpaque(opaque).takeUnretainedValue()
                    do {
                        try owner.sink(Int(row), UnsafeBufferPointer(start: indices, count: Int(count)).map(Int.init),
                                       Array(UnsafeBufferPointer(start: distances, count: Int(count))))
                        return 0
                    } catch { owner.error = error; return 1 }
                }, &raw, { _ in Task.isCancelled ? 1 : 0 }, Unmanaged.passUnretained(output).toOpaque())
            }
        }
        if let error = output.error { throw error }
        switch status {
        case 0: break
        case 2: throw VivoOmicsError.limit("HNSW distance-evaluation budget")
        case 3: throw VivoOmicsError.invalid("HNSW score file, coordinates or values")
        case 4: throw VivoOmicsError.invalid("HNSW nonfinite distance")
        case 5: throw CancellationError()
        case 7: throw VivoOmicsError.limit("HNSW allocation failed")
        default: throw VivoOmicsError.invalid("HNSW native arguments or index failure: \(status)")
        }
        try Task.checkCancellation()
        let report = VivoHNSWReport(implementation: "hnswlib-0.8.0-3f3429661187e4c24a490a0f148fc6bc89042b3d-fp64-id-space-v1", options: approximation,
            constructionDistances: Int(raw.construction_distances), queryDistances: Int(raw.query_distances), scoreReadBytes: Int(raw.score_read_bytes),
            scoreReadCalls: Int(raw.score_read_calls), scoreCacheHits: Int(raw.score_cache_hits), peakCachedScoreBytes: Int(raw.peak_cached_score_bytes),
            indexStorageBytes: Int(raw.index_storage_bytes),
            qualification: "Serial source-order HNSW construction/query with pinned seed, FP64 metric and bounded LRU score tiles. Index topology remains resident; graph residency depends on the selected output storage. Search is approximate; recall needs independent measurement. Index storage bytes exclude allocator, mutex and container overhead. No million-cell or biological qualification.")
        return report
    }
    static func run(source: URL, cells: [VivoOmicsCellIdentity], dimensions: Int, options: VivoSingleCellNeighborOptions,
                    approximation: VivoHNSWOptions) throws -> (VivoSingleCellNeighborGraph, VivoHNSWReport) {
        try options.validate()
        guard cells.count <= 4_000_000/options.neighbors else { throw VivoOmicsError.limit("resident HNSW graph-entry bound") }
        var indices: [Int] = [], distances: [Double] = []
        let report = try stream(source: source, rows: cells.count, dimensions: dimensions, options: options, approximation: approximation) { _, i, d in
            indices.append(contentsOf: i); distances.append(contentsOf: d)
        }
        let graph = try VivoSingleCellNeighbors.finish(indices: indices, distances: distances, cells: cells, dimensions: dimensions,
            options: options, distancePairs: report.constructionDistances+report.queryDistances,
            method: "approximate-HNSW-euclidean-PCA-knn-umap-fuzzy-union-v1",
            qualification: "Approximate HNSW neighbors with FP64 returned distances and the shared fuzzy graph. distancePairs counts all construction/query metric evaluations, including repeats; it is not a unique-pair count. No exact-neighbor, embedding, clustering, integration or biological qualification.")
        return (graph, report)
    }
}
