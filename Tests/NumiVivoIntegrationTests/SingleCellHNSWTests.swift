import Foundation
import Testing
import NumiVivoCore
@testable import NumiVivoKit

@Suite struct SingleCellHNSWTests {
    @Test func nativeMillionRowAdmissionAndLargeResourceBounds() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("scores.bin")
        // A sparse file is enough to exercise real index admission and budget
        // rejection above both former caps, without computing a million-cell graph.
        _ = try VivoH5ADPCA.writeMatrix([[0], [1]], columns: 1, to: source)
        let file = try FileHandle(forWritingTo: source)
        try file.truncate(atOffset: UInt64(1_000_001 * 16))
        try file.seek(toOffset: 0)
        var tile = Data()
        for row in 0..<256 {
            var record = (UInt64(row).littleEndian, Double(row).bitPattern.littleEndian)
            withUnsafeBytes(of: &record) { tile.append(contentsOf: $0) }
        }
        try file.write(contentsOf: tile); try file.close()
        var native = NVivoHNSWOptions()
        native.struct_size = UInt32(MemoryLayout<NVivoHNSWOptions>.size); native.abi_version = 1
        native.rows = 1_000_001; native.dimensions = 1; native.neighbors = 128
        native.connections = 8; native.ef_construction = 8; native.ef_search = 128
        native.seed = 7; native.maximum_distance_evaluations = 1; native.score_cache_bytes = 262_144
        var report = NVivoHNSWReport()
        func invoke(_ path: URL) -> Int32 {
            path.path.withCString {
                nvivo_omics_hnsw_neighbors_stream($0, &native, { _, _, _, _, _ in 1 }, &report, nil, nil)
            }
        }
        #expect(invoke(source) == 2) // Reaches metric budget; no row callback publishes.
        native.rows = UInt32(VivoPCAStorageLimits.maximumRows)
        native.dimensions = UInt32(VivoPCAStorageLimits.maximumColumns)
        native.maximum_distance_evaluations = NVIVO_OMICS_HNSW_MAXIMUM_DISTANCE_EVALUATIONS
        native.score_cache_bytes = NVIVO_OMICS_HNSW_MAXIMUM_CACHE_BYTES
        let missing = root.appendingPathComponent("missing.bin")
        #expect(invoke(missing) == 3) // Accepted axes/resources, then missing input.
        native.rows += 1; #expect(invoke(missing) == 1); native.rows -= 1
        native.maximum_distance_evaluations += 1
        #expect(invoke(missing) == 1); native.maximum_distance_evaluations -= 1
        native.score_cache_bytes += 1; #expect(invoke(missing) == 1)
        var options = VivoHNSWOptions()
        #expect(options.maximumDistanceEvaluations == 500_000_000 && options.scoreCacheBytes == 33_554_432)
        options.maximumDistanceEvaluations = VivoHNSWOptions.maximumSupportedDistanceEvaluations
        options.scoreCacheBytes = VivoHNSWOptions.maximumSupportedCacheBytes
        try options.validate(neighbors: 128)
        options.maximumDistanceEvaluations += 1
        #expect(throws: (any Error).self) { try options.validate(neighbors: 128) }
        options.maximumDistanceEvaluations -= 1; options.scoreCacheBytes += 1
        #expect(throws: (any Error).self) { try options.validate(neighbors: 128) }
    }
    @Test func exhaustiveSmallSearchMatchesExactAndEvictsScoreTiles() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        var scores: [[Double]] = []
        for i in 0..<600 {
            var row = [Double](repeating: 0, count: 64)
            row[0] = sin(Double(i)*0.71); row[1] = cos(Double(i)*0.43); row[2] = Double(i)/37
            scores.append(row)
        }
        let cells = scores.indices.map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        let source = root.appendingPathComponent("scores.bin")
        _ = try VivoH5ADPCA.writeMatrix(scores, columns: 64, to: source)
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 7
        var ann = VivoHNSWOptions(); ann.searchWidth = 600; ann.scoreCacheBytes = 262_144
        let expected = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        let (graph, report) = try VivoHNSWNeighbors.run(source: source, cells: cells, dimensions: 64, options: options, approximation: ann)
        #expect(graph.neighborIndices == expected.neighborIndices)
        #expect(graph.neighborDistances == expected.neighborDistances)
        #expect(graph.rowOffsets == expected.rowOffsets && graph.columnIndices == expected.columnIndices && graph.weights == expected.weights)
        #expect(graph.method.contains("approximate-HNSW"))
        #expect(report.peakCachedScoreBytes <= ann.scoreCacheBytes && report.scoreReadBytes > 600*64*16)
        #expect(report.constructionDistances > 0 && report.queryDistances > 0)
        let replay = try VivoHNSWNeighbors.run(source: source, cells: cells, dimensions: 64, options: options, approximation: ann)
        #expect(replay.0 == graph && replay.1 == report)
        ann.maximumDistanceEvaluations = 1
        #expect(throws: (any Error).self) { try VivoHNSWNeighbors.run(source: source, cells: cells, dimensions: 64, options: options, approximation: ann) }
        ann.maximumDistanceEvaluations = 500_000_000
        var bytes = try Data(contentsOf: source); bytes[0] = 1; try bytes.write(to: source)
        #expect(throws: (any Error).self) { try VivoHNSWNeighbors.run(source: source, cells: cells, dimensions: 64, options: options, approximation: ann) }
    }
    @Test func approximatePlanHasExplicitSerialAndResourceContract() throws {
        var execution = VivoPCANeighborExecution(); execution.workers = 1
        let plan = VivoPCANeighborPlan(execution: execution, approximation: .init())
        try plan.validate()
        #expect(try VivoCanonicalJSON.decode(VivoPCANeighborPlan.self, from: VivoCanonicalJSON.encode(plan)) == plan)
        #expect(throws: (any Error).self) { try VivoPCANeighborPlan(approximation: .init()).validate() }
        var invalid = VivoHNSWOptions(); invalid.scoreCacheBytes = 1
        #expect(throws: (any Error).self) { try invalid.validate(neighbors: 15) }
    }
}
