import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellHNSWTests {
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
