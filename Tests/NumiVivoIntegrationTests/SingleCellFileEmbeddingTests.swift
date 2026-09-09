import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileEmbeddingTests {
    private func fixture(_ n: Int, root: URL) throws -> (VivoSingleCellNeighborGraph, [[Double]]) {
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        var scores: [[Double]] = []
        for i in 0..<n { scores.append([Double(i % 31), sin(Double(i) * 0.71), cos(Double(i) * 0.37)]) }
        let cells = (0..<n).map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 7
        let graph = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        let writer = try VivoCountRecordWriter(root.appendingPathComponent("edges.bin"))
        for i in 0..<n { for j in graph.rowOffsets[i]..<graph.rowOffsets[i+1] {
            try writer.append(row: i, feature: graph.columnIndices[j], bits: graph.weights[j].bitPattern)
        } }
        _ = try writer.finish(); return (graph, scores)
    }
    @Test func multipleWindowsPreserveFullTrajectoryAndRemoveScratch() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let (graph, scores) = try fixture(350, root: root)
        for dimensions in [2, 3] {
            var options = VivoSingleCellEmbeddingOptions(); options.epochs = 20; options.dimensions = dimensions
            let schedule = try VivoFileEmbeddingSchedule(edges: root.appendingPathComponent("edges.bin"), entries: graph.weights.count,
                cells: graph.cells.count, options: options, scratch: root, windowBytes: 16_384)
            #expect(schedule.fileBytes > schedule.windowBytes)
            let result = try VivoSingleCellEmbedding.run(cells: graph.cells, scores: scores, schedule: schedule, options: options)
            #expect(result == (try VivoSingleCellEmbedding.run(graph, scores: scores, options: options)))
            try schedule.remove()
            #expect(try !FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".embedding-schedule-") })
        }
    }
    @Test func exhaustedBudgetAndMalformedEdgeRejectWithoutScratch() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let (graph, _) = try fixture(16, root: root)
        var options = VivoSingleCellEmbeddingOptions(); options.maximumUpdates = 1
        #expect(throws: (any Error).self) { try VivoFileEmbeddingSchedule(edges: root.appendingPathComponent("edges.bin"), entries: graph.weights.count, cells: 16, options: options, scratch: root) }
        #expect(try !FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".embedding-schedule-") })
        var bytes = try Data(contentsOf: root.appendingPathComponent("edges.bin")); bytes[4] = bytes[0]
        try bytes.write(to: root.appendingPathComponent("edges.bin")); options.maximumUpdates = 1_000_000
        #expect(throws: (any Error).self) { try VivoFileEmbeddingSchedule(edges: root.appendingPathComponent("edges.bin"), entries: graph.weights.count, cells: 16, options: options, scratch: root) }
    }
    @Test func callbackFailureReleasesScheduleAndExtendedBudgetIsExplicit() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let (graph, _) = try fixture(16, root: root)
        func fail() throws {
            let schedule = try VivoFileEmbeddingSchedule(edges: root.appendingPathComponent("edges.bin"), entries: graph.weights.count, cells: 16, options: .init(), scratch: root)
            try schedule.visit(epoch: 2) { _, _, _ in throw VivoOmicsError.invalid("injected callback failure") }
        }
        #expect(throws: (any Error).self) { try fail() }
        #expect(try !FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".embedding-schedule-") })
        var options = VivoSingleCellEmbeddingOptions(); #expect(options.maximumUpdates == 200_000_000)
        options.maximumUpdates = 20_000_000_000; try options.validate()
        options.maximumUpdates += 1; #expect(throws: (any Error).self) { try options.validate() }
    }
}
