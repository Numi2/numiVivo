import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileClusteringTests {
    @Test func bufferedRecordsPreserveBitsAcrossPagesAndRejectChangedFiles() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let path = root.appendingPathComponent("records.bin")
        let writer = try VivoCountRecordWriter(path)
        for i in 0..<600 { try writer.append(row: i, feature: 599 - i, bits: UInt64.max - UInt64(i)) }
        _ = try writer.finish()
        let reader = try VivoBufferedCountRecords(path, entries: 600)
        for i in [0, 255, 256, 599, 512, 0] {
            let record = try reader.record(i)
            #expect(record.row == i && record.feature == 599 - i && record.bits == UInt64.max - UInt64(i))
        }
        #expect(reader.bufferLoads == 4)
        #expect(reader.loadedBytes == 3 * 4_096 + 88 * 16)
        #expect(throws: (any Error).self) { try reader.record(-1) }
        #expect(throws: (any Error).self) { try reader.record(600) }
        #expect(throws: (any Error).self) { try VivoBufferedCountRecords(path, entries: 599) }
        #expect(throws: (any Error).self) { try VivoBufferedCountRecords(root, entries: 0) }
        let link = root.appendingPathComponent("linked.bin")
        try FileManager.default.createSymbolicLink(at: link, withDestinationURL: path)
        #expect(throws: (any Error).self) { try VivoBufferedCountRecords(link, entries: 600) }
        let changed = try FileHandle(forWritingTo: path)
        try changed.truncate(atOffset: 16); try changed.close()
        #expect(throws: (any Error).self) { try reader.record(599) }
    }

    private func write(_ rows: [[VivoSingleCellClustering.Edge]], root: URL) throws {
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        let edges = try VivoCountRecordWriter(root.appendingPathComponent("edges.bin")), offsets = try VivoCountRecordWriter(root.appendingPathComponent("offsets.bin"))
        for i in rows.indices {
            try offsets.append(row: i, feature: 0, bits: UInt64(edges.entries))
            for edge in rows[i] { try edges.append(row: i, feature: edge.column, bits: edge.weight.bitPattern) }
        }
        try offsets.append(row: rows.count, feature: 0, bits: UInt64(edges.entries))
        _ = try edges.finish(); _ = try offsets.finish()
    }
    @Test func fileAndResidentResultsMatchAndScratchIsRemoved() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        var scores: [[Double]] = []
        for i in 0..<90 {
            let x = Double(i % 10), y = Double(i / 10) * 30.0, z = sin(Double(i) * 0.71)
            scores.append([x, y, z])
        }
        let cells = scores.indices.map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 7
        let graph = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        let rows = cells.indices.map { i in (graph.rowOffsets[i]..<graph.rowOffsets[i+1]).map { VivoSingleCellClustering.Edge(column: graph.columnIndices[$0], weight: graph.weights[$0]) } }
        try write(rows, root: root.appendingPathComponent("input"))
        let work = VivoClusteringWork(maximumEdgeVisits: 1_000_000)
        let input = try VivoFileClusteringGraph(root: root.appendingPathComponent("input"), rows: cells.count, entries: graph.weights.count, work: work)
        let result = try VivoSingleCellClustering.run(original: input, cells: cells, options: .init(), scratch: root)
        #expect(result == (try VivoSingleCellClustering.run(graph, options: .init())))
        #expect(work.aggregatedLevels > 0 && work.edgeVisits > graph.weights.count && work.maximumEdgeMapEntries > 0)
        #expect(work.edgeBufferLoads == 0 && work.edgeBytesRead == 0)
        #expect(try !FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".clustering-level-") })
    }
    @Test func graphAboveOneWindowUsesBoundedBufferedReads() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let entries = VivoWindowedCountRecords.windowBytes / 16 + 1
        let offsets = try VivoCountRecordWriter(root.appendingPathComponent("offsets.bin"))
        for (row, value) in [0, 1, entries].enumerated() {
            try offsets.append(row: row, feature: 0, bits: UInt64(value))
        }
        _ = try offsets.finish()
        let path = root.appendingPathComponent("edges.bin")
        #expect(FileManager.default.createFile(atPath: path.path, contents: nil))
        let file = try FileHandle(forWritingTo: path)
        try file.truncate(atOffset: UInt64(entries * 16))
        try file.seek(toOffset: 0)
        try file.write(contentsOf: Data([0,0,0,0,1,0,0,0,0,0,0,0,0,0,240,63]))
        try file.close()
        let work = VivoClusteringWork()
        let graph = try VivoFileClusteringGraph(root: root, rows: 2, entries: entries, work: work)
        var edges = 0
        try graph.forEachEdge(in: 0) { edge in
            #expect(edge.column == 1 && edge.weight == 1); edges += 1
        }
        #expect(edges == 1 && work.edgeBufferLoads == 1 && work.edgeBytesRead == 4_096)
        // The sparse tail deliberately has invalid coordinates; it must not be
        // accepted merely because the byte count and backend admission succeed.
        #expect(throws: (any Error).self) { try graph.forEachEdge(in: 1) { _ in } }
    }
    @Test func diskAggregationPreservesSummationOrderAndDiagonalMass() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let rows: [[VivoSingleCellClustering.Edge]] = [
            [.init(column: 1, weight: 1e16), .init(column: 2, weight: 1), .init(column: 3, weight: 1)],
            [.init(column: 0, weight: 1e16), .init(column: 2, weight: 1), .init(column: 3, weight: 1)],
            [.init(column: 0, weight: 1), .init(column: 1, weight: 1), .init(column: 3, weight: 1)],
            [.init(column: 0, weight: 1), .init(column: 1, weight: 1), .init(column: 2, weight: 1)]]
        try write(rows, root: root.appendingPathComponent("input"))
        let file = try VivoFileClusteringGraph(root: root.appendingPathComponent("input"), rows: 4, entries: 12, work: .init())
        let expected = try VivoResidentClusteringGraph(rows).aggregate(labels: [0,1,0,1], groups: 2, scratch: nil)
        let observed = try file.aggregate(labels: [0,1,0,1], groups: 2, scratch: root)
        for row in 0..<2 {
            var indices: [Int] = [], weights: [UInt64] = [], observedIndices: [Int] = [], observedWeights: [UInt64] = []
            try expected.forEachEdge(in: row) { indices.append($0.column); weights.append($0.weight.bitPattern) }
            try observed.forEachEdge(in: row) { observedIndices.append($0.column); observedWeights.append($0.weight.bitPattern) }
            #expect(indices == observedIndices && weights == observedWeights)
            #expect(indices.contains(row))
        }
    }
    @Test func malformedCSRAndExhaustedWorkAreRejected() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let rows: [[VivoSingleCellClustering.Edge]] = [[.init(column: 1, weight: 1)], [.init(column: 0, weight: 1)]]
        try write(rows, root: root)
        let graph = try VivoFileClusteringGraph(root: root, rows: 2, entries: 2, work: .init(maximumEdgeVisits: 1))
        let cells = (0..<2).map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        #expect(throws: (any Error).self) { try VivoSingleCellClustering.run(original: graph, cells: cells, options: .init(), scratch: root) }
        var bad = try Data(contentsOf: root.appendingPathComponent("offsets.bin")); bad[8] = 1; try bad.write(to: root.appendingPathComponent("offsets.bin"))
        #expect(throws: (any Error).self) { try VivoFileClusteringGraph(root: root, rows: 2, entries: 2, work: .init()) }
        #expect(throws: (any Error).self) { try VivoPCAGraphClusteringPlan(maximumEdgeVisits: 0).validate() }
    }
}
