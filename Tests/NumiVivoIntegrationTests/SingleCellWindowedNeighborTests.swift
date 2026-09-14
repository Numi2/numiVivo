import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellWindowedNeighborTests {
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_METAL"] == "1"))
    func metalFP32TilesMatchCPUOracleOnBoundedFixture() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        var scores: [[Double]] = []
        for i in 0..<257 {
            // Integer-valued coordinates keep the membership oracle separated
            // while still exercising the FP32 upload and tiled dispatch path.
            let a = Double(i % 13), b = Double(i % 7)
            let c = Double(i % 5), d = Double(i % 3)
            scores.append([a, b, c, d])
        }
        let cells = scores.indices.map { VivoOmicsCellIdentity(sampleID: "s", barcode: "metal-c\($0)") }
        let source = root.appendingPathComponent("scores.bin")
        _ = try VivoH5ADPCA.writeMatrix(scores, columns: 4, to: source)
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 7
        let cpu = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        var execution = VivoPCANeighborExecution()
        execution.backend = .metalFP32; execution.workers = 1
        execution.queryBlockRows = 23; execution.candidateBlockRows = 31
        let metal = try VivoWindowedPCANeighbors.run(source: source, cells: cells, dimensions: 4,
            options: options, execution: execution)
        #expect(metal.method == "metal-FP32-windowed-PCA-knn-v1")
        #expect(metal.neighborIndices == cpu.neighborIndices)
        #expect(metal.rowOffsets == cpu.rowOffsets)
        #expect(metal.columnIndices == cpu.columnIndices)
        #expect(metal.weights.count == cpu.weights.count)
        for (actual, expected) in zip(metal.neighborDistances, cpu.neighborDistances) {
            #expect(abs(actual - expected) <= 1e-6)
        }
        for (actual, expected) in zip(metal.weights, cpu.weights) {
            #expect(abs(actual - expected) <= 1e-6)
        }
        #expect(metal.distancePairs == cpu.distancePairs)
    }

    @Test func parallelTilesMatchResidentGraphIncludingDistanceTies() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        var scores: [[Double]] = []
        for i in 0..<257 {
            let a = Double(i % 13), b = Double(i % 7)
            let c = Double(i % 5), d = Double(i % 3)
            scores.append([a, b, c, d])
        }
        let cells = scores.indices.map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        let url = root.appendingPathComponent("scores.bin")
        _ = try VivoH5ADPCA.writeMatrix(scores, columns: 4, to: url)
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 7
        let expected = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        for (workers, query, candidate) in [(1, 8, 17), (2, 23, 31), (4, 64, 512)] {
            var execution = VivoPCANeighborExecution(); execution.workers = workers
            execution.queryBlockRows = query; execution.candidateBlockRows = candidate
            let graph = try VivoWindowedPCANeighbors.run(source: url, cells: cells, dimensions: 4, options: options, execution: execution)
            #expect(graph == expected)
        }
        options.maximumDistancePairs = 1
        #expect(throws: (any Error).self) { try VivoWindowedPCANeighbors.run(source: url, cells: cells, dimensions: 4, options: options, execution: .init()) }
        options.maximumDistancePairs = 50_000_000
        var bytes = try Data(contentsOf: url); bytes[0] = 1; try bytes.write(to: url)
        #expect(throws: (any Error).self) { try VivoWindowedPCANeighbors.run(source: url, cells: cells, dimensions: 4, options: options, execution: .init()) }
    }
    @Test func scoreTilesCrossMappingBoundaryAndRewind() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let url = root.appendingPathComponent("scores.bin"), rows = 17_000, d = 64
        let writer = try VivoCountRecordWriter(url)
        for row in 0..<rows { for c in 0..<d { try writer.append(row: row, feature: c, bits: (Double(row)-Double(c)*0.25).bitPattern) } }
        _ = try writer.finish()
        let reader = try VivoPCAScoreReader(url, rows: rows, columns: d)
        for range in [16_370..<16_400, 0..<17, 16_995..<17_000] {
            let values = try reader.readRows(range)
            for (i, row) in range.enumerated() { for c in 0..<d { #expect(values[i*d+c] == Double(row)-Double(c)*0.25) } }
        }
        #expect(throws: (any Error).self) { try reader.readRows(0..<8_193) }
    }
}
