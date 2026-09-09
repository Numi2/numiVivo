import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellGraphStoreTests {
    @Test func binaryTransposeMatchesResidentWithTiesAndDisconnectedComponents() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        // Duplicates exercise zero rho/global-mean flooring; distant groups give
        // disconnected components and unequal incoming degrees.
        let scores: [[Double]] = [[0,0],[0,0],[0,0],[0,0],[1,0],[3,1],[100,0],[100,0],[101,0],[103,2],[106,1],[107,1]]
        let cells = scores.indices.map { VivoOmicsCellIdentity(sampleID: "s", barcode: "c\($0)") }
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 4
        let graph = try VivoSingleCellNeighbors.run(scores: scores, cells: cells, options: options)
        let writer = try VivoCountRecordWriter(root.appendingPathComponent("neighbors.bin"))
        for i in graph.neighborIndices.indices { try writer.append(row: i/options.neighbors, feature: graph.neighborIndices[i], bits: graph.neighborDistances[i].bitPattern) }
        let report = try VivoPCAGraphStore.build(root: root, rows: cells.count, dimensions: 2, options: options,
            distancePairs: graph.distancePairs, approximate: false, neighborsHash: writer.finish())
        #expect(report.connectedComponents == graph.connectedComponents && report.connectedComponents > 1)
        #expect(report.isolatedCells == graph.isolatedCells && report.connectivityEntries == graph.weights.count)
        let bandwidths = try VivoWindowedCountRecords(root.appendingPathComponent("bandwidths.bin"), entries: cells.count*3)
        for i in cells.indices {
            for (j, expected) in [graph.rhos[i], graph.sigmas[i], graph.kernelMassResiduals[i]].enumerated() {
                let record = try bandwidths.record(i*3+j)
                #expect(record.row == i && record.feature == j && record.bits == expected.bitPattern)
            }
        }
        let offsets = try VivoWindowedCountRecords(root.appendingPathComponent("offsets.bin"), entries: cells.count+1)
        let edges = try VivoWindowedCountRecords(root.appendingPathComponent("edges.bin"), entries: graph.weights.count)
        for i in 0...cells.count {
            let record = try offsets.record(i)
            #expect(record.row == i && record.feature == 0 && record.bits == UInt64(graph.rowOffsets[i]))
        }
        for i in cells.indices {
            for entry in graph.rowOffsets[i]..<graph.rowOffsets[i+1] {
                let record = try edges.record(entry)
                #expect(record.row == i && record.feature == graph.columnIndices[entry] && record.bits == graph.weights[entry].bitPattern)
            }
        }
        #expect(try !FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".graph-transpose-") })
    }
    @Test func streamingExactAndHNSWPropagateSinkFailure() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let scores = (0..<24).map { [Double($0), sin(Double($0))] }
        let source = root.appendingPathComponent("scores.bin")
        _ = try VivoH5ADPCA.writeMatrix(scores, columns: 2, to: source)
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 4
        enum SinkError: Error { case stopped }
        var rows: [Int] = []
        #expect(throws: SinkError.self) {
            try VivoHNSWNeighbors.stream(source: source, rows: scores.count, dimensions: 2, options: options, approximation: .init()) { row, _, _ in
                rows.append(row); if row == 2 { throw SinkError.stopped }
            }
        }
        #expect(rows == [0,1,2]); rows = []
        #expect(throws: SinkError.self) {
            try VivoWindowedPCANeighbors.stream(source: source, rows: scores.count, dimensions: 2, options: options, execution: .init()) { row, _, _ in
                rows.append(row); if row == 2 { throw SinkError.stopped }
            }
        }
        #expect(rows == [0,1,2])
        let plan = VivoPCANeighborPlan(storage: .binary)
        #expect(try VivoCanonicalJSON.decode(VivoPCANeighborPlan.self, from: VivoCanonicalJSON.encode(plan)) == plan)
    }
}
