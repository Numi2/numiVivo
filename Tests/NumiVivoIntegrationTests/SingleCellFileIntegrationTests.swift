import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileIntegrationTests {
    @Test func mappedRowsSurviveWindowChangesAndRejectInvalidAccess() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let matrix = try VivoIntegrationMatrix(rows: 513, columns: 88, scratch: root, windowBytes: 16_384)
        for i in 0..<513 { try matrix.setRow(i, (0..<88).map { Double(i * 88 + $0) / 17 }) }
        for i in stride(from: 512, through: 0, by: -1) {
            #expect(try matrix.row(i) == (0..<88).map { Double(i * 88 + $0) / 17 })
            #expect(try matrix.value(i, 87) == Double(i * 88 + 87) / 17)
        }
        #expect(throws: (any Error).self) { try matrix.row(-1) }
        #expect(throws: (any Error).self) { try matrix.value(513, 0) }
        #expect(throws: (any Error).self) { try matrix.value(0, 88) }
        #expect(throws: (any Error).self) { try matrix.setRow(0, [.nan]) }
        try matrix.remove(); try matrix.remove()
        #expect(throws: (any Error).self) { try matrix.row(0) }
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
    }
    @Test func fileAndResidentTrajectoriesMatchAcrossWindowsAndSeeds() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let n = 641, d = 5
        let cells = (0..<n).map { VivoOmicsCellIdentity(sampleID: "s\($0 % 3)", barcode: "c\($0)") }
        let samples = (0..<3).map { VivoOmicsSample(id: "s\($0)", biologicalReplicateID: "d\($0)", donorID: "d\($0)", condition: "shared", batchID: "b\($0)", organism: "human") }
        for seed: UInt64 in [7, 19, 41] {
            var options = VivoSingleCellIntegrationOptions(); options.clusters = 17; options.seed = seed
            var owned: [VivoIntegrationMatrix] = []
            func make(_ rows: Int, _ columns: Int) throws -> VivoIntegrationMatrix {
                let matrix = try VivoIntegrationMatrix(rows: rows, columns: columns, scratch: root, windowBytes: 16_384)
                owned.append(matrix); return matrix
            }
            let resident = try VivoIntegrationMatrix(rows: n, columns: d), file = try make(n, d)
            for i in 0..<n {
                let row = (0..<d).map { j in sin(Double(i * (j + 1)) / 17) + Double(i % 3) / Double(j + 1) }
                try resident.setRow(i, row); try file.setRow(i, row)
            }
            let a = try VivoSingleCellIntegration.run(cells: cells, x: resident, samples: samples, options: options) { try VivoIntegrationMatrix(rows: $0, columns: $1) }
            let b = try VivoSingleCellIntegration.run(cells: cells, x: file, samples: samples, options: options, matrix: make)
            #expect(try a.materialize(cells: cells, options: options) == b.materialize(cells: cells, options: options))
            for matrix in owned { try matrix.remove() }
            #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
        }
    }
    @Test func dependencyAndAdmissionChecksRemainExplicit() throws {
        #expect(throws: (any Error).self) { try VivoPCAIntegrationPlan(inputKind: .integrated).validate() }
        #expect(throws: (any Error).self) { try VivoPCANeighborPlan(inputKind: .integrated).validate() }
        var neighbors = VivoSingleCellNeighborOptions(); neighbors.representation = .integrated
        try VivoPCANeighborPlan(inputKind: .integrated, neighbors: neighbors).validate()
        #expect(throws: (any Error).self) { try VivoPCANeighborPlan(inputKind: .fitted, neighbors: neighbors).validate() }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoPCAIntegrationPlan.self, from: Data("{\"schemaVersion\":1,\"doner\":1}".utf8)) }
        var options = VivoSingleCellIntegrationOptions()
        #expect(throws: (any Error).self) { try VivoSingleCellIntegration.validateAxes(rows: 1_000_000, columns: 20, options: options) }
        options.maximumWork = 100_000_000_000
        try VivoSingleCellIntegration.validateAxes(rows: 1_000_000, columns: 20, options: options)
        #expect(throws: (any Error).self) { try VivoSingleCellIntegration.validateAxes(rows: Int.max, columns: 20, options: options) }
    }
}
