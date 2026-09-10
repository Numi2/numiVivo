import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellFileIntegrationTests {
    @Test func streamedRidgePassesPreserveScalarOrderAndInactiveLevels() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let n = 1_031, d = 5, k = 7, levels = 3
        let batch = (0..<n).map { ($0 * 7 + $0 / 11) % levels }
        let x = try VivoIntegrationMatrix(rows: n, columns: d, scratch: root, windowBytes: 16_384)
        let r = try VivoIntegrationMatrix(rows: n, columns: k, scratch: root, windowBytes: 16_384)
        let scores = try VivoIntegrationMatrix(rows: n, columns: d, scratch: root, windowBytes: 16_384)
        for i in 0..<n {
            try x.setRow(i, (0..<d).map { $0 == 0 ? -0.0 : sin(Double(i * ($0 + 1)) / 17) })
            try r.setRow(i, (0..<k).map { Double((i * 3 + $0) % 19) / 19 })
        }
        try scores.copy(from: x)
        let active = [0, 2, 5, 6]
        var expected = [Double](repeating: 0, count: k * levels * d)
        for c in active { for i in 0..<n {
            let weight = try r.value(i, c), row = try x.row(i)
            for j in 0..<d { expected[(c * levels + batch[i]) * d + j] += weight * row[j] }
        } }
        let actual = try VivoSingleCellIntegration.streamedRidgeSums(x: x, memberships: r,
            batch: batch, levels: levels, activeClusters: active)
        #expect(actual.map(\.bitPattern) == expected.map(\.bitPattern))
        var effects = [[Double]?](repeating: nil, count: k * levels)
        for c in active { for b in 0..<levels where (c + b) % 2 == 0 {
            effects[c * levels + b] = (0..<d).map { $0 == 0 ? -0.0 : Double(c + b + $0) / 13 }
        } }
        let old = try VivoIntegrationMatrix(rows: n, columns: d)
        try old.copy(from: x)
        for c in 0..<k { for b in 0..<levels {
            guard let effect = effects[c * levels + b] else { continue }
            for i in 0..<n where batch[i] == b {
                let weight = try r.value(i, c); var row = try old.row(i)
                for j in 0..<d { row[j] -= weight * effect[j] }
                try old.setRow(i, row)
            }
        } }
        try VivoSingleCellIntegration.applyStreamedRidgeEffects(scores: scores, memberships: r,
            batch: batch, levels: levels, effects: effects)
        for i in 0..<n { #expect(try scores.row(i).map(\.bitPattern) == old.row(i).map(\.bitPattern)) }
    }

    @Test func streamedRidgeRejectsInvalidTablesBeforeWrites() throws {
        let x = try VivoIntegrationMatrix(rows: 3, columns: 2)
        let r = try VivoIntegrationMatrix(rows: 3, columns: 2)
        for i in 0..<3 { try x.setRow(i, [Double(i), -0.0]); try r.setRow(i, [0.5, 0.5]) }
        #expect(throws: (any Error).self) {
            try VivoSingleCellIntegration.streamedRidgeSums(x: x, memberships: r, batch: [0, 1, 0], levels: 2, activeClusters: [0, 0])
        }
        #expect(throws: (any Error).self) {
            try VivoSingleCellIntegration.streamedRidgeSums(x: x, memberships: r, batch: [0, 2, 0], levels: 2, activeClusters: [0])
        }
        #expect(throws: (any Error).self) {
            try VivoSingleCellIntegration.applyStreamedRidgeEffects(scores: x, memberships: r,
                batch: [0, 1, 0], levels: 2, effects: [[1, 1], nil, nil, [.nan, 1]])
        }
        for i in 0..<3 { #expect(try x.row(i).map(\.bitPattern) == [Double(i), -0.0].map(\.bitPattern)) }
    }

    @Test func cancelledStreamedRidgePassDoesNotMutateScores() async throws {
        let task = Task {
            let scores = try VivoIntegrationMatrix(rows: 3, columns: 2)
            let memberships = try VivoIntegrationMatrix(rows: 3, columns: 2)
            withUnsafeCurrentTask { $0?.cancel() }
            #expect(throws: CancellationError.self) {
                try VivoSingleCellIntegration.streamedRidgeSums(x: scores, memberships: memberships,
                    batch: [0, 1, 0], levels: 2, activeClusters: [0, 1])
            }
            #expect(throws: CancellationError.self) {
                try VivoSingleCellIntegration.applyStreamedRidgeEffects(scores: scores, memberships: memberships,
                    batch: [0, 1, 0], levels: 2, effects: [[1, 1], nil, nil, nil])
            }
            for i in 0..<3 { #expect(try scores.row(i) == [0, 0]) }
        }
        try await task.value
    }

    @Test func consumedScratchPreservesOutputBitsAndFailureOwnership() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let matrix = try VivoIntegrationMatrix(rows: 513, columns: 3, scratch: root, windowBytes: 16_384)
        defer { try? matrix.remove() }
        for i in 0..<513 { try matrix.setRow(i, [Double(i) / 17, -0.0, Double.leastNonzeroMagnitude]) }
        let expected = root.appendingPathComponent("expected.bin")
        let fingerprint = try matrix.writeRecords(to: expected)
        #expect(throws: (any Error).self) { try matrix.writeRecordsAndRemove(to: expected) }
        #expect(try matrix.row(0).map(\.bitPattern) == [0.0, -0.0, Double.leastNonzeroMagnitude].map(\.bitPattern))
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).filter { $0.hasPrefix(".integration-matrix-") }.count == 1)
        let output = root.appendingPathComponent("output.bin")
        #expect(try matrix.writeRecordsAndRemove(to: output) == fingerprint)
        #expect(try Data(contentsOf: output) == Data(contentsOf: expected))
        #expect(throws: (any Error).self) { try matrix.row(0) }
        try matrix.remove()
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).sorted() == ["expected.bin", "output.bin"])
    }
    @Test func cancelledSerializationLeavesScratchForOwnerCleanup() async throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let task = Task {
            let matrix = try VivoIntegrationMatrix(rows: 1, columns: 3, scratch: root)
            defer { try? matrix.remove() }
            try matrix.setRow(0, [1, 2, 3])
            withUnsafeCurrentTask { $0?.cancel() }
            #expect(throws: CancellationError.self) { try matrix.writeRecordsAndRemove(to: root.appendingPathComponent("cancelled.bin")) }
            #expect(try matrix.row(0) == [1, 2, 3])
        }
        try await task.value
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).allSatisfy { !$0.hasPrefix(".integration-matrix-") })
    }
    @Test func completeCohortWitnessAdmissionUsesSupportedAxes() throws {
        var options = VivoSingleCellIntegrationOptions()
        options.clusters = VivoSingleCellIntegrationOptions.maximumClusters
        options.maximumWork = 100_000_000_000
        try VivoSingleCellIntegration.validateAxes(rows: 1_612_594, columns: 20, options: options)
        #expect(1_612_594 * options.clusters * 16 <= VivoIntegrationStorageLimits.maximumMembershipBytes)
        #expect(VivoIntegrationStorageLimits.maximumMembershipBytes == 3_200_000_000)
        #expect(VivoIntegrationStorageLimits.maximumAnchorBytes == 3_200_000_000)
        options.clusters += 1
        #expect(throws: (any Error).self) { try options.validate() }
        var mnn = VivoMNNIntegrationOptions()
        mnn.neighbors = VivoMNNIntegrationOptions.maximumNeighbors + 1
        #expect(throws: (any Error).self) { try mnn.validate() }

        // Exercise the real bounded snapshot reader beyond the previous ceiling
        // without allocating a resident multi-gigabyte fixture. APFS keeps the
        // zero-filled source and its immutable clone sparse.
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("memberships.bin")
        #expect(FileManager.default.createFile(atPath: source.path, contents: nil))
        let handle = try FileHandle(forWritingTo: source)
        try handle.truncate(atOffset: 1_600_000_016)
        try handle.close()
        #expect(throws: (any Error).self) {
            try VivoOmicsFileSnapshot.fingerprint(source, maximumBytes: 1_600_000_000)
        }
        let clone = root.appendingPathComponent("snapshot.bin")
        let copied = try VivoOmicsFileSnapshot.fingerprint(source, copyTo: clone,
            maximumBytes: VivoIntegrationStorageLimits.maximumMembershipBytes)
        #expect(copied == (try VivoOmicsFileSnapshot.fingerprint(clone,
            maximumBytes: VivoIntegrationStorageLimits.maximumAnchorBytes)))
        #expect((try clone.resourceValues(forKeys: [.fileSizeKey])).fileSize == 1_600_000_016)
        let oversized = try FileHandle(forWritingTo: source)
        try oversized.truncate(atOffset: UInt64(VivoIntegrationStorageLimits.maximumMembershipBytes + 16))
        try oversized.close()
        #expect(throws: (any Error).self) {
            try VivoOmicsFileSnapshot.fingerprint(source,
                maximumBytes: VivoIntegrationStorageLimits.maximumMembershipBytes)
        }
    }

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
            #expect(!resident.benefitsFromBatchedAccess)
            #expect(file.benefitsFromBatchedAccess)
            for i in 0..<n {
                let row = (0..<d).map { j in sin(Double(i * (j + 1)) / 17) + Double(i % 3) / Double(j + 1) }
                try resident.setRow(i, row); try file.setRow(i, row)
            }
            let a = try VivoSingleCellIntegration.run(cells: cells, x: resident, samples: samples, options: options) { try VivoIntegrationMatrix(rows: $0, columns: $1) }
            let b = try VivoSingleCellIntegration.run(cells: cells, x: file, samples: samples, options: options, matrix: make)
            #expect(try a.materialize(cells: cells, options: options) == b.materialize(cells: cells, options: options))
            for (left, right) in [(a.scores, b.scores), (a.memberships, b.memberships), (a.assignmentScores, b.assignmentScores)] {
                for i in 0..<n { #expect(try left.row(i).map(\.bitPattern) == right.row(i).map(\.bitPattern)) }
            }
            #expect(a.objectives.map(\.bitPattern) == b.objectives.map(\.bitPattern))
            #expect(a.relativeImprovements.map(\.bitPattern) == b.relativeImprovements.map(\.bitPattern))
            for matrix in owned { try matrix.remove() }
            #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
        }
    }
    @Test func sortedBatchesPreserveLogicalBitsAndRejectBeforeWriting() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let n = VivoIntegrationMatrix.maximumBatchRows + 9
        for scratch: URL? in [nil, root] {
            let matrix = try VivoIntegrationMatrix(rows: n, columns: 3, scratch: scratch, windowBytes: 16_384)
            #expect(matrix.benefitsFromBatchedAccess == (scratch != nil))
            func initial(_ i: Int) -> [Double] { [Double(i) / 17, i % 2 == 0 ? -0.0 : 0.0, Double.leastNonzeroMagnitude * Double(i + 1)] }
            for i in 0..<n { try matrix.setRow(i, initial(i)) }
            let order = Array((0..<n).reversed())
            for start in stride(from: 0, to: n, by: VivoIntegrationMatrix.maximumBatchRows) {
                let indices = Array(order[start..<min(n, start + VivoIntegrationMatrix.maximumBatchRows)])
                var rows = try matrix.gatherRows(indices)
                for slot in indices.indices {
                    #expect(rows[slot].map(\.bitPattern) == initial(indices[slot]).map(\.bitPattern))
                    rows[slot][0] += 1
                }
                try matrix.scatterRows(indices, rows: rows)
            }
            for i in 0..<n {
                var expected = initial(i); expected[0] += 1
                #expect(try matrix.row(i).map(\.bitPattern) == expected.map(\.bitPattern))
            }
            let before = try matrix.materialize().map { $0.map(\.bitPattern) }
            #expect(try matrix.gatherRows([]).isEmpty)
            try matrix.scatterRows([], rows: [])
            for invalid in [[-1], [n], [2, 1, 2], Array(0...VivoIntegrationMatrix.maximumBatchRows)] {
                #expect(throws: (any Error).self) { try matrix.gatherRows(invalid) }
                #expect(throws: (any Error).self) { try matrix.scatterRows(invalid, rows: invalid.map { _ in [1, 2, 3] }) }
            }
            for invalid: [[Double]] in [[], [[1, 2, 3]], [[1, 2, 3], [4, 5]], [[1, 2, 3], [.nan, 5, 6]], [[1, 2, 3], [4, .infinity, 6]]] {
                #expect(throws: (any Error).self) { try matrix.scatterRows([0, 1], rows: invalid) }
            }
            #expect(try matrix.materialize().map { $0.map(\.bitPattern) } == before)
            try matrix.remove()
            #expect(throws: (any Error).self) { try matrix.gatherRows([]) }
            #expect(throws: (any Error).self) { try matrix.scatterRows([], rows: []) }
        }
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
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
