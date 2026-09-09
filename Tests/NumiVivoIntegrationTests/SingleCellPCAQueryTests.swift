import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellPCAQueryTests {
    @Test func frozenProjectionUsesTrainingCentersAndAllGeneTotals() throws {
        let model = try VivoFrozenPCAProjection(centers: [1, 2], loadings: [[0.5, -1], [2, 0.25]], target: 10_000)
        var score = model.initial
        #expect(score == [-4.5, 0.5])
        try score.withUnsafeMutableBufferPointer { try model.add(count: 3, total: 100, selected: 0, to: $0) }
        #expect(score == [-4.5 + log1p(300) * 0.5, 0.5 - log1p(300)])
        #expect(throws: (any Error).self) {
            try score.withUnsafeMutableBufferPointer { try model.add(count: 3, total: 0, selected: 0, to: $0) }
        }
        #expect(throws: (any Error).self) { try VivoFrozenPCAProjection(centers: [.nan], loadings: [[1]], target: 10_000) }
    }
    @Test func scoreWindowsPreserveEveryValueAcrossRewindsAndPageBoundaries() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let rows = 40_000, d = 63
        let initial = (0..<d).map { Double($0) * -0.25 }
        let scores = try VivoWindowedPCAScores(root.appendingPathComponent("scratch"), rows: rows, initial: initial)
        #expect(scores.byteCount > 16 * 1_024 * 1_024)
        for row in [rows - 1, 0, 32_768, 1, rows - 1, 0] {
            try scores.withRow(row) { values in for c in 0..<d { values[c] += Double(c + 1) } }
        }
        let result = root.appendingPathComponent("scores.bin")
        let hash = try scores.write(to: result)
        #expect(try VivoH5ADCountStore.fingerprint(result) == hash)
        let records = try VivoWindowedCountRecords(result, entries: rows * d)
        for i in 0..<records.count {
            let record = try records.record(i), row = i / d, c = i % d
            let repeats = row == 0 || row == rows - 1 ? 2 : (row == 1 || row == 32_768 ? 1 : 0)
            let expected = initial[c] + Double(repeats * (c + 1))
            #expect(record.row == row && record.feature == c && Double(bitPattern: record.bits) == expected)
        }
        #expect(throws: (any Error).self) { try scores.withRow(rows) { _ in } }
    }
}
