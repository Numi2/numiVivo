import Foundation
import CryptoKit

public enum VivoAccessibilityTFIDFError: Error {
    case invalidDimensions, invalidCounts, changedReplay, nonfiniteResult
}

public struct VivoAccessibilityTFIDFReport: Codable, Sendable {
    public let cells: Int
    public let features: Int
    public let records: Int
    public let cellTotals: [UInt64]
    public let featureTotals: [UInt64]
    public let scaleFactor: Double
}

/// Signac method-1 arithmetic on positive canonical accessibility counts:
/// log1p(scale * (count / cellTotal) * (numberOfCells / featureTotal)).
/// Two source scans, O(cells + features) state, and no dense matrix. Source
/// identity, assay units, axes and transactional output belong to the caller.
public enum VivoAccessibilityTFIDF {
    private final class Signature {
        var hash = SHA256(), buffer = Data()
        func add(_ row: Int, _ column: Int, _ count: UInt64) throws {
            var coordinate = (UInt64(row) | (UInt64(column) << 32)).littleEndian
            var value = count.littleEndian
            withUnsafeBytes(of: &coordinate) { buffer.append(contentsOf: $0) }
            withUnsafeBytes(of: &value) { buffer.append(contentsOf: $0) }
            if buffer.count >= 1_048_576 {
                try Task.checkCancellation()
                hash.update(data: buffer); buffer.removeAll(keepingCapacity: true)
            }
        }
        func finish() -> SHA256.Digest { hash.update(data: buffer); return hash.finalize() }
    }
    public static func normalize(cells: Int, features: Int, maximumRecords: Int,
                                 scaleFactor: Double = 10_000,
                                 scan: (_ accept: (Int, Int, UInt64) throws -> Void) throws -> Void,
                                 emit: @escaping (Int, Int, Double) throws -> Void) throws -> VivoAccessibilityTFIDFReport {
        guard cells > 0, cells <= 10_000_000, features > 0, features <= 1_000_000,
              maximumRecords >= 0, scaleFactor.isFinite, scaleFactor > 0 else {
            throw VivoAccessibilityTFIDFError.invalidDimensions
        }
        var cellTotals = [UInt64](repeating: 0, count: cells)
        var featureTotals = [UInt64](repeating: 0, count: features)
        var previous = [Int](repeating: -1, count: cells)
        var records = 0
        let firstSignature = Signature()
        func checkedSum(_ a: UInt64, _ b: UInt64) throws -> UInt64 {
            let x = a.addingReportingOverflow(b)
            guard !x.overflow else { throw VivoAccessibilityTFIDFError.invalidCounts }
            return x.partialValue
        }
        try scan { row, column, count in
            guard row >= 0, row < cells, column >= 0, column < features,
                  column > previous[row], count > 0, records < maximumRecords else {
                throw VivoAccessibilityTFIDFError.invalidCounts
            }
            previous[row] = column; records += 1
            cellTotals[row] = try checkedSum(cellTotals[row], count)
            featureTotals[column] = try checkedSum(featureTotals[column], count)
            try firstSignature.add(row, column, count)
        }
        try Task.checkCancellation()
        var replayCells = [UInt64](repeating: 0, count: cells)
        var replayFeatures = [UInt64](repeating: 0, count: features)
        previous = [Int](repeating: -1, count: cells)
        var replayRecords = 0
        let secondSignature = Signature()
        try scan { row, column, count in
            guard row >= 0, row < cells, column >= 0, column < features,
                  column > previous[row], count > 0, replayRecords < records,
                  cellTotals[row] > 0, featureTotals[column] > 0 else {
                throw VivoAccessibilityTFIDFError.invalidCounts
            }
            previous[row] = column; replayRecords += 1
            replayCells[row] = try checkedSum(replayCells[row], count)
            replayFeatures[column] = try checkedSum(replayFeatures[column], count)
            try secondSignature.add(row, column, count)
            let value = log1p((Double(count) / Double(cellTotals[row])) *
                             (Double(cells) / Double(featureTotals[column])) * scaleFactor)
            guard value.isFinite else { throw VivoAccessibilityTFIDFError.nonfiniteResult }
            try emit(row, column, value)
        }
        guard replayRecords == records, replayCells == cellTotals, replayFeatures == featureTotals,
              firstSignature.finish() == secondSignature.finish() else {
            throw VivoAccessibilityTFIDFError.changedReplay
        }
        try Task.checkCancellation()
        return .init(cells: cells, features: features, records: records,
                     cellTotals: cellTotals, featureTotals: featureTotals, scaleFactor: scaleFactor)
    }
}
