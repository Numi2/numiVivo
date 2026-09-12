import Foundation

/// Bounded donor/group means of per-cell log1p counts per million.
/// Supply every feature's positive counts in strictly increasing (row, feature)
/// order. Totals must cover that complete feature axis, before any projection.
/// Zero-count rows contribute zero and remain in the group denominator.
public struct VivoStreamedLogCPM {
    public enum Failure: Error { case invalidPlan, invalidRecord, incompleteRow, closed }
    public struct Result: Codable, Sendable, Equatable {
        public let method: String
        public let featureIDs: [String]
        public let groupIDs: [String]
        public let cellCounts: [Int]
        public let zeroCellCounts: [Int]
        /// Group-major, on the complete supplied feature axis.
        public let means: [[Double]]
    }
    private let features: [String], groups: [String]
    private let assignments: [Int], totals: [UInt64]
    private var sizes: [Int], zeros: [Int]
    private var sums: [Double], corrections: [Double]
    private var row = -1, feature = -1
    private var rowSum: UInt64 = 0
    private var closed = false

    public init(featureIDs: [String], groupIDs: [String],
                rowGroups: [Int], rowTotals: [UInt64]) throws {
        guard !featureIDs.isEmpty, featureIDs.count <= 1_000_000,
              !groupIDs.isEmpty, groupIDs.count <= 4096,
              featureIDs.count <= 16_000_000 / groupIDs.count,
              !rowGroups.isEmpty, rowGroups.count <= 10_000_000,
              rowGroups.count == rowTotals.count,
              Set(featureIDs).count == featureIDs.count,
              Set(groupIDs).count == groupIDs.count,
              featureIDs.allSatisfy({ !$0.isEmpty }),
              groupIDs.allSatisfy({ !$0.isEmpty }),
              rowGroups.allSatisfy({ groupIDs.indices.contains($0) }) else {
            throw Failure.invalidPlan
        }
        features = featureIDs; groups = groupIDs
        assignments = rowGroups; totals = rowTotals
        sizes = Array(repeating: 0, count: groupIDs.count)
        zeros = sizes
        for i in rowGroups.indices {
            sizes[rowGroups[i]] += 1
            if rowTotals[i] == 0 { zeros[rowGroups[i]] += 1 }
        }
        guard sizes.allSatisfy({ $0 > 0 }) else { throw Failure.invalidPlan }
        sums = Array(repeating: 0, count: featureIDs.count * groupIDs.count)
        corrections = sums
    }

    public mutating func add(row nextRow: Int, feature nextFeature: Int, count: UInt64) throws {
        guard !closed else { throw Failure.closed }
        do {
            guard assignments.indices.contains(nextRow), features.indices.contains(nextFeature),
                  count > 0, nextRow > row || (nextRow == row && nextFeature > feature) else {
                throw Failure.invalidRecord
            }
            if nextRow != row {
                if row >= 0, rowSum != totals[row] { throw Failure.incompleteRow }
                for skipped in (row + 1)..<nextRow {
                    guard totals[skipped] == 0 else { throw Failure.incompleteRow }
                }
                row = nextRow; feature = -1; rowSum = 0
            }
            let (total, overflow) = rowSum.addingReportingOverflow(count)
            guard !overflow, total <= totals[row], totals[row] > 0 else { throw Failure.invalidRecord }
            rowSum = total; feature = nextFeature
            let value = log1p((Double(count) / Double(totals[row])) * 1_000_000)
            let index = assignments[row] * features.count + nextFeature
            let adjusted = value - corrections[index]
            let updated = sums[index] + adjusted
            corrections[index] = (updated - sums[index]) - adjusted
            sums[index] = updated
        } catch {
            closed = true
            throw error
        }
    }

    public mutating func finish() throws -> Result {
        guard !closed else { throw Failure.closed }
        closed = true
        if row >= 0, rowSum != totals[row] { throw Failure.incompleteRow }
        for skipped in (row + 1)..<totals.count {
            guard totals[skipped] == 0 else { throw Failure.incompleteRow }
        }
        let means = groups.indices.map { group in
            (0..<features.count).map { sums[group * features.count + $0] / Double(sizes[group]) }
        }
        return Result(method: "mean-per-cell-log1p-cpm-full-axis-v1",
                      featureIDs: features, groupIDs: groups, cellCounts: sizes,
                      zeroCellCounts: zeros, means: means)
    }
}


extension VivoStreamedLogCPM {
    private static func chunkPool<T>(_ body: () throws -> T) rethrows -> T {
        #if canImport(ObjectiveC)
        return try autoreleasepool(invoking: body)
        #else
        return try body()
        #endif
    }
    /// Consume canonical little-endian row-u32/feature-u32/count-u64 records.
    /// The caller owns the input. Incomplete or invalid streams return no result.
    public static func consume(featureIDs: [String], groupIDs: [String],
                               rowGroups: [Int], rowTotals: [UInt64],
                               read: () throws -> Data) throws -> Result {
        var accumulator = try Self(featureIDs: featureIDs, groupIDs: groupIDs,
                                   rowGroups: rowGroups, rowTotals: rowTotals)
        var buffer = Data()
        while true {
            try Task.checkCancellation()
            let ended = try chunkPool { () throws -> Bool in
                let chunk = try read()
                guard chunk.count <= 1_048_576 else { throw Failure.invalidRecord }
                if chunk.isEmpty { return true }
                buffer.append(chunk)
                let complete = buffer.count / 16 * 16
                try buffer.withUnsafeBytes { bytes in
                    for offset in stride(from: 0, to: complete, by: 16) {
                        if offset % 16_384 == 0 { try Task.checkCancellation() }
                        try accumulator.add(
                            row: Int(UInt32(littleEndian: bytes.loadUnaligned(fromByteOffset: offset, as: UInt32.self))),
                            feature: Int(UInt32(littleEndian: bytes.loadUnaligned(fromByteOffset: offset + 4, as: UInt32.self))),
                            count: UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset + 8, as: UInt64.self)))
                    }
                }
                buffer = Data(buffer.suffix(buffer.count - complete))
                return false
            }
            if ended { break }
        }
        guard buffer.isEmpty else { throw Failure.invalidRecord }
        try Task.checkCancellation()
        return try accumulator.finish()
    }
}
