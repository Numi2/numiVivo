import Foundation
import Dispatch

public struct VivoPCANeighborExecution: Codable, Sendable, Equatable {
    public var workers: Int = 4
    public var queryBlockRows: Int = 128
    public var candidateBlockRows: Int = 2_048
    public init() {}
    private enum CodingKeys: String, CodingKey { case workers, queryBlockRows, candidateBlockRows }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["workers", "queryBlockRows", "candidateBlockRows"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        workers = try c.decodeIfPresent(Int.self, forKey: .workers) ?? 4
        queryBlockRows = try c.decodeIfPresent(Int.self, forKey: .queryBlockRows) ?? 128
        candidateBlockRows = try c.decodeIfPresent(Int.self, forKey: .candidateBlockRows) ?? 2_048
    }
    public func validate() throws {
        guard (1...16).contains(workers), (1...512).contains(queryBlockRows), (1...8_192).contains(candidateBlockRows) else {
            throw VivoOmicsError.invalid("PCA neighbor worker or tile bounds")
        }
    }
}

/// One immutable score-file reader per worker. The common 16 MiB mapping window
/// and explicit tiles bound score working memory independently of cohort size.
final class VivoPCAScoreReader {
    private let records: VivoWindowedCountRecords
    let rows: Int
    let columns: Int
    init(_ url: URL, rows: Int, columns: Int) throws {
        guard (1...1_000_000).contains(rows), (1...64).contains(columns) else { throw VivoOmicsError.invalid("PCA score axes") }
        self.rows = rows; self.columns = columns
        records = try VivoWindowedCountRecords(url, entries: rows * columns)
    }
    func readRows(_ range: Range<Int>) throws -> [Double] {
        guard range.lowerBound >= 0, range.upperBound <= rows, range.count <= 8_192 else { throw VivoOmicsError.limit("PCA score read tile") }
        var values: [Double] = []; values.reserveCapacity(range.count * columns)
        for row in range {
            for c in 0..<columns {
                let r = try records.record(row * columns + c), value = Double(bitPattern: r.bits)
                guard r.row == row, r.feature == c, value.isFinite else { throw VivoOmicsError.invalid("PCA score coordinates or values") }
                values.append(value)
            }
        }
        return values
    }
}

enum VivoWindowedPCANeighbors {
    private struct Block: Sendable { let indices: [Int]; let distances: [Double] }
    /// Slots are written under the lock and read only after concurrentPerform
    /// joins. No mutable score reader, heap or numerical accumulator is shared.
    private final class Results: @unchecked Sendable {
        private let lock = NSLock()
        private var values: [Result<Block, Error>?]
        init(_ count: Int) { values = Array(repeating: nil, count: count) }
        func set(_ value: Result<Block, Error>, at index: Int) { lock.lock(); defer { lock.unlock() }; values[index] = value }
        func get(_ index: Int) throws -> Block {
            lock.lock(); defer { lock.unlock() }
            guard let value = values[index] else { throw VivoOmicsError.invalid("PCA neighbor worker did not complete") }
            return try value.get()
        }
    }
    private static func block(source: URL, range: Range<Int>, rows: Int, dimensions d: Int, neighbors k: Int, candidateRows: Int) throws -> Block {
        let reader = try VivoPCAScoreReader(source, rows: rows, columns: d)
        let query = try reader.readRows(range)
        var heaps = Array(repeating: [VivoSingleCellNeighbors.Neighbor](), count: range.count)
        for start in stride(from: 0, to: rows, by: candidateRows) {
            let end = min(rows, start + candidateRows), candidates = try reader.readRows(start..<end)
            try query.withUnsafeBufferPointer { q in try candidates.withUnsafeBufferPointer { r in
                for local in 0..<range.count {
                    let row = range.lowerBound + local, qi = local * d
                    for other in start..<end where other != row {
                        let ri = (other - start) * d
                        var distance = 0.0
                        for c in 0..<d { let delta = q[qi+c] - r[ri+c]; distance += delta * delta }
                        guard distance.isFinite else { throw VivoOmicsError.invalid("nonfinite neighbor distance") }
                        VivoSingleCellNeighbors.retain(.init(index: other, squaredDistance: distance), in: &heaps[local], capacity: k-1)
                    }
                }
            } }
        }
        var indices: [Int] = [], distances: [Double] = []
        indices.reserveCapacity(range.count*k); distances.reserveCapacity(range.count*k)
        for (local, row) in range.enumerated() {
            indices.append(row); distances.append(0)
            for neighbor in heaps[local].sorted(by: { $0.precedes($1) }) {
                indices.append(neighbor.index); distances.append(sqrt(neighbor.squaredDistance))
            }
        }
        return .init(indices: indices, distances: distances)
    }
    static func stream(source: URL, rows n: Int, dimensions: Int, options: VivoSingleCellNeighborOptions,
                       execution: VivoPCANeighborExecution, sink: (Int, [Int], [Double]) throws -> Void) throws {
        try options.validate(); try execution.validate()
        let k = options.neighbors
        guard n >= k, n <= 1_000_000, (1...64).contains(dimensions) else {
            throw VivoOmicsError.limit("PCA neighbor axes, graph-entry bound or unsupported integrated representation")
        }
        let pairs = n*(n-1)/2
        guard pairs <= options.maximumDistancePairs else {
            throw VivoOmicsError.limit("exact neighbor distance-pair budget; use the PCA bundle HNSW route for approximate search")
        }
        // Cancellation is checked between bounded parallel batches on the caller
        // task. Dispatch workers own independent data and do not inherit Tasks.
        let batch = execution.workers * execution.queryBlockRows
        for start in stride(from: 0, to: n, by: batch) {
            try Task.checkCancellation()
            let first = start, active = min(execution.workers, (n-start+execution.queryBlockRows-1)/execution.queryBlockRows)
            let results = Results(active)
            DispatchQueue.concurrentPerform(iterations: active) { worker in
                let lower = first + worker * execution.queryBlockRows, upper = min(n, lower + execution.queryBlockRows)
                results.set(Result { try block(source: source, range: lower..<upper, rows: n, dimensions: dimensions,
                    neighbors: k, candidateRows: execution.candidateBlockRows) }, at: worker)
            }
            for worker in 0..<active {
                let result = try results.get(worker)
                for local in 0..<(result.indices.count/k) {
                    try sink(first + worker*execution.queryBlockRows + local,
                             Array(result.indices[(local*k)..<((local+1)*k)]), Array(result.distances[(local*k)..<((local+1)*k)]))
                }
            }
        }
        try Task.checkCancellation()
    }
    static func run(source: URL, cells: [VivoOmicsCellIdentity], dimensions: Int, options: VivoSingleCellNeighborOptions,
                    execution: VivoPCANeighborExecution) throws -> VivoSingleCellNeighborGraph {
        try options.validate()
        guard cells.count <= 4_000_000/options.neighbors else { throw VivoOmicsError.limit("resident exact graph-entry bound") }
        var indices: [Int] = [], distances: [Double] = []
        try stream(source: source, rows: cells.count, dimensions: dimensions, options: options, execution: execution) { _, i, d in
            indices.append(contentsOf: i); distances.append(contentsOf: d)
        }
        return try VivoSingleCellNeighbors.finish(indices: indices, distances: distances, cells: cells, dimensions: dimensions, options: options, distancePairs: cells.count*(cells.count-1)/2)
    }
}
