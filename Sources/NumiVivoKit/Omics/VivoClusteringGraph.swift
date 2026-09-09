import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

final class VivoClusteringWork {
    let maximumEdgeVisits: Int
    private(set) var edgeVisits = 0
    private(set) var rowReads = 0
    private(set) var aggregatedLevels = 0
    private(set) var maximumEdgeMapEntries = 0
    init(maximumEdgeVisits: Int = Int.max) { self.maximumEdgeVisits = maximumEdgeVisits }
    func row() throws { try Task.checkCancellation(); rowReads += 1 }
    func edge() throws {
        guard edgeVisits < maximumEdgeVisits else { throw VivoOmicsError.limit("clustering edge-visit budget") }
        edgeVisits += 1
    }
    func map(_ count: Int) { maximumEdgeMapEntries = max(maximumEdgeMapEntries, count) }
    func aggregate() { aggregatedLevels += 1 }
}

protocol VivoClusteringGraph: AnyObject {
    var count: Int { get }
    var work: VivoClusteringWork { get }
    func forEachEdge(in row: Int, _ body: (VivoSingleCellClustering.Edge) throws -> Void) throws
    func aggregate(labels: [Int], groups: Int, scratch: URL?) throws -> any VivoClusteringGraph
}

final class VivoResidentClusteringGraph: VivoClusteringGraph {
    let rows: [[VivoSingleCellClustering.Edge]]
    let work: VivoClusteringWork
    var count: Int { rows.count }
    init(_ rows: [[VivoSingleCellClustering.Edge]], work: VivoClusteringWork = .init()) { self.rows = rows; self.work = work }
    func forEachEdge(in row: Int, _ body: (VivoSingleCellClustering.Edge) throws -> Void) throws {
        try work.row()
        for edge in rows[row] { try work.edge(); try body(edge) }
    }
    func aggregate(labels: [Int], groups: Int, scratch: URL?) throws -> any VivoClusteringGraph {
        var result = [[Int: Double]](repeating: [:], count: groups)
        for i in rows.indices { try forEachEdge(in: i) { edge in result[labels[i]][labels[edge.column], default: 0] += edge.weight } }
        for row in result { work.map(row.count) }; work.aggregate()
        return VivoResidentClusteringGraph(result.map { row in row.keys.sorted().map { .init(column: $0, weight: row[$0]!) } }, work: work)
    }
}

/// Private immutable CSR snapshots. Offsets are cell-scale; edge records use the
/// shared 16 MiB reader. Aggregation preserves source-row/column summation order.
final class VivoFileClusteringGraph: VivoClusteringGraph {
    let count: Int
    let work: VivoClusteringWork
    private let offsets: [Int]
    private let records: VivoWindowedCountRecords
    private let allowSelf: Bool
    private let ownedDirectory: URL?
    init(root: URL, rows: Int, entries: Int, work: VivoClusteringWork, allowSelf: Bool = false, owned: Bool = false) throws {
        guard (1...1_000_000).contains(rows), (0...254_000_000).contains(entries) else { throw VivoOmicsError.limit("clustering CSR axes") }
        count = rows; self.work = work; self.allowSelf = allowSelf; ownedDirectory = owned ? root : nil
        let source = try VivoWindowedCountRecords(root.appendingPathComponent("offsets.bin"), entries: rows+1)
        var offsets: [Int] = []; offsets.reserveCapacity(rows+1)
        for i in 0...rows {
            let record = try source.record(i)
            guard record.row == i, record.feature == 0, record.bits <= UInt64(entries),
                  i == 0 ? record.bits == 0 : record.bits >= UInt64(offsets[i-1]) else { throw VivoOmicsError.invalid("clustering CSR offsets") }
            offsets.append(Int(record.bits))
        }
        guard offsets[rows] == entries else { throw VivoOmicsError.invalid("clustering CSR terminal offset") }
        self.offsets = offsets
        records = try VivoWindowedCountRecords(root.appendingPathComponent("edges.bin"), entries: entries)
    }
    deinit { if let ownedDirectory { try? FileManager.default.removeItem(at: ownedDirectory) } }
    func forEachEdge(in row: Int, _ body: (VivoSingleCellClustering.Edge) throws -> Void) throws {
        guard row >= 0, row < count else { throw VivoOmicsError.invalid("clustering CSR row") }
        try work.row(); var previous = -1
        for i in offsets[row]..<offsets[row+1] {
            try work.edge()
            let record = try records.record(i), weight = Double(bitPattern: record.bits)
            guard record.row == row, record.feature < count, record.feature > previous,
                  allowSelf || record.feature != row, weight.isFinite, weight > 0 else { throw VivoOmicsError.invalid("clustering CSR edge") }
            previous = record.feature; try body(.init(column: record.feature, weight: weight))
        }
    }
    func aggregate(labels: [Int], groups: Int, scratch: URL?) throws -> any VivoClusteringGraph {
        guard let scratch, labels.count == count, groups > 0, groups <= count,
              labels.allSatisfy({ $0 >= 0 && $0 < groups }) else { throw VivoOmicsError.invalid("clustering aggregation axes") }
        let root = scratch.appendingPathComponent(".clustering-level-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        var retained = false
        defer { if !retained { try? FileManager.default.removeItem(at: root) } }
        var starts = [Int](repeating: 0, count: groups+1)
        for i in 0..<count { try forEachEdge(in: i) { _ in starts[labels[i]+1] += 1 } }
        for i in 0..<groups { starts[i+1] += starts[i] }
        var cursor = Array(starts.dropLast())
        let bucket = root.appendingPathComponent("bucket.bin")
        let fd = open(bucket.path, O_RDWR | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
        guard fd >= 0 else { throw VivoOmicsError.invalid("clustering bucket creation") }
        defer { _ = close(fd) }
        guard ftruncate(fd, off_t(starts[groups]*16)) == 0 else { throw VivoOmicsError.limit("clustering bucket allocation") }
        for i in 0..<count { try forEachEdge(in: i) { edge in
            let row = labels[i], column = labels[edge.column]
            var encoded = ((UInt64(row) | UInt64(column)<<32).littleEndian, edge.weight.bitPattern.littleEndian)
            try withUnsafeBytes(of: &encoded) { bytes in
                var done = 0
                while done < 16 {
                    let wrote = pwrite(fd, bytes.baseAddress!.advanced(by: done), 16-done, off_t(cursor[row]*16+done))
                    if wrote < 0 && errno == EINTR { continue }
                    guard wrote > 0 else { throw VivoOmicsError.invalid("clustering bucket write") }
                    done += wrote
                }
            }
            cursor[row] += 1
        } }
        for i in 0..<groups where cursor[i] != starts[i+1] { throw VivoOmicsError.invalid("clustering bucket cursor") }
        let input = try VivoWindowedCountRecords(bucket, entries: starts[groups])
        let edgeWriter = try VivoCountRecordWriter(root.appendingPathComponent("edges.bin"))
        let offsetWriter = try VivoCountRecordWriter(root.appendingPathComponent("offsets.bin"))
        for row in 0..<groups {
            try Task.checkCancellation()
            try offsetWriter.append(row: row, feature: 0, bits: UInt64(edgeWriter.entries))
            var sums: [Int: Double] = [:]
            for i in starts[row]..<starts[row+1] {
                try work.edge()
                let record = try input.record(i)
                guard record.row == row, record.feature < groups else { throw VivoOmicsError.invalid("clustering bucket record") }
                sums[record.feature, default: 0] += Double(bitPattern: record.bits)
            }
            work.map(sums.count)
            for column in sums.keys.sorted() {
                let weight = sums[column]!
                guard weight.isFinite, weight > 0 else { throw VivoOmicsError.invalid("clustering aggregated weight") }
                try edgeWriter.append(row: row, feature: column, bits: weight.bitPattern)
            }
        }
        try offsetWriter.append(row: groups, feature: 0, bits: UInt64(edgeWriter.entries))
        _ = try edgeWriter.finish(); _ = try offsetWriter.finish()
        try FileManager.default.removeItem(at: bucket)
        let result = try VivoFileClusteringGraph(root: root, rows: groups, entries: edgeWriter.entries, work: work, allowSelf: true, owned: true)
        work.aggregate(); retained = true; return result
    }
}
