import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Small latent axes only. File rows are padded to a power of two, so no row
/// crosses a mapping boundary. Scratch is private, mutable and never a receipt.
final class VivoIntegrationMatrix {
    static let maximumWindowBytes = 64 * 1_024 * 1_024
    static let maximumBatchRows = 8_192
    let rows: Int
    let columns: Int
    let strideBytes: Int
    let fileBytes: Int
    let windowBytes: Int
    private var values: [Double]?
    private var fd: Int32 = -1
    private var path: URL?
    private var address: UnsafeMutableRawPointer?
    private var offset = -1
    private var length = 0
    private var removed = false

    /// A single mapped window already makes shuffled access inexpensive.
    var benefitsFromBatchedAccess: Bool { fd >= 0 && fileBytes > windowBytes }

    init(rows: Int, columns: Int, scratch: URL? = nil, windowBytes: Int = maximumWindowBytes) throws {
        guard (1...VivoPCAStorageLimits.maximumRows).contains(rows), (1...128).contains(columns),
              windowBytes >= Int(getpagesize()), windowBytes <= Self.maximumWindowBytes,
              windowBytes.nonzeroBitCount == 1 else { throw VivoOmicsError.invalid("integration matrix axes or window") }
        var stride = 8
        while stride < columns * 8 { stride *= 2 }
        self.rows = rows; self.columns = columns; strideBytes = stride
        fileBytes = rows * stride; self.windowBytes = windowBytes
        if let scratch {
            let target = scratch.appendingPathComponent(".integration-matrix-" + UUID().uuidString)
            let descriptor = open(target.path, O_RDWR | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
            guard descriptor >= 0 else { throw VivoOmicsError.invalid("integration matrix creation") }
            guard ftruncate(descriptor, off_t(fileBytes)) == 0 else {
                _ = close(descriptor); try? FileManager.default.removeItem(at: target)
                throw VivoOmicsError.limit("integration matrix file size")
            }
            fd = descriptor; path = target
        } else { values = [Double](repeating: 0, count: rows * columns) }
    }
    deinit {
        if let address { _ = munmap(address, length) }
        if fd >= 0 { _ = close(fd) }
        if let path { try? FileManager.default.removeItem(at: path) }
    }
    func remove() throws {
        guard !removed else { return }
        if let address {
            guard munmap(address, length) == 0 else { throw VivoOmicsError.invalid("integration matrix unmap") }
            self.address = nil
        }
        if fd >= 0 {
            let descriptor = fd; fd = -1
            guard close(descriptor) == 0 else { throw VivoOmicsError.invalid("integration matrix close") }
        }
        if let path { try FileManager.default.removeItem(at: path); self.path = nil }
        values = nil; removed = true
    }
    private func check(_ row: Int) throws {
        guard !removed, row >= 0, row < rows else { throw VivoOmicsError.invalid("integration matrix row or closed state") }
    }
    private func pointer(_ row: Int) throws -> UnsafeMutableRawPointer {
        let byte = row * strideBytes, next = byte / windowBytes * windowBytes
        if next != offset {
            try Task.checkCancellation()
            if let address {
                guard munmap(address, length) == 0 else { throw VivoOmicsError.invalid("integration matrix unmap") }
                self.address = nil; offset = -1
            }
            length = min(windowBytes, fileBytes - next)
            let mapped = mmap(nil, length, PROT_READ | PROT_WRITE, MAP_SHARED, fd, off_t(next))
            guard mapped != MAP_FAILED, let mapped else { throw VivoOmicsError.limit("integration matrix mapping") }
            address = mapped; offset = next
        }
        return address!.advanced(by: byte - offset)
    }
    func row(_ index: Int) throws -> [Double] {
        try check(index)
        if let values { return Array(values[(index * columns)..<((index + 1) * columns)]) }
        let p = try pointer(index)
        return (0..<columns).map { Double(bitPattern: UInt64(littleEndian: p.load(fromByteOffset: $0 * 8, as: UInt64.self))) }
    }
    func value(_ row: Int, _ column: Int) throws -> Double {
        try check(row)
        guard column >= 0, column < columns else { throw VivoOmicsError.invalid("integration matrix column") }
        if let values { return values[row * columns + column] }
        return Double(bitPattern: UInt64(littleEndian: try pointer(row).load(fromByteOffset: column * 8, as: UInt64.self)))
    }
    func setRow(_ index: Int, _ row: [Double]) throws {
        try check(index)
        guard row.count == columns, row.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("integration matrix row values") }
        if values != nil {
            for j in 0..<columns { values![index * columns + j] = row[j] }
        } else {
            let p = try pointer(index)
            for j in 0..<columns { p.storeBytes(of: row[j].bitPattern.littleEndian, toByteOffset: j * 8, as: UInt64.self) }
        }
    }
    private func orderedSlots(_ indices: [Int]) throws -> [Int] {
        guard !removed, indices.count <= Self.maximumBatchRows else {
            throw VivoOmicsError.invalid("integration matrix batch size or closed state")
        }
        for index in indices { try check(index) }
        let slots = indices.indices.sorted { indices[$0] < indices[$1] }
        for (left, right) in zip(slots, slots.dropFirst()) where indices[left] == indices[right] {
            throw VivoOmicsError.invalid("integration matrix duplicate batch row")
        }
        return slots
    }
    /// Read in physical order, return in caller order. The fixed row limit bounds
    /// the value buffer independently of the cohort size or shuffled block size.
    func gatherRows(_ indices: [Int]) throws -> [[Double]] {
        let slots = try orderedSlots(indices)
        try Task.checkCancellation()
        var result = [[Double]](repeating: [], count: indices.count)
        for slot in slots { result[slot] = try row(indices[slot]) }
        return result
    }
    /// Validate the entire batch before any write. Only IO is reordered; each
    /// result row remains paired with its original logical index.
    func scatterRows(_ indices: [Int], rows: [[Double]]) throws {
        let slots = try orderedSlots(indices)
        guard rows.count == indices.count,
              rows.allSatisfy({ $0.count == columns && $0.allSatisfy(\.isFinite) }) else {
            throw VivoOmicsError.invalid("integration matrix batch values")
        }
        try Task.checkCancellation()
        for slot in slots { try setRow(indices[slot], rows[slot]) }
    }
    func copy(from source: VivoIntegrationMatrix, normalize: Bool = false) throws {
        guard rows == source.rows, columns == source.columns else { throw VivoOmicsError.invalid("integration matrix copy axes") }
        for i in 0..<rows {
            if i % 2_048 == 0 { try Task.checkCancellation() }
            var row = try source.row(i)
            if normalize {
                let norm = sqrt(row.reduce(0) { $0 + $1 * $1 })
                if norm > 0 { row = row.map { $0 / norm } }
            }
            try setRow(i, row)
        }
    }
    func materialize() throws -> [[Double]] { try (0..<rows).map { try row($0) } }
    func writeRecords(to url: URL) throws -> VivoFingerprint {
        let writer = try VivoCountRecordWriter(url)
        for i in 0..<rows {
            let row = try row(i)
            for j in 0..<columns { try writer.append(row: i, feature: j, bits: row[j].bitPattern) }
        }
        return try writer.finish()
    }
    /// Consume private scratch only after the complete output has been written
    /// and fingerprinted. The owning publication cleans up on any failure.
    func writeRecordsAndRemove(to url: URL) throws -> VivoFingerprint {
        let fingerprint = try writeRecords(to: url)
        try remove()
        return fingerprint
    }
}
