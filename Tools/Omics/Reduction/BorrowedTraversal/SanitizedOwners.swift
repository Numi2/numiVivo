import Foundation
import Darwin
enum VivoOmicsError: Error { case invalid(String), limit(String) }
enum VivoH5ADCountStore { static let maximumEntries = 2_000_000_000 }
final class VivoWindowedCountRecords {
    static let windowBytes = 16 * 1_024 * 1_024
    private let fd: Int32
    let count: Int
    private var address: UnsafeMutableRawPointer?
    private var offset = -1
    private var length = 0
    init(_ url: URL, entries: Int) throws {
        guard entries >= 0, entries <= VivoH5ADCountStore.maximumEntries else { throw VivoOmicsError.limit("mapped count entries") }
        fd = open(url.path, O_RDONLY | O_NOFOLLOW | O_NONBLOCK)
        guard fd >= 0 else { throw VivoOmicsError.invalid("cannot open count records") }
        var info = stat()
        guard fstat(fd, &info) == 0, (info.st_mode & mode_t(S_IFMT)) == mode_t(S_IFREG), info.st_size == entries * 16 else {
            _ = close(fd); throw VivoOmicsError.invalid("count record size or type")
        }
        count = entries
    }
    deinit { if let address { _ = munmap(address, length) }; _ = close(fd) }
    /// Visit immutable records in source order with one mapping check per window.
    /// The borrowed pointer never escapes the synchronous body.
    @inline(__always)
    func forEachRecord(_ body: (Int, Int, UInt64) throws -> Void) throws {
        var index = 0
        while index < count {
            try Task.checkCancellation()
            _ = try record(index)
            let pointer = address!
            let recordsInWindow = length / 16
            for localIndex in 0..<recordsInWindow {
                let local = localIndex * 16
                let packed = UInt64(littleEndian: pointer.load(fromByteOffset: local, as: UInt64.self))
                let bits = UInt64(littleEndian: pointer.load(fromByteOffset: local + 8, as: UInt64.self))
                try body(Int(packed & 0xffff_ffff), Int(packed >> 32), bits)
            }
            index += recordsInWindow
        }
    }
    func record(_ index: Int) throws -> (row: Int, feature: Int, bits: UInt64) {
        guard index >= 0, index < count else { throw VivoOmicsError.invalid("count record index") }
        let byte = index * 16, next = (byte / Self.windowBytes) * Self.windowBytes
        if next != offset {
            try Task.checkCancellation()
            if let address { _ = munmap(address, length); self.address = nil }
            length = min(Self.windowBytes, count * 16 - next)
            let mapped = mmap(nil, length, PROT_READ, MAP_PRIVATE, fd, off_t(next))
            guard mapped != MAP_FAILED, let mapped else { throw VivoOmicsError.limit("cannot map count window") }
            address = mapped; offset = next
        }
        let p = address!, local = byte - offset
        let packed = UInt64(littleEndian: p.load(fromByteOffset: local, as: UInt64.self))
        return (Int(packed & 0xffff_ffff), Int(packed >> 32), UInt64(littleEndian: p.load(fromByteOffset: local + 8, as: UInt64.self)))
    }
}

final class VivoMappedReductionEntries {
    private let records: VivoWindowedCountRecords
    private let rowCount: Int
    private let columnCount: Int
    let byteCount: Int
    let count: Int
    private(set) var visits = 0
    init(url: URL, expectedEntries: Int, rows: Int, columns: Int) throws {
        guard expectedEntries > 0, expectedEntries <= VivoH5ADCountStore.maximumEntries, rows > 0, columns > 0 else {
            throw VivoOmicsError.invalid("reduction cache dimensions")
        }
        records = try VivoWindowedCountRecords(url, entries: expectedEntries)
        rowCount = rows; columnCount = columns
        byteCount = expectedEntries * 16; count = expectedEntries
        // The owner retains an immutable private scratch file. Validate the
        // complete stream once, using the same bounded window as later passes.
        for i in 0..<count {
            let record = try records.record(i)
            let value = Double(bitPattern: record.bits)
            guard record.row < rows, record.feature < columns, value.isFinite, value > 0 else {
                throw VivoOmicsError.invalid("reduction cache record")
            }
        }
    }
    func project(_ vector: [Double], shift: Double, rows: Int) throws -> [Double] {
        guard vector.count == columnCount, rows == rowCount else {
            throw VivoOmicsError.invalid("reduction projection dimensions")
        }
        var result = [Double](repeating: -shift, count: rows)
        try result.withUnsafeMutableBufferPointer { output in
            try vector.withUnsafeBufferPointer { input in
                try records.forEachRecord { row, feature, bits in
                    output[row] += Double(bitPattern: bits) * input[feature]
                }
            }
        }
        visits += count; return result
    }
    func transpose(_ vector: [Double], initial: [Double]) throws -> [Double] {
        guard vector.count == rowCount, initial.count == columnCount else {
            throw VivoOmicsError.invalid("reduction transpose dimensions")
        }
        var result = initial
        try result.withUnsafeMutableBufferPointer { output in
            try vector.withUnsafeBufferPointer { input in
                try records.forEachRecord { row, feature, bits in
                    output[feature] += Double(bitPattern: bits) * input[row]
                }
            }
        }
        visits += count; return result
    }
}
