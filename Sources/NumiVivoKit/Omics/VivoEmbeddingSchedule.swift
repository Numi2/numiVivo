import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

protocol VivoEmbeddingSchedule: AnyObject {
    var count: Int { get }
    var discarded: Int { get }
    func visit(epoch: Int, _ body: (Int, Int, Int) throws -> Void) throws
}

private func vivoEmbeddingBudget(count: Int, bound: Double, options: VivoSingleCellEmbeddingOptions) throws {
    guard bound.isFinite, bound <= Double(options.maximumUpdates),
          Double(count) * Double(options.epochs) <= Double(options.maximumUpdates) else {
        throw VivoOmicsError.limit("UMAP update budget before optimization")
    }
}

final class VivoResidentEmbeddingSchedule: VivoEmbeddingSchedule {
    private struct Event { let head: Int; let tail: Int; let interval: Double; var positive: Double; var negative: Double }
    private var events: [Event] = []
    private let rate: Double
    let discarded: Int
    var count: Int { events.count }
    init(graph: VivoSingleCellNeighborGraph, options: VivoSingleCellEmbeddingOptions) throws {
        guard let maximum = graph.weights.max(), maximum > 0 else { throw VivoOmicsError.invalid("embedding needs a nonempty graph") }
        rate = Double(options.negativeSampleRate); var bound = 0.0
        for i in graph.cells.indices {
            for edge in graph.rowOffsets[i]..<graph.rowOffsets[i+1] where graph.weights[edge] >= maximum / Double(options.epochs) {
                let interval = maximum / graph.weights[edge]
                events.append(.init(head: i, tail: graph.columnIndices[edge], interval: interval, positive: interval, negative: interval / rate))
                bound += ceil(Double(options.epochs) / interval) * Double(1 + options.negativeSampleRate)
            }
        }
        discarded = graph.weights.count - events.count
        try vivoEmbeddingBudget(count: events.count, bound: bound, options: options)
    }
    func visit(epoch: Int, _ body: (Int, Int, Int) throws -> Void) throws {
        for i in events.indices where events[i].positive <= Double(epoch) {
            let event = events[i], negativeInterval = event.interval / rate
            let samples = max(0, Int((Double(epoch) - event.negative) / negativeInterval))
            events[i].positive += event.interval
            events[i].negative += Double(samples) * negativeInterval
            try body(event.head, event.tail, samples)
        }
    }
}

/// A private, mutable 32-byte schedule record is four little-endian UInt64 words:
/// packed head/tail u32, interval f64, next-positive f64, next-negative f64.
/// Windows are shared so updates survive unmapping between sequential epochs.
/// This is scratch state, not a restart/checkpoint artifact.
final class VivoFileEmbeddingSchedule: VivoEmbeddingSchedule {
    static let maximumWindowBytes = 16 * 1_024 * 1_024
    let count: Int
    let discarded: Int
    let fileBytes: Int
    let windowBytes: Int
    private let fd: Int32
    private let path: URL
    private let rate: Double
    private var closed = false
    private var address: UnsafeMutableRawPointer?
    private var offset = -1
    private var length = 0
    init(edges: URL, entries: Int, cells: Int, options: VivoSingleCellEmbeddingOptions, scratch: URL,
         windowBytes: Int = maximumWindowBytes) throws {
        try options.validate()
        guard entries > 0, entries <= 254_000_000, (1...1_000_000).contains(cells),
              windowBytes >= Int(getpagesize()), windowBytes <= Self.maximumWindowBytes,
              windowBytes % Int(getpagesize()) == 0 else { throw VivoOmicsError.invalid("embedding schedule axes or window") }
        let records = try VivoWindowedCountRecords(edges, entries: entries)
        var maximum = 0.0, previousRow = -1, previousColumn = -1
        for i in 0..<entries {
            let record = try records.record(i), weight = Double(bitPattern: record.bits)
            guard record.row < cells, record.feature < cells, record.row != record.feature,
                  record.row > previousRow || (record.row == previousRow && record.feature > previousColumn),
                  weight.isFinite, weight > 0 else { throw VivoOmicsError.invalid("embedding graph edge") }
            previousRow = record.row; previousColumn = record.feature; maximum = max(maximum, weight)
        }
        let path = scratch.appendingPathComponent(".embedding-schedule-" + UUID().uuidString + ".bin")
        let fd = open(path.path, O_RDWR | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
        guard fd >= 0 else { throw VivoOmicsError.invalid("embedding schedule creation") }
        var retained = false
        defer { if !retained { _ = close(fd); try? FileManager.default.removeItem(at: path) } }
        var buffer: [UInt64] = []; buffer.reserveCapacity(131_072)
        func flush() throws {
            try buffer.withUnsafeBytes { bytes in
                var done = 0
                while done < bytes.count {
                    let written = write(fd, bytes.baseAddress!.advanced(by: done), bytes.count - done)
                    if written < 0 && errno == EINTR { continue }
                    guard written > 0 else { throw VivoOmicsError.invalid("embedding schedule write") }
                    done += written
                }
            }
            buffer.removeAll(keepingCapacity: true)
        }
        let rate = Double(options.negativeSampleRate); var count = 0, bound = 0.0
        for i in 0..<entries {
            let record = try records.record(i), weight = Double(bitPattern: record.bits)
            if weight < maximum / Double(options.epochs) { continue }
            let interval = maximum / weight
            count += 1; bound += ceil(Double(options.epochs) / interval) * Double(1 + options.negativeSampleRate)
            try vivoEmbeddingBudget(count: count, bound: bound, options: options)
            buffer.append((UInt64(record.row) | UInt64(record.feature) << 32).littleEndian)
            buffer.append(interval.bitPattern.littleEndian); buffer.append(interval.bitPattern.littleEndian)
            buffer.append((interval / rate).bitPattern.littleEndian)
            if buffer.count == 131_072 { try flush() }
        }
        try flush()
        var info = stat()
        guard fstat(fd, &info) == 0, info.st_size == count * 32 else { throw VivoOmicsError.invalid("embedding schedule size") }
        self.fd = fd; self.path = path; self.rate = rate; self.windowBytes = windowBytes
        self.count = count; discarded = entries - count; fileBytes = count * 32; retained = true
    }
    deinit {
        if let address { _ = munmap(address, length) }
        if !closed { _ = close(fd) }; try? FileManager.default.removeItem(at: path)
    }
    func remove() throws {
        if let address {
            guard munmap(address, length) == 0 else { throw VivoOmicsError.invalid("embedding schedule unmap") }
            self.address = nil
        }
        if !closed { closed = true; guard close(fd) == 0 else { throw VivoOmicsError.invalid("embedding schedule close") } }
        try FileManager.default.removeItem(at: path)
    }
    func visit(epoch: Int, _ body: (Int, Int, Int) throws -> Void) throws {
        guard !closed else { throw VivoOmicsError.invalid("embedding schedule already closed") }
        for i in 0..<count {
            let byte = i * 32, next = (byte / windowBytes) * windowBytes
            if next != offset {
                try Task.checkCancellation()
                if let address { _ = munmap(address, length); self.address = nil }
                length = min(windowBytes, fileBytes - next)
                let mapped = mmap(nil, length, PROT_READ | PROT_WRITE, MAP_SHARED, fd, off_t(next))
                guard mapped != MAP_FAILED, let mapped else { throw VivoOmicsError.limit("embedding schedule mapping") }
                address = mapped; offset = next
            }
            let p = address!, local = byte - offset
            let positive = Double(bitPattern: UInt64(littleEndian: p.load(fromByteOffset: local + 16, as: UInt64.self)))
            if positive > Double(epoch) { continue }
            let packed = UInt64(littleEndian: p.load(fromByteOffset: local, as: UInt64.self))
            let interval = Double(bitPattern: UInt64(littleEndian: p.load(fromByteOffset: local + 8, as: UInt64.self)))
            let negative = Double(bitPattern: UInt64(littleEndian: p.load(fromByteOffset: local + 24, as: UInt64.self)))
            let negativeInterval = interval / rate, samples = max(0, Int((Double(epoch) - negative) / negativeInterval))
            p.storeBytes(of: (positive + interval).bitPattern.littleEndian, toByteOffset: local + 16, as: UInt64.self)
            p.storeBytes(of: (negative + Double(samples) * negativeInterval).bitPattern.littleEndian, toByteOffset: local + 24, as: UInt64.self)
            try body(Int(packed & 0xffff_ffff), Int(packed >> 32), samples)
        }
    }
}
