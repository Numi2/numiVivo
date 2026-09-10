import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// The same frozen transform serves resident reference mapping and windowed
/// query projection. Centers/loadings come exclusively from the training fit.
struct VivoFrozenPCAProjection {
    let loadings: [[Double]]
    let initial: [Double]
    let target: Double
    var components: Int { initial.count }
    init(centers: [Double], loadings: [[Double]], target: Double) throws {
        guard !centers.isEmpty, centers.count == loadings.count, let d = loadings.first?.count,
              (1...64).contains(d), target.isFinite, target > 0,
              centers.allSatisfy(\.isFinite), loadings.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }) else {
            throw VivoOmicsError.invalid("frozen PCA model dimensions or finite values")
        }
        self.loadings = loadings; self.target = target
        var shift = [Double](repeating: 0, count: d)
        for j in centers.indices { for c in 0..<d { shift[c] += centers[j] * loadings[j][c] } }
        guard shift.allSatisfy(\.isFinite) else { throw VivoOmicsError.invalid("frozen PCA shift is nonfinite") }
        initial = shift.map { -$0 }
    }
    func add(count: UInt64, total: UInt64, selected: Int, to row: UnsafeMutableBufferPointer<Double>) throws {
        guard total > 0, count > 0, count <= total, selected >= 0, selected < loadings.count, row.count == components else {
            throw VivoOmicsError.invalid("frozen PCA projection coordinate or counts")
        }
        let log = log1p(Double(count) / Double(total) * target)
        guard log.isFinite, log > 0 else { throw VivoOmicsError.invalid("frozen PCA normalization overflow/underflow") }
        for c in 0..<components { row[c] += log * loadings[selected][c] }
    }
}

/// Private row-major Float64 scratch with at most 16 MiB mapped at a time.
/// Row-aligned windows also start at page boundaries. No cells-by-PC array is
/// retained on the heap; metadata and per-cell quality remain resident.
final class VivoWindowedPCAScores {
    private let fd: Int32
    private let rows: Int
    private let components: Int
    private let windowRows: Int
    private var address: UnsafeMutableRawPointer?
    private var firstRow = -1
    private var length = 0
    let byteCount: Int
    init(_ url: URL, rows: Int, initial: [Double]) throws {
        guard (1...VivoPCAStorageLimits.maximumRows).contains(rows), (1...64).contains(initial.count), initial.allSatisfy(\.isFinite) else {
            throw VivoOmicsError.limit("windowed PCA score dimensions")
        }
        self.rows = rows; components = initial.count
        byteCount = rows * components * 8
        // macOS arm64 uses 16 KiB pages; align to the actual host page size.
        let page = Int(sysconf(Int32(_SC_PAGESIZE)))
        guard page > 0, page <= 1_048_576, page % 8 == 0 else { throw VivoOmicsError.invalid("PCA score page size") }
        let alignmentRows = page / 8
        windowRows = max(1, 16 * 1_024 * 1_024 / (components * 8 * alignmentRows)) * alignmentRows
        fd = open(url.path, O_RDWR | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
        guard fd >= 0 else { throw VivoOmicsError.invalid("cannot create PCA score scratch") }
        guard ftruncate(fd, off_t(byteCount)) == 0 else {
            _ = close(fd); throw VivoOmicsError.limit("cannot size PCA score scratch")
        }
        for row in 0..<rows { try withRow(row) { values in
            for c in 0..<components { values[c] = initial[c] }
        } }
    }
    deinit { if let address { _ = munmap(address, length) }; _ = close(fd) }
    func withRow<T>(_ row: Int, _ body: (UnsafeMutableBufferPointer<Double>) throws -> T) throws -> T {
        guard row >= 0, row < rows else { throw VivoOmicsError.invalid("PCA score row") }
        let next = row / windowRows * windowRows
        if next != firstRow {
            try Task.checkCancellation()
            if let address { _ = munmap(address, length); self.address = nil }
            length = min(windowRows, rows - next) * components * 8
            let mapped = mmap(nil, length, PROT_READ | PROT_WRITE, MAP_SHARED, fd, off_t(next * components * 8))
            guard mapped != MAP_FAILED, let mapped else { throw VivoOmicsError.limit("cannot map PCA score window") }
            address = mapped; firstRow = next
        }
        return try body(.init(start: address!.assumingMemoryBound(to: Double.self).advanced(by: (row - firstRow) * components), count: components))
    }
    func write(to url: URL) throws -> VivoFingerprint {
        let writer = try VivoCountRecordWriter(url)
        for row in 0..<rows { try withRow(row) { values in
            for c in 0..<components {
                guard values[c].isFinite else { throw VivoOmicsError.invalid("PCA query score is nonfinite") }
                try writer.append(row: row, feature: c, bits: values[c].bitPattern)
            }
        } }
        return try writer.finish()
    }
}
