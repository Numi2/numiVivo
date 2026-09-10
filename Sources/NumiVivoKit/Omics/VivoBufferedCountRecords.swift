import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Small immutable record pages for shuffled CSR row access. Sequential scans
/// can retain the mapped reader; random rows do not remap a large virtual window.
final class VivoBufferedCountRecords {
    static let bufferBytes = 4_096
    private let fd: Int32
    let count: Int
    private let buffer: UnsafeMutablePointer<UInt64>
    private var firstRecord = -1
    private var bufferedRecords = 0
    private(set) var bufferLoads = 0
    private(set) var loadedBytes = 0

    init(_ url: URL, entries: Int) throws {
        guard entries >= 0, entries <= VivoH5ADCountStore.maximumEntries else {
            throw VivoOmicsError.limit("buffered count entries")
        }
        let descriptor = open(url.path, O_RDONLY | O_NOFOLLOW | O_NONBLOCK)
        guard descriptor >= 0 else { throw VivoOmicsError.invalid("cannot open buffered count records") }
        var info = stat()
        guard fstat(descriptor, &info) == 0,
              (info.st_mode & mode_t(S_IFMT)) == mode_t(S_IFREG), info.st_size == entries * 16 else {
            _ = close(descriptor)
            throw VivoOmicsError.invalid("buffered count record size or type")
        }
        fd = descriptor; count = entries
        buffer = .allocate(capacity: Self.bufferBytes / 8)
    }
    deinit { buffer.deallocate(); _ = close(fd) }

    func record(_ index: Int) throws -> (row: Int, feature: Int, bits: UInt64) {
        guard index >= 0, index < count else { throw VivoOmicsError.invalid("buffered count record index") }
        if index < firstRecord || index >= firstRecord + bufferedRecords {
            try Task.checkCancellation()
            let pageRecords = Self.bufferBytes / 16
            let next = index / pageRecords * pageRecords
            let bytes = min(pageRecords, count - next) * 16
            firstRecord = -1; bufferedRecords = 0
            var done = 0
            while done < bytes {
                try Task.checkCancellation()
                let n = pread(fd, UnsafeMutableRawPointer(buffer).advanced(by: done), bytes - done,
                              off_t(next * 16 + done))
                if n < 0 && errno == EINTR { continue }
                guard n > 0 else { throw VivoOmicsError.invalid("buffered count record read or truncation") }
                done += n
            }
            // Logical full-page loads are deterministic despite short reads or
            // EINTR, so execution receipts can replay their exact accounting.
            bufferLoads += 1; loadedBytes += bytes
            firstRecord = next; bufferedRecords = bytes / 16
        }
        let local = (index - firstRecord) * 2
        let packed = UInt64(littleEndian: buffer[local])
        return (Int(packed & 0xffff_ffff), Int(packed >> 32), UInt64(littleEndian: buffer[local + 1]))
    }
}
