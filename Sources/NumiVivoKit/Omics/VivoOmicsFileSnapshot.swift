import Foundation
import CryptoKit
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Shared immutable snapshot I/O for native omics consumers. The caller supplies
/// its own admission bound; a fixed 1 MiB buffer hashes the exact copied bytes.
enum VivoOmicsFileSnapshot {
    /// Bounded snapshot copying, using descriptors to reject non-regular files.
    static func fingerprint(_ source: URL, copyTo destination: URL? = nil, maximumBytes: Int) throws -> VivoFingerprint {
        guard maximumBytes >= 0, maximumBytes <= 64 * 1_024 * 1_024 * 1_024 else { throw VivoOmicsError.limit("snapshot byte bound") }
        let fd = open(source.path, O_RDONLY | O_NOFOLLOW | O_NONBLOCK)
        guard fd >= 0 else { throw VivoOmicsError.invalid("cannot open omics input") }
        defer { _ = close(fd) }
        var info = stat()
        guard fstat(fd, &info) == 0, (info.st_mode & mode_t(S_IFMT)) == mode_t(S_IFREG), info.st_size >= 0,
              info.st_size <= maximumBytes else { throw VivoOmicsError.limit("omics input type or bytes") }
        var output: Int32 = -1
        if let destination {
            output = open(destination.path, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW, 0o600)
            guard output >= 0 else { throw VivoOmicsError.invalid("cannot create omics snapshot") }
        }
        defer { if output >= 0 { _ = close(output) } }
        var hash = SHA256(), bytes = 0
        var buffer = [UInt8](repeating: 0, count: 1_048_576)
        while true {
            try Task.checkCancellation()
            let count = buffer.withUnsafeMutableBytes { p in
                #if canImport(Darwin)
                Darwin.read(fd, p.baseAddress, p.count)
                #else
                Glibc.read(fd, p.baseAddress, p.count)
                #endif
            }
            if count < 0 { if errno == EINTR { continue }; throw VivoOmicsError.invalid("omics snapshot read") }
            if count == 0 { break }
            guard count <= maximumBytes - bytes else { throw VivoOmicsError.limit("omics input bytes") }
            bytes += count
            try buffer.withUnsafeBytes { full in
                let block = UnsafeRawBufferPointer(rebasing: full[..<count])
                hash.update(bufferPointer: block)
                if output >= 0 {
                    var written = 0
                    while written < count {
                        #if canImport(Darwin)
                        let n = Darwin.write(output, block.baseAddress!.advanced(by: written), count - written)
                        #else
                        let n = Glibc.write(output, block.baseAddress!.advanced(by: written), count - written)
                        #endif
                        if n < 0, errno == EINTR { continue }
                        guard n > 0 else { throw VivoOmicsError.invalid("omics snapshot write") }
                        written += n
                    }
                }
            }
        }
        guard bytes == info.st_size else { throw VivoOmicsError.invalid("omics input changed size during snapshot") }
        if output >= 0 { guard fsync(output) == 0 else { throw VivoOmicsError.invalid("omics snapshot sync") } }
        return try VivoFingerprint(bytes: Array(hash.finalize()))
    }
}
