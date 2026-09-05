import Foundation
#if canImport(Darwin)
import Darwin
#else
import Glibc
#endif

/// Explicit, caller-selected exports; not a replacement for library-relative
/// artifact-store access. Parent directories must already exist. Publication
/// never follows a final-component symlink and preserves the prior file on a
/// failed write. JSON checkpoint updates replace one complete file atomically.
public enum VivoMDAtomicFileExport {
    public static func write(_ data: Data, to url: URL, overwrite: Bool,
                             maximumBytes: UInt64 = 512 << 20) throws {
        guard url.isFileURL, UInt64(data.count) <= maximumBytes else {
            throw VivoArtifactValidationError.invalid("MD export requires a bounded file payload")
        }
        let target = url.standardizedFileURL
        let name = target.lastPathComponent
        guard !name.isEmpty, name != ".", name != "..", !name.contains("/"), !name.contains("\0") else {
            throw VivoArtifactValidationError.invalid("invalid MD export filename")
        }
        let directory = open(target.deletingLastPathComponent().path, O_RDONLY | O_DIRECTORY | O_CLOEXEC)
        guard directory >= 0 else { throw io("cannot open export directory") }
        defer { close(directory) }
        let temporary = ".numivivo-export-\(UUID().uuidString).partial"
        let descriptor = openat(directory, temporary, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC, mode_t(0o600))
        guard descriptor >= 0 else { throw io("cannot exclusively create export staging file") }
        defer { _ = unlinkat(directory, temporary, 0) }
        let handle = FileHandle(fileDescriptor: descriptor, closeOnDealloc: true)
        do {
            try handle.write(contentsOf: data)
            try handle.synchronize()
            try handle.close()
        } catch {
            try? handle.close()
            throw error
        }
        let status = overwrite
            ? renameat(directory, temporary, directory, name)
            : linkat(directory, temporary, directory, name, 0)
        guard status == 0 else { throw io("cannot publish export atomically") }
        if fsync(directory) != 0, errno != EINVAL, errno != ENOTSUP {
            throw io("export was published, but directory synchronization failed")
        }
    }

    public static func read(_ url: URL, maximumBytes: UInt64 = 256 << 20) throws -> Data {
        guard url.isFileURL, maximumBytes <= UInt64(Int.max) else {
            throw VivoArtifactValidationError.invalid("invalid MD input path or size limit")
        }
        let descriptor = open(url.path, O_RDONLY | O_NOFOLLOW | O_CLOEXEC)
        guard descriptor >= 0 else { throw io("cannot open MD input") }
        let handle = FileHandle(fileDescriptor: descriptor, closeOnDealloc: true)
        defer { try? handle.close() }
        var info = stat()
        guard fstat(descriptor, &info) == 0, (info.st_mode & S_IFMT) == S_IFREG,
              info.st_size >= 0, UInt64(info.st_size) <= maximumBytes else {
            throw VivoArtifactValidationError.invalid("MD input must be a bounded regular file")
        }
        let expected = Int(info.st_size)
        var data = Data(capacity: expected)
        while data.count < expected {
            guard let next = try handle.read(upToCount: min(expected - data.count, 1 << 20)), !next.isEmpty else {
                throw VivoArtifactValidationError.invalid("MD input changed or was truncated during read")
            }
            data.append(next)
        }
        if let extra = try handle.read(upToCount: 1), !extra.isEmpty {
            throw VivoArtifactValidationError.invalid("MD input grew during read")
        }
        return data
    }

    private static func io(_ reason: String) -> VivoMDTrajectoryArchiveError {
        .io("\(reason) (errno \(errno))")
    }
}
