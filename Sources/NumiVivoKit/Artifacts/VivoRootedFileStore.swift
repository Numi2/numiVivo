import Foundation
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#endif

/// A root-directory capability. All descendants are opened relative to pinned
/// directory descriptors with O_NOFOLLOW; checking a path and opening it later
/// through Foundation would reintroduce a symlink/rename race.
final class VivoRootedFileStore: @unchecked Sendable {
    enum Failure: Error, Sendable, CustomStringConvertible {
        case invalidPath(String)
        case missing(String)
        case exceededLimit(String)
        case io(String, Int32)
        var description: String {
            switch self {
            case .invalidPath(let path): return "Unsafe artifact-store path: \(path)"
            case .missing(let path): return "Missing artifact-store path: \(path)"
            case .exceededLimit(let path): return "Artifact-store I/O bound exceeded: \(path)"
            case .io(let operation, let code): return "Artifact-store \(operation) failed (errno \(code))"
            }
        }
    }
    let rootURL: URL
    private let root: Int32
    init(rootURL: URL, createIfNeeded: Bool) throws {
        guard rootURL.isFileURL else { throw Failure.invalidPath("root is not a file URL") }
        self.rootURL = rootURL.standardizedFileURL
        if createIfNeeded {
            try FileManager.default.createDirectory(at: self.rootURL, withIntermediateDirectories: true,
                                                     attributes: [.posixPermissions: 0o700])
        }
        let fd = open(self.rootURL.path, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW)
        guard fd >= 0 else { throw Failure.io("open root", errno) }
        root = fd
    }
    deinit { _ = close(root) }

    private func components(_ relative: String) throws -> [String] {
        let values = relative.split(separator: "/", omittingEmptySubsequences: false).map(String.init)
        guard !values.isEmpty, !relative.hasPrefix("/"), relative.utf8.count <= 4096,
              values.allSatisfy({ !$0.isEmpty && $0 != "." && $0 != ".." &&
                                  !$0.contains("\0") && !$0.contains("\\") && $0.utf8.count <= 240 }) else {
            throw Failure.invalidPath(relative)
        }
        return values
    }
    private func directory(_ parts: ArraySlice<String>, create: Bool) throws -> Int32 {
        var parent = openat(root, ".", O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW)
        guard parent >= 0 else { throw Failure.io("reopen root", errno) }
        do {
            for part in parts {
                if create, mkdirat(parent, part, 0o700) != 0, errno != EEXIST {
                    throw Failure.io("mkdirat", errno)
                }
                let next = openat(parent, part, O_RDONLY | O_DIRECTORY | O_CLOEXEC | O_NOFOLLOW)
                guard next >= 0 else {
                    if errno == ENOENT { throw Failure.missing(part) }
                    throw Failure.io("open directory", errno)
                }
                _ = close(parent)
                parent = next
            }
            return parent
        } catch { _ = close(parent); throw error }
    }
    func createDirectory(_ relative: String) throws {
        let fd = try directory(try components(relative)[...], create: true)
        _ = close(fd)
    }
    func isRegularFile(_ relative: String) throws -> Bool {
        let parts = try components(relative)
        let parent = try directory(parts.dropLast(), create: false)
        defer { _ = close(parent) }
        let fd = openat(parent, parts.last!, O_RDONLY | O_NOFOLLOW | O_CLOEXEC | O_NONBLOCK)
        guard fd >= 0 else {
            if errno == ENOENT { return false }
            throw Failure.io("open object", errno)
        }
        defer { _ = close(fd) }
        var info = stat()
        guard fstat(fd, &info) == 0 else { throw Failure.io("fstat object", errno) }
        return info.st_mode & mode_t(S_IFMT) == mode_t(S_IFREG)
    }
    func readFile(_ relative: String, maximumBytes: Int) throws -> Data {
        guard maximumBytes >= 0 else { throw Failure.exceededLimit(relative) }
        let parts = try components(relative)
        let parent = try directory(parts.dropLast(), create: false)
        defer { _ = close(parent) }
        // O_NONBLOCK prevents an attacker-supplied FIFO from blocking before
        // fstat has established that this is a regular file.
        let fd = openat(parent, parts.last!, O_RDONLY | O_NOFOLLOW | O_CLOEXEC | O_NONBLOCK)
        guard fd >= 0 else {
            if errno == ENOENT { throw Failure.missing(relative) }
            throw Failure.io("open object", errno)
        }
        defer { _ = close(fd) }
        var before = stat()
        guard fstat(fd, &before) == 0 else { throw Failure.io("fstat object", errno) }
        guard before.st_mode & mode_t(S_IFMT) == mode_t(S_IFREG) else { throw Failure.invalidPath(relative) }
        guard before.st_size >= 0, UInt64(before.st_size) <= UInt64(maximumBytes) else { throw Failure.exceededLimit(relative) }
        var result = Data(count: Int(before.st_size))
        try result.withUnsafeMutableBytes { bytes in
            var position = 0
            while position < bytes.count {
                let amount = read(fd, bytes.baseAddress!.advanced(by: position), bytes.count - position)
                if amount < 0 && errno == EINTR { continue }
                guard amount > 0 else { throw Failure.io("read object", amount < 0 ? errno : EIO) }
                position += amount
            }
        }
        var extra: UInt8 = 0
        var extraCount: Int
        repeat { extraCount = read(fd, &extra, 1) } while extraCount < 0 && errno == EINTR
        guard extraCount == 0 else { throw Failure.io("object changed during read", EIO) }
        var after = stat()
        guard fstat(fd, &after) == 0, after.st_size == before.st_size else { throw Failure.io("object changed during read", EIO) }
        return result
    }

    /// Returns false when an immutable destination already exists. Mutable
    /// publication replaces only the directory entry, never follows its target.
    @discardableResult
    func writeFile(_ data: Data, relative: String, immutable: Bool) throws -> Bool {
        let parts = try components(relative)
        let parent = try directory(parts.dropLast(), create: true)
        defer { _ = close(parent) }
        let temporary = ".nv-\(UUID().uuidString).tmp"
        let fd = openat(parent, temporary, O_WRONLY | O_CREAT | O_EXCL | O_NOFOLLOW | O_CLOEXEC, 0o600)
        guard fd >= 0 else { throw Failure.io("create temporary", errno) }
        defer { _ = close(fd); _ = unlinkat(parent, temporary, 0) }
        try data.withUnsafeBytes { bytes in
            var position = 0
            while position < bytes.count {
                let amount = write(fd, bytes.baseAddress!.advanced(by: position), bytes.count - position)
                if amount < 0 && errno == EINTR { continue }
                guard amount > 0 else { throw Failure.io("write object", amount < 0 ? errno : EIO) }
                position += amount
            }
        }
        guard fsync(fd) == 0 else { throw Failure.io("fsync object", errno) }
        if immutable {
            if linkat(parent, temporary, parent, parts.last!, 0) != 0 {
                if errno == EEXIST { return false }
                throw Failure.io("publish immutable object", errno)
            }
        } else if renameat(parent, temporary, parent, parts.last!) != 0 {
            throw Failure.io("publish reference", errno)
        }
        // Some file systems do not implement directory fsync. Atomic name
        // publication remains valid there; crash durability is filesystem-specific.
        if fsync(parent) != 0 && errno != EINVAL && errno != ENOTSUP {
            throw Failure.io("fsync parent", errno)
        }
        return true
    }
    func removeFile(_ relative: String) throws {
        let parts = try components(relative)
        let parent = try directory(parts.dropLast(), create: false)
        defer { _ = close(parent) }
        guard unlinkat(parent, parts.last!, 0) == 0 else {
            if errno == ENOENT { throw Failure.missing(relative) }
            throw Failure.io("remove reference", errno)
        }
    }
    func children(_ relative: String, maximumCount: Int) throws -> [String] {
        let fd = try directory(try components(relative)[...], create: false)
        guard let directory = fdopendir(fd) else { _ = close(fd); throw Failure.io("fdopendir", errno) }
        defer { _ = closedir(directory) }
        var names: [String] = []
        while true {
            errno = 0
            guard let entry = readdir(directory) else {
                if errno != 0 { throw Failure.io("readdir", errno) }
                break
            }
            var storage = entry.pointee.d_name
            let name = withUnsafeBytes(of: &storage) { bytes in
                String(bytes: bytes.prefix { $0 != 0 }, encoding: .utf8)
            }
            guard let name else { throw Failure.invalidPath("non-UTF8 directory entry") }
            if name == "." || name == ".." || name.hasPrefix(".nv-") { continue }
            guard names.count < maximumCount else { throw Failure.exceededLimit(relative) }
            _ = try components(name)
            names.append(name)
        }
        return names.sorted()
    }
}
