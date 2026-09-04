// Compiled with the actual VivoRootedFileStore.swift, not a mock filesystem.
import Foundation
#if canImport(Darwin)
import Darwin
#elseif canImport(Glibc)
import Glibc
#endif

@main
struct RootedFileStoreChecks {
    static func main() throws {
        let fm = FileManager.default
        let base = fm.temporaryDirectory.appendingPathComponent("numivivo-audit-" + UUID().uuidString)
        defer { try? fm.removeItem(at: base) }
        let root = base.appendingPathComponent("store")
        let outside = base.appendingPathComponent("outside")
        try fm.createDirectory(at: outside, withIntermediateDirectories: true)
        let victim = outside.appendingPathComponent("victim")
        let original = Data("untouched".utf8)
        try original.write(to: victim)
        let files = try VivoRootedFileStore(rootURL: root, createIfNeeded: true)
        var checked = 0
        func expect(_ condition: @autoclosure () throws -> Bool) throws {
            guard try condition() else { throw NSError(domain: "portable-audit", code: 1) }
            checked += 1
        }
        func rejects(_ action: () throws -> Void) throws {
            do { try action() } catch { checked += 1; return }
            throw NSError(domain: "portable-audit-expected-rejection", code: 2)
        }
        try expect(files.writeFile(Data([1,2,3]), relative: "objects/a", immutable: true))
        try expect(files.readFile("objects/a", maximumBytes: 3) == Data([1,2,3]))
        try expect(!files.writeFile(Data([9]), relative: "objects/a", immutable: true))
        try expect(files.readFile("objects/a", maximumBytes: 3) == Data([1,2,3]))
        try rejects { _ = try files.readFile("objects/a", maximumBytes: 2) }
        try rejects { _ = try files.readFile("../outside/victim", maximumBytes: 100) }
        try rejects { _ = try files.writeFile(Data(), relative: "a//b", immutable: true) }
        try fm.createSymbolicLink(atPath: root.appendingPathComponent("link").path, withDestinationPath: outside.path)
        try rejects { _ = try files.readFile("link/victim", maximumBytes: 100) }
        try rejects { _ = try files.writeFile(Data([9]), relative: "link/new", immutable: true) }
        try fm.createSymbolicLink(atPath: root.appendingPathComponent("objects/symlink").path, withDestinationPath: victim.path)
        try rejects { _ = try files.readFile("objects/symlink", maximumBytes: 100) }
        try expect(files.writeFile(Data([7]), relative: "objects/symlink", immutable: false))
        try expect(Data(contentsOf: victim) == original)
        try expect(files.readFile("objects/symlink", maximumBytes: 100) == Data([7]))
        try expect(files.writeFile(Data(), relative: "objects/empty", immutable: true))
        try expect(files.readFile("objects/empty", maximumBytes: 0).isEmpty)
        try rejects { _ = try files.children("objects", maximumCount: 1) }
        let fifo = root.appendingPathComponent("objects/pipe").path
        guard mkfifo(fifo, 0o600) == 0 else { throw NSError(domain: "mkfifo", code: Int(errno)) }
        try rejects { _ = try files.readFile("objects/pipe", maximumBytes: 100) }
        // Pinning survives a rename of the root directory; newly created content
        // cannot redirect the already-open capability into an attacker directory.
        let moved = base.appendingPathComponent("moved")
        try fm.moveItem(at: root, to: moved)
        try fm.createSymbolicLink(atPath: root.path, withDestinationPath: outside.path)
        try expect(files.writeFile(Data([8]), relative: "objects/pinned", immutable: true))
        try expect(fm.fileExists(atPath: moved.appendingPathComponent("objects/pinned").path))
        try expect(!fm.fileExists(atPath: outside.appendingPathComponent("objects/pinned").path))
        print("Rooted-file-store portable checks passed: \(checked)")
    }
}
