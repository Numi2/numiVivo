import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct RootedStreamingPublicationTests {
    @Test func streamingCopyPreservesBytesBoundsAndImmutableDestination() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoRootedFileStore(rootURL: root, createIfNeeded: false)
        let source = root.appendingPathComponent("source")
        let bytes = Data((0..<(2_097_152 + 31)).map { UInt8($0 % 251) })
        try bytes.write(to: source)
        #expect(try store.writeFile(from: source, relative: "copy", maximumBytes: bytes.count, immutable: true))
        #expect(try store.readFile("copy", maximumBytes: bytes.count) == bytes)
        #expect(try !store.writeFile(from: source, relative: "copy", maximumBytes: bytes.count, immutable: true))
        #expect(throws: (any Error).self) {
            try store.writeFile(from: source, relative: "oversized", maximumBytes: bytes.count - 1, immutable: true)
        }
        try FileManager.default.createSymbolicLink(atPath: root.appendingPathComponent("alias").path, withDestinationPath: source.path)
        #expect(throws: (any Error).self) {
            try store.writeFile(from: root.appendingPathComponent("alias"), relative: "linked", maximumBytes: bytes.count, immutable: true)
        }
        #expect(try !store.writeFile(from: source, relative: "alias", maximumBytes: bytes.count, immutable: true))
        #expect(throws: (any Error).self) {
            try store.writeFile(from: root, relative: "directory", maximumBytes: bytes.count, immutable: true)
        }
        #expect(throws: (any Error).self) {
            try store.writeFile(from: source, relative: "../escape", maximumBytes: bytes.count, immutable: true)
        }
        #expect(try Data(contentsOf: source) == bytes)
        #expect(try Set(FileManager.default.contentsOfDirectory(atPath: root.path)) == ["source", "copy", "alias"])
    }
}
