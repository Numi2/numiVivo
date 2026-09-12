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

    @Test func clonedAndStreamedPublicationsHaveIndependentOwnership() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoRootedFileStore(rootURL: root, createIfNeeded: false)
        let source = root.appendingPathComponent("source")
        let bytes = Data((0..<65_537).map { UInt8($0 % 251) })
        for preferClone in [true, false] {
            try bytes.write(to: source)
            let name = preferClone ? "cloned" : "streamed"
            #expect(try store.writeFile(from: source, relative: name, maximumBytes: bytes.count, immutable: true, preferClone: preferClone))
            let target = root.appendingPathComponent(name)
            let sourceInfo = try FileManager.default.attributesOfItem(atPath: source.path)
            let targetInfo = try FileManager.default.attributesOfItem(atPath: target.path)
            #expect(sourceInfo[.systemFileNumber] as? NSNumber != targetInfo[.systemFileNumber] as? NSNumber)
            #expect((targetInfo[.posixPermissions] as? NSNumber)?.intValue == 0o600)
            let writer = try FileHandle(forWritingTo: source)
            try writer.write(contentsOf: Data([255])); try writer.close()
            #expect(try Data(contentsOf: target) == bytes)
            let targetWriter = try FileHandle(forWritingTo: target)
            try targetWriter.seek(toOffset: 1); try targetWriter.write(contentsOf: Data([254])); try targetWriter.close()
            let changed = try Data(contentsOf: source)
            #expect(changed[0] == 255 && changed[1] == 1)
            #expect(try store.writeFile(from: source, relative: "mutable", maximumBytes: bytes.count, immutable: false, preferClone: preferClone))
            #expect(try store.readFile("mutable", maximumBytes: bytes.count) == changed)
        }
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).allSatisfy { !$0.hasPrefix(".nv-") })
    }
}
