import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct KineticsPreparedOutputTests {
    private func directory() throws -> URL {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-prepared-output-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        return root
    }

    @Test func parentRenameAndReplacementCannotRedirectPreparedOutput() throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let parent = root.appendingPathComponent("parent"), moved = root.appendingPathComponent("captured-parent")
        let destination = parent.appendingPathComponent("result.json")
        let prepared = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: destination)
        #expect(try FileManager.default.contentsOfDirectory(atPath: parent.path).isEmpty)
        try FileManager.default.moveItem(at: parent, to: moved)
        try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: false)
        let data = Data("captured parent".utf8)
        try prepared.write(data)
        #expect(try Data(contentsOf: moved.appendingPathComponent("result.json")) == data)
        #expect(try FileManager.default.contentsOfDirectory(atPath: parent.path).isEmpty)
        #expect(!FileManager.default.fileExists(atPath: destination.path))
    }

    @Test func ancestorReplacementWithSymlinkCannotRedirectPreparedOutput() throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let ancestor = root.appendingPathComponent("ancestor"), moved = root.appendingPathComponent("captured-ancestor")
        let redirected = root.appendingPathComponent("redirected")
        let destination = ancestor.appendingPathComponent("child/result.json")
        let prepared = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: destination)
        try FileManager.default.moveItem(at: ancestor, to: moved)
        try FileManager.default.createDirectory(at: redirected.appendingPathComponent("child"), withIntermediateDirectories: true)
        try FileManager.default.createSymbolicLink(at: ancestor, withDestinationURL: redirected)
        let data = Data("captured ancestor".utf8)
        try prepared.write(data)
        #expect(try Data(contentsOf: moved.appendingPathComponent("child/result.json")) == data)
        #expect(try FileManager.default.contentsOfDirectory(atPath: redirected.appendingPathComponent("child").path).isEmpty)
        #expect(!FileManager.default.fileExists(atPath: destination.path))
    }

    @Test(arguments: [false, true])
    func existingLeafRemainsNoClobberBeforeOrAfterPreparation(_ collisionAfterPreparation: Bool) throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let destination = root.appendingPathComponent("result.json"), original = Data("original".utf8)
        if !collisionAfterPreparation { try original.write(to: destination, options: .withoutOverwriting) }
        let prepared = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: destination)
        if collisionAfterPreparation { try original.write(to: destination, options: .withoutOverwriting) }
        do {
            try prepared.write(Data("replacement".utf8))
            Issue.record("prepared publication clobbered an existing leaf")
        } catch VivoKineticsError.invalid(let message) {
            #expect(message == "output already exists; explicitly request overwrite")
        }
        #expect(try Data(contentsOf: destination) == original)
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path) == ["result.json"])
    }

    @Test func symlinkLeafRemainsUntouchedAndNeverOverwritesItsTarget() throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let destination = root.appendingPathComponent("result.json"), target = root.appendingPathComponent("target.json")
        let prepared = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: destination)
        let original = Data("protected target".utf8)
        try original.write(to: target, options: .withoutOverwriting)
        try FileManager.default.createSymbolicLink(at: destination, withDestinationURL: target)
        do {
            try prepared.write(Data("replacement".utf8))
            Issue.record("prepared publication replaced an existing symlink leaf")
        } catch VivoKineticsError.invalid(_) {}
        #expect(try Data(contentsOf: target) == original)
        #expect(try FileManager.default.destinationOfSymbolicLink(atPath: destination.path) == target.path)
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).sorted() == ["result.json", "target.json"])
    }

    @Test func byteAdmissionRejectsBeforePublicationAndAllowsEmptyDocuments() throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let destination = root.appendingPathComponent("result.json")
        let prepared = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: destination, maximumBytes: 1)
        do {
            try prepared.write(Data([1, 2]))
            Issue.record("prepared publication ignored its byte admission limit")
        } catch VivoKineticsError.capacity(let message) {
            #expect(message == "document output")
        }
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
        try prepared.write(Data([1]))
        #expect(try Data(contentsOf: destination) == Data([1]))
        let empty = root.appendingPathComponent("empty.json")
        try VivoKineticsDocumentIO.prepareNoClobberOutput(to: empty, maximumBytes: 1).write(Data())
        #expect(try Data(contentsOf: empty).isEmpty)
    }

    @Test func invalidPathAndLimitAdmissionPrecedeParentCreation() throws {
        let root = try directory(); defer { try? FileManager.default.removeItem(at: root) }
        let absentParent = root.appendingPathComponent("must-not-create")
        for limit in [-1, 0, 536_870_913] {
            do {
                _ = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: absentParent.appendingPathComponent("result.json"),
                    maximumBytes: limit)
                Issue.record("invalid prepared output byte limit was accepted")
            } catch VivoKineticsError.invalid(_) {}
        }
        let invalidPaths = [try #require(URL(string: "https://example.invalid/result.json")),
            absentParent.appendingPathComponent("bad\\leaf.json"),
            absentParent.appendingPathComponent(String(repeating: "a", count: 241)),
            URL(fileURLWithPath: "/")]
        for path in invalidPaths {
            do {
                _ = try VivoKineticsDocumentIO.prepareNoClobberOutput(to: path)
                Issue.record("invalid prepared output path was accepted")
            } catch VivoKineticsError.invalid(_) {}
        }
        #expect(!FileManager.default.fileExists(atPath: absentParent.path))
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).isEmpty)
        let longest = root.appendingPathComponent(String(repeating: "a", count: 240))
        try VivoKineticsDocumentIO.prepareNoClobberOutput(to: longest).write(Data())
        #expect(try Data(contentsOf: longest).isEmpty)
    }
}
