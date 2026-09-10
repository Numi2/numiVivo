import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellSnapshotTests {
    @Test func downstreamSnapshotsAcceptTheStreamedSourceByteBudget() throws {
        let root=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root,withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let input=root.appendingPathComponent("input"),reference=input.appendingPathComponent("reference")
        try FileManager.default.createDirectory(at: reference,withIntermediateDirectories: true)
        let length: UInt64=1_073_741_824+19
        for directory in [input,reference] {
            let file=directory.appendingPathComponent("original.h5ad")
            try Data("snapshot-boundary-only-not-an-H5AD".utf8).write(to: file)
            let handle=try FileHandle(forWritingTo: file);try handle.truncate(atOffset: length);try handle.close()
            for name in ["plan.json","receipt.json","metadata.json","quality.json","model.json","scores.bin","loadings.bin","report.json"] {
                try Data().write(to: directory.appendingPathComponent(name))
            }
        }
        // Exercise the actual fixed-buffer snapshot boundaries, independently
        // of HDF5 parsing or model qualification. Sparse inputs use no dense RAM.
        let fitted=root.appendingPathComponent("fitted")
        try VivoH5ADPCAQuery.snapshotReference(reference,to: fitted)
        #expect((try FileManager.default.attributesOfItem(atPath: fitted.appendingPathComponent("original.h5ad").path)[.size] as? NSNumber)?.uint64Value == length)
        try FileManager.default.removeItem(at: fitted)
        let query=root.appendingPathComponent("query")
        try VivoPCANeighborBundle.snapshot(input,kind: .query,to: query)
        for file in [query.appendingPathComponent("original.h5ad"),query.appendingPathComponent("reference/original.h5ad")] {
            #expect((try FileManager.default.attributesOfItem(atPath: file.path)[.size] as? NSNumber)?.uint64Value == length)
        }
    }
    @Test func snapshotCopiesExactBytesAndRejectsSmallerAdmission() throws {
        let root=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root,withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source=root.appendingPathComponent("source"),copy=root.appendingPathComponent("copy")
        // Cross the fixed I/O buffer boundary; verify the independent digest.
        let bytes=Data((0..<(1_048_576+19)).map { UInt8($0%251) })
        try bytes.write(to: source)
        let id=try VivoOmicsFileSnapshot.fingerprint(source,copyTo: copy,maximumBytes: bytes.count)
        #expect(id == (try VivoCanonicalJSON.fingerprint(bytes)))
        #expect(try Data(contentsOf: copy) == bytes)
        let rejected=root.appendingPathComponent("rejected")
        #expect(throws: (any Error).self) {
            try VivoOmicsFileSnapshot.fingerprint(source,copyTo: rejected,maximumBytes: bytes.count-1)
        }
        #expect(!FileManager.default.fileExists(atPath: rejected.path))
        #expect(throws: (any Error).self) {
            try VivoOmicsFileSnapshot.fingerprint(source,copyTo: copy,maximumBytes: bytes.count)
        }
        #expect(try Data(contentsOf: copy) == bytes)
    }
    @Test func oversizedSparseSourceIsRejectedBeforeCopying() throws {
        let root=FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root,withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source=root.appendingPathComponent("source"),copy=root.appendingPathComponent("copy")
        #expect(FileManager.default.createFile(atPath: source.path,contents: Data()))
        let handle=try FileHandle(forWritingTo: source)
        // Sparse logical length consumes no 64 GiB payload on disk.
        try handle.truncate(atOffset: UInt64(VivoH5ADPseudobulk.sourceLimits.maximumInputBytes)+1)
        try handle.close()
        #expect(throws: (any Error).self) { try VivoH5ADPseudobulk.fingerprint(source,copyTo: copy) }
        #expect(!FileManager.default.fileExists(atPath: copy.path))
    }
}
