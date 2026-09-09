import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct CountStoreTests {
    @Test func windowBoundaryAndExactUInt64Records() throws {
        let dir = FileManager.default.temporaryDirectory.appendingPathComponent("count-window-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: dir, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: dir) }
        let url = dir.appendingPathComponent("records")
        let writer = try VivoCountRecordWriter(url)
        let n = VivoWindowedCountRecords.windowBytes / 16 + 1
        for i in 0..<n { try writer.append(row: i % 100, feature: i % 50, bits: UInt64.max - UInt64(i)) }
        let hash = try writer.finish()
        #expect(try VivoH5ADCountStore.fingerprint(url) == hash)
        let snapshot = dir.appendingPathComponent("snapshot")
        #expect(throws: (any Error).self) { try VivoOmicsFileSnapshot.fingerprint(url, copyTo: snapshot, maximumBytes: 10) }
        #expect(!FileManager.default.fileExists(atPath: snapshot.path))
        #expect(try VivoOmicsFileSnapshot.fingerprint(url, copyTo: snapshot, maximumBytes: n * 16) == hash)
        #expect(try VivoH5ADCountStore.fingerprint(snapshot) == hash)
        #expect(throws: (any Error).self) { try VivoOmicsFileSnapshot.fingerprint(url, copyTo: snapshot, maximumBytes: n * 16) }
        let mapped = try VivoWindowedCountRecords(snapshot, entries: n)
        for i in [0,n-2,n-1,1] {
            let value = try mapped.record(i)
            #expect(value.row == i % 100 && value.feature == i % 50 && value.bits == UInt64.max - UInt64(i))
        }
        #expect(throws: (any Error).self) { try mapped.record(n) }
        #expect(throws: (any Error).self) { try VivoWindowedCountRecords(url, entries: n-1) }
    }
    @Test func emptyRecordsAndNonRegularInputs() throws {
        let dir = FileManager.default.temporaryDirectory.appendingPathComponent("count-empty-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: dir, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: dir) }
        let url = dir.appendingPathComponent("empty"), writer = try VivoCountRecordWriter(dir.appendingPathComponent("empty"))
        _ = try writer.finish()
        let mapped = try VivoWindowedCountRecords(url, entries: 0)
        #expect(mapped.count == 0)
        #expect(throws: (any Error).self) { try mapped.record(0) }
        #expect(throws: (any Error).self) { try VivoH5ADCountStore.fingerprint(dir) }
        let link = dir.appendingPathComponent("link")
        try FileManager.default.createSymbolicLink(at: link, withDestinationURL: url)
        #expect(throws: (any Error).self) { try VivoWindowedCountRecords(link, entries: 0) }
        #expect(throws: (any Error).self) { try VivoH5ADCountStore.fingerprint(link) }
    }
}
