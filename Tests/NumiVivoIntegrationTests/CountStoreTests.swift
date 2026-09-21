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

    @Test func countStoreSnapshotPinsCountsAfterLiveMutation() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("count-store-snapshot-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = try VivoSingleCellExamples.pairedCounts()
        let h5ad = root.appendingPathComponent("source.h5ad")
        let store = root.appendingPathComponent("store")
        let implementation = try VivoFingerprint(bytes: Array(repeating: 13, count: 32))
        let mapping = VivoH5ADImportPlan(
            id: "count-store-snapshot-fixture", evidence: .synthetic,
            sourceDescription: "Synthetic count-store snapshot fixture", countUnit: .umiCount,
            matrixPath: "X", samples: source.samples, sampleColumn: "sample",
            barcodeColumn: "barcode", groupColumn: "group", featureNameColumn: "name",
            mitochondrialFeatureIDs: ["g31"])
        try VivoSingleCellH5AD.write(source, to: h5ad)
        _ = try VivoH5ADCountStore.publish(source: h5ad, plan: mapping, implementation: implementation, to: store)

        let snapshot = try VivoH5ADCountStore.openSnapshot(store, implementation: implementation)
        let original = try snapshot.records.record(0)
        var changed = (original.bits + 1).littleEndian
        let live = try FileHandle(forWritingTo: store.appendingPathComponent("counts.bin"))
        try live.seek(toOffset: 8)
        try withUnsafeBytes(of: &changed) { bytes in
            try live.write(contentsOf: Data(bytes))
        }
        try live.synchronize()
        try live.close()

        let retained = try snapshot.records.record(0)
        #expect(retained.row == original.row)
        #expect(retained.feature == original.feature)
        #expect(retained.bits == original.bits)
        #expect(throws: (any Error).self) {
            _ = try VivoH5ADCountStore.openSnapshot(store, implementation: implementation)
        }
    }
}
