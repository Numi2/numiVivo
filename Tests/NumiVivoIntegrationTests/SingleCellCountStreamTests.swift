import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellCountStreamTests {
    static let implementation = try! VivoFingerprint(bytes: Array(repeating: 3, count: 32))
    static func fixture() -> VivoSingleCellDataset {
        .init(id: "stream-control", evidence: .synthetic, sourceDescription: "Numerical control, not biological evidence", countUnit: .umiCount,
              samples: [.init(id: "s", biologicalReplicateID: "d", donorID: "d", condition: "PBS", batchID: "b", organism: "human")],
              features: [.init(id: "A", name: "A"), .init(id: "B", name: "B", mitochondrial: true)],
              cells: [.init(barcode: "c0", sampleID: "s"), .init(barcode: "empty", sampleID: "s"), .init(barcode: "c2", sampleID: "s")],
              matrix: .init(cellCount: 3, featureCount: 2, rowOffsets: [0,2,2,3], featureIndices: [0,1,1], counts: [4,2,7]))
    }
    static func plan(_ data: VivoSingleCellDataset = fixture()) -> VivoCountStreamPlan {
        let quality = try! VivoSingleCellAnalysis.quality(data)
        return .init(metadata: data.metadata, sourceDeclaration: "Synthetic stream fixture",
                     rowNonzeros: quality.map(\.detectedFeatures), rowTotals: quality.map(\.totalCounts))
    }
    static func bytes(_ records: [(UInt32, UInt32, UInt64)] = [(0,0,4),(0,1,2),(2,1,7)]) -> Data {
        var bytes = Data()
        for (row, feature, count) in records {
            var packed = (UInt64(row) | UInt64(feature) << 32).littleEndian, value = count.littleEndian
            withUnsafeBytes(of: &packed) { bytes.append(contentsOf: $0) }
            withUnsafeBytes(of: &value) { bytes.append(contentsOf: $0) }
        }
        return bytes
    }
    static func evaluate(_ bytes: Data, plan: VivoCountStreamPlan = plan(), width: Int = 17) throws -> VivoCountStreamReport {
        var position = 0
        return try VivoCountStreamPseudobulk.evaluate(plan: plan) {
            let stop = min(position + width, bytes.count)
            defer { position = stop }
            return Data(bytes[position..<stop])
        }.0
    }
    @Test func fragmentedStreamMatchesIndependentResidentOwnersAndPreservesEmptyCells() throws {
        let source = Self.fixture()
        for width in [1,7,16,17,48,1_048_576] {
            let report = try Self.evaluate(Self.bytes(), width: width)
            #expect(report.quality == (try VivoSingleCellAnalysis.quality(source)))
            #expect(report.pseudobulk == (try VivoSingleCellAnalysis.pseudobulk(source)))
            #expect(report.quality[1].totalCounts == 0 && report.canonicalNonzeros == 3)
        }
        let noIndependentTotals = VivoCountStreamPlan(metadata: source.metadata, sourceDeclaration: "No independent totals supplied", rowNonzeros: [2,0,1])
        #expect(try Self.evaluate(Self.bytes(), plan: noIndependentTotals).quality == VivoSingleCellAnalysis.quality(source))
    }
    @Test func rejectsTruncationExtraRecordsDuplicatesUnsortedAndOutOfAxis() throws {
        let corrupt = [Data(Self.bytes().dropLast()), Self.bytes() + Data([0]),
                       Self.bytes([(0,0,4),(0,0,2),(2,1,7)]), Self.bytes([(0,1,2),(0,0,4),(2,1,7)]),
                       Self.bytes([(0,0,4),(3,1,2),(2,1,7)]), Self.bytes([(0,0,4),(0,2,2),(2,1,7)]),
                       Self.bytes([(0,0,4),(0,1,0),(2,1,7)]), Self.bytes([(0,0,4),(0,1,2),(2,1,6)])]
        for data in corrupt { #expect(throws: (any Error).self) { try Self.evaluate(data) } }
        #expect(throws: (any Error).self) {
            try VivoCountStreamPseudobulk.evaluate(plan: Self.plan()) { Data(repeating: 0, count: 1_048_577) }
        }
    }
    @Test func rejectsAxisContractsOverflowAndUnknownPlanFields() throws {
        for (nonzeros, totals) in [([2,0], [UInt64(6),0]), ([2,0,1],[6,1,7]), ([-1,0,1],[6,0,7]), ([3,0,1],[6,0,7])] {
            #expect(throws: (any Error).self) {
                try VivoCountStreamPlan(metadata: Self.fixture().metadata, sourceDeclaration: "fixture", rowNonzeros: nonzeros, rowTotals: totals).validate()
            }
        }
        let overflow = VivoCountStreamPlan(metadata: Self.fixture().metadata, sourceDeclaration: "fixture", rowNonzeros: [2,0,1], rowTotals: [UInt64.max,0,7])
        #expect(throws: (any Error).self) { try Self.evaluate(Self.bytes([(0,0,UInt64.max),(0,1,1),(2,1,7)]), plan: overflow) }
        let aggregateOverflow = VivoCountStreamPlan(metadata: Self.fixture().metadata, sourceDeclaration: "fixture", rowNonzeros: [2,0,1], rowTotals: [UInt64.max,0,7])
        #expect(throws: (any Error).self) { try Self.evaluate(Self.bytes([(0,0,1),(0,1,UInt64.max-1),(2,1,7)]), plan: aggregateOverflow) }
        var json = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(Self.plan())) as? [String: Any])
        json["silentSelection"] = true
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoCountStreamPlan.self, from: JSONSerialization.data(withJSONObject: json)) }
    }
    @Test func publishesTransactionallyAndReconstructsWithReplay() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("count-stream-test-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("counts.bin"), output = root.appendingPathComponent("bundle")
        try Self.bytes().write(to: source)
        func handle() throws -> FileHandle { try FileHandle(forReadingFrom: source) }
        let input = try handle(); defer { try? input.close() }
        let receipt = try VivoCountStreamPseudobulk.publish(plan: Self.plan(), input: input, implementation: Self.implementation, to: output)
        #expect(receipt.streamBytes == 48)
        #expect(receipt.stream == (try VivoCanonicalJSON.fingerprint(Self.bytes())))
        let replay = try handle(); defer { try? replay.close() }
        #expect(try VivoCountStreamPseudobulk.verify(output, input: replay, implementation: Self.implementation) == receipt)
        // Forging the report and its receipt hash cannot evade stream reconstruction.
        let reportURL = output.appendingPathComponent("report.json"), original = try Data(contentsOf: reportURL)
        var object = try #require(JSONSerialization.jsonObject(with: original) as? [String: Any])
        object["method"] = "forged"
        let altered = try JSONSerialization.data(withJSONObject: object, options: [.sortedKeys])
        try altered.write(to: reportURL)
        let forged = try VivoCountStreamReceipt(schemaVersion: 1, encoding: receipt.encoding, streamBytes: receipt.streamBytes,
            stream: receipt.stream, plan: receipt.plan, report: VivoCanonicalJSON.fingerprint(altered), implementation: receipt.implementation)
        try VivoCanonicalJSON.encode(forged).write(to: output.appendingPathComponent("receipt.json"))
        let recheck = try handle(); defer { try? recheck.close() }
        #expect(throws: (any Error).self) { try VivoCountStreamPseudobulk.verify(output, input: recheck, implementation: Self.implementation) }
        try Self.bytes().dropLast().write(to: source)
        let broken = try handle(); defer { try? broken.close() }
        let failed = root.appendingPathComponent("failed")
        #expect(throws: (any Error).self) { try VivoCountStreamPseudobulk.publish(plan: Self.plan(), input: broken, implementation: Self.implementation, to: failed) }
        #expect(!FileManager.default.fileExists(atPath: failed.path))
        #expect(!(try FileManager.default.contentsOfDirectory(atPath: root.path)).contains { $0.hasPrefix(".numivivo-count-stream-") })
    }
}
