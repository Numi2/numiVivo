import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct MultiAssayH5MUImportTests {
    static func plan(matrix: String = "X") -> VivoH5MUMultiAssayPlan {
        let d = MultiAssayTests.fixture()
        return .init(schemaVersion: 1, id: d.id, evidence: .synthetic, sourceDescription: d.sourceDescription,
            samples: d.samples, sampleColumn: "sample", barcodeColumn: "barcode", groupColumn: nil,
            defaultObservationKind: .cell, observationKindColumn: "observation_kind", spatial: nil,
            assays: d.assays.map { .init(sourceName: $0.id, id: $0.id, kind: $0.kind, featureNamespace: $0.featureNamespace,
                countUnit: $0.countUnit, genomeAssembly: $0.genomeAssembly, matrixPath: matrix, featureIDColumn: nil, featureNameColumn: "name", peakIDConvention: nil) })
    }
    @Test func mappingRejectsUnknownOptionsAndAmbiguousCountSelection() throws {
        try VivoMultiAssayH5MUImport.validate(Self.plan())
        #expect(throws: (any Error).self) { try VivoMultiAssayH5MUImport.validate(Self.plan(matrix: "raw/X")) }
        #expect(throws: (any Error).self) { try VivoMultiAssayH5MUImport.validate(Self.plan(matrix: "layers/../X")) }
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(Self.plan())) as? [String: Any])
        object["ignoreMissingModalities"] = true
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoH5MUMultiAssayPlan.self, from: JSONSerialization.data(withJSONObject: object)) }
    }
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1")) func spatialArraysRoundTripPositionsAndMissingRows() throws {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent("h5mu-spatial-" + UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: url) }
        let fixture = MultiAssayTests.fixture(), base = Self.plan()
        try VivoMultiAssayH5MU.writeSnapshot(fixture, to: url)
        let plan = VivoH5MUMultiAssayPlan(schemaVersion: 1, id: base.id, evidence: base.evidence,
            sourceDescription: base.sourceDescription, samples: base.samples, sampleColumn: base.sampleColumn,
            barcodeColumn: base.barcodeColumn, groupColumn: nil, defaultObservationKind: .cell,
            observationKindColumn: base.observationKindColumn,
            spatial: .init(path: "obsm/spatial", frame: fixture.spatialFrames[0]), assays: base.assays)
        let result = try VivoMultiAssayH5MUImport.readSnapshot(url, plan: plan)
        #expect(result.observations == fixture.observations)
        #expect(result.spatialFrames == fixture.spatialFrames)
        #expect(result.assays == fixture.assays)
    }
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1")) func nativeReadPreservesPartialMapsKindsAndExactIntegers() throws {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent("h5mu-read-" + UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: url) }
        try VivoMultiAssayH5MU.writeSnapshot(MultiAssayTests.fixture(), to: url)
        let d = try VivoMultiAssayH5MUImport.readSnapshot(url, plan: Self.plan())
        #expect(d.assays[1].observationIndices == [2,0])
        #expect(d.observations[2].kind == .spot)
        #expect(try d.count(assayID: "rna", observationIndex: 0, featureID: "shared-id") == 9_007_199_254_740_993)
        #expect(try d.count(assayID: "protein", observationIndex: 1, featureID: "shared-id") == nil)
        #expect(try d.count(assayID: "protein", observationIndex: 0, featureID: "shared-id") == 0)
    }
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1")) func conflictingAndIncompleteObservationMapsReject() throws {
        for values: [UInt64] in [[2,0,9],[1,0,2],[2,0,0]] {
            let url = FileManager.default.temporaryDirectory.appendingPathComponent("h5mu-map-" + UUID().uuidString)
            defer { try? FileManager.default.removeItem(at: url) }
            try VivoMultiAssayH5MU.writeSnapshot(MultiAssayTests.fixture(), to: url)
            try VivoHDF5.lock.withLock {
                let h = try VivoHDF5(), f = try h.file(url.path, writable: true); defer { h.close(f, "H5Fclose") }
                let d = try h.dataset(f, "obsmap/protein"); defer { h.close(d, "H5Dclose") }
                let write: @convention(c) (Int64, Int64, Int64, Int64, Int64, UnsafeRawPointer?) -> Int32 = try h.symbol("H5Dwrite")
                try values.withUnsafeBytes { try h.check(write(d, h.native("NATIVE_ULLONG"), 0, 0, 0, $0.baseAddress), "mutate fixture map") }
            }
            #expect(throws: (any Error).self) { try VivoMultiAssayH5MUImport.readSnapshot(url, plan: Self.plan()) }
        }
    }
}
