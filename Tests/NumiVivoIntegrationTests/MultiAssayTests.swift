import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct MultiAssayTests {
    static func fixture(proteinRows: [Int] = [2,0], frameID: String = "slide", duplicateFeature: Bool = false) -> VivoMultiAssayDataset {
        let sample = VivoOmicsSample(id: "sample", biologicalReplicateID: "rep", condition: "control", batchID: "batch", organism: "NCBITaxon:9606")
        let observations: [VivoAssayObservation] = (0..<3).map { .init(identity: .init(barcode: "b\($0)", sampleID: "sample"), kind: $0 == 2 ? .spot : .cell,
            position: $0 == 2 ? .init(frameID: frameID, coordinates: [2.5,4.0]) : nil) }
        let gene = VivoAssayFeature(id: "shared-id", name: "RNA", interval: nil)
        let protein = VivoAssayFeature(id: "shared-id", name: "Antibody", interval: nil)
        let rna = VivoAssaySpace(id: "rna", kind: .rna, featureNamespace: "test-gene", countUnit: .umiCount, genomeAssembly: nil,
            sourceDescription: "Synthetic structural fixture", features: duplicateFeature ? [gene,gene] : [gene], observationIndices: [0,1,2],
            matrix: .init(cellCount: 3, featureCount: duplicateFeature ? 2 : 1, rowOffsets: [0,1,1,2], featureIndices: [0,0], counts: [9_007_199_254_740_993,4]))
        let prot = VivoAssaySpace(id: "protein", kind: .antibodyCapture, featureNamespace: "test-antibody", countUnit: .umiCount, genomeAssembly: nil,
            sourceDescription: "Synthetic structural fixture", features: [protein], observationIndices: proteinRows,
            matrix: .init(cellCount: 2, featureCount: 1, rowOffsets: [0,1,1], featureIndices: [0], counts: [7]))
        return .init(schemaVersion: 1, id: "multi", evidence: .synthetic, sourceDescription: "Synthetic structural fixture", samples: [sample],
            observations: observations, spatialFrames: [.init(id: "slide", unit: .micrometer, axes: ["x","y"], sourceDescription: "Synthetic coordinate frame")], assays: [rna,prot])
    }
    @Test func missingMeasurementsDifferFromMeasuredZerosAndNamespacesRemainSeparate() throws {
        let d = Self.fixture(); try d.validate()
        #expect(try d.count(assayID: "protein", observationIndex: 1, featureID: "shared-id") == nil)
        #expect(try d.count(assayID: "protein", observationIndex: 0, featureID: "shared-id") == 0)
        #expect(try d.count(assayID: "protein", observationIndex: 2, featureID: "shared-id") == 7)
        #expect(try d.count(assayID: "rna", observationIndex: 0, featureID: "shared-id") == 9_007_199_254_740_993)
        #expect(try VivoCanonicalJSON.decode(VivoMultiAssayDataset.self, from: VivoCanonicalJSON.encode(d)) == d)
    }
    @Test func invalidAlignmentFramesAndFeatureIdentitiesReject() throws {
        for d in [Self.fixture(proteinRows: [0,0]), Self.fixture(proteinRows: [3,0]), Self.fixture(frameID: "missing"), Self.fixture(duplicateFeature: true)] {
            #expect(throws: (any Error).self) { try d.validate() }
        }
        let bad = Data("{\"featureType\":\"Antibody Capture\",\"id\":\"adt\",\"kind\":\"antibodyCapture\",\"featureNamespace\":\"biolegend\",\"countUnit\":\"umiCount\",\"genomeAssembly\":null,\"ignoreUnknown\":true}".utf8)
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoTenXAssayMapping.self, from: bad) }
    }
    @Test func accessibilityRequiresAssemblyIntervalsAndCorrectUnits() throws {
        let original = Self.fixture()
        func make(_ assembly: String?, _ unit: VivoAssayCountUnit, _ interval: VivoGenomicInterval?) -> VivoMultiAssayDataset {
            let a = VivoAssaySpace(id: "atac", kind: .chromatinAccessibility, featureNamespace: "test-peaks", countUnit: unit, genomeAssembly: assembly,
                sourceDescription: "Synthetic peaks", features: [.init(id: "peak", name: "peak", interval: interval)], observationIndices: [1],
                matrix: .init(cellCount: 1, featureCount: 1, rowOffsets: [0,1], featureIndices: [0], counts: [2]))
            return .init(schemaVersion: 1, id: original.id, evidence: original.evidence, sourceDescription: original.sourceDescription,
                samples: original.samples, observations: original.observations, spatialFrames: original.spatialFrames, assays: original.assays + [a])
        }
        let interval = VivoGenomicInterval(contig: "chr1", start: 10, end: 20)
        try make("GRCh38", .fragmentCount, interval).validate()
        try make("GRCh38", .cutSiteCount, interval).validate()
        for d in [make(nil, .fragmentCount, interval), make("GRCh38", .umiCount, interval), make("GRCh38", .fragmentCount, nil),
                  make("GRCh38", .fragmentCount, .init(contig: "chr1", start: 20, end: 10))] {
            #expect(throws: (any Error).self) { try d.validate() }
        }
    }
    @Test func cutSiteUnitsCannotBeAppliedToRNA() throws {
        var value = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(Self.fixture())) as? [String: Any])
        var assays = try #require(value["assays"] as? [[String: Any]])
        assays[0]["countUnit"] = "cutSiteCount"; value["assays"] = assays
        let data = try VivoCanonicalJSON.decode(VivoMultiAssayDataset.self, from: JSONSerialization.data(withJSONObject: value))
        #expect(throws: (any Error).self) { try data.validate() }
        let plan = VivoTenXMultiAssayPlan(schemaVersion: 1, id: "unit", evidence: .synthetic,
            sourceDescription: "Synthetic wrong unit", sample: Self.fixture().samples[0],
            assays: [.init(featureType: "Gene Expression", id: "rna", kind: .rna,
                featureNamespace: "genes", countUnit: .cutSiteCount, genomeAssembly: nil)])
        #expect(throws: (any Error).self) { try VivoMultiAssayTenX.validate(plan) }
        let base = MultiAssayH5MUImportTests.plan()
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(base)) as? [String: Any])
        var mappings = try #require(object["assays"] as? [[String: Any]])
        mappings[0]["countUnit"] = "cutSiteCount"; object["assays"] = mappings
        let mapped = try VivoCanonicalJSON.decode(VivoH5MUMultiAssayPlan.self, from: JSONSerialization.data(withJSONObject: object))
        #expect(throws: (any Error).self) { try VivoMultiAssayH5MUImport.validate(mapped) }
    }
    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1")) func nativeH5MUMapsPreservePartialRowsAndUInt64() throws {
        let path = FileManager.default.temporaryDirectory.appendingPathComponent("multiassay-test-" + UUID().uuidString + ".h5mu")
        defer { try? FileManager.default.removeItem(at: path) }
        try VivoMultiAssayH5MU.writeSnapshot(Self.fixture(), to: path)
        try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), f = try h.file(path.path); defer { h.close(f, "H5Fclose") }
            func integers(_ name: String) throws -> [UInt64] {
                let d = try h.dataset(f, name); defer { h.close(d, "H5Dclose") }
                return try h.integers(d, maximum: 10)
            }
            #expect(try integers("obsmap/protein") == [2,0,1])
            #expect(try integers("varmap/protein") == [0,1])
            #expect(try integers("mod/rna/X/data") == [9_007_199_254_740_993,4])
        }
        #expect(throws: (any Error).self) { try VivoMultiAssayH5MU.writeSnapshot(Self.fixture(), to: path) }
    }
}
