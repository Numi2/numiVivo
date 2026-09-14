import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct QuantitativeH5MUTests {
    private static func fixture() -> VivoQuantitativeAssayDataset {
        let sample = VivoOmicsSample(id: "sample", biologicalReplicateID: "rep",
                                     condition: "control", batchID: "batch", organism: "NCBITaxon:9606")
        let observations: [VivoAssayObservation] = [
            .init(identity: .init(barcode: "b0", sampleID: "sample", group: nil), kind: .cell, position: nil),
            .init(identity: .init(barcode: "b1", sampleID: "sample", group: "treated"), kind: .cell, position: nil),
            .init(identity: .init(barcode: "b2", sampleID: "sample", group: nil), kind: .spot,
                  position: .init(frameID: "slide", coordinates: [4, 8]))]
        let protein = VivoQuantitativeAssaySpace(
            id: "protein", kind: .proteomics, featureNamespace: "uniprot", unit: "log2-intensity",
            sourceDescription: "Measured protein panel",
            features: [.init(id: "P1", name: "protein-1", interval: .init(contig: "chr1", start: 5, end: 8))], observationIndices: [0, 2],
            matrix: .init(observationCount: 2, featureCount: 1, rowOffsets: [0, 1, 2],
                          featureIndices: [0, 0], values: [0, 1.25]))
        let metabolite = VivoQuantitativeAssaySpace(
            id: "metabolite", kind: .metabolomics, featureNamespace: "hmdb", unit: "log10-concentration",
            sourceDescription: "Measured metabolite panel",
            features: [.init(id: "HMDB1", name: "metabolite-1", interval: nil)], observationIndices: [1],
            matrix: .init(observationCount: 1, featureCount: 1, rowOffsets: [0, 1],
                          featureIndices: [0], values: [-2.5]))
        return .init(id: "quantitative-h5mu", evidence: .measured,
                     sourceDescription: "Bounded quantitative H5MU fixture", samples: [sample],
                     observations: observations,
                     spatialFrames: [.init(id: "slide", unit: .pixel, axes: ["x", "y"], sourceDescription: "Measured image frame")],
                     assays: [protein, metabolite])
    }

    private static func plan(for data: VivoQuantitativeAssayDataset) -> VivoH5MUQuantitativePlan {
        .init(id: data.id, evidence: data.evidence, sourceDescription: data.sourceDescription,
              samples: data.samples, sampleColumn: "sample", barcodeColumn: "barcode", groupColumn: "group",
              defaultObservationKind: .cell, observationKindColumn: "observation_kind",
              spatial: [.init(path: "obsm/spatial", frame: data.spatialFrames[0])],
              assays: data.assays.map { assay in
                  .init(sourceName: assay.id, id: assay.id, kind: assay.kind,
                        featureNamespace: assay.featureNamespace, unit: assay.unit,
                        sourceDescription: assay.sourceDescription, matrixPath: "X",
                        featureNameColumn: "name", intervalContigColumn: "interval_contig",
                        intervalStartColumn: "interval_start", intervalEndColumn: "interval_end",
                        intervalPresentColumn: "interval_present")
              })
    }

    @Test func planRejectsUnknownFieldsAndIncompleteIntervalMappings() throws {
        let data = Self.fixture(), plan = Self.plan(for: data)
        try VivoQuantitativeH5MUImport.validate(plan)
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(plan)) as? [String: Any])
        object["ignoreUnknown"] = true
        #expect(throws: (any Error).self) {
            try VivoCanonicalJSON.decode(VivoH5MUQuantitativePlan.self, from: JSONSerialization.data(withJSONObject: object))
        }
        let incomplete = VivoH5MUQuantitativeAssayMapping(sourceName: "protein", id: "protein", kind: .proteomics,
            featureNamespace: "uniprot", unit: "log2-intensity", sourceDescription: "Measured protein panel",
            intervalContigColumn: "contig")
        let invalid = VivoH5MUQuantitativePlan(id: data.id, evidence: data.evidence,
            sourceDescription: data.sourceDescription, samples: data.samples, sampleColumn: "sample", assays: [incomplete])
        #expect(throws: (any Error).self) { try VivoQuantitativeH5MUImport.validate(invalid) }
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1"))
    func nativeH5MURoundTripPreservesRealValuesMissingZerosGroupsAndSpatialPositions() throws {
        let data = Self.fixture(), plan = Self.plan(for: data)
        let url = FileManager.default.temporaryDirectory.appendingPathComponent("quantitative-h5mu-" + UUID().uuidString + ".h5mu")
        defer { try? FileManager.default.removeItem(at: url) }
        try VivoQuantitativeH5MU.writeSnapshot(data, to: url)
        let result = try VivoQuantitativeH5MUImport.readSnapshot(url, plan: plan)
        #expect(result == data)
        #expect(try result.value(assayID: "protein", observationIndex: 0, featureID: "P1") == 0)
        #expect(try result.value(assayID: "protein", observationIndex: 1, featureID: "P1") == nil)
        #expect(try result.value(assayID: "protein", observationIndex: 2, featureID: "P1") == 1.25)
        #expect(result.observations[1].identity.group == "treated")
        #expect(result.observations[2].position == data.observations[2].position)
        try VivoHDF5.lock.withLock {
            let h = try VivoHDF5(), file = try h.file(url.path); defer { h.close(file, "H5Fclose") }
            let values = try h.dataset(file, "mod/protein/X/data"); defer { h.close(values, "H5Dclose") }
            #expect(try h.doubles(values, maximum: 4, allowInteger: false) == [0, 1.25])
        }
    }

    @Test(.enabled(if: ProcessInfo.processInfo.environment["NUMIVIVO_TEST_HDF5"] == "1"))
    func nativeH5MUImportAcceptsDenseAndCSCFloatMatrices() throws {
        let data = Self.fixture(), plan = Self.plan(for: data)
        for representation in ["dense", "csc"] {
            let url = FileManager.default.temporaryDirectory.appendingPathComponent("quantitative-\(representation)-" + UUID().uuidString + ".h5mu")
            defer { try? FileManager.default.removeItem(at: url) }
            try VivoQuantitativeH5MU.writeSnapshot(data, to: url)
            try VivoHDF5.lock.withLock {
                let h = try VivoHDF5(), file = try h.file(url.path, writable: true)
                defer { h.close(file, "H5Fclose") }
                let delete: @convention(c) (Int64, UnsafePointer<CChar>, Int64) -> Int32 = try h.symbol("H5Ldelete")
                try h.check(delete(file, "mod/protein/X", 0), "replace quantitative matrix")
                let modality = try h.object(file, "mod/protein"); defer { h.close(modality, "H5Oclose") }
                if representation == "dense" {
                    let matrix = try h.projectionDataset(modality, "X", type: h.native("NATIVE_DOUBLE"), shape: [2, 1])
                    defer { h.close(matrix, "H5Dclose") }
                    let write: @convention(c) (Int64, Int64, Int64, Int64, Int64, UnsafeRawPointer?) -> Int32 = try h.symbol("H5Dwrite")
                    let values = [0.0, 1.25]
                    try values.withUnsafeBytes { try h.check(write(matrix, h.native("NATIVE_DOUBLE"), 0, 0, 0, $0.baseAddress), "write quantitative dense fixture") }
                    try h.encoding(matrix, "array", "0.2.0")
                } else {
                    let matrix = try h.projectionGroup(modality, "X"); defer { h.close(matrix, "H5Gclose") }
                    try h.encoding(matrix, "csc_matrix", "0.1.0")
                    try h.writeIntegers(matrix, "shape", [2, 1], attribute: true)
                    try h.writeIndices(matrix, "indptr", [0, 2])
                    try h.writeIndices(matrix, "indices", [0, 1])
                    try h.writeDoubles(matrix, "data", [0, 1.25])
                }
            }
            #expect(try VivoQuantitativeH5MUImport.readSnapshot(url, plan: plan) == data)
        }
    }
}
