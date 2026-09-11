import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct MultiAssayVisiumTests {
    static let header = "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres"
    static let rows = "b2,1,2,4,17,23\nb0,1,0,0,13,19\nb1,1,0,2,15,21\noutside,0,4,4,70,90\n"
    @Test func bothPublisherFormatsJoinByIdentityAndPreserveAxes() throws {
        let legacy = try VivoMultiAssayVisium.positions(Data(Self.rows.utf8), format: .legacyHeaderlessCSV)
        let modern = try VivoMultiAssayVisium.positions(Data((Self.header + "\r\n" + Self.rows.replacingOccurrences(of: "\n", with: "\r\n")).utf8), format: .csvWithHeader)
        #expect(legacy == modern)
        #expect(legacy["outside"]?.inTissue == false)
        #expect(legacy["b2"]?.arrayRow == 2 && legacy["b2"]?.arrayColumn == 4)
        let source = MultiAssayTests.fixture()
        let frame = VivoSpatialFrame(id: "image", unit: .pixel, axes: ["x", "y"], sourceDescription: "Test pixels")
        let attached = try VivoMultiAssayVisium.attach(legacy, to: source, frame: frame)
        #expect(attached.assays == source.assays)
        #expect(attached.observations.map(\.kind) == [.spot, .spot, .spot])
        #expect(attached.observations.map { $0.position!.coordinates } == [[19,13], [21,15], [23,17]])
    }
    @Test func malformedRowsAndAmbiguousFormatsReject() {
        for text in ["", Self.header + "\n" + Self.rows, Self.rows + "b2,1,0,0,1,2\n",
                     "b,2,0,0,1,2", "b,1,0,0,NaN,2", "b,1,0,0,1.5,2", "b,1,0,0,-1,2",
                     "b,1,0,0,9007199254740993,2", "b,1,0,0,1,2,3", "b,1,0,0,1,2\n\n",
                     "b,1,0,0,+1,2", "b,1,0,0,1,", "b,1,0,0,1,2\rb,1,0,0,1,2"] {
            #expect(throws: (any Error).self) { try VivoMultiAssayVisium.positions(Data(text.utf8), format: .legacyHeaderlessCSV) }
        }
        #expect(throws: (any Error).self) { try VivoMultiAssayVisium.positions(Data(Self.rows.utf8), format: .csvWithHeader) }
        #expect(throws: (any Error).self) { try VivoMultiAssayVisium.positions(Data([0xff]), format: .legacyHeaderlessCSV) }
    }
    @Test func missingExtraOrOutsideTissueMatrixRowsReject() throws {
        let source = MultiAssayTests.fixture(), frame = source.spatialFrames[0]
        for text in [Self.rows.replacingOccurrences(of: "b1,1", with: "b1,0"),
                     Self.rows.replacingOccurrences(of: "b1,1", with: "absent,1"),
                     Self.rows.replacingOccurrences(of: "outside,0", with: "outside,1")] {
            let table = try VivoMultiAssayVisium.positions(Data(text.utf8), format: .legacyHeaderlessCSV)
            #expect(throws: (any Error).self) { try VivoMultiAssayVisium.attach(table, to: source, frame: frame) }
        }
    }
    @Test func explicitPlanRequiresPixelAxesAndRejectsUnknownFields() throws {
        let counts = VivoTenXMultiAssayPlan(schemaVersion: 1, id: "visium", evidence: .synthetic, sourceDescription: "Test input",
            sample: MultiAssayTests.fixture().samples[0], assays: [.init(featureType: "Gene Expression", id: "rna", kind: .rna,
                featureNamespace: "test-gene", countUnit: .umiCount, genomeAssembly: "GRCh38")])
        func plan(_ unit: VivoSpatialUnit, _ axes: [String]) -> VivoVisiumPlan {
            .init(schemaVersion: 1, counts: counts, positionsFormat: .legacyHeaderlessCSV,
                  frame: .init(id: "image", unit: unit, axes: axes, sourceDescription: "Full-resolution pixels"))
        }
        try VivoMultiAssayVisium.validate(plan(.pixel, ["x", "y"]))
        for p in [plan(.micrometer, ["x", "y"]), plan(.pixel, ["y", "x"]), plan(.pixel, ["x", "y", "z"])] {
            #expect(throws: (any Error).self) { try VivoMultiAssayVisium.validate(p) }
        }
        var value = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(plan(.pixel, ["x", "y"]))) as? [String: Any])
        value["inferMicrometers"] = true
        #expect(throws: (any Error).self) { try VivoCanonicalJSON.decode(VivoVisiumPlan.self, from: JSONSerialization.data(withJSONObject: value)) }
    }
}
