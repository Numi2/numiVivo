import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderPredictionTests: XCTestCase {
    func input(sample: String = "seed-1", separation: Double = 4, predictor: String = "fixture", sequence: String = "G",
               points: Int = 96) throws -> VivoBinderPredictionAnalysis.Input {
        var text = "data_test\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n"
        for (r, chain) in ["B", "T"].enumerated() {
            for (a, name) in ["N", "CA", "C", "O"].enumerated() {
                let element = name == "N" ? "N" : name == "O" ? "O" : "C"
                text += "ATOM \(r*4+a+1) \(element) \(name) GLY \(chain) 1 \(a) \(Double(r)*separation) 0 1\n"
            }
        }
        return .init(identity: .init(candidateID: "candidate", predictor: predictor, sampleID: sample, targetConstructID: "one-gly"),
            binderSequence: sequence, source: .init(candidateID: "candidate", target: "target", sourceLabel: "synthetic fixture",
                format: .mmcif, contents: text, sha256: try VivoCanonicalJSON.fingerprint(Data(text.utf8)).hex,
                targetChainSequences: ["T": "G"], interfacePlan: .init(binderChains: ["B"], targetChains: ["T"], conformerID: "model-1")),
            surfacePlan: .init(radiusProfile: "test-explicit", radiiNM: ["C": 0.17, "N": 0.155, "O": 0.152], pointsPerAtom: points))
    }
    func testSourceAnalysisAndSequenceAdmission() throws {
        let r = try VivoBinderPredictionAnalysis.analyze(input())
        XCTAssertEqual(r.contacts.count, 1)
        XCTAssertEqual(r.contacts[0].binderOffset, 0)
        XCTAssertEqual(r.contacts[0].targetOffset, 0)
        XCTAssertGreaterThan(r.features["numi.surface.buriedAreaSumNM2"]!, 0)
        XCTAssertThrowsError(try VivoBinderPredictionAnalysis.analyze(input(sequence: "A")))
    }
    func testUnknownOutcomeFieldRejected() throws {
        let data = try VivoCanonicalJSON.encode(input())
        var object = try JSONSerialization.jsonObject(with: data) as! [String: Any]
        object["outcome"] = "binder"
        XCTAssertThrowsError(try VivoBinderPredictionAnalysis.decodeInput(JSONSerialization.data(withJSONObject: object)))
        object.removeValue(forKey: "outcome")
        var source = object["source"] as! [String: Any]; source["affinity"] = 1; object["source"] = source
        XCTAssertThrowsError(try VivoBinderPredictionAnalysis.decodeInput(JSONSerialization.data(withJSONObject: object)))
    }
    func testSeparatedPartnersRetainZeroAreaAndOmitUndefinedRatio() throws {
        let r = try VivoBinderPredictionAnalysis.analyze(input(separation: 30))
        XCTAssertEqual(r.features["numi.surface.buriedAreaSumNM2"], 0)
        XCTAssertNil(r.features["numi.surface.contactResiduePairsPerNM2"])
    }
    func testMultipleSeedsAndPredictorsRemainDistinct() throws {
        var a = VivoBinderPredictionAnalysis.Accumulator()
        try a.append(input()); try a.append(input(sample: "seed-2", separation: 5))
        try a.append(input(predictor: "other"))
        let groups = try a.finish()
        XCTAssertEqual(groups.count, 2)
        let g = groups.first { $0.predictor == "fixture" }!
        XCTAssertEqual(g.sampleIDs, ["seed-1", "seed-2"])
        XCTAssertEqual(g.meanContactJaccard, 0)
        XCTAssertEqual(g.features["numi.surface.buriedAreaSumNM2"]?.availableCount, 2)
        XCTAssertGreaterThan(g.features["numi.surface.buriedAreaSumNM2"]!.populationSD, 0)
    }
    func testDuplicateAppendLeavesStateUnchanged() throws {
        var a = VivoBinderPredictionAnalysis.Accumulator(); try a.append(input())
        let before = try VivoCanonicalJSON.encode(a.finish())
        XCTAssertThrowsError(try a.append(input(separation: 10)))
        XCTAssertEqual(before, try VivoCanonicalJSON.encode(a.finish()))
    }
    func testEmptyContactsAreNotPerfectAgreement() throws {
        var a = VivoBinderPredictionAnalysis.Accumulator()
        try a.append(input(separation: 30)); try a.append(input(sample: "seed-2", separation: 35))
        let g = try a.finish()[0]
        XCTAssertNil(g.meanContactJaccard); XCTAssertEqual(g.bothEmptyContactComparisons, 1)
    }
    func testInconsistentSurfaceProfilesRejectAggregation() throws {
        var a = VivoBinderPredictionAnalysis.Accumulator()
        try a.append(input()); try a.append(input(sample: "seed-2", points: 128))
        XCTAssertThrowsError(try a.finish())
    }
    func testInputOrderDoesNotChangeSummary() throws {
        var a = VivoBinderPredictionAnalysis.Accumulator(), b = VivoBinderPredictionAnalysis.Accumulator()
        let x = try input(), y = try input(sample: "seed-2", separation: 5)
        try a.append(x); try a.append(y); try b.append(y); try b.append(x)
        XCTAssertEqual(try VivoCanonicalJSON.encode(a.finish()), try VivoCanonicalJSON.encode(b.finish()))
    }
    func testBundleRecomputesAndRejectsTampering() throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let file = root.appendingPathComponent("prediction.json"), out = root.appendingPathComponent("result")
        try VivoCanonicalJSON.encode(input()).write(to: file)
        let identity = String(repeating: "c", count: 64)
        _ = try VivoBinderBundleIO.analyzePrediction(input: file, to: out, implementationSHA256: identity)
        XCTAssertEqual(try VivoBinderBundleIO.verify(out, implementationSHA256: identity).kind, "predictionAnalysis")
        XCTAssertThrowsError(try VivoBinderBundleIO.analyzePrediction(input: file, to: out, implementationSHA256: identity))
        let report = out.appendingPathComponent("analysis.json")
        var object = try JSONSerialization.jsonObject(with: Data(contentsOf: report)) as! [String: Any]
        object["features"] = ["invented": 42]
        let data = try JSONSerialization.data(withJSONObject: object, options: [.sortedKeys]); try data.write(to: report)
        let receipt = out.appendingPathComponent("receipt.json")
        var r = try JSONSerialization.jsonObject(with: Data(contentsOf: receipt)) as! [String: Any]
        var files = r["files"] as! [String: String]; files["analysis.json"] = try VivoCanonicalJSON.fingerprint(data).hex
        r["files"] = files; try JSONSerialization.data(withJSONObject: r).write(to: receipt)
        XCTAssertThrowsError(try VivoBinderBundleIO.verify(out, implementationSHA256: identity))
    }
}
