import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderImportTests: XCTestCase {
    typealias I = VivoBinderAnthropicImport
    let sha = String(repeating: "a", count: 64)
    let header = "uuid,target,sequence,adaptyv_binding,twist_binding,ipsae_min_boltz2,adaptyv_kd_nM\n"
    func config(_ assay: I.Assay = .adaptyv) -> I.Configuration {
        I.Configuration(assay: assay, targets: ["MBP"], features: ["ipsae_min_boltz2"])
    }
    func parse(_ rows: String, _ assay: I.Assay = .adaptyv) throws -> I.Result {
        try I.parse(Data((header + rows).utf8), sourceSHA256: sha, configuration: config(assay))
    }
    func testAssaysRemainSeparate() throws {
        let rows = "a,MBP,ACDE,binder,non_binder,0.8,20\n"
        XCTAssertEqual(try parse(rows).dataset.records[0].outcome, .binder)
        XCTAssertEqual(try parse(rows, .twist).dataset.records[0].outcome, .nonBinder)
        XCTAssertEqual(try parse(rows).dataset.records[0].sourceFields["adaptyv_kd_nM"], "20")
    }
    func testUntestedNotNegativeAndUnknownPreserved() throws {
        let r = try parse("a,MBP,ACDE,not_tested,non_binder,0,\nb,MBP,ACDF,weak_signal,non_binder,,\n")
        XCTAssertEqual(r.dataset.records[0].outcome, .notTested)
        XCTAssertEqual(r.dataset.records[0].features["ipsae_min_boltz2"], 0)
        XCTAssertEqual(r.dataset.records[1].outcome, .inconclusive)
        XCTAssertNil(r.dataset.records[1].features["ipsae_min_boltz2"])
        XCTAssertEqual(r.unknownOutcomeCounts["weak_signal"], 1)
    }
    func testAssayFieldCannotBecomePredictor() {
        let c = I.Configuration(assay: .adaptyv, targets: ["MBP"], features: ["adaptyv_kd_nM"])
        XCTAssertThrowsError(try I.parse(Data(header.utf8), sourceSHA256: sha, configuration: c))
    }
    func testTargetScopeAndSourceRowWitness() throws {
        let r = try parse("skip,OTHER,ACDE,binder,binder,0.8,\na,MBP,ACDF,binder,binder,0.8,\n")
        XCTAssertEqual(r.dataset.records.count, 1); XCTAssertEqual(r.excludedTargetCounts["OTHER"], 1)
        XCTAssertEqual(r.sourceRowByID["a"], 3)
    }
    func testCSVQuotesCRLFAndUnicode() throws {
        let text = "a,b\r\n\"one, two\",\"a\"\"b\nµ\"\r\n"
        XCTAssertEqual(try I.csv(Data(text.utf8)), [["a", "b"], ["one, two", "a\"b\nµ"]])
    }
    func testCSVFinalEmptyAndBOM() throws {
        XCTAssertEqual(try I.csv(Data("\u{feff}a,b\n1,".utf8)), [["a", "b"], ["1", ""]])
    }
    func testMalformedCSVRejected() {
        for text in ["a\"b,c", "\"a\"x,b", "\"unclosed", "a,b\rc,d"] {
            XCTAssertThrowsError(try I.csv(Data(text.utf8)))
        }
        XCTAssertThrowsError(try I.csv(Data([0xff])))
    }
    func testDuplicateAndRaggedHeaderRejected() {
        XCTAssertThrowsError(try parse("a,MBP,ACDE,binder,binder,0.8\n"))
        XCTAssertThrowsError(try I.parse(Data("uuid,uuid\na,b\n".utf8), sourceSHA256: sha, configuration: config()))
    }
    func testInvalidNormalizedScoresRejected() {
        for value in ["NaN", "inf", "-0.1", "1.1", "missing"] {
            XCTAssertThrowsError(try parse("a,MBP,ACDE,binder,binder,\(value),\n"))
        }
    }
    func testSequenceIdentityAndDuplicateCandidateGuard() throws {
        let r = try parse("a,MBP,ACDE,binder,binder,0.8,\nb,MBP,ACDE,binder,binder,0.9,\n")
        XCTAssertEqual(r.dataset.records[0].leakageGroup, r.dataset.records[1].leakageGroup)
        XCTAssertThrowsError(try parse("a,MBP,ACDE,binder,binder,0.8,\na,MBP,ACDF,binder,binder,0.9,\n"))
    }
    func testMissingTargetAndUnsupportedSequenceRejected() {
        XCTAssertThrowsError(try parse("a,OTHER,ACDE,binder,binder,0.8,\n"))
        XCTAssertThrowsError(try parse("a,MBP,AXDE,binder,binder,0.8,\n"))
    }
}
