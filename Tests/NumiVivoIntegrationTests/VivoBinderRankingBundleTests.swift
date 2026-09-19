import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderRankingBundleTests: XCTestCase {
    typealias IO = VivoBinderBundleIO
    let identity = String(repeating: "c", count: 64)
    func workspace(_ body: (URL, URL, URL) throws -> Void) throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("source.csv"), config = root.appendingPathComponent("config.json")
        try Data("uuid,target,sequence,adaptyv_binding,ipsae_min_boltz2\na,T,AAAA,not_tested,0.9\nb,T,CCCC,binder,0.8\nc,T,DDDD,non_binder,0.1\n".utf8).write(to: source)
        try VivoCanonicalJSON.encode(VivoBinderAnthropicImport.Configuration(assay: .adaptyv, targets: ["T"],
            features: ["ipsae_min_boltz2"])).write(to: config)
        let imported = root.appendingPathComponent("imported"), query = root.appendingPathComponent("query")
        _ = try IO.importCSV(source: source, configuration: config, to: imported, implementationSHA256: identity)
        _ = try IO.rankingQuery(bundle: imported, to: query, implementationSHA256: identity)
        try body(root, imported, query)
    }
    func plan(_ root: URL, extra: Bool = false) throws -> URL {
        let url = root.appendingPathComponent("ranking-plan.json")
        var object: [String: Any] = ["schemaVersion": 1, "method": "fixedScore", "topK": 1,
            "components": [["feature": "ipsae_min_boltz2", "weight": 1]]]
        if extra { object["outcomes"] = ["binder"] }
        try JSONSerialization.data(withJSONObject: object).write(to: url)
        return url
    }
    func testQueryHasNoExperimentalPayload() throws {
        try workspace { _, _, query in
            let receipt = try IO.verify(query, implementationSHA256: identity)
            XCTAssertEqual(Set(receipt.files.keys), ["query.json"])
            let q = try VivoBinderRanking.decodeQuery(Data(contentsOf: query.appendingPathComponent("query.json")))
            XCTAssertEqual(q.candidates.count, 3)
        }
    }
    func testRankThenAssessAndReplay() throws {
        try workspace { root, imported, query in
            let ranked = root.appendingPathComponent("ranked"), assessed = root.appendingPathComponent("assessed")
            _ = try IO.rank(bundle: query, plan: plan(root), to: ranked, implementationSHA256: identity)
            _ = try IO.assessRanking(bundle: imported, ranking: ranked, to: assessed, implementationSHA256: identity)
            XCTAssertEqual(try IO.verify(ranked, implementationSHA256: identity).kind, "ranking")
            XCTAssertEqual(try IO.verify(assessed, implementationSHA256: identity).kind, "rankingAssessment")
            let report = try VivoCanonicalJSON.decode(VivoBinderRanking.Assessment.self,
                from: Data(contentsOf: assessed.appendingPathComponent("assessment.json")))
            XCTAssertEqual(report.targets[0].measuredBinderHits, 0)
            XCTAssertEqual(report.targets[0].selectedNonbinaryWeight, 1)
        }
    }
    func testUnknownPlanFieldsNeverPublish() throws {
        try workspace { root, _, query in
            let output = root.appendingPathComponent("bad")
            XCTAssertThrowsError(try IO.rank(bundle: query, plan: plan(root, extra: true), to: output, implementationSHA256: identity))
            XCTAssertFalse(FileManager.default.fileExists(atPath: output.path))
        }
    }
    func testWrongBundleKindAndOverwriteRejected() throws {
        try workspace { root, imported, query in
            let p = try plan(root), output = root.appendingPathComponent("ranking")
            XCTAssertThrowsError(try IO.rank(bundle: imported, plan: p, to: output, implementationSHA256: identity))
            _ = try IO.rank(bundle: query, plan: p, to: output, implementationSHA256: identity)
            XCTAssertThrowsError(try IO.rank(bundle: query, plan: p, to: output, implementationSHA256: identity))
            XCTAssertThrowsError(try IO.assessRanking(bundle: imported, ranking: query,
                to: root.appendingPathComponent("bad"), implementationSHA256: identity))
        }
    }
    func testRehashedChangedRankingStillRejected() throws {
        try workspace { root, _, query in
            let output = root.appendingPathComponent("ranking")
            _ = try IO.rank(bundle: query, plan: plan(root), to: output, implementationSHA256: identity)
            let path = output.appendingPathComponent("ranking.json")
            var object = try XCTUnwrap(JSONSerialization.jsonObject(with: Data(contentsOf: path)) as? [String: Any])
            object["methodVersion"] = "invented"
            let changed = try JSONSerialization.data(withJSONObject: object); try changed.write(to: path)
            let receiptPath = output.appendingPathComponent("receipt.json")
            var receipt = try XCTUnwrap(JSONSerialization.jsonObject(with: Data(contentsOf: receiptPath)) as? [String: Any])
            var hashes = receipt["files"] as! [String: String]
            hashes["ranking.json"] = try VivoCanonicalJSON.fingerprint(changed).hex; receipt["files"] = hashes
            try JSONSerialization.data(withJSONObject: receipt).write(to: receiptPath)
            XCTAssertThrowsError(try IO.verify(output, implementationSHA256: identity))
        }
    }
    func testNoExtraFilesInOutcomeBlindQuery() throws {
        try workspace { _, _, query in
            try Data("binder".utf8).write(to: query.appendingPathComponent("outcomes.txt"))
            XCTAssertThrowsError(try IO.verify(query, implementationSHA256: identity))
        }
    }
}
