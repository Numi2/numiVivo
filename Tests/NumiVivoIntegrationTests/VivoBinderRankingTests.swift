import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderRankingTests: XCTestCase {
    typealias R = VivoBinderRanking
    let digest = String(repeating: "a", count: 64)
    let f = "ipsae_min_boltz2", g = "ipsae_min_ptxv2"
    func query(_ values: [Double]) -> R.Query {
        .init(sourceSHA256: digest, candidates: values.enumerated().map {
            .init(id: "c\($0.offset)", target: "T", scores: [f: $0.element])
        })
    }
    func fixed(_ k: Int = 2) -> R.Plan { .init(method: .fixedScore, components: [.init(feature: f)], topK: k) }
    func dataset(_ outcomes: [VivoBinderBenchmark.Outcome]) -> VivoBinderBenchmark.Dataset {
        .init(sourceSHA256: digest, assay: "synthetic", groupingMethod: "synthetic",
            records: outcomes.enumerated().map {
                .init(id: "c\($0.offset)", target: "T", leakageGroup: "g\($0.offset)", outcome: $0.element,
                      rawOutcome: $0.element.rawValue, features: [f: 1 - Double($0.offset) / 10],
                      sourceFields: ["secret": "experimental values"])
            })
    }
    func testFixedRankingAndTies() throws {
        let report = try R.rank(query([0.2, 0.8, 0.8, 0.8]), plan: fixed())
        XCTAssertEqual(report.targets[0].candidates.map(\.id), ["c1", "c2", "c3", "c0"])
        XCTAssertEqual(report.targets[0].candidates[0].selectionWeight, 2.0 / 3.0)
        XCTAssertEqual(report.targets[0].candidates.reduce(0) { $0 + $1.selectionWeight }, 2, accuracy: 1e-12)
    }
    func testPopulationZScores() throws {
        let p = R.Plan(method: .targetStandardizedMean, components: [.init(feature: f)], topK: 2)
        let result = try R.rank(query([0, 0.5, 1]), plan: p).targets[0]
        XCTAssertEqual(result.normalization[0].mean, 0.5, accuracy: 1e-12)
        XCTAssertEqual(result.normalization[0].populationSD, sqrt(1.0 / 6.0), accuracy: 1e-12)
        XCTAssertEqual(result.candidates[0].score, sqrt(1.5), accuracy: 1e-12)
    }
    func testWeightedComponentsAndConstantHandling() throws {
        let q = R.Query(sourceSHA256: digest, candidates: [
            .init(id: "a", target: "T", scores: [f: 0, g: 0.4]),
            .init(id: "b", target: "T", scores: [f: 1, g: 0.4])])
        let p = R.Plan(method: .targetStandardizedMean,
                      components: [.init(feature: f, weight: 3), .init(feature: g)], topK: 1)
        let t = try R.rank(q, plan: p).targets[0]
        XCTAssertTrue(t.normalization[1].constant)
        XCTAssertEqual(t.candidates[0].score, 0.75, accuracy: 1e-12)
    }
    func testAllConstantScoresSelectFractionally() throws {
        let p = R.Plan(method: .targetStandardizedMean, components: [.init(feature: f)], topK: 1)
        let t = try R.rank(query([0.5, 0.5, 0.5]), plan: p).targets[0]
        XCTAssertTrue(t.candidates.allSatisfy { $0.score == 0 && $0.selectionWeight == 1.0 / 3 })
    }
    func testMissingCandidateAndEntireMissingTargetAreExplicit() throws {
        let q = R.Query(sourceSHA256: digest, candidates: [
            .init(id: "a", target: "T", scores: [f: 1]), .init(id: "b", target: "T", scores: [:]),
            .init(id: "c", target: "U", scores: [:])])
        let r = try R.rank(q, plan: fixed())
        XCTAssertEqual(r.targets[0].excluded["b"], [f])
        XCTAssertEqual(r.targets[0].effectiveK, 1)
        XCTAssertEqual(r.targets[1].inputCount, 1)
        XCTAssertEqual(r.targets[1].effectiveK, 0)
        XCTAssertTrue(r.targets[1].candidates.isEmpty)
    }
    func testProjectionStripsOutcomesAndSourceFields() throws {
        let a = try R.query(from: dataset([.binder, .nonBinder, .notTested]))
        let b = try R.query(from: dataset([.nonBinder, .binder, .expressionFailure]))
        XCTAssertEqual(try VivoCanonicalJSON.encode(a), try VivoCanonicalJSON.encode(b))
        XCTAssertEqual(try VivoCanonicalJSON.encode(R.rank(a, plan: fixed())),
                       try VivoCanonicalJSON.encode(R.rank(b, plan: fixed())))
        let text = String(decoding: try VivoCanonicalJSON.encode(a), as: UTF8.self)
        XCTAssertFalse(text.contains("outcome")); XCTAssertFalse(text.contains("secret"))
    }
    func testQueryDecodeRejectsOutcomeFields() throws {
        let data = try VivoCanonicalJSON.encode(query([0.2]))
        _ = try R.decodeQuery(data)
        var object = try XCTUnwrap(JSONSerialization.jsonObject(with: data) as? [String: Any])
        object["outcome"] = "binder"
        XCTAssertThrowsError(try R.decodeQuery(JSONSerialization.data(withJSONObject: object)))
        object.removeValue(forKey: "outcome")
        var candidates = object["candidates"] as! [[String: Any]]
        candidates[0]["outcome"] = "binder"; object["candidates"] = candidates
        XCTAssertThrowsError(try R.decodeQuery(JSONSerialization.data(withJSONObject: object)))
    }
    func testUntestedSelectionDoesNotRefillTopK() throws {
        let d = dataset([.notTested, .nonBinder, .binder])
        let q = try R.query(from: d), p = fixed(1), r = try R.rank(q, plan: p)
        let a = try R.assess(d, query: q, plan: p, ranking: r).targets[0]
        XCTAssertEqual(a.measuredBinderHits, 0)
        XCTAssertEqual(a.selectedNonbinaryWeight, 1)
        XCTAssertNil(a.precisionAmongBinarySelected)
        XCTAssertEqual(a.binaryOutcomeCount, 2)
    }
    func testEvaluationPreservesEveryOutcomeCategory() throws {
        let d = dataset([.binder, .nonBinder, .notTested, .inconclusive, .expressionFailure])
        let q = try R.query(from: d), p = fixed(10), r = try R.rank(q, plan: p)
        let a = try R.assess(d, query: q, plan: p, ranking: r).targets[0]
        XCTAssertEqual(a.selectedWeightByOutcome.count, 5)
        XCTAssertEqual(a.measuredBinderHits, 1)
        XCTAssertEqual(a.selectedNonbinaryWeight, 3)
        XCTAssertEqual(a.precisionAmongBinarySelected, 0.5)
    }
    func testTamperedRankingAndUnrelatedQueryRejected() throws {
        let d = dataset([.binder, .nonBinder]), q = try R.query(from: d), p = fixed(1)
        let r = try R.rank(q, plan: p)
        XCTAssertThrowsError(try R.assess(d, query: query([0.4, 0.3]), plan: p, ranking: r))
        let changed = try R.rank(q, plan: fixed(2))
        XCTAssertThrowsError(try R.assess(d, query: q, plan: p, ranking: changed))
    }
    func testTargetsNormalizeSeparately() throws {
        let q = R.Query(sourceSHA256: digest, candidates: [
            .init(id: "a", target: "T", scores: [f: 0]), .init(id: "b", target: "T", scores: [f: 0.2]),
            .init(id: "c", target: "U", scores: [f: 0.7]), .init(id: "d", target: "U", scores: [f: 1])])
        let p = R.Plan(method: .targetStandardizedMean, components: [.init(feature: f)], topK: 1)
        for t in try R.rank(q, plan: p).targets { XCTAssertEqual(t.candidates[0].score, 1, accuracy: 1e-12) }
    }
    func testQueryOrderDoesNotChangeScoresOrSelection() throws {
        let q = query([0.1, 0.8, 0.4, 0.8])
        let r = R.Query(sourceSHA256: digest, candidates: Array(q.candidates.reversed()))
        XCTAssertEqual(try VivoCanonicalJSON.encode(R.rank(q, plan: fixed()).targets),
                       try VivoCanonicalJSON.encode(R.rank(r, plan: fixed()).targets))
    }
    func testInvalidPlanAndScoresReject() throws {
        for value in [Double.nan, Double.infinity, -0.1, 1.1] {
            XCTAssertThrowsError(try R.rank(query([value]), plan: fixed()))
        }
        XCTAssertThrowsError(try R.rank(query([0.5]), plan: fixed(0)))
        XCTAssertThrowsError(try R.rank(query([0.5]), plan: .init(method: .fixedScore, components: [])))
        XCTAssertThrowsError(try R.rank(query([0.5]), plan: .init(method: .fixedScore, components: [.init(feature: f, weight: 2)])))
        XCTAssertThrowsError(try R.rank(query([0.5]), plan: .init(method: .targetStandardizedMean, components: [.init(feature: f, weight: -1)])))
        XCTAssertThrowsError(try R.rank(query([0.5]), plan: .init(method: .targetStandardizedMean, components: [.init(feature: f), .init(feature: f)])))
        XCTAssertThrowsError(try R.rank(.init(sourceSHA256: "bad", candidates: query([0.5]).candidates), plan: fixed()))
        XCTAssertThrowsError(try R.rank(.init(sourceSHA256: digest, candidates: [.init(id: "x", target: "T", scores: ["measured_affinity": 1])]), plan: fixed()))
    }
    func testDuplicateCandidateRejected() throws {
        let c = query([0.5]).candidates[0]
        XCTAssertThrowsError(try R.rank(.init(sourceSHA256: digest, candidates: [c, c]), plan: fixed()))
    }
}
