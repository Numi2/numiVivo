import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderBenchmarkTests: XCTestCase {
    typealias B = VivoBinderBenchmark
    let sha = String(repeating: "a", count: 64)
    func fixture() -> B.Dataset {
        let rows = (0..<12).map { i in
            B.Record(id: "c\(i)", target: i < 8 ? "train" : "test", leakageGroup: "g\(i)",
                outcome: i % 2 == 0 ? .nonBinder : .binder, rawOutcome: i % 2 == 0 ? "non_binder" : "binder",
                features: ["score": i % 2 == 0 ? 0.1 : 0.9, "extra": Double(i % 3)])
        }
        return B.Dataset(sourceSHA256: sha, assay: "one-vendor:binding", groupingMethod: "synthetic-groups", records: rows)
    }
    func plan(_ features: [String] = ["score"], k: Int = 2) -> B.Plan {
        B.Plan(sourceSHA256: sha, trainingTargets: ["train"], testTargets: ["test"],
               baselineFeature: "score", modelFeatures: features, topK: k)
    }
    func replace(_ dataset: B.Dataset, rows: [B.Record]) -> B.Dataset {
        B.Dataset(sourceSHA256: dataset.sourceSHA256, assay: dataset.assay, groupingMethod: dataset.groupingMethod, records: rows)
    }
    func testPerfectRanking() throws {
        let m = try B.metrics(labels: [1, 0, 1, 0], scores: [0.9, 0.2, 0.8, 0.1], topK: 2)
        XCTAssertEqual(m.auroc!, 1); XCTAssertEqual(m.averagePrecision!, 1)
        XCTAssertEqual(m.precisionAtK, 1); XCTAssertNil(m.brier)
    }
    func testTiesDoNotDependOnInputOrder() throws {
        let a = try B.metrics(labels: [1, 0, 0, 1], scores: [1, 1, 1, 1], topK: 1)
        let b = try B.metrics(labels: [0, 1, 1, 0], scores: [1, 1, 1, 1], topK: 1)
        XCTAssertEqual(a.precisionAtK, 0.5); XCTAssertEqual(a.auroc!, 0.5)
        XCTAssertEqual(a.averagePrecision!, 0.5); XCTAssertEqual(a.expectedHitsAtK, b.expectedHitsAtK)
    }
    func testBoundaryTieFraction() throws {
        let m = try B.metrics(labels: [1, 1, 0, 0], scores: [3, 2, 2, 1], topK: 2)
        XCTAssertEqual(m.expectedHitsAtK, 1.5); XCTAssertEqual(m.auroc!, 0.875)
    }
    func testSingleClassUnavailableAUROC() throws {
        let m = try B.metrics(labels: [0, 0], scores: [0.5, 0.4], topK: 9)
        XCTAssertNil(m.auroc); XCTAssertNil(m.averagePrecision); XCTAssertEqual(m.effectiveK, 2)
    }
    func testBrierAndProbabilityBounds() throws {
        let m = try B.metrics(labels: [1, 0], scores: [0.5, 0.5], topK: 1, probabilities: true)
        XCTAssertEqual(m.brier!, 0.25)
        XCTAssertThrowsError(try B.metrics(labels: [1], scores: [2], topK: 1, probabilities: true))
    }
    func testBadMetricsRejected() {
        XCTAssertThrowsError(try B.metrics(labels: [], scores: [], topK: 1))
        XCTAssertThrowsError(try B.metrics(labels: [0], scores: [.nan], topK: 1))
        XCTAssertThrowsError(try B.metrics(labels: [2], scores: [0.5], topK: 1))
        XCTAssertThrowsError(try B.metrics(labels: [1], scores: [0.5], topK: 0))
    }
    func testNativeFit() throws {
        let r = try B.evaluate(fixture(), plan: plan())
        XCTAssertEqual(r.trainingIDs.count, 8); XCTAssertEqual(r.targets[0].learned.auroc!, 1)
        XCTAssertLessThanOrEqual(r.model.gradientInfinityNorm, 1e-7)
        XCTAssertEqual(r.model.means[0], 0.5, accuracy: 1e-12)
        XCTAssertEqual(r.evidenceStatus, "retrospective-development-only")
    }
    func testHeldOutLabelsCannotChangeModel() throws {
        let d = fixture(), before = try B.evaluate(d, plan: plan())
        let changed = d.records.map { r in
            B.Record(id: r.id, target: r.target, leakageGroup: r.leakageGroup,
                outcome: r.target == "test" ? .nonBinder : r.outcome,
                rawOutcome: r.rawOutcome, features: r.features)
        }
        let after = try B.evaluate(replace(d, rows: changed), plan: plan())
        XCTAssertEqual(before.model.coefficients, after.model.coefficients)
        XCTAssertEqual(before.model.means, after.model.means)
        XCTAssertEqual(before.targets[0].modelProbabilities, after.targets[0].modelProbabilities)
    }
    func testInputOrderDeterministic() throws {
        let d = fixture(), a = try B.evaluate(d, plan: plan(["score", "extra"]))
        let b = try B.evaluate(replace(d, rows: d.records.reversed()), plan: plan(["score", "extra"]))
        XCTAssertEqual(a.model.coefficients, b.model.coefficients); XCTAssertEqual(a.trainingIDs, b.trainingIDs)
    }
    func testTargetOverlapRejected() {
        let p = B.Plan(sourceSHA256: sha, trainingTargets: ["train"], testTargets: ["train"],
                       baselineFeature: "score", modelFeatures: ["score"])
        XCTAssertThrowsError(try B.evaluate(fixture(), plan: p))
    }
    func testSourceMismatchRejected() {
        let p = B.Plan(sourceSHA256: String(repeating: "b", count: 64), trainingTargets: ["train"],
                       testTargets: ["test"], baselineFeature: "score", modelFeatures: ["score"])
        XCTAssertThrowsError(try B.evaluate(fixture(), plan: p))
    }
    func testUnknownTargetRejected() {
        let p = B.Plan(sourceSHA256: sha, trainingTargets: ["train"], testTargets: ["absent"],
                       baselineFeature: "score", modelFeatures: ["score"])
        XCTAssertThrowsError(try B.evaluate(fixture(), plan: p))
    }
    func testGroupLeakagePurged() throws {
        let d = fixture(), first = d.records[0]
        var rows = d.records
        rows[0] = B.Record(id: first.id, target: first.target, leakageGroup: "g8", outcome: first.outcome,
                           rawOutcome: first.rawOutcome, features: first.features)
        let r = try B.evaluate(replace(d, rows: rows), plan: plan())
        XCTAssertEqual(r.excluded["c0"], "training group overlaps test"); XCTAssertEqual(r.trainingIDs.count, 7)
    }
    func testMissingAndUntestedStayExcluded() throws {
        let d = fixture(); var rows = d.records
        rows[8] = B.Record(id: "c8", target: "test", leakageGroup: "g8", outcome: .notTested,
                           rawOutcome: "not_tested", features: ["score": 1])
        rows[9] = B.Record(id: "c9", target: "test", leakageGroup: "g9", outcome: .binder,
                           rawOutcome: "binder", features: [:])
        let r = try B.evaluate(replace(d, rows: rows), plan: plan())
        XCTAssertEqual(r.targets[0].baseline.count, 2); XCTAssertEqual(r.targets[0].learned.count, 2)
        XCTAssertEqual(r.excluded.count, 2)
    }
    func testDuplicateIdentityRejected() {
        let d = fixture()
        XCTAssertThrowsError(try B.validate(replace(d, rows: d.records + [d.records[0]])))
    }
    func testUnknownSchemaRejectedAfterDecode() throws {
        let data = try JSONEncoder().encode(fixture())
        var json = try XCTUnwrap(JSONSerialization.jsonObject(with: data) as? [String: Any])
        json["schemaVersion"] = 2
        let decoded = try JSONDecoder().decode(B.Dataset.self, from: JSONSerialization.data(withJSONObject: json))
        XCTAssertThrowsError(try B.evaluate(decoded, plan: plan()))
    }
    func testBadFeatureRejected() {
        let d = fixture(); var rows = d.records
        rows[0] = B.Record(id: "c0", target: "train", leakageGroup: "g0", outcome: .binder,
                           rawOutcome: "binder", features: ["score": .infinity])
        XCTAssertThrowsError(try B.validate(replace(d, rows: rows)))
    }
    func testSingleClassTrainingRejected() {
        let d = fixture(), rows = d.records.map { r in
            B.Record(id: r.id, target: r.target, leakageGroup: r.leakageGroup, outcome: .nonBinder,
                     rawOutcome: "non_binder", features: r.features)
        }
        XCTAssertThrowsError(try B.evaluate(replace(d, rows: rows), plan: plan()))
    }
}
