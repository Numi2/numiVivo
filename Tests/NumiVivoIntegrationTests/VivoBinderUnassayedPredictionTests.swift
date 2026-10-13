import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderUnassayedPredictionTests: XCTestCase {
    typealias B = VivoBinderBenchmark
    let sha = String(repeating: "a", count: 64)
    func fixture() -> B.Dataset {
        let rows = (0..<8).map { i in B.Record(id: "t\(i)", target: "train", leakageGroup: "g\(i)",
            outcome: i.isMultiple(of: 2) ? .nonBinder : .binder,
            rawOutcome: i.isMultiple(of: 2) ? "negative" : "positive", features: ["score": Double(i % 2)]) }
        return B.Dataset(sourceSHA256: sha, assay: "fixture", groupingMethod: "fixture", records: rows)
    }
    func plan() -> B.Plan { .init(sourceSHA256: sha, trainingTargets: ["train"], testTargets: ["test"],
        baselineFeature: "score", modelFeatures: ["score"]) }
    func query(_ id: String = "q", group: String = "heldout", features: [String: Double] = ["score": 0.75]) -> B.QueryCandidate {
        .init(id: id, target: "test", leakageGroup: group, features: features)
    }
    func testUnassayedScoresMatchSharedFitterWithoutRetrospectiveOutcomes() throws {
        let training = fixture(), q = query()
        let prediction = try B.predict(training: training, query: [q], plan: plan())
        let heldout = B.Record(id: q.id, target: q.target, leakageGroup: q.leakageGroup,
            outcome: .binder, rawOutcome: "heldout", features: q.features)
        let dataset = B.Dataset(sourceSHA256: sha, assay: "fixture", groupingMethod: "fixture",
            records: training.records + [heldout])
        let retrospective = try B.evaluate(dataset, plan: plan())
        XCTAssertEqual(prediction.predictions[0].score, retrospective.targets[0].modelProbabilities[0])
        XCTAssertEqual(prediction.trainingIDs, retrospective.trainingIDs)
        XCTAssertEqual(prediction.trainingSHA256, try VivoBinderRanking.fingerprint(training))
        XCTAssertEqual(prediction.evidenceStatus, "unassayed-query-development-only")
    }
    func testQueryOrderDoesNotChangeIdentityOrPredictions() throws {
        let q = [query("z"), query("a", features: ["score": 0.25])]
        let a = try B.predict(training: fixture(), query: q, plan: plan())
        let b = try B.predict(training: fixture(), query: Array(q.reversed()), plan: plan())
        XCTAssertEqual(a.querySHA256, b.querySHA256)
        XCTAssertEqual(try VivoCanonicalJSON.encode(a), try VivoCanonicalJSON.encode(b))
        XCTAssertEqual(a.predictions.map(\.id), ["a", "z"])
    }
    func testUnsupportedQueryStillPurgesOverlappingTrainingGroup() throws {
        let r = try B.predict(training: fixture(), query: [query(group: "g0", features: [:])], plan: plan())
        XCTAssertFalse(r.trainingIDs.contains("t0"))
        XCTAssertEqual(r.excluded["t0"], "training group overlaps test")
        XCTAssertEqual(r.excluded["q"], "missing features:score")
        XCTAssertTrue(r.predictions.isEmpty)
    }
    func testDuplicateCrossSplitIdentityNonfiniteFeaturesAndUnknownTargetsReject() throws {
        XCTAssertThrowsError(try B.predict(training: fixture(), query: [query("t0")], plan: plan()))
        XCTAssertThrowsError(try B.predict(training: fixture(), query: [query(), query()], plan: plan()))
        XCTAssertThrowsError(try B.predict(training: fixture(), query: [query(features: ["score": .nan])], plan: plan()))
        XCTAssertThrowsError(try B.predict(training: fixture(), query: [], plan: plan()))
        let q = B.QueryCandidate(id: "q", target: "unknown", leakageGroup: "q", features: ["score": 1])
        XCTAssertThrowsError(try B.predict(training: fixture(), query: [q], plan: plan()))
    }
    func testTrainingCannotContainHeldoutTargetAndChangedQueryHasNewDigest() throws {
        let training = B.Dataset(sourceSHA256: sha, assay: "fixture", groupingMethod: "fixture",
            records: fixture().records + [.init(id: "q", target: "test", leakageGroup: "q", outcome: .binder,
                rawOutcome: "heldout", features: ["score": 1])])
        XCTAssertThrowsError(try B.predict(training: training, query: [query()], plan: plan()))
        let a = try B.predict(training: fixture(), query: [query()], plan: plan())
        let b = try B.predict(training: fixture(), query: [query(features: ["score": 0.25])], plan: plan())
        XCTAssertNotEqual(a.querySHA256, b.querySHA256)
        XCTAssertEqual(a.trainingSHA256, b.trainingSHA256)
    }
}
