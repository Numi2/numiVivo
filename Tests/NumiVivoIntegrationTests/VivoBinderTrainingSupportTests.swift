import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderTrainingSupportTests: XCTestCase {
    typealias B = VivoBinderBenchmark
    typealias S = VivoBinderTrainingSupport
    let digest = String(repeating: "a", count: 64)
    let feature = "ipsae_min_boltz2"
    func dataset() -> B.Dataset {
        .init(sourceSHA256: digest, assay: "synthetic", groupingMethod: "synthetic exact groups",
            records: (0..<12).map { i in
                .init(id: "c\(i)", target: i < 4 ? "A" : i < 8 ? "B" : "TEST",
                    leakageGroup: "g\(i)", outcome: i % 2 == 0 ? .binder : .nonBinder,
                    rawOutcome: "synthetic", features: [feature: i % 2 == 0 ? 0.8 : 0.2])
            })
    }
    func plan() -> B.Plan {
        .init(sourceSHA256: digest, trainingTargets: ["A", "B"], testTargets: ["TEST"],
              baselineFeature: feature, modelFeatures: [feature])
    }
    func policy(positives: Int = 2, targets: Int = 2) -> S.Policy {
        .init(identifier: "synthetic-policy", rationale: "test fixture, not a calibrated biological criterion",
              minimumPositiveGroups: positives, minimumNegativeGroups: 2,
              minimumTrainingTargets: targets, minimumTargetsWithBothClasses: targets, maximumFeatures: 2)
    }
    func changing(_ dataset: B.Dataset, _ change: (B.Record) -> B.Record) -> B.Dataset {
        .init(sourceSHA256: dataset.sourceSHA256, assay: dataset.assay, groupingMethod: dataset.groupingMethod,
              records: dataset.records.map(change))
    }
    func testSupportCountsAndLegacyFitRemainIdentical() throws {
        let d = dataset(), p = plan()
        let e = try S.evaluate(d, plan: p, policy: policy())
        XCTAssertEqual(e.disposition, "experimentalFitOnly")
        XCTAssertEqual(e.support.overall.positiveGroups, 4)
        XCTAssertEqual(e.support.overall.negativeGroups, 4)
        XCTAssertEqual(e.support.byTarget.count, 2)
        XCTAssertEqual(try VivoCanonicalJSON.encode(XCTUnwrap(e.fittedReport)),
                       try VivoCanonicalJSON.encode(B.evaluate(d, plan: p)))
    }
    func testUnsupportedFitNeverRuns() throws {
        let e = try S.evaluate(dataset(), plan: plan(), policy: policy(positives: 5))
        XCTAssertEqual(e.disposition, "retainFixedRanking")
        XCTAssertNil(e.fittedReport)
        XCTAssertFalse(e.support.eligibleForExperimentalFit)
        XCTAssertEqual(e.support.unmetRequirements, ["insufficient exclusively positive training groups"])
    }
    func testSingleClassGetsSupportReportInsteadOfInvalidFit() throws {
        let d = changing(dataset()) { r in
            .init(id: r.id, target: r.target, leakageGroup: r.leakageGroup,
                  outcome: .nonBinder, rawOutcome: r.rawOutcome, features: r.features)
        }
        let e = try S.evaluate(d, plan: plan(), policy: policy())
        XCTAssertNil(e.fittedReport); XCTAssertEqual(e.support.overall.positiveGroups, 0)
    }
    func testHeldOutLabelsDoNotAffectSupport() throws {
        let d = dataset()
        let changed = changing(d) { r in
            .init(id: r.id, target: r.target, leakageGroup: r.leakageGroup,
                  outcome: r.target == "TEST" ? .notTested : r.outcome,
                  rawOutcome: r.rawOutcome, features: r.features)
        }
        XCTAssertEqual(try VivoCanonicalJSON.encode(S.assess(d, plan: plan(), policy: policy())),
                       try VivoCanonicalJSON.encode(S.assess(changed, plan: plan(), policy: policy())))
    }
    func testRepeatedAndConflictingGroupsAreNotExtraSupport() throws {
        let d = changing(dataset()) { r in
            .init(id: r.id, target: r.target, leakageGroup: r.target == "TEST" ? r.leakageGroup : "same",
                  outcome: r.outcome, rawOutcome: r.rawOutcome, features: r.features)
        }
        let s = try S.assess(d, plan: plan(), policy: policy())
        XCTAssertEqual(s.overall.mixedGroups, 1)
        XCTAssertEqual(s.overall.positiveGroups, 0)
        XCTAssertEqual(s.overall.negativeGroups, 0)
        XCTAssertFalse(s.eligibleForExperimentalFit)
    }
    func testLeakagePurgeAndMissingFeaturesCountBeforeSupport() throws {
        let d = changing(dataset()) { r in
            .init(id: r.id, target: r.target, leakageGroup: r.id == "c0" ? "g8" : r.leakageGroup,
                  outcome: r.outcome, rawOutcome: r.rawOutcome, features: r.id == "c2" ? [:] : r.features)
        }
        let s = try S.assess(d, plan: plan(), policy: policy())
        XCTAssertEqual(s.overall.positiveRows, 2)
        XCTAssertEqual(s.trainingExclusions["c0"], "training group overlaps test")
        XCTAssertEqual(s.trainingExclusions["c2"], "missing features:\(feature)")
        XCTAssertFalse(s.trainingExclusions.keys.contains("c8"))
        XCTAssertFalse(s.eligibleForExperimentalFit)
    }
    func testTargetAndFeaturePoliciesAreExplicit() throws {
        let s = try S.assess(dataset(), plan: plan(), policy: policy(targets: 3))
        XCTAssertEqual(s.unmetRequirements.count, 2)
        XCTAssertThrowsError(try S.assess(dataset(), plan: plan(), policy: policy(positives: 0)))
    }
    func testStrictPolicyAndPlanValidation() throws {
        let d = try VivoCanonicalJSON.encode(policy())
        _ = try S.decodePolicy(d)
        var object = try XCTUnwrap(JSONSerialization.jsonObject(with: d) as? [String: Any])
        object["testPositives"] = 4
        XCTAssertThrowsError(try S.decodePolicy(JSONSerialization.data(withJSONObject: object)))
        let overlap = B.Plan(sourceSHA256: digest, trainingTargets: ["A", "TEST"], testTargets: ["TEST"],
                             baselineFeature: feature, modelFeatures: [feature])
        XCTAssertThrowsError(try S.assess(dataset(), plan: overlap, policy: policy()))
    }
}
