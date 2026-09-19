import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderSupportBundleTests: XCTestCase {
    let identity = String(repeating: "c", count: 64)
    func workspace(_ body: (URL, URL, URL, URL) throws -> Void) throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("source.csv"), config = root.appendingPathComponent("config.json")
        var csv = "uuid,target,sequence,adaptyv_binding,ipsae_min_boltz2\n"
        let letters = Array("ACDEFGHIKLMN")
        for i in 0..<12 {
            csv += "c\(i),\(i < 4 ? "A" : i < 8 ? "B" : "TEST"),\(String(repeating: String(letters[i]), count: 4)),\(i % 2 == 0 ? "binder" : "non_binder"),\(i % 2 == 0 ? "0.8" : "0.2")\n"
        }
        try Data(csv.utf8).write(to: source)
        try VivoCanonicalJSON.encode(VivoBinderAnthropicImport.Configuration(assay: .adaptyv, targets: ["A", "B", "TEST"],
            features: ["ipsae_min_boltz2"])).write(to: config)
        let bundle = root.appendingPathComponent("input")
        _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: bundle, implementationSHA256: identity)
        let plan = root.appendingPathComponent("plan.json"), policy = root.appendingPathComponent("policy.json")
        try VivoCanonicalJSON.encode(VivoBinderBenchmark.Plan(sourceSHA256: VivoCanonicalJSON.fingerprint(Data(csv.utf8)).hex,
            trainingTargets: ["A", "B"], testTargets: ["TEST"], baselineFeature: "ipsae_min_boltz2",
            modelFeatures: ["ipsae_min_boltz2"])).write(to: plan)
        try writePolicy(policy, minimum: 2)
        try body(root, bundle, plan, policy)
    }
    func writePolicy(_ url: URL, minimum: Int) throws {
        try VivoCanonicalJSON.encode(VivoBinderTrainingSupport.Policy(identifier: "fixture", rationale: "synthetic test only",
            minimumPositiveGroups: minimum, minimumNegativeGroups: 2, minimumTrainingTargets: 2,
            minimumTargetsWithBothClasses: 2, maximumFeatures: 1)).write(to: url)
    }
    func testSupportedAndDeclinedBundlesReplay() throws {
        try workspace { root, bundle, plan, policy in
            for minimum in [2, 100] {
                try writePolicy(policy, minimum: minimum)
                let out = root.appendingPathComponent("output-\(minimum)")
                _ = try VivoBinderBundleIO.evaluateSupported(bundle: bundle, plan: plan, policy: policy, to: out, implementationSHA256: identity)
                XCTAssertEqual(try VivoBinderBundleIO.verify(out, implementationSHA256: identity).kind, "supportedEvaluation")
                let result = try VivoCanonicalJSON.decode(VivoBinderTrainingSupport.Evaluation.self,
                    from: Data(contentsOf: out.appendingPathComponent("report.json")))
                XCTAssertEqual(result.support.eligibleForExperimentalFit, minimum == 2)
                XCTAssertEqual(result.fittedReport == nil, minimum != 2)
            }
        }
    }
    func testRehashedPolicyChangeCannotRetainOldReport() throws {
        try workspace { root, bundle, plan, policy in
            let out = root.appendingPathComponent("output")
            _ = try VivoBinderBundleIO.evaluateSupported(bundle: bundle, plan: plan, policy: policy, to: out, implementationSHA256: identity)
            let changed = out.appendingPathComponent("support-policy.json")
            try writePolicy(changed, minimum: 100)
            let receiptURL = out.appendingPathComponent("receipt.json")
            var receipt = try XCTUnwrap(JSONSerialization.jsonObject(with: Data(contentsOf: receiptURL)) as? [String: Any])
            var hashes = receipt["files"] as! [String: String]
            hashes["support-policy.json"] = try VivoCanonicalJSON.fingerprint(Data(contentsOf: changed)).hex
            receipt["files"] = hashes
            try JSONSerialization.data(withJSONObject: receipt).write(to: receiptURL)
            XCTAssertThrowsError(try VivoBinderBundleIO.verify(out, implementationSHA256: identity))
        }
    }
}
