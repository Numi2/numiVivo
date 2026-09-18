import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderBundleTests: XCTestCase {
    let identity = String(repeating: "c", count: 64)
    func workspace(_ body: (URL, URL, URL) throws -> Void) throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        try FileManager.default.createDirectory(at: root, withIntermediateDirectories: true)
        defer { try? FileManager.default.removeItem(at: root) }
        let source = root.appendingPathComponent("input.csv"), config = root.appendingPathComponent("config.json")
        var csv = "uuid,target,sequence,adaptyv_binding,ipsae_min_boltz2\n"
        let sequences = ["ACDA", "ACDC", "ACDD", "ACDE", "ACDF", "ACDG", "ACDH", "ACDI", "ACDK", "ACDL", "ACDM", "ACDN"]
        for i in 0..<12 {
            csv += "c\(i),\(i < 8 ? "BBF-14" : "MBP"),\(sequences[i]),\(i % 2 == 0 ? "non_binder" : "binder"),\(i % 2 == 0 ? "0.1" : "0.9")\n"
        }
        try Data(csv.utf8).write(to: source)
        try VivoCanonicalJSON.encode(VivoBinderAnthropicImport.Configuration(assay: .adaptyv,
            targets: ["BBF-14", "MBP"], features: ["ipsae_min_boltz2"])).write(to: config)
        try body(root, source, config)
    }
    func testImportReplayAndHash() throws {
        try workspace { root, source, config in
            let out = root.appendingPathComponent("bundle")
            let r = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: out, implementationSHA256: identity)
            XCTAssertEqual(r.files["source.csv"], try VivoCanonicalJSON.fingerprint(Data(contentsOf: source)).hex)
            XCTAssertEqual(try VivoBinderBundleIO.verify(out, implementationSHA256: identity).kind, "import")
        }
    }
    func testEvaluationReplay() throws {
        try workspace { root, source, config in
            let input = root.appendingPathComponent("bundle"), output = root.appendingPathComponent("result")
            _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: input, implementationSHA256: identity)
            let plan = VivoBinderBenchmark.Plan(sourceSHA256: try VivoCanonicalJSON.fingerprint(Data(contentsOf: source)).hex,
                trainingTargets: ["BBF-14"], testTargets: ["MBP"], baselineFeature: "ipsae_min_boltz2",
                modelFeatures: ["ipsae_min_boltz2"], topK: 2)
            let file = root.appendingPathComponent("plan.json"); try VivoCanonicalJSON.encode(plan).write(to: file)
            _ = try VivoBinderBundleIO.evaluate(bundle: input, plan: file, to: output, implementationSHA256: identity)
            XCTAssertEqual(try VivoBinderBundleIO.verify(output, implementationSHA256: identity).kind, "evaluation")
        }
    }
    func testNoOverwriteOrPartialOutput() throws {
        try workspace { root, source, config in
            let output = root.appendingPathComponent("bundle")
            _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: output, implementationSHA256: identity)
            XCTAssertThrowsError(try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: output, implementationSHA256: identity))
            try Data("invalid".utf8).write(to: source)
            let failed = root.appendingPathComponent("failed")
            XCTAssertThrowsError(try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: failed, implementationSHA256: identity))
            XCTAssertFalse(FileManager.default.fileExists(atPath: failed.path))
            XCTAssertFalse(try FileManager.default.contentsOfDirectory(atPath: root.path).contains { $0.hasPrefix(".numivivo-binder-") })
        }
    }
    func testTamperingAndWrongExecutableRejected() throws {
        try workspace { root, source, config in
            let output = root.appendingPathComponent("bundle")
            _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: output, implementationSHA256: identity)
            XCTAssertThrowsError(try VivoBinderBundleIO.verify(output, implementationSHA256: String(repeating: "d", count: 64)))
            try Data("changed".utf8).write(to: output.appendingPathComponent("source.csv"))
            XCTAssertThrowsError(try VivoBinderBundleIO.verify(output, implementationSHA256: identity))
        }
    }
    func testChangedImportedRecordEvenWithUpdatedHashRejected() throws {
        try workspace { root, source, config in
            let output = root.appendingPathComponent("bundle")
            _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: output, implementationSHA256: identity)
            let artifact = output.appendingPathComponent("imported.json")
            let changed = Data(String(decoding: try Data(contentsOf: artifact), as: UTF8.self)
                .replacingOccurrences(of: "\"rawOutcome\":\"binder\"", with: "\"rawOutcome\":\"invented\"").utf8)
            try changed.write(to: artifact)
            let receiptURL = output.appendingPathComponent("receipt.json")
            var receipt = try XCTUnwrap(JSONSerialization.jsonObject(with: Data(contentsOf: receiptURL)) as? [String: Any])
            var hashes = try XCTUnwrap(receipt["files"] as? [String: String])
            hashes["imported.json"] = try VivoCanonicalJSON.fingerprint(changed).hex; receipt["files"] = hashes
            try JSONSerialization.data(withJSONObject: receipt).write(to: receiptURL)
            XCTAssertThrowsError(try VivoBinderBundleIO.verify(output, implementationSHA256: identity))
        }
    }
    func testUnexpectedFilesRejected() throws {
        try workspace { root, source, config in
            let output = root.appendingPathComponent("bundle")
            _ = try VivoBinderBundleIO.importCSV(source: source, configuration: config, to: output, implementationSHA256: identity)
            try Data().write(to: output.appendingPathComponent("unexpected"))
            XCTAssertThrowsError(try VivoBinderBundleIO.verify(output, implementationSHA256: identity))
        }
    }
}
