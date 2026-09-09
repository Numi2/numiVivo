import Foundation
import XCTest
@testable import NumiVivoKit

@MainActor
final class VivoNeoantigenIntegrationTests: XCTestCase {
    func testMetadataAndMissingEvidence() throws {
        let example = try VivoNeoantigenExample.inputs()
        let implementation = try VivoCanonicalJSON.fingerprint(Data("test-implementation".utf8))
        let result = try VivoNeoantigenWorkbench.analyze(manifest: example.manifest, tsv: example.tsv, implementationSHA256: implementation.hex)
        XCTAssertEqual(result.state, .researchReview)
        XCTAssertEqual(result.candidates.count, 3)
        XCTAssertNil(result.candidates[1].geneExpression)
        XCTAssertEqual(result.candidates[2].geneExpression, 0)
        XCTAssertEqual(result.candidates[0].sourceFields["Evaluation"], "Accept")
        XCTAssertTrue(result.candidates[0].evidenceGaps.contains { $0.contains("not been experimentally") })
        var mismatch = example.manifest
        mismatch.samples[0].assembly = "GRCh37"
        let blocked = try VivoNeoantigenWorkbench.analyze(manifest: mismatch, tsv: example.tsv, implementationSHA256: implementation.hex)
        XCTAssertEqual(blocked.state, .blocked)
        XCTAssertTrue(blocked.candidates.isEmpty)
    }

    func testArtifactReconstructionAndCorrectlyHashedForgery() async throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root)
        let implementation = try VivoCanonicalJSON.fingerprint(Data("test-implementation".utf8))
        let example = try VivoNeoantigenExample.inputs()
        let result = try await VivoNeoantigenArtifacts.publish(manifestData: VivoCanonicalJSON.encode(example.manifest),
            tsv: example.tsv, implementation: implementation, store: store)
        let verified = try await VivoNeoantigenArtifacts.verify(result.receipt, implementation: implementation, store: store)
        XCTAssertEqual(verified, result.report)
        var forged = result.report
        forged.candidates[1].geneExpression = 99
        let saved = try await store.put(data: VivoCanonicalJSON.encode(forged), kind: "neoantigen.report", mediaType: "application/json")
        var receipt = result.receipt; receipt.report = saved.fingerprint
        do {
            _ = try await VivoNeoantigenArtifacts.verify(receipt, implementation: implementation, store: store)
            XCTFail("A correctly hashed but reconstructed-result mismatch must be rejected.")
        } catch let error as VivoNeoantigenError {
            XCTAssertTrue(error.description.contains("does not reconstruct"))
        }
    }

    func testBoundResearchReviewAndNoOverwrite() async throws {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent(UUID().uuidString)
        defer { try? FileManager.default.removeItem(at: root) }
        let store = try VivoArtifactStore(rootURL: root)
        let implementation = try VivoCanonicalJSON.fingerprint(Data("test-implementation".utf8))
        let example = try VivoNeoantigenExample.inputs()
        let result = try await VivoNeoantigenArtifacts.publish(manifestData: VivoCanonicalJSON.encode(example.manifest),
            tsv: example.tsv, implementation: implementation, store: store)
        var review = VivoNeoantigenReview(schema: "numivivo.org/neoantigen-review/v1", reportSHA256: result.receipt.report.hex,
            reviewerID: "research-reviewer", reviewedAt: Date(timeIntervalSince1970: 1_700_000_000),
            decisions: [.init(candidateID: result.report.candidates[0].id, disposition: .deferReview,
                              rationale: "Await independent experimental recognition evidence.")])
        let first = try await VivoNeoantigenArtifacts.record(review, receipt: result.receipt, implementation: implementation, store: store)
        review.decisions[0].rationale = "Additional evidence is still required."
        let second = try await VivoNeoantigenArtifacts.record(review, receipt: result.receipt, implementation: implementation, store: store)
        XCTAssertNotEqual(first.fingerprint, second.fingerprint)
        let original = try await store.data(for: first.fingerprint)
        XCTAssertTrue(String(decoding: original, as: UTF8.self).contains("independent experimental"))
        review.reportSHA256 = String(repeating: "0", count: 64)
        do {
            _ = try await VivoNeoantigenArtifacts.record(review, receipt: result.receipt, implementation: implementation, store: store)
            XCTFail("Review transplantation must be rejected.")
        } catch is VivoNeoantigenError { }
    }
}
