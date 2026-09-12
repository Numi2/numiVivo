import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoArc2026Tests: XCTestCase {
    private func query(phase: String = "validation") -> VivoArc2026Query {
        let genes = (0..<18_533).map { "G\($0)" }, targets = (0..<300).map { "G\($0)" }
        let inspection = VivoArc2026CountInspection(cells: 1, entries: 1, totalCounts: 1,
            maximumCellTotal: 1, cellsPerLabel: ["non-targeting": 1])
        return .init(schemaVersion: 1, phase: phase, evaluatorCommit: VivoArc2026.evaluatorCommit,
            contexts: ["A", "B", "C"].map {
                .init(id: $0, perturbationColumn: "perturbation", genes: genes, targets: targets,
                      controlsSHA256: String(repeating: "a", count: 64), targetsSHA256: String(repeating: "b", count: 64), inspection: inspection)
            })
    }
    func testExact2026IdentityAndCompletePhase() throws {
        XCTAssertEqual(VivoArc2026.evaluatorVersion, "0.16.0")
        XCTAssertEqual(VivoArc2026.ruleVersion, 3)
        XCTAssertEqual(VivoArc2026.metrics.count, 6)
        XCTAssertEqual(Set(VivoArc2026.metrics).count, 6)
        try query().validate(); try query(phase: "final-test").validate()
        XCTAssertThrowsError(try query(phase: "training").validate())
        let q = query()
        XCTAssertThrowsError(try VivoArc2026Query(schemaVersion: 1, phase: q.phase,
            evaluatorCommit: q.evaluatorCommit, contexts: Array(q.contexts.prefix(2))).validate())
    }
    func testCountBudgetAndUnconstrainedPredictionCellCounts() throws {
        // Three cells, not 400: 400 is the reference preparation, not an output requirement.
        var a = try VivoArc2026CountAccumulator(labels: ["non-targeting", "G0", "G0"], featureCount: 3,
            expectedLabels: ["non-targeting", "G0"])
        try a.append(row: 2, feature: 1, count: 7) // CSC-like non-row traversal is valid.
        try a.append(row: 0, feature: 0, count: 1_000_000)
        XCTAssertThrowsError(try a.append(row: 0, feature: 1, count: 1))
        XCTAssertThrowsError(try a.append(row: 0, feature: 2, count: UInt64.max))
        let r = try a.finish()
        XCTAssertEqual(r.cells, 3); XCTAssertEqual(r.entries, 2)
        XCTAssertEqual(r.totalCounts, 1_000_007); XCTAssertEqual(r.cellsPerLabel["G0"], 2)
    }
    func testMissingAndExtraLabelsAreRejected() throws {
        XCTAssertThrowsError(try VivoArc2026CountAccumulator(labels: ["G0"], featureCount: 2,
            expectedLabels: ["non-targeting", "G0"]))
        XCTAssertThrowsError(try VivoArc2026CountAccumulator(labels: ["non-targeting", "G1"], featureCount: 2,
            expectedLabels: ["non-targeting", "G0"]))
        XCTAssertThrowsError(try VivoArc2026CountAccumulator(labels: [], featureCount: 2,
            expectedLabels: ["non-targeting"]))
    }
    func testUnknownRequiredSemanticsAreNotIgnored() throws {
        let json = Data(#"{"schemaVersion":1,"phase":"validation","contexts":[],"referenceOutcomes":"forbidden"}"#.utf8)
        XCTAssertThrowsError(try JSONDecoder().decode(VivoArc2026PreparePlan.self, from: json))
        XCTAssertFalse(VivoArc2026.identifier("../A"))
        XCTAssertFalse(VivoArc2026.sha256(String(repeating: "z", count: 64)))
    }
    func testCommonOrderedGeneAxisIsRequired() throws {
        let q = query(), c = q.contexts[0]
        var genes = c.genes; genes.swapAt(0, 1)
        let altered = VivoArc2026QueryContext(id: c.id, perturbationColumn: c.perturbationColumn,
            genes: genes, targets: c.targets, controlsSHA256: c.controlsSHA256,
            targetsSHA256: c.targetsSHA256, inspection: c.inspection)
        XCTAssertThrowsError(try VivoArc2026Query(schemaVersion: 1, phase: q.phase,
            evaluatorCommit: q.evaluatorCommit, contexts: [altered] + Array(q.contexts.dropFirst())).validate())
    }
}
