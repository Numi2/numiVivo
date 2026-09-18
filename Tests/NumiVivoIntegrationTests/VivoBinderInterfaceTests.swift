import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderInterfaceTests: XCTestCase {
    private func fixture(distance: Double = 0.4) -> VivoMolecularStructure {
        let names = ["N", "CA", "C", "O"], symbols = ["N", "C", "C", "O"]
        var atoms: [VivoMolecularAtom] = [], points: [VivoVector3D] = []
        for r in 0..<2 {
            for a in 0..<4 {
                atoms.append(.init(index: UInt32(atoms.count), name: names[a],
                    element: VivoElement.from(symbol: symbols[a])!, residueIndex: UInt32(r)))
                points.append(.init(Double(a) * 0.1, r == 0 ? 0 : distance, 0))
            }
        }
        return .init(identifier: "synthetic-two-glycines", atoms: atoms,
            residues: [.init(index: 0, name: "GLY", chainIndex: 0, atomIndices: [0,1,2,3]),
                       .init(index: 1, name: "GLY", chainIndex: 1, atomIndices: [4,5,6,7])],
            chains: [.init(index: 0, identifier: "binder", residueIndices: [0]),
                     .init(index: 1, identifier: "target", residueIndices: [1])],
            conformers: [.init(identifier: "prediction-0", positionsNM: points)])
    }
    private var plan: VivoMolecularInterface.Plan {
        .init(binderChains: ["binder"], targetChains: ["target"], conformerID: "prediction-0")
    }
    private func row(_ id: String = "a", target: String = "T", sequence: String = "G") -> VivoBinderBenchmark.Record {
        .init(id: id, target: target, leakageGroup: id, outcome: .binder, rawOutcome: "binder",
              features: ["ipsae_min_boltz2": 0.6], sourceFields: ["sequence": sequence])
    }
    private func data(_ rows: [VivoBinderBenchmark.Record]) -> VivoBinderBenchmark.Dataset {
        .init(sourceSHA256: String(repeating: "a", count: 64), assay: "synthetic", groupingMethod: "synthetic IDs", records: rows)
    }
    private func input(_ structure: VivoMolecularStructure, id: String = "a", target: String = "T") -> VivoBinderStructuralFeatures.Input {
        .init(observations: [.init(candidateID: id, target: target, sourceLabel: "synthetic geometry fixture",
                                  structure: structure, interfacePlan: plan)])
    }

    func testAnalyticalGeometry() throws {
        let result = try VivoMolecularInterface.analyze(fixture(), plan: plan)
        XCTAssertEqual(result.evaluatedAtomPairs, 16)
        XCTAssertEqual(result.contactAtomPairs, 14)
        XCTAssertEqual(result.shortDistanceAtomPairs, 0)
        XCTAssertEqual(result.residueContacts.count, 1)
        XCTAssertEqual(result.minimumDistanceNM, 0.4, accuracy: 1e-14)
        XCTAssertEqual(result.geometricCentroidDistanceNM, 0.4, accuracy: 1e-14)
        XCTAssertEqual(result.binderGeometricRadiusOfGyrationNM, sqrt(0.0125), accuracy: 1e-14)
        XCTAssertEqual(result.binderContactAtomIndices, [0,1,2,3])
        XCTAssertEqual(result.targetContactAtomIndices, [4,5,6,7])
    }
    func testRotationTranslationInvariantObservations() throws {
        let original = try VivoMolecularInterface.analyze(fixture(), plan: plan)
        var moved = fixture()
        moved.conformers[0].positionsNM = moved.conformers[0].positionsNM.map { .init(-$0.y + 10, $0.z - 3, $0.x + 7) }
        let transformed = try VivoMolecularInterface.analyze(moved, plan: plan)
        XCTAssertEqual(original.contactAtomPairs, transformed.contactAtomPairs)
        XCTAssertEqual(original.minimumDistanceNM, transformed.minimumDistanceNM, accuracy: 1e-12)
        XCTAssertEqual(original.binderGeometricRadiusOfGyrationNM, transformed.binderGeometricRadiusOfGyrationNM, accuracy: 1e-12)
        XCTAssertNotEqual(original.structureSHA256, transformed.structureSHA256)
    }
    func testSeparatedPartnersAreZeroContactsNotMissing() throws {
        let report = try VivoMolecularInterface.analyze(fixture(distance: 1), plan: plan)
        XCTAssertEqual(report.contactAtomPairs, 0)
        XCTAssertTrue(report.residueContacts.isEmpty)
        XCTAssertEqual(report.minimumDistanceNM, 1, accuracy: 1e-14)
    }
    func testCutoffSemanticsAndShortDistances() throws {
        let boundary = try VivoMolecularInterface.analyze(fixture(distance: 0.45), plan: plan)
        XCTAssertEqual(boundary.contactAtomPairs, 4)
        let shortBoundary = try VivoMolecularInterface.analyze(fixture(distance: 0.2), plan: plan)
        XCTAssertEqual(shortBoundary.shortDistanceAtomPairs, 0)
        let close = try VivoMolecularInterface.analyze(fixture(distance: 0.1), plan: plan)
        XCTAssertEqual(close.shortDistanceAtomPairs, 10)
    }
    func testRejectsPeriodicAndCovalentInterfaces() throws {
        var periodic = fixture()
        periodic.periodicCell = .init(a: .init(2,0,0), b: .init(0,2,0), c: .init(0,0,2))
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(periodic, plan: plan))
        var bonded = fixture(); bonded.bonds = [.init(atomA: 0, atomB: 4)]
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(bonded, plan: plan))
    }
    func testRejectsAmbiguousOccupancyAndIdentity() throws {
        var structure = fixture(); structure.atoms[0].alternateLocation = "A"
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(structure, plan: plan))
        structure = fixture(); structure.atoms[0].occupancy = 0.5
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(structure, plan: plan))
        structure = fixture(); structure.atoms[0].name = "CA"
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(structure, plan: plan))
    }
    func testRejectsBadSelectionsAndBudget() throws {
        for invalid in [
            VivoMolecularInterface.Plan(binderChains: ["binder"], targetChains: ["binder"], conformerID: "prediction-0"),
            .init(binderChains: ["binder"], targetChains: ["absent"], conformerID: "prediction-0"),
            .init(binderChains: ["binder"], targetChains: ["target"], conformerID: "absent"),
            .init(binderChains: ["binder"], targetChains: ["target"], conformerID: "prediction-0", maximumPairEvaluations: 15),
            .init(binderChains: ["binder"], targetChains: ["target"], conformerID: "prediction-0", shortDistanceNM: 1)
        ] { XCTAssertThrowsError(try VivoMolecularInterface.analyze(fixture(), plan: invalid)) }
    }
    func testRejectsNonfiniteGeometryAndBadSchema() throws {
        var structure = fixture(); structure.conformers[0].positionsNM[0].x = .nan
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(structure, plan: plan))
        structure = fixture(); structure.schemaVersion = 99
        XCTAssertThrowsError(try VivoMolecularInterface.analyze(structure, plan: plan))
    }
    func testGeometryJoinsByExactCandidateAndSequence() throws {
        let augmented = try VivoBinderStructuralFeatures.augment(data([row(), row("b")]), input: input(fixture()))
        XCTAssertEqual(augmented.dataset.records[0].features["numi.interface.contactAtomPairs"], 14)
        XCTAssertNil(augmented.dataset.records[1].features["numi.interface.contactAtomPairs"])
        XCTAssertEqual(augmented.unavailableCandidateIDs, ["b"])
        XCTAssertEqual(augmented.dataset.records[0].outcome, .binder)
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row(sequence: "A")]), input: input(fixture())))
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row()]), input: input(fixture(), target: "OTHER")))
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row()]), input: input(fixture(), id: "unknown")))
    }
    func testMissingHeavyAtomsAreNotRepaired() throws {
        var structure = fixture(); structure.residues[0].name = "ALA"
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row(sequence: "A")]), input: input(structure)))
        structure = fixture(); structure.atoms[1].element = VivoElement.from(symbol: "Ca")!
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row()]), input: input(structure)))
    }
    func testDuplicateObservationsAndFeatureSpoofingRejected() throws {
        let observation = input(fixture()).observations[0]
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row()]), input: .init(observations: [observation, observation])))
        let spoof = VivoBinderBenchmark.Record(id: "a", target: "T", leakageGroup: "a", outcome: .binder,
            rawOutcome: "binder", features: ["numi.interface.contactAtomPairs": 100], sourceFields: ["sequence": "G"])
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([spoof]), input: input(fixture())))
    }
    func testMatchedScoreOnlyKeepsExactPopulations() throws {
        let geometry = "numi.interface.minimumDistanceNM"
        var rows: [VivoBinderBenchmark.Record] = []
        for i in 0..<10 {
            var features = ["ipsae_min_boltz2": Double(i % 4) / 4]
            if i != 9 { features[geometry] = Double(i % 3) / 5 }
            rows.append(.init(id: "row-\(i)", target: i < 6 ? "TRAIN" : "TEST", leakageGroup: "g-\(i)",
                outcome: i % 2 == 0 ? .binder : .nonBinder, rawOutcome: i % 2 == 0 ? "binder" : "non_binder", features: features))
        }
        let dataset = data(rows), request = VivoBinderBenchmark.Plan(sourceSHA256: dataset.sourceSHA256,
            trainingTargets: ["TRAIN"], testTargets: ["TEST"], baselineFeature: "ipsae_min_boltz2",
            modelFeatures: ["ipsae_min_boltz2", geometry])
        let full = try VivoBinderBenchmark.evaluate(dataset, plan: request)
        let control = try VivoBinderStructuralFeatures.matchedScoreOnly(dataset, plan: request)
        XCTAssertEqual(full.trainingIDs, control.trainingIDs)
        XCTAssertEqual(full.targets.map(\.candidateIDs), control.targets.map(\.candidateIDs))
        XCTAssertEqual(control.model.featureNames, ["ipsae_min_boltz2"])
        XCTAssertEqual(control.targets[0].candidateIDs.count, 3)
        XCTAssertTrue(control.excluded.keys.contains("row-9"))
    }
    func testCandidateSpecificCutoffsRejected() throws {
        let a = input(fixture()).observations[0]
        let b = VivoBinderStructuralFeatures.Observation(candidateID: "b", target: "T", sourceLabel: "synthetic",
            structure: fixture(), interfacePlan: .init(binderChains: ["binder"], targetChains: ["target"],
                conformerID: "prediction-0", contactDistanceNM: 0.6))
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row(), row("b")]),
            input: .init(observations: [a, b])))
    }
    func testInconsistentTargetSequenceRejected() throws {
        let a = input(fixture()).observations[0]
        var structure = fixture()
        structure.residues[1].name = "ALA"; structure.residues[1].atomIndices.append(8)
        structure.atoms.append(.init(index: 8, name: "CB", element: VivoElement.from(symbol: "C")!, residueIndex: 1))
        structure.conformers[0].positionsNM.append(.init(0.1, 0.5, 0))
        let b = VivoBinderStructuralFeatures.Observation(candidateID: "b", target: "T", sourceLabel: "synthetic",
            structure: structure, interfacePlan: plan)
        XCTAssertThrowsError(try VivoBinderStructuralFeatures.augment(data([row(), row("b")]),
            input: .init(observations: [a, b])))
    }
    func testGeometryReportRoundTripIsDeterministic() throws {
        let first = try VivoMolecularInterface.analyze(fixture(), plan: plan)
        let bytes = try VivoCanonicalJSON.encode(first)
        XCTAssertEqual(first, try VivoCanonicalJSON.decode(VivoMolecularInterface.Report.self, from: bytes))
        XCTAssertEqual(bytes, try VivoCanonicalJSON.encode(VivoMolecularInterface.analyze(fixture(), plan: plan)))
    }
}
