import Foundation
import XCTest
@testable import NumiVivoKit

final class ProteinStressObservableTests: XCTestCase {
    private let masses = [12.0, 16, 12, 1]
    private let coordinates: [VivoVector3D] = [.init(0, 0, 0), .init(0.1, 0.2, 0), .init(0.8, 0.3, 0.2), .init(1, 0.1, 0.4)]
    private func definition(axis: VivoVector3D? = nil) -> VivoProteinPullDefinition {
        .init(reference: .init(particles: [0, 1]), moving: .init(particles: [2, 3]),
              projectionAxis: axis, stiffnessKJPerMolNM2: 450)
    }
    private func evaluation(_ positions: [VivoVector3D], axis: VivoVector3D? = nil, reference: Double = 1.2) throws -> VivoProteinPullEvaluation {
        try VivoProteinPullMath.evaluate(definition(axis: axis), wholePositionsNM: positions, massesDa: masses, referenceNM: reference)
    }
    func testRadialForceMatchesIndependentEnergyDifference() throws {
        let e = try evaluation(coordinates), h = 1e-6
        for i in coordinates.indices {
            for axis in 0..<3 {
                var plus = coordinates, minus = coordinates
                switch axis {
                case 0: plus[i].x += h; minus[i].x -= h
                case 1: plus[i].y += h; minus[i].y -= h
                default: plus[i].z += h; minus[i].z -= h
                }
                let derivative = try (evaluation(plus).energyKJPerMol - evaluation(minus).energyKJPerMol) / (2 * h)
                let actual = [e.forcesKJPerMolNM[i].x, e.forcesKJPerMolNM[i].y, e.forcesKJPerMolNM[i].z][axis]
                XCTAssertEqual(actual, -derivative, accuracy: 1e-7)
            }
        }
        XCTAssertLessThan(e.forcesKJPerMolNM.reduce(.zero, +).norm, 1e-12)
        XCTAssertEqual(e.forcesKJPerMolNM[0].norm / e.forcesKJPerMolNM[1].norm, 12 / 16, accuracy: 1e-12)
    }
    func testProjectionGradientAndAxisNormalization() throws {
        let axis = VivoVector3D(2, -1, 3), e = try evaluation(coordinates, axis: axis)
        XCTAssertEqual(e.energyKJPerMol, try evaluation(coordinates, axis: axis * 10).energyKJPerMol, accuracy: 1e-12)
        let h = 1e-6
        var p = coordinates, m = coordinates; p[3].y += h; m[3].y -= h
        let d = try (evaluation(p, axis: axis).energyKJPerMol - evaluation(m, axis: axis).energyKJPerMol) / (2 * h)
        XCTAssertEqual(e.forcesKJPerMolNM[3].y, -d, accuracy: 1e-7)
        XCTAssertGreaterThan(try evaluation(coordinates, axis: axis, reference: -0.5).energyKJPerMol, 0)
    }
    func testUnitConversionAndSimpleSpring() throws {
        let pull = VivoProteinPullDefinition(reference: .init(particles: [0]), moving: .init(particles: [1]), stiffnessKJPerMolNM2: 100)
        let e = try VivoProteinPullMath.evaluate(pull, wholePositionsNM: [.zero, .init(1, 0, 0)], massesDa: [1, 1], referenceNM: 1.1)
        XCTAssertEqual(e.energyKJPerMol, 0.5, accuracy: 1e-12)
        XCTAssertEqual(e.tensileForcePN, 16.605390671738467, accuracy: 1e-11)
        XCTAssertEqual(e.forcesKJPerMolNM[1].x, 10, accuracy: 1e-12)
    }
    func testForceAndEnergyAreTranslationInvariant() throws {
        let translated = coordinates.map { $0 + .init(2, -3, 1) }
        let a = try evaluation(coordinates), b = try evaluation(translated)
        XCTAssertEqual(a.coordinateNM, b.coordinateNM, accuracy: 1e-14)
        XCTAssertEqual(a.energyKJPerMol, b.energyKJPerMol, accuracy: 1e-12)
    }
    func testSwitchingWorkHasCorrectSignAndTelescopes() throws {
        let q = 0.7, k = 450.0
        let w = try VivoProteinPullMath.switchingWork(coordinateNM: q, from: 0.8, to: 1.2, stiffnessKJPerMolNM2: k)
        XCTAssertEqual(w, 0.5 * k * (pow(q - 1.2, 2) - pow(q - 0.8, 2)), accuracy: 1e-12)
        let a = try VivoProteinPullMath.switchingWork(coordinateNM: q, from: 0.8, to: 1, stiffnessKJPerMolNM2: k)
        let b = try VivoProteinPullMath.switchingWork(coordinateNM: q, from: 1, to: 1.2, stiffnessKJPerMolNM2: k)
        XCTAssertEqual(a + b, w, accuracy: 1e-12)
        XCTAssertEqual(try VivoProteinPullMath.switchingWork(coordinateNM: q, from: 1.2, to: 0.8, stiffnessKJPerMolNM2: k), -w, accuracy: 1e-12)
    }
    func testMalformedPullsAndSingularGeometryAreRejected() throws {
        let duplicate = VivoProteinPullDefinition(reference: .init(particles: [0, 0]), moving: .init(particles: [2]), stiffnessKJPerMolNM2: 1)
        let overlap = VivoProteinPullDefinition(reference: .init(particles: [0]), moving: .init(particles: [0]), stiffnessKJPerMolNM2: 1)
        XCTAssertThrowsError(try duplicate.validate(massesDa: masses))
        XCTAssertThrowsError(try overlap.validate(massesDa: masses))
        XCTAssertThrowsError(try definition(axis: .zero).validate(massesDa: masses))
        XCTAssertThrowsError(try definition().validate(massesDa: [0, 0, 1, 1]))
        XCTAssertThrowsError(try evaluation([.zero, .zero, .zero, .zero]))
        XCTAssertThrowsError(try evaluation(coordinates, reference: -.infinity))
        XCTAssertThrowsError(try evaluation(coordinates, reference: -1))
    }
    func testHydrogenBondUsesDHAAngleNotCovalentBondCount() throws {
        let bond = VivoProteinHydrogenBond(donor: 0, hydrogen: 1, acceptor: 2)
        let straight = try VivoProteinHydrogenBondAnalysis.evaluate([bond], positionsNM: [.zero, .init(0.1, 0, 0), .init(0.28, 0, 0)])[0]
        XCTAssertTrue(straight.present); XCTAssertEqual(straight.dhaAngleDegrees, 180, accuracy: 1e-12)
        let bent = try VivoProteinHydrogenBondAnalysis.evaluate([bond], positionsNM: [.zero, .init(0.1, 0, 0), .init(0.1, 0.18, 0)])[0]
        XCTAssertFalse(bent.present); XCTAssertEqual(bent.dhaAngleDegrees, 90, accuracy: 1e-12)
        let far = try VivoProteinHydrogenBondAnalysis.evaluate([bond], positionsNM: [.zero, .init(0.1, 0, 0), .init(0.5, 0, 0)])[0]
        XCTAssertFalse(far.present)
        XCTAssertThrowsError(try VivoProteinHydrogenBondAnalysis.evaluate([bond, bond], positionsNM: [.zero, .init(0.1, 0, 0), .init(0.28, 0, 0)]))
        XCTAssertThrowsError(try VivoProteinHydrogenBondAnalysis.evaluate([bond], positionsNM: [.zero, .zero, .init(0.28, 0, 0)]))
    }
    func testIrregularTimeOccupancyAndReformation() throws {
        let p = try VivoProteinHydrogenBondAnalysis.persistence(timesPS: [0, 1, 4], present: [[true, false], [false, true], [true, true]])
        XCTAssertEqual(p[0].timeWeightedOccupancy, 0.25, accuracy: 1e-14)
        XCTAssertEqual(p[0].firstObservedLossTimePS, 1); XCTAssertEqual(p[0].observedReformations, 1)
        XCTAssertEqual(p[1].timeWeightedOccupancy, 0.75, accuracy: 1e-14)
        XCTAssertFalse(p[1].initiallyPresent); XCTAssertNil(p[1].firstObservedLossTimePS)
        XCTAssertThrowsError(try VivoProteinHydrogenBondAnalysis.persistence(timesPS: [1, 1], present: [[true], [false]]))
        XCTAssertThrowsError(try VivoProteinHydrogenBondAnalysis.persistence(timesPS: [0, 1], present: [[true], []]))
    }
    func testStructuralRetentionIsInvariantAndMissingContactsStayMissing() throws {
        let x: [VivoVector3D] = [.zero, .init(1, 0, 0)]
        let contact = VivoProteinNativeContact(a: 0, b: 1, referenceDistanceNM: 1)
        let a = try VivoProteinStructuralAnalysis.evaluate(positionsNM: x, massesDa: [1, 1], selection: [0, 1], contacts: [contact])
        XCTAssertEqual(a.radiusOfGyrationNM, 0.5, accuracy: 1e-14)
        XCTAssertEqual(a.retainedContactFraction, 1); XCTAssertEqual(a.contactDistanceRMSErrorNM, 0)
        let rotated: [VivoVector3D] = [.init(2, 3, 4), .init(2, 4, 4)]
        XCTAssertEqual(a, try VivoProteinStructuralAnalysis.evaluate(positionsNM: rotated, massesDa: [1, 1], selection: [0, 1], contacts: [contact]))
        let missing = try VivoProteinStructuralAnalysis.evaluate(positionsNM: x, massesDa: [1, 1], selection: [0, 1], contacts: [])
        XCTAssertNil(missing.retainedContactFraction); XCTAssertNil(missing.contactDistanceRMSErrorNM)
        let unfolded = try VivoProteinStructuralAnalysis.evaluate(positionsNM: [.zero, .init(2, 0, 0)], massesDa: [1, 1], selection: [0, 1], contacts: [contact])
        XCTAssertEqual(unfolded.retainedContactFraction, 0)
        XCTAssertThrowsError(try VivoProteinStructuralAnalysis.evaluate(positionsNM: x, massesDa: [1, 1], selection: [0, 1], contacts: [contact, .init(a: 1, b: 0, referenceDistanceNM: 1)]))
    }
    func testReplicaStatisticsDoNotFabricateSingleRunUncertainty() throws {
        let single = try VivoProteinReplicaStatistics.calculate([4])
        XCTAssertNil(single.sampleStandardDeviation); XCTAssertNil(single.standardErrorOfMean)
        let multiple = try VivoProteinReplicaStatistics.calculate([1, 2, 3])
        XCTAssertEqual(multiple.mean, 2); XCTAssertEqual(multiple.sampleStandardDeviation, 1)
        XCTAssertEqual(multiple.standardErrorOfMean!, 1 / sqrt(3), accuracy: 1e-14)
        XCTAssertThrowsError(try VivoProteinReplicaStatistics.calculate([]))
        XCTAssertThrowsError(try VivoProteinReplicaStatistics.calculate([.nan]))
    }
    func testPublicModelsRoundTrip() throws {
        let d = definition(axis: .init(1, 2, 3)), encoder = JSONEncoder()
        XCTAssertEqual(d, try JSONDecoder().decode(VivoProteinPullDefinition.self, from: encoder.encode(d)))
    }
}
