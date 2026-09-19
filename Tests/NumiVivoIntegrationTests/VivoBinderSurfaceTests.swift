import Foundation
import XCTest
@testable import NumiVivoKit

final class VivoBinderSurfaceTests: XCTestCase {
    let ip = VivoMolecularInterface.Plan(binderChains: ["A"], targetChains: ["B"], conformerID: "m")
    func plan(_ points: Int = 960, tests: Int = 1_000_000) -> VivoMolecularInterface.SurfacePlan {
        .init(radiusProfile: "explicit-test-carbon", radiiNM: ["C": 0.17], pointsPerAtom: points,
              maximumPointTests: tests)
    }
    func spheres(_ d: Double, offset: VivoVector3D = .zero) -> VivoMolecularStructure {
        .init(identifier: "two-spheres", atoms: [
            .init(index: 0, name: "CA", element: .from(symbol: "C")!, residueIndex: 0),
            .init(index: 1, name: "CA", element: .from(symbol: "C")!, residueIndex: 1)],
            residues: [.init(index: 0, name: "GLY", chainIndex: 0, atomIndices: [0]),
                       .init(index: 1, name: "GLY", chainIndex: 1, atomIndices: [1])],
            chains: [.init(index: 0, identifier: "A", residueIndices: [0]),
                     .init(index: 1, identifier: "B", residueIndices: [1])],
            conformers: [.init(identifier: "m", positionsNM: [offset, offset + .init(0, 0, d)])])
    }
    func testSeparatedSphereAnalyticalAreaAndUnits() throws {
        let r = try VivoMolecularInterface.analyzeSurface(spheres(1), interfacePlan: ip, plan: plan())
        XCTAssertEqual(r.binderIsolatedAreaNM2, 4 * .pi * 0.31 * 0.31, accuracy: 1e-14)
        XCTAssertEqual(r.buriedAreaSumNM2, 0)
        XCTAssertEqual(r.neighborPairs, 0)
        XCTAssertEqual(r.pointTests, 0)
    }
    func testOverlappingSpheresAnalyticalCap() throws {
        let r = try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip, plan: plan(4096))
        let cap = 2 * Double.pi * 0.31 * (0.31 - 0.15)
        XCTAssertEqual(r.binderBuriedAreaNM2, cap, accuracy: 0.001)
        XCTAssertEqual(r.targetBuriedAreaNM2, cap, accuracy: 0.001)
        XCTAssertEqual(r.halfBuriedAreaSumNM2, cap, accuracy: 0.001)
        XCTAssertEqual(r.binderIsolatedAreaNM2 + r.targetIsolatedAreaNM2 - r.complexAreaNM2,
                       r.buriedAreaSumNM2, accuracy: 1e-14)
        XCTAssertTrue(r.atoms.allSatisfy { $0.buriedAreaNM2 >= 0 })
    }
    func testOverlapIsVDWNotProbeExpanded() throws {
        let r = try VivoMolecularInterface.analyzeSurface(spheres(0.2), interfacePlan: ip, plan: plan())
        XCTAssertEqual(r.maximumCrossPartnerOverlapNM, 0.14, accuracy: 1e-14)
        XCTAssertEqual(r.crossPartnerOverlapPairs, 1)
    }
    func testSharedAdmissionAndMissingRadii() throws {
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(spheres(0), interfacePlan: ip, plan: plan()))
        var s = spheres(0.3); s.atoms[0].occupancy = 0.5
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(s, interfacePlan: ip, plan: plan()))
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip,
            plan: .init(radiusProfile: "missing", radiiNM: ["O": 0.15])))
    }
    func testResolutionAndBudgetsReject() throws {
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip, plan: plan(10)))
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip, plan: plan(tests: 1)))
        XCTAssertThrowsError(try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip,
            plan: .init(radiusProfile: "bad", radiiNM: ["C": -.infinity])))
    }
    func testTranslationAndReplay() throws {
        let a = try VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip, plan: plan())
        let b = try VivoMolecularInterface.analyzeSurface(spheres(0.3, offset: .init(10, -20, 30)), interfacePlan: ip, plan: plan())
        XCTAssertEqual(a.buriedAreaSumNM2, b.buriedAreaSumNM2, accuracy: 1e-14)
        let bytes = try VivoCanonicalJSON.encode(a)
        XCTAssertEqual(try VivoCanonicalJSON.decode(VivoMolecularInterface.SurfaceReport.self, from: bytes), a)
    }
    func testPolarLabelIsNotInferredFromElements() throws {
        let data = try VivoCanonicalJSON.encode(VivoMolecularInterface.analyzeSurface(spheres(0.3), interfacePlan: ip, plan: plan()))
        XCTAssertFalse(String(decoding: data, as: UTF8.self).contains("polarArea"))
    }
}
