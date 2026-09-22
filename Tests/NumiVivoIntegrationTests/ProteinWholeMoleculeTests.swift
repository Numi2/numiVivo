import XCTest
@testable import NumiVivoKit

final class ProteinWholeMoleculeTests: XCTestCase {
    private func periodic(_ x: VivoVector3D) -> VivoVector3D {
        .init(x.x - x.x.rounded(), x.y - x.y.rounded(), x.z - x.z.rounded())
    }
    func testWrappedChainRetainsExtensionBeyondHalfBox() throws {
        let layout = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 1, 2, 3], parents: [.max, 0, 1, 2],
            closureEdges: [[0, 1], [1, 2], [2, 3]])
        let p: [VivoVector3D] = [.init(0.8, 0, 0), .init(0.1, 0, 0), .init(0.4, 0, 0), .init(0.7, 0, 0)]
        let whole = try layout.reconstruct(positionsNM: p, minimumImage: periodic)
        XCTAssertEqual((whole[3] - whole[0]).norm, 0.9, accuracy: 1e-12)
        XCTAssertEqual(periodic(p[3] - p[0]).norm, 0.1, accuracy: 1e-12)
        let pull = VivoProteinPullDefinition(reference: .init(particles: [0]), moving: .init(particles: [3]), stiffnessKJPerMolNM2: 10)
        let e = try VivoProteinPullMath.evaluate(pull, wholePositionsNM: whole, massesDa: [1, 1, 1, 1], referenceNM: 1)
        XCTAssertEqual(e.coordinateNM, 0.9, accuracy: 1e-12)
    }
    func testPeriodicWindingIsRejected() throws {
        let layout = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 1, 2], parents: [.max, 0, 1],
            closureEdges: [[0, 1], [1, 2], [0, 2]])
        let p: [VivoVector3D] = [.zero, .init(0.4, 0, 0), .init(0.8, 0, 0)]
        XCTAssertThrowsError(try layout.reconstruct(positionsNM: p, minimumImage: periodic))
    }
    func testMalformedTraversalIsRejected() {
        let forward = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 1, 2], parents: [.max, 2, 0], closureEdges: [])
        let duplicate = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 1, 1], parents: [.max, 0, 0], closureEdges: [])
        let edge = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 1], parents: [.max, 0], closureEdges: [[0, 2]])
        XCTAssertThrowsError(try forward.validate(particleCount: 3))
        XCTAssertThrowsError(try duplicate.validate(particleCount: 3))
        XCTAssertThrowsError(try edge.validate(particleCount: 3))
    }
    func testFiniteGeometryAndOtherParticlesRemainUnchanged() throws {
        let layout = VivoProteinWholeMoleculeLayout(orderedMembers: [0, 2], parents: [.max, 0], closureEdges: [[0, 2]])
        let p: [VivoVector3D] = [.zero, .init(100, 20, 1), .init(0.25, 0, 0)]
        XCTAssertEqual(try layout.reconstruct(positionsNM: p, minimumImage: { $0 }), p)
    }
}
