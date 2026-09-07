import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct QMMMReactionCoordinateTests {
    private func system() throws -> VivoClassicalSystem {
        let structure=try VivoCanonicalJSON.fingerprint(Data("shared-transfer-structure".utf8))
        return VivoClassicalSystem(identifier:"shared-transfer",structureFingerprint:structure,particles:[
            .init(index:0,atomIndex:0,typeIdentifier:"D",massDa:14,chargeE:0,sigmaNM:0,epsilonKJPerMol:0),
            .init(index:1,atomIndex:1,typeIdentifier:"H",massDa:1,chargeE:0,sigmaNM:0,epsilonKJPerMol:0),
            .init(index:2,atomIndex:2,typeIdentifier:"A",massDa:16,chargeE:0,sigmaNM:0,epsilonKJPerMol:0)
        ])
    }

    @Test func sharedHydrogenDistanceDifferenceAccumulatesJacobian() throws {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"proton-transfer",kind:.distanceDifference,
            atomIndices:[0,1,1,2])
        try coordinate.validate()
        let resolved=try VivoQMMMResolvedCoordinate(source:coordinate,system:system())
        let geometry=try VivoMDCandidateGeometry(particlePositionsNM:[.init(0,0,0),.init(0.10,0,0),.init(0.25,0,0)],periodicCell:nil)
        let value=try resolved.evaluate(geometry)
        #expect(abs(value.valueNM + 0.05)<1e-12)
        #expect(value.gradients.count==3)
        #expect(abs(value.gradients[0]!.x + 1)<1e-12)
        #expect(abs(value.gradients[1]!.x - 2)<1e-12)
        #expect(abs(value.gradients[2]!.x + 1)<1e-12)
    }

    @Test func sharedCoordinateJacobianMatchesFiniteDifference() throws {
        let coordinate=VivoQMMMReactionCoordinate(identifier:"bent-proton-transfer",kind:.distanceDifference,
            atomIndices:[0,1,1,2])
        let resolved=try VivoQMMMResolvedCoordinate(source:coordinate,system:system())
        var positions=[VivoVector3D(0,0,0),VivoVector3D(0.10,0.03,0),VivoVector3D(0.22,-0.02,0)]
        let base=try resolved.evaluate(VivoMDCandidateGeometry(particlePositionsNM:positions,periodicCell:nil))
        let h=1e-7
        for particle in 0..<3 {
            positions[particle].x += h
            let plus=try resolved.evaluate(VivoMDCandidateGeometry(particlePositionsNM:positions,periodicCell:nil)).valueNM
            positions[particle].x -= 2*h
            let minus=try resolved.evaluate(VivoMDCandidateGeometry(particlePositionsNM:positions,periodicCell:nil)).valueNM
            positions[particle].x += h
            let numeric=(plus-minus)/(2*h)
            #expect(abs(numeric-base.gradients[UInt32(particle)]!.x)<1e-7)
        }
    }
}
