import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMPathCoordinateTests {
    private func fingerprint(_ text: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }
    private func system() throws -> VivoClassicalSystem {
        let fp = try fingerprint("pathcv-structure")
        let particles = (0..<4).map { index in
            VivoClassicalParticle(index: UInt32(index), atomIndex: UInt32(index),
                typeIdentifier: "X", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
        }
        return VivoClassicalSystem(identifier: "pathcv", structureFingerprint: fp, particles: particles)
    }
    private func path() -> VivoQMMMPathCoordinate {
        .init(identifier: "two-distance-path",
            components: [
                .init(identifier: "forming", kind: .distance, atomIndices: [0,1]),
                .init(identifier: "breaking", kind: .distance, atomIndices: [2,3])
            ],
            scalesNM: [0.1,0.1],
            nodes: [
                .init(identifier: "r", valuesNM: [0.1,0.2]),
                .init(identifier: "ts", valuesNM: [0.2,0.3]),
                .init(identifier: "p", valuesNM: [0.3,0.4])
            ], lambda: 4)
    }
    private func geometry(_ x1: Double = 0.2, _ x3: Double = 0.3) throws -> VivoMDCandidateGeometry {
        try .init(particlePositionsNM: [.zero,.init(x1,0,0),.init(0,1,0),.init(x3,1,0)], periodicCell: nil)
    }

    @Test func middleNodeHasHalfProgressAndSymmetricDistanceGradient() throws {
        let resolved = try VivoQMMMResolvedPathCoordinate(source: path(), system: system())
        let result = try resolved.evaluate(geometry())
        #expect(abs(result.progress - 0.5) < 1e-12)
        #expect(result.componentValuesNM == [0.2,0.3])
        #expect(abs(result.normalizedNodeWeights.reduce(0,+) - 1) < 1e-12)
        #expect(abs(result.distanceGradients[1]?.x ?? 1) < 1e-10)
        #expect(abs(result.distanceGradients[3]?.x ?? 1) < 1e-10)
        #expect((result.progressGradients[1]?.norm ?? 0) > 0)
    }

    @Test func analyticGradientsMatchCentralDifferences() throws {
        let resolved = try VivoQMMMResolvedPathCoordinate(source: path(), system: system())
        let center = try geometry(0.215,0.285)
        let value = try resolved.evaluate(center)
        let h = 1e-6
        let plusX = try resolved.evaluate(geometry(0.215+h,0.285))
        let minusX = try resolved.evaluate(geometry(0.215-h,0.285))
        let dsDx = (plusX.progress-minusX.progress)/(2*h)
        let dzDx = (plusX.distanceFromPath-minusX.distanceFromPath)/(2*h)
        #expect(abs(dsDx - (value.progressGradients[1]?.x ?? .nan)) < 1e-5)
        #expect(abs(dzDx - (value.distanceGradients[1]?.x ?? .nan)) < 1e-5)

        let plusY = try resolved.evaluate(geometry(0.215,0.285+h))
        let minusY = try resolved.evaluate(geometry(0.215,0.285-h))
        let dsDy = (plusY.progress-minusY.progress)/(2*h)
        let dzDy = (plusY.distanceFromPath-minusY.distanceFromPath)/(2*h)
        #expect(abs(dsDy - (value.progressGradients[3]?.x ?? .nan)) < 1e-5)
        #expect(abs(dzDy - (value.distanceGradients[3]?.x ?? .nan)) < 1e-5)
    }

    @Test func malformedPathIsRejected() throws {
        var duplicate = path()
        duplicate.nodes[1].valuesNM = duplicate.nodes[0].valuesNM
        #expect(throws: (any Error).self) { try duplicate.validate() }
        var badScale = path(); badScale.scalesNM[0] = 0
        #expect(throws: (any Error).self) { try badScale.validate() }
    }
}
