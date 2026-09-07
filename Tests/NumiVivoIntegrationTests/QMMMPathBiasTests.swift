import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMPathBiasTests {
    private func fingerprint(_ text: String) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(Data(text.utf8))
    }
    private func system() throws -> VivoClassicalSystem {
        let fp = try fingerprint("path-bias-structure")
        let particles = (0..<4).map { index in
            VivoClassicalParticle(index: UInt32(index), atomIndex: UInt32(index),
                typeIdentifier: "X", massDa: 12, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
        }
        return VivoClassicalSystem(identifier: "path-bias", structureFingerprint: fp, particles: particles)
    }
    private func path() -> VivoQMMMPathCoordinate {
        .init(identifier: "two-distance-path",
            components: [
                .init(identifier: "a", kind: .distance, atomIndices: [0,1]),
                .init(identifier: "b", kind: .distance, atomIndices: [2,3])
            ], scalesNM: [0.1,0.1], nodes: [
                .init(identifier: "r", valuesNM: [0.1,0.2]),
                .init(identifier: "m", valuesNM: [0.2,0.3]),
                .init(identifier: "p", valuesNM: [0.3,0.4])
            ], lambda: 3)
    }
    private func geometry(_ x1: Double, _ x3: Double) throws -> VivoMDCandidateGeometry {
        try .init(particlePositionsNM: [.zero,.init(x1,0,0),.init(0,1,0),.init(x3,1,0)], periodicCell: nil)
    }
    private func base(_ system: VivoClassicalSystem) throws -> VivoMDCandidateForceProvider {
        let id = try fingerprint("zero-path-bias-provider"), systemID = try system.fingerprint()
        return try .init(fingerprint: id, retainedSystemFingerprint: systemID,
            boundary: .finiteCluster, supportsCellMoves: false,
            maximumAcceptedResidual: 1e-8, molecularConnectivitySystem: system) { geometry in
            try .init(providerFingerprint: id, geometry: geometry, additionalEnergyKJPerMol: 0,
                physicalParticleForcesKJPerMolNM: Array(repeating: .zero, count: system.particles.count),
                derivativeMethod: "synthetic zero provider", convergenceResidual: 0, requiredResidual: 1e-8)
        }
    }

    @Test func analyticBiasForceMatchesEnergyDifference() async throws {
        let system = try system()
        let provider = try VivoQMMMPathBias.provider(base: base(system), system: system, path: path(),
            bias: .init(progressCenter: 0.35, progressForceConstantKJPerMol: 12,
                        distanceCenter: 0, distanceForceConstantKJPerMol: 7))
        let center = try geometry(0.218,0.287)
        let evaluation = try await provider.evaluate(center)
        let h = 1e-6
        let plus = try await provider.evaluate(geometry(0.218+h,0.287))
        let minus = try await provider.evaluate(geometry(0.218-h,0.287))
        let derivative = (plus.additionalEnergyKJPerMol-minus.additionalEnergyKJPerMol)/(2*h)
        #expect(abs(evaluation.physicalParticleForcesKJPerMolNM[1].x + derivative) < 2e-4)
        #expect(evaluation.additionalEnergyKJPerMol > 0)
        #expect(evaluation.additionalAffineStrainDerivativeKJPerMol == nil)
    }

    @Test func disabledBothTermsIsRejected() throws {
        let invalid = VivoQMMMPathBiasConfiguration(progressCenter: 0.5,
            progressForceConstantKJPerMol: 0, distanceCenter: 0, distanceForceConstantKJPerMol: 0)
        #expect(throws: (any Error).self) { try invalid.validate() }
    }
}
