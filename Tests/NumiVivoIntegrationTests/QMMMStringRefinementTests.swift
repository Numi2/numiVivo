import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct QMMMStringRefinementTests {
    private func path() -> VivoQMMMPathCoordinate {
        .init(identifier: "string-fixture",
            components: [
                .init(identifier: "q1", kind: .distance, atomIndices: [0,1]),
                .init(identifier: "q2", kind: .distance, atomIndices: [2,3])
            ], scalesNM: [0.1,0.1], nodes: [
                .init(identifier: "n0", valuesNM: [0.1,0.1]),
                .init(identifier: "n1", valuesNM: [0.2,0.2]),
                .init(identifier: "n2", valuesNM: [0.3,0.3]),
                .init(identifier: "n3", valuesNM: [0.4,0.4])
            ], lambda: 4)
    }

    @Test func purelyTangentialForcesLeaveUniformStraightStringUnchanged() throws {
        let source = path()
        let force = 10 / sqrt(2.0)
        let observations = ["n1","n2"].map {
            VivoQMMMStringNodeForce(nodeIdentifier: $0,
                meanForceKJPerMolNM: [force,force],
                standardErrorKJPerMolNM: [0.01,0.01], effectiveSamples: 200)
        }
        let result = try VivoQMMMStringRefinement.refine(path: source, nodeForces: observations,
            configuration: .init(mobilityPerKJPerMol: 0.01, smoothing: 0,
                                 maximumScaledNodeDisplacement: 1,
                                 minimumEffectiveSamples: 50,
                                 perpendicularForceToleranceKJPerMol: 1e-10))
        for (a,b) in zip(result.refinedPath.nodes, source.nodes) {
            for (x,y) in zip(a.valuesNM,b.valuesNM) { #expect(abs(x-y) < 1e-12) }
        }
        #expect(result.converged)
        #expect(result.maximumPerpendicularForceNormKJPerMol < 1e-12)
    }

    @Test func perpendicularForceMovesInteriorButKeepsEndpointsFixedAndReparameterizes() throws {
        let source = path()
        let observations = [
            VivoQMMMStringNodeForce(nodeIdentifier: "n1", meanForceKJPerMolNM: [-20,20],
                standardErrorKJPerMolNM: [0.1,0.1], effectiveSamples: 500),
            VivoQMMMStringNodeForce(nodeIdentifier: "n2", meanForceKJPerMolNM: [-20,20],
                standardErrorKJPerMolNM: [0.1,0.1], effectiveSamples: 500)
        ]
        let result = try VivoQMMMStringRefinement.refine(path: source, nodeForces: observations,
            configuration: .init(mobilityPerKJPerMol: 0.02, smoothing: 0,
                                 maximumScaledNodeDisplacement: 0.5,
                                 minimumEffectiveSamples: 50,
                                 perpendicularForceToleranceKJPerMol: 0.1))
        #expect(result.refinedPath.nodes.first == source.nodes.first)
        #expect(result.refinedPath.nodes.last == source.nodes.last)
        #expect(result.refinedPath.nodes[1].valuesNM != source.nodes[1].valuesNM)
        #expect(result.refinedPath.nodes[2].valuesNM != source.nodes[2].valuesNM)
        #expect(!result.converged)
        #expect(result.maximumPerpendicularForceNormKJPerMol > 0.1)
        let nodes = result.refinedPath.nodes.map(\.valuesNM)
        func scaledDistance(_ a:[Double],_ b:[Double])->Double {
            sqrt(zip(a,b).enumerated().reduce(0.0) { partial,item in
                let d=(item.element.0-item.element.1)/source.scalesNM[item.offset]
                return partial+d*d
            })
        }
        let lengths = (1..<nodes.count).map { scaledDistance(nodes[$0-1],nodes[$0]) }
        #expect((lengths.max() ?? 0) - (lengths.min() ?? 0) < 1e-10)
    }

    @Test func missingOrUndersampledNodeForceIsRejected() throws {
        let source = path()
        let insufficient = [
            VivoQMMMStringNodeForce(nodeIdentifier: "n1", meanForceKJPerMolNM: [1,0],
                standardErrorKJPerMolNM: [1,1], effectiveSamples: 2),
            VivoQMMMStringNodeForce(nodeIdentifier: "n2", meanForceKJPerMolNM: [1,0],
                standardErrorKJPerMolNM: [1,1], effectiveSamples: 100)
        ]
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMStringRefinement.refine(path: source, nodeForces: insufficient)
        }
        #expect(throws: (any Error).self) {
            _ = try VivoQMMMStringRefinement.refine(path: source, nodeForces: Array(insufficient.prefix(1)))
        }
    }
}
