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
                standardErrorKJPerMolNM: [0,0], effectiveSamples: 200)
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
        // Equal distances ALONG the provisional polyline do not imply equal
        // Euclidean chords between resampled points around a corner.
        let outer = sqrt(0.96*0.96 + 1.04*1.04), middle = sqrt(2.0)
        let spacing = (2*outer + middle)/3
        let t = spacing/outer
        let expected1 = [0.1*(1+0.96*t), 0.1*(1+1.04*t)]
        let u = (2*spacing - outer - middle)/outer
        let expected2 = [0.1*(2.96+1.04*u), 0.1*(3.04+0.96*u)]
        for d in 0..<2 {
            #expect(abs(result.refinedPath.nodes[1].valuesNM[d] - expected1[d]) < 1e-12)
            #expect(abs(result.refinedPath.nodes[2].valuesNM[d] - expected2[d]) < 1e-12)
        }
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
    @Test func zeroMeanForceWithLargeUncertaintyDoesNotDeclareConvergence() throws {
        let source = path()
        let observations = ["n1","n2"].map {
            VivoQMMMStringNodeForce(nodeIdentifier: $0, meanForceKJPerMolNM: [0,0],
                standardErrorKJPerMolNM: [10,10], effectiveSamples: 500)
        }
        let result = try VivoQMMMStringRefinement.refine(path: source, nodeForces: observations,
            configuration: .init(smoothing: 0, perpendicularForceToleranceKJPerMol: 0.1))
        #expect(result.maximumPerpendicularForceNormKJPerMol == 0)
        #expect(result.maximumUncertaintyGuardedPerpendicularForceNormKJPerMol > 2.8)
        #expect(!result.converged)
        #expect(abs(result.diagnostics[0].perpendicularForceStandardDeviationUpperBoundKJPerMol - sqrt(2)) < 1e-12)
    }

    @Test func reparameterizationMotionIsPartOfConvergenceAndOldSchemaIsRejected() throws {
        var source = path(); source.nodes[1].valuesNM = [0.12,0.12]
        let observations = ["n1","n2"].map {
            VivoQMMMStringNodeForce(nodeIdentifier: $0, meanForceKJPerMolNM: [0,0],
                standardErrorKJPerMolNM: [0,0], effectiveSamples: 500)
        }
        let result = try VivoQMMMStringRefinement.refine(path: source, nodeForces: observations,
            configuration: .init(smoothing: 0))
        #expect(result.maximumUncertaintyGuardedPerpendicularForceNormKJPerMol == 0)
        #expect(result.maximumScaledDisplacement > 0.5 && !result.converged)
        var old = VivoQMMMStringRefinementConfiguration(); old.schema = "numivivo.org/qmmm-string-refinement/v1"
        #expect(throws: (any Error).self) { try old.validate() }
    }

}
