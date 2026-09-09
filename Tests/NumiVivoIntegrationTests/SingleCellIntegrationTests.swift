import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellIntegrationTests {
    @Test func ridgeMatchesIndependentDenseNormalSolve() throws {
        let fit = try VivoSingleCellIntegration.ridgeFit(masses: [2,3,5],sums: [[4,2],[9,-3],[0,10]],ridge: 0.7)
        // NumPy linalg.solve on the full 4x4 intercept/dummy normal matrix.
        let expected = [[1.6114967462039047,0.6934924078091108],[0.28778018799710753,0.2270426608821401],[1.1258134490238607,-1.3731019522776575],[-1.413593637020969,1.146059291395517]]
        for (actual,row) in zip([fit.intercept]+fit.effects,expected) { for (a,b) in zip(actual,row) { #expect(abs(a-b)<1e-12) } }
        #expect(fit.residual<1e-12)
        let constant = try VivoSingleCellIntegration.ridgeFit(masses: [2,3,5],sums: [[14],[21],[35]],ridge: 1)
        #expect(abs(constant.intercept[0]-7)<1e-12)
        #expect(constant.effects.allSatisfy { abs($0[0])<1e-12 })
    }
    static func fixture() -> VivoSingleCellReductionResult {
        let cells: [VivoOmicsCellIdentity] = (0..<40).map { .init(sampleID: "s\($0%2)",barcode: "c\($0)") }
        let scores: [[Double]] = (0..<40).map { i in
            let t = Double(i)
            return [cos(t/5)+Double(i%2),sin(t/5),t/40]
        }
        return .init(method: "numerical-fixture",options: .init(),features: [],selectedFeatureIndices: [],
            cells: cells,scores: scores,
            loadings: [],explainedVariance: [],explainedVarianceRatio: [],relativeResiduals: [],maximumLoadingOrthogonalityError: 0,basisSize: 0,qualification: "numerical fixture only")
    }
    @Test func correctionIsReproducibleAndRejectsUnavailableMetadata() throws {
        let input = Self.fixture()
        let samples = (0..<2).map { VivoOmicsSample(id: "s\($0)",biologicalReplicateID: "d\($0)",donorID: "d\($0)",condition: "c",batchID: "unreported",organism: "human") }
        var options = VivoSingleCellIntegrationOptions(); options.clusters=3
        let result = try VivoSingleCellIntegration.run(input,samples: samples,options: options)
        #expect(result == (try VivoSingleCellIntegration.run(input,samples: samples,options: options)))
        #expect(result.cells == input.cells && result.scores != input.scores)
        #expect(result.maximumRidgeResidual<1e-10)
        #expect(result.objectives.count == result.relativeImprovements.count+1)
        let confounded = (0..<2).map { VivoOmicsSample(id: "s\($0)",biologicalReplicateID: "d\($0)",donorID: "d\($0)",condition: "c\($0)",batchID: "unreported",organism: "human") }
        #expect(throws: (any Error).self) { try VivoSingleCellIntegration.run(input,samples: confounded,options: options) }
        options.covariate = .batch
        #expect(throws: (any Error).self) { try VivoSingleCellIntegration.run(input,samples: samples,options: options) }
        options.covariate = .donor; options.maximumWork=1
        #expect(throws: (any Error).self) { try VivoSingleCellIntegration.run(input,samples: samples,options: options) }
    }
    @Test func representationDependenciesAreExplicit() throws {
        #expect(throws: (any Error).self) { try VivoSingleCellAnalysisPlan(id: "missing-pca",integration: .init()).validate() }
        var neighbor = VivoSingleCellNeighborOptions(); neighbor.representation = .integrated
        #expect(throws: (any Error).self) { try VivoSingleCellAnalysisPlan(id: "missing-integration",reduction: .init(),neighbors: neighbor).validate() }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoSingleCellIntegrationOptions.self,from: Data("{\"lamda\":1}".utf8)) }
    }
}
