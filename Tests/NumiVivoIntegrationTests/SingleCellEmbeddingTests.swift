import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellEmbeddingTests {
    @Test func curveFitMatchesIndependentDefaultAndScaleCases() throws {
        let fit = try VivoSingleCellEmbedding.fitCurve(spread: 1,minimumDistance: 0.1)
        #expect(abs(fit.a-1.5769434605754993) < 2e-5)
        #expect(abs(fit.b-0.8950608781680347) < 2e-5)
        for (spread,distance) in [(0.01,0.0),(10.0,10.0),(1.0,0.0)] {
            let fit = try VivoSingleCellEmbedding.fitCurve(spread: spread,minimumDistance: distance)
            #expect(fit.a.isFinite && fit.a>0 && fit.b.isFinite && fit.b>0 && fit.error.isFinite)
        }
    }
    @Test func attractionAndRegularizedRepulsionMatchFiniteDifferences() {
        let a=1.5769434605754993,b=0.8950608781680347,x=0.7,y = -0.2,h=1e-6
        let squared=x*x+y*y
        func attractiveLoss(_ x: Double) -> Double { log1p(a*pow(x*x+y*y,b)) }
        let expected = -(attractiveLoss(x+h)-attractiveLoss(x-h))/(2*h)
        #expect(abs(VivoSingleCellEmbedding.attractiveCoefficient(squaredDistance: squared,a: a,b: b)*x-expected) < 1e-8)
        func repulsiveLoss(_ x: Double) -> Double { -log(1-1/(1+a*pow(x*x+y*y,b))) }
        let raw = -(repulsiveLoss(x+h)-repulsiveLoss(x-h))/(2*h)
        let regularized = raw*squared/(0.001+squared)
        #expect(abs(VivoSingleCellEmbedding.repulsiveCoefficient(squaredDistance: squared,a: a,b: b,strength: 1)*x-regularized) < 1e-8)
    }
    @Test func sparseOptimizationIsDeterministicAndBudgeted() throws {
        let scores = (0..<12).map { i in [cos(Double(i)/2),sin(Double(i)/2),Double(i)/10] }
        var neighbor = VivoSingleCellNeighborOptions(); neighbor.neighbors=5
        let graph = try VivoSingleCellNeighbors.run(scores: scores,cells: (0..<12).map { .init(sampleID: "s",barcode: "c\($0)") },options: neighbor)
        var options = VivoSingleCellEmbeddingOptions(); options.epochs=20
        let result = try VivoSingleCellEmbedding.run(graph,scores: scores,options: options)
        #expect(result == (try VivoSingleCellEmbedding.run(graph,scores: scores,options: options)))
        #expect(result.coordinates != result.initialCoordinates)
        #expect(result.coordinates.allSatisfy { $0.count==2 && $0.allSatisfy(\.isFinite) })
        #expect(result.attractiveUpdates>0 && result.negativeSamples>0)
        #expect(result.edgeVisits==result.retainedDirectedEdges*options.epochs)
        options.maximumUpdates=1
        #expect(throws: (any Error).self) { try VivoSingleCellEmbedding.run(graph,scores: scores,options: options) }
        #expect(throws: (any Error).self) { try VivoSingleCellAnalysisPlan(id: "missing-graph",embedding: .init()).validate() }
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellEmbeddingOptions.self,from: Data("{\"minimumDistnace\":0.1}".utf8))
        }
    }
}
