import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNeighborTests {
    static func cells(_ n: Int) -> [VivoOmicsCellIdentity] {
        (0..<n).map { .init(sampleID: "s",barcode: "c\($0)") }
    }
    @Test func tiesAndDuplicateCellsPreserveSelfAndDeterminism() throws {
        let scores = [[0.0,0],[0,0],[1,0],[-1,0],[0,1],[0,-1]]
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 3
        let graph = try VivoSingleCellNeighbors.run(scores: scores,cells: Self.cells(6),options: options)
        #expect(graph == (try VivoSingleCellNeighbors.run(scores: scores,cells: Self.cells(6),options: options)))
        #expect(Array(graph.neighborIndices[0..<3]) == [0,1,2])
        #expect(Array(graph.neighborIndices[3..<6]) == [1,0,2])
        #expect(graph.neighborDistances[0..<3].elementsEqual([0,0,1]))
        #expect(graph.connectedComponents == 1 && graph.isolatedCells == 0)
        #expect(graph.distancePairs == 15)
        for i in 0..<6 {
            for edge in graph.rowOffsets[i]..<graph.rowOffsets[i+1] {
                let j = graph.columnIndices[edge]
                let reverse = try #require((graph.rowOffsets[j]..<graph.rowOffsets[j+1]).first { graph.columnIndices[$0] == i })
                #expect(graph.weights[edge] == graph.weights[reverse])
                #expect(i != j && graph.weights[edge] > 0 && graph.weights[edge] <= 1)
            }
        }
        let duplicate = try VivoSingleCellNeighbors.run(scores: Array(repeating: [0.0],count: 6),cells: Self.cells(6),options: options)
        #expect(duplicate.weights.allSatisfy { $0 == 1 })
        #expect(duplicate.kernelMassResiduals.allSatisfy { abs($0-(2-log2(3.0))) < 1e-12 })
    }
    @Test func heapRecoversKnownNearestNeighborsAndDisconnectedGraph() throws {
        let scores = [[0.0],[1],[3],[100],[101],[103]]
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 2
        let graph = try VivoSingleCellNeighbors.run(scores: scores,cells: Self.cells(6),options: options)
        #expect(graph.neighborIndices == [0,1,1,0,2,1,3,4,4,3,5,4])
        #expect(graph.neighborDistances == [0,1,0,1,0,2,0,1,0,1,0,2])
        #expect(graph.connectedComponents == 2 && graph.isolatedCells == 0)
    }
    @Test func invalidInputsAndWorkBudgetFailBeforePublication() throws {
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 2; options.maximumDistancePairs = 1
        #expect(throws: (any Error).self) {
            try VivoSingleCellNeighbors.run(scores: [[0],[1],[2]],cells: Self.cells(3),options: options)
        }
        options.maximumDistancePairs = 10
        #expect(throws: (any Error).self) {
            try VivoSingleCellNeighbors.run(scores: [[0],[Double.nan]],cells: Self.cells(2),options: options)
        }
        #expect(throws: (any Error).self) { try VivoSingleCellAnalysisPlan(id: "missing-pca",neighbors: options).validate() }
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellNeighborOptions.self,from: Data("{\"neigbors\":15}".utf8))
        }
    }
}
