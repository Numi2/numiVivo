import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellClusteringTests {
    static func graph() throws -> VivoSingleCellNeighborGraph {
        var options = VivoSingleCellNeighborOptions(); options.neighbors = 3
        return try VivoSingleCellNeighbors.run(scores: [[0],[1],[2],[100],[101],[102]],
            cells: (0..<6).map { .init(sampleID: "s",barcode: "c\($0)") },options: options)
    }
    @Test func separatedComponentsHaveAnalyticalModularity() throws {
        let graph = try Self.graph(), result = try VivoSingleCellClustering.run(graph,options: .init())
        #expect(result.labels == [0,0,0,1,1,1])
        #expect(result.clusterSizes == [3,3])
        #expect(abs(result.modularity-0.5) < 1e-12)
        #expect(result.disconnectedCommunities == 0)
        #expect(result == (try VivoSingleCellClustering.run(graph,options: .init())))
        #expect(zip(result.levels,result.levels.dropFirst()).allSatisfy { $1.modularity+1e-12 >= $0.modularity })
    }
    @Test func resolutionAndInsufficientIterationLimitsAreExplicit() throws {
        let graph = try Self.graph()
        var options = VivoSingleCellClusteringOptions(); options.resolution = 100
        let result = try VivoSingleCellClustering.run(graph,options: options)
        #expect(result.labels == Array(0..<6))
        options = .init(); options.maximumSweeps = 1
        #expect(throws: (any Error).self) { try VivoSingleCellClustering.run(graph,options: options) }
        options = .init(); options.maximumLevels = 1
        #expect(throws: (any Error).self) { try VivoSingleCellClustering.run(graph,options: options) }
        #expect(throws: (any Error).self) { try VivoSingleCellAnalysisPlan(id: "missing-graph",clustering: .init()).validate() }
        #expect(throws: (any Error).self) {
            try JSONDecoder().decode(VivoSingleCellClusteringOptions.self,from: Data("{\"resoluton\":1}".utf8))
        }
    }
    @Test func aggregationKeepsInternalEdgesInObjective() throws {
        let original: [[VivoSingleCellClustering.Edge]] = [
            [.init(column: 1,weight: 2),.init(column: 2,weight: 1)],
            [.init(column: 0,weight: 2),.init(column: 2,weight: 1)],
            [.init(column: 0,weight: 1),.init(column: 1,weight: 1)]]
        let aggregate: [[VivoSingleCellClustering.Edge]] = [
            [.init(column: 0,weight: 4),.init(column: 1,weight: 2)],
            [.init(column: 0,weight: 2)]]
        let q = VivoSingleCellClustering.modularity(original,labels: [0,0,1],resolution: 1)
        #expect(abs(q - (-0.125)) < 1e-12)
        #expect(q == VivoSingleCellClustering.modularity(aggregate,labels: [0,1],resolution: 1))
    }
}
