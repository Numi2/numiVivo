import Foundation

public struct VivoWeightedNeighborResult: Codable, Sendable {
    public let method: String
    public let modalityWeights: VivoModalityWeightResult
    public let graph: VivoSingleCellNeighborGraph
}

/// Exact all-candidate two-modality neighbor selection; input conventions are
/// inherited from VivoModalityWeights. Distances are directed before fuzzy union.
public enum VivoWeightedNeighbors {
    public static func fit(first: [[Double]], second: [[Double]], cells: [VivoOmicsCellIdentity],
                           neighborsIncludingSelf k: Int = 15) throws -> VivoWeightedNeighborResult {
        let prepared = try VivoModalityWeights.fitWithGraphs(first: first, second: second, cells: cells, neighborsIncludingSelf: k)
        let weights = prepared.result, graphs = prepared.graphs
        let n = cells.count, matrices = [first, second]
        var heaps = [[VivoSingleCellNeighbors.Neighbor]](repeating: [], count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            for j in (i+1)..<n {
                var distances = [Double](repeating: 0, count: 2)
                for m in 0..<2 {
                    var energy = 0.0
                    for f in matrices[m][i].indices {
                        let delta = matrices[m][i][f]-matrices[m][j][f]; energy += delta*delta
                    }
                    distances[m] = sqrt(energy)
                }
                for (row, other) in [(i,j),(j,i)] {
                    var dissimilarity = 0.0
                    for m in 0..<2 {
                        let distance = max(0, distances[m]-graphs[m].neighborDistances[row*k+1])
                        dissimilarity += -expm1(-distance/weights.bandwidths[m][row])*weights.weights[m][row]
                    }
                    // Stable 1-exp(-x), avoiding cancellation near identical neighbors.
                    let squared = min(1, max(0, dissimilarity/2))
                    guard squared.isFinite else { throw VivoOmicsError.invalid("weighted neighbor affinity") }
                    VivoSingleCellNeighbors.retain(.init(index: other, squaredDistance: squared), in: &heaps[row], capacity: k-1)
                }
            }
        }
        var indices = [Int](), distances = [Double]()
        for i in 0..<n {
            indices.append(i); distances.append(0)
            for item in heaps[i].sorted(by: { $0.precedes($1) }) {
                indices.append(item.index); distances.append(sqrt(item.squaredDistance))
            }
        }
        let graph = try VivoSingleCellNeighbors.finish(indices: indices, distances: distances, cells: cells,
            dimensions: first[0].count+second[0].count, options: graphs[0].options, distancePairs: graphs[0].distancePairs,
            method: "exact-all-candidate-directed-modality-affinity-fuzzy-union-v1",
            qualification: "Exact paired modality affinity neighbors followed by UMAP fuzzy union; not Seurat approximate-search or SNN-graph parity, biological preservation, or prediction")
        return .init(method: "paired-two-modality-weighted-neighbors-v1", modalityWeights: weights, graph: graph)
    }
}
