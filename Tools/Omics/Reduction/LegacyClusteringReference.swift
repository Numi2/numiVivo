// Frozen resident oracle from 8635cfb734f8ffb6562cba083a0f4347ece35811. Not a production implementation.
import Foundation
import NumiVivoKit

public struct VivoSingleCellClusteringLevel: Codable, Sendable, Equatable {
    public let vertices: Int
    public let communities: Int
    public let sweeps: Int
    public let moves: Int
    public let modularity: Double
}
public struct VivoSingleCellClusteringResult: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoSingleCellClusteringOptions
    public let cells: [VivoOmicsCellIdentity]
    public let labels: [Int]
    public let clusterSizes: [Int]
    public let modularity: Double
    public let levels: [VivoSingleCellClusteringLevel]
    public let disconnectedCommunities: Int
    public let termination: String
    public let qualification: String
}

enum VivoLegacyClusteringReference {
    struct Edge { let column: Int; let weight: Double }
    static func modularity(_ rows: [[Edge]], labels: [Int], resolution: Double) -> Double {
        var volume = [Double](repeating: 0,count: rows.count), internalWeight = 0.0, total = 0.0
        for i in rows.indices { for edge in rows[i] {
            total += edge.weight; volume[labels[i]] += edge.weight
            if labels[i] == labels[edge.column] { internalWeight += edge.weight }
        } }
        if total == 0 { return 0 }
        return internalWeight/total-resolution*volume.reduce(0) { $0+pow($1/total,2) }
    }
    static func run(_ graph: VivoSingleCellNeighborGraph,options: VivoSingleCellClusteringOptions) throws -> VivoSingleCellClusteringResult {
        try options.validate()
        let n = graph.cells.count
        var rows = (0..<n).map { i in
            (graph.rowOffsets[i]..<graph.rowOffsets[i+1]).map { Edge(column: graph.columnIndices[$0],weight: graph.weights[$0]) }
        }
        let original = rows
        var labels = Array(0..<n), levels: [VivoSingleCellClusteringLevel] = [], state = options.seed
        var previous = modularity(rows,labels: labels,resolution: options.resolution), termination: String?
        func random() -> UInt64 {
            state &+= 0x9E3779B97F4A7C15; var z = state
            z = (z^(z>>30)) &* 0xBF58476D1CE4E5B9
            z = (z^(z>>27)) &* 0x94D049BB133111EB
            return z^(z>>31)
        }
        for _ in 0..<options.maximumLevels {
            try Task.checkCancellation()
            let count = rows.count, degrees = rows.map { $0.reduce(0) { $0+$1.weight } }
            let total = degrees.reduce(0,+)
            if total == 0 { termination = "edgeless graph"; break }
            var community = Array(0..<count), totals = degrees, order = Array(0..<count)
            for i in stride(from: count-1,through: 1,by: -1) { order.swapAt(i,Int(random()%UInt64(i+1))) }
            var completed = false, moves = 0, sweeps = 0
            for _ in 0..<options.maximumSweeps {
                try Task.checkCancellation(); sweeps += 1; var changed = 0
                for i in order {
                    let old = community[i], degree = degrees[i]
                    var incident: [Int: Double] = [:]
                    for edge in rows[i] where edge.column != i { incident[community[edge.column],default: 0] += edge.weight }
                    totals[old] -= degree
                    let removal = -(incident[old] ?? 0)+options.resolution*degree*totals[old]/total
                    var best = old, bestGain = 1e-12
                    for candidate in incident.keys.sorted() where candidate != old {
                        let gain = 2*(removal+(incident[candidate] ?? 0)-options.resolution*degree*totals[candidate]/total)/total
                        if gain > bestGain { best = candidate; bestGain = gain }
                    }
                    totals[best] += degree
                    if best != old { community[i] = best; changed += 1 }
                }
                moves += changed
                if changed == 0 { completed = true; break }
            }
            guard completed else { throw VivoOmicsError.invalid("Louvain sweep limit before local convergence") }
            // Canonical community IDs follow the first original cell, not hash order.
            let sourceLabels = labels.map { community[$0] }
            var mapping: [Int: Int] = [:]
            for value in sourceLabels where mapping[value] == nil { mapping[value] = mapping.count }
            labels = sourceLabels.map { mapping[$0]! }
            let compact = community.map { mapping[$0]! }, groups = mapping.count
            let objective = modularity(original,labels: labels,resolution: options.resolution)
            guard objective.isFinite,objective+1e-10 >= previous else { throw VivoOmicsError.invalid("Louvain objective decreased") }
            levels.append(.init(vertices: count,communities: groups,sweeps: sweeps,moves: moves,modularity: objective))
            if groups == count || groups == 1 || objective-previous <= options.levelTolerance {
                termination = groups == count ? "no aggregate moves" : (groups == 1 ? "one community" : "level objective tolerance")
                previous = objective; break
            }
            previous = objective
            var aggregated = [[Int: Double]](repeating: [:],count: groups)
            for i in rows.indices { for edge in rows[i] {
                aggregated[compact[i]][compact[edge.column],default: 0] += edge.weight
            } }
            rows = aggregated.map { row in row.keys.sorted().map { Edge(column: $0,weight: row[$0]!) } }
        }
        guard let termination else { throw VivoOmicsError.invalid("Louvain level limit before convergence") }
        let clusters = (labels.max() ?? -1)+1
        var sizes = [Int](repeating: 0,count: clusters), components = sizes, visited = [Bool](repeating: false,count: n)
        for label in labels { sizes[label] += 1 }
        for i in 0..<n where !visited[i] {
            let label = labels[i]; components[label] += 1; visited[i] = true
            var queue = [i], cursor = 0
            while cursor < queue.count {
                let row = queue[cursor]; cursor += 1
                for edge in original[row] where labels[edge.column] == label && !visited[edge.column] {
                    visited[edge.column] = true; queue.append(edge.column)
                }
            }
        }
        return .init(method: "seeded-multilevel-Louvain-fuzzy-modularity-v1",options: options,cells: graph.cells,labels: labels,
            clusterSizes: sizes,modularity: previous,levels: levels,disconnectedCommunities: components.filter { $0>1 }.count,
            termination: termination,qualification: "Descriptive Louvain communities; connectivity diagnostics retained, not Leiden refinement, cell types, donor integration or biological validation")
    }
}
