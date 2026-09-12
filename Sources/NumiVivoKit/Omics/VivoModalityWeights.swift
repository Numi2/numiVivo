import Foundation

public struct VivoModalityWeightResult: Codable, Sendable {
    public let method: String
    public let cells: [VivoOmicsCellIdentity]
    public let neighborsIncludingSelf: Int
    public let bandwidths: [[Double]]
    public let withinDistances: [[Double]]
    public let crossDistances: [[Double]]
    public let scores: [[Double]]
    public let weights: [[Double]]
}

/// Paired row-L2-normalized representations. Exact kNN, SNN far-neighbor
/// bandwidth, no smoothing, cross constant 1e-4, scores clipped to [0,200].
/// Weights describe embedding predictability, not biological validity.
public enum VivoModalityWeights {
    public static func fit(first: [[Double]], second: [[Double]], cells: [VivoOmicsCellIdentity],
                           neighborsIncludingSelf k: Int = 15) throws -> VivoModalityWeightResult {
        try fitWithGraphs(first: first, second: second, cells: cells, neighborsIncludingSelf: k).result
    }
    static func fitWithGraphs(first: [[Double]], second: [[Double]], cells: [VivoOmicsCellIdentity],
                              neighborsIncludingSelf k: Int) throws -> (result: VivoModalityWeightResult, graphs: [VivoSingleCellNeighborGraph]) {
        let n = cells.count
        guard n >= k, Set(cells).count == n, cells.allSatisfy({ vivoOmicsID($0.sampleID) && vivoOmicsID($0.barcode) }) else {
            throw VivoOmicsError.invalid("modality cell identities")
        }
        let matrices = [first, second]
        for x in matrices {
            let d = x.first?.count ?? 0
            guard x.count == n, (1...64).contains(d), x.allSatisfy({ row in
                row.count == d && row.allSatisfy(\.isFinite) &&
                abs(row.reduce(0) { $0 + $1*$1 } - 1) <= 1e-8
            }) else { throw VivoOmicsError.invalid("modality requires paired finite unit rows") }
        }
        var options = VivoSingleCellNeighborOptions(); options.neighbors = k
        try options.validate()
        let graphs = try matrices.map { try VivoSingleCellNeighbors.run(scores: $0, cells: cells, options: options) }
        var widths = [[Double]](), within = [[Double]](), cross = [[Double]](), scores = [[Double]]()
        for m in 0..<2 {
            let x = matrices[m], graph = graphs[m], other = graphs[1-m]
            // Inverted membership avoids materializing a cells-by-cells SNN.
            var inverted = [[Int]](repeating: [], count: n)
            for i in 0..<n { for t in 0..<k { inverted[graph.neighborIndices[i*k+t]].append(i) } }
            var width = [Double](), wd = [Double](), cd = [Double](), score = [Double]()
            for i in 0..<n {
                try Task.checkCancellation()
                let nearest = graph.neighborDistances[i*k+1]
                var shared = [Int: Int]()
                for t in 0..<k {
                    for j in inverted[graph.neighborIndices[i*k+t]] { shared[j, default: 0] += 1 }
                }
                let ranked = shared.values.sorted(), count = min(k, ranked.count)
                guard count > 0 else { throw VivoOmicsError.invalid("empty modality SNN") }
                let cutoff = ranked[count-1]
                var far = [Double]()
                for (j, overlap) in shared where overlap <= cutoff {
                    var energy = 0.0
                    for f in x[i].indices { let delta = x[j][f]-x[i][f]; energy += delta*delta }
                    far.append(max(0, sqrt(energy)-nearest))
                }
                far.sort(by: >)
                let sigma = far.prefix(count).reduce(0, +)/Double(count)
                guard sigma.isFinite, sigma > 0 else { throw VivoOmicsError.invalid("zero modality bandwidth") }
                func error(_ indices: [Int]) -> Double {
                    var energy = 0.0
                    for f in x[i].indices {
                        var mean = 0.0
                        for t in 1..<k { mean += x[indices[i*k+t]][f] }
                        let delta = x[i][f]-mean/Double(k-1); energy += delta*delta
                    }
                    return max(0, sqrt(energy)-nearest)
                }
                let w = error(graph.neighborIndices), c = error(other.neighborIndices)
                let value = min(200, max(0, exp(-w/sigma)/(exp(-c/sigma)+1e-4)))
                guard value.isFinite else { throw VivoOmicsError.invalid("nonfinite modality score") }
                width.append(sigma); wd.append(w); cd.append(c); score.append(value)
            }
            widths.append(width); within.append(wd); cross.append(cd); scores.append(score)
        }
        var weights = [[Double]](repeating: [Double](repeating: 0, count: n), count: 2)
        for i in 0..<n {
            let shift = max(scores[0][i], scores[1][i])
            let a = exp(scores[0][i]-shift), b = exp(scores[1][i]-shift)
            weights[0][i] = a/(a+b); weights[1][i] = b/(a+b)
        }
        return (result: .init(method: "paired-two-modality-snn-predictability-v1", cells: cells,
            neighborsIncludingSelf: k, bandwidths: widths, withinDistances: within,
            crossDistances: cross, scores: scores, weights: weights), graphs: graphs)
    }
}
