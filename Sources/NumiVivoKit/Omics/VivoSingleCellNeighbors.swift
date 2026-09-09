import Foundation

public struct VivoSingleCellNeighborOptions: Codable, Sendable, Equatable {
    public enum Representation: String, Codable, Sendable { case pca, integrated }
    public var representation: Representation? = nil
    /// Includes the cell itself in slot zero, as in UMAP's reference convention.
    public var neighbors: Int = 15
    public var localConnectivity: Double = 1
    public var maximumDistancePairs: Int = 50_000_000
    public init() {}
    private enum CodingKeys: String, CodingKey { case neighbors, localConnectivity, maximumDistancePairs, representation }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["neighbors", "localConnectivity", "maximumDistancePairs", "representation"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        representation = try c.decodeIfPresent(Representation.self,forKey: .representation)
        neighbors = try c.decodeIfPresent(Int.self, forKey: .neighbors) ?? 15
        localConnectivity = try c.decodeIfPresent(Double.self, forKey: .localConnectivity) ?? 1
        maximumDistancePairs = try c.decodeIfPresent(Int.self, forKey: .maximumDistancePairs) ?? 50_000_000
    }
    public func validate() throws {
        guard (2...128).contains(neighbors), localConnectivity.isFinite,
              localConnectivity >= 0, localConnectivity <= Double(neighbors-1),
              (1...500_000_000).contains(maximumDistancePairs) else {
            throw VivoOmicsError.invalid("neighbor options")
        }
    }
}
public struct VivoSingleCellNeighborGraph: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoSingleCellNeighborOptions
    public let cells: [VivoOmicsCellIdentity]
    public let dimensions: Int
    /// Cell-major fixed-width kNN arrays. Self is first, then distance/index order.
    public let neighborIndices: [Int]
    public let neighborDistances: [Double]
    public let rhos: [Double]
    public let sigmas: [Double]
    /// Absolute difference from log2(k), after the minimum bandwidth floor.
    /// Duplicate points can make the target mass unattainable; retain this fact.
    public let kernelMassResiduals: [Double]
    /// Symmetric fuzzy-union connectivity CSR, sorted by column, no self edges.
    public let rowOffsets: [Int]
    public let columnIndices: [Int]
    public let weights: [Double]
    public let connectedComponents: Int
    public let isolatedCells: Int
    public let distancePairs: Int
    public let qualification: String
}

enum VivoSingleCellNeighbors {
    struct Neighbor {
        let index: Int
        let squaredDistance: Double
        func precedes(_ other: Self) -> Bool {
            squaredDistance == other.squaredDistance ? index < other.index : squaredDistance < other.squaredDistance
        }
    }
    /// A bounded max heap retains exact nearest neighbors without a pair matrix.
    static func retain(_ value: Neighbor, in heap: inout [Neighbor], capacity: Int) {
        if heap.count < capacity {
            heap.append(value); var i = heap.count-1
            while i > 0 {
                let p = (i-1)/2
                if !heap[p].precedes(heap[i]) { break }
                heap.swapAt(i,p); i = p
            }
        } else if value.precedes(heap[0]) {
            heap[0] = value; var i = 0
            while 2*i+1 < heap.count {
                var child = 2*i+1
                if child+1 < heap.count, heap[child].precedes(heap[child+1]) { child += 1 }
                if !heap[i].precedes(heap[child]) { break }
                heap.swapAt(i,child); i = child
            }
        }
    }
    static func run(scores: [[Double]], cells: [VivoOmicsCellIdentity], options: VivoSingleCellNeighborOptions) throws -> VivoSingleCellNeighborGraph {
        try options.validate()
        let n = scores.count, k = options.neighbors, d = scores.first?.count ?? 0
        guard n == cells.count, n >= k, (1...64).contains(d),
              scores.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }) else {
            throw VivoOmicsError.invalid("neighbor scores or cell axes")
        }
        // Divide before multiplying to avoid overflow on untrusted dimensions.
        let pairProduct = (n % 2 == 0 ? n/2 : n).multipliedReportingOverflow(by: n % 2 == 0 ? n-1 : (n-1)/2)
        guard !pairProduct.overflow, pairProduct.partialValue <= options.maximumDistancePairs else {
            throw VivoOmicsError.limit("exact neighbor distance-pair budget; approximate search is not yet implemented")
        }
        var heaps = [[Neighbor]](repeating: [],count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            for j in (i+1)..<n {
                var distance = 0.0
                for f in 0..<d { let delta = scores[i][f]-scores[j][f]; distance += delta*delta }
                guard distance.isFinite else { throw VivoOmicsError.invalid("nonfinite neighbor distance") }
                retain(.init(index: j,squaredDistance: distance),in: &heaps[i],capacity: k-1)
                retain(.init(index: i,squaredDistance: distance),in: &heaps[j],capacity: k-1)
            }
        }
        var indices: [Int] = [], distances: [Double] = []
        indices.reserveCapacity(n*k); distances.reserveCapacity(n*k)
        for i in 0..<n {
            indices.append(i); distances.append(0)
            for neighbor in heaps[i].sorted(by: { $0.precedes($1) }) {
                indices.append(neighbor.index); distances.append(sqrt(neighbor.squaredDistance))
            }
        }
        return try finish(indices: indices, distances: distances, cells: cells, dimensions: d, options: options, distancePairs: pairProduct.partialValue)
    }
    /// Shared UMAP-compatible graph construction for resident and file-backed
    /// exact searches. Search execution settings do not alter this calculation.
    static func finish(indices: [Int], distances: [Double], cells: [VivoOmicsCellIdentity], dimensions d: Int,
                       options: VivoSingleCellNeighborOptions, distancePairs: Int) throws -> VivoSingleCellNeighborGraph {
        let n = cells.count, k = options.neighbors
        guard indices.count == n*k, distances.count == n*k,
              indices.allSatisfy({ $0 >= 0 && $0 < n }), distances.allSatisfy({ $0.isFinite && $0 >= 0 }),
              (0..<n).allSatisfy({ indices[$0*k] == $0 && distances[$0*k] == 0 }) else {
            throw VivoOmicsError.invalid("neighbor graph arrays")
        }
        let globalMean = distances.reduce(0,+)/Double(distances.count), target = log2(Double(k))
        var rhos = [Double](repeating: 0,count: n), sigmas = rhos, residuals = rhos
        var directed = [[Int: Double]](repeating: [:],count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            let row = Array(distances[(i*k)..<((i+1)*k)]), positive = row.filter { $0 > 0 }
            let local = options.localConnectivity, integer = Int(floor(local)), fraction = local-Double(integer)
            var rho = 0.0
            if Double(positive.count) >= local {
                if integer > 0 {
                    rho = positive[integer-1]
                    if fraction > 1e-5 { rho += fraction*(positive[integer]-positive[integer-1]) }
                } else if let first = positive.first { rho = fraction*first }
            } else { rho = positive.last ?? 0 }
            func mass(_ sigma: Double) -> Double {
                row.dropFirst().reduce(0) { $0 + ($1 <= rho ? 1 : exp(-($1-rho)/sigma)) }
            }
            var low = 0.0, high = Double.infinity, sigma = 1.0
            for _ in 0..<64 {
                let sum = mass(sigma)
                if abs(sum-target) < 1e-5 { break }
                if sum > target { high = sigma; sigma = (low+high)/2 }
                else { low = sigma; sigma = high.isInfinite ? sigma*2 : (low+high)/2 }
            }
            let mean = rho > 0 ? row.reduce(0,+)/Double(k) : globalMean
            sigma = max(sigma, 1e-3*mean)
            guard rho.isFinite, sigma.isFinite, sigma > 0 else { throw VivoOmicsError.invalid("nonfinite graph bandwidth") }
            rhos[i] = rho; sigmas[i] = sigma; residuals[i] = abs(mass(sigma)-target)
            for slot in 1..<k {
                let delta = row[slot]-rho
                let weight = delta <= 0 ? 1 : exp(-delta/sigma)
                if weight > 0 { directed[i][indices[i*k+slot]] = weight }
            }
        }
        var symmetric = [[Int: Double]](repeating: [:],count: n)
        for i in 0..<n {
            try Task.checkCancellation()
            for (j,forward) in directed[i] {
                let reverse = directed[j][i] ?? 0
                let weight = min(1,max(0,forward+reverse-forward*reverse))
                symmetric[i][j] = weight; symmetric[j][i] = weight
            }
        }
        var offsets = [0], columns: [Int] = [], weights: [Double] = []
        for row in symmetric {
            for j in row.keys.sorted() { columns.append(j); weights.append(row[j]!) }
            offsets.append(columns.count)
        }
        var visited = [Bool](repeating: false,count: n), components = 0
        for i in 0..<n where !visited[i] {
            components += 1; var queue = [i], cursor = 0; visited[i] = true
            while cursor < queue.count {
                let row = queue[cursor]; cursor += 1
                for edge in offsets[row]..<offsets[row+1] {
                    let j = columns[edge]
                    if !visited[j] { visited[j] = true; queue.append(j) }
                }
            }
        }
        return .init(method: options.representation == .integrated ? "exact-euclidean-integrated-knn-umap-fuzzy-union-v1" : "exact-euclidean-PCA-knn-umap-fuzzy-union-v1",options: options,cells: cells,dimensions: d,
            neighborIndices: indices,neighborDistances: distances,rhos: rhos,sigmas: sigmas,kernelMassResiduals: residuals,
            rowOffsets: offsets,columnIndices: columns,weights: weights,connectedComponents: components,
            isolatedCells: (0..<n).filter { offsets[$0] == offsets[$0+1] }.count,distancePairs: distancePairs,
            qualification: "Exact neighbors in the declared representation and UMAP-compatible fuzzy graph; not an embedding, clustering, integration or biological qualification")
    }
}
