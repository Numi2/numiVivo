import Foundation

public struct VivoSingleCellIntegrationOptions: Codable, Sendable, Equatable {
    public enum Covariate: String, Codable, Sendable { case donor, batch }
    public var covariate: Covariate = .donor
    public var clusters: Int = 88
    public var diversity: Double = 2
    public var ridge: Double = 1
    public var temperature: Double = 0.1
    public var maximumIterations: Int = 10
    public var relativeTolerance: Double = 0.01
    public var seed: UInt64 = 7
    public var maximumWork: Int = 200_000_000
    public init() {}
    private enum CodingKeys: String, CodingKey { case covariate, clusters, diversity, ridge, temperature, maximumIterations, relativeTolerance, seed, maximumWork }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["covariate","clusters","diversity","ridge","temperature","maximumIterations","relativeTolerance","seed","maximumWork"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        covariate = try c.decodeIfPresent(Covariate.self,forKey: .covariate) ?? .donor
        clusters = try c.decodeIfPresent(Int.self,forKey: .clusters) ?? 88
        diversity = try c.decodeIfPresent(Double.self,forKey: .diversity) ?? 2
        ridge = try c.decodeIfPresent(Double.self,forKey: .ridge) ?? 1
        temperature = try c.decodeIfPresent(Double.self,forKey: .temperature) ?? 0.1
        maximumIterations = try c.decodeIfPresent(Int.self,forKey: .maximumIterations) ?? 10
        relativeTolerance = try c.decodeIfPresent(Double.self,forKey: .relativeTolerance) ?? 0.01
        seed = try c.decodeIfPresent(UInt64.self,forKey: .seed) ?? 7
        maximumWork = try c.decodeIfPresent(Int.self,forKey: .maximumWork) ?? 200_000_000
    }
    public func validate() throws {
        guard (2...100).contains(clusters), diversity.isFinite, (0...10).contains(diversity),
              ridge.isFinite, (0.001...100).contains(ridge), temperature.isFinite, (0.01...1).contains(temperature),
              (2...100).contains(maximumIterations), relativeTolerance.isFinite, (1e-8...0.05).contains(relativeTolerance),
              (1...1_000_000_000).contains(maximumWork) else { throw VivoOmicsError.invalid("integration options") }
    }
}

public struct VivoSingleCellIntegrationResult: Codable, Sendable, Equatable {
    public let method: String
    public let options: VivoSingleCellIntegrationOptions
    public let cells: [VivoOmicsCellIdentity]
    public let levels: [String]
    public let cellLevels: [Int]
    public let scores: [[Double]]
    /// Final soft memberships and pre-correction cosine-space inputs permit an
    /// independent objective and weighted ridge reconstruction.
    public let memberships: [[Double]]
    public let assignmentScores: [[Double]]
    public let assignmentCenters: [[Double]]
    public let objectives: [Double]
    public let relativeImprovements: [Double]
    public let stoppingReason: String
    public let maximumRidgeResidual: Double
    public let qualification: String
}

enum VivoSingleCellIntegration {
    /// Solve an intercept plus categorical ridge model through its Schur complement.
    /// The stable complement is sum(mass * lambda / (mass + lambda)); no subtraction
    /// of nearly equal total masses. Intercept is unpenalized and never removed.
    static func ridgeFit(masses: [Double], sums: [[Double]], ridge: Double) throws -> (intercept: [Double], effects: [[Double]], residual: Double) {
        let d = sums.first?.count ?? 0
        guard d > 0, masses.count == sums.count, masses.count >= 2, ridge.isFinite, ridge > 0,
              masses.allSatisfy({ $0.isFinite && $0 >= 0 }), sums.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration ridge inputs") }
        let denominator = masses.reduce(0) { $0 + $1 * ridge / ($1 + ridge) }
        guard denominator > 1e-14 else { throw VivoOmicsError.invalid("integration cluster has no effective mass") }
        var intercept = [Double](repeating: 0,count: d)
        for b in masses.indices { for j in 0..<d { intercept[j] += sums[b][j] * ridge / (masses[b] + ridge) / denominator } }
        var effects = sums, residual = 0.0
        for b in masses.indices { for j in 0..<d {
            effects[b][j] = (sums[b][j] - masses[b] * intercept[j]) / (masses[b] + ridge)
            let error = masses[b] * intercept[j] + (masses[b]+ridge) * effects[b][j] - sums[b][j]
            residual = max(residual,abs(error)/(1+abs(sums[b][j])))
        } }
        for j in 0..<d {
            var error = 0.0, scale = 1.0
            for b in masses.indices { error += masses[b]*(intercept[j]+effects[b][j])-sums[b][j]; scale += abs(sums[b][j]) }
            residual = max(residual,abs(error)/scale)
        }
        guard residual < 1e-10, intercept.allSatisfy(\.isFinite), effects.allSatisfy({ $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration ridge normal-equation residual") }
        return (intercept,effects,residual)
    }

    static func run(_ reduction: VivoSingleCellReductionResult, samples: [VivoOmicsSample], options o: VivoSingleCellIntegrationOptions) throws -> VivoSingleCellIntegrationResult {
        try o.validate()
        let x = reduction.scores, n = x.count, d = x.first?.count ?? 0, k = o.clusters
        guard n >= k, n == reduction.cells.count, (1...64).contains(d), x.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration PCA axes") }
        // Bounds memberships and arithmetic before allocating N*K, not N*genes.
        var work = n
        for factor in [k,d,o.maximumIterations+10] {
            let product = work.multipliedReportingOverflow(by: factor)
            guard !product.overflow, product.partialValue <= o.maximumWork else { throw VivoOmicsError.limit("integration work budget") }
            work = product.partialValue
        }
        guard Set(samples.map(\.id)).count == samples.count else { throw VivoOmicsError.invalid("duplicate integration sample identities") }
        let lookup = Dictionary(uniqueKeysWithValues: samples.map { ($0.id,$0) })
        let names = try reduction.cells.map { cell -> String in
            guard let sample = lookup[cell.sampleID], let value = o.covariate == .donor ? sample.donorID : sample.batchID,
                  !value.isEmpty, value != "unreported" else { throw VivoOmicsError.invalid("integration requires known selected covariate") }
            return value
        }
        let levels = Set(names).sorted(), bCount = levels.count
        guard (2...128).contains(bCount) else { throw VivoOmicsError.invalid("integration requires 2...128 covariate levels") }
        let ids = Dictionary(uniqueKeysWithValues: levels.enumerated().map { ($0.element,$0.offset) })
        let batch = names.map { ids[$0]! }
        // Reject disconnected designs in which condition and covariate cannot be
        // separated. This eligibility check does not feed labels into correction.
        var conditionLevels: [String:Set<Int>] = [:]
        for (i,cell) in reduction.cells.enumerated() { conditionLevels[lookup[cell.sampleID]!.condition,default: []].insert(batch[i]) }
        var reachable: Set<Int> = [0]
        for _ in 0..<bCount {
            let before = reachable.count
            for group in conditionLevels.values where !reachable.isDisjoint(with: group) { reachable.formUnion(group) }
            if reachable.count == before { break }
        }
        guard reachable.count == bCount else { throw VivoOmicsError.invalid("integration covariate is confounded with condition; no connected shared-condition design") }
        var sizes = [Double](repeating: 0,count: bCount)
        for b in batch { sizes[b] += 1 }
        var state = o.seed
        func uniform() -> Double {
            state &+= 0x9E3779B97F4A7C15; var z = state
            z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9; z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
            return Double((z ^ (z >> 31)) >> 11) / 9007199254740992.0
        }
        func normalized(_ rows: [[Double]]) -> [[Double]] {
            rows.map { row in let norm = sqrt(row.reduce(0) { $0+$1*$1 }); return norm > 0 ? row.map { $0/norm } : row }
        }
        func distance(_ a: [Double],_ b: [Double]) -> Double { zip(a,b).reduce(0) { $0+($1.0-$1.1)*($1.0-$1.1) } }
        let unit = normalized(x)
        // Independent seeded kmeans++ and Lloyd initialization. Reference uses a
        // different random generator/initializer; coordinate identity is not claimed.
        var centers = [unit[min(n-1,Int(uniform()*Double(n)))]]
        var nearest = [Double](repeating: .infinity,count: n)
        while centers.count < k {
            for i in 0..<n { nearest[i] = min(nearest[i],distance(unit[i],centers.last!)) }
            let total = nearest.reduce(0,+)
            guard total > 1e-14 else { throw VivoOmicsError.invalid("integration cluster count exceeds distinct PCA directions") }
            var threshold = uniform()*total, selected = n-1
            for i in 0..<n { threshold -= nearest[i]; if threshold <= 0 { selected = i; break } }
            centers.append(unit[selected])
        }
        for _ in 0..<10 {
            var sums = [[Double]](repeating: [Double](repeating: 0,count: d),count: k), counts = [Int](repeating: 0,count: k)
            for i in 0..<n {
                let c = (0..<k).min { distance(unit[i],centers[$0]) < distance(unit[i],centers[$1]) }!
                counts[c] += 1; for j in 0..<d { sums[c][j] += unit[i][j] }
            }
            for c in 0..<k where counts[c] > 0 { centers[c] = sums[c].map { $0/Double(counts[c]) } }
        }
        centers = normalized(centers)
        var scores = x, r = [Double](repeating: 0,count: n*k), distances = r
        var observed = [Double](repeating: 0,count: k*bCount), masses = [Double](repeating: 0,count: k)
        func setProbabilities(_ i: Int,_ diversity: Bool) {
            var logits = [Double](repeating: 0,count: k)
            for c in 0..<k {
                let expected = max(0,masses[c])*sizes[batch[i]]/Double(n), obs = max(0,observed[c*bCount+batch[i]])
                logits[c] = -distances[i*k+c]/o.temperature + (diversity ? o.diversity*log((2*expected+1)/(obs+expected+1)) : 0)
            }
            let peak = logits.max()!, weights = logits.map { exp($0-peak) }, total = weights.reduce(0,+)
            for c in 0..<k { r[i*k+c] = weights[c]/total }
        }
        func add(_ i: Int,_ sign: Double) { for c in 0..<k { let value = sign*r[i*k+c]; masses[c] += value; observed[c*bCount+batch[i]] += value } }
        func objective() -> Double {
            var result = 0.0
            for i in r.indices where r[i] > 0 { result += r[i]*distances[i]+o.temperature*r[i]*log(r[i]) }
            for c in 0..<k { for b in 0..<bCount {
                let obs = max(0,observed[c*bCount+b]), expected = max(0,masses[c])*sizes[b]/Double(n)
                result += o.temperature*o.diversity*obs*log((obs+expected+1)/(2*expected+1))
            } }
            return result*2000/Double(n)
        }
        var objectives: [Double] = [], improvements: [Double] = [], maxResidual = 0.0, stopping = "iteration-limit"
        var assignmentScores: [[Double]] = [], assignmentCenters: [[Double]] = []
        for iteration in 0..<o.maximumIterations {
            try Task.checkCancellation()
            let current = normalized(scores)
            assignmentScores = current; assignmentCenters = centers
            observed = observed.map { _ in 0 }; masses = masses.map { _ in 0 }
            for i in 0..<n {
                for c in 0..<k { distances[i*k+c] = 2*(1-zip(current[i],centers[c]).reduce(0) { $0+$1.0*$1.1 }) }
                setProbabilities(i,false); add(i,1)
            }
            if iteration == 0 { objectives.append(objective()) }
            // Four block-coordinate sweeps; block fraction .05, as the pinned reference.
            for _ in 0..<4 {
                var order = Array(0..<n)
                for i in stride(from: n-1,through: 1,by: -1) { order.swapAt(i,min(i,Int(uniform()*Double(i+1)))) }
                let blockSize = max(1,n/20)
                for block in 0..<20 {
                    let start = block*blockSize, end = block == 19 ? n : min(n,start+blockSize)
                    if start >= n { break }
                    for t in start..<end { add(order[t],-1) }
                    for t in start..<end { setProbabilities(order[t],true) }
                    for t in start..<end { add(order[t],1) }
                }
            }
            let value = objective()
            guard value.isFinite else { throw VivoOmicsError.invalid("integration objective nonfinite") }
            let previous = objectives.last!, improvement = (previous-value)/max(abs(previous),1e-12)
            objectives.append(value); improvements.append(improvement)
            // Every correction is fitted to and subtracted from original PCA.
            scores = x
            for c in 0..<k {
                let active = (0..<bCount).filter { observed[c*bCount+$0]/sizes[$0] > 1e-5 }
                if active.count < 2 { continue }
                var sums = [[Double]](repeating: [Double](repeating: 0,count: d),count: bCount)
                for i in 0..<n { for j in 0..<d { sums[batch[i]][j] += r[i*k+c]*x[i][j] } }
                let fit = try ridgeFit(masses: active.map { max(0,observed[c*bCount+$0]) },sums: active.map { sums[$0] },ridge: o.ridge)
                maxResidual = max(maxResidual,fit.residual); centers[c] = fit.intercept
                for (slot,b) in active.enumerated() { for i in 0..<n where batch[i] == b { for j in 0..<d { scores[i][j] -= r[i*k+c]*fit.effects[slot][j] } } }
            }
            centers = normalized(centers)
            guard scores.allSatisfy({ $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration corrected scores nonfinite") }
            if improvement < 0 { stopping = "objective-increase"; break }
            if improvement < o.relativeTolerance { stopping = "relative-objective-tolerance"; break }
        }
        return .init(method: "diversity-soft-clustering-categorical-ridge-Double-v1",options: o,cells: reduction.cells,levels: levels,cellLevels: batch,scores: scores,
            memberships: (0..<n).map { Array(r[($0*k)..<(($0+1)*k)]) },assignmentScores: assignmentScores,assignmentCenters: assignmentCenters,
            objectives: objectives,relativeImprovements: improvements,stoppingReason: stopping,maximumRidgeResidual: maxResidual,
            qualification: "Transductive single-covariate PCA correction with Harmony2 objective and intercept centers; independent initialization. Original counts/PCA preserved. Stopping diagnostics are not biological preservation or prospective prediction evidence.")
    }
}
