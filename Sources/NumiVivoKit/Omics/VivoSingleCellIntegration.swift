import Foundation

public struct VivoSingleCellIntegrationOptions: Codable, Sendable, Equatable {
    static let maximumClusters = 100
    public enum Covariate: String, Codable, Sendable, CaseIterable { case donor, batch }
    public enum RidgeScaling: String, Codable, Sendable { case expectedClusterBatchMass }
    /// Nil preserves the original fixed penalty; otherwise ridge is the expected-mass coefficient.
    public var ridgeScaling: RidgeScaling? = nil
    /// Optional, complete observation-sample map for an explicitly protected
    /// categorical group. Nil preserves historical encoding and behavior.
    public var protectedSampleGroups: [String: String]? = nil
    public var covariate: Covariate = .donor
    /// Optional joint correction factors. Omission preserves the historical
    /// single-covariate representation and output bytes. When supplied, use
    /// both donor and batch exactly once; the correction is additive across
    /// their categorical effects rather than an interaction-level table.
    public var covariates: [Covariate]? = nil
    public var clusters: Int = 88
    public var diversity: Double = 2
    public var ridge: Double = 1
    public var temperature: Double = 0.1
    public var maximumIterations: Int = 10
    public var relativeTolerance: Double = 0.01
    public var seed: UInt64 = 7
    public var maximumWork: Int = 200_000_000
    public init() {}
    private enum CodingKeys: String, CodingKey { case protectedSampleGroups, covariate, covariates, clusters, diversity, ridge, ridgeScaling, temperature, maximumIterations, relativeTolerance, seed, maximumWork }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["protectedSampleGroups","covariate","covariates","clusters","diversity","ridge","ridgeScaling","temperature","maximumIterations","relativeTolerance","seed","maximumWork"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        protectedSampleGroups = try c.decodeIfPresent([String: String].self, forKey: .protectedSampleGroups)
        covariate = try c.decodeIfPresent(Covariate.self,forKey: .covariate) ?? .donor
        covariates = try c.decodeIfPresent([Covariate].self, forKey: .covariates)
        clusters = try c.decodeIfPresent(Int.self,forKey: .clusters) ?? 88
        diversity = try c.decodeIfPresent(Double.self,forKey: .diversity) ?? 2
        ridge = try c.decodeIfPresent(Double.self,forKey: .ridge) ?? 1
        ridgeScaling = try c.decodeIfPresent(RidgeScaling.self, forKey: .ridgeScaling)
        temperature = try c.decodeIfPresent(Double.self,forKey: .temperature) ?? 0.1
        maximumIterations = try c.decodeIfPresent(Int.self,forKey: .maximumIterations) ?? 10
        relativeTolerance = try c.decodeIfPresent(Double.self,forKey: .relativeTolerance) ?? 0.01
        seed = try c.decodeIfPresent(UInt64.self,forKey: .seed) ?? 7
        maximumWork = try c.decodeIfPresent(Int.self,forKey: .maximumWork) ?? 200_000_000
    }
    public func validate() throws {
        try VivoIntegrationProtection.validate(protectedSampleGroups)
        if let covariates {
            guard covariates.count == 2, Set(covariates).count == covariates.count,
                  Set(covariates) == Set(Covariate.allCases) else {
                throw VivoOmicsError.invalid("joint integration requires donor and batch exactly once")
            }
        }
        guard (2...Self.maximumClusters).contains(clusters), diversity.isFinite, (0...10).contains(diversity),
              ridge.isFinite, (0.001...100).contains(ridge), temperature.isFinite, (0.01...1).contains(temperature),
              (2...100).contains(maximumIterations), relativeTolerance.isFinite, (1e-8...0.05).contains(relativeTolerance),
              (1...100_000_000_000).contains(maximumWork) else { throw VivoOmicsError.invalid("integration options") }
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
    /// Final cluster-by-level penalties; absent for the original fixed-ridge mode.
    public var ridgePenalties: [[Double]]? = nil
    /// Factor-major levels and cell assignments for the opt-in joint path.
    /// Nil preserves the historical single-covariate result representation.
    public var factorLevels: [[String]]? = nil
    public var cellFactorLevels: [[Int]]? = nil
    public var factorRidgePenalties: [[[Double]]]? = nil
}


/// Categorical design admission only. Connectivity does not qualify preservation
/// of expression programs, cell identities, or any undeclared biological factor.
enum VivoIntegrationProtection {
    static func validate(_ groups: [String: String]?) throws {
        guard let groups else { return }
        guard !groups.isEmpty, groups.count <= 4_096,
              groups.allSatisfy({ vivoOmicsID($0.key) && vivoOmicsID($0.value) && $0.value != "unreported" }),
              try VivoCanonicalJSON.encode(groups).count <= 32_768 else {
            throw VivoOmicsError.invalid("protected integration sample groups or byte budget")
        }
    }
    static func check(_ groups: [String: String]?, cells: [VivoOmicsCellIdentity],
                      samples: [String: VivoOmicsSample],
                      covariate: VivoSingleCellIntegrationOptions.Covariate,
                      levelIDs: [String: Int]) throws {
        guard let groups else { return }
        try validate(groups)
        let observed = Set(cells.lazy.map(\.sampleID))
        guard Set(groups.keys) == observed else {
            throw VivoOmicsError.invalid("protected integration groups must match every observed sample exactly")
        }
        var strata: [[String]: Set<Int>] = [:]
        for id in observed.sorted() {
            try Task.checkCancellation()
            guard let sample = samples[id], let group = groups[id],
                  let name = covariate == .donor ? sample.donorID : sample.batchID,
                  let level = levelIDs[name] else {
                throw VivoOmicsError.invalid("protected integration sample or covariate identity")
            }
            // An array key keeps condition/group pairs distinct even if labels
            // contain punctuation. Protect their joint strata, not just marginals.
            strata[[sample.condition, group], default: []].insert(level)
        }
        var reachable: Set<Int> = [0]
        for _ in 0..<levelIDs.count {
            let before = reachable.count
            for levels in strata.values where !reachable.isDisjoint(with: levels) {
                reachable.formUnion(levels)
            }
            if reachable.count == before { break }
        }
        guard reachable.count == levelIDs.count else {
            throw VivoOmicsError.invalid("integration covariate is confounded with protected condition/group strata")
        }
    }
    static func check(_ groups: [String: String]?, cells: [VivoOmicsCellIdentity],
                      samples: [String: VivoOmicsSample],
                      covariates: [VivoSingleCellIntegrationOptions.Covariate],
                      levelIDs: [[String: Int]]) throws {
        guard let groups else { return }
        try validate(groups)
        let observed = Set(cells.lazy.map(\.sampleID))
        guard Set(groups.keys) == observed, levelIDs.count == covariates.count else {
            throw VivoOmicsError.invalid("protected integration groups or joint covariate axes")
        }
        for factor in covariates.indices {
            var strata: [[String]: Set<Int>] = [:]
            for id in observed.sorted() {
                try Task.checkCancellation()
                guard let sample = samples[id], let group = groups[id],
                      let name = covariates[factor] == .donor ? sample.donorID : sample.batchID,
                      let level = levelIDs[factor][name] else {
                    throw VivoOmicsError.invalid("protected integration sample or joint covariate identity")
                }
                strata[[sample.condition, group], default: []].insert(level)
            }
            var reachable: Set<Int> = [0]
            for _ in 0..<levelIDs[factor].count {
                let before = reachable.count
                for levels in strata.values where !reachable.isDisjoint(with: levels) { reachable.formUnion(levels) }
                if reachable.count == before { break }
            }
            guard reachable.count == levelIDs[factor].count else {
                throw VivoOmicsError.invalid("joint integration covariate is confounded with protected condition/group strata")
            }
        }
    }
}

enum VivoSingleCellIntegration {
    /// One sequential read of the latent and membership matrices. Each scalar
    /// accumulator still sees cells in exactly the original ascending order.
    /// The table is bounded by clusters * levels * components, never cell count.
    static func streamedRidgeSums(x: VivoIntegrationMatrix, memberships: VivoIntegrationMatrix,
                                 batch: [Int], levels: Int, activeClusters: [Int]) throws -> [Double] {
        let n = x.rows, d = x.columns, k = memberships.columns
        guard memberships.rows == n, batch.count == n, (1...64).contains(d),
              (2...100).contains(k), (2...128).contains(levels),
              batch.allSatisfy({ (0..<levels).contains($0) }),
              activeClusters.allSatisfy({ (0..<k).contains($0) }),
              Set(activeClusters).count == activeClusters.count else { throw VivoOmicsError.invalid("streamed ridge axes") }
        var sums = [Double](repeating: 0, count: k * levels * d)
        for i in 0..<n {
            if i % 2_048 == 0 { try Task.checkCancellation() }
            let row = try x.row(i), weights = try memberships.row(i)
            for c in activeClusters {
                let weight = weights[c], offset = (c * levels + batch[i]) * d
                for j in 0..<d { sums[offset + j] += weight * row[j] }
            }
        }
        return sums
    }

    /// Each cell receives effects in ascending cluster order, exactly as in the
    /// cluster-major implementation. Nil entries skip inactive levels entirely.
    static func applyStreamedRidgeEffects(scores: VivoIntegrationMatrix, memberships: VivoIntegrationMatrix,
                                         batch: [Int], levels: Int, effects: [[Double]?]) throws {
        let n = scores.rows, d = scores.columns, k = memberships.columns
        guard memberships.rows == n, batch.count == n, (1...64).contains(d),
              (2...100).contains(k), (2...128).contains(levels), effects.count == k * levels,
              batch.allSatisfy({ (0..<levels).contains($0) }),
              effects.allSatisfy({ $0 == nil || ($0!.count == d && $0!.allSatisfy(\.isFinite)) }) else {
            throw VivoOmicsError.invalid("streamed ridge effect axes or values")
        }
        for i in 0..<n {
            if i % 2_048 == 0 { try Task.checkCancellation() }
            let weights = try memberships.row(i)
            var row = try scores.row(i)
            for c in 0..<k {
                guard let effect = effects[c * levels + batch[i]] else { continue }
                for j in 0..<d { row[j] -= weights[c] * effect[j] }
            }
            try scores.setRow(i, row)
        }
    }

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

    /// Categorical ridge with independently specified positive penalties. The
    /// unpenalized intercept uses the same positive Schur complement as fixed ridge.
    static func ridgeFit(masses: [Double], sums: [[Double]], penalties: [Double]) throws -> (intercept: [Double], effects: [[Double]], residual: Double) {
        let d = sums.first?.count ?? 0
        guard d > 0, masses.count >= 2, masses.count == sums.count, masses.count == penalties.count,
              masses.allSatisfy({ $0.isFinite && $0 >= 0 }), penalties.allSatisfy({ $0.isFinite && $0 > 0 }),
              sums.allSatisfy({ $0.count == d && $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration adaptive ridge inputs") }
        var denominator = 0.0
        for b in masses.indices { denominator += masses[b] * penalties[b] / (masses[b] + penalties[b]) }
        guard denominator > 1e-14 else { throw VivoOmicsError.invalid("integration cluster has no effective mass") }
        var intercept = [Double](repeating: 0, count: d)
        for b in masses.indices { for j in 0..<d { intercept[j] += sums[b][j] * penalties[b] / (masses[b] + penalties[b]) / denominator } }
        var effects = sums, residual = 0.0
        for b in masses.indices { for j in 0..<d {
            effects[b][j] = (sums[b][j] - masses[b] * intercept[j]) / (masses[b] + penalties[b])
            let error = masses[b] * intercept[j] + (masses[b] + penalties[b]) * effects[b][j] - sums[b][j]
            residual = max(residual, abs(error) / (1 + abs(sums[b][j])))
        } }
        for j in 0..<d {
            var error = 0.0, scale = 1.0
            for b in masses.indices { error += masses[b] * (intercept[j] + effects[b][j]) - sums[b][j]; scale += abs(sums[b][j]) }
            residual = max(residual, abs(error) / scale)
        }
        guard residual < 1e-10, intercept.allSatisfy(\.isFinite), effects.allSatisfy({ $0.allSatisfy(\.isFinite) }) else { throw VivoOmicsError.invalid("integration adaptive ridge normal-equation residual") }
        return (intercept, effects, residual)
    }

    /// Solve an additive intercept plus multiple categorical ridge model. The
    /// cross-factor masses are retained explicitly so a donor effect cannot be
    /// mistaken for a batch effect when the two factors are correlated. Ridge
    /// penalties make the bounded normal system positive definite; the solved
    /// coefficients and normalized normal-equation residual remain inspectable.
    static func additiveRidgeFit(totalMass: Double, totalSums: [Double], masses: [[Double]],
                                 sums: [[[Double]]], pairMasses: [[Double]],
                                 penalties: [[Double]]) throws -> (intercept: [Double], effects: [[[Double]]], residual: Double) {
        let factorCount = masses.count, d = totalSums.count
        guard factorCount >= 2, factorCount == sums.count, factorCount == penalties.count,
              pairMasses.count == factorCount * (factorCount - 1) / 2,
              totalMass.isFinite, totalMass > 0, totalSums.allSatisfy(\.isFinite) else {
            throw VivoOmicsError.invalid("joint integration ridge inputs")
        }
        let coefficientCounts = masses.map(\.count)
        let sumsShapeValid = sums.indices.allSatisfy { factor in
            sums[factor].count == coefficientCounts[factor] &&
            sums[factor].allSatisfy { $0.count == d && $0.allSatisfy(\.isFinite) }
        }
        let penaltiesShapeValid = penalties.indices.allSatisfy { factor in
            penalties[factor].count == coefficientCounts[factor] &&
            penalties[factor].allSatisfy { $0.isFinite && $0 > 0 }
        }
        guard coefficientCounts.allSatisfy({ $0 >= 2 }),
              masses.allSatisfy({ $0.allSatisfy({ $0.isFinite && $0 > 0 }) }),
              penalties.count == coefficientCounts.count,
              sumsShapeValid,
              penaltiesShapeValid else {
            throw VivoOmicsError.invalid("joint integration ridge factor axes")
        }
        let p = 1 + coefficientCounts.reduce(0, +)
        guard p <= 256 else { throw VivoOmicsError.limit("joint integration ridge coefficient budget") }
        var offsets = [Int](), next = 1
        for count in coefficientCounts { offsets.append(next); next += count }
        var normal = [[Double]](repeating: [Double](repeating: 0, count: p), count: p)
        normal[0][0] = totalMass
        var rhs = [[Double]](repeating: [Double](repeating: 0, count: d), count: p)
        rhs[0] = totalSums
        for factor in 0..<factorCount {
            for level in masses[factor].indices {
                let index = offsets[factor] + level, mass = masses[factor][level]
                normal[0][index] = mass; normal[index][0] = mass
                normal[index][index] = mass + penalties[factor][level]
                rhs[index] = sums[factor][level]
            }
        }
        var pair = 0
        for left in 0..<factorCount {
            for right in (left + 1)..<factorCount {
                let leftCount = coefficientCounts[left], rightCount = coefficientCounts[right]
                let table = pairMasses[pair]; pair += 1
                guard table.count == leftCount * rightCount,
                      table.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
                    throw VivoOmicsError.invalid("joint integration ridge cross-factor masses")
                }
                for i in 0..<leftCount { for j in 0..<rightCount {
                    let value = table[i * rightCount + j]
                    normal[offsets[left] + i][offsets[right] + j] = value
                    normal[offsets[right] + j][offsets[left] + i] = value
                } }
            }
        }
        // Cholesky is used only for this small (at most 256-square) positive
        // ridge system; it avoids a per-component QR allocation in the cell
        // loop while retaining an explicit residual check below.
        var lower = [[Double]](repeating: [Double](repeating: 0, count: p), count: p)
        let scale = max(1, totalMass)
        for i in 0..<p {
            for j in 0...i {
                var value = normal[i][j]
                if j > 0 { for q in 0..<j { value -= lower[i][q] * lower[j][q] } }
                if i == j {
                    guard value.isFinite, value > 1e-14 * scale else {
                        throw VivoOmicsError.invalid("joint integration ridge system is not positive definite")
                    }
                    lower[i][j] = sqrt(value)
                } else {
                    lower[i][j] = value / lower[j][j]
                }
            }
        }
        var solution = [[Double]](repeating: [Double](repeating: 0, count: d), count: p)
        for j in 0..<d {
            var forward = [Double](repeating: 0, count: p)
            for i in 0..<p {
                var value = rhs[i][j]
                if i > 0 { for q in 0..<i { value -= lower[i][q] * forward[q] } }
                forward[i] = value / lower[i][i]
            }
            for i in stride(from: p - 1, through: 0, by: -1) {
                var value = forward[i]
                if i + 1 < p { for q in (i + 1)..<p { value -= lower[q][i] * solution[q][j] } }
                solution[i][j] = value / lower[i][i]
            }
        }
        var residual = 0.0
        for i in 0..<p {
            for j in 0..<d {
                var value = 0.0
                for q in 0..<p { value += normal[i][q] * solution[q][j] }
                let error = abs(value - rhs[i][j]) / (1 + abs(rhs[i][j]))
                guard error.isFinite else { throw VivoOmicsError.invalid("joint integration ridge residual overflow") }
                residual = max(residual, error)
            }
        }
        guard residual < 1e-9, solution.allSatisfy({ $0.allSatisfy(\.isFinite) }) else {
            throw VivoOmicsError.invalid("joint integration ridge normal-equation residual")
        }
        var effects: [[[Double]]] = []
        for factor in 0..<factorCount {
            effects.append((0..<coefficientCounts[factor]).map { solution[offsets[factor] + $0] })
        }
        return (solution[0], effects, residual)
    }

    static func run(_ reduction: VivoSingleCellReductionResult, samples: [VivoOmicsSample], options o: VivoSingleCellIntegrationOptions) throws -> VivoSingleCellIntegrationResult {
        let n = reduction.scores.count, d = reduction.scores.first?.count ?? 0
        try validateAxes(rows: n, columns: d, options: o)
        guard n == reduction.cells.count else { throw VivoOmicsError.invalid("integration PCA axes") }
        let x = try VivoIntegrationMatrix(rows: n, columns: d)
        for i in 0..<n { try x.setRow(i, reduction.scores[i]) }
        let result = try run(cells: reduction.cells, x: x, samples: samples, options: o) { rows, columns in
            try VivoIntegrationMatrix(rows: rows, columns: columns)
        }
        return try result.materialize(cells: reduction.cells, options: o)
    }
    static func validateAxes(rows n: Int, columns d: Int, options o: VivoSingleCellIntegrationOptions) throws {
        try o.validate()
        guard (o.clusters...VivoPCAStorageLimits.maximumRows).contains(n), (1...64).contains(d) else { throw VivoOmicsError.invalid("integration PCA axes") }
        // An explicit work index, not an instruction count. Check before allocation.
        var work = n
        for factor in [o.clusters, d, o.maximumIterations + 10] {
            let product = work.multipliedReportingOverflow(by: factor)
            guard !product.overflow, product.partialValue <= o.maximumWork else { throw VivoOmicsError.limit("integration work budget") }
            work = product.partialValue
        }
    }

    /// Opt-in additive correction for both donor and batch. The historical
    /// single-factor path below is intentionally left byte-for-byte stable;
    /// this route has its own explicit factor axes and normal-equation witness.
    static func runJoint(cells: [VivoOmicsCellIdentity], x: VivoIntegrationMatrix, samples: [VivoOmicsSample],
                         options o: VivoSingleCellIntegrationOptions, matrix: (Int, Int) throws -> VivoIntegrationMatrix) throws -> VivoIntegrationSolution {
        guard let selected = o.covariates, selected.count == 2 else {
            throw VivoOmicsError.invalid("joint integration requires donor and batch covariates")
        }
        let n = x.rows, d = x.columns, k = o.clusters
        try validateAxes(rows: n, columns: d, options: o)
        guard cells.count == n, Set(samples.map(\.id)).count == samples.count else {
            throw VivoOmicsError.invalid("joint integration PCA identities or samples")
        }
        let lookup = Dictionary(uniqueKeysWithValues: samples.map { ($0.id, $0) })
        let names: [[String]] = try selected.map { covariate in
            try cells.map { cell in
                guard let sample = lookup[cell.sampleID],
                      let value = covariate == .donor ? sample.donorID : sample.batchID,
                      !value.isEmpty, value != "unreported" else {
                    throw VivoOmicsError.invalid("joint integration requires known selected covariates")
                }
                return value
            }
        }
        let factorLevels = names.map { Set($0).sorted() }
        guard factorLevels.allSatisfy({ (2...128).contains($0.count) }) else {
            throw VivoOmicsError.invalid("joint integration requires 2...128 levels per covariate")
        }
        let factorIDs = factorLevels.indices.map { factor in
            let ids = Dictionary(uniqueKeysWithValues: factorLevels[factor].enumerated().map { ($0.element, $0.offset) })
            return names[factor].map { ids[$0]! }
        }
        let levelCounts = factorLevels.map(\.count), totalLevels = levelCounts.reduce(0, +)
        var work = n
        for factor in [k, d, o.maximumIterations + 10, max(1, totalLevels)] {
            let product = work.multipliedReportingOverflow(by: factor)
            guard !product.overflow, product.partialValue <= o.maximumWork else {
                throw VivoOmicsError.limit("joint integration work budget")
            }
            work = product.partialValue
        }
        // Every factor must be connected through shared condition labels. A
        // pairwise factor graph additionally rejects perfect donor/batch
        // confounding, for which additive effects have no data identification.
        for factor in factorIDs.indices {
            var conditionLevels: [String: Set<Int>] = [:]
            for i in cells.indices {
                conditionLevels[lookup[cells[i].sampleID]!.condition, default: []].insert(factorIDs[factor][i])
            }
            var reachable: Set<Int> = [0]
            for _ in 0..<levelCounts[factor] {
                let before = reachable.count
                for levels in conditionLevels.values where !reachable.isDisjoint(with: levels) { reachable.formUnion(levels) }
                if reachable.count == before { break }
            }
            guard reachable.count == levelCounts[factor] else {
                throw VivoOmicsError.invalid("joint integration covariate is confounded with condition")
            }
        }
        for left in factorIDs.indices {
            for right in (left + 1)..<factorIDs.count {
                var adjacency = [[Int]](repeating: [], count: levelCounts[left] + levelCounts[right])
                for i in cells.indices {
                    let a = factorIDs[left][i], b = factorIDs[right][i] + levelCounts[left]
                    adjacency[a].append(b); adjacency[b].append(a)
                }
                var reachable: Set<Int> = [0]
                for _ in adjacency.indices {
                    let before = reachable.count
                    for node in Array(reachable).sorted() { reachable.formUnion(adjacency[node]) }
                    if reachable.count == before { break }
                }
                guard reachable.count == adjacency.count else {
                    throw VivoOmicsError.invalid("joint integration covariates are perfectly confounded")
                }
            }
        }
        try VivoIntegrationProtection.check(o.protectedSampleGroups, cells: cells, samples: lookup,
            covariates: selected, levelIDs: factorLevels.indices.map { factor in
                Dictionary(uniqueKeysWithValues: factorLevels[factor].enumerated().map { ($0.element, $0.offset) })
            })
        let factorSizes = factorIDs.map { ids in
            var result = [Double](repeating: 0, count: (ids.max() ?? -1) + 1)
            for id in ids { result[id] += 1 }
            return result
        }
        guard factorSizes.count == factorLevels.count else { throw VivoOmicsError.invalid("joint integration factor sizes") }
        var state = o.seed
        func uniform() -> Double {
            state &+= 0x9E3779B97F4A7C15; var z = state
            z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9; z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
            return Double((z ^ (z >> 31)) >> 11) / 9007199254740992.0
        }
        func normalized(_ rows: [[Double]]) -> [[Double]] {
            rows.map { row in let norm = sqrt(row.reduce(0) { $0 + $1 * $1 }); return norm > 0 ? row.map { $0 / norm } : row }
        }
        func distance(_ a: [Double], _ b: [Double]) -> Double { zip(a, b).reduce(0) { $0 + ($1.0 - $1.1) * ($1.0 - $1.1) } }
        let unit = try matrix(n, d); try unit.copy(from: x, normalize: true)
        var centers = [try unit.row(min(n - 1, Int(uniform() * Double(n))))]
        var nearest = [Double](repeating: .infinity, count: n)
        while centers.count < k {
            for i in 0..<n { nearest[i] = min(nearest[i], distance(try unit.row(i), centers.last!)) }
            let total = nearest.reduce(0, +)
            guard total > 1e-14 else { throw VivoOmicsError.invalid("joint integration clusters exceed distinct PCA directions") }
            var threshold = uniform() * total, chosen = n - 1
            for i in 0..<n { threshold -= nearest[i]; if threshold <= 0 { chosen = i; break } }
            centers.append(try unit.row(chosen))
        }
        for _ in 0..<10 {
            var sums = [[Double]](repeating: [Double](repeating: 0, count: d), count: k), counts = [Int](repeating: 0, count: k)
            for i in 0..<n {
                let row = try unit.row(i), c = (0..<k).min { distance(row, centers[$0]) < distance(row, centers[$1]) }!
                counts[c] += 1; for j in 0..<d { sums[c][j] += row[j] }
            }
            for c in 0..<k where counts[c] > 0 { centers[c] = sums[c].map { $0 / Double(counts[c]) } }
        }
        centers = normalized(centers)
        let scores = try matrix(n, d), r = try matrix(n, k), distances = try matrix(n, k), assignmentScores = try matrix(n, d)
        try scores.copy(from: x)
        var observed = factorLevels.indices.map { factor in [Double](repeating: 0, count: k * levelCounts[factor]) }
        var masses = [Double](repeating: 0, count: k)
        func probabilities(_ i: Int, _ row: [Double], _ diversity: Bool) -> [Double] {
            var logits = [Double](repeating: 0, count: k)
            for c in 0..<k {
                logits[c] = -row[c] / o.temperature
                if diversity {
                    for factor in factorIDs.indices {
                        let level = factorIDs[factor][i]
                        let expected = max(0, masses[c]) * factorSizes[factor][level] / Double(n)
                        let value = max(0, observed[factor][c * levelCounts[factor] + level])
                        logits[c] += o.diversity * log((2 * expected + 1) / (value + expected + 1))
                    }
                }
            }
            let peak = logits.max()!, weights = logits.map { exp($0 - peak) }, total = weights.reduce(0, +)
            return weights.map { $0 / total }
        }
        func addRow(_ i: Int, _ row: [Double], _ sign: Double) {
            for c in 0..<k {
                let value = sign * row[c]; masses[c] += value
                for factor in factorIDs.indices {
                    observed[factor][c * levelCounts[factor] + factorIDs[factor][i]] += value
                }
            }
        }
        func objective() throws -> Double {
            var result = 0.0
            for i in 0..<n {
                let row = try r.row(i), distance = try distances.row(i)
                for c in 0..<k where row[c] > 0 { result += row[c] * distance[c] + o.temperature * row[c] * log(row[c]) }
            }
            for factor in factorIDs.indices {
                for c in 0..<k { for level in 0..<levelCounts[factor] {
                    let value = max(0, observed[factor][c * levelCounts[factor] + level])
                    let expected = max(0, masses[c]) * factorSizes[factor][level] / Double(n)
                    result += o.temperature * o.diversity * value * log((value + expected + 1) / (2 * expected + 1))
                } }
            }
            return result * 2000 / Double(n)
        }
        var objectives: [Double] = [], improvements: [Double] = [], maxResidual = 0.0, stopping = "iteration-limit"
        var assignmentCenters: [[Double]] = []
        var factorPenalties: [[[Double]]]? = o.ridgeScaling == nil ? nil : factorLevels.map { levels in
            Array(repeating: [Double](repeating: 0, count: levels.count), count: k)
        }
        var factorEffects = factorLevels.map { levels in Array<[Double]?>(repeating: nil, count: k * levels.count) }
        for iteration in 0..<o.maximumIterations {
            try Task.checkCancellation(); try assignmentScores.copy(from: scores, normalize: true); assignmentCenters = centers
            observed = factorLevels.indices.map { factor in [Double](repeating: 0, count: k * levelCounts[factor]) }
            masses = [Double](repeating: 0, count: k)
            for i in 0..<n {
                let current = try assignmentScores.row(i)
                try distances.setRow(i, (0..<k).map { c in 2 * (1 - zip(current, centers[c]).reduce(0) { $0 + $1.0 * $1.1 }) })
                try r.setRow(i, probabilities(i, distances.row(i), false)); addRow(i, try r.row(i), 1)
            }
            if iteration == 0 { objectives.append(try objective()) }
            for _ in 0..<4 {
                var order = Array(0..<n)
                for i in stride(from: n - 1, through: 1, by: -1) { order.swapAt(i, min(i, Int(uniform() * Double(i + 1)))) }
                let blockSize = max(1, n / 20)
                for block in 0..<20 {
                    let start = block * blockSize, end = block == 19 ? n : min(n, start + blockSize)
                    if start >= n { break }
                    for t in start..<end { let i = order[t]; addRow(i, try r.row(i), -1) }
                    for t in start..<end { let i = order[t]; try r.setRow(i, probabilities(i, try distances.row(i), true)) }
                    for t in start..<end { let i = order[t]; addRow(i, try r.row(i), 1) }
                }
            }
            let value = try objective(); guard value.isFinite else { throw VivoOmicsError.invalid("joint integration objective nonfinite") }
            let previous = objectives.last!, improvement = (previous - value) / max(abs(previous), 1e-12)
            objectives.append(value); improvements.append(improvement)
            try scores.copy(from: x)
            factorEffects = factorLevels.map { levels in Array<[Double]?>(repeating: nil, count: k * levels.count) }
            for c in 0..<k {
                let active = factorIDs.indices.map { factor in
                    (0..<levelCounts[factor]).filter { observed[factor][c * levelCounts[factor] + $0] / factorSizes[factor][$0] > 1e-5 }
                }
                guard active.allSatisfy({ $0.count >= 2 }) else { continue }
                let activeSets = active.map { Set($0) }
                var factorMasses = active.map { [Double](repeating: 0, count: $0.count) }
                var factorSums = active.map { levels in levels.map { _ in [Double](repeating: 0, count: d) } }
                var pairMasses = [[Double]]()
                for left in factorIDs.indices { for right in (left + 1)..<factorIDs.count {
                    pairMasses.append([Double](repeating: 0, count: active[left].count * active[right].count))
                } }
                var totalMass = 0.0, totalSums = [Double](repeating: 0, count: d), pairIndex = 0
                for i in 0..<n where factorIDs.indices.allSatisfy({ activeSets[$0].contains(factorIDs[$0][i]) }) {
                    let weight = try r.value(i, c); guard weight > 0 else { continue }
                    let row = try x.row(i); totalMass += weight
                    for j in 0..<d { totalSums[j] += weight * row[j] }
                    for factor in factorIDs.indices {
                        let slot = active[factor].firstIndex(of: factorIDs[factor][i])!
                        factorMasses[factor][slot] += weight
                        for j in 0..<d { factorSums[factor][slot][j] += weight * row[j] }
                    }
                    pairIndex = 0
                    for left in factorIDs.indices { for right in (left + 1)..<factorIDs.count {
                        let a = active[left].firstIndex(of: factorIDs[left][i])!, b = active[right].firstIndex(of: factorIDs[right][i])!
                        pairMasses[pairIndex][a * active[right].count + b] += weight; pairIndex += 1
                    } }
                }
                guard totalMass > 1e-14 else { continue }
                // A marginally active level can lose all support after the
                // complete factor-product filter. Leave that cluster
                // uncorrected instead of constructing a singular design.
                guard factorMasses.allSatisfy({ $0.allSatisfy({ $0 > 1e-14 }) }) else { continue }
                var penalties = active.map { levels in [Double](repeating: o.ridge, count: levels.count) }
                if o.ridgeScaling == .expectedClusterBatchMass {
                    for factor in factorIDs.indices { for slot in active[factor].indices {
                        let level = active[factor][slot]
                        penalties[factor][slot] = o.ridge * max(0, masses[c]) * factorSizes[factor][level] / Double(n)
                    } }
                }
                if factorPenalties != nil {
                    for factor in factorIDs.indices { for slot in active[factor].indices {
                        factorPenalties![factor][c][active[factor][slot]] = penalties[factor][slot]
                    } }
                }
                let fit = try additiveRidgeFit(totalMass: totalMass, totalSums: totalSums, masses: factorMasses,
                    sums: factorSums, pairMasses: pairMasses, penalties: penalties)
                maxResidual = max(maxResidual, fit.residual); centers[c] = fit.intercept
                for factor in factorIDs.indices { for slot in active[factor].indices {
                    factorEffects[factor][c * levelCounts[factor] + active[factor][slot]] = fit.effects[factor][slot]
                } }
            }
            for i in 0..<n {
                var row = try scores.row(i)
                for c in 0..<k {
                    let weight = try r.value(i, c); if weight == 0 { continue }
                    for factor in factorIDs.indices {
                        if let effect = factorEffects[factor][c * levelCounts[factor] + factorIDs[factor][i]] {
                            for j in 0..<d { row[j] -= weight * effect[j] }
                        }
                    }
                }
                try scores.setRow(i, row)
            }
            centers = normalized(centers)
            if improvement < 0 { stopping = "objective-increase"; break }
            if improvement < o.relativeTolerance { stopping = "relative-objective-tolerance"; break }
        }
        return .init(levels: factorLevels[0], cellLevels: factorIDs[0], scores: scores, memberships: r,
            assignmentScores: assignmentScores, assignmentCenters: assignmentCenters, objectives: objectives,
            relativeImprovements: improvements, stoppingReason: stopping, maximumRidgeResidual: maxResidual,
            ridgePenalties: factorPenalties?[0], factorLevels: factorLevels, cellFactorLevels: factorIDs,
            factorRidgePenalties: factorPenalties)
    }
    static func run(cells: [VivoOmicsCellIdentity], x: VivoIntegrationMatrix, samples: [VivoOmicsSample],
                    options o: VivoSingleCellIntegrationOptions, matrix: (Int, Int) throws -> VivoIntegrationMatrix) throws -> VivoIntegrationSolution {
        let n = x.rows, d = x.columns, k = o.clusters
        try validateAxes(rows: n, columns: d, options: o)
        guard cells.count == n else { throw VivoOmicsError.invalid("integration PCA identities") }
        if o.covariates != nil {
            return try runJoint(cells: cells, x: x, samples: samples, options: o, matrix: matrix)
        }
        guard Set(samples.map(\.id)).count == samples.count else { throw VivoOmicsError.invalid("duplicate integration sample identities") }
        let lookup = Dictionary(uniqueKeysWithValues: samples.map { ($0.id,$0) })
        let names = try cells.map { cell -> String in
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
        for (i,cell) in cells.enumerated() { conditionLevels[lookup[cell.sampleID]!.condition,default: []].insert(batch[i]) }
        var reachable: Set<Int> = [0]
        for _ in 0..<bCount {
            let before = reachable.count
            for group in conditionLevels.values where !reachable.isDisjoint(with: group) { reachable.formUnion(group) }
            if reachable.count == before { break }
        }
        guard reachable.count == bCount else { throw VivoOmicsError.invalid("integration covariate is confounded with condition; no connected shared-condition design") }
        try VivoIntegrationProtection.check(o.protectedSampleGroups, cells: cells,
            samples: lookup, covariate: o.covariate, levelIDs: ids)
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
        let unit = try matrix(n, d)
        try unit.copy(from: x, normalize: true)
        // Independent seeded kmeans++ and Lloyd initialization. Reference uses a
        // different random generator/initializer; coordinate identity is not claimed.
        var centers = [try unit.row(min(n-1,Int(uniform()*Double(n))))]
        var nearest = [Double](repeating: .infinity,count: n)
        while centers.count < k {
            for i in 0..<n { nearest[i] = min(nearest[i],distance(try unit.row(i),centers.last!)) }
            let total = nearest.reduce(0,+)
            guard total > 1e-14 else { throw VivoOmicsError.invalid("integration cluster count exceeds distinct PCA directions") }
            var threshold = uniform()*total, selected = n-1
            for i in 0..<n { threshold -= nearest[i]; if threshold <= 0 { selected = i; break } }
            centers.append(try unit.row(selected))
        }
        for _ in 0..<10 {
            var sums = [[Double]](repeating: [Double](repeating: 0,count: d),count: k), counts = [Int](repeating: 0,count: k)
            for i in 0..<n {
                let row = try unit.row(i)
                let c = (0..<k).min { distance(row,centers[$0]) < distance(row,centers[$1]) }!
                counts[c] += 1; for j in 0..<d { sums[c][j] += row[j] }
            }
            for c in 0..<k where counts[c] > 0 { centers[c] = sums[c].map { $0/Double(counts[c]) } }
        }
        centers = normalized(centers)
        let scores = try matrix(n, d), r = try matrix(n, k), distances = try matrix(n, k)
        let assignmentScores = try matrix(n, d)
        try scores.copy(from: x)
        var observed = [Double](repeating: 0,count: k*bCount), masses = [Double](repeating: 0,count: k)
        func probabilities(_ i: Int, _ row: [Double], _ diversity: Bool) -> [Double] {
            var logits = [Double](repeating: 0,count: k)
            for c in 0..<k {
                let expected = max(0,masses[c])*sizes[batch[i]]/Double(n), obs = max(0,observed[c*bCount+batch[i]])
                logits[c] = -row[c]/o.temperature + (diversity ? o.diversity*log((2*expected+1)/(obs+expected+1)) : 0)
            }
            let peak = logits.max()!, weights = logits.map { exp($0-peak) }, total = weights.reduce(0,+)
            return weights.map { $0 / total }
        }
        func setProbabilities(_ i: Int,_ diversity: Bool) throws {
            try r.setRow(i, probabilities(i, distances.row(i), diversity))
        }
        func addRow(_ i: Int, _ row: [Double], _ sign: Double) {
            for c in 0..<k { let value = sign*row[c]; masses[c] += value; observed[c*bCount+batch[i]] += value }
        }
        func add(_ i: Int,_ sign: Double) throws { try addRow(i, r.row(i), sign) }
        func objective() throws -> Double {
            var result = 0.0
            for i in 0..<n {
                let row = try r.row(i), distance = try distances.row(i)
                for c in 0..<k where row[c] > 0 { result += row[c]*distance[c]+o.temperature*row[c]*log(row[c]) }
            }
            for c in 0..<k { for b in 0..<bCount {
                let obs = max(0,observed[c*bCount+b]), expected = max(0,masses[c])*sizes[b]/Double(n)
                result += o.temperature*o.diversity*obs*log((obs+expected+1)/(2*expected+1))
            } }
            return result*2000/Double(n)
        }
        var objectives: [Double] = [], improvements: [Double] = [], maxResidual = 0.0, stopping = "iteration-limit"
        var assignmentCenters: [[Double]] = []
        var ridgePenalties: [[Double]]? = o.ridgeScaling == nil ? nil : [[Double]](repeating: [Double](repeating: 0, count: bCount), count: k)
        for iteration in 0..<o.maximumIterations {
            try Task.checkCancellation()
            try assignmentScores.copy(from: scores, normalize: true)
            assignmentCenters = centers
            observed = observed.map { _ in 0 }; masses = masses.map { _ in 0 }
            for i in 0..<n {
                let current = try assignmentScores.row(i)
                try distances.setRow(i, (0..<k).map { c in 2*(1-zip(current,centers[c]).reduce(0) { $0+$1.0*$1.1 }) })
                try setProbabilities(i,false); try add(i,1)
            }
            if iteration == 0 { objectives.append(try objective()) }
            // Four block-coordinate sweeps; block fraction .05, as the pinned reference.
            for _ in 0..<4 {
                var order = Array(0..<n)
                for i in stride(from: n-1,through: 1,by: -1) { order.swapAt(i,min(i,Int(uniform()*Double(i+1)))) }
                let blockSize = max(1,n/20)
                for block in 0..<20 {
                    let start = block*blockSize, end = block == 19 ? n : min(n,start+blockSize)
                    if start >= n { break }
                    if r.benefitsFromBatchedAccess || distances.benefitsFromBatchedAccess {
                        // Keep all three whole-block phases separate. Tiling the
                        // phases together would change the diversity statistics.
                        let tileRows = VivoIntegrationMatrix.maximumBatchRows
                        for tileStart in stride(from: start, to: end, by: tileRows) {
                            let indices = Array(order[tileStart..<min(end, tileStart + tileRows)])
                            let rows = try r.gatherRows(indices)
                            for slot in indices.indices { addRow(indices[slot], rows[slot], -1) }
                        }
                        for tileStart in stride(from: start, to: end, by: tileRows) {
                            let indices = Array(order[tileStart..<min(end, tileStart + tileRows)])
                            var rows = try distances.gatherRows(indices)
                            for slot in indices.indices { rows[slot] = probabilities(indices[slot], rows[slot], true) }
                            try r.scatterRows(indices, rows: rows)
                        }
                        for tileStart in stride(from: start, to: end, by: tileRows) {
                            let indices = Array(order[tileStart..<min(end, tileStart + tileRows)])
                            let rows = try r.gatherRows(indices)
                            for slot in indices.indices { addRow(indices[slot], rows[slot], 1) }
                        }
                    } else {
                        for t in start..<end { try add(order[t],-1) }
                        for t in start..<end { try setProbabilities(order[t],true) }
                        for t in start..<end { try add(order[t],1) }
                    }
                }
            }
            let value = try objective()
            guard value.isFinite else { throw VivoOmicsError.invalid("integration objective nonfinite") }
            let previous = objectives.last!, improvement = (previous-value)/max(abs(previous),1e-12)
            objectives.append(value); improvements.append(improvement)
            // Every correction is fitted to and subtracted from original PCA.
            try scores.copy(from: x)
            let streamRidge = x.benefitsFromBatchedAccess || r.benefitsFromBatchedAccess || scores.benefitsFromBatchedAccess
            let activeLevels = (0..<k).map { c in (0..<bCount).filter { observed[c*bCount+$0]/sizes[$0] > 1e-5 } }
            let streamedSums = try streamRidge ? streamedRidgeSums(x: x, memberships: r, batch: batch,
                levels: bCount, activeClusters: (0..<k).filter { activeLevels[$0].count >= 2 }) : nil
            var streamedEffects = streamRidge ? [[Double]?](repeating: nil, count: k * bCount) : []
            for c in 0..<k {
                if o.ridgeScaling == .expectedClusterBatchMass {
                    ridgePenalties![c] = sizes.map { size in
                        let expected = max(0, masses[c]) * size / Double(n)
                        return o.ridge * expected
                    }
                }
                let active = activeLevels[c]
                if active.count < 2 { continue }
                var sums = [[Double]](repeating: [Double](repeating: 0,count: d),count: bCount)
                if let streamedSums {
                    for b in 0..<bCount {
                        let offset = (c * bCount + b) * d
                        sums[b] = Array(streamedSums[offset..<(offset + d)])
                    }
                } else {
                    for i in 0..<n {
                        let weight = try r.value(i, c), row = try x.row(i)
                        for j in 0..<d { sums[batch[i]][j] += weight*row[j] }
                    }
                }
                let fit: (intercept: [Double], effects: [[Double]], residual: Double)
                if let ridgePenalties {
                    fit = try ridgeFit(masses: active.map { max(0, observed[c*bCount+$0]) }, sums: active.map { sums[$0] }, penalties: active.map { ridgePenalties[c][$0] })
                } else {
                    fit = try ridgeFit(masses: active.map { max(0,observed[c*bCount+$0]) },sums: active.map { sums[$0] },ridge: o.ridge)
                }
                maxResidual = max(maxResidual,fit.residual); centers[c] = fit.intercept
                if streamRidge {
                    for (slot, b) in active.enumerated() { streamedEffects[c * bCount + b] = fit.effects[slot] }
                } else {
                    for (slot,b) in active.enumerated() { for i in 0..<n where batch[i] == b {
                        let weight = try r.value(i, c); var row = try scores.row(i)
                        for j in 0..<d { row[j] -= weight*fit.effects[slot][j] }
                        try scores.setRow(i, row)
                    } }
                }
            }
            if streamRidge { try applyStreamedRidgeEffects(scores: scores, memberships: r,
                batch: batch, levels: bCount, effects: streamedEffects) }
            centers = normalized(centers)
            // Each matrix row write rejects nonfinite corrected values.
            if improvement < 0 { stopping = "objective-increase"; break }
            if improvement < o.relativeTolerance { stopping = "relative-objective-tolerance"; break }
        }
        return .init(levels: levels, cellLevels: batch, scores: scores, memberships: r,
            assignmentScores: assignmentScores, assignmentCenters: assignmentCenters, objectives: objectives,
            relativeImprovements: improvements, stoppingReason: stopping, maximumRidgeResidual: maxResidual,
            ridgePenalties: ridgePenalties, factorLevels: nil, cellFactorLevels: nil, factorRidgePenalties: nil)
    }
}

struct VivoIntegrationSolution {
    static func method(options: VivoSingleCellIntegrationOptions) -> String {
        if options.covariates != nil {
            return options.ridgeScaling == nil ? "diversity-soft-clustering-multi-categorical-ridge-Double-v1" : "diversity-soft-clustering-multi-categorical-expected-mass-ridge-Double-v1"
        }
        return options.ridgeScaling == nil ? "diversity-soft-clustering-categorical-ridge-Double-v1" : "diversity-soft-clustering-expected-mass-ridge-Double-v1"
    }
    static func qualification(options: VivoSingleCellIntegrationOptions) -> String {
        let scope = options.covariates == nil ? "single-covariate" : "additive multi-covariate"
        return "Transductive \(scope) PCA correction with Harmony2 objective and intercept centers; independent initialization. Original counts/PCA preserved. Stopping diagnostics are not biological preservation or prospective prediction evidence."
    }
    let levels: [String]
    let cellLevels: [Int]
    let scores: VivoIntegrationMatrix
    let memberships: VivoIntegrationMatrix
    let assignmentScores: VivoIntegrationMatrix
    let assignmentCenters: [[Double]]
    let objectives: [Double]
    let relativeImprovements: [Double]
    let stoppingReason: String
    let maximumRidgeResidual: Double
    let ridgePenalties: [[Double]]?
    let factorLevels: [[String]]?
    let cellFactorLevels: [[Int]]?
    let factorRidgePenalties: [[[Double]]]?
    func materialize(cells: [VivoOmicsCellIdentity], options: VivoSingleCellIntegrationOptions) throws -> VivoSingleCellIntegrationResult {
        try .init(method: Self.method(options: options), options: options, cells: cells, levels: levels, cellLevels: cellLevels,
            scores: scores.materialize(), memberships: memberships.materialize(), assignmentScores: assignmentScores.materialize(),
            assignmentCenters: assignmentCenters, objectives: objectives, relativeImprovements: relativeImprovements,
            stoppingReason: stoppingReason, maximumRidgeResidual: maximumRidgeResidual, qualification: Self.qualification(options: options),
            ridgePenalties: ridgePenalties, factorLevels: factorLevels, cellFactorLevels: cellFactorLevels,
            factorRidgePenalties: factorRidgePenalties)
    }
}
