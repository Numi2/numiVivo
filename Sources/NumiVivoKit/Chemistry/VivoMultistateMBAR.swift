import Foundation

public struct VivoMBARReducedPotentialSample: Codable, Sendable, Equatable {
    public var identifier: String
    public var originStateIndex: Int
    /// Dimensionless reduced potentials u_k(x)=beta*U_k(x) for every state.
    public var reducedPotentials: [Double]

    public init(identifier: String, originStateIndex: Int, reducedPotentials: [Double]) {
        self.identifier = identifier
        self.originStateIndex = originStateIndex
        self.reducedPotentials = reducedPotentials
    }
}

public struct VivoMBARConfiguration: Codable, Sendable, Equatable {
    public var residualTolerance: Double
    public var maximumIterations: Int
    /// Minimum of mean Metropolis acceptance in both directions after correcting
    /// the cross-state energy difference by the converged MBAR free-energy offset.
    public var minimumBidirectionalAcceptanceOverlap: Double
    public var minimumTargetEffectiveSamples: Double
    public var bootstrapReplicates: Int
    public var bootstrapSeed: UInt64
    public var maximumWorkElements: Int

    public init(residualTolerance: Double = 1e-10,
                maximumIterations: Int = 20_000,
                minimumBidirectionalAcceptanceOverlap: Double = 0.01,
                minimumTargetEffectiveSamples: Double = 20,
                bootstrapReplicates: Int = 128,
                bootstrapSeed: UInt64 = 0x4D424152,
                maximumWorkElements: Int = 100_000_000) {
        self.residualTolerance = residualTolerance
        self.maximumIterations = maximumIterations
        self.minimumBidirectionalAcceptanceOverlap = minimumBidirectionalAcceptanceOverlap
        self.minimumTargetEffectiveSamples = minimumTargetEffectiveSamples
        self.bootstrapReplicates = bootstrapReplicates
        self.bootstrapSeed = bootstrapSeed
        self.maximumWorkElements = maximumWorkElements
    }

    public func validate() throws {
        guard residualTolerance.isFinite, residualTolerance > 0, residualTolerance <= 1e-3,
              (1...1_000_000).contains(maximumIterations),
              minimumBidirectionalAcceptanceOverlap.isFinite,
              minimumBidirectionalAcceptanceOverlap > 0,
              minimumBidirectionalAcceptanceOverlap <= 1,
              minimumTargetEffectiveSamples.isFinite,
              minimumTargetEffectiveSamples >= 2,
              (0...4096).contains(bootstrapReplicates),
              maximumWorkElements > 0 else {
            throw VivoChemistryError.invalid("MBAR numerical configuration")
        }
    }
}

public struct VivoMBARPairwiseOverlap: Codable, Sendable, Equatable {
    public let stateA: Int
    public let stateB: Int
    /// Normalized importance-weight ESS fraction. Useful for diagnosing weight
    /// concentration, but not sufficient by itself to establish phase-space overlap.
    public let reweightingESSFractionAToB: Double
    public let reweightingESSFractionBToA: Double
    /// Mean free-energy-corrected Metropolis acceptance from configurations drawn
    /// in each originating state. A constant energy offset therefore preserves
    /// overlap, while mutually inaccessible wells have near-zero acceptance.
    public let meanAcceptanceAToB: Double
    public let meanAcceptanceBToA: Double
    public let symmetricAcceptanceOverlap: Double
    public let connected: Bool
}

public struct VivoMBARStateEstimate: Codable, Sendable, Equatable {
    public let stateIndex: Int
    /// Dimensionless free energy relative to state zero.
    public let relativeFreeEnergy: Double
    public let bootstrapStandardDeviation: Double?
    /// Effective sample count under the final mixture-reweighting weights.
    public let effectiveSamples: Double
}

public struct VivoMBARResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/multistate-mbar-result/v1"
    public let schema: String
    public let stateCount: Int
    public let sampleCount: Int
    public let originCounts: [Int]
    public let estimates: [VivoMBARStateEstimate]
    public let pairwiseOverlap: [VivoMBARPairwiseOverlap]
    public let iterations: Int
    public let maximumResidual: Double
    public let converged: Bool
    public let overlapGraphConnected: Bool
    public let issues: [String]
    public let interpretation: String
}

/// General unbinned MBAR over dimensionless reduced potentials. The caller owns
/// sample generation and decorrelation. This solver never interprets correlated
/// trajectory frames as independent evidence and does not infer missing states.
public enum VivoMultistateMBAR {
    public static let interpretation = "Log-domain multistate Bennett acceptance ratio over explicitly supplied decorrelated reduced-potential samples. Acceptance requires a connected free-energy-corrected bidirectional acceptance-overlap graph and sufficient target-state effective sample counts. Stratified bootstrap resamples within originating states and quantifies finite independent-sample dispersion only."

    public static func solve(samples: [VivoMBARReducedPotentialSample],
                             stateCount: Int,
                             configuration cfg: VivoMBARConfiguration = .init()) throws -> VivoMBARResult {
        try cfg.validate()
        guard stateCount >= 1, stateCount <= 4096,
              !samples.isEmpty, samples.count <= 10_000_000,
              Set(samples.map(\.identifier)).count == samples.count else {
            throw VivoChemistryError.invalid("MBAR state/sample identity or capacity")
        }
        let (matrixWork, matrixOverflow) = samples.count.multipliedReportingOverflow(by: stateCount)
        guard !matrixOverflow, matrixWork <= cfg.maximumWorkElements else {
            throw VivoChemistryError.resourceLimit("MBAR reduced-potential matrix exceeds work capacity")
        }
        var counts = [Int](repeating: 0, count: stateCount)
        for sample in samples {
            guard sample.originStateIndex >= 0, sample.originStateIndex < stateCount,
                  sample.reducedPotentials.count == stateCount,
                  sample.reducedPotentials.allSatisfy({ $0.isFinite && abs($0) <= 1e12 }) else {
                throw VivoChemistryError.invalid("MBAR sample dimensions, origin or reduced potential")
            }
            counts[sample.originStateIndex] += 1
        }
        guard counts.allSatisfy({ $0 >= 2 }) else {
            throw VivoChemistryError.invalid("MBAR requires at least two decorrelated samples from every state")
        }
        let core = try solveCore(samples: samples, stateCount: stateCount, counts: counts,
                                 tolerance: cfg.residualTolerance, maximumIterations: cfg.maximumIterations)
        let effective = targetEffectiveSamples(samples: samples, counts: counts, freeEnergies: core.freeEnergies)
        let overlap = pairwiseOverlap(samples: samples, counts: counts, freeEnergies: core.freeEnergies,
                                      threshold: cfg.minimumBidirectionalAcceptanceOverlap)
        let connected = graphConnected(stateCount: stateCount, overlap: overlap)
        var issues: [String] = []
        if !core.converged { issues.append("MBAR fixed-point residual did not converge") }
        if !connected { issues.append("pairwise bidirectional acceptance-overlap graph is disconnected") }
        for (index, value) in effective.enumerated() where value < cfg.minimumTargetEffectiveSamples {
            issues.append("state \(index) effective sample count \(value) is below \(cfg.minimumTargetEffectiveSamples)")
        }
        var bootstrapSD: [Double?] = [Double?](repeating: nil, count: stateCount)
        if cfg.bootstrapReplicates > 1, core.converged, connected {
            let (bootstrapWork, overflow) = matrixWork.multipliedReportingOverflow(by: cfg.bootstrapReplicates)
            guard !overflow, bootstrapWork <= cfg.maximumWorkElements else {
                throw VivoChemistryError.resourceLimit("MBAR bootstrap exceeds work capacity")
            }
            bootstrapSD = try bootstrap(samples: samples, stateCount: stateCount, counts: counts,
                configuration: cfg)
        }
        let estimates = (0..<stateCount).map { index in
            VivoMBARStateEstimate(stateIndex: index,
                relativeFreeEnergy: core.freeEnergies[index],
                bootstrapStandardDeviation: bootstrapSD[index],
                effectiveSamples: effective[index])
        }
        return .init(schema: VivoMBARResult.schema, stateCount: stateCount,
            sampleCount: samples.count, originCounts: counts, estimates: estimates,
            pairwiseOverlap: overlap, iterations: core.iterations,
            maximumResidual: core.residual, converged: issues.isEmpty,
            overlapGraphConnected: connected, issues: issues, interpretation: interpretation)
    }

    private struct CoreResult {
        let freeEnergies: [Double]
        let iterations: Int
        let residual: Double
        let converged: Bool
    }

    private static func solveCore(samples: [VivoMBARReducedPotentialSample], stateCount: Int,
                                  counts: [Int], tolerance: Double,
                                  maximumIterations: Int) throws -> CoreResult {
        let logCounts = counts.map { log(Double($0)) }
        var f = [Double](repeating: 0, count: stateCount)
        var residual = Double.infinity
        for iteration in 1...maximumIterations {
            var logDenominators = [Double](repeating: 0, count: samples.count)
            for n in samples.indices {
                var terms = [Double](repeating: 0, count: stateCount)
                for k in 0..<stateCount { terms[k] = logCounts[k] + f[k] - samples[n].reducedPotentials[k] }
                logDenominators[n] = try logSumExp(terms)
            }
            var next = [Double](repeating: 0, count: stateCount)
            for k in 0..<stateCount {
                var terms = [Double](repeating: 0, count: samples.count)
                for n in samples.indices { terms[n] = -samples[n].reducedPotentials[k] - logDenominators[n] }
                next[k] = -(try logSumExp(terms))
            }
            let reference = next[0]
            for k in next.indices { next[k] -= reference }
            residual = zip(next, f).reduce(0.0) { max($0, abs($1.0 - $1.1)) }
            guard residual.isFinite, next.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.convergence("nonfinite MBAR fixed-point iteration")
            }
            f = next
            if residual <= tolerance { return .init(freeEnergies: f, iterations: iteration, residual: residual, converged: true) }
        }
        return .init(freeEnergies: f, iterations: maximumIterations, residual: residual, converged: false)
    }

    private static func targetEffectiveSamples(samples: [VivoMBARReducedPotentialSample],
                                               counts: [Int], freeEnergies f: [Double]) -> [Double] {
        let kCount = counts.count, logCounts = counts.map { log(Double($0)) }
        var sumSquares = [Double](repeating: 0, count: kCount)
        for sample in samples {
            var maxTerm = -Double.infinity
            for k in 0..<kCount { maxTerm = max(maxTerm, logCounts[k] + f[k] - sample.reducedPotentials[k]) }
            var denom = 0.0
            for k in 0..<kCount { denom += exp(logCounts[k] + f[k] - sample.reducedPotentials[k] - maxTerm) }
            let logDen = maxTerm + log(denom)
            for k in 0..<kCount {
                let w = exp(f[k] - sample.reducedPotentials[k] - logDen)
                sumSquares[k] += w * w
            }
        }
        return sumSquares.map { $0 > 0 ? 1 / $0 : 0 }
    }

    private static func pairwiseOverlap(samples: [VivoMBARReducedPotentialSample], counts: [Int],
                                        freeEnergies f: [Double], threshold: Double) -> [VivoMBARPairwiseOverlap] {
        let kCount = counts.count
        func essFraction(from: Int, to: Int) -> Double {
            var logs: [Double] = []
            logs.reserveCapacity(counts[from])
            for sample in samples where sample.originStateIndex == from {
                logs.append(-(sample.reducedPotentials[to] - sample.reducedPotentials[from]))
            }
            guard let maximum = logs.max(), maximum.isFinite else { return 0 }
            var sum = 0.0, sum2 = 0.0
            for value in logs {
                let w = exp(value - maximum); sum += w; sum2 += w * w
            }
            guard sum2 > 0 else { return 0 }
            return min(1, max(0, (sum * sum / sum2) / Double(logs.count)))
        }
        func meanAcceptance(from: Int, to: Int) -> Double {
            let deltaF = f[to] - f[from]
            var sum = 0.0, count = 0
            for sample in samples where sample.originStateIndex == from {
                let corrected = (sample.reducedPotentials[to] - sample.reducedPotentials[from]) - deltaF
                let probability = corrected <= 0 ? 1.0 : exp(-min(corrected, 745))
                sum += probability; count += 1
            }
            return count > 0 ? sum / Double(count) : 0
        }
        var result: [VivoMBARPairwiseOverlap] = []
        for a in 0..<kCount { for b in (a + 1)..<kCount {
            let essAB = essFraction(from: a, to: b), essBA = essFraction(from: b, to: a)
            let acceptanceAB = meanAcceptance(from: a, to: b)
            let acceptanceBA = meanAcceptance(from: b, to: a)
            let symmetric = min(acceptanceAB, acceptanceBA)
            result.append(.init(stateA: a, stateB: b,
                reweightingESSFractionAToB: essAB, reweightingESSFractionBToA: essBA,
                meanAcceptanceAToB: acceptanceAB, meanAcceptanceBToA: acceptanceBA,
                symmetricAcceptanceOverlap: symmetric, connected: symmetric >= threshold))
        } }
        return result
    }

    private static func graphConnected(stateCount: Int, overlap: [VivoMBARPairwiseOverlap]) -> Bool {
        var adjacency = [[Int]](repeating: [], count: stateCount)
        for edge in overlap where edge.connected {
            adjacency[edge.stateA].append(edge.stateB); adjacency[edge.stateB].append(edge.stateA)
        }
        var visited = Set([0]), queue = [0], cursor = 0
        while cursor < queue.count {
            let node = queue[cursor]; cursor += 1
            for neighbor in adjacency[node] where !visited.contains(neighbor) {
                visited.insert(neighbor); queue.append(neighbor)
            }
        }
        return visited.count == stateCount
    }

    private static func bootstrap(samples: [VivoMBARReducedPotentialSample], stateCount: Int,
                                  counts: [Int], configuration cfg: VivoMBARConfiguration) throws -> [Double?] {
        var grouped = [[VivoMBARReducedPotentialSample]](repeating: [], count: stateCount)
        for sample in samples { grouped[sample.originStateIndex].append(sample) }
        var random = VivoSplitMix64(state: cfg.bootstrapSeed)
        var values = [[Double]](repeating: [], count: stateCount)
        for replicate in 0..<cfg.bootstrapReplicates {
            var resampled: [VivoMBARReducedPotentialSample] = []
            resampled.reserveCapacity(samples.count)
            for origin in 0..<stateCount {
                for draw in 0..<counts[origin] {
                    let chosen = Int(random.next() % UInt64(grouped[origin].count))
                    var sample = grouped[origin][chosen]
                    sample.identifier = "bootstrap-\(replicate)-\(origin)-\(draw)"
                    resampled.append(sample)
                }
            }
            let core = try solveCore(samples: resampled, stateCount: stateCount, counts: counts,
                                     tolerance: cfg.residualTolerance, maximumIterations: cfg.maximumIterations)
            guard core.converged else { throw VivoChemistryError.convergence("MBAR bootstrap replicate did not converge") }
            for k in 0..<stateCount { values[k].append(core.freeEnergies[k]) }
        }
        return values.map { series in
            guard series.count > 1 else { return nil }
            let mean = series.reduce(0, +) / Double(series.count)
            let variance = series.reduce(0) { $0 + ($1 - mean) * ($1 - mean) } / Double(series.count - 1)
            return sqrt(max(0, variance))
        }
    }

    private static func logSumExp(_ values: [Double]) throws -> Double {
        guard !values.isEmpty, let maximum = values.max(), maximum.isFinite else {
            throw VivoChemistryError.invalid("MBAR log-sum-exp input")
        }
        let scaled = values.reduce(0.0) { $0 + exp($1 - maximum) }
        guard scaled.isFinite, scaled > 0 else { throw VivoChemistryError.convergence("MBAR log-sum-exp underflow") }
        return maximum + log(scaled)
    }
}