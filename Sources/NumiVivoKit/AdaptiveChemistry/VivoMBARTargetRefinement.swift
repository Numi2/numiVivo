import Foundation

/// One physical configuration, evaluated under every sampled and target
/// Hamiltonian. Reduced potentials are dimensionless and use one common measure.
/// Block IDs declare independent units; naming them does not prove independence.
public struct VivoMBARTargetSample: Codable, Sendable, Equatable {
    public var identifier: String
    public var originStateIndex: Int
    public var independentBlockIdentifier: String
    public var coordinate: Double
    public var sampledReducedPotentials: [Double]
    public var targetReducedPotentials: [Double]
    public init(identifier: String, originStateIndex: Int, independentBlockIdentifier: String,
                coordinate: Double, sampledReducedPotentials: [Double], targetReducedPotentials: [Double]) {
        self.identifier = identifier; self.originStateIndex = originStateIndex
        self.independentBlockIdentifier = independentBlockIdentifier; self.coordinate = coordinate
        self.sampledReducedPotentials = sampledReducedPotentials; self.targetReducedPotentials = targetReducedPotentials
    }
}
public struct VivoMBARTargetConfiguration: Codable, Sendable, Equatable {
    public var sourceSolver: VivoMBARConfiguration
    public var bootstrapReplicates: Int
    public var bootstrapSeed: UInt64
    public var minimumIndependentBlocksPerOrigin: Int
    public var minimumBinEffectiveSamples: Double
    public var maximumBinWeight: Double
    public var maximumPrimitiveElements: Int
    public init(sourceSolver: VivoMBARConfiguration = .init(bootstrapReplicates: 0),
                bootstrapReplicates: Int = 64, bootstrapSeed: UInt64 = 0x544152474554,
                minimumIndependentBlocksPerOrigin: Int = 4, minimumBinEffectiveSamples: Double = 10,
                maximumBinWeight: Double = 0.25, maximumPrimitiveElements: Int = 200_000_000) {
        self.sourceSolver = sourceSolver; self.bootstrapReplicates = bootstrapReplicates
        self.bootstrapSeed = bootstrapSeed; self.minimumIndependentBlocksPerOrigin = minimumIndependentBlocksPerOrigin
        self.minimumBinEffectiveSamples = minimumBinEffectiveSamples; self.maximumBinWeight = maximumBinWeight
        self.maximumPrimitiveElements = maximumPrimitiveElements
    }
    public func validate() throws {
        try sourceSolver.validate()
        guard sourceSolver.bootstrapReplicates == 0,
              bootstrapReplicates == 0 || (32...1024).contains(bootstrapReplicates),
              (2...4096).contains(minimumIndependentBlocksPerOrigin),
              minimumBinEffectiveSamples.isFinite, minimumBinEffectiveSamples >= 2,
              maximumBinWeight.isFinite, maximumBinWeight > 0, maximumBinWeight <= 1,
              maximumPrimitiveElements > 0 else {
            throw VivoChemistryError.invalid("target MBAR block-bootstrap or work policy")
        }
    }
}
public struct VivoMBARTargetRefinementRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mbar-target-refinement/v1"
    public var schema: String
    public var sampledStateIdentifiers: [String]
    public var targetIdentifiers: [String]
    /// A match requires equality of the entire reduced-potential column. An
    /// unmatched target is exploratory even when its importance weights look good.
    public var matchingSampledStateIndices: [Int?]
    public var referenceTargetIndex: Int
    public var binEdges: [Double]
    public var coordinateUnit: String
    public var referenceBinIndex: Int
    public var samples: [VivoMBARTargetSample]
    public var configuration: VivoMBARTargetConfiguration
    public init(sampledStateIdentifiers: [String], targetIdentifiers: [String],
                matchingSampledStateIndices: [Int?], referenceTargetIndex: Int = 0,
                binEdges: [Double], coordinateUnit: String, referenceBinIndex: Int = 0,
                samples: [VivoMBARTargetSample], configuration: VivoMBARTargetConfiguration = .init()) {
        schema = Self.schema; self.sampledStateIdentifiers = sampledStateIdentifiers
        self.targetIdentifiers = targetIdentifiers; self.matchingSampledStateIndices = matchingSampledStateIndices
        self.referenceTargetIndex = referenceTargetIndex; self.binEdges = binEdges; self.coordinateUnit = coordinateUnit
        self.referenceBinIndex = referenceBinIndex; self.samples = samples; self.configuration = configuration
    }
    public func validate() throws {
        try configuration.validate()
        let k = sampledStateIdentifiers.count, t = targetIdentifiers.count
        func validNames(_ names: [String]) -> Bool {
            Set(names).count == names.count && names.allSatisfy { !$0.isEmpty && $0.utf8.count <= 512 }
        }
        guard schema == Self.schema, (1...128).contains(k), (2...16).contains(t),
              validNames(sampledStateIdentifiers), validNames(targetIdentifiers),
              matchingSampledStateIndices.count == t, (0..<t).contains(referenceTargetIndex),
              (3...65).contains(binEdges.count), binEdges.allSatisfy(\.isFinite),
              zip(binEdges,binEdges.dropFirst()).allSatisfy({ $0 < $1 }),
              !coordinateUnit.isEmpty, coordinateUnit.utf8.count <= 128,
              (0..<(binEdges.count-1)).contains(referenceBinIndex),
              (2...1_000_000).contains(samples.count), validNames(samples.map(\.identifier)) else {
            throw VivoChemistryError.invalid("target MBAR identities, bins or dimensions")
        }
        let elements = samples.count.multipliedReportingOverflow(by: k+t)
        guard !elements.overflow, elements.partialValue <= configuration.maximumPrimitiveElements,
              elements.partialValue <= configuration.sourceSolver.maximumWorkElements else {
            throw VivoChemistryError.resourceLimit("target MBAR input matrix")
        }
        for match in matchingSampledStateIndices.compactMap({ $0 }) where !(0..<k).contains(match) {
            throw VivoChemistryError.invalid("target MBAR sampled-state match")
        }
        var origins: [String: Int] = [:]
        for sample in samples {
            guard (0..<k).contains(sample.originStateIndex), sample.coordinate.isFinite,
                  sample.coordinate >= binEdges[0], sample.coordinate <= binEdges.last!,
                  !sample.independentBlockIdentifier.isEmpty, sample.independentBlockIdentifier.utf8.count <= 512,
                  sample.sampledReducedPotentials.count == k, sample.targetReducedPotentials.count == t,
                  (sample.sampledReducedPotentials + sample.targetReducedPotentials).allSatisfy({ $0.isFinite && abs($0) <= 1e12 }) else {
                throw VivoChemistryError.invalid("target MBAR sample, coordinate coverage or potential matrix")
            }
            if let prior = origins[sample.independentBlockIdentifier], prior != sample.originStateIndex {
                throw VivoChemistryError.invalid("block crosses independent origins; supply joint-exchange block analysis instead")
            }
            origins[sample.independentBlockIdentifier] = sample.originStateIndex
            for target in 0..<t {
                if let origin = matchingSampledStateIndices[target],
                   abs(sample.targetReducedPotentials[target]-sample.sampledReducedPotentials[origin]) > 1e-10 {
                    throw VivoChemistryError.invalid("declared sampled target does not match its physical potential column")
                }
            }
        }
        // Equal block lengths within an origin preserve its sampling mixture
        // under resampling. Trailing partial blocks must be handled explicitly.
        for origin in 0..<k {
            let grouped = Dictionary(grouping: samples.filter { $0.originStateIndex == origin },by: \.independentBlockIdentifier)
            guard !grouped.isEmpty, Set(grouped.values.map(\.count)).count == 1,
                  grouped.values.reduce(0,{ $0+$1.count }) >= 2 else {
                throw VivoChemistryError.invalid("target MBAR needs equal nonempty blocks and two samples per origin")
            }
        }
    }
}
public struct VivoMBARTargetBin: Codable, Sendable, Equatable {
    public let index: Int
    public let probability: Double
    public let effectiveSamples: Double
    public let maximumConditionalWeight: Double
    public let relativeFreeEnergy: Double?
}
public struct VivoMBARTargetProfile: Codable, Sendable, Equatable {
    public let identifier: String
    public let relativeFreeEnergy: Double
    public let effectiveSamples: Double
    public let sampledTarget: Bool
    public let bins: [VivoMBARTargetBin]
}
public struct VivoMBARTargetCorrection: Codable, Sendable, Equatable {
    public let targetIdentifier: String
    public let binIndex: Int
    /// Difference between target and reference bin-relative free energies, in kT.
    public let deltaRelativeFreeEnergy: Double?
    public let bootstrapStandardDeviation: Double?
}
public enum VivoMBARTargetQualification: String, Codable, Sendable {
    case sampledTargetsPassDiagnostics, exploratoryUnsampledTarget, insufficientEvidence
}
public struct VivoMBARTargetRefinementResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/mbar-target-refinement-result/v1"
    public let schema: String
    public let source: VivoMBARResult
    public let profiles: [VivoMBARTargetProfile]
    public let corrections: [VivoMBARTargetCorrection]
    /// Corrections order is explicit above. All targets and bins share each
    /// bootstrap draw; covariance is not discarded before forming contrasts.
    public let correctionCovariance: [[Double]]?
    public let independentBlockCounts: [Int]
    public let completedBootstrapReplicates: Int
    public let chargedPrimitiveElements: Int
    public let qualification: VivoMBARTargetQualification
    public let issues: [String]
    public let interpretation: String
}

/// Reuses the native MBAR estimating-equation solver. This supplies equilibrium
/// bin-probability corrections, not continuous PMF derivatives, free energies of
/// a dividing surface, a transmission coefficient, or a kinetic rate.
public enum VivoMBARTargetRefinement {
    public static let interpretation = "Common-measure equilibrium reweighting with shared stratified block-bootstrap covariance. Every target must be independently sampled, connected by the native overlap graph, and pass per-bin weight diagnostics for sampled-target qualification. Declared block independence and unobserved phase-space support are not established by this calculation. No dynamical transmission or rate is transferred between Hamiltonians."
    private static func lse(_ values: [Double]) -> Double? {
        guard let maximum = values.max(), maximum.isFinite else { return nil }
        return maximum + log(values.reduce(0) { $0 + exp($1-maximum) })
    }
    private struct Estimate {
        let source: VivoMBARResult
        let profiles: [VivoMBARTargetProfile]
        let corrections: [Double?]
    }
    public static func calculate(_ request: VivoMBARTargetRefinementRequest) throws -> VivoMBARTargetRefinementResult {
        try request.validate()
        let cfg = request.configuration, k = request.sampledStateIdentifiers.count
        let t = request.targetIdentifiers.count, bins = request.binEdges.count-1
        var work = 0
        func charge(_ amount: Int) throws {
            guard amount >= 0, amount <= cfg.maximumPrimitiveElements-work else {
                throw VivoChemistryError.resourceLimit("target MBAR aggregate primitive work")
            }
            work += amount
        }
        func estimate(_ samples: [VivoMBARTargetSample]) throws -> Estimate {
            let n = samples.count
            let fixed = n*(20*k*k+8*t*bins+20)
            try charge(fixed)
            let perIteration = max(1,8*n*k)
            var solver = cfg.sourceSolver
            solver.maximumIterations = min(solver.maximumIterations,(cfg.maximumPrimitiveElements-work)/perIteration)
            guard solver.maximumIterations > 0 else { throw VivoChemistryError.resourceLimit("target MBAR iteration allowance") }
            let source = try VivoMultistateMBAR.solve(samples: samples.map {
                .init(identifier: $0.identifier,originStateIndex: $0.originStateIndex,reducedPotentials: $0.sampledReducedPotentials)
            },stateCount: k,configuration: solver)
            try charge(source.iterations*perIteration)
            let f = source.estimates.map(\.relativeFreeEnergy), logCounts = source.originCounts.map { log(Double($0)) }
            let denominator = samples.map { sample in
                lse((0..<k).map { logCounts[$0]+f[$0]-sample.sampledReducedPotentials[$0] })!
            }
            let membership = samples.map { sample -> Int in
                if sample.coordinate == request.binEdges.last! { return bins-1 }
                return (0..<bins).first { sample.coordinate >= request.binEdges[$0] && sample.coordinate < request.binEdges[$0+1] }!
            }
            var profiles: [VivoMBARTargetProfile] = []
            for target in 0..<t {
                let logs = samples.indices.map { -samples[$0].targetReducedPotentials[target]-denominator[$0] }
                let logZ = lse(logs)!, weights = logs.map { exp($0-logZ) }
                let binZ = (0..<bins).map { bin in lse(logs.indices.filter { membership[$0] == bin }.map { logs[$0] }) }
                var observations: [VivoMBARTargetBin] = []
                for bin in 0..<bins {
                    let z = binZ[bin]
                    let conditional = z.map { z in logs.indices.filter { membership[$0] == bin }.map { exp(logs[$0]-z) } } ?? []
                    let sum2 = conditional.reduce(0) { $0+$1*$1 }
                    let relative: Double?
                    if let z, let ref = binZ[request.referenceBinIndex] { relative = ref-z } else { relative = nil }
                    observations.append(.init(index: bin,probability: z.map { exp($0-logZ) } ?? 0,
                        effectiveSamples: sum2 > 0 ? 1/sum2 : 0,maximumConditionalWeight: conditional.max() ?? 1,
                        relativeFreeEnergy: relative))
                }
                profiles.append(.init(identifier: request.targetIdentifiers[target],relativeFreeEnergy: -logZ,
                    effectiveSamples: 1/weights.reduce(0) { $0+$1*$1 },sampledTarget: request.matchingSampledStateIndices[target] != nil,
                    bins: observations))
            }
            var differences: [Double?] = []
            for target in 0..<t where target != request.referenceTargetIndex {
                for bin in 0..<bins {
                    if let a = profiles[target].bins[bin].relativeFreeEnergy,
                       let b = profiles[request.referenceTargetIndex].bins[bin].relativeFreeEnergy { differences.append(a-b) }
                    else { differences.append(nil) }
                }
            }
            return .init(source: source,profiles: profiles,corrections: differences)
        }
        let base = try estimate(request.samples)
        let grouped = (0..<k).map { origin in
            let blocks = Dictionary(grouping: request.samples.filter { $0.originStateIndex == origin },by: \.independentBlockIdentifier)
            return blocks.keys.sorted().map { blocks[$0]! }
        }
        let counts = grouped.map(\.count)
        var issues = base.source.issues
        for origin in 0..<k where counts[origin] < cfg.minimumIndependentBlocksPerOrigin {
            issues.append("origin \(origin) has too few declared independent blocks")
        }
        for profile in base.profiles {
            for bin in profile.bins where bin.effectiveSamples < cfg.minimumBinEffectiveSamples || bin.maximumConditionalWeight > cfg.maximumBinWeight {
                issues.append("\(profile.identifier) bin \(bin.index) lacks adequate effective support")
            }
        }
        if cfg.bootstrapReplicates == 0 { issues.append("paired block-bootstrap uncertainty was not requested") }
        let d = base.corrections.count
        guard d*d <= cfg.maximumPrimitiveElements else { throw VivoChemistryError.resourceLimit("correction covariance capacity") }
        var mean = [Double](repeating: 0,count: d), m2 = [[Double]](repeating: [Double](repeating: 0,count: d),count: d)
        var completed = 0, bootstrapFailed = false
        if cfg.bootstrapReplicates > 0 && base.source.converged && base.corrections.allSatisfy({ $0 != nil }) {
            var random = VivoSplitMix64(state: cfg.bootstrapSeed)
            for replicate in 0..<cfg.bootstrapReplicates {
                try charge(d*d+request.samples.count)
                var resampled: [VivoMBARTargetSample] = []
                for origin in 0..<k {
                    for draw in grouped[origin].indices {
                        let chosen = Int(random.next() % UInt64(grouped[origin].count))
                        for (offset,original) in grouped[origin][chosen].enumerated() {
                            var sample = original
                            sample.identifier = "resample-\(replicate)-\(origin)-\(draw)-\(offset)"
                            resampled.append(sample)
                        }
                    }
                }
                let trial = try estimate(resampled)
                guard trial.source.converged, trial.corrections.allSatisfy({ $0 != nil }) else {
                    bootstrapFailed = true; issues.append("bootstrap replicate \(replicate) lacks source convergence or bin support"); break
                }
                let values = trial.corrections.map { $0! }; completed += 1
                let delta = zip(values,mean).map { $0-$1 }
                for i in 0..<d { mean[i] += delta[i]/Double(completed) }
                for i in 0..<d { for j in 0..<d { m2[i][j] += delta[i]*(values[j]-mean[j]) } }
            }
        }
        let covariance: [[Double]]?
        if completed == cfg.bootstrapReplicates && completed >= 32 && !bootstrapFailed {
            let denominator = 2 * Double(completed - 1)
            var symmetric = [[Double]](repeating: [Double](repeating: 0, count: d), count: d)
            for i in 0..<d {
                for j in 0..<d {
                    symmetric[i][j] = (m2[i][j] + m2[j][i]) / denominator
                }
            }
            covariance = symmetric
        } else { covariance = nil }
        if covariance == nil && cfg.bootstrapReplicates > 0 { issues.append("complete shared-bootstrap covariance is unavailable") }
        var corrections: [VivoMBARTargetCorrection] = [], cursor = 0
        for target in 0..<t where target != request.referenceTargetIndex {
            for bin in 0..<bins {
                corrections.append(.init(targetIdentifier: request.targetIdentifiers[target],binIndex: bin,
                    deltaRelativeFreeEnergy: base.corrections[cursor],
                    bootstrapStandardDeviation: covariance.map { sqrt(max(0,$0[cursor][cursor])) }))
                cursor += 1
            }
        }
        let sampled = request.matchingSampledStateIndices.allSatisfy { $0 != nil }
        let qualification: VivoMBARTargetQualification = !issues.isEmpty ? .insufficientEvidence
            : (sampled ? .sampledTargetsPassDiagnostics : .exploratoryUnsampledTarget)
        if !sampled { issues.append("unsampled targets remain exploratory; weight diagnostics do not establish missing target support") }
        return .init(schema: VivoMBARTargetRefinementResult.schema,source: base.source,profiles: base.profiles,
            corrections: corrections,correctionCovariance: covariance,independentBlockCounts: counts,
            completedBootstrapReplicates: completed,chargedPrimitiveElements: work,qualification: qualification,
            issues: issues,interpretation: interpretation)
    }
    public static func validate(_ result: VivoMBARTargetRefinementResult,request: VivoMBARTargetRefinementRequest) throws {
        guard result == (try calculate(request)) else { throw VivoChemistryError.invalid("target MBAR evidence does not reconstruct") }
    }
}
