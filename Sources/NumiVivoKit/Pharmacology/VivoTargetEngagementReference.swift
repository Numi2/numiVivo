import Foundation

public struct VivoTargetEngagementSample: Codable, Equatable, Sendable {
    public let timeSeconds: Double
    public let unboundDrugM: Double
    public let fractions: VivoTargetFractions
    public let totalTargetM: Double
    public var drugOccupancy: Double { fractions.drugOccupancy }
    public var covalentOccupancy: Double { fractions.covalent / fractions.total }
}

public struct VivoTargetEngagementNumerics: Codable, Equatable, Sendable {
    public let maximumMatrixSquarings: Int
    public let maximumPropagations: Int
    public let fractionConservationTolerance: Double
    public init(maximumMatrixSquarings: Int = 48, maximumPropagations: Int = 1_048_576,
                fractionConservationTolerance: Double = 1e-8) {
        self.maximumMatrixSquarings = maximumMatrixSquarings
        self.maximumPropagations = maximumPropagations
        self.fractionConservationTolerance = fractionConservationTolerance
    }
    public func validate() throws {
        guard (0...60).contains(maximumMatrixSquarings), (1...1_048_576).contains(maximumPropagations),
              fractionConservationTolerance.isFinite, fractionConservationTolerance >= 1e-14,
              fractionConservationTolerance <= 1e-3 else {
            throw VivoKineticsError.invalid("reference numerical policy")
        }
    }
}

public struct VivoTargetEngagementResult: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let backend: String
    public let samples: [VivoTargetEngagementSample]
    public let propagationCount: Int
    public let maximumScalingDepth: Int
    public let maximumFractionMassError: Double
    public let containsAssumedParameters: Bool
    public let hasUnknownParameterUncertainty: Bool
    public let limitations: [String]
}

/// Deterministic FP64 reference for externally prescribed free-drug exposure.
/// It solves the affine four-state binding/turnover system by a nonnegative
/// uniformization series with scaling and squaring. It does not integrate a
/// second physiology model, infer clinical outcomes, or estimate rate constants.
public enum VivoTargetEngagementReference {
    public static func run(_ experiment: VivoTargetEngagementExperiment,
                           numerics: VivoTargetEngagementNumerics = .init()) throws -> VivoTargetEngagementResult {
        try experiment.validate(); try numerics.validate()
        try experiment.initial.validate(hasCompetitor: experiment.kinetics.competitor != nil,
                                        tolerance: numerics.fractionConservationTolerance)
        let model = experiment.kinetics, knots = experiment.exposure.knots
        var state = experiment.initial.values + [1.0]
        var time = 0.0, knotIndex = 0, propagationCount = 0, maximumDepth = 0
        var maximumError = abs(experiment.initial.total - 1)
        var samples: [VivoTargetEngagementSample] = []
        samples.reserveCapacity(experiment.sampleTimesSeconds.count)
        for sampleTime in experiment.sampleTimesSeconds {
            while time < sampleTime {
                while knotIndex + 1 < knots.count, knots[knotIndex + 1].timeSeconds <= time { knotIndex += 1 }
                guard knotIndex + 1 < knots.count else { throw VivoKineticsError.invalid("exposure support exhausted") }
                let boundary = min(sampleTime, knots[knotIndex + 1].timeSeconds)
                let dt = boundary - time
                guard dt.isFinite, dt > 0, propagationCount < numerics.maximumPropagations else {
                    throw VivoKineticsError.capacity("time progress or propagation budget")
                }
                let (transition, depth) = try propagator(model, drugM: knots[knotIndex].unboundDrugM,
                                                        dt: dt, maximumDepth: numerics.maximumMatrixSquarings)
                state = multiply(transition, vector: state)
                // The last coordinate represents a constant source, not mass.
                guard state.allSatisfy({ $0.isFinite && $0 >= 0 }), state[4] == 1 else {
                    throw VivoKineticsError.numerical("nonfinite/negative state or corrupted source coordinate")
                }
                let total = state.prefix(4).reduce(0, +)
                let error = abs(total - 1)
                maximumError = max(maximumError, error)
                guard error <= numerics.fractionConservationTolerance else {
                    throw VivoKineticsError.numerical("target balance failed; no clipping or renormalization was applied")
                }
                maximumDepth = max(maximumDepth, depth); propagationCount += 1
                time = boundary
            }
            while knotIndex + 1 < knots.count, knots[knotIndex + 1].timeSeconds <= sampleTime { knotIndex += 1 }
            let fractions = VivoTargetFractions(free: state[0], reversible: state[1], covalent: state[2], competitor: state[3])
            let totalTarget = fractions.total * model.baselineTarget.value
            guard totalTarget.isFinite, totalTarget > 0 else { throw VivoKineticsError.numerical("target amount overflow/underflow") }
            samples.append(.init(timeSeconds: sampleTime, unboundDrugM: knots[knotIndex].unboundDrugM,
                                 fractions: fractions, totalTargetM: totalTarget))
        }
        return .init(schemaVersion: 1, backend: "numivivo.target-engagement.fp64-uniformization.v1",
                     samples: samples, propagationCount: propagationCount, maximumScalingDepth: maximumDepth,
                     maximumFractionMassError: maximumError, containsAssumedParameters: model.containsAssumptions,
                     hasUnknownParameterUncertainty: model.hasUnknownParameterUncertainty,
                     limitations: [
                        "Conditional model output, not evidence of biological efficacy or clinical safety.",
                        "Exposure is externally maintained UNBOUND concentration; ligand depletion and finite drug mass balance are not represented.",
                        "All target states share one loss rate; synthesis equals loss rate times baseline target abundance.",
                        "A single constant-concentration reversible competitor is supported; metabolite binding and state-dependent turnover are not represented.",
                        "Exposure is right-continuous piecewise constant. Errors in exposure reconstruction and model parameters are not numerical solver error.",
                        "Parameter uncertainty is recorded but this deterministic trajectory is not a predictive interval."
                     ])
    }

    private static let dimension = 5
    private static func identity() -> [Double] {
        var result = [Double](repeating: 0, count: dimension * dimension)
        for i in 0..<dimension { result[i * dimension + i] = 1 }
        return result
    }
    private static func multiply(_ a: [Double], _ b: [Double]) -> [Double] {
        var result = [Double](repeating: 0, count: dimension * dimension)
        for i in 0..<dimension {
            for k in 0..<dimension where a[i * dimension + k] != 0 {
                for j in 0..<dimension {
                    result[i * dimension + j] += a[i * dimension + k] * b[k * dimension + j]
                }
            }
        }
        return result
    }
    private static func multiply(_ a: [Double], vector: [Double]) -> [Double] {
        (0..<dimension).map { i in
            (0..<dimension).reduce(0) { $0 + a[i * dimension + $1] * vector[$1] }
        }
    }

    private static func propagator(_ model: VivoCovalentKineticPack, drugM: Double,
                                   dt: Double, maximumDepth: Int) throws -> ([Double], Int) {
        let a = model.association.value * drugM
        let b = model.dissociation.value, k = model.inactivation.value, d = model.targetTurnover.value
        let c = (model.competitor?.association.value ?? 0) * (model.competitor?.unboundConcentration.value ?? 0)
        let q = model.competitor?.dissociation.value ?? 0
        let rate = max(a + c + d, b + k + d, q + d, d)
        guard rate > 0 else { return (identity(), 0) }
        let scaledTime = rate * dt
        guard scaledTime.isFinite else { throw VivoKineticsError.numerical("rate-times-interval overflow") }
        let depthValue = scaledTime > 0.5 ? ceil(log2(scaledTime) + 1) : 0
        guard depthValue <= Double(maximumDepth) else {
            throw VivoKineticsError.capacity("stiff propagation exceeds declared matrix scaling budget")
        }
        let depth = Int(depthValue)
        let mu = rate * (dt / pow(2, Double(depth)))
        guard mu.isFinite, mu >= 0, mu <= 0.50000000000001 else {
            throw VivoKineticsError.numerical("invalid scaled propagation interval")
        }
        // P = I + A/rate in augmented coordinates [free, reversible, covalent,
        // competitor, 1]. P is nonnegative; its maximum column sum is <= 2.
        var p = identity()
        p[0] -= (a + c + d) / rate; p[1] = b / rate; p[3] = q / rate; p[4] = d / rate
        p[5] = a / rate; p[6] -= (b + k + d) / rate
        p[11] = k / rate; p[12] -= d / rate
        p[15] = c / rate; p[18] -= (q + d) / rate
        guard p.allSatisfy({ $0.isFinite && $0 >= 0 }) else { throw VivoKineticsError.numerical("invalid nonnegative transition matrix") }
        var power = identity(), sum = identity(), weight = 1.0
        // mu * ||P||_1 <= 1. The 32-term truncation is below FP64 rounding
        // before squaring; floating point balance is checked at every boundary.
        for order in 1...32 {
            power = multiply(power, p)
            weight *= mu / Double(order)
            for i in sum.indices { sum[i] += weight * power[i] }
        }
        let prefactor = exp(-mu)
        for i in sum.indices { sum[i] *= prefactor }
        // Exact known affine coordinate, independent of truncation roundoff.
        for j in 0..<dimension { sum[20 + j] = j == 4 ? 1 : 0 }
        for _ in 0..<depth { sum = multiply(sum, sum) }
        return (sum, depth)
    }
}
