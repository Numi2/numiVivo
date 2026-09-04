import Foundation

/// Bounded relaxed fixed-point update. Residual stopping is component-relative
/// and therefore invariant to a nonzero multiplicative change of channel units.
/// It is not a biological accuracy or nonlinear convergence certificate.
struct VivoCouplingIterationController: Sendable {
    enum Failure: Error { case invalidInput }
    private let initialRelaxation: Double
    private var relaxation: Double
    private var scales: [Double]?
    private var previousResidual: [Double]?

    init(relaxation: Double) throws {
        guard relaxation.isFinite, relaxation > 0, relaxation <= 1 else { throw Failure.invalidInput }
        initialRelaxation = relaxation
        self.relaxation = relaxation
    }
    static func residual(candidate: [Float], reference: [Float]) throws -> Double {
        guard candidate.count == reference.count else { throw Failure.invalidInput }
        var maximum = 0.0
        for (a, b) in zip(candidate, reference) {
            guard a.isFinite, b.isFinite else { throw Failure.invalidInput }
            let scale = max(abs(Double(a)), abs(Double(b)))
            let value = scale == 0 ? 0 : abs(Double(a) - Double(b)) / scale
            maximum = max(maximum, value)
        }
        return maximum
    }
    mutating func next(current: [Float], proposed: [Float]) throws -> [Float] {
        guard current.count == proposed.count,
              current.allSatisfy(\.isFinite), proposed.allSatisfy(\.isFinite) else { throw Failure.invalidInput }
        if scales == nil {
            scales = zip(current, proposed).map {
                let value = max(abs(Double($0.0)), abs(Double($0.1)))
                return value == 0 ? 1 : value
            }
        }
        guard let scales, scales.count == current.count else { throw Failure.invalidInput }
        let residual = current.indices.map { (Double(proposed[$0]) - Double(current[$0])) / scales[$0] }
        if let previousResidual {
            var numerator = 0.0, denominator = 0.0
            for i in residual.indices {
                let delta = residual[i] - previousResidual[i]
                numerator += previousResidual[i] * delta
                denominator += delta * delta
            }
            if denominator.isFinite, denominator > Double.leastNormalMagnitude, numerator.isFinite {
                let estimate = -relaxation * numerator / denominator
                // No extrapolation beyond the proposed state. Invalid or
                // nonpositive acceleration resets to the caller's safe factor.
                relaxation = estimate.isFinite && estimate > 0 ? min(1, max(0.05, estimate)) : initialRelaxation
            }
        }
        previousResidual = residual
        return try current.indices.map {
            let value = Double(current[$0]) + relaxation * (Double(proposed[$0]) - Double(current[$0]))
            guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else { throw Failure.invalidInput }
            return Float(value)
        }
    }
}
