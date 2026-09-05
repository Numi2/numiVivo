import Foundation

/// Censoring describes information actually observed, not a substituted value.
/// An interval is a reported interval for one observation, NOT a cohort range.
public enum VivoAssaySupport: Codable, Equatable, Sendable {
    case exact
    case below(limit: Double)
    case above(limit: Double)
    case interval(lower: Double, upper: Double)

    public func validate() throws {
        switch self {
        case .exact: break
        case .below(let limit), .above(let limit):
            guard limit.isFinite else { throw VivoPosteriorError.invalid("nonfinite censoring limit") }
        case .interval(let lower, let upper):
            guard lower.isFinite, upper.isFinite, lower < upper, (upper - lower).isFinite else {
                throw VivoPosteriorError.invalid("invalid individual-observation interval")
            }
        }
    }
}

public struct VivoGaussianAssayNoise: Codable, Equatable, Sendable {
    public let scale: Double
    public let additionalStandardDeviation: Double
    public let bias: Double
    public let correlatedFraction: Double
    public let correlationTimeSeconds: Double
    public init(scale: Double = 1, additionalStandardDeviation: Double = 0, bias: Double = 0,
                correlatedFraction: Double = 0, correlationTimeSeconds: Double = 1) {
        self.scale = scale; self.additionalStandardDeviation = additionalStandardDeviation; self.bias = bias
        self.correlatedFraction = correlatedFraction; self.correlationTimeSeconds = correlationTimeSeconds
    }
    public func validate() throws {
        guard scale.isFinite, scale > 0, additionalStandardDeviation.isFinite, additionalStandardDeviation >= 0,
              bias.isFinite, correlatedFraction.isFinite, correlatedFraction >= 0, correlatedFraction < 1,
              correlationTimeSeconds.isFinite, correlationTimeSeconds > 0 else {
            throw VivoPosteriorError.invalid("assay noise scale, bias or correlation policy")
        }
    }
    /// Missing reported precision is usable only with explicit positive residual
    /// SD. That SD is a model assumption/parameter, not a fabricated measurement.
    public func standardDeviation(reported: Double?) throws -> Double {
        try validate()
        if let reported, !reported.isFinite || reported <= 0 {
            throw VivoPosteriorError.invalid("reported assay SD must be positive or unknown")
        }
        let sd = hypot((reported ?? 0) * scale, additionalStandardDeviation)
        guard sd.isFinite, sd > 0 else { throw VivoPosteriorError.invalid("assay uncertainty is unresolved") }
        return sd
    }
}

public struct VivoGaussianAssayDatum: Codable, Equatable, Sendable {
    public let timeSeconds: Double
    public let value: Double
    public let reportedStandardDeviation: Double?
    public let support: VivoAssaySupport
    public init(timeSeconds: Double, value: Double, reportedStandardDeviation: Double?, support: VivoAssaySupport = .exact) {
        self.timeSeconds = timeSeconds; self.value = value
        self.reportedStandardDeviation = reportedStandardDeviation; self.support = support
    }
}

/// Gaussian density for exact observations, probability for censored observations.
/// For exact correlated observations R_ij = rho exp(-|t_i-t_j|/tau), R_ii = 1.
/// Different cases are independent. Correlated censoring is rejected explicitly;
/// multiplying marginal censoring probabilities would be the wrong likelihood.
public enum VivoGaussianAssay {
    public static func logLikelihood(predictions: [Double], observations: [VivoGaussianAssayDatum],
                                     noise: VivoGaussianAssayNoise = .init()) throws -> Double {
        try noise.validate()
        guard !predictions.isEmpty, predictions.count == observations.count, predictions.count <= 4096,
              predictions.allSatisfy(\.isFinite) else { throw VivoPosteriorError.invalid("assay data dimensions") }
        var sd: [Double] = []
        for observation in observations {
            try observation.support.validate()
            guard observation.value.isFinite, observation.timeSeconds.isFinite, observation.timeSeconds >= 0 else {
                throw VivoPosteriorError.invalid("assay observation value or clock") }
            sd.append(try noise.standardDeviation(reported: observation.reportedStandardDeviation))
        }
        if noise.correlatedFraction == 0 {
            var sum = 0.0, correction = 0.0
            for i in predictions.indices {
                let value = try logContribution(mean: predictions[i] + noise.bias, sd: sd[i],
                                                value: observations[i].value, support: observations[i].support)
                let adjusted = value - correction, next = sum + adjusted
                correction = (next - sum) - adjusted; sum = next
            }
            guard sum.isFinite else { throw VivoPosteriorError.numerical("assay likelihood overflow") }
            return sum
        }
        let n = predictions.count
        guard n <= 256, observations.allSatisfy({ $0.support == .exact }) else {
            throw VivoPosteriorError.invalid("correlated Gaussian assays require <=256 exact observations; correlated censoring is unsupported")
        }
        var factor = [Double](repeating: 0, count: n * n)
        var residual = [Double](repeating: 0, count: n), logDetHalf = 0.0
        for i in 0..<n {
            for j in 0...i {
                var value = i == j ? 1 : noise.correlatedFraction * exp(-abs(observations[i].timeSeconds - observations[j].timeSeconds) / noise.correlationTimeSeconds)
                for k in 0..<j { value -= factor[i * n + k] * factor[j * n + k] }
                if i == j {
                    guard value.isFinite, value > 0 else { throw VivoPosteriorError.numerical("assay correlation factorization failed") }
                    factor[i * n + j] = sqrt(value)
                    logDetHalf += log(factor[i * n + j]) + log(sd[i])
                } else { factor[i * n + j] = value / factor[j * n + j] }
            }
            var z = (observations[i].value - predictions[i] - noise.bias) / sd[i]
            for j in 0..<i { z -= factor[i * n + j] * residual[j] }
            residual[i] = z / factor[i * n + i]
        }
        let result = -0.5 * residual.reduce(0) { $0 + $1 * $1 } - logDetHalf - 0.5 * Double(n) * log(2 * Double.pi)
        guard result.isFinite else { throw VivoPosteriorError.numerical("correlated assay likelihood overflow") }
        return result
    }

    public static func logContribution(mean: Double, sd: Double, value: Double,
                                       support: VivoAssaySupport = .exact) throws -> Double {
        try support.validate()
        guard mean.isFinite, sd.isFinite, sd > 0, value.isFinite else { throw VivoPosteriorError.invalid("Gaussian contribution") }
        let result: Double
        switch support {
        case .exact:
            let z = (value - mean) / sd
            result = -0.5 * z * z - log(sd) - 0.5 * log(2 * Double.pi)
        case .below(let limit): result = try logCDF((limit - mean) / sd)
        case .above(let limit): result = try logCDF((mean - limit) / sd)
        case .interval(let lower, let upper): result = try logInterval((lower - mean) / sd, (upper - mean) / sd)
        }
        guard result.isFinite else { throw VivoPosteriorError.numerical("Gaussian contribution exceeds finite support") }
        return result
    }

    /// Stable normal log CDF, including tails where erfc underflows. Mills'
    /// asymptotic series is truncated before its terms grow. No probability floor.
    public static func logCDF(_ z: Double) throws -> Double {
        guard z.isFinite, abs(z) <= 1e150 else { throw VivoPosteriorError.numerical("normal variate outside supported range") }
        if z >= 0 { return log1p(-0.5 * erfc(z / sqrt(2))) }
        if z > -10 { return log(0.5 * erfc(-z / sqrt(2))) }
        let inverseSquare = (1 / z) * (1 / z)
        var term = 1.0, sum = 1.0
        for i in 1...100 {
            let next = -term * Double(2 * i - 1) * inverseSquare
            if abs(next) >= abs(term) { break }
            sum += next; term = next
            if abs(term) < abs(sum) * 1e-17 { break }
        }
        return -0.5 * z * z - log(-z) - 0.5 * log(2 * Double.pi) + log(sum)
    }

    private static func logInterval(_ lower: Double, _ upper: Double) throws -> Double {
        guard lower.isFinite, upper.isFinite, lower < upper, (upper - lower).isFinite else {
            throw VivoPosteriorError.numerical("standardized censoring interval is unrepresentable")
        }
        let width = upper - lower, middle = lower + width * 0.5
        // Narrow intervals need direct integration: subtracting nearly identical
        // CDFs loses significant digits even when the probability is representable.
        if width * (1 + max(abs(lower), abs(upper))) < 1e-3 {
            let nodes = [0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363]
            let weights = [0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763]
            var integral = 0.0
            for i in 0..<4 {
                let dx = 0.5 * width * nodes[i]
                integral += weights[i] * (exp(-middle * dx - 0.5 * dx * dx) + exp(middle * dx - 0.5 * dx * dx))
            }
            return log(width) - 0.5 * middle * middle - 0.5 * log(2 * Double.pi) + log(0.5 * integral)
        }
        let high: Double, low: Double
        if lower > 0 { high = try logCDF(-lower); low = try logCDF(-upper) }
        else { high = try logCDF(upper); low = try logCDF(lower) }
        guard high > low else { throw VivoPosteriorError.numerical("censoring probability lost numerical resolution") }
        return high + log(-expm1(low - high))
    }

    public static func mixtureQuantile(means: [Double], standardDeviations: [Double], probability: Double) throws -> Double {
        guard !means.isEmpty, means.count == standardDeviations.count, means.count <= 8192,
              probability.isFinite, probability >= 0.0001, probability <= 0.9999,
              means.allSatisfy(\.isFinite), standardDeviations.allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoPosteriorError.invalid("Gaussian mixture dimensions, scales or probability")
        }
        var lower = zip(means, standardDeviations).map { $0 - 12 * $1 }.min()!
        var upper = zip(means, standardDeviations).map { $0 + 12 * $1 }.max()!
        guard lower.isFinite, upper.isFinite, lower < upper, (upper - lower).isFinite else {
            throw VivoPosteriorError.numerical("Gaussian mixture interval exceeds numerical resolution")
        }
        for _ in 0..<96 {
            let middle = lower + 0.5 * (upper - lower)
            if middle == lower || middle == upper { break }
            var probabilityAtMiddle = 0.0
            for i in means.indices {
                probabilityAtMiddle += 0.5 * erfc(-(middle - means[i]) / standardDeviations[i] / sqrt(2)) / Double(means.count)
            }
            if probabilityAtMiddle < probability { lower = middle } else { upper = middle }
        }
        return lower + 0.5 * (upper - lower)
    }
}
