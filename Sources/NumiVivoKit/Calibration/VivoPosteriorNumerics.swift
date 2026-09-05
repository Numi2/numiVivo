import Foundation

/// Small FP64 statistical operations. Dense covariance belongs here; the costly
/// forward calculations remain behind the batch evaluator. No GPU capability is
/// assumed and no posterior convergence claim follows from finite arithmetic.
internal enum VivoPosteriorNumerics {
    static func openUnit(_ random: inout VivoSplitMix64) -> Double {
        // 52-bit midpoint grid is strictly inside (0,1), including at the ends.
        (Double(random.next() >> 12) + 0.5) * 0x1.0p-52
    }

    static func normalizedWeights(logLikelihoods: [Double], delta: Double) throws -> ([Double], Double) {
        guard !logLikelihoods.isEmpty, delta.isFinite, delta >= 0,
              logLikelihoods.allSatisfy({ $0.isFinite && abs($0) <= 1e250 }),
              let maximum = logLikelihoods.max() else { throw VivoPosteriorError.invalid("tempering weights") }
        let unnormalized = logLikelihoods.map { exp(delta * ($0 - maximum)) }
        let total = unnormalized.reduce(0, +)
        guard total.isFinite, total > 0 else { throw VivoPosteriorError.numerical("zero tempering weight") }
        let weights = unnormalized.map { $0 / total }
        let squares = weights.reduce(0) { $0 + $1 * $1 }
        guard squares.isFinite, squares > 0 else { throw VivoPosteriorError.numerical("invalid weight concentration") }
        return (weights, 1 / squares)
    }

    static func nextTemperature(logLikelihoods: [Double], beta: Double,
                                configuration: VivoSMCConfiguration) throws -> Double {
        guard beta.isFinite, beta >= 0, beta < 1 else { throw VivoPosteriorError.invalid("tempering clock") }
        let target = configuration.targetESSFraction * Double(logLikelihoods.count)
        if try normalizedWeights(logLikelihoods: logLikelihoods, delta: 1 - beta).1 >= target { return 1 }
        var lower = 0.0, upper = 1 - beta
        for _ in 0..<64 {
            let middle = lower + (upper - lower) * 0.5
            if try normalizedWeights(logLikelihoods: logLikelihoods, delta: middle).1 >= target { lower = middle }
            else { upper = middle }
        }
        let next = beta + lower
        guard next > beta, next - beta >= configuration.minimumTemperatureIncrement else {
            throw VivoPosteriorError.budget("likelihood concentration requires a smaller tempering increment than declared")
        }
        return next
    }

    static func proposalCholesky(particles: [VivoPosteriorParticle], weights: [Double],
                                configuration: VivoSMCConfiguration) throws -> [Double] {
        guard let first = particles.first, weights.count == particles.count else {
            throw VivoPosteriorError.invalid("proposal population")
        }
        let d = first.coordinates.count
        guard (1...32).contains(d), particles.allSatisfy({ $0.coordinates.count == d }),
              weights.allSatisfy({ $0.isFinite && $0 >= 0 }), abs(weights.reduce(0, +) - 1) < 1e-10 else {
            throw VivoPosteriorError.invalid("proposal covariance support")
        }
        var mean = [Double](repeating: 0, count: d)
        for (p, w) in zip(particles, weights) {
            for j in 0..<d { mean[j] += w * p.coordinates[j] }
        }
        var covariance = [Double](repeating: 0, count: d * d)
        for (p, w) in zip(particles, weights) {
            for i in 0..<d {
                for j in 0...i {
                    covariance[i * d + j] += w * (p.coordinates[i] - mean[i]) * (p.coordinates[j] - mean[j])
                }
            }
        }
        let scale = configuration.proposalScale * 2.38 / sqrt(Double(d))
        for i in 0..<d {
            covariance[i * d + i] += configuration.covarianceRegularization
            for j in 0...i {
                covariance[i * d + j] *= scale * scale
                covariance[j * d + i] = covariance[i * d + j]
            }
        }
        var factor = [Double](repeating: 0, count: d * d)
        for i in 0..<d {
            for j in 0...i {
                var value = covariance[i * d + j]
                for k in 0..<j { value -= factor[i * d + k] * factor[j * d + k] }
                if i == j {
                    guard value.isFinite, value > 0 else { throw VivoPosteriorError.numerical("proposal covariance is not positive definite") }
                    factor[i * d + j] = sqrt(value)
                } else {
                    factor[i * d + j] = value / factor[j * d + j]
                }
            }
        }
        guard factor.allSatisfy(\.isFinite) else { throw VivoPosteriorError.numerical("proposal factor overflow") }
        return factor
    }

    static func resample(_ particles: [VivoPosteriorParticle], weights: [Double],
                         random: inout VivoSplitMix64) -> [VivoPosteriorParticle] {
        let n = particles.count, first = openUnit(&random) / Double(particles.count)
        var index = 0, cumulative = weights[0]
        return (0..<n).map { draw in
            let threshold = first + Double(draw) / Double(n)
            while index + 1 < n, threshold >= cumulative {
                index += 1; cumulative += weights[index]
            }
            return particles[index]
        }
    }

    static func quantile(sorted: [Double], probability: Double) -> Double {
        let position = probability * Double(sorted.count - 1)
        let lower = Int(floor(position)), upper = min(lower + 1, sorted.count - 1)
        return sorted[lower] + (position - Double(lower)) * (sorted[upper] - sorted[lower])
    }

    /// Quantile of an equally weighted Gaussian mixture. Numerical bisection,
    /// not random measurement-noise draws. Observation SD is conditional input.
    static func gaussianMixtureQuantile(means: [Double], sd: Double, probability: Double) throws -> Double {
        guard let minimum = means.min(), let maximum = means.max(), means.allSatisfy(\.isFinite),
              sd.isFinite, sd > 0, probability > 0, probability < 1 else {
            throw VivoPosteriorError.invalid("predictive mixture")
        }
        var lower = minimum - 12 * sd, upper = maximum + 12 * sd
        guard lower.isFinite, upper.isFinite, lower < upper else {
            throw VivoPosteriorError.numerical("predictive interval cannot be represented")
        }
        let denominator = Double(means.count)
        for _ in 0..<80 {
            let middle = lower + 0.5 * (upper - lower)
            var cdf = 0.0
            for mean in means {
                cdf += (0.5 * erfc(-(middle - mean) / sd / sqrt(2))) / denominator
            }
            if cdf < probability { lower = middle } else { upper = middle }
        }
        return lower + (upper - lower) * 0.5
    }
}
