import Foundation

public enum VivoOmicsStatisticsError: Error, LocalizedError, Sendable {
    case invalid(String)
    public var errorDescription: String? {
        switch self { case .invalid(let message): return "omics statistics: \(message)" }
    }
}

/// FP64 numerical authority for bounded sample-by-covariate linear models.
/// Twice-reorthogonalized QR avoids forming X'X. Rank failures are errors,
/// never silently regularized coefficients or automatic dropped covariates.
public struct VivoOmicsQR: Sendable {
    public let observationCount: Int
    public let coefficientCount: Int
    private let design: [[Double]]
    private let q: [[Double]]
    private let r: [[Double]]

    public init(design: [[Double]], relativeRankTolerance: Double = 1e-10) throws {
        let n = design.count, p = design.first?.count ?? 0
        guard n <= 512, p > 0, p <= 128, n > p,
              relativeRankTolerance.isFinite, relativeRankTolerance > 0,
              relativeRankTolerance < 1,
              design.allSatisfy({ $0.count == p && $0.allSatisfy(\.isFinite) }) else {
            throw VivoOmicsStatisticsError.invalid("design shape, finite entries or residual degrees of freedom")
        }
        var columns: [[Double]] = [], triangular = [[Double]](repeating: [Double](repeating: 0, count: p), count: p)
        for j in 0..<p {
            var v = design.map { $0[j] }
            let scale = sqrt(Self.dot(v, v))
            guard scale.isFinite, scale > 0 else { throw VivoOmicsStatisticsError.invalid("zero or overflowing design column \(j)") }
            for _ in 0..<2 {
                for i in 0..<j {
                    let projection = Self.dot(columns[i], v)
                    triangular[i][j] += projection
                    for row in 0..<n { v[row] -= projection * columns[i][row] }
                }
            }
            let length = sqrt(Self.dot(v, v))
            guard length.isFinite, length > relativeRankTolerance * scale else {
                throw VivoOmicsStatisticsError.invalid("rank-deficient or confounded design at column \(j)")
            }
            triangular[j][j] = length
            columns.append(v.map { $0 / length })
        }
        observationCount = n; coefficientCount = p
        self.design = design; q = columns; r = triangular
    }
    private static func dot(_ a: [Double], _ b: [Double]) -> Double {
        var total = 0.0, correction = 0.0
        for i in a.indices {
            let term = a[i] * b[i] - correction, next = total + term
            correction = (next - total) - term; total = next
        }
        return total
    }
    public func fit(_ response: [Double], contrast: [Double]) throws -> VivoOmicsLinearFit {
        guard response.count == observationCount, contrast.count == coefficientCount,
              response.allSatisfy(\.isFinite), contrast.allSatisfy(\.isFinite) else {
            throw VivoOmicsStatisticsError.invalid("response or contrast shape/value")
        }
        var beta = q.map { Self.dot($0, response) }
        for i in stride(from: coefficientCount - 1, through: 0, by: -1) {
            for j in (i + 1)..<coefficientCount { beta[i] -= r[i][j] * beta[j] }
            beta[i] /= r[i][i]
        }
        // c'(X'X)^(-1)c = ||R^(-T)c||^2, without explicitly inverting X'X.
        var z = contrast
        for i in 0..<coefficientCount {
            for j in 0..<i { z[i] -= r[j][i] * z[j] }
            z[i] /= r[i][i]
        }
        let residuals = (0..<observationCount).map { response[$0] - Self.dot(design[$0], beta) }
        let df = observationCount - coefficientCount
        let variance = Self.dot(residuals, residuals) / Double(df)
        let effect = Self.dot(contrast, beta), scale = Self.dot(z, z)
        guard beta.allSatisfy(\.isFinite), variance.isFinite, variance >= 0,
              effect.isFinite, scale.isFinite, scale > 0 else {
            throw VivoOmicsStatisticsError.invalid("nonfinite fit or null contrast")
        }
        return .init(coefficients: beta, effect: effect, residualVariance: variance,
                     contrastVarianceScale: scale, residualDegreesOfFreedom: df)
    }
}
public struct VivoOmicsLinearFit: Codable, Sendable, Equatable {
    public let coefficients: [Double]
    public let effect: Double
    public let residualVariance: Double
    public let contrastVarianceScale: Double
    public let residualDegreesOfFreedom: Int
}
public struct VivoOmicsVariancePrior: Codable, Sendable, Equatable {
    public let variance: Double
    public let degreesOfFreedom: Double
    public let fittedFeatures: Int
    public let omittedZeroVariances: Int
    public let upperBoundReached: Bool
}

public enum VivoOmicsLinearStatistics {
    public static func median(_ values: [Double]) throws -> Double {
        guard !values.isEmpty, values.allSatisfy(\.isFinite) else {
            throw VivoOmicsStatisticsError.invalid("median requires finite nonempty values")
        }
        let ordered = values.sorted(), n = ordered.count
        return n % 2 == 1 ? ordered[n / 2] : ordered[n / 2 - 1] / 2 + ordered[n / 2] / 2
    }
    private static func digamma(_ input: Double) -> Double {
        var x = input, result = 0.0
        while x < 12 { result -= 1 / x; x += 1 }
        let t = 1 / (x * x)
        return result + log(x) - 0.5 / x - t * (1.0 / 12 - t * (1.0 / 120 - t * (1.0 / 252 - t / 240)))
    }
    private static func trigamma(_ input: Double) -> Double {
        var x = input, result = 0.0
        while x < 12 { result += 1 / (x * x); x += 1 }
        let a = 1 / x, t = a * a
        return result + a + t / 2 + a * t * (1.0 / 6 - t * (1.0 / 30 - t * (1.0 / 42 - t / 30)))
    }
    /// Moment-fitted scaled inverse-chi-square prior for equal residual df.
    /// This is the untrended, non-robust model, not voom or a limma reimplementation.
    public static func variancePrior(_ variances: [Double], residualDF: Int,
                                     minimumFeatures: Int = 20) throws -> VivoOmicsVariancePrior {
        guard residualDF > 0, minimumFeatures >= 3, variances.allSatisfy({ $0.isFinite && $0 >= 0 }) else {
            throw VivoOmicsStatisticsError.invalid("variance-prior inputs")
        }
        let positive = variances.filter { $0 > 0 }
        guard positive.count >= minimumFeatures else {
            throw VivoOmicsStatisticsError.invalid("empirical Bayes requires at least \(minimumFeatures) positive residual variances")
        }
        var mean = 0.0, m2 = 0.0
        for (i, value) in positive.enumerated() {
            let x = log(value), delta = x - mean
            mean += delta / Double(i + 1); m2 += delta * (x - mean)
        }
        let excess = m2 / Double(positive.count - 1) - trigamma(Double(residualDF) / 2)
        let maximumDF = 1_000_000.0
        var priorDF = maximumDF
        if excess > trigamma(maximumDF / 2) {
            var low = 1e-6, high = maximumDF / 2
            for _ in 0..<100 {
                let middle = (low + high) / 2
                if trigamma(middle) > excess { low = middle } else { high = middle }
            }
            priorDF = low + high
        }
        let d = Double(residualDF) / 2, d0 = priorDF / 2
        let scale = exp(mean - digamma(d) + log(d) + digamma(d0) - log(d0))
        guard scale.isFinite, scale > 0, priorDF.isFinite, priorDF > 0 else {
            throw VivoOmicsStatisticsError.invalid("variance-prior fit did not produce a finite positive scale")
        }
        return .init(variance: scale, degreesOfFreedom: priorDF, fittedFeatures: positive.count,
                     omittedZeroVariances: variances.count - positive.count, upperBoundReached: priorDF == maximumDF)
    }
    private static func betaFraction(_ a: Double, _ b: Double, _ x: Double) throws -> Double {
        let tiny = 1e-300
        func nonzero(_ value: Double) -> Double { abs(value) < tiny ? (value < 0 ? -tiny : tiny) : value }
        var c = 1.0, d = 1 / nonzero(1 - (a + b) * x / (a + 1)), h = d
        for index in 1...512 {
            let m = Double(index), twice = 2 * m
            var aa = m * (b - m) * x / ((a - 1 + twice) * (a + twice))
            d = 1 / nonzero(1 + aa * d); c = nonzero(1 + aa / c); h *= d * c
            aa = -(a + m) * (a + b + m) * x / ((a + twice) * (a + 1 + twice))
            d = 1 / nonzero(1 + aa * d); c = nonzero(1 + aa / c)
            let change = d * c; h *= change
            if abs(change - 1) < 2e-14 { return h }
        }
        throw VivoOmicsStatisticsError.invalid("incomplete beta continued fraction did not converge")
    }
    public static func studentTwoSidedP(t: Double, degreesOfFreedom df: Double) throws -> Double {
        guard t.isFinite, df.isFinite, df > 0 else { throw VivoOmicsStatisticsError.invalid("Student t inputs") }
        if t == 0 { return 1 }
        let x = 1 / (1 + (abs(t) / sqrt(df)) * (abs(t) / sqrt(df)))
        if x == 0 { return 0 }
        if x == 1 { return 1 }
        let a = df / 2, b = 0.5
        let prefactor = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log1p(-x))
        let value: Double
        if x < (a + 1) / (a + b + 2) { value = try prefactor * betaFraction(a, b, x) / a }
        else { value = try 1 - prefactor * betaFraction(b, a, 1 - x) / b }
        guard value.isFinite, value >= -1e-12, value <= 1 + 1e-12 else {
            throw VivoOmicsStatisticsError.invalid("invalid Student t probability")
        }
        return min(1, max(0, value))
    }
    public static func studentCriticalValue(degreesOfFreedom df: Double, coverage: Double = 0.95) throws -> Double {
        guard coverage.isFinite, coverage > 0, coverage < 1 else { throw VivoOmicsStatisticsError.invalid("interval coverage") }
        let alpha = 1 - coverage
        var low = 0.0, high = 1.0
        while try studentTwoSidedP(t: high, degreesOfFreedom: df) > alpha {
            high *= 2
            guard high <= 1e12 else { throw VivoOmicsStatisticsError.invalid("Student quantile exceeds bounded search") }
        }
        for _ in 0..<80 {
            let middle = (low + high) / 2
            if try studentTwoSidedP(t: middle, degreesOfFreedom: df) > alpha { low = middle } else { high = middle }
        }
        return (low + high) / 2
    }
    /// Family consists of the supplied tested hypotheses only; nil/filtered
    /// hypotheses must not be assigned fabricated p-values by the caller.
    public static func benjaminiHochberg(_ probabilities: [Double]) throws -> [Double] {
        guard probabilities.allSatisfy({ $0.isFinite && (0...1).contains($0) }) else {
            throw VivoOmicsStatisticsError.invalid("FDR inputs")
        }
        if probabilities.isEmpty { return [] }
        let order = probabilities.indices.sorted {
            probabilities[$0] == probabilities[$1] ? $0 < $1 : probabilities[$0] < probabilities[$1]
        }
        var result = [Double](repeating: 1, count: probabilities.count), running = 1.0
        for i in stride(from: order.count - 1, through: 0, by: -1) {
            running = min(running, probabilities[order[i]] * Double(order.count) / Double(i + 1))
            result[order[i]] = running
        }
        return result
    }
}
