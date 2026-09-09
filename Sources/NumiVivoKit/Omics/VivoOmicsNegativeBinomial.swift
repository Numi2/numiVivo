import Foundation

/// A conditional NB2 GLM fit: Var(Y) = mu + dispersion * mu^2.
/// Coefficients/effects use natural logs. No Wald p-value is reported until
/// dispersion estimation and its cohort-level diagnostics are qualified.
public struct VivoOmicsNBFit: Codable, Sendable, Equatable {
    public let coefficients: [Double]
    public let means: [Double]
    public let effect: Double?
    public let standardError: Double?
    public let logLikelihood: Double
    public let coxReidLogLikelihood: Double?
    public let positiveCountDesignRankDeficient: Bool
    public let dispersion: Double
    public let iterations: Int
    public let converged: Bool
    public let maximumScaledScore: Double
    public let pearsonResiduals: [Double]
    public let leverage: [Double]
    /// Unavailable on deficient support or unit-leverage observations.
    public let cooksDistances: [Double]?
}
public struct VivoOmicsNBDispersionFit: Codable, Sendable, Equatable {
    public let fit: VivoOmicsNBFit
    public let objective: Double
    public let lowerBoundary: Bool
    public let upperBoundary: Bool
    public let profileEvaluations: Int
    public let logPriorMean: Double?
    public let logPriorVariance: Double?
}

public enum VivoOmicsNegativeBinomial {
    /// Stable NB log mass, including the near-Poisson small-dispersion limit.
    public static func logMass(count: UInt64, mean: Double, dispersion: Double) throws -> Double {
        guard count <= 9_007_199_254_740_992, mean.isFinite, mean > 0,
              dispersion.isFinite, (1e-8...100).contains(dispersion) else {
            throw VivoOmicsStatisticsError.invalid("NB count, mean or dispersion domain")
        }
        let value = mass(Double(count), mean, dispersion)
        guard value.isFinite else { throw VivoOmicsStatisticsError.invalid("NB log mass overflow") }
        return value
    }
    private static func mass(_ y: Double, _ mu: Double, _ a: Double) -> Double {
        let r = 1 / a
        let ratio: Double
        if r >= 8 {
            // Stirling difference avoids cancellation of two huge log-Gammas.
            func correction(_ z: Double) -> Double {
                let v = 1 / z, v2 = v * v
                return v * (1/12.0 + v2 * (-1/360.0 + v2 * (1/1260.0 + v2 * (-1/1680.0 + v2/1188.0))))
            }
            ratio = (r + y - 0.5) * log1p(y/r) - y + correction(r+y) - correction(r)
        } else { ratio = lgamma(r+y) - lgamma(r) - y * log(r) }
        return ratio - lgamma(y+1) + y * log(mu) - (y+r) * log1p(a*mu)
    }
    /// Call only after response/design dimensions have been validated.
    static func positiveSupportIsRankDeficient(counts: [UInt64],design: [[Double]]) -> Bool {
        let p = design[0].count
        let support = counts.indices.filter { counts[$0] > 0 }.map { design[$0] }
        let rows = support.count == p ? support + [support[0]] : support
        return support.count < p || (try? VivoOmicsQR(design: rows)) == nil
    }
    public static func fit(counts: [UInt64], design: [[Double]], offsets: [Double],
                           contrast: [Double], dispersion: Double, maximumIterations: Int = 100) throws -> VivoOmicsNBFit {
        let base = try VivoOmicsQR(design: design)
        let n = base.observationCount, p = base.coefficientCount
        guard counts.count == n, offsets.count == n, contrast.count == p,
              offsets.allSatisfy(\.isFinite), contrast.allSatisfy(\.isFinite),
              counts.allSatisfy({ $0 <= 9_007_199_254_740_992 }), counts.contains(where: { $0 > 0 }),
              dispersion.isFinite, (1e-8...100).contains(dispersion), (1...1000).contains(maximumIterations) else {
            throw VivoOmicsStatisticsError.invalid("NB response, offsets, contrast or iteration domain")
        }
        let deficient = positiveSupportIsRankDeficient(counts: counts,design: design)
        let y = counts.map { Double($0) }
        // Pseudocount is initialization only; the fitted likelihood uses raw y.
        var beta = try base.fit(y.indices.map { log(y[$0]+0.5)-offsets[$0] }, contrast: contrast).coefficients
        func means(_ coefficients: [Double]) -> [Double]? {
            let eta = (0..<n).map { i in offsets[i] + zip(design[i],coefficients).reduce(0.0) { $0 + $1.0 * $1.1 } }
            guard eta.allSatisfy({ $0.isFinite && (-700...700).contains($0) }) else { return nil }
            return eta.map(exp)
        }
        guard var mu = means(beta) else { throw VivoOmicsStatisticsError.invalid("NB initialization overflows") }
        func likelihood(_ m: [Double]) -> Double { (0..<n).reduce(0) { $0 + mass(y[$1],m[$1],dispersion) } }
        func weighted(_ m: [Double]) throws -> (VivoOmicsQR, [Double]) {
            let roots = m.map { sqrt($0 / (1 + dispersion * $0)) }
            return (try VivoOmicsQR(design: (0..<n).map { i in design[i].map { $0 * roots[i] } }), roots)
        }
        func score(_ m: [Double]) -> Double {
            (0..<p).map { j in
                var value = 0.0, info = 0.0
                for i in 0..<n {
                    value += design[i][j] * (y[i]-m[i]) / (1+dispersion*m[i])
                    info += design[i][j]*design[i][j]*m[i]/(1+dispersion*m[i])
                }
                return abs(value) / sqrt(info)
            }.max() ?? .infinity
        }
        var ll = likelihood(mu), iterations = 0, converged = false
        guard ll.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite NB likelihood") }
        for iteration in 1...maximumIterations {
            try Task.checkCancellation()
            iterations = iteration
            if score(mu) <= 1e-7 { converged = true; break }
            // Observed-information Newton steps converge at high dispersion
            // where Fisher scoring can stall. Final covariance still uses
            // expected Fisher information, not this optimization curvature.
            let curvature = (0..<n).map { i in (1+dispersion*y[i])*mu[i]/pow(1+dispersion*mu[i],2) }
            let roots = curvature.map(sqrt)
            let qr = try VivoOmicsQR(design: (0..<n).map { i in design[i].map { $0*roots[i] } })
            let response = (0..<n).map { i in
                let score = (y[i]-mu[i])/(1+dispersion*mu[i])
                return (log(mu[i])-offsets[i]+score/curvature[i])*roots[i]
            }
            let proposed = try qr.fit(response, contrast: contrast).coefficients
            var fraction = 1.0, accepted = false
            for _ in 0..<30 {
                let candidate = zip(beta,proposed).map { $0 + fraction * ($1-$0) }
                if let next = means(candidate) {
                    let nextLL = likelihood(next)
                    if nextLL.isFinite && nextLL >= ll - 1e-10 * max(1,abs(ll)) {
                        beta = candidate; mu = next; ll = nextLL; accepted = true; break
                    }
                }
                fraction *= 0.5
            }
            if !accepted { break }
        }
        let scaledScore = score(mu)
        converged = converged || scaledScore <= 1e-7
        let (qr,roots) = try weighted(mu)
        let covariance = try qr.fit([Double](repeating: 0,count: n),contrast: contrast).contrastVarianceScale
        let h = qr.leverage
        let residuals = (0..<n).map { (y[$0]-mu[$0]) * roots[$0] / mu[$0] }
        let rawCooks = (0..<n).map { residuals[$0]*residuals[$0]*h[$0] / (Double(p)*pow(1-h[$0],2)) }
        let cooks = !deficient && h.allSatisfy({ $0 < 1 }) && rawCooks.allSatisfy(\.isFinite) ? rawCooks : nil
        guard covariance.isFinite, scaledScore.isFinite else {
            throw VivoOmicsStatisticsError.invalid("nonfinite NB information or influence diagnostics")
        }
        return .init(coefficients: beta, means: mu, effect: deficient ? nil : zip(beta,contrast).reduce(0) { $0+$1.0*$1.1 },
            standardError: deficient ? nil : sqrt(covariance), logLikelihood: ll,
            coxReidLogLikelihood: deficient ? nil : ll-0.5*qr.logInformationDeterminant,
            positiveCountDesignRankDeficient: deficient, dispersion: dispersion,
            iterations: iterations, converged: converged, maximumScaledScore: scaledScore,
            pearsonResiduals: residuals, leverage: h, cooksDistances: cooks)
    }
    /// Bounded gene-wise adjusted-profile or supplied log-normal MAP dispersion.
    /// The prior is explicit; this routine does not learn a cohort dispersion trend.
    public static func estimateDispersion(counts: [UInt64], design: [[Double]], offsets: [Double],
        contrast: [Double], lower: Double = 1e-8, upper: Double = 100,
        logPriorMean: Double? = nil, logPriorVariance: Double? = nil) throws -> VivoOmicsNBDispersionFit {
        guard lower.isFinite, upper.isFinite, lower >= 1e-8, upper <= 100, lower < upper,
              (logPriorMean == nil) == (logPriorVariance == nil),
              logPriorMean.map(\.isFinite) ?? true,
              logPriorVariance.map({ $0.isFinite && $0 > 0 }) ?? true else {
            throw VivoOmicsStatisticsError.invalid("NB profile bounds or prior")
        }
        var evaluations = 0
        func evaluate(_ t: Double) throws -> (VivoOmicsNBFit, Double) {
            let fit = try fit(counts: counts, design: design, offsets: offsets, contrast: contrast, dispersion: min(upper,max(lower,exp(t))))
            guard fit.converged else { throw VivoOmicsStatisticsError.invalid("NB dispersion profile contains an unconverged coefficient fit at alpha=\(fit.dispersion), scaledScore=\(fit.maximumScaledScore), iterations=\(fit.iterations)") }
            guard let adjusted = fit.coxReidLogLikelihood else {
                throw VivoOmicsStatisticsError.invalid("NB positive-count support is rank deficient; dispersion profile is not identified")
            }
            evaluations += 1
            let penalty = logPriorMean.map { pow(t-$0,2)/(2*logPriorVariance!) } ?? 0
            return (fit,adjusted-penalty)
        }
        // Coarse global bracket, followed by golden-section refinement. Boundary
        // estimates are retained and labeled, never interpreted as interior optima.
        let lo = log(lower), hi = log(upper), grid = (0...24).map { lo+(hi-lo)*Double($0)/24 }
        var samples: [(VivoOmicsNBFit,Double)] = []
        for t in grid { samples.append(try evaluate(t)) }
        let best = samples.indices.max { samples[$0].1 < samples[$1].1 }!
        var result = samples[best]
        var left = grid[max(0,best-1)], right = grid[min(24,best+1)]
        let ratio = (sqrt(5.0)-1)/2
        var x = right-ratio*(right-left), z = left+ratio*(right-left)
        var fx = try evaluate(x), fz = try evaluate(z)
        while right-left > 1e-6 {
            if fx.1 > result.1 { result = fx }; if fz.1 > result.1 { result = fz }
            if fx.1 > fz.1 {
                right = z; z = x; fz = fx; x = right-ratio*(right-left); fx = try evaluate(x)
            } else {
                left = x; x = z; fx = fz; z = left+ratio*(right-left); fz = try evaluate(z)
            }
        }
        if fx.1 > result.1 { result = fx }; if fz.1 > result.1 { result = fz }
        // Values within 0.01% on the log scale are conservatively labeled
        // boundary estimates: likelihood rounding can shift a near-Poisson
        // optimum by more than the 1e-6 scalar-search stopping width.
        return .init(fit: result.0, objective: result.1,
            lowerBoundary: abs(log(result.0.dispersion)-lo) < 1e-4,
            upperBoundary: abs(log(result.0.dispersion)-hi) < 1e-4,
            profileEvaluations: evaluations, logPriorMean: logPriorMean, logPriorVariance: logPriorVariance)
    }
}
