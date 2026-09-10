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

/// Conditional count-likelihood MAP for c' beta ~ Normal(0, priorSD^2).
/// Effects and SDs use natural logs. The Laplace SD conditions on the supplied
/// dispersion and prior; it is not a Wald SE or calibrated posterior coverage.
public struct VivoOmicsNBContrastMAPFit: Codable, Sendable, Equatable {
    public let coefficients: [Double]
    public let means: [Double]
    public let effect: Double
    public let posteriorStandardDeviation: Double
    public let priorStandardDeviation: Double
    public let dispersion: Double
    public let logLikelihood: Double
    /// Log likelihood minus quadratic prior penalty, omitting prior constants.
    public let objective: Double
    public let iterations: Int
    public let converged: Bool
    public let maximumScaledScore: Double
}

/// Tests c' beta = 0 with dispersion fixed. Null coefficients use the original
/// coordinates. Failed null optimization retains diagnostics but no probability.
public struct VivoOmicsNBLikelihoodRatioFit: Codable, Sendable, Equatable {
    public let nullCoefficients: [Double]
    public let nullMeans: [Double]
    public let nullLogLikelihood: Double
    public let nullIterations: Int
    public let nullConverged: Bool
    public let nullMaximumScaledScore: Double
    public let rawStatistic: Double
    public let statistic: Double?
    public let pValue: Double?
    public let error: String?
    public let degreesOfFreedom: Int
}

public enum VivoOmicsNegativeBinomial {
    /// Fixed-dispersion NB likelihood ratio; asymptotic chi-square with one DF.
    /// This is not a quasi-likelihood test or dispersion-uncertainty adjustment.
    public static func fitContrastLikelihoodRatio(counts: [UInt64], design: [[Double]], offsets: [Double],
        contrast: [Double], dispersion: Double) throws -> VivoOmicsNBLikelihoodRatioFit {
        let full = try fit(counts: counts,design: design,offsets: offsets,contrast: contrast,dispersion: dispersion)
        return try contrastLikelihoodRatio(counts: counts,design: design,offsets: offsets,contrast: contrast,full: full)
    }
    // Internal reuse only: full must be the unpenalized fit to these exact inputs.
    static func contrastLikelihoodRatio(counts: [UInt64], design: [[Double]], offsets: [Double],
        contrast: [Double], full: VivoOmicsNBFit) throws -> VivoOmicsNBLikelihoodRatioFit {
        guard full.converged, !full.positiveCountDesignRankDeficient,
              let pivot = contrast.indices.max(by: { abs(contrast[$0]) < abs(contrast[$1]) }),
              contrast[pivot] != 0 else {
            throw VivoOmicsStatisticsError.invalid("NB likelihood ratio requires an identified converged full fit and nonzero contrast")
        }
        let free = contrast.indices.filter { $0 != pivot }
        var beta = [Double](repeating: 0,count: contrast.count)
        let mu: [Double], ll: Double, iterations: Int, converged: Bool, scaledScore: Double
        if free.isEmpty {
            guard offsets.allSatisfy({ (-700...700).contains($0) }) else {
                throw VivoOmicsStatisticsError.invalid("NB constrained offsets overflow")
            }
            mu = offsets.map(exp)
            ll = counts.indices.reduce(0) { $0 + mass(Double(counts[$1]),mu[$1],full.dispersion) }
            iterations = 0; converged = true; scaledScore = 0
        } else {
            // Eliminate the largest contrast coefficient. Ratios stay <= 1;
            // no special case assumes that treatment is a particular column.
            let reduced = design.map { row in free.map { row[$0]-row[pivot]*(contrast[$0]/contrast[pivot]) } }
            // The QR fit requires a nonzero vector for its information report.
            // Select a nuisance coefficient solely for that discarded report;
            // the null constraint is already encoded by the reduced design.
            let nuisanceContrast = [1.0] + [Double](repeating: 0,count: free.count-1)
            let null = try fit(counts: counts,design: reduced,offsets: offsets,
                contrast: nuisanceContrast,dispersion: full.dispersion)
            for (j,column) in free.enumerated() { beta[column] = null.coefficients[j] }
            beta[pivot] = -free.reduce(0) { $0 + beta[$1]*(contrast[$1]/contrast[pivot]) }
            mu = null.means; ll = null.logLikelihood; iterations = null.iterations
            converged = null.converged && !null.positiveCountDesignRankDeficient
            scaledScore = null.maximumScaledScore
        }
        // Cancel count-only log-Gamma terms before summing. Subtracting complete
        // log likelihoods loses small LR differences at large counts.
        let a = full.dispersion
        var difference = 0.0, compensation = 0.0
        for i in counts.indices {
            let delta = full.means[i]-mu[i], relative = delta/mu[i]
            let logMeanRatio = abs(relative) < 0.5 ? log1p(relative) : log(full.means[i])-log(mu[i])
            let ratio = a*delta/(1+a*mu[i])
            let logSizeRatio = abs(ratio) < 0.5 ? log1p(ratio) : log1p(a*full.means[i])-log1p(a*mu[i])
            let term = Double(counts[i])*logMeanRatio-(Double(counts[i])+1/a)*logSizeRatio
            let corrected = term-compensation, next = difference+corrected
            compensation = (next-difference)-corrected; difference = next
        }
        let raw = 2*difference
        guard raw.isFinite, ll.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite NB likelihood ratio") }
        let failure: String? = !converged ? "NB likelihood-ratio null fit did not converge" :
            raw < -1e-7 ? "NB constrained likelihood exceeds the full fit beyond numerical tolerance" : nil
        let statistic = failure == nil ? max(0,raw) : nil
        return .init(nullCoefficients: beta,nullMeans: mu,nullLogLikelihood: ll,nullIterations: iterations,
            nullConverged: converged,nullMaximumScaledScore: scaledScore,rawStatistic: raw,
            statistic: statistic,pValue: statistic.map { erfc(sqrt($0/2)) },error: failure,degreesOfFreedom: 1)
    }
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
    /// Shrinks the requested contrast while jointly refitting nuisance coefficients.
    /// A prior never rescues rank-deficient support or an unconverged initial MLE.
    /// The explicit prior SD is in natural-log effect units; no prior is learned.
    public static func fitContrastMAP(counts: [UInt64], design: [[Double]], offsets: [Double],
        contrast: [Double], dispersion: Double, priorStandardDeviation: Double,
        maximumIterations: Int = 100) throws -> VivoOmicsNBContrastMAPFit {
        guard priorStandardDeviation.isFinite, (1e-6...1e6).contains(priorStandardDeviation),
              (1...1000).contains(maximumIterations) else {
            throw VivoOmicsStatisticsError.invalid("NB contrast prior SD or iteration domain")
        }
        let initial = try fit(counts: counts, design: design, offsets: offsets, contrast: contrast, dispersion: dispersion)
        guard initial.converged, !initial.positiveCountDesignRankDeficient else {
            throw VivoOmicsStatisticsError.invalid("NB contrast MAP requires a converged identified unpenalized fit")
        }
        // Spell out numeric conversion: Double.init also has a UInt64
        // bitPattern overload, which is not a count-to-FP64 conversion.
        let n = counts.count, p = contrast.count, y = counts.map { Double($0) }
        let priorRow = contrast.map { $0 / priorStandardDeviation }
        let precision = 1 / (priorStandardDeviation * priorStandardDeviation)
        func effect(_ beta: [Double]) -> Double { zip(beta,contrast).reduce(0) { $0+$1.0*$1.1 } }
        func likelihood(_ mu: [Double]) -> Double { y.indices.reduce(0) { $0+mass(y[$1],mu[$1],dispersion) } }
        func means(_ beta: [Double]) -> [Double]? {
            let eta = y.indices.map { i in offsets[i]+zip(design[i],beta).reduce(0.0) { $0+$1.0*$1.1 } }
            guard eta.allSatisfy({ $0.isFinite && (-700...700).contains($0) }) else { return nil }
            return eta.map(exp)
        }
        func curvature(_ mu: [Double]) -> [Double] {
            y.indices.map { (1+dispersion*y[$0])*mu[$0]/pow(1+dispersion*mu[$0],2) }
        }
        func score(_ beta: [Double], _ mu: [Double]) -> Double {
            let weights = curvature(mu), e = effect(beta)
            return (0..<p).map { j in
                var value = -contrast[j]*e*precision, info = priorRow[j]*priorRow[j]
                for i in 0..<n {
                    value += design[i][j]*(y[i]-mu[i])/(1+dispersion*mu[i])
                    info += design[i][j]*design[i][j]*weights[i]
                }
                return abs(value)/sqrt(info)
            }.max() ?? .infinity
        }
        func information(_ roots: [Double]) throws -> VivoOmicsQR {
            try VivoOmicsQR(design: y.indices.map { i in design[i].map { $0*roots[i] } } + [priorRow],
                            maximumObservations: 513)
        }
        var beta = initial.coefficients, mu = initial.means, ll = initial.logLikelihood
        var objective = ll-pow(effect(beta),2)*precision/2, iterations = 0
        for iteration in 1...maximumIterations {
            try Task.checkCancellation()
            iterations = iteration
            if score(beta,mu) <= 1e-7 { break }
            let weights = curvature(mu), roots = weights.map(sqrt)
            let qr = try information(roots)
            let response = y.indices.map { i in
                (log(mu[i])-offsets[i]+(y[i]-mu[i])/(1+dispersion*mu[i])/weights[i])*roots[i]
            } + [0.0]
            let proposed = try qr.fit(response,contrast: contrast).coefficients
            var fraction = 1.0, accepted = false
            for _ in 0..<30 {
                let candidate = zip(beta,proposed).map { $0+fraction*($1-$0) }
                if let next = means(candidate) {
                    let nextLL = likelihood(next), nextObjective = nextLL-pow(effect(candidate),2)*precision/2
                    if nextObjective.isFinite && nextObjective >= objective-1e-10*max(1,abs(objective)) {
                        beta = candidate; mu = next; ll = nextLL; objective = nextObjective; accepted = true; break
                    }
                }
                fraction *= 0.5
            }
            if !accepted { break }
        }
        let scaledScore = score(beta,mu), qr = try information(curvature(mu).map(sqrt))
        let variance = try qr.fit([Double](repeating: 0,count: n+1),contrast: contrast).contrastVarianceScale
        guard scaledScore.isFinite, objective.isFinite, variance.isFinite, variance > 0 else {
            throw VivoOmicsStatisticsError.invalid("nonfinite NB contrast MAP diagnostics")
        }
        return .init(coefficients: beta,means: mu,effect: effect(beta),posteriorStandardDeviation: sqrt(variance),
            priorStandardDeviation: priorStandardDeviation,dispersion: dispersion,logLikelihood: ll,objective: objective,
            iterations: iterations,converged: scaledScore <= 1e-7,maximumScaledScore: scaledScore)
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
