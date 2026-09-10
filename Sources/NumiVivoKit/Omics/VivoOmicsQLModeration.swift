import Foundation

public struct VivoOmicsQLProfilePoint: Codable, Sendable, Equatable {
    public let parameter: Double
    public let shape: Double
    public let objective: Double
}
public struct VivoOmicsQLPriorProfile: Codable, Sendable, Equatable {
    public let shape: Double
    public let objective: Double
    public let lowerBoundary: Bool
    public let upperBoundary: Bool
    public let informativeWeight: Double
    public let priorWeights: [Double]
    public let logVarianceTrend: [Double]
    public let smoother: VivoOmicsPrecisionLowessFit?
    public let evaluations: [VivoOmicsQLProfilePoint]
}
public struct VivoOmicsQLModeratedVariances: Codable, Sendable, Equatable {
    public let priorScales: [Double]
    public let commonPriorDegreesOfFreedom: Double
    public let priorDegreesOfFreedom: [Double]
    public let posteriorVariances: [Double]
    public let informativeIndices: [Int]
    public let excludedLowDFIndices: [Int]
    public let flooredVarianceIndices: [Int]
    public let profiles: [VivoOmicsQLPriorProfile]
    public let screeningRightProbabilities: [Double]?
    public let screeningFDRWeights: [Double]?
    public let notOutlierProbabilities: [Double]?
    public let outlierDegreesOfFreedom: Double?
    public let featureEvaluations: Int
}

/// Unequal-DF scaled-F prior and robust posterior variance moderation.
/// This is not a cohort hypothesis test or a claim of calibrated error rates.
public enum VivoOmicsQLModeration {
    static func quantile(_ values: [Double], probability: Double) -> Double {
        let sorted = values.sorted(), position = probability*Double(values.count-1)
        let lower = Int(floor(position)), fraction = position-Double(lower)
        return (1-fraction)*sorted[lower]+fraction*sorted[min(lower+1,sorted.count-1)]
    }
    public static func fit(variances: [Double], degreesOfFreedom: [Double], abundance: [Double]? = nil,
        robust: Bool = true, maximumProfileEvaluations: Int = 128,
        maximumFeatureEvaluations: Int = 100_000_000,
        maximumNeighborhoodVisits: Int = 100_000_000) throws -> VivoOmicsQLModeratedVariances {
        let n = variances.count
        guard (3...100_000).contains(n), degreesOfFreedom.count == n,
              variances.allSatisfy({ $0.isFinite && (0...1e100).contains($0) }),
              degreesOfFreedom.allSatisfy({ $0.isFinite && (0...1e6).contains($0) }),
              abundance == nil || (abundance!.count == n && abundance!.allSatisfy(\.isFinite)),
              (4...1024).contains(maximumProfileEvaluations), maximumFeatureEvaluations > 0,
              maximumNeighborhoodVisits > 0 else {
            throw VivoOmicsStatisticsError.invalid("QL moderation family, variance, DF, abundance or work bounds")
        }
        let excluded = variances.indices.filter { degreesOfFreedom[$0] < 0.01 }
        let informative = variances.indices.filter { variances[$0] > 0 && degreesOfFreedom[$0] >= 0.01 }
        guard informative.count >= 3 else { throw VivoOmicsStatisticsError.invalid("QL prior requires three positive informative variances") }
        let median = quantile(informative.map { variances[$0] },probability: 0.5), logFloor = log(median)+log(1e-12)
        let logs = variances.map { $0 > 0 ? max(logFloor,log($0)) : logFloor }
        let floored = variances.indices.filter { variances[$0] == 0 || log(variances[$0]) < logFloor }
        // Dummy positive DF only evaluates ignored rows; their likelihood weights
        // remain zero in both profiles and the posterior uses the original DF.
        let shape1 = degreesOfFreedom.map { ($0 < 0.01 ? 1 : $0)/2 }
        let mask = degreesOfFreedom.map { $0 < 0.01 ? 0.0 : 1.0 }
        let corrected = try zip(logs,shape1).map { $0+(try VivoOmicsQLSpecialFunctions.logMinusDigamma($1)) }
        let precision = try shape1.map { 1/(try VivoOmicsQLSpecialFunctions.trigamma($0)) }
        let span = min(1,0.3+0.7*pow(500/Double(n),1.0/3))
        var featureEvaluations = 0
        func profile(_ weights: [Double]) throws -> VivoOmicsQLPriorProfile {
            let totalWeight = weights.reduce(0,+)
            guard totalWeight > 0 else { throw VivoOmicsStatisticsError.invalid("QL robust profile has no informative weight") }
            let w = zip(precision,weights).map(*)
            let trend: [Double], smoother: VivoOmicsPrecisionLowessFit?
            if let abundance {
                let normalization = quantile(w,probability: 0.75)
                guard normalization > 0 else { throw VivoOmicsStatisticsError.invalid("QL trend precision upper quartile is zero") }
                let fit = try VivoOmicsPrecisionLowess.fit(x: abundance,y: corrected,weights: w.map { $0/normalization },
                    span: span,maximumNeighborhoodVisits: maximumNeighborhoodVisits)
                smoother = fit; trend = fit.fitted
            } else {
                let total = w.reduce(0,+), mean = zip(w,corrected).reduce(0.0) { $0+$1.0*$1.1/total }
                trend = [Double](repeating: mean,count: n); smoother = nil
            }
            var points: [VivoOmicsQLProfilePoint] = []
            func evaluate(_ parameter: Double) throws -> VivoOmicsQLProfilePoint {
                try Task.checkCancellation()
                guard points.count < maximumProfileEvaluations, n <= maximumFeatureEvaluations-featureEvaluations else {
                    throw VivoOmicsStatisticsError.invalid("QL prior exhausted profile or feature evaluations")
                }
                featureEvaluations += n
                let shape = parameter/(1-parameter), h = try VivoOmicsQLSpecialFunctions.logMinusDigamma(shape)
                let logShape = log(shape)
                var sum = 0.0, correction = 0.0
                for i in 0..<n where weights[i] > 0 {
                    let a = shape1[i], logScaledPrior = logShape+trend[i]-h
                    let logRatio = log(a)+logs[i]-logScaledPrior
                    let gamma = try VivoOmicsQLSpecialFunctions.logGammaIncrement(base: shape,increment: a)
                    let term = -2*weights[i]*(gamma-a*logScaledPrior-(a+shape)*VivoOmicsQLSpecialFunctions.softplus(logRatio))
                    let adjusted = term-correction, next = sum+adjusted
                    correction = (next-sum)-adjusted; sum = next
                }
                guard sum.isFinite else { throw VivoOmicsStatisticsError.invalid("nonfinite QL prior profile objective") }
                let point = VivoOmicsQLProfilePoint(parameter: parameter,shape: shape,objective: sum)
                points.append(point); return point
            }
            let lowerEndpoint = try evaluate(0.5), upperEndpoint = try evaluate(0.9998)
            var low = 0.5, high = 0.9998
            let ratio = (sqrt(5)-1)/2
            var left = high-ratio*(high-low), right = low+ratio*(high-low)
            var first = try evaluate(left), second = try evaluate(right)
            while high-low > 1e-10 {
                if first.objective <= second.objective {
                    high = right; right = left; second = first
                    left = high-ratio*(high-low); first = try evaluate(left)
                } else {
                    low = left; left = right; first = second
                    right = low+ratio*(high-low); second = try evaluate(right)
                }
            }
            _ = try evaluate(low+(high-low)/2)
            let best = points.min { $0.objective < $1.objective }!
            return .init(shape: best.shape,objective: best.objective,
                lowerBoundary: best.parameter == lowerEndpoint.parameter,upperBoundary: best.parameter == upperEndpoint.parameter,
                informativeWeight: totalWeight,priorWeights: weights,logVarianceTrend: trend,smoother: smoother,evaluations: points)
        }
        func scales(_ profile: VivoOmicsQLPriorProfile) throws -> [Double] {
            let h = try VivoOmicsQLSpecialFunctions.logMinusDigamma(profile.shape)
            let result = profile.logVarianceTrend.map { exp($0-h) }
            guard result.allSatisfy({ $0.isFinite && $0 > 0 }) else {
                throw VivoOmicsStatisticsError.invalid("QL prior scales exceed positive finite arithmetic")
            }
            return result
        }
        var profiles = [try profile(mask)], priorScales = try scales(profiles[0])
        var commonDF = 2*profiles[0].shape, priorDF = [Double](repeating: commonDF,count: n)
        var rightProbabilities: [Double]?, screeningWeights: [Double]?, notOutlier: [Double]?, outlierDF: Double?
        if robust {
            var logRight = [Double](repeating: 0,count: n), twoSided = logRight
            let logF = variances.indices.map { variances[$0] > 0 ? log(variances[$0])-log(priorScales[$0]) : -.infinity }
            for i in 0..<n {
                try Task.checkCancellation()
                if variances[i] == 0 { logRight[i] = 0; twoSided[i] = 0 }
                else {
                    let tails = try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: logF[i],numeratorDF: 2*shape1[i],denominatorDF: commonDF)
                    logRight[i] = tails.upper; twoSided[i] = min(1,2*exp(min(tails.lower,tails.upper)))
                }
            }
            let right = logRight.map(exp)
            let weights = try VivoOmicsLinearStatistics.benjaminiHochberg(twoSided).map { $0 > 0.3 ? 1 : $0 }
            rightProbabilities = right; screeningWeights = weights
            if weights.contains(where: { $0 < 1 }) {
                profiles.append(try profile(zip(mask,weights).map(*)))
                priorScales = try scales(profiles[1]); commonDF = 2*profiles[1].shape
                priorDF = [Double](repeating: commonDF,count: n)
                let order = logF.indices.sorted { logF[$0] == logF[$1] ? $0 < $1 : logF[$0] < logF[$1] }
                var ranks = [Double](repeating: 0,count: n), start = 0
                while start < n {
                    var end = start+1
                    while end < n, logF[order[end]] == logF[order[start]] { end += 1 }
                    let rank = Double(start+1+end)/2
                    for j in start..<end { ranks[order[j]] = rank }
                    start = end
                }
                let probability = right.indices.map { min(1,right[$0]/((Double(n)-ranks[$0]+0.5)/Double(n))) }
                notOutlier = probability
                if probability.contains(where: { $0 < 1 }) {
                    let worst = right.indices.min { right[$0] < right[$1] }!
                    let outlier: Double
                    if right[worst] == 0 { outlier = 0 }
                    else {
                        let first = log(0.5)/logRight[worst]*commonDF
                        let tail = try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: logF[worst],
                            numeratorDF: 2*shape1[worst],denominatorDF: first).upper
                        outlier = log(0.5)/tail*first
                    }
                    guard outlier.isFinite, outlier >= 0 else { throw VivoOmicsStatisticsError.invalid("invalid QL outlier prior DF") }
                    outlierDF = outlier
                    priorDF = probability.map { $0*commonDF+(1-$0)*outlier }
                    let sorted = right.indices.sorted { right[$0] == right[$1] ? $0 < $1 : right[$0] < right[$1] }
                    var prefix = 0.0, minimum = Double.infinity, minimumIndex = 0
                    for j in 0..<n {
                        prefix += priorDF[sorted[j]]
                        let mean = prefix/Double(j+1)
                        if mean < minimum { minimum = mean; minimumIndex = j }
                    }
                    for j in 0...minimumIndex { priorDF[sorted[j]] = minimum }
                    var previous = 0.0
                    for j in sorted { previous = max(previous,priorDF[j]); priorDF[j] = previous }
                }
            }
        }
        var posterior = [Double](repeating: 0,count: n)
        for i in 0..<n {
            let denominator = degreesOfFreedom[i]+priorDF[i]
            guard denominator > 0 else { throw VivoOmicsStatisticsError.invalid("QL posterior has neither residual nor prior information") }
            posterior[i] = degreesOfFreedom[i]/denominator*variances[i]+priorDF[i]/denominator*priorScales[i]
        }
        guard posterior.allSatisfy({ $0.isFinite && $0 >= 0 }) else { throw VivoOmicsStatisticsError.invalid("nonfinite QL posterior variance") }
        return .init(priorScales: priorScales,commonPriorDegreesOfFreedom: commonDF,priorDegreesOfFreedom: priorDF,
            posteriorVariances: posterior,informativeIndices: informative,excludedLowDFIndices: excluded,
            flooredVarianceIndices: floored,profiles: profiles,screeningRightProbabilities: rightProbabilities,
            screeningFDRWeights: screeningWeights,notOutlierProbabilities: notOutlier,
            outlierDegreesOfFreedom: outlierDF,featureEvaluations: featureEvaluations)
    }
    public static func fit(from native: VivoOmicsNBQLNativeAbundanceFit, robust: Bool = true) throws -> VivoOmicsQLModeratedVariances {
        guard native.completed, native.failures.isEmpty, let global = native.globalFit, global.completed,
              global.adjustedResiduals.count == native.abundanceFits.count,
              global.adjustedResiduals.allSatisfy({ $0?.quasiDispersion != nil }),
              native.abundanceFits.allSatisfy({ $0 != nil }) else {
            throw VivoOmicsStatisticsError.invalid("QL moderation requires complete identified native abundance and residual inputs")
        }
        return try fit(variances: global.adjustedResiduals.map { $0!.quasiDispersion! },
            degreesOfFreedom: global.adjustedResiduals.map { $0!.degreesOfFreedom },
            abundance: native.abundanceFits.map { $0!.log2CountsPerMillion },robust: robust)
    }
}
