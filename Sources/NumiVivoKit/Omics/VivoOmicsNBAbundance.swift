import Foundation

public struct VivoOmicsNBAbundanceFit: Codable, Sendable, Equatable {
    public let log2CountsPerMillion: Double
    public let logProportion: Double
    public let priorCount: Double
    public let dispersion: Double
    public let scoreEvaluations: Int
    public let logProportionBracketWidth: Double
    public let scaledScore: Double
}

/// Intercept-only NB abundance with library-scaled prior counts.
/// The augmented counts are used only for this descriptive covariate.
public enum VivoOmicsNBAbundance {
    public static func fit(counts: [UInt64], offsets: [Double], dispersion: Double,
        priorCount: Double = 2, maximumScoreEvaluations: Int = 128) throws -> VivoOmicsNBAbundanceFit {
        let n = counts.count
        guard (1...512).contains(n), offsets.count == n,
              counts.allSatisfy({ $0 <= 9_007_199_254_740_992 }),
              offsets.allSatisfy({ $0.isFinite && (-700...700).contains($0) }),
              dispersion.isFinite, dispersion == 0 || (1e-8...100).contains(dispersion),
              priorCount.isFinite, (0...1e8).contains(priorCount),
              (1...1024).contains(maximumScoreEvaluations) else {
            throw VivoOmicsStatisticsError.invalid("NB abundance counts, offsets, dispersion, prior or work bound")
        }
        let center = offsets.max()!
        let relativeLibraries = offsets.map { exp($0-center) }
        let meanRelativeLibrary = relativeLibraries.reduce(0,+)/Double(n)
        let logMeanLibrary = center+log(meanRelativeLibrary)
        let logAdjustment: Double
        if priorCount == 0 { logAdjustment = 0 }
        else {
            let t = log(2*priorCount)-logMeanLibrary
            logAdjustment = max(0,t)+log1p(exp(-abs(t)))
        }
        let centeredOffsets = offsets.map { $0-center }
        var augmented = [Double](); augmented.reserveCapacity(n)
        for i in counts.indices {
            let prior = priorCount*(relativeLibraries[i]/meanRelativeLibrary)
            guard priorCount == 0 || prior > 0 else {
                throw VivoOmicsStatisticsError.invalid("NB abundance scaled prior underflow")
            }
            augmented.append(Double(counts[i])+prior)
        }
        let total = augmented.reduce(0,+)
        guard total > 0, total.isFinite else {
            throw VivoOmicsStatisticsError.invalid("NB abundance has no finite log proportion at zero augmented count")
        }
        // The common prior offset adjustment cancels from the centered geometry.
        // Restore it only when translating the fitted centered intercept.
        let initial = log(total)-log(relativeLibraries.reduce(0,+))
        var evaluations = 0
        func evaluate(_ coefficient: Double) throws -> (score: Double, information: Double) {
            try Task.checkCancellation()
            guard evaluations < maximumScoreEvaluations else {
                throw VivoOmicsStatisticsError.invalid("NB abundance exhausted \(maximumScoreEvaluations) score evaluations")
            }
            evaluations += 1
            var score = 0.0, information = 0.0
            for i in counts.indices {
                let eta = coefficient+centeredOffsets[i]
                if dispersion == 0 {
                    let mean = exp(eta)
                    score += augmented[i]-mean; information += mean
                } else {
                    let t = log(dispersion)+eta
                    let small = exp(-abs(t)), denominator = 1+small
                    let inverse = t >= 0 ? small/denominator : 1/denominator
                    let saturated = t >= 0 ? 1/denominator : small/denominator
                    let weight = saturated/dispersion
                    score += augmented[i]*inverse-weight; information += weight
                }
            }
            guard !score.isNaN, !information.isNaN else {
                throw VivoOmicsStatisticsError.invalid("nonfinite NB abundance score geometry")
            }
            return (score,information)
        }
        var low = initial, high = initial
        let start = try evaluate(initial)
        var answer = initial, answerScore = start, width = 0.0
        if start.score != 0 {
            var step = 1.0
            if start.score > 0 {
                repeat { high = initial+step; step *= 2 } while try evaluate(high).score > 0
            } else {
                repeat { low = initial-step; step *= 2 } while try evaluate(low).score < 0
            }
            while high-low > 2e-12 {
                let middle = low+(high-low)/2
                guard middle > low, middle < high else { break }
                let value = try evaluate(middle)
                if value.score == 0 { low = middle; high = middle; break }
                if value.score > 0 { low = middle } else { high = middle }
            }
            answer = low+(high-low)/2; width = high-low
            answerScore = try evaluate(answer)
        }
        let logProportion = (answer-center)-logAdjustment
        let logCPM = (logProportion+log(1_000_000))/log(2)
        let scaledScore = abs(answerScore.score)/sqrt(answerScore.information)
        guard logCPM.isFinite, logProportion.isFinite, scaledScore.isFinite else {
            throw VivoOmicsStatisticsError.invalid("NB abundance output exceeds finite arithmetic")
        }
        return .init(log2CountsPerMillion: logCPM,logProportion: logProportion,
            priorCount: priorCount,dispersion: dispersion,scoreEvaluations: evaluations,
            logProportionBracketWidth: width,scaledScore: scaledScore)
    }
}

public struct VivoOmicsNBQLNativeAbundanceFit: Codable, Sendable, Equatable {
    public let completed: Bool
    public let abundanceFits: [VivoOmicsNBAbundanceFit?]
    public let globalFit: VivoOmicsNBQLGlobalFit?
    public let failures: [VivoOmicsNBQLFailure]
}

extension VivoOmicsNBQLGlobalScale {
    /// Estimate native abundance before the existing two-update global scale fit.
    /// Original counts enter the GLMs unchanged; prior counts affect abundance only.
    public static func fitEstimatingAbundance(counts: [[UInt64]], design: [[Double]], offsets: [Double],
        contrast: [Double], trendDispersions: [Double], abundancePriorCount: Double = 2,
        maximumAbundanceScoreEvaluations: Int = 128,
        momentMethod: VivoOmicsNBMomentMethod = .adaptive, maximumFitIterations: Int = 100,
        maximumMomentEvaluations: Int = 1_000_000,
        maximumNeighborhoodVisits: Int = 100_000_000) throws -> VivoOmicsNBQLNativeAbundanceFit {
        let qr = try VivoOmicsQR(design: design), n = qr.observationCount
        guard counts.count >= 3, counts.count <= 100_000, counts.count <= 1_000_000/n,
              counts.allSatisfy({ $0.count == n }), trendDispersions.count == counts.count else {
            throw VivoOmicsStatisticsError.invalid("native-abundance QL family dimensions")
        }
        var abundance = [VivoOmicsNBAbundanceFit?](repeating: nil,count: counts.count)
        var failures: [VivoOmicsNBQLFailure] = []
        for g in counts.indices {
            try Task.checkCancellation()
            do {
                abundance[g] = try VivoOmicsNBAbundance.fit(counts: counts[g],offsets: offsets,
                    dispersion: trendDispersions[g],priorCount: abundancePriorCount,
                    maximumScoreEvaluations: maximumAbundanceScoreEvaluations)
            } catch {
                if error is CancellationError { throw error }
                failures.append(.init(featureIndex: g,stage: "abundance",message: error.localizedDescription))
            }
        }
        guard failures.isEmpty else {
            return .init(completed: false,abundanceFits: abundance,globalFit: nil,failures: failures)
        }
        let global = try fit(counts: counts,design: design,offsets: offsets,contrast: contrast,
            trendDispersions: trendDispersions,abundanceCovariates: abundance.map { $0!.log2CountsPerMillion },
            momentMethod: momentMethod,maximumFitIterations: maximumFitIterations,
            maximumMomentEvaluations: maximumMomentEvaluations,maximumNeighborhoodVisits: maximumNeighborhoodVisits)
        return .init(completed: global.completed,abundanceFits: abundance,globalFit: global,failures: global.failures)
    }
}
