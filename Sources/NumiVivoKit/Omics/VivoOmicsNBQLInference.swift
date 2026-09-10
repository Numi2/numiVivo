import Foundation

public struct VivoOmicsNBQLTest: Codable, Sendable, Equatable {
    public let likelihoodRatio: VivoOmicsNBLikelihoodRatioFit
    /// nil when a positive F statistic exceeds finite Double representation.
    public let fStatistic: Double?
    /// nil only for the exactly zero statistic (whose p-value is one).
    public let logFStatistic: Double?
    public let numeratorDegreesOfFreedom: Int
    public let denominatorDegreesOfFreedom: Double
    public let logPValue: Double
    public let pValue: Double
    public var adjustedPValue: Double?
}
public struct VivoOmicsNBQLTestFamily: Codable, Sendable, Equatable {
    public let completed: Bool
    public let tests: [VivoOmicsNBQLTest?]
    public let testedIndices: [Int]
    public let excludedInfluentialIndices: [Int]
    public let ordinaryResidualDFCap: Double
    public let failedNullFits: [VivoOmicsNBLikelihoodRatioFit?]
    public let failures: [VivoOmicsNBQLFailure]
    public let poissonBound: String
}
public struct VivoOmicsNBQLInferenceFit: Codable, Sendable, Equatable {
    public let completed: Bool
    public let native: VivoOmicsNBQLNativeAbundanceFit
    public let moderation: VivoOmicsQLModeratedVariances?
    public let inference: VivoOmicsNBQLTestFamily?
    public let failures: [VivoOmicsNBQLFailure]
}

/// Modern adjusted-residual QL inference at a supplied common/trended NB dispersion.
/// Dispersion must not be a per-gene MAP substitute for the cohort trend.
public enum VivoOmicsNBQLInference {
    public static func fit(counts: [[UInt64]], design: [[Double]], offsets: [Double],
        contrast: [Double], trendDispersions: [Double], maximumCooksDistance: Double? = nil,
        maximumFitIterations: Int = 100, maximumNullIterations: Int = 100,
        maximumMomentEvaluations: Int = 1_000_000,
        maximumProfileEvaluations: Int = 128) throws -> VivoOmicsNBQLInferenceFit {
        guard (1...1000).contains(maximumNullIterations),
              contrast.allSatisfy(\.isFinite), contrast.contains(where: { $0 != 0 }),
              maximumCooksDistance.map({ $0.isFinite && $0 > 0 }) ?? true else {
            throw VivoOmicsStatisticsError.invalid("QL contrast, null-fit or influence bound")
        }
        let native = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: counts,design: design,
            offsets: offsets,contrast: contrast,trendDispersions: trendDispersions,
            maximumFitIterations: maximumFitIterations,maximumMomentEvaluations: maximumMomentEvaluations)
        guard native.completed, let global = native.globalFit else {
            return .init(completed: false,native: native,moderation: nil,inference: nil,failures: native.failures)
        }
        let moderation: VivoOmicsQLModeratedVariances
        do {
            guard global.adjustedResiduals.allSatisfy({ $0?.quasiDispersion != nil }),
                  native.abundanceFits.allSatisfy({ $0 != nil }) else {
                throw VivoOmicsStatisticsError.invalid("QL moderation requires available residual dispersion and abundance for every row")
            }
            moderation = try VivoOmicsQLModeration.fit(
                variances: global.adjustedResiduals.map { $0!.quasiDispersion! },
                degreesOfFreedom: global.adjustedResiduals.map { $0!.degreesOfFreedom },
                abundance: native.abundanceFits.map { $0!.log2CountsPerMillion },
                maximumProfileEvaluations: maximumProfileEvaluations)
        } catch is CancellationError { throw CancellationError() }
        catch {
            return .init(completed: false,native: native,moderation: nil,inference: nil,
                failures: [.init(featureIndex: nil,stage: "moderation",message: String(describing: error))])
        }
        let inference = try testConditional(counts: counts,design: design,offsets: offsets,contrast: contrast,
            fullFits: global.refittedFits.map { $0! },
            adjustedResidualDF: global.adjustedResiduals.map { $0!.degreesOfFreedom },moderation: moderation,
            maximumCooksDistance: maximumCooksDistance,maximumNullIterations: maximumNullIterations)
        return .init(completed: inference.completed,native: native,moderation: moderation,
                     inference: inference,failures: inference.failures)
    }

    /// Internal reuse requires fits/moderation from these exact identified inputs.
    /// The public entry point computes them under one owner; reference harnesses
    /// bind cached inputs to the original count and native-output hashes.
    static func testConditional(counts: [[UInt64]], design: [[Double]], offsets: [Double],
        contrast: [Double], fullFits: [VivoOmicsNBFit], adjustedResidualDF: [Double],
        moderation: VivoOmicsQLModeratedVariances, maximumCooksDistance: Double? = nil,
        maximumNullIterations: Int = 100) throws -> VivoOmicsNBQLTestFamily {
        let qr = try VivoOmicsQR(design: design), n = qr.observationCount, genes = counts.count
        guard (3...100_000).contains(genes), genes <= 1_000_000/n,
              counts.allSatisfy({ $0.count == n }), fullFits.count == genes,
              adjustedResidualDF.count == genes, moderation.posteriorVariances.count == genes,
              moderation.priorDegreesOfFreedom.count == genes,
              offsets.count == n, offsets.allSatisfy(\.isFinite),
              contrast.count == qr.coefficientCount, contrast.allSatisfy(\.isFinite), contrast.contains(where: { $0 != 0 }),
              adjustedResidualDF.allSatisfy({ $0.isFinite && $0 >= 0 }),
              moderation.posteriorVariances.allSatisfy({ $0.isFinite && $0 > 0 }),
              moderation.priorDegreesOfFreedom.allSatisfy({ $0.isFinite && $0 >= 0 }),
              fullFits.allSatisfy({ $0.converged && !$0.positiveCountDesignRankDeficient && $0.means.count == n }),
              (1...1000).contains(maximumNullIterations),
              maximumCooksDistance.map({ $0.isFinite && $0 > 0 }) ?? true else {
            throw VivoOmicsStatisticsError.invalid("QL inference family, fits, moments, posterior or work bounds")
        }
        let cap = Double(genes*(n-qr.coefficientCount))
        var tests = [VivoOmicsNBQLTest?](repeating: nil,count: genes)
        var failedNullFits = [VivoOmicsNBLikelihoodRatioFit?](repeating: nil,count: genes)
        var tested: [Int] = [], excluded: [Int] = [], failures: [VivoOmicsNBQLFailure] = []
        for g in 0..<genes {
            try Task.checkCancellation()
            do {
                if let threshold = maximumCooksDistance {
                    guard let cooks = fullFits[g].cooksDistances else {
                        throw VivoOmicsStatisticsError.invalid("QL influence is unavailable")
                    }
                    if cooks.contains(where: { $0 > threshold }) { excluded.append(g); continue }
                }
                let lr = try VivoOmicsNegativeBinomial.contrastLikelihoodRatio(counts: counts[g],design: design,
                    offsets: offsets,contrast: contrast,full: fullFits[g],maximumNullIterations: maximumNullIterations)
                failedNullFits[g] = lr
                guard let statistic = lr.statistic else {
                    failures.append(.init(featureIndex: g,stage: "nullFit",message: lr.error ?? "Unavailable constrained fit"))
                    // Keep failed null fits separately from valid test probabilities.
                    failedNullFits[g] = lr
                    continue
                }
                let denominatorDF = min(cap,adjustedResidualDF[g]+moderation.priorDegreesOfFreedom[g])
                guard denominatorDF > 0 else { throw VivoOmicsStatisticsError.invalid("QL denominator has no residual or prior information") }
                let logF: Double?, value: Double?, logP: Double
                if statistic == 0 { logF = nil; value = 0; logP = 0 }
                else {
                    let lf = log(statistic)-log(moderation.posteriorVariances[g])
                    logF = lf
                    let f = exp(lf); value = f.isFinite ? f : nil
                    logP = try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: lf,numeratorDF: 1,denominatorDF: denominatorDF).upper
                }
                tests[g] = .init(likelihoodRatio: lr,fStatistic: value,logFStatistic: logF,numeratorDegreesOfFreedom: 1,
                    denominatorDegreesOfFreedom: denominatorDF,logPValue: logP,pValue: exp(logP),adjustedPValue: nil)
                tested.append(g)
                failedNullFits[g] = nil
            } catch is CancellationError { throw CancellationError() }
            catch { failures.append(.init(featureIndex: g,stage: "hypothesis",message: String(describing: error))) }
        }
        let adjusted = try VivoOmicsLinearStatistics.benjaminiHochberg(tested.map { tests[$0]!.pValue })
        for (i,g) in tested.enumerated() { tests[g]!.adjustedPValue = adjusted[i] }
        return .init(completed: failures.isEmpty,tests: tests,testedIndices: tested,excludedInfluentialIndices: excluded,
            ordinaryResidualDFCap: cap,failedNullFits: failedNullFits,failures: failures,poissonBound: "not-applicable-to-modern-adjusted-QL")
    }
}
