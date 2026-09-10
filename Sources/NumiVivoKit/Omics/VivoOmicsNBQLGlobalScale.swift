import Foundation

public struct VivoOmicsNBQLFailure: Codable, Sendable, Equatable {
    public let featureIndex: Int?
    public let stage: String
    public let message: String
}
public struct VivoOmicsNBQLScaleUpdate: Codable, Sendable, Equatable {
    public let inputScale: Double
    public let outputScale: Double?
    public let quasiDispersions: [Double?]
    public let residualDegreesOfFreedom: [Double?]
    public let residualDeviances: [Double?]
    public let eligibleIndices: [Int]
    public let omittedIndices: [Int]
    public let smoother: VivoOmicsLowessFit?
    public let quarterRootQuantile90: Double?
    public let momentEvaluations: Int
}
public struct VivoOmicsNBQLGlobalFit: Codable, Sendable, Equatable {
    public let completed: Bool
    public let averageQuasiDispersion: Double?
    public let initialFits: [VivoOmicsNBFit?]
    public let updates: [VivoOmicsNBQLScaleUpdate]
    public let refittedFits: [VivoOmicsNBFit?]
    public let adjustedResiduals: [VivoOmicsNBAdjustedResiduals?]
    public let failures: [VivoOmicsNBQLFailure]
}

/// Two global scale updates at fixed initial means, followed by one NB refit.
/// Abundance coordinates and trend dispersions are explicit caller inputs.
/// This stage provides no QL prior, posterior moderation or hypothesis test.
public enum VivoOmicsNBQLGlobalScale {
    static func quantile90(_ values: [Double]) throws -> Double {
        guard values.count >= 2, values.allSatisfy(\.isFinite) else {
            throw VivoOmicsStatisticsError.invalid("QL scale quantile requires finite values")
        }
        let sorted = values.sorted(), index = 0.9*Double(values.count-1), low = Int(floor(index))
        let fraction = index-Double(low)
        return (1-fraction)*sorted[low]+fraction*sorted[min(low+1,sorted.count-1)]
    }
    public static func fit(counts: [[UInt64]], design: [[Double]], offsets: [Double],
        contrast: [Double], trendDispersions: [Double], abundanceCovariates: [Double],
        momentMethod: VivoOmicsNBMomentMethod = .adaptive, maximumFitIterations: Int = 100,
        maximumMomentEvaluations: Int = 1_000_000, maximumNeighborhoodVisits: Int = 100_000_000) throws -> VivoOmicsNBQLGlobalFit {
        let qr = try VivoOmicsQR(design: design), n = qr.observationCount, genes = counts.count
        guard genes >= 3, genes <= 100_000, genes <= 1_000_000/n,
              counts.allSatisfy({ $0.count == n }), offsets.count == n, offsets.allSatisfy(\.isFinite),
              contrast.count == qr.coefficientCount, contrast.allSatisfy(\.isFinite),
              trendDispersions.count == genes, abundanceCovariates.count == genes,
              abundanceCovariates.allSatisfy(\.isFinite),
              trendDispersions.allSatisfy({ $0.isFinite && (1e-8...100).contains($0) }),
              (1...1000).contains(maximumFitIterations), (1...10_000_000).contains(maximumMomentEvaluations),
              maximumNeighborhoodVisits > 0 else {
            throw VivoOmicsStatisticsError.invalid("global QL family dimensions, covariates, dispersions or work limits")
        }
        var initial = [VivoOmicsNBFit?](repeating: nil,count: genes)
        var refitted = initial, residuals = [VivoOmicsNBAdjustedResiduals?](repeating: nil,count: genes)
        var updates: [VivoOmicsNBQLScaleUpdate] = [], failures: [VivoOmicsNBQLFailure] = []
        var globalScale: Double?
        func result() -> VivoOmicsNBQLGlobalFit {
            .init(completed: failures.isEmpty && updates.count == 2 && residuals.allSatisfy({ $0 != nil }),
                averageQuasiDispersion: globalScale,initialFits: initial,updates: updates,
                refittedFits: refitted,adjustedResiduals: residuals,failures: failures)
        }
        func record(_ error: Error, index: Int?, stage: String) throws {
            if error is CancellationError { throw error }
            failures.append(.init(featureIndex: index,stage: stage,message: error.localizedDescription))
        }
        for g in 0..<genes {
            try Task.checkCancellation()
            do {
                let fit = try VivoOmicsNegativeBinomial.fit(counts: counts[g],design: design,offsets: offsets,
                    contrast: contrast,dispersion: trendDispersions[g],maximumIterations: maximumFitIterations)
                initial[g] = fit
                if !fit.converged || fit.positiveCountDesignRankDeficient {
                    failures.append(.init(featureIndex: g,stage: "initialFit",message: "Unconverged or unidentified NB fit; diagnostics retained"))
                }
            } catch { try record(error,index: g,stage: "initialFit") }
        }
        if !failures.isEmpty { return result() }
        var scale = 1.0
        for update in 1...2 {
            var quasi = [Double?](repeating: nil,count: genes), df = quasi, deviance = quasi
            var eligible: [Int] = [], omitted: [Int] = [], work = 0
            for g in 0..<genes {
                try Task.checkCancellation()
                do {
                    let adjusted = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: counts[g],means: initial[g]!.means,
                        design: design,dispersion: trendDispersions[g],averageQuasiDispersion: scale,
                        maximumTerms: maximumMomentEvaluations,method: momentMethod)
                    quasi[g] = adjusted.quasiDispersion; df[g] = adjusted.degreesOfFreedom; deviance[g] = adjusted.deviance
                    work += adjusted.moments.reduce(0) { $0+$1.evaluatedCounts }
                    if adjusted.quasiDispersion != nil { eligible.append(g) } else { omitted.append(g) }
                } catch { try record(error,index: g,stage: "scaleUpdate\(update)") }
            }
            var smoother: VivoOmicsLowessFit?, quantile: Double?, nextScale: Double?
            if failures.isEmpty {
                do {
                    guard eligible.count >= 3 else { throw VivoOmicsStatisticsError.invalid("global QL scale requires at least three informative residuals") }
                    let fitted = try VivoOmicsRobustLowess.fit(x: eligible.map { abundanceCovariates[$0] },
                        y: eligible.map { sqrt(sqrt(quasi[$0]!)) },maximumNeighborhoodVisits: maximumNeighborhoodVisits)
                    smoother = fitted
                    guard fitted.emptyWeightAnchors == 0 else { throw VivoOmicsStatisticsError.invalid("global QL smoother has empty weighted neighborhoods") }
                    let q = try quantile90(fitted.fitted), candidate = pow(max(1,q),4)
                    quantile = q
                    guard candidate.isFinite else { throw VivoOmicsStatisticsError.invalid("global QL scale overflow") }
                    nextScale = candidate
                } catch { try record(error,index: nil,stage: "scaleUpdate\(update)") }
            }
            updates.append(.init(inputScale: scale,outputScale: nextScale,quasiDispersions: quasi,
                residualDegreesOfFreedom: df,residualDeviances: deviance,eligibleIndices: eligible,
                omittedIndices: omitted,smoother: smoother,quarterRootQuantile90: quantile,momentEvaluations: work))
            if !failures.isEmpty { return result() }
            scale = nextScale!
        }
        globalScale = scale
        for g in 0..<genes {
            try Task.checkCancellation()
            do {
                let fit = try VivoOmicsNegativeBinomial.fit(counts: counts[g],design: design,offsets: offsets,
                    contrast: contrast,dispersion: trendDispersions[g]/scale,maximumIterations: maximumFitIterations)
                refitted[g] = fit
                guard fit.converged, !fit.positiveCountDesignRankDeficient else {
                    failures.append(.init(featureIndex: g,stage: "refit",message: "Unconverged or unidentified scaled NB fit; diagnostics retained"));continue
                }
                residuals[g] = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: counts[g],means: fit.means,
                    design: design,dispersion: trendDispersions[g],averageQuasiDispersion: scale,
                    maximumTerms: maximumMomentEvaluations,method: momentMethod)
            } catch { try record(error,index: g,stage: "refitAndResiduals") }
        }
        return result()
    }
}
