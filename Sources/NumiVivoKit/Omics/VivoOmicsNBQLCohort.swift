import Foundation

public struct VivoOmicsNBQLCohortDiagnostics: Codable, Sendable, Equatable {
    /// Maps local QL prior/test rows back to the original feature table.
    public let featureIndices: [Int]
    public let averageQLDispersion: Double?
    public let abundanceFits: [VivoOmicsNBAbundanceFit?]
    public let moderation: VivoOmicsQLModeratedVariances?
    public let inference: VivoOmicsNBQLTestFamily?
    public let failures: [VivoOmicsNBQLFailure]
    public let failedUpstreamFit: VivoOmicsNBQLNativeAbundanceFit?
}

extension VivoOmicsNBCohort {
    static func evaluateAdjustedQL(metadata: VivoSingleCellCountMetadata,entries: [[(row: Int,count: UInt64)]],
        design: VivoOmicsDesignMatrix,request: VivoOmicsExpressionContrast,
        totals: [UInt64],means: [Double],profiles: [VivoOmicsNBDispersionFit?],
        diagnostics originalDiagnostics: [VivoOmicsNBFeatureDiagnostics],
        statuses originalStatuses: [VivoOmicsExpressionStatus],trend: VivoOmicsNBTrend) throws -> VivoOmicsExpressionResult {
        let options = request.negativeBinomialOptions ?? .init()
        var diagnostics = originalDiagnostics, statuses = originalStatuses
        var genes: [Int] = [], dispersions: [Double] = [], counts: [[UInt64]] = []
        for g in entries.indices where profiles[g] != nil {
            try Task.checkCancellation()
            let target = trend.dispersion(mean: means[g])
            guard target.isFinite, (1e-8...100).contains(target) else {
                statuses[g] = .numericalFailure; diagnostics[g].error = "QL trend prediction outside dispersion bounds"
                continue
            }
            diagnostics[g].trendDispersion = target
            genes.append(g); dispersions.append(target)
            var y = [UInt64](repeating: 0,count: design.rows.count)
            for entry in entries[g] { y[entry.row] = entry.count }
            counts.append(y)
        }
        var attempt: VivoOmicsNBQLInferenceFit?, ownerFailures: [VivoOmicsNBQLFailure] = []
        do {
            attempt = try VivoOmicsNBQLInference.fit(counts: counts,design: design.rows,
                offsets: design.sizeFactorValues.map(log),contrast: design.contrast,
                trendDispersions: dispersions,maximumCooksDistance: options.maximumCooksDistance)
        } catch is CancellationError { throw CancellationError() }
        catch { ownerFailures.append(.init(featureIndex: nil,stage: "family",message: String(describing: error))) }
        let failures = (attempt?.failures ?? [])+ownerFailures
        var local = [Int:Int]()
        for (i,g) in genes.enumerated() { local[g] = i }
        let excluded = Set(attempt?.inference?.excludedInfluentialIndices ?? [])
        var features: [VivoOmicsExpressionFeature] = []
        for g in entries.indices {
            var full: VivoOmicsNBFit?, test: VivoOmicsNBQLTest?, residual: Double?, posterior: Double?
            if let i = local[g] {
                full = attempt?.native.globalFit?.refittedFits[i]
                diagnostics[g].finalFit = full; diagnostics[g].finalDispersion = full?.dispersion
                residual = attempt?.native.globalFit?.adjustedResiduals[i]?.quasiDispersion
                posterior = attempt?.moderation?.posteriorVariances[i]
                test = attempt?.inference?.tests[i]
                if excluded.contains(i) { statuses[g] = .influentialObservation }
                else if let test {
                    statuses[g] = .tested; diagnostics[g].likelihoodRatioFit = test.likelihoodRatio
                } else {
                    statuses[g] = .numericalFailure
                    diagnostics[g].likelihoodRatioFit = attempt?.inference?.failedNullFits[i]
                    let messages = failures.filter { $0.featureIndex == nil || $0.featureIndex == i }.map { $0.stage+": "+$0.message }
                    diagnostics[g].error = messages.isEmpty ? "QL family prerequisite unavailable; see retained family failures" : messages.joined(separator: "; ")
                }
            }
            features.append(.init(featureIndex: g,featureID: metadata.features[g].id,status: statuses[g],
                totalCounts: totals[g],expressingPseudobulks: entries[g].count,meanNormalizedCount: means[g],
                log2FoldChange: full?.effect.map { $0/log(2) },residualVariance: residual,posteriorVariance: posterior,
                standardError: nil,tStatistic: nil,degreesOfFreedom: test?.denominatorDegreesOfFreedom,
                intervalLower: nil,intervalUpper: nil,pValue: test?.pValue,adjustedPValue: test?.adjustedPValue,
                fStatistic: test?.fStatistic))
        }
        let ql = VivoOmicsNBQLCohortDiagnostics(featureIndices: genes,
            averageQLDispersion: attempt?.native.globalFit?.averageQuasiDispersion,
            abundanceFits: attempt?.native.abundanceFits ?? [],moderation: attempt?.moderation,
            inference: attempt?.inference,failures: failures,
            failedUpstreamFit: attempt?.native.completed == false ? attempt?.native : nil)
        return .init(method: "donor-aware-NB2-trend-adjusted-QL-v1",request: request,evidence: metadata.evidence,
            design: design,variancePrior: nil,features: features,
            testedFeatures: attempt?.inference?.testedIndices.count ?? 0,
            multiplicityScope: "BH across available one-contrast adjusted QL tests; prior uses the complete identified common-design family, including any subsequently influence-excluded rows; no cross-contrast or selection-adjusted calibration claim",
            negativeBinomial: .init(trend: trend,features: diagnostics,
                qualification: "Experimental native adjusted-residual QL: common/trended dispersion, native abundance and global refit, robust unequal-DF prior, constrained LR/posterior F test; denominator DF capped by complete-family ordinary residual DF; Poisson bound applies only to legacy QL and is disabled here; no QL effect intervals or biological calibration claim",
                quasiLikelihood: ql))
    }
}
