import Foundation

public enum VivoOmicsNBZeroDonorPolicy: String, Codable, Sendable { case activeDonorProfile }
public enum VivoOmicsNBTrendMethod: String, Codable, Sendable { case parametric, mean, gammaParametric }
public struct VivoOmicsNBCohortOptions: Codable, Sendable, Equatable {
    public var zeroTotalDonorPolicy: VivoOmicsNBZeroDonorPolicy?
    public var trend: VivoOmicsNBTrendMethod = .parametric
    public var minimumTrendGenes: Int = 20
    public var minimumPriorVariance: Double = 0.25
    public var outlierStandardDeviations: Double = 2
    /// nil reports influence without excluding observations or genes.
    public var maximumCooksDistance: Double?
    /// Optional zero-centered normal prior on the requested log2 contrast.
    /// Fixed by the caller; never learned from the test outcomes.
    public var effectPriorStandardDeviationLog2: Double?
    public init() {}
    private enum CodingKeys: String, CodingKey {
        case zeroTotalDonorPolicy, trend, minimumTrendGenes, minimumPriorVariance, outlierStandardDeviations, maximumCooksDistance, effectPriorStandardDeviationLog2
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["zeroTotalDonorPolicy","trend", "minimumTrendGenes", "minimumPriorVariance", "outlierStandardDeviations", "maximumCooksDistance", "effectPriorStandardDeviationLog2"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        zeroTotalDonorPolicy = try c.decodeIfPresent(VivoOmicsNBZeroDonorPolicy.self,forKey: .zeroTotalDonorPolicy)
        trend = try c.decodeIfPresent(VivoOmicsNBTrendMethod.self, forKey: .trend) ?? .parametric
        minimumTrendGenes = try c.decodeIfPresent(Int.self, forKey: .minimumTrendGenes) ?? 20
        minimumPriorVariance = try c.decodeIfPresent(Double.self, forKey: .minimumPriorVariance) ?? 0.25
        outlierStandardDeviations = try c.decodeIfPresent(Double.self, forKey: .outlierStandardDeviations) ?? 2
        maximumCooksDistance = try c.decodeIfPresent(Double.self, forKey: .maximumCooksDistance)
        effectPriorStandardDeviationLog2 = try c.decodeIfPresent(Double.self, forKey: .effectPriorStandardDeviationLog2)
    }
    public func validate() throws {
        guard (20...100_000).contains(minimumTrendGenes), minimumPriorVariance.isFinite,
              (0.01...10).contains(minimumPriorVariance), outlierStandardDeviations.isFinite,
              (1...10).contains(outlierStandardDeviations),
              maximumCooksDistance.map({ $0.isFinite && $0 > 0 }) ?? true,
              effectPriorStandardDeviationLog2.map({ $0.isFinite && (0.01...100).contains($0) }) ?? true else {
            throw VivoOmicsError.invalid("NB trend, prior or influence options")
        }
    }
}
public struct VivoOmicsNBTrend: Codable, Sendable, Equatable {
    public let method: VivoOmicsNBTrendMethod
    public let intercept: Double
    public let inverseMeanCoefficient: Double
    public let referenceFeatureIndices: [Int]
    public let robustLogResidualVariance: Double
    public let samplingLogVariance: Double
    public let priorLogVariance: Double
    public let priorVarianceFloorReached: Bool
    public let iterations: Int
    /// Gamma trend rejection subset; prior variance still uses the full interior cohort.
    public var trendFitFeatureIndices: [Int]? = nil
    public func dispersion(mean: Double) -> Double { intercept + inverseMeanCoefficient / mean }
}
public struct VivoOmicsNBFeatureDiagnostics: Codable, Sendable, Equatable {
    public let featureIndex: Int
    public var supportResolution: VivoOmicsNBSupportResolution?
    public var geneWiseDispersion: Double?
    public var geneWiseLowerBoundary: Bool?
    public var geneWiseUpperBoundary: Bool?
    public var trendDispersion: Double?
    public var dispersionOutlier: Bool?
    public var finalDispersion: Double?
    public var finalFit: VivoOmicsNBFit?
    public var error: String?
    /// Raw natural-log MAP/Laplace diagnostics, including an unconverged attempt.
    public var effectShrinkageFit: VivoOmicsNBContrastMAPFit?
    public var effectShrinkageError: String?
}
public struct VivoOmicsNBEffectShrinkageSummary: Codable, Sendable, Equatable {
    public let method: String
    public let priorStandardDeviationLog2: Double
    public let eligibleFeatures: Int
    public let convergedFeatures: Int
    public let failedFeatures: Int
    public let qualification: String
}
public struct VivoOmicsNBCohortDiagnostics: Codable, Sendable, Equatable {
    public let trend: VivoOmicsNBTrend
    public let features: [VivoOmicsNBFeatureDiagnostics]
    public let qualification: String
    public var effectShrinkage: VivoOmicsNBEffectShrinkageSummary? = nil
}

public enum VivoOmicsNBCohort {
    /// Robust log-residual fit of alpha = a0 + a1 / mean. A named mean-only
    /// alternative is explicit; failed parametric fits do not switch methods.
    public static func fitTrend(means: [Double], dispersions: [Double], featureIndices: [Int],
                                residualDF: Int, options: VivoOmicsNBCohortOptions) throws -> VivoOmicsNBTrend {
        try options.validate()
        guard means.count == dispersions.count, means.count == featureIndices.count,
              means.count >= options.minimumTrendGenes, residualDF > 0,
              means.allSatisfy({ $0.isFinite && $0 > 0 }),
              dispersions.allSatisfy({ $0.isFinite && $0 > 1e-8 && $0 < 100 }) else {
            throw VivoOmicsError.invalid("NB trend needs sufficient interior gene-wise estimates")
        }
        let logs = dispersions.map(log), scale = try VivoOmicsLinearStatistics.median(means)
        let x = means.map { scale / $0 }
        var a = exp(try VivoOmicsLinearStatistics.median(logs)), b = 0.0, iterations = 0
        var gammaReferences: [Int]?
        if options.trend == .gammaParametric {
            var selected=Array(x.indices), outerConverged=false
            b=a*0.1
            for _ in 0..<50 {
                let oldA=a,oldB=b
                func loss(_ aa: Double,_ bb: Double) -> Double {
                    selected.reduce(0) { sum,i in
                        let mu=aa+bb*x[i]
                        return sum+dispersions[i]/mu+log(mu)
                    }
                }
                var converged=false
                for _ in 0..<200 {
                    iterations+=1
                    var h00=0.0,h01=0.0,h11=0.0,g0=0.0,g1=0.0
                    for i in selected {
                        let mu=a+b*x[i],r=dispersions[i]/mu-1,d0=1/mu,d1=x[i]/mu
                        h00+=d0*d0;h01+=d0*d1;h11+=d1*d1;g0+=d0*r;g1+=d1*r
                    }
                    let determinant=h00*h11-h01*h01
                    guard determinant.isFinite,determinant>1e-12*h00*h11 else {
                        throw VivoOmicsError.invalid("NB Gamma trend is not identifiable")
                    }
                    let da=(h11*g0-h01*g1)/determinant,db=(h00*g1-h01*g0)/determinant,before=loss(a,b)
                    var fraction=1.0,accepted=false
                    for _ in 0..<50 {
                        let aa=a+fraction*da,bb=b+fraction*db
                        if aa>1e-12,bb>1e-12 {
                            let after=loss(aa,bb)
                            if after.isFinite,after<=before+1e-12*max(1,abs(before)) {
                                let change=max(abs(aa/a-1),abs(bb/b-1))
                                a=aa;b=bb;accepted=true
                                converged=change<1e-8 && max(abs(g0)/sqrt(h00),abs(g1)/sqrt(h11))<1e-6
                                break
                            }
                        }
                        fraction*=0.5
                    }
                    guard accepted else { throw VivoOmicsError.invalid("NB Gamma trend failed to improve") }
                    if converged { break }
                }
                guard converged,a>1e-10,b/scale>1e-10 else { throw VivoOmicsError.invalid("NB Gamma trend did not converge to positive coefficients") }
                let retained=selected.filter { i in
                    let ratio=dispersions[i]/(a+b*x[i]);return ratio>=1e-4 && ratio<15
                }
                guard retained.count>=options.minimumTrendGenes else { throw VivoOmicsError.invalid("NB Gamma trend has too few retained genes") }
                let unchanged=retained==selected
                selected=retained
                if unchanged,pow(log(a/oldA),2)+pow(log(b/oldB),2)<1e-6 { outerConverged=true;break }
            }
            guard outerConverged else { throw VivoOmicsError.invalid("NB Gamma trend rejection did not converge") }
            gammaReferences=selected.map { featureIndices[$0] }
        }
        if options.trend == .parametric {
            b = a * 0.1
            func objective(_ a: Double, _ b: Double) -> Double {
                zip(x,logs).reduce(0) { sum, pair in
                    let prediction = a + b * pair.0
                    if prediction <= 0 { return .infinity }
                    let r = abs(pair.1-log(prediction))
                    return sum + (r <= 1.345 ? r*r/2 : 1.345*(r-1.345/2))
                }
            }
            var converged = false
            for iteration in 1...200 {
                iterations = iteration
                var h00 = 0.0, h01 = 0.0, h11 = 0.0, g0 = 0.0, g1 = 0.0
                for i in x.indices {
                    let mu = a+b*x[i], r = logs[i]-log(mu), w = min(1,1.345/max(abs(r),1e-100))
                    let d0 = 1/mu, d1 = x[i]/mu
                    h00 += w*d0*d0; h01 += w*d0*d1; h11 += w*d1*d1
                    g0 += w*d0*r; g1 += w*d1*r
                }
                let det = h00*h11-h01*h01
                guard det.isFinite, det > 1e-12*h00*h11 else { throw VivoOmicsError.invalid("NB parametric trend is not identifiable") }
                let da = (h11*g0-h01*g1)/det, db = (h00*g1-h01*g0)/det
                let before = objective(a,b)
                var fraction = 1.0, accepted = false
                for _ in 0..<40 {
                    let aa = max(0,a+fraction*da), bb = max(0,b+fraction*db)
                    let after = objective(aa,bb)
                    if after.isFinite && after <= before {
                        let change = x.map { abs((aa+bb*$0)/(a+b*$0)-1) }.max()!
                        a = aa; b = bb; accepted = true
                        if change < 1e-7 { converged = true }
                        break
                    }
                    fraction *= 0.5
                }
                if converged { break }
                if !accepted { throw VivoOmicsError.invalid("NB parametric trend step did not improve") }
            }
            guard converged else { throw VivoOmicsError.invalid("NB parametric trend did not converge") }
        }
        let residuals = x.indices.map { logs[$0]-log(a+b*x[$0]) }
        let center = try VivoOmicsLinearStatistics.median(residuals)
        let mad = try VivoOmicsLinearStatistics.median(residuals.map { abs($0-center) }) / 0.6744897501960817
        let sampling = VivoOmicsLinearStatistics.trigamma(Double(residualDF)/2)
        let excess = mad*mad-sampling
        return .init(method: options.trend, intercept: a, inverseMeanCoefficient: b*scale,
            referenceFeatureIndices: featureIndices, robustLogResidualVariance: mad*mad,
            samplingLogVariance: sampling, priorLogVariance: max(options.minimumPriorVariance,excess),
            priorVarianceFloorReached: excess <= options.minimumPriorVariance, iterations: iterations,trendFitFeatureIndices: gammaReferences)
    }
    static func evaluate(metadata: VivoSingleCellCountMetadata, entries: [[(row: Int,count: UInt64)]],
                         design: VivoOmicsDesignMatrix, request: VivoOmicsExpressionContrast) throws -> VivoOmicsExpressionResult {
        let options = request.negativeBinomialOptions ?? .init()
        guard options.zeroTotalDonorPolicy == nil || request.design == .pairedDonors else {
            throw VivoOmicsError.invalid("NB zero-total donor policy requires paired donors")
        }
        let n = design.rows.count, offsets = design.sizeFactorValues.map(log)
        let count = entries.count
        var totals = [UInt64](repeating: 0,count: count), means = [Double](repeating: 0,count: count)
        var profiles = [VivoOmicsNBDispersionFit?](repeating: nil,count: count)
        var diagnostics = entries.indices.map { VivoOmicsNBFeatureDiagnostics(featureIndex: $0) }
        var statuses = [VivoOmicsExpressionStatus](repeating: .filteredLowExpression,count: count)
        func response(_ gene: Int) -> [UInt64] {
            var y = [UInt64](repeating: 0,count: n)
            for entry in entries[gene] { y[entry.row] = entry.count }
            return y
        }
        for gene in entries.indices {
            try Task.checkCancellation()
            for entry in entries[gene] {
                totals[gene] = try vivoOmicsSum(totals[gene],entry.count)
                means[gene] += Double(entry.count)/design.sizeFactorValues[entry.row]/Double(n)
            }
            if totals[gene] < request.minimumFeatureCounts || entries[gene].count < request.minimumExpressingPseudobulks { continue }
            let y = response(gene)
            do {
                if options.zeroTotalDonorPolicy != nil {
                    diagnostics[gene].supportResolution = try VivoOmicsNBSupport.resolve(counts: y,design: design,request: request)
                    if let resolution=diagnostics[gene].supportResolution,resolution.outcome != .ready {
                        statuses[gene] = resolution.outcome == .insufficientReplication ? .insufficientActiveDonors : .rankDeficientSupport
                        diagnostics[gene].error = "Active-donor profile unavailable: " + resolution.outcome.rawValue
                        continue
                    }
                }
                let resolution=diagnostics[gene].supportResolution
                let rows=resolution?.retainedObservationIndices ?? Array(0..<n)
                let matrix=resolution?.rows ?? design.rows, contrast=resolution?.contrast ?? design.contrast
                let counts=rows.map { y[$0] }, localOffsets=rows.map { offsets[$0] }
                if VivoOmicsNegativeBinomial.positiveSupportIsRankDeficient(counts: counts,design: matrix) {
                    statuses[gene] = .rankDeficientSupport
                    diagnostics[gene].error = "Positive-count support is rank deficient; no inferential fit"
                    continue
                }
                let profile = try VivoOmicsNegativeBinomial.estimateDispersion(counts: counts,design: matrix,offsets: localOffsets,contrast: contrast)
                profiles[gene] = profile
                diagnostics[gene].geneWiseDispersion = profile.fit.dispersion
                diagnostics[gene].geneWiseLowerBoundary = profile.lowerBoundary
                diagnostics[gene].geneWiseUpperBoundary = profile.upperBoundary
                statuses[gene] = .tested
            } catch is CancellationError { throw CancellationError() }
            catch { statuses[gene] = .numericalFailure; diagnostics[gene].error = error.localizedDescription }
        }
        // Keep the original full-design reference cohort and its sampling variance.
        // Gene-specific profiles borrow this prior; they do not mix residual DFs
        // into the full-cohort prior-variance estimate.
        let reference = entries.indices.filter { diagnostics[$0].supportResolution == nil && (profiles[$0].map { !$0.lowerBoundary && !$0.upperBoundary } ?? false) }
        let trend = try fitTrend(means: reference.map { means[$0] },dispersions: reference.map { profiles[$0]!.fit.dispersion },
                                 featureIndices: reference,residualDF: design.residualDegreesOfFreedom,options: options)
        var low = 0.0, high = 10.0
        for _ in 0..<80 {
            let mid = (low+high)/2
            if erfc(mid/sqrt(2)) > 1-request.intervalCoverage { low = mid } else { high = mid }
        }
        let critical = (low+high)/2
        var features: [VivoOmicsExpressionFeature] = [], tested: [Int] = [], probabilities: [Double] = []
        for gene in entries.indices {
            try Task.checkCancellation()
            var result = VivoOmicsExpressionFeature(featureIndex: gene,featureID: metadata.features[gene].id,status: statuses[gene],
                totalCounts: totals[gene],expressingPseudobulks: entries[gene].count,meanNormalizedCount: means[gene],
                log2FoldChange: nil,residualVariance: nil,posteriorVariance: nil,standardError: nil,tStatistic: nil,
                degreesOfFreedom: nil,intervalLower: nil,intervalUpper: nil,pValue: nil,adjustedPValue: nil)
            if let profile = profiles[gene] {
                do {
                    let resolution=diagnostics[gene].supportResolution
                    let rows=resolution?.retainedObservationIndices ?? Array(0..<n)
                    let y=response(gene)
                    let target = trend.dispersion(mean: resolution?.meanNormalizedCount ?? means[gene])
                    guard target.isFinite, (1e-8...100).contains(target) else { throw VivoOmicsError.invalid("NB trend prediction outside dispersion bounds") }
                    diagnostics[gene].trendDispersion = target
                    let outlier = log(profile.fit.dispersion/target) > options.outlierStandardDeviations * sqrt(trend.robustLogResidualVariance)
                    diagnostics[gene].dispersionOutlier = outlier
                    let final = outlier ? profile : try VivoOmicsNegativeBinomial.estimateDispersion(counts: rows.map { y[$0] },design: resolution?.rows ?? design.rows,
                        offsets: rows.map { offsets[$0] },contrast: resolution?.contrast ?? design.contrast,logPriorMean: log(target),logPriorVariance: trend.priorLogVariance)
                    diagnostics[gene].finalDispersion = final.fit.dispersion
                    diagnostics[gene].finalFit = final.fit
                    guard let cooks = final.fit.cooksDistances else { throw VivoOmicsError.invalid("NB influence unavailable on unit-leverage design") }
                    let status: VivoOmicsExpressionStatus
                    if final.lowerBoundary || final.upperBoundary { status = .dispersionBoundary }
                    else if let threshold = options.maximumCooksDistance, cooks.contains(where: { $0 > threshold }) { status = .influentialObservation }
                    else { status = .tested }
                    guard let effect = final.fit.effect, let error = final.fit.standardError, error > 0 else { throw VivoOmicsError.invalid("NB final fit lacks identified effect/information") }
                    let z = effect/error, probability = erfc(abs(z)/sqrt(2))
                    result = .init(featureIndex: gene,featureID: metadata.features[gene].id,status: status,totalCounts: totals[gene],
                        expressingPseudobulks: entries[gene].count,meanNormalizedCount: means[gene],
                        log2FoldChange: effect/log(2),residualVariance: nil,posteriorVariance: nil,standardError: error/log(2),tStatistic: nil,
                        degreesOfFreedom: nil,intervalLower: status == .tested ? (effect-critical*error)/log(2) : nil,
                        intervalUpper: status == .tested ? (effect+critical*error)/log(2) : nil,pValue: status == .tested ? probability : nil,
                        adjustedPValue: nil,zStatistic: status == .tested ? z : nil)
                    if status == .tested { tested.append(gene); probabilities.append(probability) }
                } catch is CancellationError { throw CancellationError() }
                catch {
                    diagnostics[gene].error = error.localizedDescription
                    result = .init(featureIndex: gene,featureID: metadata.features[gene].id,status: .numericalFailure,totalCounts: totals[gene],
                        expressingPseudobulks: entries[gene].count,meanNormalizedCount: means[gene],log2FoldChange: nil,residualVariance: nil,
                        posteriorVariance: nil,standardError: nil,tStatistic: nil,degreesOfFreedom: nil,intervalLower: nil,intervalUpper: nil,pValue: nil,adjustedPValue: nil)
                }
            }
            // Keep shrinkage failures outside the Wald inference path. Only
            // tested genes are eligible, using exactly their retained donors.
            if result.status == .tested, let priorSD = options.effectPriorStandardDeviationLog2 {
                do {
                    let resolution = diagnostics[gene].supportResolution
                    let rows = resolution?.retainedObservationIndices ?? Array(0..<n), y = response(gene)
                    let map = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: rows.map { y[$0] },
                        design: resolution?.rows ?? design.rows, offsets: rows.map { offsets[$0] },
                        contrast: resolution?.contrast ?? design.contrast, dispersion: diagnostics[gene].finalDispersion!,
                        priorStandardDeviation: priorSD*log(2))
                    diagnostics[gene].effectShrinkageFit = map
                    if !map.converged { diagnostics[gene].effectShrinkageError = "NB contrast MAP did not converge" }
                } catch is CancellationError { throw CancellationError() }
                catch { diagnostics[gene].effectShrinkageError = error.localizedDescription }
            }
            features.append(result)
        }
        let adjusted = try VivoOmicsLinearStatistics.benjaminiHochberg(probabilities)
        for (i,gene) in tested.enumerated() { features[gene].adjustedPValue = adjusted[i] }
        return .init(method: options.zeroTotalDonorPolicy == nil ? "donor-aware-NB2-adjusted-profile-log-prior-Wald-v1" : "donor-aware-NB2-active-donor-adjusted-profile-log-prior-Wald-v1",request: request,evidence: metadata.evidence,design: design,
            variancePrior: nil,features: features,testedFeatures: tested.count,
            multiplicityScope: "BH across available NB Wald tests within this contrast; no selection-adjusted or cross-contrast calibration claim",
            negativeBinomial: .init(trend: trend,features: diagnostics,
                qualification: "Experimental NB cohort method; asymptotic Wald intervals, heuristic influence gate if requested, no count replacement; multi-study calibration remains open",
                effectShrinkage: options.effectPriorStandardDeviationLog2.map { priorSD in
                    let completed = diagnostics.filter { $0.effectShrinkageFit?.converged == true }.count
                    return .init(method: "NB2-contrast-normal-prior-count-likelihood-MAP-Laplace-v1",
                        priorStandardDeviationLog2: priorSD,eligibleFeatures: tested.count,convergedFeatures: completed,
                        failedFeatures: tested.count-completed,
                        qualification: "Explicit fixed prior; nuisance coefficients jointly refitted; natural-log fit diagnostics; Laplace SD conditional on fixed dispersion and prior; original Wald tests and BH family unchanged; posterior coverage and biological calibration unqualified")
                }))
    }
}
