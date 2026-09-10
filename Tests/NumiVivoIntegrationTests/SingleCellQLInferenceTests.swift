import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellQLInferenceTests {
    static let design = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
    static let counts: [[UInt64]] = (1...12).map { g in
        let n = UInt64(g); return [n,0,5*n,2*n,1,10*n]
    }
    static func fit(nullIterations: Int = 100,profileEvaluations: Int = 128,cooks: Double? = nil,
                    contrast: [Double] = [0,1]) throws -> VivoOmicsNBQLInferenceFit {
        try VivoOmicsNBQLInference.fit(counts: counts,design: design,
            offsets: [Double](repeating: log(1000),count: 6),contrast: contrast,
            trendDispersions: [Double](repeating: 0.05,count: 12),maximumCooksDistance: cooks,
            maximumNullIterations: nullIterations,maximumProfileEvaluations: profileEvaluations)
    }
    @Test func completeNativeChainOwnsFAndMultiplicity() throws {
        let fit = try Self.fit(), inference = try #require(fit.inference), prior = try #require(fit.moderation)
        #expect(fit.completed && inference.completed && inference.testedIndices.count == 12)
        #expect(inference.poissonBound == "not-applicable-to-modern-adjusted-QL")
        #expect(inference.ordinaryResidualDFCap == 48)
        let values = try inference.tests.map { try #require($0) }
        let bh = try VivoOmicsLinearStatistics.benjaminiHochberg(values.map(\.pValue))
        for i in values.indices {
            let t = values[i], residual = fit.native.globalFit!.adjustedResiduals[i]!
            #expect(abs(t.fStatistic!-t.likelihoodRatio.statistic!/prior.posteriorVariances[i]) < 1e-12)
            #expect(t.denominatorDegreesOfFreedom == min(48,residual.degreesOfFreedom+prior.priorDegreesOfFreedom[i]))
            #expect(t.adjustedPValue == bh[i])
            #expect(t.likelihoodRatio.nullCoefficients[1] == 0)
            #expect(t.likelihoodRatio.nullMaximumScaledScore <= 1e-7)
        }
    }
    @Test func contrastScalingPreservesHypothesis() throws {
        let original = try Self.fit(), reversed = try Self.fit(contrast: [0,-3])
        for i in Self.counts.indices {
            #expect(abs(original.inference!.tests[i]!.pValue-reversed.inference!.tests[i]!.pValue) < 1e-12)
            #expect(original.inference!.tests[i]!.likelihoodRatio.nullMeans == reversed.inference!.tests[i]!.likelihoodRatio.nullMeans)
        }
        #expect(throws: (any Error).self) { try Self.fit(contrast: [0,0]) }
    }
    @Test func failedNullFitsRemainAvailableWithoutProbabilities() throws {
        let fit = try Self.fit(nullIterations: 1), inference = try #require(fit.inference)
        #expect(!fit.completed && !inference.failures.isEmpty)
        #expect(fit.native.completed && fit.moderation != nil)
        for failure in inference.failures {
            let i = try #require(failure.featureIndex)
            #expect(inference.tests[i] == nil)
            #expect(inference.failedNullFits[i] != nil)
        }
    }
    @Test func failedModerationRetainsNativePrerequisites() throws {
        let fit = try Self.fit(profileEvaluations: 4)
        #expect(!fit.completed && fit.native.completed)
        #expect(fit.moderation == nil && fit.inference == nil)
        #expect(fit.failures.contains { $0.stage == "moderation" })
    }
    @Test func influenceExclusionPreservesPriorFamily() throws {
        let original = try Self.fit(), excluded = try Self.fit(cooks: 1e-15)
        #expect(excluded.completed && excluded.inference!.testedIndices.isEmpty)
        #expect(excluded.inference!.excludedInfluentialIndices.count == Self.counts.count)
        #expect(excluded.moderation == original.moderation)
        #expect(excluded.inference!.ordinaryResidualDFCap == original.inference!.ordinaryResidualDFCap)
        #expect(excluded.inference!.tests.allSatisfy { $0 == nil })
    }
    @Test func expressionOptionUsesExistingDesignAndReportsQLStatistics() throws {
        let data = try SingleCellNBCohortTests.fixture()
        var request = VivoOmicsExpressionContrast(id: "ql",controlCondition: "ctrl",treatmentCondition: "stim",design: .pairedDonors)
        request.model = .negativeBinomial; request.minimumCellsPerPseudobulk = 1
        var options = VivoOmicsNBCohortOptions(); options.trend = .mean
        request.negativeBinomialOptions = options
        let baseline = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        options.testMethod = .quasiLikelihoodAdjusted; request.negativeBinomialOptions = options
        let ql = try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        #expect(ql.design == baseline.design && ql.negativeBinomial?.trend == baseline.negativeBinomial?.trend)
        #expect(ql.method == "donor-aware-NB2-trend-adjusted-QL-v1")
        #expect(ql.negativeBinomial!.quasiLikelihood!.inference!.completed)
        #expect(ql.testedFeatures > 80)
        #expect(ql.features.allSatisfy { $0.zStatistic == nil && $0.tStatistic == nil && $0.intervalLower == nil && $0.standardError == nil })
        let family = ql.negativeBinomial!.quasiLikelihood!
        for (i,g) in family.featureIndices.enumerated() {
            let test = try #require(family.inference!.tests[i])
            #expect(ql.features[g].fStatistic == test.fStatistic)
            #expect(ql.features[g].pValue == test.pValue && ql.features[g].adjustedPValue == test.adjustedPValue)
        }
        #expect(baseline.negativeBinomial?.quasiLikelihood == nil)
        #expect(try JSONDecoder().decode(VivoOmicsNBCohortOptions.self,from: JSONEncoder().encode(options)) == options)
    }
    @Test func unsupportedQLCombinationsAreExplicit() throws {
        var options = VivoOmicsNBCohortOptions(); options.testMethod = .quasiLikelihoodAdjusted
        options.zeroTotalDonorPolicy = .activeDonorProfile
        #expect(throws: (any Error).self) { try options.validate() }
        options.zeroTotalDonorPolicy = nil; options.effectPriorStandardDeviationLog2 = 1
        #expect(throws: (any Error).self) { try options.validate() }
        options.effectPriorStandardDeviationLog2 = nil; options.effectPriorEstimation = .weightedUpperQuantile
        #expect(throws: (any Error).self) { try options.validate() }
    }
}
