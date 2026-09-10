import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNegativeBinomialTests {
    @Test func unitDevianceResolvesSaturationPoissonAndExtremeMeans() throws {
        // Independent 100-digit saturated likelihood differences at exact
        // binary64 inputs, including the cancellation-prone large-count case.
        let cases: [(UInt64,Double,Double,Double)] = [
            (0,2,0.1,3.6464311358790925069),
            (7,3,0,3.8621700454208505919),
            (100,100.00000001,0.1,9.0908976751672006888e-20),
            (9_007_199_254_740_992,9_007_199_254_740_991,100,1.232595164407831127e-34),
            (1_000_000_000,1_000_000_001,1e-8,9.0909090793388428147e-11),
            (1,1e-200,0.1,918.93721324192312672),
            (1,1e200,100,9.1902407036527832081)]
        for (y,mu,a,expected) in cases {
            let value = try VivoOmicsNegativeBinomial.unitDeviance(count: y,mean: mu,dispersion: a)
            #expect(abs(value/expected-1) < 2e-12)
        }
        #expect(try VivoOmicsNegativeBinomial.unitDeviance(count: 0,mean: 0,dispersion: 0) == 0)
        #expect(try VivoOmicsNegativeBinomial.unitDeviance(count: 3,mean: 3,dispersion: 0.2) == 0)
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.unitDeviance(count: 1,mean: 0,dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.unitDeviance(count: 9_007_199_254_740_993,mean: 2,dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.unitDeviance(count: 1,mean: .nan,dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.unitDeviance(count: 1,mean: 2,dispersion: -1) }
    }
    @Test func likelihoodRatioMatchesAnalyticGroupMeansAndReparameterization() throws {
        let y: [UInt64] = [11,20,12,55,43,61]
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]], offsets = Array(repeating: 0.0,count: 6)
        let fit = try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: y,design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1)
        let groupMeans = Array(repeating: 43.0/3,count: 3)+Array(repeating: 53.0,count: 3)
        var expected = 0.0
        for i in y.indices {
            expected += 2 * (Double(y[i])*log(groupMeans[i]/(202.0/6)) - (Double(y[i])+10)*log((1+0.1*groupMeans[i])/(1+0.1*202/6)))
        }
        #expect(fit.error == nil && fit.degreesOfFreedom == 1)
        #expect(abs(try #require(fit.statistic)-expected) < 1e-9)
        #expect(fit.nullMeans.allSatisfy { abs($0-202.0/6) < 1e-6 })
        let rebased = try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: y,design: x.map { [1-$0[1],$0[1]] },offsets: offsets.map { $0+2 },contrast: [-3,3],dispersion: 0.1)
        #expect(abs(try #require(rebased.statistic)-expected) < 1e-9)
        #expect(abs(rebased.nullCoefficients[0]-rebased.nullCoefficients[1]) < 1e-12)
        #expect(abs(try #require(fit.pValue)-erfc(sqrt(expected/2))) < 1e-12)
    }
    @Test func likelihoodRatioPreservesNullAndSupportBoundaries() throws {
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]], offsets = Array(repeating: 0.0,count: 6)
        let null = try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: [10,12,14,10,12,14],design: x,offsets: offsets,contrast: [0,1],dispersion: 0.2)
        #expect(abs(try #require(null.statistic)) < 1e-10)
        for c in [[0.0,0],[.nan,1]] {
            #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: [1,2,3,4,5,6],design: x,offsets: offsets,contrast: c,dispersion: 0.1) }
        }
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: [1,2,3,0,0,0],design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1) }
        let scalar = try VivoOmicsNegativeBinomial.fitContrastLikelihoodRatio(counts: [1,1,1],design: [[1],[1],[1]],offsets: [0,0,0],contrast: [2],dispersion: 0.2)
        #expect(scalar.nullCoefficients == [0] && scalar.nullMeans == [1,1,1])
        let scalarStatistic = try #require(scalar.statistic)
        #expect(scalar.nullIterations == 0 && abs(scalarStatistic) < 1e-10)
    }
    // Independent mpmath 1.3.0, 80-digit loggamma calculation.
    @Test func logMassMatchesHighPrecisionReference() throws {
        let cases: [(UInt64,Double,Double,Double)] = [
            (0, 1.0, 0.2, -0.9116077839697732),
            (12, 12.0, 0.2, -2.7919657092583265),
            (3, 0.2, 0.01, -6.7961145226363335),
            (1000, 800.0, 0.1, -6.954183343333024),
            (1000000, 1000000.0, 1e-08, -7.831669060954978),
            (0, 100.0, 1e-08, -99.99995000003334),
            (8, 11.0, 100, -6.73606751589651),
            (12, 12.0, 0.125, -2.632725009712935),
        ]
        for (y,mu,a,expected) in cases {
            #expect(abs(try VivoOmicsNegativeBinomial.logMass(count: y,mean: mu,dispersion: a)-expected) < 1e-8)
        }
    }
    @Test func interceptHasAnalyticMeanAndInformation() throws {
        let y: [UInt64] = [2,5,10,12,20,23]
        let result = try VivoOmicsNegativeBinomial.fit(counts: y,design: y.map { _ in [1.0] },
            offsets: y.map { _ in 0 },contrast: [1],dispersion: 0.2)
        #expect(result.converged)
        #expect(abs(result.coefficients[0]-log(12)) < 1e-8)
        let expectedError = sqrt(3.4 / 72.0)
        let errorDifference = (result.standardError ?? Double.nan) - expectedError
        #expect(abs(errorDifference) < 1e-9)
        #expect(abs(result.leverage.reduce(0,+)-1) < 1e-12)
    }
    @Test func offsetReparameterizationPreservesFit() throws {
        let y: [UInt64] = [11,20,12,55,43,61]
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        let a = try VivoOmicsNegativeBinomial.fit(counts: y,design: x,offsets: [0,0.2,0.1,0.3,0.1,0.2],contrast: [0,1],dispersion: 0.1)
        let b = try VivoOmicsNegativeBinomial.fit(counts: y,design: x,offsets: [2,2.2,2.1,2.3,2.1,2.2],contrast: [0,1],dispersion: 0.1)
        #expect(a.converged && b.converged)
        #expect(abs((a.effect ?? .nan)-(b.effect ?? .nan)) < 1e-9)
        #expect(abs(a.logLikelihood-b.logLikelihood) < 1e-9)
        #expect(abs(a.coefficients[0]-b.coefficients[0]-2) < 1e-9)
    }
    @Test func invalidInputsAndNonconvergenceAreNotSuccessfulFits() throws {
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fit(counts: [0,0,0,0,0,0],design: x,offsets: Array(repeating: 0,count: 6),contrast: [0,1],dispersion: 0.1) }
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fit(counts: [1,2,3,4,5,6],design: x.map { $0+[1] },offsets: Array(repeating: 0,count: 6),contrast: [0,1,0],dispersion: 0.1) }
        let unfinished = try VivoOmicsNegativeBinomial.fit(counts: [1,2,1000,3,4000,5],design: x,
            offsets: Array(repeating: 0,count: 6),contrast: [0,1],dispersion: 0.1,maximumIterations: 1)
        #expect(!unfinished.converged)
        #expect(unfinished.maximumScaledScore > 1e-7)
    }
    @Test func zeroOnlyCovariateSupportDoesNotProduceInference() throws {
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        let y: [UInt64] = [10,12,9,0,0,0]
        let value = try VivoOmicsNegativeBinomial.fit(counts: y,design: x,offsets: Array(repeating: 0,count: 6),contrast: [0,1],dispersion: 0.15)
        #expect(value.positiveCountDesignRankDeficient)
        #expect(value.effect == nil && value.standardError == nil && value.coxReidLogLikelihood == nil)
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.estimateDispersion(counts: y,design: x,offsets: Array(repeating: 0,count: 6),contrast: [0,1]) }
    }
    @Test func profileAndExplicitPriorAreDistinct() throws {
        let y: [UInt64] = [2,5,10,12,20,23]
        let x = y.map { _ in [1.0] }, offsets = y.map { _ in 0.0 }
        let mle = try VivoOmicsNegativeBinomial.estimateDispersion(counts: y,design: x,offsets: offsets,contrast: [1])
        let map = try VivoOmicsNegativeBinomial.estimateDispersion(counts: y,design: x,offsets: offsets,contrast: [1],logPriorMean: log(0.1),logPriorVariance: 0.01)
        #expect(mle.fit.converged && map.fit.converged)
        #expect(abs(log(map.fit.dispersion/0.1)) < abs(log(mle.fit.dispersion/0.1)))
        #expect(!mle.lowerBoundary && !mle.upperBoundary)
    }
    @Test func contrastMAPRefitsNuisanceAndPreservesParameterization() throws {
        let y: [UInt64] = [11,20,12,55,43,61]
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]], offsets = [0.0,0.2,0.1,0.3,0.1,0.2]
        let mle = try VivoOmicsNegativeBinomial.fit(counts: y,design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1)
        let map = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y,design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1,priorStandardDeviation: 0.4)
        let cellMeans = x.map { [1-$0[1],$0[1]] }
        let rebased = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y,design: cellMeans,offsets: offsets.map { $0+2 },contrast: [-1,1],dispersion: 0.1,priorStandardDeviation: 0.4)
        #expect(map.converged && rebased.converged)
        #expect(map.effect > 0 && map.effect < mle.effect!)
        #expect(abs(map.coefficients[0]-mle.coefficients[0]) > 0.01)
        #expect(abs(map.effect-rebased.effect) < 1e-8)
        #expect(abs(map.posteriorStandardDeviation-rebased.posteriorStandardDeviation) < 1e-8)
        #expect(zip(map.means,rebased.means).allSatisfy { abs($0-$1) < 1e-6 })
        #expect(map.objective >= mle.logLikelihood-pow(mle.effect!,2)/(2*0.4*0.4))
        let weak = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y,design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1,priorStandardDeviation: 1e6)
        #expect(weak.converged && abs(weak.effect-mle.effect!) < 1e-7)
    }
    @Test func contrastMAPMatchesHighPrecisionScalarCountLikelihood() throws {
        // mpmath 1.3.0, 70 digits: solve 6*(12-exp(b))/(1+0.2*exp(b))-b/0.4^2=0.
        // Guards numerical UInt64 conversion as well as the posterior curvature.
        let fit = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: Array(repeating: 12,count: 6),
            design: Array(repeating: [1.0],count: 6),offsets: Array(repeating: 0,count: 6),
            contrast: [1],dispersion: 0.2,priorStandardDeviation: 0.4)
        #expect(fit.converged)
        #expect(abs(fit.effect-1.9568087071818467) < 1e-8)
        #expect(fit.means.allSatisfy { abs($0-7.07670714626965) < 1e-7 })
        #expect(abs(fit.posteriorStandardDeviation-0.17961700266278638) < 1e-9)
    }
    @Test func contrastMAPPreservesSupportAndReportsIterationExhaustion() throws {
        let x = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]], offsets = Array(repeating: 0.0,count: 6)
        for sd in [0.0,-1,.infinity,.nan,1e-7,1e7] {
            #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fitContrastMAP(counts: [1,2,3,4,5,6],design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1,priorStandardDeviation: sd) }
        }
        for y: [UInt64] in [[0,0,0,0,0,0],[10,12,9,0,0,0]] {
            #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y,design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1,priorStandardDeviation: 1) }
        }
        let unfinished = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: [1,2,1000,3,4000,5],design: x,offsets: offsets,contrast: [0,1],dispersion: 0.1,priorStandardDeviation: 0.01,maximumIterations: 1)
        #expect(!unfinished.converged && unfinished.maximumScaledScore > 1e-7)
    }
    @Test func contrastMAPSupportsFullObservationLimitWithoutAddingReplication() throws {
        let y = (0..<512).map { UInt64(5+$0%7) }, x = y.map { _ in [1.0] }
        let fit = try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y,design: x,offsets: y.map { _ in 0 },contrast: [1],dispersion: 0.2,priorStandardDeviation: 1)
        #expect(fit.converged && fit.means.count == 512)
        #expect(throws: (any Error).self) { try VivoOmicsNegativeBinomial.fitContrastMAP(counts: y+[5],design: x+[[1]],offsets: Array(repeating: 0,count: 513),contrast: [1],dispersion: 0.2,priorStandardDeviation: 1) }
    }
}
