import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNegativeBinomialTests {
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
}
