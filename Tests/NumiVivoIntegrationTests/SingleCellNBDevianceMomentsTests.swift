import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBDevianceMomentsTests {
    @Test func momentsResolveSmallCountsAndPoissonWithTailBounds() throws {
        let cases: [(Double,Double,Double,Double)] = [
            (0.01,4,0.07248941841715136,0.2934516848743014),
            (1,0.1,1.1263775802197873,1.259389087010258),
            (10,0,1.0188285396938748,2.0876874939573415)]
        for (mu,a,mean,variance) in cases {
            let m = try VivoOmicsNBResidualAdjustment.moments(mean: mu,dispersion: a)
            #expect(abs(m.mean/mean-1) < 2e-9)
            #expect(abs(m.variance/variance-1) < 2e-9)
            #expect(m.omittedProbabilityBound <= 1e-10)
            #expect(m.meanTruncationBound <= 1e-10*m.mean)
            #expect(m.varianceTruncationBound <= 1e-10*m.variance)
            let tighter = try VivoOmicsNBResidualAdjustment.moments(mean: mu,dispersion: a,relativeTolerance: 1e-12)
            #expect(abs(tighter.mean-m.mean) <= m.meanTruncationBound+1e-13)
            #expect(abs(tighter.variance-m.variance) <= m.varianceTruncationBound+1e-13)
            #expect(tighter.evaluatedCounts >= m.evaluatedCounts)
        }
    }
    @Test func zeroInvalidAndExhaustedMomentsAreExplicit() throws {
        let zero = try VivoOmicsNBResidualAdjustment.moments(mean: 0,dispersion: 0.2)
        #expect(zero.mean == 0 && zero.variance == 0 && zero.degreesOfFreedom == 0)
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.moments(mean: 10,dispersion: 0.2,maximumTerms: 1) }
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.moments(mean: -1,dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.moments(mean: .infinity,dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.moments(mean: 1,dispersion: 101) }
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.moments(mean: 1,dispersion: 0.2,relativeTolerance: 0) }
    }
    @Test func adjustedResidualsPreserveWeightedGeometryAndUnitLeverageExclusion() throws {
        let y: [UInt64] = [5,7,4,8,10,13], mu = [5.0,6,5,9,11,12]
        let design = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        let fit = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: y,means: mu,design: design,dispersion: 0.2,averageQuasiDispersion: 1.3)
        let rebased = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: y,means: mu,design: design.map { [1-$0[1],$0[1]] },dispersion: 0.2,averageQuasiDispersion: 1.3)
        #expect(abs(fit.deviance-rebased.deviance) < 1e-12)
        #expect(abs(fit.degreesOfFreedom-rebased.degreesOfFreedom) < 1e-12)
        #expect(abs(fit.leverage.reduce(0,+)-2) < 1e-12)
        for i in y.indices {
            let raw = try VivoOmicsNegativeBinomial.unitDeviance(count: y[i],mean: mu[i],dispersion: 0.2/1.3)
            #expect(abs(fit.unitDeviance[i]-raw*fit.moments[i].devianceScale) < 1e-12)
        }
        let unique = [[1.0,1],[1,0],[1,0],[1,0],[1,0],[1,0]]
        let limited = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: y,means: mu,design: unique,dispersion: 0.2,averageQuasiDispersion: 1.3)
        #expect(limited.unitDeviance[0] == 0 && limited.unitDegreesOfFreedom[0] == 0)
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: y,means: mu,design: design,dispersion: 0.2,averageQuasiDispersion: 0) }
        #expect(throws: (any Error).self) { try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: [9_007_199_254_740_993]+Array(y.dropFirst()),means: mu,design: unique,dispersion: 0.2,averageQuasiDispersion: 1) }
        #expect(throws: (any Error).self) {
            try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: [0,0,0],means: [1e-200,1e-200,1e-200],design: [[1],[1],[1]],dispersion: 0,averageQuasiDispersion: 1e200)
        }
    }
}
