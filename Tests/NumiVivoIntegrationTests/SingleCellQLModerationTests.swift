import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellQLModerationTests {
    @Test func specialFunctionsRespectRecurrencesAndTailSymmetry() throws {
        for x in [0.005, 0.5, 1, 10, 1000] {
            let h = try VivoOmicsQLSpecialFunctions.logMinusDigamma(x)
            let next = try VivoOmicsQLSpecialFunctions.logMinusDigamma(x+1)
            #expect(abs(h-next-(1/x-log1p(1/x))) < 1e-10)
            let t = try VivoOmicsQLSpecialFunctions.trigamma(x)
            let tn = try VivoOmicsQLSpecialFunctions.trigamma(x+1)
            #expect(abs((t-tn)-1/(x*x)) < 1e-8)
            let increment = try VivoOmicsQLSpecialFunctions.logGammaIncrement(base: x,increment: 1)
            #expect(abs(increment-log(x)) < 1e-10)
        }
        for value in [-100.0, -4, 0, 4, 100] {
            let f = try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: value,numeratorDF: 3,denominatorDF: 18)
            let inverse = try VivoOmicsQLSpecialFunctions.logFTails(logStatistic: -value,numeratorDF: 18,denominatorDF: 3)
            #expect(abs(f.lower-inverse.upper) < 1e-10)
            #expect(abs(exp(f.lower)+exp(f.upper)-1) < 1e-12)
        }
    }
    @Test func precisionWeightsChangeTrendAndLinearFunctionsRemainExact() throws {
        let x = (0..<30).map(Double.init)
        let weights = (0..<30).map { $0 == 15 ? 50.0 : 0.1 }
        let line = try VivoOmicsPrecisionLowess.fit(x: x,y: x.map { 2+3*$0 },weights: weights,span: 0.6)
        for i in x.indices { #expect(abs(line.fitted[i]-(2+3*x[i])) < 1e-9) }
        let y = x.map { sin($0) }
        let weighted = try VivoOmicsPrecisionLowess.fit(x: x,y: y,weights: weights,span: 0.6)
        let equal = try VivoOmicsPrecisionLowess.fit(x: x,y: y,weights: x.map { _ in 1 },span: 0.6)
        #expect(abs(weighted.fitted[15]-y[15]) < abs(equal.fitted[15]-y[15]))
        #expect(throws: (any Error).self) {
            try VivoOmicsPrecisionLowess.fit(x: x,y: y,weights: weights,span: 0.6,maximumNeighborhoodVisits: 1)
        }
        #expect(throws: (any Error).self) {
            try VivoOmicsPrecisionLowess.fit(x: [0,1,2],y: [1e308,-1e308,1e308],weights: [1,2,3],span: 0.6)
        }
    }
    @Test func constantVariancesChooseExplicitUpperEndpoint() throws {
        let fit = try VivoOmicsQLModeration.fit(variances: [Double](repeating: 1,count: 30),
            degreesOfFreedom: [Double](repeating: 5,count: 30),robust: false)
        #expect(fit.profiles[0].upperBoundary)
        #expect(abs(fit.commonPriorDegreesOfFreedom-9998) < 1e-6)
        #expect(fit.profiles[0].evaluations.contains { $0.parameter == 0.5 })
        #expect(fit.profiles[0].evaluations.contains { $0.parameter == 0.9998 })
        let heterogeneous = try VivoOmicsQLModeration.fit(
            variances: (0..<30).map { $0%2 == 0 ? 1e-8 : 1e8 },
            degreesOfFreedom: [Double](repeating: 5,count: 30),robust: false)
        #expect(heterogeneous.profiles[0].lowerBoundary)
        #expect(heterogeneous.commonPriorDegreesOfFreedom == 2)
    }
    @Test func priorAndPosteriorRespectVarianceScaleAndPermutation() throws {
        let x = (0..<40).map { 0.5+Double(($0*17)%13)/8 }
        let df = (0..<40).map { 2+Double($0%5) }
        let fit = try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: df,robust: false)
        let scaled = try VivoOmicsQLModeration.fit(variances: x.map { $0*100 },degreesOfFreedom: df,robust: false)
        let reversed = try VivoOmicsQLModeration.fit(variances: x.reversed(),degreesOfFreedom: df.reversed(),robust: false)
        for i in x.indices {
            #expect(abs(scaled.posteriorVariances[i]/100-fit.posteriorVariances[i]) < 1e-5)
            #expect(abs(reversed.posteriorVariances[x.count-1-i]-fit.posteriorVariances[i]) < 1e-8)
            let d = fit.priorDegreesOfFreedom[i], s = fit.priorScales[i]
            #expect(abs(fit.posteriorVariances[i]-(df[i]*x[i]+d*s)/(df[i]+d)) < 1e-12)
            #expect(fit.posteriorVariances[i] >= min(x[i],s))
            #expect(fit.posteriorVariances[i] <= max(x[i],s))
        }
    }
    @Test func robustOutlierReceivesLessPriorInformation() throws {
        var x = (0..<100).map { 0.8+Double($0%11)/25 }
        x[99] = 1e6
        let df = [Double](repeating: 10,count: 100)
        let fit = try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: df)
        #expect(fit.profiles.count == 2)
        #expect(fit.screeningFDRWeights![99] < 0.01)
        #expect(fit.priorDegreesOfFreedom[99] < fit.commonPriorDegreesOfFreedom/10)
        #expect(fit.posteriorVariances[99] > x[99]/2)
        #expect(fit.outlierDegreesOfFreedom != nil)
    }
    @Test func excludedDFFloorAndOriginalPosteriorRemainIdentified() throws {
        var x = (0..<100).map { 0.8+Double($0%11)/25 }
        var df = [Double](repeating: 10,count: 100)
        x[0] = 0; df[1] = 0; x[1] = 1e10; x[99] = 1e6
        let fit = try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: df)
        #expect(fit.excludedLowDFIndices == [1])
        #expect(fit.flooredVarianceIndices.contains(0))
        #expect(!fit.informativeIndices.contains(0) && !fit.informativeIndices.contains(1))
        #expect(fit.profiles.allSatisfy { $0.priorWeights[1] == 0 })
        #expect(fit.posteriorVariances[1] == fit.priorScales[1])
    }
    @Test func workAndInputFailuresDoNotReturnFabricatedPosterior() throws {
        let x = [0.5,1,2,0.8,1.5], df = [Double](repeating: 5,count: 5)
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: df,maximumProfileEvaluations: 4) }
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: df,maximumFeatureEvaluations: 5) }
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(variances: [0,0,1,2,0],degreesOfFreedom: df) }
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(variances: x,degreesOfFreedom: [0,0,0,1,2]) }
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(variances: [.nan,1,2,3,4],degreesOfFreedom: df) }
    }
    @Test func nativeAbundanceIntegrationRequiresCompleteInputs() throws {
        let counts: [[UInt64]] = (1...12).map { g in
            let n = UInt64(g); return [n,0,5*n,2*n,1,10*n]
        }
        let design = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        let fit = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: counts,design: design,
            offsets: [Double](repeating: log(1000),count: 6),contrast: [0,1],
            trendDispersions: [Double](repeating: 0.05,count: 12))
        let moderated = try VivoOmicsQLModeration.fit(from: fit)
        #expect(moderated.posteriorVariances.count == counts.count)
        var zero = counts; zero[2] = [UInt64](repeating: 0,count: 6)
        let failed = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: zero,design: design,
            offsets: [Double](repeating: log(1000),count: 6),contrast: [0,1],
            trendDispersions: [Double](repeating: 0.05,count: 12))
        #expect(throws: (any Error).self) { try VivoOmicsQLModeration.fit(from: failed) }
    }
}
