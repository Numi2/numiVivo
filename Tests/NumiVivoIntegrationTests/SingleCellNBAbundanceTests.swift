import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBAbundanceTests {
    @Test func equalLibrariesHaveClosedFormForPoissonAndNB() throws {
        let counts: [UInt64] = [0,2,7,11]
        for phi in [0.0,1e-8,0.2,100] {
            for prior in [0.0,0.5,2,10] {
                let fit = try VivoOmicsNBAbundance.fit(counts: counts,
                    offsets: [Double](repeating: log(1000),count: 4),dispersion: phi,priorCount: prior)
                let expected = log2(1_000_000*(5+prior)/(1000+2*prior))
                #expect(abs(fit.log2CountsPerMillion-expected) < 1e-10)
                #expect(fit.logProportionBracketWidth <= 2e-12)
            }
        }
        let zero = try VivoOmicsNBAbundance.fit(counts: [0,0,0],offsets: [0,log(10),log(100)],dispersion: 0.2)
        #expect(abs(zero.log2CountsPerMillion-log2(2_000_000.0/41)) < 1e-10)
    }
    @Test func offsetShiftAndPermutationRespectAbundanceDefinition() throws {
        let counts: [UInt64] = [0,3,7,0], offsets = [0.0,log(10),log(100),log(1000)]
        let native = try VivoOmicsNBAbundance.fit(counts: counts,offsets: offsets,dispersion: 0.2,priorCount: 0)
        let shifted = try VivoOmicsNBAbundance.fit(counts: counts,offsets: offsets.map { $0+3 },dispersion: 0.2,priorCount: 0)
        #expect(abs(shifted.log2CountsPerMillion-native.log2CountsPerMillion+3/log(2)) < 1e-10)
        let permuted = try VivoOmicsNBAbundance.fit(counts: counts.reversed(),offsets: offsets.reversed(),dispersion: 0.2,priorCount: 0)
        #expect(abs(permuted.log2CountsPerMillion-native.log2CountsPerMillion) < 1e-10)
        #expect(native.scaledScore < 1e-7)
    }
    @Test func workAndRepresentabilityFailuresRemainExplicit() throws {
        #expect(throws: (any Error).self) { try VivoOmicsNBAbundance.fit(counts: [0,3,7,0],offsets: [0,1,2,3],dispersion: 0.2,maximumScoreEvaluations: 1) }
        #expect(throws: (any Error).self) { try VivoOmicsNBAbundance.fit(counts: [0,0],offsets: [0,0],dispersion: 0.2,priorCount: 0) }
        #expect(throws: (any Error).self) { try VivoOmicsNBAbundance.fit(counts: [UInt64.max],offsets: [0],dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNBAbundance.fit(counts: [1],offsets: [.nan],dispersion: 0.2) }
        #expect(throws: (any Error).self) { try VivoOmicsNBAbundance.fit(counts: [1,1],offsets: [-700,700],dispersion: 0.2) }
    }
    @Test func integratedScaleUsesOriginalCountsAndRetainsBothFailureStages() throws {
        let counts: [[UInt64]] = (1...12).map { g in
            let n = UInt64(g)
            return [n,0,5*n,2*n,1,10*n]
        }
        let design = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]]
        let offsets = [Double](repeating: log(1000),count: 6)
        let phi = [Double](repeating: 0.05,count: 12)
        let fit = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: counts,design: design,
            offsets: offsets,contrast: [0,1],trendDispersions: phi)
        #expect(fit.completed && fit.failures.isEmpty)
        let first = try VivoOmicsNegativeBinomial.fit(counts: counts[0],design: design,offsets: offsets,
            contrast: [0,1],dispersion: phi[0])
        #expect(fit.globalFit!.initialFits[0] == first)
        let limited = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: counts,design: design,
            offsets: [0,1,2,0,1,2],contrast: [0,1],trendDispersions: phi,maximumAbundanceScoreEvaluations: 1)
        #expect(!limited.completed && limited.globalFit == nil && !limited.failures.isEmpty)
        #expect(limited.failures.allSatisfy { $0.stage == "abundance" })
        var zero = counts; zero[2] = [UInt64](repeating: 0,count: 6)
        let failed = try VivoOmicsNBQLGlobalScale.fitEstimatingAbundance(counts: zero,design: design,
            offsets: offsets,contrast: [0,1],trendDispersions: phi)
        #expect(!failed.completed && failed.abundanceFits[2] != nil && failed.globalFit != nil)
        #expect(failed.failures.contains { $0.featureIndex == 2 && $0.stage == "initialFit" })
    }
}
