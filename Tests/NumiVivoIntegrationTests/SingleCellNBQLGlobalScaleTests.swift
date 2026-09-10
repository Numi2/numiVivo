import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBQLGlobalScaleTests {
    @Test func smootherPreservesLinesOrderingAndExplicitWorkLimit() throws {
        let x = (0..<100).map { Double($0)/10 }, y = x.map { 2+3*$0 }
        let fitted = try VivoOmicsRobustLowess.fit(x: x,y: y)
        for i in x.indices { #expect(abs(fitted.fitted[i]-y[i]) < 1e-10) }
        let reversed = try VivoOmicsRobustLowess.fit(x: x.reversed(),y: y.reversed())
        for i in x.indices { #expect(abs(reversed.fitted[99-i]-fitted.fitted[i]) < 1e-10) }
        let constant = try VivoOmicsRobustLowess.fit(x: [1,1,1],y: [3,3,3])
        #expect(constant.fitted == [3,3,3] && constant.emptyWeightAnchors == 0)
        #expect(throws: (any Error).self) { try VivoOmicsRobustLowess.fit(x: x,y: y,maximumNeighborhoodVisits: 1) }
        #expect(throws: (any Error).self) { try VivoOmicsRobustLowess.fit(x: [0,.infinity],y: [1,2]) }
    }
    @Test func perfectCountFitsRetainScaleFloorAndTwoUpdates() throws {
        let counts = (1...12).map { [UInt64](repeating: UInt64($0),count: 4) }
        let result = try VivoOmicsNBQLGlobalScale.fit(counts: counts,design: [[1],[1],[1],[1]],offsets: [0,0,0,0],
            contrast: [1],trendDispersions: [Double](repeating: 0.1,count: 12),abundanceCovariates: (1...12).map(Double.init))
        #expect(result.completed && result.failures.isEmpty)
        #expect(result.averageQuasiDispersion == 1)
        #expect(result.updates.count == 2)
        #expect(result.updates[0].inputScale == 1 && result.updates[1].inputScale == result.updates[0].outputScale)
        #expect(result.updates.allSatisfy { $0.omittedIndices.isEmpty })
        #expect(result.refittedFits.allSatisfy { $0?.dispersion == 0.1 })
    }
    @Test func heterogeneousCountsUseFixedInitialMeansAndExposeFamilyFailures() throws {
        let counts: [[UInt64]] = (1...24).map { (g: Int) -> [UInt64] in
            let value = UInt64(g)
            return [value,0,5*value,2*value,1,10*value]
        }
        let design = [[1.0,0],[1,0],[1,0],[1,1],[1,1],[1,1]], abundance = (1...24).map { log(Double($0)) }
        let dispersion = [Double](repeating: 0.05,count: 24)
        let result = try VivoOmicsNBQLGlobalScale.fit(counts: counts,design: design,offsets: [0,0,0,0,0,0],contrast: [0,1],trendDispersions: dispersion,abundanceCovariates: abundance)
        #expect(result.completed)
        let scale = try #require(result.averageQuasiDispersion)
        #expect(scale > 1)
        for update in result.updates {
            let first = try VivoOmicsNBResidualAdjustment.adjustedResiduals(counts: counts[0],means: result.initialFits[0]!.means,design: design,dispersion: 0.05,averageQuasiDispersion: update.inputScale,method: .adaptive)
            #expect(abs(first.deviance-update.residualDeviances[0]!) < 1e-10)
        }
        #expect(result.refittedFits.allSatisfy { $0!.converged && abs($0!.dispersion-0.05/scale) < 1e-14 })
        var invalid = counts;invalid[3] = [0,0,0,0,0,0]
        let failed = try VivoOmicsNBQLGlobalScale.fit(counts: invalid,design: design,offsets: [0,0,0,0,0,0],contrast: [0,1],trendDispersions: dispersion,abundanceCovariates: abundance)
        #expect(!failed.completed && failed.averageQuasiDispersion == nil && failed.updates.isEmpty)
        #expect(failed.failures.contains { $0.featureIndex == 3 && $0.stage == "initialFit" })
        #expect(failed.initialFits.count == 24)
        let limited = try VivoOmicsNBQLGlobalScale.fit(counts: counts,design: design,offsets: [0,0,0,0,0,0],contrast: [0,1],trendDispersions: dispersion,abundanceCovariates: abundance,maximumMomentEvaluations: 1)
        #expect(!limited.completed && !limited.failures.isEmpty && limited.refittedFits.allSatisfy { $0 == nil })
    }
}
