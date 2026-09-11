import Foundation
import Testing
@testable import NumiVivoKit
@Suite struct SingleCellCountCalibrationTests {
    static func group(_ rows: [[UInt64]]) throws -> VivoCellCountMoments {
        let a=try VivoCellCountMomentAccumulator(featureCount: rows[0].count)
        for row in rows {
            let ids=row.indices.filter { row[$0]>0 }
            try a.appendCell(featureIndices: ids,counts: ids.map { row[$0] })
        }
        return try a.finish()
    }
    @Test func sparseZerosAndUnequalDepthsMatchDirectMoments() throws {
        let rows: [[UInt64]]=[[1,9,0],[0,20,0],[6,24,0]]
        let g=try Self.group(rows)
        #expect(g.cells==3 && g.libraryCounts==60 && g.counts==[7,53,0])
        let z=[100_000.0,0,200_000.0]
        let v=z.reduce(0) { $0+pow($1-100_000,2) }/2
        #expect(g.meanCPM[0]==100_000 && abs(g.sampleVarianceCPM[0]-v)<1e-5)
        #expect(abs(g.distinctCellRateProductCPM2[0]-2e10/3)<1e-5)
        #expect(g.positiveCells==[2,3,0] && g.sampleVarianceCPM[2]==0)
    }
    @Test func invalidCellDoesNotPartiallyMutateAccumulator() throws {
        let a=try VivoCellCountMomentAccumulator(featureCount: 2)
        try a.appendCell(featureIndices: [0,1],counts: [1,9])
        #expect(throws: (any Error).self) { try a.appendCell(featureIndices: [1,0],counts: [2,8]) }
        try a.appendCell(featureIndices: [0,1],counts: [1,9])
        #expect(try a.finish()==Self.group([[1,9],[1,9]]))
    }
    @Test func calibrationKeepsDonorVariationOutOfCellDispersion() throws {
        let groups=try [1,4,8].map { v in try Self.group(Array(repeating: [UInt64(v),UInt64(10-v),0],count: 10)) }
        let model=try VivoCountObservationCalibration.fit(groups,featureIDs: ["gene","other","zero"],donorIDs: ["a","b","c"],conditionID: "control",trainingSource: VivoFingerprint(bytes: Array(repeating: 0,count: 32)))
        let f=model.features[0]
        #expect(f.cellDispersion==0 && f.poissonBoundary && f.latentDonorRateVarianceCPM!>0)
        #expect(f.status=="availableConditionalMomentCalibration")
        #expect(model.features[2].status=="insufficientWithinDonorCountPairs")
        let p=try VivoCountObservationCalibration.posterior(counts: [0,0],libraryCounts: [10,10],featureID: "gene",queryDonorID: "query",conditionID: "control",model: model)
        #expect(p.varianceCPM>0)
        #expect(throws: (any Error).self) { try VivoCountObservationCalibration.posterior(counts: [0,0],libraryCounts: [10,10],featureID: "gene",queryDonorID: "a",conditionID: "control",model: model) }
        #expect(throws: (any Error).self) { try VivoCountObservationCalibration.posterior(counts: [0,0],libraryCounts: [10,10],featureID: "gene",queryDonorID: "query",conditionID: "treated",model: model) }
    }
    @Test func zeroLatentVarianceStillHasTrainingMeanUncertainty() throws {
        let g=try Self.group([[1,9],[1,9],[1,9]])
        let m=try VivoCountObservationCalibration.fit([g,g,g],featureIDs: ["gene","other"],donorIDs: ["a","b","c"],conditionID: "control",trainingSource: VivoFingerprint(bytes: Array(repeating: 0,count: 32)))
        let f=m.features[0]
        #expect(f.latentDonorRateVarianceCPM==0 && f.latentVarianceBoundary)
        #expect(f.newDonorRateVarianceCPM!>0 && f.gammaPriorShape != nil)
        #expect(abs(f.newDonorRateVarianceCPM!-f.meanDonorMeasurementVarianceCPM!/3)<1e-5)
    }
}
