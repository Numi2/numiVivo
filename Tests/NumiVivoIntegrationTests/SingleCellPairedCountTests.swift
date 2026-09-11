import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellPairedCountTests {
    let source=try! VivoFingerprint(bytes: Array(repeating: 7,count: 32))
    func group(_ donor: Int,_ condition: String,_ rows: [[UInt64]]) throws -> VivoPairedCountGroup {
        let a=try VivoCellCountMomentAccumulator(featureCount: 2)
        for row in rows {
            let ids=row.indices.filter { row[$0]>0 }
            try a.appendCell(featureIndices: ids,counts: ids.map { row[$0] })
        }
        return .init(donorID: "d\(donor)",conditionID: condition,moments: try a.finish())
    }
    func analyze(_ groups: [VivoPairedCountGroup]) throws -> VivoPairedCountReport {
        try VivoPairedCountMoments.analyze(groups: groups,featureIDs: ["g","background"],
            controlConditionID: "control",treatedConditionID: "treated",trainingSource: source)
    }
    @Test func pairingUsesIdentityAndRejectsMissingOrDuplicateStrata() throws {
        var groups: [VivoPairedCountGroup]=[]
        for i in 0..<4 {
            groups.append(try group(i,"control",[[UInt64(10+i*20),90],[UInt64(12+i*20),88]]))
            groups.append(try group(i,"treated",[[UInt64(80-i*20),20],[UInt64(78-i*20),22]]))
        }
        #expect(try analyze(groups)==analyze(groups.reversed()))
        #expect(throws: (any Error).self) { try analyze(Array(groups.dropLast())) }
        var duplicated=groups;duplicated[duplicated.count-1]=groups[1]
        #expect(throws: (any Error).self) { try analyze(duplicated) }
    }
    @Test func diagonalClippingDoesNotRepairAnImpossibleJointCovariance() throws {
        var groups: [VivoPairedCountGroup]=[]
        for (i,x) in [10,30,60,90].enumerated() {
            let row=[UInt64(x),UInt64(100-x)]
            groups.append(try group(i,"control",[row,row]))
            groups.append(try group(i,"treated",[row,row]))
        }
        let f=try analyze(groups).features[0]
        #expect(f.rawLatentControlVarianceCPM2!>0)
        #expect(f.rawCovarianceStatus=="indefinite")
        #expect(f.marginalClippedCovarianceStatus=="indefinite")
        #expect(f.observedResponseVarianceCPM2==0)
        #expect(f.rawLatentResponseVarianceCPM2!<0)
    }
    @Test func negativeAssociationIsRetainedAndEndpointUsesFullLibraries() throws {
        var groups: [VivoPairedCountGroup]=[]
        for (i,x) in [10,30,60,90].enumerated() {
            groups.append(try group(i,"control",[[UInt64(x),UInt64(100-x)],[UInt64(x),UInt64(300-x)]]))
            groups.append(try group(i,"treated",[[UInt64(100-x),UInt64(x)],[UInt64(100-x),UInt64(x)]]))
        }
        let f=try analyze(groups).features[0]
        #expect(f.observedCovarianceCPM2<0)
        #expect(abs(f.meanCellRateResponseCPM-208333.33333333334)<1e-8)
        #expect(abs(f.meanPseudobulkResponseCPM-287500)<1e-8)
        #expect(abs(f.maximumAbsoluteEndpointDifferenceCPM-150000)<1e-8)
        #expect(abs(f.observedResponseVarianceCPM2-(f.observedControlVarianceCPM2+f.observedTreatedVarianceCPM2-2*f.observedCovarianceCPM2))<1e-3)
    }
    @Test func noCountSupportRemainsUnavailable() throws {
        var groups: [VivoPairedCountGroup]=[]
        for i in 0..<3 { for condition in ["control","treated"] { groups.append(try group(i,condition,[[0,100],[0,100]])) } }
        let f=try analyze(groups).features[0]
        #expect(f.rawCovarianceStatus=="unavailableCellDispersion")
        #expect(f.rawMinimumEigenvalueCPM2==nil)
        #expect(f.controlCalibrationStatus=="insufficientWithinDonorCountPairs")
    }
    @Test func independentPairedRatesHavePositiveJointVariance() throws {
        var groups: [VivoPairedCountGroup]=[]
        for (i,pair) in [(10,10),(10,90),(90,10),(90,90)].enumerated() {
            for (condition,x) in [("control",pair.0),("treated",pair.1)] {
                let row=[UInt64(x),UInt64(100-x)]
                groups.append(try group(i,condition,[row,row]))
            }
        }
        let f=try analyze(groups).features[0]
        #expect(f.observedCovarianceCPM2==0)
        #expect(f.rawCovarianceStatus=="positiveSemidefinite")
        #expect(abs(f.rawMinimumEigenvalueCPM2!-210833333333.33334)<1e-3)
        #expect(f.meanLog1pCellRateResponse==f.meanLog1pPseudobulkResponse)
    }
}
