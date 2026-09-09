import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellNBSupportTests {
    static func request() -> VivoOmicsExpressionContrast {
        var r=VivoOmicsExpressionContrast(id: "boundary",controlCondition: "ctrl",treatmentCondition: "stim",design: .pairedDonors)
        r.model = .negativeBinomial;r.minimumCellsPerPseudobulk=1
        var options=VivoOmicsNBCohortOptions();options.zeroTotalDonorPolicy = .activeDonorProfile;options.trend = .mean
        r.negativeBinomialOptions=options;return r
    }
    static func design(batch: Bool = false) throws -> (VivoOmicsDesignMatrix,VivoOmicsExpressionContrast) {
        var request=Self.request();request.adjustForBatch=batch
        let observations=(0..<8).map { row in
            VivoPseudobulkGroup(biologicalReplicateID: "d\(row/2)",donorID: "d\(row/2)",condition: row%2==0 ? "ctrl" : "stim",organism: "fixture",
                cellGroup: nil,sampleIDs: ["s\(row)"],batchIDs: [batch && row>1 && row%2==1 ? "B" : "A"],sourceCellIndices: [row])
        }
        let b=try VivoPseudobulkDifferentialExpression.makeDesign(observations: observations,request: request)
        return (.init(columnNames: b.columnNames,rows: b.rows,contrast: b.contrast,sourcePseudobulkIndices: Array(0..<8),observations: observations,
            excludedSmallPseudobulkIndices: [],controlReplicates: 4,treatmentReplicates: 4,residualDegreesOfFreedom: 8-b.columnNames.count,
            sizeFactorValues: Array(repeating: 1,count: 8),libraryCounts: Array(repeating: 1000,count: 8),referenceFeatureIndices: []),request)
    }
    @Test func baselineDonorRemovalRebasesAndMatchesAnalyticInformationLimit() throws {
        let (d,r)=try Self.design(),y: [UInt64]=[0,0,20,40,30,60,50,100]
        let resolved=try VivoOmicsNBSupport.resolve(counts: y,design: d,request: r)
        let s=try #require(resolved)
        #expect(s.outcome == .ready && s.excludedObservationIndices==[0,1])
        #expect(s.columnNames==["intercept","treatment-minus-control","donor:d2","donor:d3"])
        #expect(s.residualDegreesOfFreedom==2 && s.activeDonorIDs==["d1","d2","d3"])
        let fit=try VivoOmicsNegativeBinomial.fit(counts: Array(y.dropFirst(2)),design: s.rows!,offsets: Array(repeating: 0,count: 6),contrast: s.contrast!,dispersion: 0.2)
        #expect(fit.converged)
        let effectError=abs(fit.effect!-log(2))
        // The fitter stops at maximum scaled score <= 1e-7; compare the
        // analytic effect in standard-error units at that numerical scale.
        #expect(effectError/fit.standardError!<1e-7)
        #expect(fit.maximumScaledScore<=1e-7)
        let information=[20.0,30,50].reduce(0.0) { sum,m in
            let w0=m/(1+0.2*m),w1=2*m/(1+0.4*m);return sum+w0*w1/(w0+w1)
        }
        #expect(abs(fit.standardError!-sqrt(1/information))<1e-10)
        // Zero-total pair likelihood and Schur-complement treatment information
        // vanish as its nuisance rate approaches zero. No finite zero-donor
        // coefficient, full-design determinant or sampling calibration is claimed.
        for t in [10.0,20,30] {
            let mu=exp(-t),w0=mu/(1+0.2*mu),w1=2*mu/(1+0.4*mu)
            let zeroLL=try VivoOmicsNegativeBinomial.logMass(count: 0,mean: mu,dispersion: 0.2)+VivoOmicsNegativeBinomial.logMass(count: 0,mean: 2*mu,dispersion: 0.2)
            #expect(abs(zeroLL)<=3*mu)
            #expect(w0*w1/(w0+w1)<mu)
        }
    }
    @Test func replicationSeparationAndBatchConfoundingRemainRejected() throws {
        let (d,r)=try Self.design()
        #expect(try VivoOmicsNBSupport.resolve(counts: [0,0,0,0,5,10,8,16],design: d,request: r)?.outcome == .insufficientReplication)
        #expect(try VivoOmicsNBSupport.resolve(counts: [0,0,0,10,0,16,0,12],design: d,request: r)?.outcome == .rankDeficientPositiveSupport)
        #expect(try VivoOmicsNBSupport.resolve(counts: [0,1,5,10,8,16,6,12],design: d,request: r) == nil)
        let (batch,br)=try Self.design(batch: true)
        #expect(try VivoOmicsNBSupport.resolve(counts: [0,0,5,10,8,16,6,12],design: batch,request: br)?.outcome == .rankDeficientDesign)
        var wrong=r;wrong.design = .independentReplicates
        #expect(throws: (any Error).self) { try wrong.validate() }
    }
    @Test func cohortPreservesPriorAndFullSupportFitsWhileRecordingNewDesign() throws {
        let source=try SingleCellNBCohortTests.fixture();var offsets=[0],indices: [Int]=[],counts: [UInt64]=[]
        for row in source.cells.indices {
            for k in source.matrix.rowOffsets[row]..<source.matrix.rowOffsets[row+1] {
                let gene=source.matrix.featureIndices[k]
                if gene==90 && row<2 || gene==91 && row<8 || gene==92 && row%2==0 { continue }
                indices.append(gene);counts.append(source.matrix.counts[k])
            }
            offsets.append(counts.count)
        }
        let data=VivoSingleCellDataset(id: source.id,evidence: source.evidence,sourceDescription: source.sourceDescription,countUnit: source.countUnit,
            samples: source.samples,features: source.features,cells: source.cells,
            matrix: .init(cellCount: source.cells.count,featureCount: source.features.count,rowOffsets: offsets,featureIndices: indices,counts: counts))
        let request=Self.request();var baseline=request;baseline.negativeBinomialOptions!.zeroTotalDonorPolicy=nil
        let old=try VivoPseudobulkDifferentialExpression.run(data,contrast: baseline)
        let new=try VivoPseudobulkDifferentialExpression.run(data,contrast: request)
        #expect(old.features[90].status == .rankDeficientSupport)
        #expect(new.negativeBinomial!.features[90].supportResolution?.outcome == .ready)
        #expect(new.negativeBinomial!.features[90].finalFit != nil)
        #expect(new.features[91].status == .insufficientActiveDonors && new.features[91].pValue == nil)
        #expect(new.features[92].status == .rankDeficientSupport && new.features[92].pValue == nil)
        #expect(old.negativeBinomial!.trend==new.negativeBinomial!.trend)
        #expect(old.negativeBinomial!.features[0]==new.negativeBinomial!.features[0])
        #expect(old.features[0].pValue==new.features[0].pValue)
    }
}
