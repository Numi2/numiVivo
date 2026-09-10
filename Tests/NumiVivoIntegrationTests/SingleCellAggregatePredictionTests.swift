import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellAggregatePredictionTests {
    static func inputs() throws -> (VivoH5ADPseudobulkReport,VivoPerturbationPlan,VivoPerturbationAggregateFold,VivoPerturbationAggregateBatchPlan,Int) {
        let (report,original)=try SingleCellPerturbationTests.fixture(),groups=report.pseudobulk.groups
        let donor=try #require(groups.first?.donorID)
        let train=groups.indices.filter { groups[$0].donorID != donor }
        let query=groups.indices.filter { groups[$0].donorID==donor && groups[$0].condition==original.controlCondition }
        let target=try #require(groups.firstIndex { $0.donorID==donor && $0.condition==original.treatmentCondition })
        let fold=VivoPerturbationAggregateFold(id: "held-out",perturbationID: original.perturbationID,
            controlCondition: original.controlCondition,treatmentCondition: original.treatmentCondition,
            cellGroup: "synthetic-population",trainingGroupIndices: train,queryGroupIndices: query)
        let batch=try VivoPerturbationAggregateBatchPlan(sourceReport: VivoFingerprint(bytes: Array(repeating: 0,count: 32)),
            featureNamespace: "synthetic",provenance: "Synthetic isolation regression, not biological evidence",folds: [fold])
        return (report,original,fold,batch,target)
    }
    static func changedCounts(_ bulk: VivoPseudobulkCounts,row: Int) -> VivoPseudobulkCounts {
        var counts=bulk.matrix.counts;counts[bulk.matrix.rowOffsets[row]]+=100_000
        return .init(method: bulk.method,countUnit: bulk.countUnit,groups: bulk.groups,featureIDs: bulk.featureIDs,
            matrix: .init(cellCount: bulk.matrix.cellCount,featureCount: bulk.matrix.featureCount,rowOffsets: bulk.matrix.rowOffsets,
                          featureIndices: bulk.matrix.featureIndices,counts: counts))
    }
    @Test func heldOutTreatedCountsCannotChangeModelOrPrediction() throws {
        let (report,original,fold,batch,target)=try Self.inputs()
        let a=try VivoPerturbationAggregateBatch.evaluateFold(report.pseudobulk,mapping: original.mapping,batch: batch,fold: fold)
        let b=try VivoPerturbationAggregateBatch.evaluateFold(Self.changedCounts(report.pseudobulk,row: target),mapping: original.mapping,batch: batch,fold: fold)
        #expect(a.0.status=="completed" && a.0==b.0 && a.1==b.1 && a.2==b.2)
        let model=try VivoCanonicalJSON.decode(VivoPerturbationModel.self,from: #require(a.1))
        #expect(model.trainingDonors.count==5)
        #expect(!model.trainingDonors.contains(report.pseudobulk.groups[target].donorID!))
        let (training,query)=try VivoPerturbationAggregateBatch.isolatedInputs(report.pseudobulk,fold: fold)
        for (local,source) in fold.trainingGroupIndices.enumerated() {
            #expect(training.groups[local].sourceCellIndices==report.pseudobulk.groups[source].sourceCellIndices)
            let range=report.pseudobulk.matrix.rowOffsets[source]..<report.pseudobulk.matrix.rowOffsets[source+1]
            let projected=training.matrix.rowOffsets[local]..<training.matrix.rowOffsets[local+1]
            #expect(Array(training.matrix.counts[projected])==Array(report.pseudobulk.matrix.counts[range]))
        }
        #expect(query.groups.count==1 && query.groups[0].sourceCellIndices==report.pseudobulk.groups[fold.queryGroupIndices[0]].sourceCellIndices)
    }
    @Test func queryChangesPredictionButCannotRefitTraining() throws {
        let (report,original,fold,batch,_)=try Self.inputs()
        let a=try VivoPerturbationAggregateBatch.evaluateFold(report.pseudobulk,mapping: original.mapping,batch: batch,fold: fold)
        let b=try VivoPerturbationAggregateBatch.evaluateFold(Self.changedCounts(report.pseudobulk,row: fold.queryGroupIndices[0]),mapping: original.mapping,batch: batch,fold: fold)
        #expect(a.0.status=="completed" && b.0.status=="completed")
        #expect(a.1==b.1 && a.0.trainingAggregate==b.0.trainingAggregate)
        #expect(a.2 != b.2 && a.0.queryAggregate != b.0.queryAggregate)
        let c=try VivoPerturbationAggregateBatch.evaluateFold(Self.changedCounts(report.pseudobulk,row: fold.trainingGroupIndices[0]),mapping: original.mapping,batch: batch,fold: fold)
        #expect(c.0.status=="completed" && a.1 != c.1)
    }
    @Test func treatedQueriesAndDonorOverlapFailBeforeModelFitting() throws {
        let (report,original,fold,batch,target)=try Self.inputs()
        let treated=VivoPerturbationAggregateFold(id: "treated-query",perturbationID: fold.perturbationID,controlCondition: fold.controlCondition,
            treatmentCondition: fold.treatmentCondition,cellGroup: fold.cellGroup,trainingGroupIndices: fold.trainingGroupIndices,queryGroupIndices: [target])
        let bad=try VivoPerturbationAggregateBatch.evaluateFold(report.pseudobulk,mapping: original.mapping,batch: batch,fold: treated)
        #expect(bad.0.status=="failed" && bad.0.trainingAggregate==nil && bad.1==nil && bad.2==nil)
        let overlap=VivoPerturbationAggregateFold(id: "overlap",perturbationID: fold.perturbationID,controlCondition: fold.controlCondition,
            treatmentCondition: fold.treatmentCondition,cellGroup: fold.cellGroup,trainingGroupIndices: fold.trainingGroupIndices+[target],queryGroupIndices: fold.queryGroupIndices)
        let donor=try VivoPerturbationAggregateBatch.evaluateFold(report.pseudobulk,mapping: original.mapping,batch: batch,fold: overlap)
        #expect(donor.0.status=="failed" && donor.0.trainingAggregate==nil && donor.1==nil && donor.2==nil)
    }
    @Test func relabelingCannotHideSharedSourceCellsAndUnknownSettingsFail() throws {
        let (report,original,fold,batch,_)=try Self.inputs();var groups=report.pseudobulk.groups
        let row=fold.queryGroupIndices[0],g=groups[row]
        groups[row] = .init(biologicalReplicateID: g.biologicalReplicateID,donorID: g.donorID,condition: g.condition,organism: g.organism,
            cellGroup: g.cellGroup,sampleIDs: g.sampleIDs,batchIDs: g.batchIDs,sourceCellIndices: groups[fold.trainingGroupIndices[0]].sourceCellIndices)
        let source=VivoPseudobulkCounts(method: report.pseudobulk.method,countUnit: report.pseudobulk.countUnit,
            groups: groups,featureIDs: report.pseudobulk.featureIDs,matrix: report.pseudobulk.matrix)
        let result=try VivoPerturbationAggregateBatch.evaluateFold(source,mapping: original.mapping,batch: batch,fold: fold)
        #expect(result.0.status=="failed" && result.0.trainingAggregate==nil)
        var json=try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(batch)) as? [String:Any])
        json["useHeldOutTreated"]=true
        #expect(throws: (any Error).self) {
            try VivoCanonicalJSON.decode(VivoPerturbationAggregateBatchPlan.self,from: JSONSerialization.data(withJSONObject: json))
        }
    }
}
