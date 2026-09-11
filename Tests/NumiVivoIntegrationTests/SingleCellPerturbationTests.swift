import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellPerturbationTests {
    static func fixture(constantControls: Bool = false) throws -> (VivoH5ADPseudobulkReport,VivoPerturbationPlan) {
        let original=try VivoSingleCellExamples.pairedCounts()
        var counts=original.matrix.counts
        if constantControls {
            for row in original.cells.indices where original.cells[row].sampleID.hasSuffix("-0") {
                for k in original.matrix.rowOffsets[row]..<original.matrix.rowOffsets[row+1] { counts[k]=100 }
            }
        }
        let matrix=VivoSparseCounts(cellCount: original.cells.count,featureCount: original.features.count,rowOffsets: original.matrix.rowOffsets,featureIndices: original.matrix.featureIndices,counts: counts)
        let d=VivoSingleCellDataset(id: original.id,evidence: original.evidence,sourceDescription: original.sourceDescription,countUnit: original.countUnit,samples: original.samples,features: original.features,cells: original.cells,matrix: matrix)
        let metadata=VivoSingleCellCountMetadata(id: d.id,evidence: d.evidence,sourceDescription: d.sourceDescription,countUnit: d.countUnit,samples: d.samples,features: d.features,cells: d.cells)
        let report=try VivoH5ADPseudobulkReport(schemaVersion: 1,method: "synthetic-control",metadata: metadata,quality: VivoSingleCellAnalysis.quality(d),pseudobulk: VivoSingleCellAnalysis.pseudobulk(d),canonicalNonzeros: d.matrix.counts.count,hdf5Version: "fixture",contrasts: [])
        let mapping=VivoH5ADImportPlan(id: "fixture",evidence: .synthetic,sourceDescription: "control fixture",countUnit: .umiCount,matrixPath: "X",samples: d.samples,sampleColumn: "sample",groupColumn: "group")
        let plan=VivoPerturbationPlan(mapping: mapping,featureNamespace: "synthetic",perturbationID: "synthetic-treatment",controlCondition: "control",treatmentCondition: "treated",provenance: "Synthetic numerical control, not biological qualification")
        return (report,plan)
    }
    @Test func pairedModelUsesTrainingDataAndReplays() throws {
        let (report,plan)=try Self.fixture(),source=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
        let a=try VivoPerturbation.model(report,plan: plan,source: source),b=try VivoPerturbation.model(report,plan: plan,source: source)
        #expect(a==b && a.trainingDonors.count==6 && a.meanResponse.count==32)
        #expect(a.maximumSolveResidual<1e-10)
        #expect(a.meanResponse[0]>0 && a.meanResponse[1]>0 && a.meanResponse[2]>0)
        let (_,logs,_)=try VivoPerturbation.dense(report.pseudobulk)
        for j in a.featureIDs.indices {
            var values: [Double]=[]
            for donor in a.trainingDonors {
                let rows=report.pseudobulk.groups
                let control=try #require(rows.firstIndex { $0.donorID==donor && $0.condition=="control" })
                let treated=try #require(rows.firstIndex { $0.donorID==donor && $0.condition=="treated" })
                values.append(logs[treated][j]-logs[control][j])
            }
            #expect(abs(values.reduce(0,+)/6-a.meanResponse[j])<1e-12)
        }
    }
    @Test func constantControlsDoNotInventContextVariation() throws {
        let (report,plan)=try Self.fixture(constantControls: true)
        let result=try VivoPerturbation.model(report,plan: plan,source: VivoFingerprint(bytes: Array(repeating: 0,count: 32)))
        #expect(result.contextScales.allSatisfy { $0==1 })
        #expect(result.contexts.allSatisfy { $0.allSatisfy { $0==0 } })
        #expect(result.maximumSolveResidual<1e-12)
    }
    @Test func unknownPredictionSettingsAreRejected() throws {
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoPerturbationQueryPlan.self,from: Data("{\"schemaVersion\":1,\"useTreated\":true}".utf8)) }
        let (_,plan)=try Self.fixture()
        let wrong=VivoPerturbationPlan(mapping: plan.mapping,featureNamespace: plan.featureNamespace,perturbationID: plan.perturbationID,controlCondition: "control",treatmentCondition: "control",provenance: plan.provenance)
        #expect(throws: (any Error).self) { try wrong.validate() }
    }
    @Test func explicitPanelKeepsFullLibraryDenominators() throws {
        let (report,original)=try Self.fixture(),source=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
        let ids=Array(report.pseudobulk.featureIDs.prefix(3).reversed())
        let plan=VivoPerturbationPlan(mapping: original.mapping,featureNamespace: original.featureNamespace,
            perturbationID: original.perturbationID,controlCondition: original.controlCondition,
            treatmentCondition: original.treatmentCondition,provenance: original.provenance,responseFeatureIDs: ids)
        let full=try VivoPerturbation.model(report,plan: original,source: source)
        let panel=try VivoPerturbation.model(report,plan: plan,source: source)
        #expect(panel.featureIDs==ids)
        #expect(panel.meanResponse==Array(full.meanResponse.prefix(3).reversed()))
        #expect(try VivoCanonicalJSON.decode(VivoPerturbationPlan.self,from: VivoCanonicalJSON.encode(plan))==plan)
        let group=report.pseudobulk.groups.first { $0.condition==original.controlCondition }!
        let queryGroup=VivoPseudobulkGroup(biologicalReplicateID: "query",donorID: "query",condition: group.condition,
            organism: group.organism,cellGroup: group.cellGroup,sampleIDs: ["query"],batchIDs: [],sourceCellIndices: [0])
        let queryPlan=VivoPerturbationQueryPlan(mapping: .init(id: "query",evidence: .synthetic,sourceDescription: "query control",
            countUnit: .umiCount,matrixPath: "X",samples: [],sampleColumn: "sample"),
            featureNamespace: original.featureNamespace,perturbationID: original.perturbationID)
        func query(_ features: [String],_ counts: [UInt64]) -> VivoPseudobulkCounts {
            .init(method: "synthetic-query",countUnit: .umiCount,groups: [queryGroup],featureIDs: features,
                matrix: .init(cellCount: 1,featureCount: features.count,rowOffsets: [0,counts.count],featureIndices: Array(counts.indices),counts: counts))
        }
        let a=try VivoPerturbation.evaluate(query(ids+["extra"],[10,20,30,40]),plan: queryPlan,model: panel)
        let b=try VivoPerturbation.evaluate(query(ids+["extra"],[10,20,30,940]),plan: queryPlan,model: panel)
        #expect(a.predictions[0].libraryCounts==100 && b.predictions[0].libraryCounts==1000)
        #expect(abs(a.predictions[0].control[0]-log1p(100_000.0))<1e-12)
        #expect(abs(b.predictions[0].control[0]-log1p(10_000.0))<1e-12)
        let reordered=try VivoPerturbation.evaluate(query(["extra"]+ids.reversed(),[40,30,20,10]),plan: queryPlan,model: panel)
        #expect(reordered==a)
        #expect(throws: (any Error).self) { try VivoPerturbation.evaluate(query(Array(ids.dropLast()),[10,20]),plan: queryPlan,model: panel) }
        #expect(throws: (any Error).self) { try VivoPerturbation.evaluate(query(ids+["extra"],[10,20,30,40]),plan: queryPlan,model: full) }
    }
    @Test func responsePanelRejectsInvalidAndUnmeasuredFeatures() throws {
        let (report,original)=try Self.fixture()
        for ids in [[],["missing"],[report.pseudobulk.featureIDs[0],report.pseudobulk.featureIDs[0]]] {
            let plan=VivoPerturbationPlan(mapping: original.mapping,featureNamespace: original.featureNamespace,
                perturbationID: original.perturbationID,controlCondition: original.controlCondition,
                treatmentCondition: original.treatmentCondition,provenance: original.provenance,responseFeatureIDs: ids)
            #expect(throws: (any Error).self) { try VivoPerturbation.model(report,plan: plan,source: VivoFingerprint(bytes: Array(repeating: 0,count: 32))) }
        }
        // An omitted panel must retain the historical serialized plan shape.
        let object=try JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(original)) as! [String: Any]
        #expect(object["responseFeatureIDs"]==nil)
    }
}
