import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellCompositionTests {
    static let context=VivoCompositionContext(id: "synthetic-pool",organism: "NCBITaxon:9606",featureNamespace: "synthetic",perturbationNamespace: "synthetic-target",countUnit: .umiCount)
    static func fixture(pair: [UInt64] = [9,99,9,0]) throws -> VivoCompositionTraining {
        let rows: [[UInt64]]=[[10,30,0,0],[20,20,1,0],[30,10,0,0],pair]
        let names=["control","a","b","pair"]
        let samples=names.map { VivoOmicsSample(id: $0,biologicalReplicateID: "pool",condition: $0,batchID: "unreported",organism: context.organism) }
        var offsets=[0],indices: [Int]=[],values: [UInt64]=[]
        for row in rows {
            for j in row.indices where row[j]>0 { indices.append(j);values.append(row[j]) };offsets.append(values.count)
        }
        let d=VivoSingleCellDataset(id: "synthetic",evidence: .synthetic,sourceDescription: "Numerical control, not biological evidence",countUnit: .umiCount,
            samples: samples,features: (0..<4).map { .init(id: "g\($0)",name: "g\($0)",mitochondrial: false) },
            cells: names.map { .init(barcode: $0,sampleID: $0,group: "synthetic") },
            matrix: .init(cellCount: 4,featureCount: 4,rowOffsets: offsets,featureIndices: indices,counts: values))
        let metadata=VivoSingleCellCountMetadata(id: d.id,evidence: d.evidence,sourceDescription: d.sourceDescription,countUnit: d.countUnit,samples: d.samples,features: d.features,cells: d.cells)
        let report=try VivoH5ADPseudobulkReport(schemaVersion: 1,method: "synthetic-control",metadata: metadata,quality: VivoSingleCellAnalysis.quality(d),
            pseudobulk: VivoSingleCellAnalysis.pseudobulk(d),canonicalNonzeros: d.matrix.counts.count,hdf5Version: "fixture",contrasts: [])
        let plan=VivoCompositionPreparation(context: context,controlCondition: "control",targets: [.init(id: "b",condition: "b"),.init(id: "a",condition: "a")],provenance: "Synthetic numerical control")
        let hash=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
        return try VivoComposition.select(report,plan: plan,source: hash,sourceReport: hash)
    }
    @Test func trainingExcludesPairedOutcomesAndCompositionReplays() throws {
        let training=try Self.fixture(),altered=try Self.fixture(pair: [999,1,99,2])
        #expect(training==altered && training.matrix.cellCount==3)
        let model=try VivoComposition.fit(training)
        #expect(model.targetIDs==["a","b"] && model.controlCPM==[250000,750000,0,0])
        let plan=VivoCompositionQueryPlan(context: Self.context,queries: [.init(id: "ab",targets: ["a","b"]),.init(id: "ba",targets: ["b","a"])])
        let prediction=try VivoComposition.predict(model,plan: plan)
        #expect(prediction.methods.count==6)
        for method in prediction.methods {
            #expect(method.expression[method.queryRows[0]]==method.expression[method.queryRows[1]])
            #expect(method.diagnostics[0]==method.diagnostics[1])
        }
        let cpm=try #require(prediction.methods.first { $0.method == .additiveCPMResponse })
        #expect(cpm.diagnostics[0].clippedFeatures==1 && cpm.expression[0][1]==0)
        #expect(abs(cpm.expression[0][0]-log1p(model.targetCPM[0][0]+model.targetCPM[1][0]-model.controlCPM[0]))<1e-12)
        let averaged=try #require(prediction.methods.first { $0.method == .meanConstituentLogResponse })
        #expect(abs(averaged.expression[0][0]-(log1p(model.targetCPM[0][0])+log1p(model.targetCPM[1][0]))/2)<1e-12)
    }
    @Test func unknownSettingsAndTargetsAreRejectedBeforePrediction() throws {
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoCompositionQuery.self,from: Data("{\"id\":\"q\",\"targets\":[\"a\",\"b\"],\"pairedCounts\":[1]}".utf8)) }
        let model=try VivoComposition.fit(Self.fixture())
        for targets in [["a","unknown"],["a","a"],["a"],["a","b","c"]] {
            let plan=VivoCompositionQueryPlan(context: Self.context,queries: [.init(id: "q",targets: targets)])
            #expect(throws: (any Error).self) { try VivoComposition.predict(model,plan: plan) }
        }
        let wrong=VivoCompositionContext(id: "other-pool",organism: Self.context.organism,featureNamespace: Self.context.featureNamespace,perturbationNamespace: Self.context.perturbationNamespace,countUnit: .umiCount)
        #expect(throws: (any Error).self) { try VivoComposition.predict(model,plan: .init(context: wrong,queries: [.init(id: "q",targets: ["a","b"])])) }
    }
    @Test func queryStorageBudgetRejectsBeforeAllocatingOutput() throws {
        let original=try Self.fixture()
        let training=VivoCompositionTraining(selection: original.selection,source: original.source,sourceReport: original.sourceReport,evidence: .synthetic,
            featureIDs: (0..<20_000).map { "g\($0)" },matrix: .init(cellCount: 3,featureCount: 20_000,rowOffsets: [0,1,2,3],featureIndices: [0,0,0],counts: [1,1,1]))
        let model=try VivoComposition.fit(training)
        let plan=VivoCompositionQueryPlan(context: Self.context,queries: (0..<256).map { .init(id: "q\($0)",targets: ["a","b"]) })
        #expect(throws: (any Error).self) { try VivoComposition.predict(model,plan: plan) }
    }
    @Test func emptyInexactAndOverflowingLibrariesAreRejected() throws {
        let original=try Self.fixture()
        let matrices: [VivoSparseCounts]=[
            .init(cellCount: 3,featureCount: 4,rowOffsets: [0,0,1,2],featureIndices: [0,0],counts: [1,1]),
            .init(cellCount: 3,featureCount: 4,rowOffsets: [0,1,2,3],featureIndices: [0,0,0],counts: [9_007_199_254_740_993,1,1]),
            .init(cellCount: 3,featureCount: 4,rowOffsets: [0,2,3,4],featureIndices: [0,1,0,0],counts: [UInt64.max,1,1,1])
        ]
        for matrix in matrices {
            let training=VivoCompositionTraining(selection: original.selection,source: original.source,sourceReport: original.sourceReport,evidence: .synthetic,featureIDs: original.featureIDs,matrix: matrix)
            #expect(throws: (any Error).self) { try VivoComposition.fit(training) }
        }
    }
    @Test func bundlesReconstructAndRejectRehashedModelTampering() throws {
        let directory=FileManager.default.temporaryDirectory.appendingPathComponent("composition-test-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: directory,withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: directory) }
        let source=directory.appendingPathComponent("training.json"),model=directory.appendingPathComponent("model"),output=directory.appendingPathComponent("prediction")
        try VivoCanonicalJSON.encode(Self.fixture()).write(to: source)
        let hash=try VivoFingerprint(bytes: Array(repeating: 0,count: 32))
        let receipt=try VivoComposition.fit(source: source,implementation: hash,to: model)
        _=try VivoComposition.verifyModel(model,implementation: hash)
        _=try VivoComposition.predict(reference: model,plan: .init(context: Self.context,queries: [.init(id: "q",targets: ["a","b"])]),implementation: hash,to: output)
        _=try VivoComposition.verifyPrediction(output,implementation: hash)
        #expect(throws: (any Error).self) { try VivoComposition.fit(source: source,implementation: hash,to: model) }
        let path=model.appendingPathComponent("model.json")
        var object=try #require(JSONSerialization.jsonObject(with: Data(contentsOf: path)) as? [String:Any])
        var cpm=try #require(object["controlCPM"] as? [Double]);cpm[0]+=1;object["controlCPM"]=cpm
        let changed=try JSONSerialization.data(withJSONObject: object,options: .sortedKeys);try changed.write(to: path)
        let updated=try VivoCompositionReceipt(schemaVersion: 1,training: receipt.training,model: VivoCanonicalJSON.fingerprint(changed),query: nil,result: nil,implementation: hash)
        try VivoCanonicalJSON.encode(updated).write(to: model.appendingPathComponent("receipt.json"))
        #expect(throws: (any Error).self) { try VivoComposition.verifyModel(model,implementation: hash) }
    }
}
