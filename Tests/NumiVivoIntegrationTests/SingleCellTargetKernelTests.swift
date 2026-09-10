import Foundation
import Testing
@testable import NumiVivoKit

@Suite struct SingleCellTargetKernelTests {
    static let context=SingleCellCompositionTests.context
    static var hash: VivoFingerprint { get throws { try .init(bytes: Array(repeating: 0,count: 32)) } }
    static let descriptors: [VivoTargetDescriptor]=[
        .init(targetID: "a",featureID: "gene-a",status: .available,terms: ["x","y"]),
        .init(targetID: "b",featureID: "gene-b",status: .available,terms: ["y","z"])]
    static func plan(_ targets: [VivoTargetDescriptor] = descriptors,regularization: Double = 1,work: Int = 200_000_000) throws -> VivoTargetKernelPlan {
        try .init(context: context,descriptorNamespace: "synthetic-terms",source: hash,provenance: "Synthetic numerical control",targets: targets,regularization: regularization,maximumWork: work)
    }
    static func query(_ descriptor: VivoTargetDescriptor = .init(targetID: "held",featureID: "gene-held",status: .available,terms: ["x","y"]),work: Int = 200_000_000) throws -> VivoTargetKernelQueryPlan {
        try .init(context: context,descriptorNamespace: "synthetic-terms",source: hash,queries: [.init(id: "q",descriptor: descriptor)],maximumWork: work)
    }
    @Test func analyticalWeightsAndMatchedShuffle() throws {
        let model=try VivoTargetKernel.fit(SingleCellCompositionTests.fixture(),plan: Self.plan())
        let report=try VivoTargetKernel.predict(model,plan: Self.query()),q=try #require(report.queries.first)
        let weights=try #require(q.weights),expression=try #require(q.expression),shuffled=try #require(q.shuffledExpression)
        // With Jaccard=1/3, lambda=1 and query equal to a, intercept-constrained weights are 7/10 and 3/10.
        #expect(abs(weights[0]-0.7)<1e-14 && abs(weights[1]-0.3)<1e-14)
        for j in model.base.featureIDs.indices {
            let a=log1p(model.base.targetCPM[0][j]),b=log1p(model.base.targetCPM[1][j])
            #expect(abs(expression[j]-(0.7*a+0.3*b))<1e-12)
            #expect(abs(shuffled[j]-(0.3*a+0.7*b))<1e-12)
            #expect(abs(report.baselines[1].expression[j]-(a+b)/2)<1e-12)
        }
        #expect(q.status=="predicted" && q.maximumSimilarity==1 && q.diagnostic?.renormalized==false)
        #expect(model.maximumInverseResidual<1e-14)
        let repeatModel=try VivoTargetKernel.fit(SingleCellCompositionTests.fixture(pair: [999,1,99,2]),plan: Self.plan(Self.descriptors.reversed()))
        #expect(model==repeatModel)
    }
    @Test func missingAndDisjointDescriptorsRemainUnsupported() throws {
        let model=try VivoTargetKernel.fit(SingleCellCompositionTests.fixture(),plan: Self.plan())
        for (descriptor,status) in [
            (VivoTargetDescriptor(targetID: "held",featureID: "gene-held",status: .available,terms: ["unrelated"]),"noSharedTerms"),
            (.init(targetID: "held",featureID: "gene-held",status: .noData,terms: []),"noData"),
            (.init(targetID: "held",featureID: nil,status: .unresolvedIdentity,terms: []),"unresolvedIdentity")] {
            let report=try VivoTargetKernel.predict(model,plan: Self.query(descriptor)),q=try #require(report.queries.first)
            #expect(q.status==status && q.expression==nil && q.shuffledExpression==nil && q.weights==nil && q.diagnostic==nil)
            #expect(report.baselines.count==3)
        }
    }
    @Test func descriptorsIdentityAndBudgetsReject() throws {
        let training=try SingleCellCompositionTests.fixture(),model=try VivoTargetKernel.fit(training,plan: Self.plan())
        for bad in [Self.descriptors[0],.init(targetID: "new-alias",featureID: "gene-b",status: .available,terms: ["x"]),
                    .init(targetID: "held",featureID: "new",status: .available,terms: ["x","x"]),
                    .init(targetID: "held",featureID: nil,status: .available,terms: ["x"]),
                    .init(targetID: "held",featureID: "new",status: .noData,terms: ["x"])] {
            #expect(throws: (any Error).self) { try VivoTargetKernel.predict(model,plan: Self.query(bad)) }
        }
        for regularization in [-1.0,0,Double.nan,Double.infinity,1001] {
            #expect(throws: (any Error).self) { try VivoTargetKernel.fit(training,plan: Self.plan(regularization: regularization)) }
        }
        for targets in [[Self.descriptors[0]],Self.descriptors+[.init(targetID: "extra",featureID: "gene-extra",status: .available,terms: ["x"])],
                        [Self.descriptors[0],.init(targetID: "b",featureID: "gene-a",status: .available,terms: ["z"])],
                        [Self.descriptors[0],.init(targetID: "b",featureID: "gene-b",status: .noData,terms: [])]] {
            #expect(throws: (any Error).self) { try VivoTargetKernel.fit(training,plan: Self.plan(targets)) }
        }
        #expect(throws: (any Error).self) { try VivoTargetKernel.fit(training,plan: Self.plan(work: 1)) }
        #expect(throws: (any Error).self) { try VivoTargetKernel.predict(model,plan: Self.query(work: 1)) }
        let wrong=try VivoTargetKernelQueryPlan(context: Self.context,descriptorNamespace: "other",source: Self.hash,queries: Self.query().queries)
        #expect(throws: (any Error).self) { try VivoTargetKernel.predict(model,plan: wrong) }
        let invalid=VivoCompositionTraining(selection: training.selection,source: training.source,sourceReport: training.sourceReport,evidence: .synthetic,
            featureIDs: training.featureIDs,matrix: .init(cellCount: 3,featureCount: 4,rowOffsets: [0,0,1,2],featureIndices: [0,0],counts: [1,1]))
        #expect(throws: (any Error).self) { try VivoTargetKernel.fit(invalid,plan: Self.plan()) }
    }
    @Test func unknownFieldsRejectAtEveryInputBoundary() throws {
        func withUnknown<T: Encodable>(_ value: T) throws -> Data {
            var json=try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(value)) as? [String:Any])
            json["heldOutExpression"]=[1,2];return try JSONSerialization.data(withJSONObject: json,options: .sortedKeys)
        }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoTargetDescriptor.self,from: withUnknown(Self.descriptors[0])) }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoTargetKernelPlan.self,from: withUnknown(Self.plan())) }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoTargetKernelQueryPlan.self,from: withUnknown(Self.query())) }
        #expect(throws: (any Error).self) { try JSONDecoder().decode(VivoTargetKernelQuery.self,from: withUnknown(Self.query().queries[0])) }
    }
    @Test func bundlesReconstructAndRejectRehashedChanges() throws {
        let root=FileManager.default.temporaryDirectory.appendingPathComponent("target-kernel-test-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: root,withIntermediateDirectories: false)
        defer { try? FileManager.default.removeItem(at: root) }
        let source=root.appendingPathComponent("training.json"),model=root.appendingPathComponent("model"),output=root.appendingPathComponent("prediction")
        try VivoCanonicalJSON.encode(SingleCellCompositionTests.fixture()).write(to: source)
        let hash=try Self.hash
        _=try VivoTargetKernel.fit(source: source,plan: Self.plan(),implementation: hash,to: model)
        _=try VivoTargetKernel.verifyModel(model,implementation: hash)
        _=try VivoTargetKernel.predict(reference: model,plan: Self.query(),implementation: hash,to: output)
        _=try VivoTargetKernel.verifyPrediction(output,implementation: hash)
        #expect(throws: (any Error).self) { try VivoTargetKernel.fit(source: source,plan: Self.plan(),implementation: hash,to: model) }
        #expect(throws: (any Error).self) { try VivoTargetKernel.predict(reference: model,plan: Self.query(),implementation: hash,to: output) }
        let other=try VivoFingerprint(bytes: Array(repeating: 1,count: 32))
        #expect(throws: (any Error).self) { try VivoTargetKernel.verifyModel(model,implementation: other) }
        for name in ["model","plan","query","report"] {
            let isPrediction=name=="query" || name=="report",copy=root.appendingPathComponent("tampered-"+name)
            try FileManager.default.copyItem(at: isPrediction ? output : model,to: copy)
            let path=copy.appendingPathComponent(name+".json")
            var object=try #require(JSONSerialization.jsonObject(with: Data(contentsOf: path)) as? [String:Any])
            switch name {
            case "model": var inverse=try #require(object["inverse"] as? [Double]);inverse[0]+=0.01;object["inverse"]=inverse
            case "plan": object["regularization"]=2
            case "query": var queries=try #require(object["queries"] as? [[String:Any]]);queries[0]["id"]="altered";object["queries"]=queries
            default: object["qualification"]="altered"
            }
            let changed=try JSONSerialization.data(withJSONObject: object,options: .sortedKeys);try changed.write(to: path)
            let receiptPath=copy.appendingPathComponent("receipt.json"),r=try VivoCanonicalJSON.decode(VivoTargetKernelReceipt.self,from: Data(contentsOf: receiptPath)),fingerprint=try VivoCanonicalJSON.fingerprint(changed)
            let updated=VivoTargetKernelReceipt(schemaVersion: 1,training: r.training,plan: name=="plan" ? fingerprint : r.plan,model: name=="model" ? fingerprint : r.model,
                query: name=="query" ? fingerprint : r.query,result: name=="report" ? fingerprint : r.result,implementation: hash)
            try VivoCanonicalJSON.encode(updated).write(to: receiptPath)
            if isPrediction { #expect(throws: (any Error).self) { try VivoTargetKernel.verifyPrediction(copy,implementation: hash) } }
            else { #expect(throws: (any Error).self) { try VivoTargetKernel.verifyModel(copy,implementation: hash) } }
        }
        #expect(throws: (any Error).self) { try VivoTargetKernel.predict(reference: model,plan: Self.query(Self.descriptors[0]),implementation: hash,to: root.appendingPathComponent("rejected")) }
        #expect(!FileManager.default.fileExists(atPath: root.appendingPathComponent("rejected").path))
        #expect(try FileManager.default.contentsOfDirectory(atPath: root.path).allSatisfy { !$0.hasPrefix(".numivivo-target-kernel-") })
    }
}
