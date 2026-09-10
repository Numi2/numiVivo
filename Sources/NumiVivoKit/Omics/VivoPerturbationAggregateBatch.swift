import Foundation

public struct VivoPerturbationAggregateFold: Codable, Sendable, Equatable {
    public let id: String
    public let perturbationID: String
    public let controlCondition: String
    public let treatmentCondition: String
    public let cellGroup: String
    public let trainingGroupIndices: [Int]
    public let queryGroupIndices: [Int]
    private enum CodingKeys: String, CodingKey {
        case id, perturbationID, controlCondition, treatmentCondition, cellGroup, trainingGroupIndices, queryGroupIndices
    }
    public init(id: String, perturbationID: String, controlCondition: String, treatmentCondition: String,
                cellGroup: String, trainingGroupIndices: [Int], queryGroupIndices: [Int]) {
        self.id=id;self.perturbationID=perturbationID;self.controlCondition=controlCondition
        self.treatmentCondition=treatmentCondition;self.cellGroup=cellGroup
        self.trainingGroupIndices=trainingGroupIndices;self.queryGroupIndices=queryGroupIndices
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["id","perturbationID","controlCondition","treatmentCondition","cellGroup","trainingGroupIndices","queryGroupIndices"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        id=try c.decode(String.self,forKey: .id);perturbationID=try c.decode(String.self,forKey: .perturbationID)
        controlCondition=try c.decode(String.self,forKey: .controlCondition);treatmentCondition=try c.decode(String.self,forKey: .treatmentCondition)
        cellGroup=try c.decode(String.self,forKey: .cellGroup)
        trainingGroupIndices=try c.decode([Int].self,forKey: .trainingGroupIndices);queryGroupIndices=try c.decode([Int].self,forKey: .queryGroupIndices)
    }
    func validate() throws {
        guard [id,perturbationID,controlCondition,treatmentCondition,cellGroup].allSatisfy(vivoOmicsID),
              controlCondition != treatmentCondition,
              (4...128).contains(trainingGroupIndices.count),(1...128).contains(queryGroupIndices.count),
              [trainingGroupIndices,queryGroupIndices].allSatisfy({ $0.allSatisfy { $0>=0 } && Set($0).count==$0.count }),
              Set(trainingGroupIndices).isDisjoint(with: queryGroupIndices) else {
            throw VivoOmicsError.invalid("aggregate prediction fold identities, conditions or group selection")
        }
    }
}

public struct VivoPerturbationAggregateBatchPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let sourceReport: VivoFingerprint
    public let featureNamespace: String
    public let provenance: String
    public let folds: [VivoPerturbationAggregateFold]
    private enum CodingKeys: String, CodingKey { case schemaVersion, sourceReport, featureNamespace, provenance, folds }
    public init(sourceReport: VivoFingerprint,featureNamespace: String,provenance: String,folds: [VivoPerturbationAggregateFold]) {
        schemaVersion=1;self.sourceReport=sourceReport;self.featureNamespace=featureNamespace;self.provenance=provenance;self.folds=folds
    }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder,allowed: ["schemaVersion","sourceReport","featureNamespace","provenance","folds"])
        let c=try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion=try c.decode(Int.self,forKey: .schemaVersion);sourceReport=try c.decode(VivoFingerprint.self,forKey: .sourceReport)
        featureNamespace=try c.decode(String.self,forKey: .featureNamespace);provenance=try c.decode(String.self,forKey: .provenance)
        folds=try c.decode([VivoPerturbationAggregateFold].self,forKey: .folds)
    }
    func validate() throws {
        guard schemaVersion==1,vivoOmicsID(featureNamespace),!provenance.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              provenance.utf8.count<=16_384,(1...128).contains(folds.count),Set(folds.map(\.id)).count==folds.count else {
            throw VivoOmicsError.invalid("aggregate prediction batch schema, provenance or fold identities")
        }
        for fold in folds { try fold.validate() }
    }
}

public struct VivoPerturbationAggregateFoldReceipt: Codable, Sendable, Equatable {
    public let id: String
    public let status: String
    public let trainingAggregate: VivoFingerprint?
    public let queryAggregate: VivoFingerprint?
    public let model: VivoFingerprint?
    public let prediction: VivoFingerprint?
    public let failure: String?
}
public struct VivoPerturbationAggregateBatchReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let source: VivoFingerprint
    public let sourceReport: VivoFingerprint
    public let sourceReceipt: VivoFingerprint
    public let plan: VivoFingerprint
    public let implementation: VivoFingerprint
    public let folds: [VivoPerturbationAggregateFoldReceipt]
}

/// Reconstruct the source once, then isolate each fold before entering the model.
/// Dense arrays contain donor aggregates only; source cell membership is retained.
public enum VivoPerturbationAggregateBatch {
    static let method="verified-pseudobulk-donor-response-batch-v1"
    static func project(_ source: VivoPseudobulkCounts,indices: [Int],cellGroup: String) throws -> VivoPseudobulkCounts {
        guard !indices.isEmpty,indices.count<=128,Set(indices).count==indices.count,
              indices.allSatisfy({ source.groups.indices.contains($0) }) else {
            throw VivoOmicsError.invalid("aggregate prediction source group selection")
        }
        var groups: [VivoPseudobulkGroup]=[],offsets=[0],columns: [Int]=[],counts: [UInt64]=[]
        for row in indices {
            let g=source.groups[row]
            groups.append(.init(biologicalReplicateID: g.biologicalReplicateID,donorID: g.donorID,condition: g.condition,
                organism: g.organism,cellGroup: cellGroup,sampleIDs: g.sampleIDs,batchIDs: g.batchIDs,sourceCellIndices: g.sourceCellIndices))
            let range=source.matrix.rowOffsets[row]..<source.matrix.rowOffsets[row+1]
            columns.append(contentsOf: source.matrix.featureIndices[range]);counts.append(contentsOf: source.matrix.counts[range]);offsets.append(counts.count)
        }
        return .init(method: "explicit-source-group-projection-v1",countUnit: source.countUnit,groups: groups,featureIDs: source.featureIDs,
            matrix: .init(cellCount: groups.count,featureCount: source.featureIDs.count,rowOffsets: offsets,featureIndices: columns,counts: counts))
    }
    static func isolatedInputs(_ source: VivoPseudobulkCounts,fold: VivoPerturbationAggregateFold) throws -> (VivoPseudobulkCounts,VivoPseudobulkCounts) {
        try fold.validate()
        let training=try project(source,indices: fold.trainingGroupIndices,cellGroup: fold.cellGroup)
        let query=try project(source,indices: fold.queryGroupIndices,cellGroup: fold.cellGroup)
        guard training.groups.allSatisfy({ $0.donorID != nil && ($0.condition==fold.controlCondition || $0.condition==fold.treatmentCondition) }),
              query.groups.allSatisfy({ $0.donorID != nil && $0.condition==fold.controlCondition }),
              Set(training.groups.compactMap(\.donorID)).isDisjoint(with: query.groups.compactMap(\.donorID)),
              Set(training.groups.flatMap(\.sampleIDs)).isDisjoint(with: query.groups.flatMap(\.sampleIDs)) else {
            throw VivoOmicsError.invalid("aggregate prediction requires disjoint training/query donors and control-only queries")
        }
        // Reject shared source cells even if the caller's sample/donor labels differ.
        let trainCells=Set(training.groups.flatMap(\.sourceCellIndices))
        guard query.groups.allSatisfy({ trainCells.isDisjoint(with: $0.sourceCellIndices) }) else {
            throw VivoOmicsError.invalid("aggregate prediction training/query source cells overlap")
        }
        return (training,query)
    }
    private static func mapping(_ original: VivoH5ADImportPlan,bulk: VivoPseudobulkCounts,id: String,provenance: String) throws -> VivoH5ADImportPlan {
        let ids=Set(bulk.groups.flatMap(\.sampleIDs)),samples=original.samples.filter { ids.contains($0.id) }
        guard samples.count==ids.count else { throw VivoOmicsError.invalid("aggregate prediction sample dictionary differs") }
        return .init(id: id,evidence: original.evidence,
            sourceDescription: "Aggregate projection; original source mapping retained for provenance. Cell group is explicitly remapped by the batch fold. "+provenance,
            countUnit: original.countUnit,matrixPath: original.matrixPath,samples: samples,sampleColumn: original.sampleColumn,
            barcodeColumn: original.barcodeColumn,groupColumn: original.groupColumn,featureIDColumn: original.featureIDColumn,
            featureNameColumn: original.featureNameColumn,mitochondrialFeatureIDs: original.mitochondrialFeatureIDs)
    }
    static func evaluateFold(_ source: VivoPseudobulkCounts,mapping original: VivoH5ADImportPlan,
                             batch: VivoPerturbationAggregateBatchPlan,fold: VivoPerturbationAggregateFold) throws -> (VivoPerturbationAggregateFoldReceipt,Data?,Data?) {
        var trainingID: VivoFingerprint?,queryID: VivoFingerprint?
        do {
            let (training,query)=try isolatedInputs(source,fold: fold)
            trainingID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(training))
            queryID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(query))
            let fit=try VivoPerturbationPlan(mapping: mapping(original,bulk: training,id: fold.id+"-training",provenance: batch.provenance),
                featureNamespace: batch.featureNamespace,perturbationID: fold.perturbationID,controlCondition: fold.controlCondition,
                treatmentCondition: fold.treatmentCondition,provenance: batch.provenance)
            let queryPlan=try VivoPerturbationQueryPlan(mapping: mapping(original,bulk: query,id: fold.id+"-query",provenance: batch.provenance),
                featureNamespace: batch.featureNamespace,perturbationID: fold.perturbationID)
            // Neither function receives the complete source, a target URL, or held-out treated rows.
            let model=try VivoPerturbation.model(training,plan: fit,source: trainingID!)
            let prediction=try VivoPerturbation.evaluate(query,plan: queryPlan,model: model)
            let modelBytes=try VivoCanonicalJSON.encode(model),predictionBytes=try VivoCanonicalJSON.encode(prediction)
            guard modelBytes.count<=67_108_864,predictionBytes.count<=67_108_864 else { throw VivoOmicsError.limit("aggregate prediction fold document") }
            return (.init(id: fold.id,status: "completed",trainingAggregate: trainingID,queryAggregate: queryID,
                model: try VivoCanonicalJSON.fingerprint(modelBytes),prediction: try VivoCanonicalJSON.fingerprint(predictionBytes),failure: nil),modelBytes,predictionBytes)
        } catch is CancellationError {
            throw CancellationError()
        } catch {
            return (.init(id: fold.id,status: "failed",trainingAggregate: trainingID,queryAggregate: queryID,
                          model: nil,prediction: nil,failure: String(describing: error)),nil,nil)
        }
    }
    private static func read(_ directory: URL,_ name: String,maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name),maximumBytes: maximum)
    }
    private static func staging(_ parent: URL) throws -> URL {
        let url=parent.appendingPathComponent(".numivivo-aggregate-prediction-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: url,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700]);return url
    }
    private static func copySource(_ source: URL,to destination: URL) throws {
        try FileManager.default.createDirectory(at: destination,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        for (name,limit) in [("original.h5ad",VivoH5ADPseudobulk.sourceLimits.maximumInputBytes),("plan.json",2_097_152),
                             ("report.json",536_870_912),("receipt.json",65_536)] {
            _=try VivoOmicsFileSnapshot.fingerprint(source.appendingPathComponent(name),copyTo: destination.appendingPathComponent(name),maximumBytes: limit)
        }
    }
    private static func source(_ directory: URL,plan: VivoPerturbationAggregateBatchPlan) throws -> (VivoH5ADPseudobulkReport,VivoH5ADPseudobulkPlan,VivoH5ADPseudobulkReceipt,VivoFingerprint) {
        let bytes=try read(directory,"receipt.json",maximum: 65_536)
        let receipt=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkReceipt.self,from: bytes)
        guard receipt.report==plan.sourceReport else { throw VivoOmicsError.invalid("aggregate prediction source report differs from frozen selection") }
        // Keep the publisher's receipt unchanged. Current code requalifies every
        // raw count and reconstructs the entire report, even across binaries.
        let report=try VivoH5ADPseudobulk.verify(directory,implementation: receipt.implementation)
        let mapping=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkPlan.self,from: read(directory,"plan.json",maximum: 2_097_152))
        return (report,mapping,receipt,try VivoCanonicalJSON.fingerprint(bytes))
    }
    public static func publish(sourceBundle: URL,plan: VivoPerturbationAggregateBatchPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoPerturbationAggregateBatchReceipt {
        try plan.validate()
        let planBytes=try VivoCanonicalJSON.encode(plan)
        guard planBytes.count<=2_097_152 else { throw VivoOmicsError.limit("aggregate prediction plan") }
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("aggregate prediction output already exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let sourceDirectory=temporary.appendingPathComponent("source");try copySource(sourceBundle,to: sourceDirectory)
        let (report,sourcePlan,sourceReceipt,sourceReceiptID)=try source(sourceDirectory,plan: plan)
        let output=temporary.appendingPathComponent("folds");try FileManager.default.createDirectory(at: output,withIntermediateDirectories: false)
        var receipts: [VivoPerturbationAggregateFoldReceipt]=[]
        for (index,fold) in plan.folds.enumerated() {
            try Task.checkCancellation()
            let (receipt,model,prediction)=try evaluateFold(report.pseudobulk,mapping: sourcePlan.mapping,batch: plan,fold: fold)
            let folder=output.appendingPathComponent(String(format: "%03d",index));try FileManager.default.createDirectory(at: folder,withIntermediateDirectories: false)
            try VivoCanonicalJSON.encode(receipt).write(to: folder.appendingPathComponent("receipt.json"),options: .withoutOverwriting)
            if let model,let prediction {
                try model.write(to: folder.appendingPathComponent("model.json"),options: .withoutOverwriting)
                try prediction.write(to: folder.appendingPathComponent("prediction.json"),options: .withoutOverwriting)
            }
            receipts.append(receipt)
        }
        let result=try VivoPerturbationAggregateBatchReceipt(schemaVersion: 1,method: method,source: sourceReceipt.source,sourceReport: sourceReceipt.report,
            sourceReceipt: sourceReceiptID,plan: VivoCanonicalJSON.fingerprint(planBytes),implementation: implementation,folds: receipts)
        try planBytes.write(to: temporary.appendingPathComponent("plan.json"),options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(result).write(to: temporary.appendingPathComponent("receipt.json"),options: .withoutOverwriting)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return result
    }
    public static func verify(_ directory: URL,implementation: VivoFingerprint) throws -> VivoPerturbationAggregateBatchReceipt {
        let receipt=try VivoCanonicalJSON.decode(VivoPerturbationAggregateBatchReceipt.self,from: read(directory,"receipt.json",maximum: 1_048_576))
        let planBytes=try read(directory,"plan.json",maximum: 2_097_152)
        guard receipt.schemaVersion==1,receipt.method==method,receipt.implementation==implementation,
              try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan else { throw VivoOmicsError.invalid("aggregate prediction receipt or plan changed") }
        let plan=try VivoCanonicalJSON.decode(VivoPerturbationAggregateBatchPlan.self,from: planBytes);try plan.validate()
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let copied=temporary.appendingPathComponent("source");try copySource(directory.appendingPathComponent("source"),to: copied)
        let (report,sourcePlan,sourceReceipt,sourceReceiptID)=try source(copied,plan: plan)
        guard receipt.source==sourceReceipt.source,receipt.sourceReport==sourceReceipt.report,receipt.sourceReceipt==sourceReceiptID,
              receipt.folds.count==plan.folds.count else { throw VivoOmicsError.invalid("aggregate prediction source or fold inventory changed") }
        for (index,fold) in plan.folds.enumerated() {
            try Task.checkCancellation()
            let (rebuilt,model,prediction)=try evaluateFold(report.pseudobulk,mapping: sourcePlan.mapping,batch: plan,fold: fold)
            let folder=directory.appendingPathComponent("folds").appendingPathComponent(String(format: "%03d",index))
            guard rebuilt==receipt.folds[index],try VivoCanonicalJSON.encode(rebuilt)==read(folder,"receipt.json",maximum: 65_536) else {
                throw VivoOmicsError.invalid("aggregate prediction fold receipt does not reconstruct")
            }
            if let model,let prediction {
                guard try model==read(folder,"model.json",maximum: 67_108_864),
                      try prediction==read(folder,"prediction.json",maximum: 67_108_864) else {
                    throw VivoOmicsError.invalid("aggregate prediction model or query does not reconstruct")
                }
            } else {
                guard !FileManager.default.fileExists(atPath: folder.appendingPathComponent("model.json").path),
                      !FileManager.default.fileExists(atPath: folder.appendingPathComponent("prediction.json").path) else {
                    throw VivoOmicsError.invalid("failed aggregate prediction fold contains fabricated output")
                }
            }
        }
        return receipt
    }
}
