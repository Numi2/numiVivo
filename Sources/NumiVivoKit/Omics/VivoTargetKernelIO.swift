import Foundation

public struct VivoTargetKernelReceipt: Codable,Sendable,Equatable {
    public let schemaVersion: Int
    public let training: VivoFingerprint
    public let plan: VivoFingerprint
    public let model: VivoFingerprint
    public let query: VivoFingerprint?
    public let result: VivoFingerprint?
    public let implementation: VivoFingerprint
}
extension VivoTargetKernel {
    private static func read(_ root: URL,_ name: String,maximum: Int = 67_108_864) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(root.appendingPathComponent(name),maximumBytes: maximum)
    }
    private static func write<T: Encodable>(_ value: T,_ root: URL,_ name: String,maximum: Int) throws -> Data {
        let data=try VivoCanonicalJSON.encode(value)
        guard data.count<=maximum else { throw VivoOmicsError.limit("target-kernel document: "+name) }
        try data.write(to: root.appendingPathComponent(name),options: .withoutOverwriting);return data
    }
    private static func staging(_ parent: URL) throws -> URL {
        let root=parent.appendingPathComponent(".numivivo-target-kernel-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: root,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700]);return root
    }
    private static func copyModel(_ source: URL,_ destination: URL) throws {
        try FileManager.default.createDirectory(at: destination,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        for (name,limit) in [("training.json",67_108_864),("plan.json",2_097_152),("model.json",134_217_728),("receipt.json",65_536)] {
            try read(source,name,maximum: limit).write(to: destination.appendingPathComponent(name),options: .withoutOverwriting)
        }
    }
    public static func fit(source: URL,plan: VivoTargetKernelPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoTargetKernelReceipt {
        try VivoH5ADCountStore.requireNew(destination)
        let bytes=try VivoSingleCellCampaignIO.readDocument(source,maximumBytes: 67_108_864)
        let training=try VivoCanonicalJSON.decode(VivoCompositionTraining.self,from: bytes)
        let model=try fit(training,plan: plan)
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        try bytes.write(to: temporary.appendingPathComponent("training.json"),options: .withoutOverwriting)
        let planBytes=try write(plan,temporary,"plan.json",maximum: 2_097_152)
        let modelBytes=try write(model,temporary,"model.json",maximum: 134_217_728)
        let receipt=try VivoTargetKernelReceipt(schemaVersion: 1,training: VivoCanonicalJSON.fingerprint(bytes),plan: VivoCanonicalJSON.fingerprint(planBytes),
            model: VivoCanonicalJSON.fingerprint(modelBytes),query: nil,result: nil,implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    private static func reconstruct(_ root: URL,implementation: VivoFingerprint) throws -> VivoTargetKernelModel {
        let receiptBytes=try read(root,"receipt.json",maximum: 65_536)
        let receipt=try VivoCanonicalJSON.decode(VivoTargetKernelReceipt.self,from: receiptBytes)
        let training=try read(root,"training.json"),plan=try read(root,"plan.json",maximum: 2_097_152),model=try read(root,"model.json",maximum: 134_217_728)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,receipt.query==nil,receipt.result==nil,
              try VivoCanonicalJSON.encode(receipt)==receiptBytes,
              try VivoCanonicalJSON.fingerprint(training)==receipt.training,try VivoCanonicalJSON.fingerprint(plan)==receipt.plan,
              try VivoCanonicalJSON.fingerprint(model)==receipt.model else { throw VivoOmicsError.invalid("target-kernel model fingerprints") }
        let result=try fit(VivoCanonicalJSON.decode(VivoCompositionTraining.self,from: training),plan: VivoCanonicalJSON.decode(VivoTargetKernelPlan.self,from: plan))
        guard try VivoCanonicalJSON.encode(result)==model else { throw VivoOmicsError.invalid("target-kernel model reconstruction") }
        return result
    }
    public static func verifyModel(_ directory: URL,implementation: VivoFingerprint) throws -> VivoTargetKernelModel {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let reference=temporary.appendingPathComponent("reference");try copyModel(directory,reference)
        return try reconstruct(reference,implementation: implementation)
    }
    public static func predict(reference: URL,plan: VivoTargetKernelQueryPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoTargetKernelReceipt {
        try VivoH5ADCountStore.requireNew(destination)
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyModel(reference,copy)
        let model=try reconstruct(copy,implementation: implementation),result=try predict(model,plan: plan)
        let queryBytes=try write(plan,temporary,"query.json",maximum: 2_097_152),resultBytes=try write(result,temporary,"report.json",maximum: 536_870_912)
        let parent=try VivoCanonicalJSON.decode(VivoTargetKernelReceipt.self,from: read(copy,"receipt.json",maximum: 65_536))
        let receipt=try VivoTargetKernelReceipt(schemaVersion: 1,training: parent.training,plan: parent.plan,model: parent.model,
            query: VivoCanonicalJSON.fingerprint(queryBytes),result: VivoCanonicalJSON.fingerprint(resultBytes),implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    public static func verifyPrediction(_ directory: URL,implementation: VivoFingerprint) throws -> VivoTargetKernelReport {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyModel(directory.appendingPathComponent("reference"),copy)
        let receiptBytes=try read(directory,"receipt.json",maximum: 65_536),query=try read(directory,"query.json",maximum: 2_097_152),report=try read(directory,"report.json",maximum: 536_870_912)
        let receipt=try VivoCanonicalJSON.decode(VivoTargetKernelReceipt.self,from: receiptBytes)
        let parent=try VivoCanonicalJSON.decode(VivoTargetKernelReceipt.self,from: read(copy,"receipt.json",maximum: 65_536))
        guard receipt.schemaVersion==1,receipt.implementation==implementation,try VivoCanonicalJSON.encode(receipt)==receiptBytes,
              receipt.training==parent.training,receipt.plan==parent.plan,receipt.model==parent.model,
              try VivoCanonicalJSON.fingerprint(query)==receipt.query,try VivoCanonicalJSON.fingerprint(report)==receipt.result else { throw VivoOmicsError.invalid("target-kernel prediction fingerprints") }
        let model=try reconstruct(copy,implementation: implementation)
        let result=try predict(model,plan: VivoCanonicalJSON.decode(VivoTargetKernelQueryPlan.self,from: query))
        guard try VivoCanonicalJSON.encode(result)==report else { throw VivoOmicsError.invalid("target-kernel prediction reconstruction") }
        return result
    }
}
