import Foundation

public struct VivoCompositionReceipt: Codable,Sendable,Equatable {
    public let schemaVersion: Int
    public let training: VivoFingerprint
    public let model: VivoFingerprint
    public let query: VivoFingerprint?
    public let result: VivoFingerprint?
    public let implementation: VivoFingerprint
}

extension VivoComposition {
    private static func read(_ directory: URL,_ name: String,maximum: Int = 67_108_864) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name),maximumBytes: maximum)
    }
    private static func write<T: Encodable>(_ value: T,_ directory: URL,_ name: String,maximum: Int = 67_108_864) throws -> Data {
        let data=try VivoCanonicalJSON.encode(value)
        guard data.count<=maximum else { throw VivoOmicsError.limit("composition encoded document: "+name) }
        try data.write(to: directory.appendingPathComponent(name),options: .withoutOverwriting);return data
    }
    private static func staging(_ parent: URL) throws -> URL {
        let url=parent.appendingPathComponent(".numivivo-composition-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: url,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700]);return url
    }
    private static func copyModel(_ source: URL,_ destination: URL) throws {
        try FileManager.default.createDirectory(at: destination,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        for (name,limit) in [("training.json",67_108_864),("model.json",134_217_728),("receipt.json",65_536)] {
            try read(source,name,maximum: limit).write(to: destination.appendingPathComponent(name),options: .withoutOverwriting)
        }
    }
    /// Preparation may inspect the full source, but the resulting file contains
    /// only declared control/single conditions. It is the sole fitting input.
    public static func prepare(from source: URL,plan: VivoCompositionPreparation,implementation: VivoFingerprint,to destination: URL) throws {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("composition training output exists") }
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        for (name,limit) in [("plan.json",2_097_152),("report.json",536_870_912),("receipt.json",65_536)] {
            try read(source,name,maximum: limit).write(to: temporary.appendingPathComponent(name),options: .withoutOverwriting)
        }
        _=try VivoH5ADPseudobulk.fingerprint(source.appendingPathComponent("original.h5ad"),copyTo: temporary.appendingPathComponent("original.h5ad"))
        let report=try VivoH5ADPseudobulk.verify(temporary,implementation: implementation)
        let receipt=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkReceipt.self,from: read(temporary,"receipt.json",maximum: 65_536))
        let training=try select(report,plan: plan,source: receipt.source,sourceReport: receipt.report)
        let bytes=try VivoCanonicalJSON.encode(training)
        guard bytes.count<=67_108_864 else { throw VivoOmicsError.limit("composition training exceeds 64 MiB") }
        try VivoSingleCellH5AD.publish(bytes,to: destination)
    }
    public static func fit(source: URL,implementation: VivoFingerprint,to destination: URL) throws -> VivoCompositionReceipt {
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("composition model output exists") }
        let bytes=try VivoSingleCellCampaignIO.readDocument(source,maximumBytes: 67_108_864)
        let training=try VivoCanonicalJSON.decode(VivoCompositionTraining.self,from: bytes)
        let model=try fit(training)
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        try bytes.write(to: temporary.appendingPathComponent("training.json"),options: .withoutOverwriting)
        let modelBytes=try write(model,temporary,"model.json",maximum: 134_217_728)
        let receipt=try VivoCompositionReceipt(schemaVersion: 1,training: VivoCanonicalJSON.fingerprint(bytes),model: VivoCanonicalJSON.fingerprint(modelBytes),query: nil,result: nil,implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    private static func reconstructModel(_ directory: URL,implementation: VivoFingerprint) throws -> VivoCompositionModel {
        let receipt=try VivoCanonicalJSON.decode(VivoCompositionReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let trainingBytes=try read(directory,"training.json"),modelBytes=try read(directory,"model.json",maximum: 134_217_728)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,receipt.query==nil,receipt.result==nil,
              try VivoCanonicalJSON.fingerprint(trainingBytes)==receipt.training,try VivoCanonicalJSON.fingerprint(modelBytes)==receipt.model else {
            throw VivoOmicsError.invalid("composition model hashes or implementation mismatch")
        }
        let training=try VivoCanonicalJSON.decode(VivoCompositionTraining.self,from: trainingBytes)
        let model=try fit(training)
        guard try VivoCanonicalJSON.encode(model)==modelBytes else { throw VivoOmicsError.invalid("composition model does not reconstruct") }
        return model
    }
    public static func verifyModel(_ directory: URL,implementation: VivoFingerprint) throws -> VivoCompositionModel {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let reference=temporary.appendingPathComponent("reference");try copyModel(directory,reference)
        return try reconstructModel(reference,implementation: implementation)
    }
    public static func predict(reference: URL,plan: VivoCompositionQueryPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoCompositionReceipt {
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("composition prediction output exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyModel(reference,copy)
        let model=try reconstructModel(copy,implementation: implementation)
        let result=try predict(model,plan: plan)
        let queryBytes=try write(plan,temporary,"query.json",maximum: 2_097_152)
        let reportBytes=try write(result,temporary,"report.json",maximum: 536_870_912)
        let receipt=try VivoCompositionReceipt(schemaVersion: 1,training: VivoCanonicalJSON.fingerprint(read(copy,"training.json")),
            model: VivoCanonicalJSON.fingerprint(read(copy,"model.json",maximum: 134_217_728)),query: VivoCanonicalJSON.fingerprint(queryBytes),
            result: VivoCanonicalJSON.fingerprint(reportBytes),implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    public static func verifyPrediction(_ directory: URL,implementation: VivoFingerprint) throws -> VivoCompositionReport {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyModel(directory.appendingPathComponent("reference"),copy)
        let receipt=try VivoCanonicalJSON.decode(VivoCompositionReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let queryBytes=try read(directory,"query.json",maximum: 2_097_152),reportBytes=try read(directory,"report.json",maximum: 536_870_912)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,
              try VivoCanonicalJSON.fingerprint(read(copy,"training.json"))==receipt.training,
              try VivoCanonicalJSON.fingerprint(read(copy,"model.json",maximum: 134_217_728))==receipt.model,
              try VivoCanonicalJSON.fingerprint(queryBytes)==receipt.query,try VivoCanonicalJSON.fingerprint(reportBytes)==receipt.result else {
            throw VivoOmicsError.invalid("composition prediction hashes or implementation mismatch")
        }
        let model=try reconstructModel(copy,implementation: implementation)
        let plan=try VivoCanonicalJSON.decode(VivoCompositionQueryPlan.self,from: queryBytes)
        let result=try predict(model,plan: plan)
        guard try VivoCanonicalJSON.encode(result)==reportBytes else { throw VivoOmicsError.invalid("composition prediction does not reconstruct") }
        return result
    }
}
