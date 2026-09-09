import Foundation

public struct VivoSingleCellReferenceReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let result: VivoFingerprint
    public let reference: VivoFingerprint?
    public let implementation: VivoFingerprint
}

extension VivoSingleCellReference {
    private static func read(_ directory: URL,_ name: String,maximum: Int = 67_108_864) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name),maximumBytes: maximum)
    }
    private static func write<T: Encodable>(_ value: T,_ directory: URL,_ name: String,maximum: Int = 67_108_864) throws -> Data {
        let bytes=try VivoCanonicalJSON.encode(value)
        guard bytes.count<=maximum else { throw VivoOmicsError.limit("reference document size: "+name) }
        try bytes.write(to: directory.appendingPathComponent(name),options: .withoutOverwriting);return bytes
    }
    private static func staging(_ parent: URL) throws -> URL {
        let url=parent.appendingPathComponent(".numivivo-reference-"+UUID().uuidString)
        try FileManager.default.createDirectory(at: url,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700]);return url
    }
    private static func copyReference(_ source: URL,to destination: URL) throws {
        try FileManager.default.createDirectory(at: destination,withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        try FileManager.default.createDirectory(at: destination.appendingPathComponent("training"),withIntermediateDirectories: false,attributes: [.posixPermissions: 0o700])
        for (name,limit) in [("plan.json",2_097_152),("model.json",67_108_864),("receipt.json",65_536),
                            ("training/plan.json",2_097_152),("training/report.json",536_870_912),("training/receipt.json",65_536)] {
            try read(source,name,maximum: limit).write(to: destination.appendingPathComponent(name),options: .withoutOverwriting)
        }
        _=try VivoH5ADPseudobulk.fingerprint(source.appendingPathComponent("training/original.h5ad"),copyTo: destination.appendingPathComponent("training/original.h5ad"))
    }
    private static func reconstructReference(_ directory: URL,implementation: VivoFingerprint) throws -> VivoSingleCellReferenceModel {
        let receipt=try VivoCanonicalJSON.decode(VivoSingleCellReferenceReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let planBytes=try read(directory,"plan.json",maximum: 2_097_152),modelBytes=try read(directory,"model.json")
        guard receipt.schemaVersion==1,receipt.implementation==implementation,receipt.reference==nil,
              try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan,try VivoCanonicalJSON.fingerprint(modelBytes)==receipt.result else {
            throw VivoOmicsError.invalid("reference bundle hash, schema or implementation mismatch")
        }
        let plan=try VivoCanonicalJSON.decode(VivoSingleCellReferencePlan.self,from: planBytes)
        try plan.validate()
        let trainingPlan=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkPlan.self,from: read(directory,"training/plan.json",maximum: 2_097_152))
        guard trainingPlan==plan.trainingPlan else { throw VivoOmicsError.invalid("reference training plan mismatch") }
        let training=try VivoH5ADPseudobulk.verify(directory.appendingPathComponent("training"),implementation: implementation)
        let source=try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("training/original.h5ad"))
        guard source==receipt.source else { throw VivoOmicsError.invalid("reference training source mismatch") }
        let rebuilt=try model(training,plan: plan,source: source)
        guard try VivoCanonicalJSON.encode(rebuilt)==modelBytes else { throw VivoOmicsError.invalid("reference model does not reconstruct") }
        return rebuilt
    }
    public static func fit(source: URL,plan: VivoSingleCellReferencePlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoSingleCellReferenceReceipt {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("reference output already exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let training=temporary.appendingPathComponent("training")
        let receipt=try VivoH5ADPseudobulk.publish(source: source,plan: plan.trainingPlan,implementation: implementation,to: training)
        let report=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkReport.self,from: read(training,"report.json",maximum: 536_870_912))
        let model=try model(report,plan: plan,source: receipt.source)
        let planBytes=try write(plan,temporary,"plan.json",maximum: 2_097_152),modelBytes=try write(model,temporary,"model.json")
        let result=try VivoSingleCellReferenceReceipt(schemaVersion: 1,source: receipt.source,plan: VivoCanonicalJSON.fingerprint(planBytes),
            result: VivoCanonicalJSON.fingerprint(modelBytes),reference: nil,implementation: implementation)
        _=try write(result,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return result
    }
    public static func verifyReference(_ directory: URL,implementation: VivoFingerprint) throws -> VivoSingleCellReferenceModel {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyReference(directory,to: copy)
        return try reconstructReference(copy,implementation: implementation)
    }
    public static func map(source: URL,plan: VivoSingleCellReferenceQueryPlan,reference: URL,implementation: VivoFingerprint,to destination: URL) throws -> VivoSingleCellReferenceReceipt {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("reference query output already exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyReference(reference,to: copy)
        let model=try reconstructReference(copy,implementation: implementation)
        let snapshot=temporary.appendingPathComponent("original.h5ad"),sourceID=try VivoH5ADPseudobulk.fingerprint(source,copyTo: snapshot)
        let report=try evaluate(snapshot: snapshot,plan: plan,model: model)
        let planBytes=try write(plan,temporary,"plan.json",maximum: 2_097_152),reportBytes=try write(report,temporary,"report.json")
        let receipt=try VivoSingleCellReferenceReceipt(schemaVersion: 1,source: sourceID,plan: VivoCanonicalJSON.fingerprint(planBytes),
            result: VivoCanonicalJSON.fingerprint(reportBytes),reference: VivoCanonicalJSON.fingerprint(read(copy,"model.json")),implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    public static func verifyMapping(_ directory: URL,implementation: VivoFingerprint) throws -> VivoSingleCellReferenceReport {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let receipt=try VivoCanonicalJSON.decode(VivoSingleCellReferenceReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let planBytes=try read(directory,"plan.json",maximum: 2_097_152),reportBytes=try read(directory,"report.json")
        let copy=temporary.appendingPathComponent("reference");try copyReference(directory.appendingPathComponent("reference"),to: copy)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,
              try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan,try VivoCanonicalJSON.fingerprint(reportBytes)==receipt.result,
              try VivoCanonicalJSON.fingerprint(read(copy,"model.json"))==receipt.reference else { throw VivoOmicsError.invalid("reference query bundle hashes or implementation mismatch") }
        let model=try reconstructReference(copy,implementation: implementation)
        let snapshot=temporary.appendingPathComponent("original.h5ad")
        guard try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("original.h5ad"),copyTo: snapshot)==receipt.source else {
            throw VivoOmicsError.invalid("reference query source changed")
        }
        let plan=try VivoCanonicalJSON.decode(VivoSingleCellReferenceQueryPlan.self,from: planBytes)
        let result=try evaluate(snapshot: snapshot,plan: plan,model: model)
        guard try VivoCanonicalJSON.encode(result)==reportBytes else { throw VivoOmicsError.invalid("reference query result does not reconstruct") }
        return result
    }
}
