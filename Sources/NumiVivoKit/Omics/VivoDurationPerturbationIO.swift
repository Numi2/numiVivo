import Foundation

extension VivoDurationPerturbation {
    private static func read(_ directory: URL, _ name: String, maximum: Int = 67_108_864) throws -> Data {
        try VivoPerturbation.read(directory, name, maximum: maximum)
    }
    private static func write<T: Encodable>(_ value: T, _ directory: URL, _ name: String, maximum: Int = 67_108_864) throws -> Data {
        try VivoPerturbation.write(value, directory, name, maximum: maximum)
    }
    private static func staging(_ parent: URL) throws -> URL { try VivoPerturbation.staging(parent) }
    private static func copyReference(_ source: URL, to destination: URL) throws { try VivoPerturbation.copyReference(source, to: destination) }
    private static func evaluate(snapshot: URL, plan: VivoDurationQueryPlan, model: VivoDurationModel) throws -> VivoDurationReport {
        try plan.validate()
        let report = try VivoH5ADPseudobulk.evaluateSnapshot(snapshot, plan: .init(mapping: plan.mapping))
        return try evaluate(report.pseudobulk, plan: plan, model: model)
    }
    private static func reconstructReference(_ directory: URL,implementation: VivoFingerprint) throws -> VivoDurationModel {
        let receipt=try VivoCanonicalJSON.decode(VivoPerturbationReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let planBytes=try read(directory,"plan.json",maximum: 2_097_152),modelBytes=try read(directory,"model.json")
        guard receipt.schemaVersion==1,receipt.implementation==implementation,receipt.reference==nil,
              try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan,try VivoCanonicalJSON.fingerprint(modelBytes)==receipt.result else {
            throw VivoOmicsError.invalid("duration perturbation bundle hash, schema or implementation mismatch")
        }
        let plan=try VivoCanonicalJSON.decode(VivoDurationPlan.self,from: planBytes)
        try plan.validate()
        let trainingPlan=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkPlan.self,from: read(directory,"training/plan.json",maximum: 2_097_152))
        guard trainingPlan == .init(mapping: plan.mapping) else { throw VivoOmicsError.invalid("duration perturbation training plan mismatch") }
        let training=try VivoH5ADPseudobulk.verify(directory.appendingPathComponent("training"),implementation: implementation)
        let source=try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("training/original.h5ad"))
        guard source==receipt.source else { throw VivoOmicsError.invalid("duration perturbation training source mismatch") }
        let rebuilt=try model(training.pseudobulk,plan: plan,source: source)
        guard try VivoCanonicalJSON.encode(rebuilt)==modelBytes else { throw VivoOmicsError.invalid("duration perturbation model does not reconstruct") }
        return rebuilt
    }
    public static func fit(source: URL,plan: VivoDurationPlan,implementation: VivoFingerprint,to destination: URL) throws -> VivoPerturbationReceipt {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("reference output already exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let training=temporary.appendingPathComponent("training")
        let receipt=try VivoH5ADPseudobulk.publish(source: source,plan: .init(mapping: plan.mapping),implementation: implementation,to: training)
        let report=try VivoCanonicalJSON.decode(VivoH5ADPseudobulkReport.self,from: read(training,"report.json",maximum: 536_870_912))
        let model=try model(report.pseudobulk,plan: plan,source: receipt.source)
        let planBytes=try write(plan,temporary,"plan.json",maximum: 2_097_152),modelBytes=try write(model,temporary,"model.json")
        let result=try VivoPerturbationReceipt(schemaVersion: 1,source: receipt.source,plan: VivoCanonicalJSON.fingerprint(planBytes),
            result: VivoCanonicalJSON.fingerprint(modelBytes),reference: nil,implementation: implementation)
        _=try write(result,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return result
    }
    public static func verifyModel(_ directory: URL,implementation: VivoFingerprint) throws -> VivoDurationModel {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyReference(directory,to: copy)
        return try reconstructReference(copy,implementation: implementation)
    }
    public static func map(source: URL,plan: VivoDurationQueryPlan,reference: URL,implementation: VivoFingerprint,to destination: URL) throws -> VivoPerturbationReceipt {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("duration perturbation query output already exists") }
        let temporary=try staging(destination.deletingLastPathComponent());defer { try? FileManager.default.removeItem(at: temporary) }
        let copy=temporary.appendingPathComponent("reference");try copyReference(reference,to: copy)
        let model=try reconstructReference(copy,implementation: implementation)
        let snapshot=temporary.appendingPathComponent("original.h5ad"),sourceID=try VivoH5ADPseudobulk.fingerprint(source,copyTo: snapshot)
        let report=try evaluate(snapshot: snapshot,plan: plan,model: model)
        let planBytes=try write(plan,temporary,"plan.json",maximum: 2_097_152),reportBytes=try write(report,temporary,"report.json")
        let receipt=try VivoPerturbationReceipt(schemaVersion: 1,source: sourceID,plan: VivoCanonicalJSON.fingerprint(planBytes),
            result: VivoCanonicalJSON.fingerprint(reportBytes),reference: VivoCanonicalJSON.fingerprint(read(copy,"model.json")),implementation: implementation)
        _=try write(receipt,temporary,"receipt.json",maximum: 65_536)
        try Task.checkCancellation();try FileManager.default.moveItem(at: temporary,to: destination);return receipt
    }
    public static func verifyPrediction(_ directory: URL,implementation: VivoFingerprint) throws -> VivoDurationReport {
        let temporary=try staging(FileManager.default.temporaryDirectory);defer { try? FileManager.default.removeItem(at: temporary) }
        let receipt=try VivoCanonicalJSON.decode(VivoPerturbationReceipt.self,from: read(directory,"receipt.json",maximum: 65_536))
        let planBytes=try read(directory,"plan.json",maximum: 2_097_152),reportBytes=try read(directory,"report.json")
        let copy=temporary.appendingPathComponent("reference");try copyReference(directory.appendingPathComponent("reference"),to: copy)
        guard receipt.schemaVersion==1,receipt.implementation==implementation,
              try VivoCanonicalJSON.fingerprint(planBytes)==receipt.plan,try VivoCanonicalJSON.fingerprint(reportBytes)==receipt.result,
              try VivoCanonicalJSON.fingerprint(read(copy,"model.json"))==receipt.reference else { throw VivoOmicsError.invalid("duration perturbation query bundle hashes or implementation mismatch") }
        let model=try reconstructReference(copy,implementation: implementation)
        let snapshot=temporary.appendingPathComponent("original.h5ad")
        guard try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("original.h5ad"),copyTo: snapshot)==receipt.source else {
            throw VivoOmicsError.invalid("duration perturbation query source changed")
        }
        let plan=try VivoCanonicalJSON.decode(VivoDurationQueryPlan.self,from: planBytes)
        let result=try evaluate(snapshot: snapshot,plan: plan,model: model)
        guard try VivoCanonicalJSON.encode(result)==reportBytes else { throw VivoOmicsError.invalid("duration perturbation query result does not reconstruct") }
        return result
    }
}
