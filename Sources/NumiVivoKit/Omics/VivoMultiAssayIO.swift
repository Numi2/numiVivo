import Foundation

public struct VivoMultiAssayReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let dataset: VivoFingerprint
    public let h5mu: VivoFingerprint
    public let implementation: VivoFingerprint
    public let sourceFormat: String?
}

public enum VivoMultiAssayIO {
    private static func read(_ directory: URL, _ name: String, maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name), maximumBytes: maximum)
    }
    private static func staging(_ parent: URL) throws -> URL {
        let path = parent.appendingPathComponent(".numivivo-multiassay-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: path, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        return path
    }
    public static func importTenX(source: URL, plan: VivoTenXMultiAssayPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoMultiAssayReceipt {
        try VivoMultiAssayTenX.validate(plan)
        let bytes = try VivoCanonicalJSON.encode(plan)
        guard bytes.count <= 131_072 else { throw VivoOmicsError.limit("10x plan bytes") }
        return try publish(source: source, planBytes: bytes, sourceFormat: "10x", implementation: implementation, to: destination) {
            try VivoMultiAssayTenX.readSnapshot($0, plan: plan)
        }
    }
    public static func importH5MU(source: URL, plan: VivoH5MUMultiAssayPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoMultiAssayReceipt {
        try VivoMultiAssayH5MUImport.validate(plan)
        return try publish(source: source, planBytes: VivoCanonicalJSON.encode(plan), sourceFormat: "h5mu", implementation: implementation, to: destination) {
            try VivoMultiAssayH5MUImport.readSnapshot($0, plan: plan)
        }
    }
    private static func publish(source: URL, planBytes: Data, sourceFormat: String, implementation: VivoFingerprint, to destination: URL,
                                readDataset: (URL) throws -> VivoMultiAssayDataset) throws -> VivoMultiAssayReceipt {
        guard destination.isFileURL, !FileManager.default.fileExists(atPath: destination.path),
              (try? FileManager.default.attributesOfItem(atPath: destination.path)) == nil else { throw VivoOmicsError.invalid("multi-assay destination exists or is not local") }
        let temporary = try staging(destination.deletingLastPathComponent())
        defer { try? FileManager.default.removeItem(at: temporary) }
        let sourceHash = try VivoH5ADPseudobulk.fingerprint(source, copyTo: temporary.appendingPathComponent("original.h5"))
        guard planBytes.count <= 2_097_152 else { throw VivoOmicsError.limit("multi-assay plan bytes") }
        try planBytes.write(to: temporary.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        let dataset = try readDataset(temporary.appendingPathComponent("original.h5"))
        let bytes = try VivoCanonicalJSON.encode(dataset)
        guard bytes.count <= 536_870_912 else { throw VivoOmicsError.limit("multi-assay encoded dataset") }
        try bytes.write(to: temporary.appendingPathComponent("dataset.json"), options: .withoutOverwriting)
        let export = temporary.appendingPathComponent("dataset.h5mu")
        try VivoMultiAssayH5MU.writeSnapshot(dataset, to: export)
        let receipt = try VivoMultiAssayReceipt(schemaVersion: 1, source: sourceHash,
            plan: VivoCanonicalJSON.fingerprint(planBytes), dataset: VivoCanonicalJSON.fingerprint(bytes),
            h5mu: VivoH5ADPseudobulk.fingerprint(export), implementation: implementation, sourceFormat: sourceFormat)
        try VivoCanonicalJSON.encode(receipt).write(to: temporary.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> VivoMultiAssayReceipt {
        let receipt = try VivoCanonicalJSON.decode(VivoMultiAssayReceipt.self, from: read(directory, "receipt.json", maximum: 65_536))
        guard receipt.schemaVersion == 1, receipt.implementation == implementation else { throw VivoOmicsError.invalid("multi-assay receipt implementation") }
        let temp = try staging(FileManager.default.temporaryDirectory); defer { try? FileManager.default.removeItem(at: temp) }
        let source = try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("original.h5"), copyTo: temp.appendingPathComponent("original.h5"))
        let planBytes = try read(directory, "plan.json", maximum: 2_097_152)
        let datasetBytes = try read(directory, "dataset.json", maximum: 536_870_912)
        guard source == receipt.source, try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan,
              try VivoCanonicalJSON.fingerprint(datasetBytes) == receipt.dataset,
              try VivoH5ADPseudobulk.fingerprint(directory.appendingPathComponent("dataset.h5mu")) == receipt.h5mu else { throw VivoOmicsError.invalid("multi-assay artifact fingerprint") }
        let rebuilt: VivoMultiAssayDataset
        switch receipt.sourceFormat {
        case nil, "10x":
            let plan = try VivoCanonicalJSON.decode(VivoTenXMultiAssayPlan.self, from: planBytes)
            rebuilt = try VivoMultiAssayTenX.readSnapshot(temp.appendingPathComponent("original.h5"), plan: plan)
        case "h5mu":
            let plan = try VivoCanonicalJSON.decode(VivoH5MUMultiAssayPlan.self, from: planBytes)
            rebuilt = try VivoMultiAssayH5MUImport.readSnapshot(temp.appendingPathComponent("original.h5"), plan: plan)
        default: throw VivoOmicsError.invalid("unknown multi-assay source format")
        }
        guard try VivoCanonicalJSON.encode(rebuilt) == datasetBytes else { throw VivoOmicsError.invalid("multi-assay source reconstruction differs") }
        try VivoMultiAssayH5MU.writeSnapshot(rebuilt, to: temp.appendingPathComponent("reconstructed.h5mu"))
        guard try VivoH5ADPseudobulk.fingerprint(temp.appendingPathComponent("reconstructed.h5mu")) == receipt.h5mu else { throw VivoOmicsError.invalid("multi-assay H5MU reconstruction differs") }
        return receipt
    }
}
