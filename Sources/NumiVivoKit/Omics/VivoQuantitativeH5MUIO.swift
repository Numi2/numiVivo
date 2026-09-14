import Foundation

public struct VivoQuantitativeH5MUReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let source: VivoFingerprint
    public let plan: VivoFingerprint
    public let dataset: VivoFingerprint
    public let h5mu: VivoFingerprint
    public let implementation: VivoFingerprint

    public init(schemaVersion: Int = 1, source: VivoFingerprint, plan: VivoFingerprint,
                dataset: VivoFingerprint, h5mu: VivoFingerprint, implementation: VivoFingerprint) {
        self.schemaVersion = schemaVersion; self.source = source; self.plan = plan
        self.dataset = dataset; self.h5mu = h5mu; self.implementation = implementation
    }

    private enum CodingKeys: String, CodingKey { case schemaVersion, source, plan, dataset, h5mu, implementation }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "source", "plan", "dataset", "h5mu", "implementation"])
        let v = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try v.decode(Int.self, forKey: .schemaVersion)
        source = try v.decode(VivoFingerprint.self, forKey: .source)
        plan = try v.decode(VivoFingerprint.self, forKey: .plan)
        dataset = try v.decode(VivoFingerprint.self, forKey: .dataset)
        h5mu = try v.decode(VivoFingerprint.self, forKey: .h5mu)
        implementation = try v.decode(VivoFingerprint.self, forKey: .implementation)
    }
}

/// Transactional, source-bound packaging for quantitative H5MU.  The source
/// file remains available for replay while the normalized H5MU export and
/// canonical JSON dataset are fingerprinted independently.
public enum VivoQuantitativeH5MUIO {
    private static let maximumSourceBytes = VivoQuantitativeAssayDataset.limits.maximumInputBytes
    private static let maximumPlanBytes = 2_097_152
    private static let maximumDatasetBytes = 536_870_912
    private static let maximumReceiptBytes = 65_536

    private static func staging(_ parent: URL) throws -> URL {
        let path = parent.appendingPathComponent(".numivivo-quantitative-h5mu-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: path, withIntermediateDirectories: false,
                                                 attributes: [.posixPermissions: 0o700])
        return path
    }

    private static func read(_ root: URL, _ name: String, maximum: Int) throws -> Data {
        try VivoSingleCellCampaignIO.readDocument(root.appendingPathComponent(name), maximumBytes: maximum)
    }

    private static func fingerprint(_ source: URL, copyTo destination: URL? = nil) throws -> VivoFingerprint {
        try VivoOmicsFileSnapshot.fingerprint(source, copyTo: destination, maximumBytes: maximumSourceBytes)
    }

    private static func requireNew(_ destination: URL) throws {
        guard destination.isFileURL,
              (try? FileManager.default.attributesOfItem(atPath: destination.path)) == nil else {
            throw VivoOmicsError.invalid("quantitative H5MU bundle destination exists or is not local")
        }
    }

    public static func importH5MU(source: URL, plan: VivoH5MUQuantitativePlan,
                                  implementation: VivoFingerprint, to destination: URL) throws -> VivoQuantitativeH5MUReceipt {
        try VivoQuantitativeH5MUImport.validate(plan)
        try requireNew(destination)
        let planBytes = try VivoCanonicalJSON.encode(plan)
        guard planBytes.count <= maximumPlanBytes else { throw VivoOmicsError.limit("quantitative H5MU plan bytes") }
        let temporary = try staging(destination.deletingLastPathComponent())
        defer { try? FileManager.default.removeItem(at: temporary) }

        let sourceSnapshot = temporary.appendingPathComponent("original.h5mu")
        let sourceHash = try fingerprint(source, copyTo: sourceSnapshot)
        let dataset = try VivoQuantitativeH5MUImport.readSnapshot(sourceSnapshot, plan: plan)
        let datasetBytes = try VivoCanonicalJSON.encode(dataset)
        guard datasetBytes.count <= maximumDatasetBytes else { throw VivoOmicsError.limit("quantitative H5MU dataset bytes") }
        try planBytes.write(to: temporary.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        try datasetBytes.write(to: temporary.appendingPathComponent("dataset.json"), options: .withoutOverwriting)

        let export = temporary.appendingPathComponent("dataset.h5mu")
        try VivoQuantitativeH5MU.writeSnapshot(dataset, to: export)
        let h5muHash = try fingerprint(export)
        let receipt = VivoQuantitativeH5MUReceipt(source: sourceHash,
            plan: try VivoCanonicalJSON.fingerprint(planBytes), dataset: try VivoCanonicalJSON.fingerprint(datasetBytes),
            h5mu: h5muHash, implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: temporary.appendingPathComponent("receipt.json"),
                                                    options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }

    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> VivoQuantitativeH5MUReceipt {
        let receiptBytes = try read(directory, "receipt.json", maximum: maximumReceiptBytes)
        let receipt = try VivoCanonicalJSON.decode(VivoQuantitativeH5MUReceipt.self, from: receiptBytes)
        guard receipt.schemaVersion == 1, try VivoCanonicalJSON.encode(receipt) == receiptBytes,
              receipt.implementation == implementation else {
            throw VivoOmicsError.invalid("quantitative H5MU receipt")
        }
        let temporary = try staging(FileManager.default.temporaryDirectory)
        defer { try? FileManager.default.removeItem(at: temporary) }
        let planBytes = try read(directory, "plan.json", maximum: maximumPlanBytes)
        let datasetBytes = try read(directory, "dataset.json", maximum: maximumDatasetBytes)
        guard try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan,
              try VivoCanonicalJSON.fingerprint(datasetBytes) == receipt.dataset else {
            throw VivoOmicsError.invalid("quantitative H5MU metadata fingerprint")
        }
        let source = try fingerprint(directory.appendingPathComponent("original.h5mu"),
                                      copyTo: temporary.appendingPathComponent("original.h5mu"))
        guard source == receipt.source,
              try fingerprint(directory.appendingPathComponent("dataset.h5mu")) == receipt.h5mu else {
            throw VivoOmicsError.invalid("quantitative H5MU artifact fingerprint")
        }
        let plan = try VivoCanonicalJSON.decode(VivoH5MUQuantitativePlan.self, from: planBytes)
        let rebuilt = try VivoQuantitativeH5MUImport.readSnapshot(temporary.appendingPathComponent("original.h5mu"), plan: plan)
        guard try VivoCanonicalJSON.encode(rebuilt) == datasetBytes else {
            throw VivoOmicsError.invalid("quantitative H5MU source reconstruction differs")
        }
        let reconstructed = temporary.appendingPathComponent("reconstructed.h5mu")
        try VivoQuantitativeH5MU.writeSnapshot(rebuilt, to: reconstructed)
        guard try fingerprint(reconstructed) == receipt.h5mu else {
            throw VivoOmicsError.invalid("quantitative H5MU reconstruction differs")
        }
        return receipt
    }
}
