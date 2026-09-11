import Foundation

/// Explicit source and execution identities allow analysis of an earlier qualified
/// count bundle without attributing that count execution to the new analysis binary.
public struct VivoFileExpressionPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let sourceReceipt: VivoFingerprint
    public let sourceImplementation: VivoFingerprint
    public let contrast: VivoOmicsExpressionContrast
    public init(sourceReceipt: VivoFingerprint, sourceImplementation: VivoFingerprint, contrast: VivoOmicsExpressionContrast) {
        schemaVersion = 1; self.sourceReceipt = sourceReceipt; self.sourceImplementation = sourceImplementation; self.contrast = contrast
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, sourceReceipt, sourceImplementation, contrast }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "sourceReceipt", "sourceImplementation", "contrast"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion); sourceReceipt = try c.decode(VivoFingerprint.self, forKey: .sourceReceipt)
        sourceImplementation = try c.decode(VivoFingerprint.self, forKey: .sourceImplementation); contrast = try c.decode(VivoOmicsExpressionContrast.self, forKey: .contrast)
        try validate()
    }
    public func validate() throws {
        guard schemaVersion == 1, sourceReceipt.bytes.count == 32, sourceImplementation.bytes.count == 32 else { throw VivoOmicsError.invalid("file expression source identity") }
        try contrast.validate()
    }
}
public struct VivoFileExpressionReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let sourceReceipt: VivoFingerprint
    public let sourceImplementation: VivoFingerprint
    public let plan: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

public enum VivoFileExpression {
    static let maximumReportBytes = 268_435_456
    /// Validate all QC/membership rows and group library totals. This does not
    /// reconstruct feature-level counts from an unavailable original count stream.
    static func evaluate(_ source: VivoFileCountSnapshot, plan: VivoFileExpressionPlan) throws -> VivoOmicsExpressionResult {
        try plan.validate()
        guard try source.receiptFingerprint() == plan.sourceReceipt, source.receipt.implementation == plan.sourceImplementation else {
            throw VivoOmicsError.invalid("file expression input differs from declared count source")
        }
        let groups = source.report.groups
        var cells = [Int](repeating: 0, count: groups.count), totals = [UInt64](repeating: 0, count: groups.count)
        for row in 0..<source.axis.header.cellCount {
            if row % 4096 == 0 { try Task.checkCancellation() }
            try vivoAxisPool {
                let item = try source.quality(row); cells[item.group] += 1
                totals[item.group] = try vivoOmicsSum(totals[item.group], item.quality.totalCounts)
            }
        }
        let counts = source.report.matrix
        for group in groups.indices {
            var total: UInt64 = 0
            for k in counts.rowOffsets[group]..<counts.rowOffsets[group + 1] { total = try vivoOmicsSum(total, counts.counts[k]) }
            guard cells[group] == groups[group].sourceCellCount, totals[group] == total else {
                throw VivoOmicsError.invalid("file expression group membership or library total")
            }
        }
        let observations = groups.enumerated().map { VivoOmicsDesignObservation($0.element, receipt: plan.sourceReceipt, index: $0.offset) }
        return try VivoPseudobulkDifferentialExpression.evaluate(metadata: source.axis.header.metadata, observations: observations, counts: counts, contrast: plan.contrast)
    }
    public static func run(source directory: URL, plan: VivoFileExpressionPlan) throws -> VivoOmicsExpressionResult {
        try plan.validate()
        return try evaluate(VivoFileCountSnapshot.open(directory, implementation: plan.sourceImplementation), plan: plan)
    }
    public static func publish(source directory: URL, plan: VivoFileExpressionPlan, implementation: VivoFingerprint, to destination: URL) throws -> VivoFileExpressionReceipt {
        try plan.validate()
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("file expression output exists") }
        let source = try VivoFileCountSnapshot.open(directory, implementation: plan.sourceImplementation)
        let result = try evaluate(source, plan: plan)
        let planBytes = try VivoCanonicalJSON.encode(plan)
        guard planBytes.count <= 131_072 else { throw VivoOmicsError.limit("file expression document size") }
        let staging = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-file-expression-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: staging, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        try source.copy(to: staging.appendingPathComponent("source"))
        let reportURL = staging.appendingPathComponent("report.json")
        // The private staging directory is owned exclusively by this publication.
        guard FileManager.default.createFile(atPath: reportURL.path, contents: nil, attributes: [.posixPermissions: 0o600]) else {
            throw VivoOmicsError.invalid("cannot create expression report")
        }
        let reportFile = try FileHandle(forWritingTo: reportURL)
        defer { try? reportFile.close() }
        let reportHash = try VivoExpressionReportJSON.fingerprint(result, maximumBytes: maximumReportBytes) { try reportFile.write(contentsOf: $0) }
        try reportFile.synchronize(); try reportFile.close()
        try planBytes.write(to: staging.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        let receipt = try VivoFileExpressionReceipt(schemaVersion: 1, method: "file-membership-pseudobulk-expression-v1", sourceReceipt: plan.sourceReceipt,
            sourceImplementation: plan.sourceImplementation, plan: VivoCanonicalJSON.fingerprint(planBytes), report: reportHash, implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation(); try FileManager.default.moveItem(at: staging, to: destination); return receipt
    }
    public static func verify(_ directory: URL, implementation: VivoFingerprint) throws -> VivoOmicsExpressionResult {
        let snapshot = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-file-expression-read-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: snapshot, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: snapshot) }
        func copy(_ name: String, _ limit: Int) throws -> VivoFingerprint {
            try VivoOmicsFileSnapshot.fingerprint(directory.appendingPathComponent(name), copyTo: snapshot.appendingPathComponent(name), maximumBytes: limit)
        }
        _ = try copy("receipt.json", 65_536)
        let receipt = try VivoCanonicalJSON.decode(VivoFileExpressionReceipt.self, from: Data(contentsOf: snapshot.appendingPathComponent("receipt.json")))
        guard receipt.schemaVersion == 1, receipt.method == "file-membership-pseudobulk-expression-v1", receipt.implementation == implementation,
              try copy("plan.json", 131_072) == receipt.plan, try copy("report.json", maximumReportBytes) == receipt.report else {
            throw VivoOmicsError.invalid("file expression receipt or document fingerprint")
        }
        let plan = try VivoCanonicalJSON.decode(VivoFileExpressionPlan.self, from: Data(contentsOf: snapshot.appendingPathComponent("plan.json")))
        guard plan.sourceReceipt == receipt.sourceReceipt, plan.sourceImplementation == receipt.sourceImplementation else { throw VivoOmicsError.invalid("file expression source binding") }
        let result = try run(source: directory.appendingPathComponent("source"), plan: plan)
        guard try VivoExpressionReportJSON.fingerprint(result, maximumBytes: maximumReportBytes) == receipt.report else { throw VivoOmicsError.invalid("file expression does not reconstruct") }
        return result
    }
}
