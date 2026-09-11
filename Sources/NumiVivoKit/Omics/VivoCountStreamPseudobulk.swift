import Foundation
import CryptoKit

/// Admission contract for a canonical raw-count stream. The source declaration
/// describes upstream provenance; it is not a fingerprint of unread source bytes.
public struct VivoCountStreamPlan: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let metadata: VivoSingleCellCountMetadata
    public let sourceDeclaration: String
    public let rowNonzeros: [Int]
    /// Optional independent totals for the exact retained feature axis. Historical
    /// pre-filter QC totals must not be relabeled as totals of this matrix.
    public let rowTotals: [UInt64]?
    public init(metadata: VivoSingleCellCountMetadata, sourceDeclaration: String,
                rowNonzeros: [Int], rowTotals: [UInt64]? = nil) {
        schemaVersion = 1; self.metadata = metadata; self.sourceDeclaration = sourceDeclaration
        self.rowNonzeros = rowNonzeros; self.rowTotals = rowTotals
    }
    private enum CodingKeys: String, CodingKey { case schemaVersion, metadata, sourceDeclaration, rowNonzeros, rowTotals }
    public init(from decoder: Decoder) throws {
        try vivoOmicsRejectUnknownKeys(decoder, allowed: ["schemaVersion", "metadata", "sourceDeclaration", "rowNonzeros", "rowTotals"])
        let c = try decoder.container(keyedBy: CodingKeys.self)
        schemaVersion = try c.decode(Int.self, forKey: .schemaVersion)
        metadata = try c.decode(VivoSingleCellCountMetadata.self, forKey: .metadata)
        sourceDeclaration = try c.decode(String.self, forKey: .sourceDeclaration)
        rowNonzeros = try c.decode([Int].self, forKey: .rowNonzeros)
        rowTotals = try c.decodeIfPresent([UInt64].self, forKey: .rowTotals)
    }
    public func validate() throws -> Int {
        try metadata.validate(limits: VivoH5ADPseudobulk.sourceLimits)
        guard schemaVersion == 1, !sourceDeclaration.isEmpty, sourceDeclaration.utf8.count <= 16_384,
              rowNonzeros.count == metadata.cells.count, rowTotals.map({ $0.count == metadata.cells.count }) ?? true else {
            throw VivoOmicsError.invalid("count stream plan schema, provenance or row axes")
        }
        var entries = 0
        for i in rowNonzeros.indices {
            let n = rowNonzeros[i]
            guard n >= 0, n <= metadata.features.count,
                  n <= VivoH5ADPseudobulk.sourceLimits.maximumNonzeros - entries else {
                throw VivoOmicsError.invalid("count stream row cardinality, total or entry limit")
            }
            if let totals = rowTotals {
                guard UInt64(n) <= totals[i], (n == 0) == (totals[i] == 0) else {
                    throw VivoOmicsError.invalid("count stream row total inconsistent with nonzeros")
                }
            }
            entries += n
        }
        return entries
    }
}

public struct VivoCountStreamReport: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let method: String
    public let canonicalNonzeros: Int
    public let quality: [VivoCellQuality]
    public let pseudobulk: VivoPseudobulkCounts
}

public struct VivoCountStreamReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let encoding: String
    public let streamBytes: Int
    public let stream: VivoFingerprint
    public let plan: VivoFingerprint
    public let report: VivoFingerprint
    public let implementation: VivoFingerprint
}

/// Matrix entries occupy a bounded buffer. Metadata, per-cell QC and bounded
/// donor aggregates remain resident. Replay requires the original count stream.
public enum VivoCountStreamPseudobulk {
    public static let encoding = "row-u32-feature-u32-count-u64-little-endian-v1"
    public static let maximumChunkBytes = 1_048_576
    public static let maximumPlanBytes = 268_435_456
    static let maximumReportBytes = 536_870_912

    static func evaluate(plan: VivoCountStreamPlan, read: () throws -> Data) throws -> (VivoCountStreamReport, VivoFingerprint, Int) {
        let expected = try plan.validate(), accumulator = try VivoH5ADPseudobulk.Accumulator(plan.metadata)
        var buffer = Data(), hash = SHA256(), entries = 0, previousRow = -1, previousFeature = -1
        while true {
            try Task.checkCancellation()
            let chunk = try read()
            guard chunk.count <= maximumChunkBytes else { throw VivoOmicsError.limit("count stream read exceeds 1 MiB") }
            if chunk.isEmpty { break }
            guard chunk.count <= expected * 16 - entries * 16 - buffer.count else {
                throw VivoOmicsError.invalid("count stream exceeds declared entry count")
            }
            hash.update(data: chunk); buffer.append(chunk)
            let complete = buffer.count / 16 * 16
            try buffer.withUnsafeBytes { bytes in
                for offset in stride(from: 0, to: complete, by: 16) {
                    let packed = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset, as: UInt64.self))
                    let row = Int(packed & 0xffff_ffff), feature = Int(packed >> 32)
                    let count = UInt64(littleEndian: bytes.loadUnaligned(fromByteOffset: offset + 8, as: UInt64.self))
                    guard row < plan.metadata.cells.count, feature < plan.metadata.features.count, count > 0,
                          row > previousRow || (row == previousRow && feature > previousFeature) else {
                        throw VivoOmicsError.invalid("count stream requires ordered unique positive in-axis records")
                    }
                    try accumulator.add(row: row, feature: feature, count: count)
                    previousRow = row; previousFeature = feature; entries += 1
                }
            }
            buffer = Data(buffer.suffix(buffer.count - complete))
        }
        guard buffer.isEmpty, entries == expected else { throw VivoOmicsError.invalid("count stream truncated or missing entries") }
        let accumulated = try accumulator.finish(version: "", contrasts: [])
        for i in accumulated.quality.indices {
            guard plan.rowTotals.map({ accumulated.quality[i].totalCounts == $0[i] }) ?? true,
                  accumulated.quality[i].detectedFeatures == plan.rowNonzeros[i] else {
                throw VivoOmicsError.invalid("count stream differs from declared row total or nonzeros")
            }
        }
        let report = VivoCountStreamReport(schemaVersion: 1, method: "canonical-count-stream-pseudobulk-v1",
            canonicalNonzeros: entries, quality: accumulated.quality, pseudobulk: accumulated.pseudobulk)
        return (report, try VivoFingerprint(bytes: Array(hash.finalize())), entries * 16)
    }
    public static func publish(plan: VivoCountStreamPlan, input: FileHandle, implementation: VivoFingerprint,
                               to destination: URL) throws -> VivoCountStreamReceipt {
        _ = try plan.validate()
        let planBytes = try VivoCanonicalJSON.encode(plan)
        guard planBytes.count <= maximumPlanBytes else { throw VivoOmicsError.limit("count stream plan exceeds 256 MiB") }
        guard !FileManager.default.fileExists(atPath: destination.path) else { throw VivoOmicsError.invalid("count stream output already exists") }
        let staging = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-count-stream-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: staging, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: staging) }
        let (report, stream, bytes) = try evaluate(plan: plan) { try input.read(upToCount: maximumChunkBytes) ?? Data() }
        let reportBytes = try VivoCanonicalJSON.encode(report)
        guard reportBytes.count <= maximumReportBytes else { throw VivoOmicsError.limit("count stream report exceeds 512 MiB") }
        let receipt = try VivoCountStreamReceipt(schemaVersion: 1, encoding: encoding, streamBytes: bytes, stream: stream,
            plan: VivoCanonicalJSON.fingerprint(planBytes), report: VivoCanonicalJSON.fingerprint(reportBytes), implementation: implementation)
        try planBytes.write(to: staging.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        try reportBytes.write(to: staging.appendingPathComponent("report.json"), options: .withoutOverwriting)
        try VivoCanonicalJSON.encode(receipt).write(to: staging.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: staging, to: destination)
        return receipt
    }
    /// The caller supplies the complete stream again. A matching receipt alone
    /// never establishes reconstruction or upstream HDF5 validity.
    public static func verify(_ directory: URL, input: FileHandle, implementation: VivoFingerprint) throws -> VivoCountStreamReceipt {
        func read(_ name: String, _ limit: Int) throws -> Data {
            try VivoSingleCellCampaignIO.readDocument(directory.appendingPathComponent(name), maximumBytes: limit)
        }
        let receipt = try VivoCanonicalJSON.decode(VivoCountStreamReceipt.self, from: read("receipt.json", 65_536))
        let planBytes = try read("plan.json", maximumPlanBytes), reportBytes = try read("report.json", maximumReportBytes)
        guard receipt.schemaVersion == 1, receipt.encoding == encoding, receipt.implementation == implementation,
              try VivoCanonicalJSON.fingerprint(planBytes) == receipt.plan,
              try VivoCanonicalJSON.fingerprint(reportBytes) == receipt.report else {
            throw VivoOmicsError.invalid("count stream receipt schema, implementation or hash mismatch")
        }
        let plan = try VivoCanonicalJSON.decode(VivoCountStreamPlan.self, from: planBytes)
        let (report, stream, bytes) = try evaluate(plan: plan) { try input.read(upToCount: maximumChunkBytes) ?? Data() }
        guard stream == receipt.stream, bytes == receipt.streamBytes,
              try VivoCanonicalJSON.encode(report) == reportBytes else {
            throw VivoOmicsError.invalid("count stream report does not reconstruct")
        }
        return receipt
    }
}
