import Foundation
import CryptoKit

public struct VivoCellTypistStreamReceipt: Codable, Sendable {
    public let schemaVersion: Int
    public let cells: Int
    public let records: Int
    public let modelSHA256: String
    public let planSHA256: String
    public let streamSHA256: String
    public let resultsSHA256: String
    public let implementation: VivoFingerprint
}

public enum VivoCellTypistReferenceIO {
    private struct Row: Encodable {
        let sourceRow: Int
        let label: String
        let decisions: [Double]
        let probabilities: [Double]
    }
    private static func hex(_ digest: SHA256.Digest) -> String {
        digest.map { String(format: "%02x", $0) }.joined()
    }

    /// Output remains in a private staging directory until every row, exact
    /// input digest and output write succeeds. The embedded plan binds row IDs.
    public static func publish(modelData: Data, plan: VivoCountStreamPlan,
                               expectedStreamSHA256: String, implementation: VivoFingerprint,
                               read: () throws -> Data, to destination: URL) throws -> VivoCellTypistStreamReceipt {
        let expectedRecords = try plan.validate()
        guard modelData.count <= 67_108_864, expectedStreamSHA256.count == 64,
              expectedStreamSHA256.allSatisfy({ "0123456789abcdef".contains($0) }),
              destination.isFileURL, !FileManager.default.fileExists(atPath: destination.path) else {
            throw VivoOmicsError.invalid("CellTypist model size, stream digest or output path")
        }
        let model = try VivoCanonicalJSON.decode(VivoCellTypistReferenceModel.self, from: modelData)
        let predictor = try VivoCellTypistReference(model: model, sourceFeatures: plan.metadata.features.map(\.id))
        let planData = try VivoCanonicalJSON.encode(plan)
        let temporary = destination.deletingLastPathComponent().appendingPathComponent(".numivivo-celltypist-" + UUID().uuidString)
        try FileManager.default.createDirectory(at: temporary, withIntermediateDirectories: false, attributes: [.posixPermissions: 0o700])
        defer { try? FileManager.default.removeItem(at: temporary) }
        try modelData.write(to: temporary.appendingPathComponent("model.json"), options: .withoutOverwriting)
        try planData.write(to: temporary.appendingPathComponent("plan.json"), options: .withoutOverwriting)
        let output = temporary.appendingPathComponent("results.jsonl")
        try Data().write(to: output, options: .withoutOverwriting)
        let handle = try FileHandle(forWritingTo: output)
        defer { try? handle.close() }
        var inputHash = SHA256(), resultHash = SHA256()
        let records = try predictor.predictStream(cellCount: plan.metadata.cells.count, maximumRecords: expectedRecords, read: {
            let data = try read(); inputHash.update(data: data); return data
        }, validateRow: { row, indices, counts in
            guard indices.count == plan.rowNonzeros[row] else { throw VivoOmicsError.invalid("CellTypist stream row cardinality") }
            if let totals = plan.rowTotals {
                var sum: UInt64 = 0
                for count in counts {
                    let next = sum.addingReportingOverflow(count)
                    guard !next.overflow else { throw VivoOmicsError.invalid("CellTypist row overflow") }
                    sum = next.partialValue
                }
                guard sum == totals[row] else { throw VivoOmicsError.invalid("CellTypist stream row total") }
            }
        }, emit: { index, result in
            var bytes = try VivoCanonicalJSON.encode(Row(sourceRow: index, label: result.label,
                decisions: result.decisions, probabilities: result.probabilities))
            bytes.append(10); try handle.write(contentsOf: bytes); resultHash.update(data: bytes)
        })
        let streamHash = hex(inputHash.finalize())
        guard records == expectedRecords, streamHash == expectedStreamSHA256 else {
            throw VivoOmicsError.invalid("CellTypist stream differs from frozen input")
        }
        try handle.synchronize(); try handle.close()
        let receipt = VivoCellTypistStreamReceipt(schemaVersion: 1, cells: plan.metadata.cells.count,
            records: records, modelSHA256: hex(SHA256.hash(data: modelData)),
            planSHA256: hex(SHA256.hash(data: planData)), streamSHA256: streamHash,
            resultsSHA256: hex(resultHash.finalize()), implementation: implementation)
        try VivoCanonicalJSON.encode(receipt).write(to: temporary.appendingPathComponent("receipt.json"), options: .withoutOverwriting)
        try Task.checkCancellation()
        try FileManager.default.moveItem(at: temporary, to: destination)
        return receipt
    }
}
