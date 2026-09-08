import Foundation

public struct VivoSingleCellRunReceipt: Codable, Sendable, Equatable {
    public let schemaVersion: Int
    public let input: VivoFingerprint
    public let result: VivoFingerprint
    public let implementation: VivoFingerprint
    public init(input: VivoFingerprint, result: VivoFingerprint, implementation: VivoFingerprint) {
        self.schemaVersion = 1; self.input = input; self.result = result; self.implementation = implementation
    }
}
private struct VivoSingleCellStoredResult: Codable, Sendable, Equatable {
    let schemaVersion: Int
    let input: VivoFingerprint
    let implementation: VivoFingerprint
    let report: VivoSingleCellReport
}

/// Uses the common immutable artifact store. No named reference is changed;
/// interrupted publication may leave reusable objects, never a success receipt.
public enum VivoSingleCellArtifacts {
    public static let inputKind = "vivo.singlecell-input-bundle-v1"
    public static let resultKind = "vivo.singlecell-stored-result-v1"
    public static func publish(input: VivoSingleCellInputBundle, implementation: VivoFingerprint,
                               store: VivoArtifactStore) async throws -> VivoSingleCellRunReceipt {
        let report = try VivoSingleCellCampaign.evaluate(input)
        let source = try VivoCanonicalJSON.encode(input)
        let sourceID = try VivoCanonicalJSON.fingerprint(source)
        let record = VivoSingleCellStoredResult(schemaVersion: 1, input: sourceID,
                                               implementation: implementation, report: report)
        let result = try VivoCanonicalJSON.encode(record)
        try Task.checkCancellation()
        let savedInput = try await store.put(data: source, kind: inputKind, mediaType: "application/json")
        guard savedInput.fingerprint == sourceID else { throw VivoOmicsError.invalid("stored source identity differs") }
        try Task.checkCancellation()
        let savedResult = try await store.put(data: result, kind: resultKind, mediaType: "application/json")
        try Task.checkCancellation()
        return .init(input: sourceID, result: savedResult.fingerprint, implementation: implementation)
    }
    public static func verify(receipt: VivoSingleCellRunReceipt, implementation: VivoFingerprint,
                              store: VivoArtifactStore) async throws -> VivoSingleCellReport {
        guard receipt.schemaVersion == 1, receipt.implementation == implementation else {
            throw VivoOmicsError.invalid("unsupported receipt or changed executable/OS; do not silently reinterpret another implementation")
        }
        let sourceDescriptor = try await store.descriptor(for: receipt.input)
        let resultDescriptor = try await store.descriptor(for: receipt.result)
        guard sourceDescriptor.kind == inputKind, resultDescriptor.kind == resultKind,
              sourceDescriptor.mediaType == "application/json", resultDescriptor.mediaType == "application/json" else {
            throw VivoOmicsError.invalid("receipt artifact kind mismatch")
        }
        let source = try await store.data(for: receipt.input, maximumBytes: 128 * 1_024 * 1_024)
        let result = try await store.data(for: receipt.result, maximumBytes: 512 * 1_024 * 1_024)
        let input = try VivoCanonicalJSON.decode(VivoSingleCellInputBundle.self, from: source)
        let record = try VivoCanonicalJSON.decode(VivoSingleCellStoredResult.self, from: result)
        guard record.schemaVersion == 1, record.input == receipt.input, record.implementation == implementation else {
            throw VivoOmicsError.invalid("result is not bound to the receipt source and implementation")
        }
        try Task.checkCancellation()
        try VivoSingleCellCampaign.verify(record.report, input: input)
        return record.report
    }
}
