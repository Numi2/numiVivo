import Foundation

public struct VivoNeoantigenReceipt: Codable, Sendable, Equatable {
    public var schema: String
    public var manifest: VivoFingerprint
    public var sourceTSV: VivoFingerprint
    public var report: VivoFingerprint
    public var implementation: VivoFingerprint
}

/// Uses the existing immutable artifact store, without publishing mutable refs.
/// Raw manifest and TSV bytes are preserved; verification reconstructs the report.
public enum VivoNeoantigenArtifacts {
    public static func publish(manifestData: Data, tsv: Data, implementation: VivoFingerprint,
                               store: VivoArtifactStore) async throws -> (receipt: VivoNeoantigenReceipt, report: VivoNeoantigenReport) {
        guard manifestData.count <= 128 * 1024 else { throw VivoNeoantigenError.invalid("Case manifest exceeds 128 KiB.") }
        let manifest = try VivoCanonicalJSON.decode(VivoNeoantigenCase.self, from: manifestData)
        let report = try VivoNeoantigenWorkbench.analyze(manifest: manifest, tsv: tsv, implementationSHA256: implementation.hex)
        let input = try await store.put(data: manifestData, kind: "neoantigen.case", mediaType: "application/json")
        let source = try await store.put(data: tsv, kind: "neoantigen.source", mediaType: "text/tab-separated-values")
        let result = try await store.put(data: VivoCanonicalJSON.encode(report), kind: "neoantigen.report", mediaType: "application/json")
        let receipt = VivoNeoantigenReceipt(schema: "numivivo.org/neoantigen-receipt/v1", manifest: input.fingerprint,
            sourceTSV: source.fingerprint, report: result.fingerprint, implementation: implementation)
        _ = try await store.put(data: VivoCanonicalJSON.encode(receipt), kind: "neoantigen.receipt", mediaType: "application/json")
        return (receipt, report)
    }

    public static func verify(_ receipt: VivoNeoantigenReceipt, implementation: VivoFingerprint,
                              store: VivoArtifactStore) async throws -> VivoNeoantigenReport {
        guard receipt.schema == "numivivo.org/neoantigen-receipt/v1", receipt.implementation == implementation else {
            throw VivoNeoantigenError.invalid("Receipt schema or implementation identity differs; use the originating binary or create a new import.")
        }
        let manifestData = try await store.data(for: receipt.manifest, maximumBytes: 128 * 1024)
        let tsv = try await store.data(for: receipt.sourceTSV, maximumBytes: VivoNeoantigenWorkbench.maximumTSVBytes)
        let saved = try await store.data(for: receipt.report, maximumBytes: 128 * 1024 * 1024)
        let manifest = try VivoCanonicalJSON.decode(VivoNeoantigenCase.self, from: manifestData)
        let regenerated = try VivoNeoantigenWorkbench.analyze(manifest: manifest, tsv: tsv, implementationSHA256: implementation.hex)
        guard try VivoCanonicalJSON.encode(regenerated) == saved else {
            throw VivoNeoantigenError.invalid("Archived report does not reconstruct from the archived inputs and implementation.")
        }
        return regenerated
    }

    public static func record(_ review: VivoNeoantigenReview, receipt: VivoNeoantigenReceipt,
                              implementation: VivoFingerprint, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        let report = try await verify(receipt, implementation: implementation, store: store)
        try VivoNeoantigenWorkbench.validate(review, against: report)
        return try await store.put(data: VivoCanonicalJSON.encode(review), kind: "neoantigen.research-review", mediaType: "application/json")
    }
}
