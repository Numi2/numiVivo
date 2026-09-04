import CryptoKit
import Foundation

public enum VivoArtifactKind: String, Codable, Sendable, CaseIterable {
    case program, hostContext, experiment, coupling, calibration, result, safetyCase, checkpoint, evidenceSnapshot
}
public enum VivoEvidenceClass: String, Codable, Sendable, CaseIterable {
    case observed, derived, calibrated, inferred, assumed, hypothetical
}
public struct VivoEvidenceReference: Codable, Sendable, Equatable, Hashable {
    public var classification: VivoEvidenceClass
    public var sourceURI: String?
    public var datasetIdentifier: String?
    public var citation: String?
    public var context: String?
    public var note: String?
    public var snapshotFingerprint: VivoFingerprint?
    public var snapshotEntryIdentifier: String?
    public init(classification: VivoEvidenceClass, sourceURI: String? = nil,
                datasetIdentifier: String? = nil, citation: String? = nil,
                context: String? = nil, note: String? = nil,
                snapshotFingerprint: VivoFingerprint? = nil, snapshotEntryIdentifier: String? = nil) {
        self.classification = classification
        self.sourceURI = sourceURI
        self.datasetIdentifier = datasetIdentifier
        self.citation = citation
        self.context = context
        self.note = note
        self.snapshotFingerprint = snapshotFingerprint
        self.snapshotEntryIdentifier = snapshotEntryIdentifier
    }
}
public struct VivoArtifactMetadata: Codable, Sendable, Equatable, Hashable {
    public var identifier: String
    public var version: String
    public var namespace: String
    public var title: String
    public var summary: String?
    public var creators: [String]
    public var labels: [String: String]
    public init(identifier: String, version: String,
                namespace: String = "https://numivivo.org/artifacts", title: String,
                summary: String? = nil, creators: [String] = [], labels: [String: String] = [:]) {
        self.identifier = identifier
        self.version = version
        self.namespace = namespace
        self.title = title
        self.summary = summary
        self.creators = creators
        self.labels = labels
    }
}
public struct VivoProvenance: Codable, Sendable, Equatable, Hashable {
    public var createdAt: Date
    public var createdBy: String
    public var softwareVersion: String
    public var sourceArtifacts: [VivoFingerprint]
    public var evidence: [VivoEvidenceReference]
    public var environment: [String: String]
    public init(createdAt: Date = Date(), createdBy: String, softwareVersion: String,
                sourceArtifacts: [VivoFingerprint] = [], evidence: [VivoEvidenceReference] = [],
                environment: [String: String] = [:]) {
        self.createdAt = createdAt
        self.createdBy = createdBy
        self.softwareVersion = softwareVersion
        self.sourceArtifacts = sourceArtifacts
        self.evidence = evidence
        self.environment = environment
    }
}
public struct VivoArtifactSignature: Codable, Sendable, Equatable, Hashable {
    public enum Algorithm: String, Codable, Sendable { case ed25519 }
    public let algorithm: Algorithm
    public let publicKey: Data
    public let signature: Data
    public init(algorithm: Algorithm = .ed25519, publicKey: Data, signature: Data) {
        self.algorithm = algorithm
        self.publicKey = publicKey
        self.signature = signature
    }
}
public struct VivoArtifactEnvelope<Payload: Codable & Sendable>: Codable, Sendable {
    public var apiVersion: String
    public var kind: VivoArtifactKind
    public var metadata: VivoArtifactMetadata
    public var provenance: VivoProvenance
    public var payload: Payload
    public var signature: VivoArtifactSignature?
    public init(apiVersion: String = "numivivo.org/v1alpha1", kind: VivoArtifactKind,
                metadata: VivoArtifactMetadata, provenance: VivoProvenance,
                payload: Payload, signature: VivoArtifactSignature? = nil) {
        self.apiVersion = apiVersion
        self.kind = kind
        self.metadata = metadata
        self.provenance = provenance
        self.payload = payload
        self.signature = signature
    }
    public func canonicalData(includeSignature: Bool = true) throws -> Data {
        var copy = self
        if !includeSignature { copy.signature = nil }
        return try VivoCanonicalJSON.encode(copy)
    }
    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(canonicalData(includeSignature: false))
    }
    public mutating func sign(using privateKey: Curve25519.Signing.PrivateKey) throws {
        let bytes = try canonicalData(includeSignature: false)
        signature = VivoArtifactSignature(publicKey: privateKey.publicKey.rawRepresentation,
                                           signature: try privateKey.signature(for: bytes))
    }
    public func verifySignature() throws -> Bool {
        guard let signature, signature.algorithm == .ed25519 else { return false }
        let publicKey = try Curve25519.Signing.PublicKey(rawRepresentation: signature.publicKey)
        return try publicKey.isValidSignature(signature.signature, for: canonicalData(includeSignature: false))
    }
}
public enum VivoCanonicalJSON {
    public static func encoder() -> JSONEncoder {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        encoder.dateEncodingStrategy = .iso8601
        encoder.dataEncodingStrategy = .base64
        encoder.nonConformingFloatEncodingStrategy = .throw
        return encoder
    }
    public static func decoder() -> JSONDecoder {
        let decoder = JSONDecoder()
        decoder.dateDecodingStrategy = .iso8601
        decoder.dataDecodingStrategy = .base64
        decoder.nonConformingFloatDecodingStrategy = .throw
        return decoder
    }
    public static func encode<T: Encodable>(_ value: T) throws -> Data { try encoder().encode(value) }
    public static func decode<T: Decodable>(_ type: T.Type, from data: Data) throws -> T {
        try decoder().decode(type, from: data)
    }
    public static func fingerprint(_ data: Data) throws -> VivoFingerprint {
        try VivoFingerprint(bytes: Array(SHA256.hash(data: data)))
    }
}

// A generic helper type must live at file scope, not inside generic append().
// Its field names retain the existing digest-chain encoding and hash material.
private struct VivoDigestChainInput<Record: Encodable>: Encodable {
    let previous: VivoFingerprint
    let index: UInt64
    let record: Record
}
public struct VivoDigestChain: Codable, Sendable, Equatable {
    public let genesis: VivoFingerprint
    public private(set) var head: VivoFingerprint
    public private(set) var length: UInt64
    public init(genesis: VivoFingerprint) {
        self.genesis = genesis
        head = genesis
        length = 0
    }
    @discardableResult
    public mutating func append<T: Encodable>(_ record: T) throws -> VivoFingerprint {
        guard length < UInt64.max else {
            throw VivoArtifactValidationError.invalid("digest-chain index overflow")
        }
        let data = try VivoCanonicalJSON.encode(VivoDigestChainInput(previous: head, index: length, record: record))
        let next = try VivoCanonicalJSON.fingerprint(data)
        head = next
        length += 1
        return next
    }
}
