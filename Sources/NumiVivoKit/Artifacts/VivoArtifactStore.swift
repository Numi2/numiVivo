import Foundation

public struct VivoStoredArtifact: Codable, Sendable, Equatable, Hashable {
    public let fingerprint: VivoFingerprint
    public let kind: String
    public let mediaType: String
    public let byteCount: UInt64
    public let objectPath: String
    public let createdAt: Date
    public let attributes: [String: String]
    public init(fingerprint: VivoFingerprint, kind: String, mediaType: String, byteCount: UInt64,
                objectPath: String, createdAt: Date, attributes: [String: String]) {
        self.fingerprint = fingerprint; self.kind = kind; self.mediaType = mediaType
        self.byteCount = byteCount; self.objectPath = objectPath; self.createdAt = createdAt; self.attributes = attributes
    }
}
public struct VivoArtifactReference: Codable, Sendable, Equatable, Hashable {
    public let name: String
    public let artifact: VivoStoredArtifact
    public let updatedAt: Date
    public init(name: String, artifact: VivoStoredArtifact, updatedAt: Date = Date()) {
        self.name = name; self.artifact = artifact; self.updatedAt = updatedAt
    }
}
public enum VivoArtifactStoreError: Error, Sendable, CustomStringConvertible {
    case invalidRoot(String), invalidDescriptor(String)
    case objectMissing(VivoFingerprint), integrityFailure(VivoFingerprint)
    case referenceMissing(String), io(String)
    public var description: String {
        switch self {
        case .invalidRoot(let value): return "Invalid artifact root: \(value)"
        case .invalidDescriptor(let value): return "Invalid artifact descriptor: \(value)"
        case .objectMissing(let value): return "Artifact object \(value.hex) is missing."
        case .integrityFailure(let value): return "Artifact \(value.hex) failed integrity verification."
        case .referenceMissing(let value): return "Artifact reference '\(value)' does not exist."
        case .io(let value): return "Artifact I/O failed: \(value)"
        }
    }
}
public struct VivoArtifactStoreLimits: Sendable {
    public let maximumObjectBytes: Int
    public let maximumMetadataBytes: Int
    public let maximumListEntries: Int
    public init(maximumObjectBytes: Int = 512 * 1024 * 1024, maximumMetadataBytes: Int = 1024 * 1024,
                maximumListEntries: Int = 100_000) {
        self.maximumObjectBytes = maximumObjectBytes
        self.maximumMetadataBytes = maximumMetadataBytes
        self.maximumListEntries = maximumListEntries
    }
}

private struct VivoReferenceDeletion: Codable {
    let schema: String
    let name: String
    let deletedAt: Date
}

/// Objects and descriptors are immutable, no-clobber entries. Named references
/// are mutable pointers. A digest establishes integrity, not authorship or trust.
public actor VivoArtifactStore {
    public let rootURL: URL
    private let files: VivoRootedFileStore
    private let limits: VivoArtifactStoreLimits
    public init(rootURL: URL, createIfNeeded: Bool = true, limits: VivoArtifactStoreLimits = .init()) throws {
        guard limits.maximumObjectBytes > 0, limits.maximumMetadataBytes > 0, limits.maximumListEntries > 0 else {
            throw VivoArtifactStoreError.invalidRoot("I/O limits must be positive")
        }
        files = try VivoRootedFileStore(rootURL: rootURL, createIfNeeded: createIfNeeded)
        self.rootURL = files.rootURL
        self.limits = limits
        if createIfNeeded {
            for path in ["objects/sha256", "descriptors/sha256", "refs/v2", "bundles"] { try files.createDirectory(path) }
        }
    }
    @discardableResult
    public func put(data: Data, kind: String, mediaType: String, attributes: [String: String] = [:]) throws -> VivoStoredArtifact {
        try validateMetadata(kind: kind, mediaType: mediaType, attributes: attributes)
        guard data.count <= limits.maximumObjectBytes else { throw VivoArtifactStoreError.invalidDescriptor("object exceeds configured size limit") }
        let fingerprint = try VivoCanonicalJSON.fingerprint(data)
        let relative = objectRelativePath(fingerprint)
        if try !files.writeFile(data, relative: relative, immutable: true) {
            let existing = try self.data(for: fingerprint)
            guard existing.count == data.count else { throw VivoArtifactStoreError.integrityFailure(fingerprint) }
        }
        let descriptor = VivoStoredArtifact(fingerprint: fingerprint, kind: kind, mediaType: mediaType,
                                             byteCount: UInt64(data.count), objectPath: relative, createdAt: Date(), attributes: attributes)
        let encoded = try VivoCanonicalJSON.encode(descriptor)
        guard encoded.count <= limits.maximumMetadataBytes else { throw VivoArtifactStoreError.invalidDescriptor("metadata exceeds size limit") }
        if try files.writeFile(encoded, relative: descriptorPath(fingerprint), immutable: true) {
            // Canonical JSON encodes dates at second precision. Return exactly
            // the persisted metadata, not an in-memory subsecond variant.
            return try VivoCanonicalJSON.decode(VivoStoredArtifact.self, from: encoded)
        }
        let existing = try self.descriptor(for: fingerprint)
        guard existing.byteCount == UInt64(data.count), existing.kind == kind,
              existing.mediaType == mediaType, existing.attributes == attributes else {
            throw VivoArtifactStoreError.invalidDescriptor("immutable content already has different metadata; encode a new provenance envelope")
        }
        return existing
    }
    @discardableResult
    public func put<Payload: Codable & Sendable>(envelope: VivoArtifactEnvelope<Payload>,
                                                 mediaType: String = "application/vnd.numivivo.artifact+json",
                                                 attributes: [String: String] = [:]) throws -> VivoStoredArtifact {
        try put(data: envelope.canonicalData(), kind: envelope.kind.rawValue, mediaType: mediaType,
                attributes: attributes.merging(["identifier": envelope.metadata.identifier, "version": envelope.metadata.version],
                                               uniquingKeysWith: { current, _ in current }))
    }
    public func contains(_ fingerprint: VivoFingerprint) -> Bool {
        (try? files.isRegularFile(objectRelativePath(fingerprint))) == true
    }
    public func descriptor(for fingerprint: VivoFingerprint) throws -> VivoStoredArtifact {
        do {
            let bytes = try files.readFile(descriptorPath(fingerprint), maximumBytes: limits.maximumMetadataBytes)
            let value = try VivoCanonicalJSON.decode(VivoStoredArtifact.self, from: bytes)
            try validateMetadata(kind: value.kind, mediaType: value.mediaType, attributes: value.attributes)
            guard value.fingerprint == fingerprint, value.objectPath == objectRelativePath(fingerprint),
                  value.byteCount <= UInt64(limits.maximumObjectBytes), value.createdAt.timeIntervalSince1970.isFinite else {
                throw VivoArtifactStoreError.invalidDescriptor("descriptor does not match its content-addressed path")
            }
            return value
        } catch VivoRootedFileStore.Failure.missing(_) { throw VivoArtifactStoreError.objectMissing(fingerprint) }
    }
    public func data(for fingerprint: VivoFingerprint, verify: Bool = true) throws -> Data {
        do {
            let data = try files.readFile(objectRelativePath(fingerprint), maximumBytes: limits.maximumObjectBytes)
            if verify, try VivoCanonicalJSON.fingerprint(data) != fingerprint { throw VivoArtifactStoreError.integrityFailure(fingerprint) }
            return data
        } catch VivoRootedFileStore.Failure.missing(_) { throw VivoArtifactStoreError.objectMissing(fingerprint) }
    }
    public func verify(_ fingerprint: VivoFingerprint) throws -> Bool {
        let descriptor = try descriptor(for: fingerprint)
        return try descriptor.byteCount == UInt64(data(for: fingerprint, verify: true).count)
    }
    public func list(kind: String? = nil) throws -> [VivoStoredArtifact] {
        var values: [VivoStoredArtifact] = []
        let base = "descriptors/sha256"
        let first: [String]
        do { first = try files.children(base, maximumCount: 256) }
        catch VivoRootedFileStore.Failure.missing(_) { return [] }
        for prefix in first {
            guard hex(prefix, count: 2) else { throw VivoArtifactStoreError.invalidDescriptor("unexpected descriptor shard") }
            for second in try files.children(base + "/" + prefix, maximumCount: 256) {
                guard hex(second, count: 2) else { throw VivoArtifactStoreError.invalidDescriptor("unexpected descriptor shard") }
                let directory = base + "/" + prefix + "/" + second
                let names = try files.children(directory, maximumCount: limits.maximumListEntries - values.count)
                for name in names {
                    guard name.hasSuffix(".json") else { throw VivoArtifactStoreError.invalidDescriptor("unexpected descriptor filename") }
                    let digest = String(name.dropLast(5))
                    guard hex(digest, count: 64), digest.hasPrefix(prefix + second) else {
                        throw VivoArtifactStoreError.invalidDescriptor("descriptor filename does not match its shard")
                    }
                    let fingerprint = try fingerprintFromHex(digest)
                    values.append(try descriptor(for: fingerprint))
                }
            }
        }
        return values.filter { kind == nil || $0.kind == kind }.sorted {
            $0.createdAt == $1.createdAt ? $0.fingerprint.hex < $1.fingerprint.hex : $0.createdAt < $1.createdAt
        }
    }
    public func materialize(_ fingerprint: VivoFingerprint, to destination: URL, verify: Bool = true) throws {
        guard destination.isFileURL else { throw VivoArtifactStoreError.io("destination must be a file URL") }
        let data = try self.data(for: fingerprint, verify: verify)
        let output = try VivoRootedFileStore(rootURL: destination.deletingLastPathComponent(), createIfNeeded: true)
        try output.writeFile(data, relative: destination.lastPathComponent, immutable: false)
    }
    @discardableResult
    public func setReference(_ name: String, to artifact: VivoStoredArtifact) throws -> VivoArtifactReference {
        let path = try referencePath(name)
        let persisted = try descriptor(for: artifact.fingerprint)
        guard persisted == artifact else { throw VivoArtifactStoreError.invalidDescriptor("reference supplied forged or stale object metadata") }
        guard try verify(persisted.fingerprint) else { throw VivoArtifactStoreError.integrityFailure(persisted.fingerprint) }
        let reference = VivoArtifactReference(name: name, artifact: persisted)
        let data = try VivoCanonicalJSON.encode(reference)
        guard data.count <= limits.maximumMetadataBytes else { throw VivoArtifactStoreError.invalidDescriptor("reference exceeds size limit") }
        try files.writeFile(data, relative: path, immutable: false)
        return try VivoCanonicalJSON.decode(VivoArtifactReference.self, from: data)
    }
    public func reference(_ name: String) throws -> VivoArtifactReference {
        let path = try referencePath(name)
        let bytes: Data
        do { bytes = try files.readFile(path, maximumBytes: limits.maximumMetadataBytes) }
        catch VivoRootedFileStore.Failure.missing(_) {
            // Read-only compatibility with v1 filenames. A subsequent set writes
            // the unambiguous hashed v2 filename; no implicit data mutation here.
            do { bytes = try files.readFile(legacyReferencePath(name), maximumBytes: limits.maximumMetadataBytes) }
            catch VivoRootedFileStore.Failure.missing(_) { throw VivoArtifactStoreError.referenceMissing(name) }
        }
        if let deleted = try? VivoCanonicalJSON.decode(VivoReferenceDeletion.self, from: bytes) {
            guard deleted.schema == "numivivo.org/reference-deletion/v1",
                  Array(deleted.name.utf8) == Array(name.utf8) else {
                throw VivoArtifactStoreError.invalidDescriptor("invalid reference tombstone")
            }
            throw VivoArtifactStoreError.referenceMissing(name)
        }
        let value = try VivoCanonicalJSON.decode(VivoArtifactReference.self, from: bytes)
        let persisted = try descriptor(for: value.artifact.fingerprint)
        guard Array(value.name.utf8) == Array(name.utf8), value.artifact == persisted,
              value.updatedAt.timeIntervalSince1970.isFinite, try verify(persisted.fingerprint) else {
            throw VivoArtifactStoreError.invalidDescriptor("reference name, persisted descriptor or object integrity mismatch")
        }
        return value
    }
    public func removeReference(_ name: String) throws {
        _ = try reference(name)
        // A v2 tombstone prevents a still-present v1 filename from resurrecting
        // a deleted reference. setReference atomically replaces the tombstone.
        let deletion = VivoReferenceDeletion(schema: "numivivo.org/reference-deletion/v1", name: name, deletedAt: Date())
        try files.writeFile(VivoCanonicalJSON.encode(deletion), relative: referencePath(name), immutable: false)
    }
    private func validateMetadata(kind: String, mediaType: String, attributes: [String: String]) throws {
        guard !kind.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty, kind.utf8.count <= 256,
              !mediaType.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty, mediaType.utf8.count <= 256,
              attributes.count <= 128, attributes.allSatisfy({ $0.key.utf8.count <= 128 && $0.value.utf8.count <= 4096 }) else {
            throw VivoArtifactStoreError.invalidDescriptor("empty or oversized object metadata")
        }
    }
    private func objectRelativePath(_ fingerprint: VivoFingerprint) -> String {
        let hex = fingerprint.hex
        return "objects/sha256/\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))/\(hex)"
    }
    private func descriptorPath(_ fingerprint: VivoFingerprint) -> String {
        let hex = fingerprint.hex
        return "descriptors/sha256/\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))/\(hex).json"
    }
    private func referencePath(_ name: String) throws -> String {
        guard !name.isEmpty, name == name.trimmingCharacters(in: .whitespacesAndNewlines), name.utf8.count <= 240,
              !name.contains("/"), !name.contains("\\"), !name.contains("\0"), name != ".", name != ".." else {
            throw VivoArtifactStoreError.invalidDescriptor("reference name is unsafe or noncanonical")
        }
        // Hash exact UTF-8 bytes: case-insensitive and Unicode-normalizing file
        // systems cannot alias two distinct logical names into one directory entry.
        return try "refs/v2/" + VivoCanonicalJSON.fingerprint(Data(name.utf8)).hex + ".json"
    }
    private func legacyReferencePath(_ name: String) throws -> String {
        _ = try referencePath(name)
        guard let encoded = name.addingPercentEncoding(withAllowedCharacters: .alphanumerics), encoded.utf8.count <= 230 else {
            throw VivoArtifactStoreError.referenceMissing(name)
        }
        return "refs/" + encoded + ".json"
    }
    private func hex(_ value: String, count: Int) -> Bool {
        value.utf8.count == count && value.utf8.allSatisfy { (48...57).contains($0) || (97...102).contains($0) }
    }
    private func fingerprintFromHex(_ value: String) throws -> VivoFingerprint {
        let bytes = Array(value.utf8)
        func digit(_ byte: UInt8) -> UInt8 { byte <= 57 ? byte - 48 : byte - 97 + 10 }
        return try VivoFingerprint(bytes: stride(from: 0, to: bytes.count, by: 2).map { digit(bytes[$0]) * 16 + digit(bytes[$0 + 1]) })
    }
}
