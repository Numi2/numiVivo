import Foundation

public struct VivoStoredArtifact: Codable, Sendable, Equatable, Hashable {
    public let fingerprint: VivoFingerprint
    public let kind: String
    public let mediaType: String
    public let byteCount: UInt64
    public let objectPath: String
    public let createdAt: Date
    public let attributes: [String: String]

    public init(
        fingerprint: VivoFingerprint,
        kind: String,
        mediaType: String,
        byteCount: UInt64,
        objectPath: String,
        createdAt: Date,
        attributes: [String: String]
    ) {
        self.fingerprint = fingerprint
        self.kind = kind
        self.mediaType = mediaType
        self.byteCount = byteCount
        self.objectPath = objectPath
        self.createdAt = createdAt
        self.attributes = attributes
    }
}

public struct VivoArtifactReference: Codable, Sendable, Equatable, Hashable {
    public let name: String
    public let artifact: VivoStoredArtifact
    public let updatedAt: Date

    public init(name: String, artifact: VivoStoredArtifact, updatedAt: Date = Date()) {
        self.name = name
        self.artifact = artifact
        self.updatedAt = updatedAt
    }
}

public enum VivoArtifactStoreError: Error, Sendable, CustomStringConvertible {
    case invalidRoot(String)
    case invalidDescriptor(String)
    case objectMissing(VivoFingerprint)
    case integrityFailure(VivoFingerprint)
    case referenceMissing(String)
    case io(String)

    public var description: String {
        switch self {
        case .invalidRoot(let message): "Invalid artifact store root: \(message)"
        case .invalidDescriptor(let message): "Invalid artifact descriptor: \(message)"
        case .objectMissing(let fingerprint): "Artifact object \(fingerprint.hex) is missing."
        case .integrityFailure(let fingerprint): "Artifact object \(fingerprint.hex) failed integrity verification."
        case .referenceMissing(let name): "Artifact reference '\(name)' does not exist."
        case .io(let message): "Artifact store I/O failed: \(message)"
        }
    }
}

public actor VivoArtifactStore {
    public let rootURL: URL
    private let fileManager: FileManager

    public init(rootURL: URL, createIfNeeded: Bool = true) throws {
        guard rootURL.isFileURL else {
            throw VivoArtifactStoreError.invalidRoot("root must be a file URL")
        }
        self.rootURL = rootURL.standardizedFileURL
        self.fileManager = .default
        if createIfNeeded {
            try Self.createLayout(at: self.rootURL, fileManager: self.fileManager)
        } else {
            var isDirectory: ObjCBool = false
            guard self.fileManager.fileExists(atPath: self.rootURL.path, isDirectory: &isDirectory),
                  isDirectory.boolValue else {
                throw VivoArtifactStoreError.invalidRoot("directory does not exist")
            }
        }
    }

    @discardableResult
    public func put(
        data: Data,
        kind: String,
        mediaType: String,
        attributes: [String: String] = [:]
    ) throws -> VivoStoredArtifact {
        guard !kind.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !mediaType.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoArtifactStoreError.invalidDescriptor("kind and mediaType are required")
        }
        let fingerprint = try VivoCanonicalJSON.fingerprint(data)
        let relativeObjectPath = objectRelativePath(fingerprint)
        let objectURL = rootURL.appendingPathComponent(relativeObjectPath)
        try fileManager.createDirectory(
            at: objectURL.deletingLastPathComponent(),
            withIntermediateDirectories: true
        )

        if fileManager.fileExists(atPath: objectURL.path) {
            let existing = try Data(contentsOf: objectURL, options: [.mappedIfSafe])
            guard try VivoCanonicalJSON.fingerprint(existing) == fingerprint else {
                throw VivoArtifactStoreError.integrityFailure(fingerprint)
            }
        } else {
            try atomicWrite(data, to: objectURL)
        }

        let descriptorURL = descriptorURL(for: fingerprint)
        if fileManager.fileExists(atPath: descriptorURL.path) {
            let descriptor = try decodeDescriptor(at: descriptorURL)
            guard descriptor.fingerprint == fingerprint,
                  descriptor.byteCount == UInt64(data.count),
                  descriptor.objectPath == relativeObjectPath else {
                throw VivoArtifactStoreError.invalidDescriptor("persisted descriptor disagrees with its object")
            }
            return descriptor
        }

        let descriptor = VivoStoredArtifact(
            fingerprint: fingerprint,
            kind: kind,
            mediaType: mediaType,
            byteCount: UInt64(data.count),
            objectPath: relativeObjectPath,
            createdAt: Date(),
            attributes: attributes
        )
        let descriptorData = try VivoCanonicalJSON.encode(descriptor)
        try fileManager.createDirectory(
            at: descriptorURL.deletingLastPathComponent(),
            withIntermediateDirectories: true
        )
        try atomicWrite(descriptorData, to: descriptorURL)
        return descriptor
    }

    @discardableResult
    public func put<Payload: Codable & Sendable>(
        envelope: VivoArtifactEnvelope<Payload>,
        mediaType: String = "application/vnd.numivivo.artifact+json",
        attributes: [String: String] = [:]
    ) throws -> VivoStoredArtifact {
        try put(
            data: envelope.canonicalData(),
            kind: envelope.kind.rawValue,
            mediaType: mediaType,
            attributes: attributes.merging([
                "identifier": envelope.metadata.identifier,
                "version": envelope.metadata.version
            ], uniquingKeysWith: { current, _ in current })
        )
    }

    public func contains(_ fingerprint: VivoFingerprint) -> Bool {
        fileManager.fileExists(atPath: objectURL(for: fingerprint).path)
    }

    public func descriptor(for fingerprint: VivoFingerprint) throws -> VivoStoredArtifact {
        let url = descriptorURL(for: fingerprint)
        guard fileManager.fileExists(atPath: url.path) else {
            throw VivoArtifactStoreError.objectMissing(fingerprint)
        }
        return try decodeDescriptor(at: url)
    }

    public func data(for fingerprint: VivoFingerprint, verify: Bool = true) throws -> Data {
        let url = objectURL(for: fingerprint)
        guard fileManager.fileExists(atPath: url.path) else {
            throw VivoArtifactStoreError.objectMissing(fingerprint)
        }
        let data = try Data(contentsOf: url, options: [.mappedIfSafe])
        if verify, try VivoCanonicalJSON.fingerprint(data) != fingerprint {
            throw VivoArtifactStoreError.integrityFailure(fingerprint)
        }
        return data
    }

    public func verify(_ fingerprint: VivoFingerprint) throws -> Bool {
        let descriptor = try descriptor(for: fingerprint)
        let data = try self.data(for: fingerprint, verify: true)
        return descriptor.byteCount == UInt64(data.count) &&
               descriptor.objectPath == objectRelativePath(fingerprint)
    }

    public func list(kind: String? = nil) throws -> [VivoStoredArtifact] {
        let directory = rootURL.appendingPathComponent("descriptors/sha256", isDirectory: true)
        guard fileManager.fileExists(atPath: directory.path) else { return [] }
        guard let enumerator = fileManager.enumerator(
            at: directory,
            includingPropertiesForKeys: [.isRegularFileKey],
            options: [.skipsHiddenFiles]
        ) else { return [] }
        var artifacts: [VivoStoredArtifact] = []
        for case let url as URL in enumerator where url.pathExtension == "json" {
            let value = try decodeDescriptor(at: url)
            if kind == nil || value.kind == kind { artifacts.append(value) }
        }
        return artifacts.sorted {
            if $0.createdAt != $1.createdAt { return $0.createdAt < $1.createdAt }
            return $0.fingerprint.hex < $1.fingerprint.hex
        }
    }

    public func materialize(
        _ fingerprint: VivoFingerprint,
        to destination: URL,
        verify: Bool = true
    ) throws {
        guard destination.isFileURL else {
            throw VivoArtifactStoreError.io("materialization destination must be a file URL")
        }
        let data = try self.data(for: fingerprint, verify: verify)
        try fileManager.createDirectory(
            at: destination.deletingLastPathComponent(),
            withIntermediateDirectories: true
        )
        try atomicWrite(data, to: destination)
    }

    @discardableResult
    public func setReference(
        _ name: String,
        to artifact: VivoStoredArtifact
    ) throws -> VivoArtifactReference {
        let component = try safeReferenceComponent(name)
        guard try verify(artifact.fingerprint) else {
            throw VivoArtifactStoreError.integrityFailure(artifact.fingerprint)
        }
        let reference = VivoArtifactReference(name: name, artifact: artifact)
        let data = try VivoCanonicalJSON.encode(reference)
        let url = rootURL.appendingPathComponent("refs/\(component).json")
        try fileManager.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        try atomicWrite(data, to: url)
        return reference
    }

    public func reference(_ name: String) throws -> VivoArtifactReference {
        let component = try safeReferenceComponent(name)
        let url = rootURL.appendingPathComponent("refs/\(component).json")
        guard fileManager.fileExists(atPath: url.path) else {
            throw VivoArtifactStoreError.referenceMissing(name)
        }
        do {
            let reference = try VivoCanonicalJSON.decode(
                VivoArtifactReference.self,
                from: Data(contentsOf: url)
            )
            guard reference.name == name else {
                throw VivoArtifactStoreError.invalidDescriptor("reference name does not match its path")
            }
            return reference
        } catch let error as VivoArtifactStoreError {
            throw error
        } catch {
            throw VivoArtifactStoreError.io(String(describing: error))
        }
    }

    public func removeReference(_ name: String) throws {
        let component = try safeReferenceComponent(name)
        let url = rootURL.appendingPathComponent("refs/\(component).json")
        guard fileManager.fileExists(atPath: url.path) else {
            throw VivoArtifactStoreError.referenceMissing(name)
        }
        try fileManager.removeItem(at: url)
    }

    private static func createLayout(at root: URL, fileManager: FileManager) throws {
        do {
            try fileManager.createDirectory(at: root, withIntermediateDirectories: true)
            for path in ["objects/sha256", "descriptors/sha256", "refs", "bundles"] {
                try fileManager.createDirectory(
                    at: root.appendingPathComponent(path, isDirectory: true),
                    withIntermediateDirectories: true
                )
            }
        } catch {
            throw VivoArtifactStoreError.io(String(describing: error))
        }
    }

    private func objectRelativePath(_ fingerprint: VivoFingerprint) -> String {
        let hex = fingerprint.hex
        return "objects/sha256/\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))/\(hex)"
    }

    private func objectURL(for fingerprint: VivoFingerprint) -> URL {
        rootURL.appendingPathComponent(objectRelativePath(fingerprint))
    }

    private func descriptorURL(for fingerprint: VivoFingerprint) -> URL {
        let hex = fingerprint.hex
        return rootURL.appendingPathComponent(
            "descriptors/sha256/\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))/\(hex).json"
        )
    }

    private func decodeDescriptor(at url: URL) throws -> VivoStoredArtifact {
        do {
            return try VivoCanonicalJSON.decode(
                VivoStoredArtifact.self,
                from: Data(contentsOf: url, options: [.mappedIfSafe])
            )
        } catch {
            throw VivoArtifactStoreError.invalidDescriptor(String(describing: error))
        }
    }

    private func safeReferenceComponent(_ name: String) throws -> String {
        let trimmed = name.trimmingCharacters(in: .whitespacesAndNewlines)
        guard !trimmed.isEmpty,
              trimmed.utf8.count <= 240,
              !trimmed.contains("/"),
              !trimmed.contains("\\"),
              trimmed != ".",
              trimmed != ".." else {
            throw VivoArtifactStoreError.invalidDescriptor("reference name is unsafe")
        }
        return trimmed.addingPercentEncoding(withAllowedCharacters: .alphanumerics) ?? UUID().uuidString
    }

    private func atomicWrite(_ data: Data, to destination: URL) throws {
        let temporary = destination.deletingLastPathComponent().appendingPathComponent(
            ".\(destination.lastPathComponent).\(UUID().uuidString).tmp"
        )
        do {
            try data.write(to: temporary, options: [.atomic])
            if fileManager.fileExists(atPath: destination.path) {
                _ = try fileManager.replaceItemAt(destination, withItemAt: temporary)
            } else {
                try fileManager.moveItem(at: temporary, to: destination)
            }
        } catch {
            try? fileManager.removeItem(at: temporary)
            throw VivoArtifactStoreError.io(String(describing: error))
        }
    }
}
