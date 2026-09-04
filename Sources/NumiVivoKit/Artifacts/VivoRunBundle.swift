import Foundation

public struct VivoRunBundleEntry: Codable, Sendable, Equatable, Hashable {
    public let path: String
    public let role: String
    public let mediaType: String
    public let byteCount: UInt64
    public let fingerprint: VivoFingerprint
    public let sourceArtifacts: [VivoFingerprint]

    public init(
        path: String,
        role: String,
        mediaType: String,
        byteCount: UInt64,
        fingerprint: VivoFingerprint,
        sourceArtifacts: [VivoFingerprint] = []
    ) {
        self.path = path
        self.role = role
        self.mediaType = mediaType
        self.byteCount = byteCount
        self.fingerprint = fingerprint
        self.sourceArtifacts = sourceArtifacts
    }
}

public struct VivoRunBundleManifest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/run-bundle/v1"

    public let schema: String
    public let runIdentifier: String
    public let createdAt: Date
    public let softwareVersion: String
    public let programFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint?
    public let experimentFingerprint: VivoFingerprint?
    public let couplingFingerprint: VivoFingerprint?
    public let entries: [VivoRunBundleEntry]
    public let labels: [String: String]

    public init(
        runIdentifier: String,
        createdAt: Date = Date(),
        softwareVersion: String,
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint? = nil,
        experimentFingerprint: VivoFingerprint? = nil,
        couplingFingerprint: VivoFingerprint? = nil,
        entries: [VivoRunBundleEntry],
        labels: [String: String] = [:]
    ) {
        self.schema = Self.schema
        self.runIdentifier = runIdentifier
        self.createdAt = createdAt
        self.softwareVersion = softwareVersion
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.entries = entries
        self.labels = labels
    }
}

public struct VivoRunBundleSource: Sendable {
    public let path: String
    public let role: String
    public let mediaType: String
    public let data: Data
    public let sourceArtifacts: [VivoFingerprint]

    public init(
        path: String,
        role: String,
        mediaType: String,
        data: Data,
        sourceArtifacts: [VivoFingerprint] = []
    ) {
        self.path = path
        self.role = role
        self.mediaType = mediaType
        self.data = data
        self.sourceArtifacts = sourceArtifacts
    }
}

public struct VivoRunBundlePlan: Sendable {
    public let manifest: VivoRunBundleManifest
    public let manifestData: Data
    public let files: [String: Data]
    public let fingerprint: VivoFingerprint

    public init(
        manifest: VivoRunBundleManifest,
        manifestData: Data,
        files: [String: Data],
        fingerprint: VivoFingerprint
    ) {
        self.manifest = manifest
        self.manifestData = manifestData
        self.files = files
        self.fingerprint = fingerprint
    }

    public func write(to destination: URL) throws {
        guard destination.isFileURL else {
            throw VivoArtifactValidationError.invalid("run-bundle destination must be a file URL")
        }
        let manager = FileManager.default
        let parent = destination.deletingLastPathComponent()
        try manager.createDirectory(at: parent, withIntermediateDirectories: true)
        let staging = parent.appendingPathComponent(
            ".\(destination.lastPathComponent).\(UUID().uuidString).staging",
            isDirectory: true
        )
        try manager.createDirectory(at: staging, withIntermediateDirectories: true)
        do {
            for (path, data) in files.sorted(by: { $0.key < $1.key }) {
                let target = staging.appendingPathComponent(path)
                try manager.createDirectory(
                    at: target.deletingLastPathComponent(),
                    withIntermediateDirectories: true
                )
                try data.write(to: target, options: [.atomic])
            }
            let completion = try VivoCanonicalJSON.encode([
                "schema": "numivivo.org/run-bundle-completion/v1",
                "bundleFingerprint": fingerprint.hex,
                "manifestFingerprint": try VivoCanonicalJSON.fingerprint(manifestData).hex
            ])
            try completion.write(
                to: staging.appendingPathComponent("COMPLETE.json"),
                options: [.atomic]
            )

            if manager.fileExists(atPath: destination.path) {
                let backup = parent.appendingPathComponent(
                    ".\(destination.lastPathComponent).\(UUID().uuidString).backup"
                )
                try manager.moveItem(at: destination, to: backup)
                do {
                    try manager.moveItem(at: staging, to: destination)
                    try? manager.removeItem(at: backup)
                } catch {
                    try? manager.moveItem(at: backup, to: destination)
                    throw error
                }
            } else {
                try manager.moveItem(at: staging, to: destination)
            }
        } catch {
            try? manager.removeItem(at: staging)
            throw error
        }
    }
}

public struct VivoRunBundleBuilder: Sendable {
    public init() {}

    public func build(
        runIdentifier: String,
        softwareVersion: String,
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint? = nil,
        experimentFingerprint: VivoFingerprint? = nil,
        couplingFingerprint: VivoFingerprint? = nil,
        sources: [VivoRunBundleSource],
        labels: [String: String] = [:],
        createdAt: Date = Date()
    ) throws -> VivoRunBundlePlan {
        guard !runIdentifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !softwareVersion.isEmpty else {
            throw VivoArtifactValidationError.invalid("runIdentifier and softwareVersion are required")
        }
        var paths = Set<String>()
        var files: [String: Data] = [:]
        var entries: [VivoRunBundleEntry] = []
        entries.reserveCapacity(sources.count)

        for source in sources {
            try validate(path: source.path)
            guard paths.insert(source.path).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate run-bundle path \(source.path)")
            }
            guard !source.role.isEmpty, !source.mediaType.isEmpty else {
                throw VivoArtifactValidationError.invalid("run-bundle source role and mediaType are required")
            }
            let fingerprint = try VivoCanonicalJSON.fingerprint(source.data)
            files[source.path] = source.data
            entries.append(.init(
                path: source.path,
                role: source.role,
                mediaType: source.mediaType,
                byteCount: UInt64(source.data.count),
                fingerprint: fingerprint,
                sourceArtifacts: source.sourceArtifacts
            ))
        }

        let manifest = VivoRunBundleManifest(
            runIdentifier: runIdentifier,
            createdAt: createdAt,
            softwareVersion: softwareVersion,
            programFingerprint: programFingerprint,
            hostContextFingerprint: hostContextFingerprint,
            experimentFingerprint: experimentFingerprint,
            couplingFingerprint: couplingFingerprint,
            entries: entries.sorted { $0.path < $1.path },
            labels: labels
        )
        let manifestData = try VivoCanonicalJSON.encode(manifest)
        files["manifest.json"] = manifestData

        struct BundleDigest: Codable {
            let manifestFingerprint: VivoFingerprint
            let entries: [VivoRunBundleEntry]
        }
        let fingerprint = try VivoCanonicalJSON.fingerprint(
            VivoCanonicalJSON.encode(BundleDigest(
                manifestFingerprint: try VivoCanonicalJSON.fingerprint(manifestData),
                entries: manifest.entries
            ))
        )
        return VivoRunBundlePlan(
            manifest: manifest,
            manifestData: manifestData,
            files: files,
            fingerprint: fingerprint
        )
    }

    private func validate(path: String) throws {
        let components = path.split(separator: "/", omittingEmptySubsequences: false)
        guard !path.isEmpty,
              !path.hasPrefix("/"),
              !path.hasSuffix("/"),
              !path.contains("\\"),
              !components.contains(""),
              !components.contains("."),
              !components.contains(".."),
              path != "manifest.json",
              path != "COMPLETE.json" else {
            throw VivoArtifactValidationError.invalid("unsafe or reserved run-bundle path \(path)")
        }
    }
}
