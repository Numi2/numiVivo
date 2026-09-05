import Foundation

public enum VivoEvidenceDataFormat: String, Codable, Sendable, CaseIterable {
    case csv
    case tsv
    case json
    case jsonLines
    case arrow
    case parquet
    case hdf5
    case binary
    case other
}

public struct VivoEvidenceTableDescriptor: Codable, Sendable, Equatable, Hashable {
    public var format: VivoEvidenceDataFormat
    public var hasHeader: Bool?
    public var columns: [String]
    public var semanticBindings: [String: String]
    public var rowCount: UInt64?
    public var notes: [String: String]

    public init(
        format: VivoEvidenceDataFormat,
        hasHeader: Bool? = nil,
        columns: [String] = [],
        semanticBindings: [String: String] = [:],
        rowCount: UInt64? = nil,
        notes: [String: String] = [:]
    ) {
        self.format = format
        self.hasHeader = hasHeader
        self.columns = columns
        self.semanticBindings = semanticBindings
        self.rowCount = rowCount
        self.notes = notes
    }

    public func validate(label: String) throws {
        guard columns.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 512 }),
              Set(columns).count == columns.count,
              semanticBindings.keys.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }),
              semanticBindings.values.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 512 }),
              rowCount.map({ $0 <= 10_000_000_000_000 }) ?? true else {
            throw VivoArtifactValidationError.invalid(
                "\(label) table descriptor is invalid"
            )
        }
        if !columns.isEmpty {
            let columnSet = Set(columns)
            guard semanticBindings.values.allSatisfy(columnSet.contains) else {
                throw VivoArtifactValidationError.invalid(
                    "\(label) semantic binding references an undeclared column"
                )
            }
        }
    }
}

public struct VivoEvidenceSnapshotEntry: Codable, Sendable, Equatable, Hashable {
    public let identifier: String
    public let classification: VivoEvidenceClass
    public let objectFingerprint: VivoFingerprint
    public let byteCount: UInt64
    public let mediaType: String
    public let originalFileName: String?
    public let sourceURI: String?
    public let datasetIdentifier: String?
    public let citation: String?
    public let license: String?
    public let table: VivoEvidenceTableDescriptor?
    public let attributes: [String: String]

    public init(
        identifier: String,
        classification: VivoEvidenceClass,
        objectFingerprint: VivoFingerprint,
        byteCount: UInt64,
        mediaType: String,
        originalFileName: String? = nil,
        sourceURI: String? = nil,
        datasetIdentifier: String? = nil,
        citation: String? = nil,
        license: String? = nil,
        table: VivoEvidenceTableDescriptor? = nil,
        attributes: [String: String] = [:]
    ) {
        self.identifier = identifier
        self.classification = classification
        self.objectFingerprint = objectFingerprint
        self.byteCount = byteCount
        self.mediaType = mediaType
        self.originalFileName = originalFileName
        self.sourceURI = sourceURI
        self.datasetIdentifier = datasetIdentifier
        self.citation = citation
        self.license = license
        self.table = table
        self.attributes = attributes
    }

    public func validate(label: String) throws {
        guard !identifier.isEmpty,
              identifier.utf8.count <= 256,
              byteCount > 0,
              !mediaType.isEmpty,
              mediaType.utf8.count <= 256,
              originalFileName.map({ !$0.isEmpty && $0.utf8.count <= 1_024 }) ?? true,
              sourceURI.map({ !$0.isEmpty && $0.utf8.count <= 8_192 }) ?? true,
              datasetIdentifier.map({ !$0.isEmpty && $0.utf8.count <= 1_024 }) ?? true,
              citation.map({ !$0.isEmpty && $0.utf8.count <= 8_192 }) ?? true,
              license.map({ !$0.isEmpty && $0.utf8.count <= 1_024 }) ?? true,
              attributes.count <= 4_096,
              attributes.keys.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }),
              attributes.values.allSatisfy({ $0.utf8.count <= 8_192 }) else {
            throw VivoArtifactValidationError.invalid(
                "\(label) evidence entry metadata is invalid"
            )
        }
        try table?.validate(label: "\(label).table")
    }
}

public struct VivoEvidenceSnapshot: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/evidence-snapshot/v1"

    public let schema: String
    public let identifier: String
    public let version: String
    public let title: String
    public let entries: [VivoEvidenceSnapshotEntry]
    public let labels: [String: String]
    public let fingerprint: VivoFingerprint

    public init(
        identifier: String,
        version: String,
        title: String,
        entries: [VivoEvidenceSnapshotEntry],
        labels: [String: String] = [:]
    ) throws {
        let orderedEntries = entries.sorted { $0.identifier < $1.identifier }
        let unsigned = VivoEvidenceSnapshot(
            schema: Self.schema,
            identifier: identifier,
            version: version,
            title: title,
            entries: orderedEntries,
            labels: labels,
            fingerprint: try VivoFingerprint(bytes: [UInt8](repeating: 0, count: 32))
        )
        try unsigned.validate(ignoreFingerprint: true)
        let digest = try VivoCanonicalJSON.fingerprint(
            VivoCanonicalJSON.encode(VivoEvidenceSnapshotIdentity(snapshot: unsigned))
        )
        self = VivoEvidenceSnapshot(
            schema: Self.schema,
            identifier: identifier,
            version: version,
            title: title,
            entries: orderedEntries,
            labels: labels,
            fingerprint: digest
        )
    }

    private init(
        schema: String,
        identifier: String,
        version: String,
        title: String,
        entries: [VivoEvidenceSnapshotEntry],
        labels: [String: String],
        fingerprint: VivoFingerprint
    ) {
        self.schema = schema
        self.identifier = identifier
        self.version = version
        self.title = title
        self.entries = entries
        self.labels = labels
        self.fingerprint = fingerprint
    }

    public func validate() throws {
        try validate(ignoreFingerprint: false)
    }

    public func entry(_ identifier: String) -> VivoEvidenceSnapshotEntry? {
        entries.first(where: { $0.identifier == identifier })
    }

    public func reference(
        entry identifier: String,
        context: String? = nil,
        note: String? = nil
    ) throws -> VivoEvidenceReference {
        guard let entry = entry(identifier) else {
            throw VivoArtifactValidationError.unresolved(
                "evidence snapshot entry \(identifier) does not exist"
            )
        }
        return VivoEvidenceReference(
            classification: entry.classification,
            sourceURI: entry.sourceURI,
            datasetIdentifier: entry.datasetIdentifier ?? entry.identifier,
            citation: entry.citation,
            context: context,
            note: note,
            snapshotFingerprint: fingerprint,
            snapshotEntryIdentifier: entry.identifier
        )
    }

    private func validate(ignoreFingerprint: Bool) throws {
        guard schema == Self.schema,
              !identifier.isEmpty,
              identifier.utf8.count <= 256,
              !version.isEmpty,
              version.utf8.count <= 128,
              !title.isEmpty,
              title.utf8.count <= 1_024,
              !entries.isEmpty,
              entries.count <= 1_000_000,
              labels.count <= 4_096,
              labels.keys.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }),
              labels.values.allSatisfy({ $0.utf8.count <= 8_192 }) else {
            throw VivoArtifactValidationError.invalid(
                "evidence snapshot metadata or entry count is invalid"
            )
        }
        guard Set(entries.map(\.identifier)).count == entries.count else {
            throw VivoArtifactValidationError.invalid(
                "evidence snapshot entry identifiers must be unique"
            )
        }
        for (index, entry) in entries.enumerated() {
            try entry.validate(label: "evidence.entries[\(index)]")
        }
        if !ignoreFingerprint {
            let expected = try VivoCanonicalJSON.fingerprint(
                VivoCanonicalJSON.encode(VivoEvidenceSnapshotIdentity(snapshot: self))
            )
            guard expected == fingerprint else {
                throw VivoArtifactValidationError.invalid(
                    "evidence snapshot fingerprint does not match its canonical identity"
                )
            }
        }
    }
}

private struct VivoEvidenceSnapshotIdentity: Codable {
    let schema: String
    let identifier: String
    let version: String
    let title: String
    let entries: [VivoEvidenceSnapshotEntry]
    let labels: [String: String]

    init(snapshot: VivoEvidenceSnapshot) {
        schema = snapshot.schema
        identifier = snapshot.identifier
        version = snapshot.version
        title = snapshot.title
        entries = snapshot.entries
        labels = snapshot.labels
    }
}

public struct VivoEvidenceSnapshotInput: Sendable {
    public var identifier: String
    public var classification: VivoEvidenceClass
    public var data: Data
    public var mediaType: String
    public var originalFileName: String?
    public var sourceURI: String?
    public var datasetIdentifier: String?
    public var citation: String?
    public var license: String?
    public var table: VivoEvidenceTableDescriptor?
    public var attributes: [String: String]

    public init(
        identifier: String,
        classification: VivoEvidenceClass,
        data: Data,
        mediaType: String,
        originalFileName: String? = nil,
        sourceURI: String? = nil,
        datasetIdentifier: String? = nil,
        citation: String? = nil,
        license: String? = nil,
        table: VivoEvidenceTableDescriptor? = nil,
        attributes: [String: String] = [:]
    ) {
        self.identifier = identifier
        self.classification = classification
        self.data = data
        self.mediaType = mediaType
        self.originalFileName = originalFileName
        self.sourceURI = sourceURI
        self.datasetIdentifier = datasetIdentifier
        self.citation = citation
        self.license = license
        self.table = table
        self.attributes = attributes
    }
}

public actor VivoEvidenceSnapshotBuilder {
    public struct Limits: Sendable, Equatable {
        public var maximumEntries: Int
        public var maximumBytesPerEntry: Int
        public var maximumTotalBytes: UInt64

        public init(
            maximumEntries: Int = 100_000,
            maximumBytesPerEntry: Int = 1_073_741_824,
            maximumTotalBytes: UInt64 = 8_589_934_592
        ) {
            self.maximumEntries = maximumEntries
            self.maximumBytesPerEntry = maximumBytesPerEntry
            self.maximumTotalBytes = maximumTotalBytes
        }
    }

    private let store: VivoArtifactStore
    private let limits: Limits

    public init(
        store: VivoArtifactStore,
        limits: Limits = .init()
    ) {
        self.store = store
        self.limits = limits
    }

    @discardableResult
    public func build(
        identifier: String,
        version: String,
        title: String,
        inputs: [VivoEvidenceSnapshotInput],
        labels: [String: String] = [:]
    ) async throws -> (snapshot: VivoEvidenceSnapshot, stored: VivoStoredArtifact) {
        guard limits.maximumEntries > 0,
              limits.maximumBytesPerEntry > 0,
              limits.maximumTotalBytes > 0,
              !inputs.isEmpty,
              inputs.count <= limits.maximumEntries,
              Set(inputs.map(\.identifier)).count == inputs.count else {
            throw VivoArtifactValidationError.invalid(
                "evidence snapshot inputs or builder limits are invalid"
            )
        }

        var totalBytes: UInt64 = 0
        var entries: [VivoEvidenceSnapshotEntry] = []
        entries.reserveCapacity(inputs.count)
        for input in inputs.sorted(by: { $0.identifier < $1.identifier }) {
            guard !input.data.isEmpty,
                  input.data.count <= limits.maximumBytesPerEntry else {
                throw VivoArtifactValidationError.invalid(
                    "evidence input \(input.identifier) is empty or exceeds its byte limit"
                )
            }
            let addition = totalBytes.addingReportingOverflow(UInt64(input.data.count))
            guard !addition.overflow,
                  addition.partialValue <= limits.maximumTotalBytes else {
                throw VivoArtifactValidationError.invalid(
                    "evidence snapshot exceeds the configured total byte limit"
                )
            }
            totalBytes = addition.partialValue

            let storedObject = try await store.put(
                data: input.data,
                kind: "evidence-object",
                mediaType: input.mediaType,
                attributes: input.attributes.merging([
                    "evidenceIdentifier": input.identifier,
                    "classification": input.classification.rawValue
                ], uniquingKeysWith: { current, _ in current })
            )
            let entry = VivoEvidenceSnapshotEntry(
                identifier: input.identifier,
                classification: input.classification,
                objectFingerprint: storedObject.fingerprint,
                byteCount: UInt64(input.data.count),
                mediaType: input.mediaType,
                originalFileName: input.originalFileName,
                sourceURI: input.sourceURI,
                datasetIdentifier: input.datasetIdentifier,
                citation: input.citation,
                license: input.license,
                table: input.table,
                attributes: input.attributes
            )
            try entry.validate(label: "evidence.\(input.identifier)")
            entries.append(entry)
        }

        let snapshot = try VivoEvidenceSnapshot(
            identifier: identifier,
            version: version,
            title: title,
            entries: entries,
            labels: labels
        )
        let snapshotData = try VivoCanonicalJSON.encode(snapshot)
        let storedSnapshot = try await store.put(
            data: snapshotData,
            kind: VivoArtifactKind.evidenceSnapshot.rawValue,
            mediaType: "application/vnd.numivivo.evidence-snapshot+json",
            attributes: [
                "identifier": identifier,
                "version": version,
                "entryCount": String(entries.count),
                "payloadBytes": String(totalBytes)
            ]
        )
        guard storedSnapshot.fingerprint == (try VivoCanonicalJSON.fingerprint(snapshotData)) else {
            throw VivoArtifactStoreError.integrityFailure(storedSnapshot.fingerprint)
        }
        return (snapshot, storedSnapshot)
    }

    public func verify(
        _ snapshot: VivoEvidenceSnapshot,
        verifyObjects: Bool = true
    ) async throws -> Bool {
        try snapshot.validate()
        guard verifyObjects else { return true }
        for entry in snapshot.entries {
            let descriptor = try await store.descriptor(for: entry.objectFingerprint)
            guard descriptor.fingerprint == entry.objectFingerprint,
                  descriptor.byteCount == entry.byteCount,
                  descriptor.mediaType == entry.mediaType,
                  try await store.verify(entry.objectFingerprint) else {
                return false
            }
        }
        return true
    }
}
