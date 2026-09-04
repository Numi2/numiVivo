import CryptoKit
import Foundation

public struct VivoArtifactLoadLimits: Sendable, Equatable {
    public var maximumBytes: Int
    public var maximumDepth: Int
    public var maximumNodes: Int
    public var maximumStringBytes: Int
    public var maximumArrayElements: Int
    public var maximumObjectMembers: Int
    public var requireRegularFile: Bool
    public var rejectSymbolicLinks: Bool

    public init(
        maximumBytes: Int = 64 * 1_024 * 1_024,
        maximumDepth: Int = 128,
        maximumNodes: Int = 1_000_000,
        maximumStringBytes: Int = 4 * 1_024 * 1_024,
        maximumArrayElements: Int = 1_000_000,
        maximumObjectMembers: Int = 250_000,
        requireRegularFile: Bool = true,
        rejectSymbolicLinks: Bool = true
    ) {
        self.maximumBytes = maximumBytes
        self.maximumDepth = maximumDepth
        self.maximumNodes = maximumNodes
        self.maximumStringBytes = maximumStringBytes
        self.maximumArrayElements = maximumArrayElements
        self.maximumObjectMembers = maximumObjectMembers
        self.requireRegularFile = requireRegularFile
        self.rejectSymbolicLinks = rejectSymbolicLinks
    }

    public func validate() throws {
        guard maximumBytes > 0,
              maximumDepth > 0,
              maximumNodes > 0,
              maximumStringBytes > 0,
              maximumArrayElements > 0,
              maximumObjectMembers > 0 else {
            throw VivoArtifactLoadError.invalidLimits
        }
    }
}

public struct VivoLoadedJSONArtifact<Value: Sendable>: Sendable {
    public let value: Value
    public let sourceFingerprint: String
    public let sourceBytes: Int

    public init(value: Value, sourceFingerprint: String, sourceBytes: Int) {
        self.value = value
        self.sourceFingerprint = sourceFingerprint
        self.sourceBytes = sourceBytes
    }
}

public enum VivoArtifactLoadError: Error, LocalizedError, Sendable {
    case invalidLimits
    case fileMetadata(String)
    case fileTooLarge(actual: Int, maximum: Int)
    case nonRegularFile
    case symbolicLinkRejected
    case invalidJSON(String)
    case structureLimit(String)
    case decoding(String)
    case unsupportedSchema(UInt32)
    case fingerprintMismatch(subject: String, expected: String, actual: String)
    case invalidFingerprint(String)
    case invalidCampaign(String)
    case campaignLedgerMismatch(job: UInt64, reason: String)
    case typeValidation(String)

    public var errorDescription: String? {
        switch self {
        case .invalidLimits:
            return "artifact load limits are invalid"
        case .fileMetadata(let reason):
            return "artifact file metadata is unavailable: \(reason)"
        case .fileTooLarge(let actual, let maximum):
            return "artifact has \(actual) bytes; the limit is \(maximum)"
        case .nonRegularFile:
            return "artifact URL does not reference a regular file"
        case .symbolicLinkRejected:
            return "artifact symbolic links are rejected by the load policy"
        case .invalidJSON(let reason):
            return "artifact JSON is invalid: \(reason)"
        case .structureLimit(let reason):
            return "artifact JSON exceeds a structural limit: \(reason)"
        case .decoding(let reason):
            return "artifact decoding failed: \(reason)"
        case .unsupportedSchema(let version):
            return "artifact schema version \(version) is unsupported"
        case .fingerprintMismatch(let subject, let expected, let actual):
            return "artifact fingerprint mismatch for \(subject): expected \(expected), found \(actual)"
        case .invalidFingerprint(let subject):
            return "artifact contains an invalid SHA-256 fingerprint for \(subject)"
        case .invalidCampaign(let reason):
            return "campaign manifest is invalid: \(reason)"
        case .campaignLedgerMismatch(let job, let reason):
            return "campaign ledger mismatch at job \(job): \(reason)"
        case .typeValidation(let reason):
            return "artifact semantic validation failed: \(reason)"
        }
    }
}

public enum VivoValidatedArtifactLoader {
    public static func data(
        at url: URL,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> Data {
        try limits.validate()
        let values: URLResourceValues
        do {
            values = try url.resourceValues(forKeys: [
                .fileSizeKey,
                .isRegularFileKey,
                .isSymbolicLinkKey
            ])
        } catch {
            throw VivoArtifactLoadError.fileMetadata(error.localizedDescription)
        }
        if limits.rejectSymbolicLinks, values.isSymbolicLink == true {
            throw VivoArtifactLoadError.symbolicLinkRejected
        }
        if limits.requireRegularFile, values.isRegularFile != true {
            throw VivoArtifactLoadError.nonRegularFile
        }
        if let size = values.fileSize, size > limits.maximumBytes {
            throw VivoArtifactLoadError.fileTooLarge(
                actual: size,
                maximum: limits.maximumBytes
            )
        }
        let data: Data
        do {
            data = try Data(contentsOf: url, options: [.mappedIfSafe])
        } catch {
            throw VivoArtifactLoadError.fileMetadata(error.localizedDescription)
        }
        guard data.count <= limits.maximumBytes else {
            throw VivoArtifactLoadError.fileTooLarge(
                actual: data.count,
                maximum: limits.maximumBytes
            )
        }
        return data
    }

    public static func decode<T: Decodable & Sendable>(
        _ type: T.Type,
        from data: Data,
        limits: VivoArtifactLoadLimits = .init(),
        validate: @Sendable (T) throws -> Void = { _ in }
    ) throws -> VivoLoadedJSONArtifact<T> {
        try limits.validate()
        guard data.count <= limits.maximumBytes else {
            throw VivoArtifactLoadError.fileTooLarge(
                actual: data.count,
                maximum: limits.maximumBytes
            )
        }
        try validateStructure(data, limits: limits)

        let value: T
        do {
            value = try JSONDecoder().decode(T.self, from: data)
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }
        do {
            try validate(value)
        } catch let error as VivoArtifactLoadError {
            throw error
        } catch {
            throw VivoArtifactLoadError.typeValidation(error.localizedDescription)
        }
        return .init(
            value: value,
            sourceFingerprint: hex(SHA256.hash(data: data)),
            sourceBytes: data.count
        )
    }

    public static func decode<T: Decodable & Sendable>(
        _ type: T.Type,
        at url: URL,
        limits: VivoArtifactLoadLimits = .init(),
        validate: @Sendable (T) throws -> Void = { _ in }
    ) throws -> VivoLoadedJSONArtifact<T> {
        try decode(
            type,
            from: data(at: url, limits: limits),
            limits: limits,
            validate: validate
        )
    }

    public static func partitionModel(
        from data: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoLoadedJSONArtifact<VivoPhysiologicalPartitionModel> {
        try decode(VivoPhysiologicalPartitionModel.self, from: data, limits: limits) { model in
            try model.validate()
            let rebuilt = try VivoPhysiologicalPartitionModel(
                name: model.name,
                compartments: model.compartments,
                analytes: model.analytes,
                edges: model.edges
            )
            try verifyOptionalGeneratedFingerprint(
                subject: "physiological partition model",
                supplied: model.fingerprint,
                generated: rebuilt.fingerprint
            )
        }
    }

    public static func populationModel(
        from data: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoLoadedJSONArtifact<VivoPopulationModel> {
        try decode(VivoPopulationModel.self, from: data, limits: limits) { model in
            try model.validate()
            if !model.fingerprint.isEmpty {
                guard isFingerprint(model.fingerprint) else {
                    throw VivoArtifactLoadError.invalidFingerprint("population model")
                }
            }
        }
    }

    public static func surrogateContract(
        from data: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoLoadedJSONArtifact<VivoSurrogateContract> {
        try decode(VivoSurrogateContract.self, from: data, limits: limits) { contract in
            try contract.validate()
            let rebuilt = try VivoSurrogateContract(
                id: contract.id,
                modelFingerprint: contract.modelFingerprint,
                trainingDataFingerprint: contract.trainingDataFingerprint,
                mechanismPackFingerprint: contract.mechanismPackFingerprint,
                hostContextFingerprint: contract.hostContextFingerprint,
                inputs: contract.inputs,
                outputs: contract.outputs,
                uncertaintyKind: contract.uncertaintyKind,
                maximumNormalizedUncertainty: contract.maximumNormalizedUncertainty,
                maximumNormalizedExtrapolation: contract.maximumNormalizedExtrapolation,
                maximumConsecutiveAcceptedSteps: contract.maximumConsecutiveAcceptedSteps,
                mandatoryAuthorityInterval: contract.mandatoryAuthorityInterval
            )
            try verifyOptionalGeneratedFingerprint(
                subject: "surrogate contract",
                supplied: contract.fingerprint,
                generated: rebuilt.fingerprint
            )
        }
    }

    public static func campaignDefinition(
        from data: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoLoadedJSONArtifact<VivoCampaignDefinition> {
        try decode(VivoCampaignDefinition.self, from: data, limits: limits) { definition in
            try definition.validate()
        }
    }

    public static func campaignManifest(
        from data: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoLoadedJSONArtifact<VivoCampaignManifest> {
        try decode(VivoCampaignManifest.self, from: data, limits: limits) { manifest in
            try validateCampaignManifest(manifest)
        }
    }

    public static func checkpoint(
        from data: Data,
        verificationKey: Curve25519.Signing.PublicKey? = nil,
        requireSignature: Bool = false,
        limits: VivoCheckpointLimits = .init()
    ) throws -> VivoDecodedCheckpoint {
        try VivoCheckpointCodec.decode(
            data,
            verificationKey: verificationKey,
            requireSignature: requireSignature,
            limits: limits
        )
    }

    public static func validateCampaignManifest(_ manifest: VivoCampaignManifest) throws {
        guard manifest.schemaVersion == 1 else {
            throw VivoArtifactLoadError.unsupportedSchema(manifest.schemaVersion)
        }
        guard isFingerprint(manifest.definitionFingerprint),
              isFingerprint(manifest.ledgerHead),
              isFingerprint(manifest.fingerprint) else {
            throw VivoArtifactLoadError.invalidFingerprint("campaign manifest")
        }
        try manifest.artifacts.validate()
        guard !manifest.candidates.isEmpty,
              !manifest.jobs.isEmpty,
              manifest.candidates.count <= Int(UInt32.max),
              manifest.candidates.map(\.index) == Array(0..<UInt32(manifest.candidates.count)),
              Set(manifest.candidates.map(\.fingerprint)).count == manifest.candidates.count,
              manifest.jobs.map(\.index) == Array(0..<UInt64(manifest.jobs.count)) else {
            throw VivoArtifactLoadError.invalidCampaign(
                "candidate and job indices must be unique contiguous ranges"
            )
        }

        var candidateByIndex: [UInt32: VivoCampaignCandidate] = [:]
        for candidate in manifest.candidates {
            guard isFingerprint(candidate.fingerprint),
                  !candidate.assignments.isEmpty,
                  Set(candidate.assignments.map(\.parameterID)).count == candidate.assignments.count,
                  candidate.assignments.allSatisfy({
                      !$0.parameterID.isEmpty && !$0.unit.isEmpty && $0.value.isFinite
                  }) else {
                throw VivoArtifactLoadError.invalidCampaign(
                    "candidate \(candidate.index) is malformed"
                )
            }
            let canonicalAssignments = candidate.assignments.sorted {
                $0.parameterID < $1.parameterID
            }
            let actual = try canonicalFingerprint(canonicalAssignments)
            guard actual == candidate.fingerprint else {
                throw VivoArtifactLoadError.fingerprintMismatch(
                    subject: "campaign candidate \(candidate.index)",
                    expected: candidate.fingerprint,
                    actual: actual
                )
            }
            candidateByIndex[candidate.index] = candidate
        }

        var previous = String(repeating: "0", count: 64)
        for job in manifest.jobs {
            guard let candidate = candidateByIndex[job.candidateIndex] else {
                throw VivoArtifactLoadError.campaignLedgerMismatch(
                    job: job.index,
                    reason: "candidate index is unknown"
                )
            }
            guard job.candidateFingerprint == candidate.fingerprint,
                  job.previousJobDigest == previous,
                  isFingerprint(job.digest) else {
                throw VivoArtifactLoadError.campaignLedgerMismatch(
                    job: job.index,
                    reason: "candidate or previous digest does not match"
                )
            }
            let digest = campaignJobDigest(
                previous: previous,
                index: job.index,
                candidate: job.candidateFingerprint,
                replicate: job.replicateIndex,
                seed: job.seed
            )
            guard digest == job.digest else {
                throw VivoArtifactLoadError.campaignLedgerMismatch(
                    job: job.index,
                    reason: "job digest does not verify"
                )
            }
            previous = digest
        }
        guard previous == manifest.ledgerHead else {
            throw VivoArtifactLoadError.invalidCampaign("ledgerHead does not match the final job")
        }

        let unsigned = VivoCampaignManifest(
            schemaVersion: manifest.schemaVersion,
            definitionFingerprint: manifest.definitionFingerprint,
            name: manifest.name,
            artifacts: manifest.artifacts,
            sampling: manifest.sampling,
            candidates: manifest.candidates,
            jobs: manifest.jobs,
            tags: manifest.tags,
            ledgerHead: manifest.ledgerHead,
            fingerprint: ""
        )
        let actualFingerprint = try canonicalFingerprint(unsigned)
        guard actualFingerprint == manifest.fingerprint else {
            throw VivoArtifactLoadError.fingerprintMismatch(
                subject: "campaign manifest",
                expected: manifest.fingerprint,
                actual: actualFingerprint
            )
        }
    }

    private static func validateStructure(
        _ data: Data,
        limits: VivoArtifactLoadLimits
    ) throws {
        let root: Any
        do {
            root = try JSONSerialization.jsonObject(
                with: data,
                options: [.fragmentsAllowed]
            )
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }

        var nodes = 0
        func visit(_ value: Any, depth: Int) throws {
            guard depth <= limits.maximumDepth else {
                throw VivoArtifactLoadError.structureLimit("depth")
            }
            nodes += 1
            guard nodes <= limits.maximumNodes else {
                throw VivoArtifactLoadError.structureLimit("node count")
            }

            switch value {
            case let string as String:
                guard string.utf8.count <= limits.maximumStringBytes else {
                    throw VivoArtifactLoadError.structureLimit("string bytes")
                }
            case let array as [Any]:
                guard array.count <= limits.maximumArrayElements else {
                    throw VivoArtifactLoadError.structureLimit("array elements")
                }
                for child in array { try visit(child, depth: depth + 1) }
            case let object as [String: Any]:
                guard object.count <= limits.maximumObjectMembers else {
                    throw VivoArtifactLoadError.structureLimit("object members")
                }
                for (key, child) in object {
                    guard key.utf8.count <= limits.maximumStringBytes else {
                        throw VivoArtifactLoadError.structureLimit("object-key bytes")
                    }
                    try visit(child, depth: depth + 1)
                }
            case is NSNull, is NSNumber:
                break
            default:
                throw VivoArtifactLoadError.invalidJSON(
                    "unsupported Foundation JSON node \(type(of: value))"
                )
            }
        }
        try visit(root, depth: 0)
    }

    private static func verifyOptionalGeneratedFingerprint(
        subject: String,
        supplied: String,
        generated: String
    ) throws {
        guard supplied.isEmpty || isFingerprint(supplied) else {
            throw VivoArtifactLoadError.invalidFingerprint(subject)
        }
        if !supplied.isEmpty, supplied != generated {
            throw VivoArtifactLoadError.fingerprintMismatch(
                subject: subject,
                expected: supplied,
                actual: generated
            )
        }
    }

    private static func canonicalFingerprint<T: Encodable>(_ value: T) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        do {
            return hex(SHA256.hash(data: try encoder.encode(value)))
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }
    }

    private static func campaignJobDigest(
        previous: String,
        index: UInt64,
        candidate: String,
        replicate: UInt32,
        seed: UInt64
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-campaign-job-v1".utf8))
        hasher.update(data: Data(previous.utf8))
        update(&hasher, index)
        hasher.update(data: Data(candidate.utf8))
        update(&hasher, replicate)
        update(&hasher, seed)
        return hex(hasher.finalize())
    }

    private static func update<T: FixedWidthInteger>(_ hasher: inout SHA256, _ value: T) {
        var littleEndian = value.littleEndian
        Swift.withUnsafeBytes(of: &littleEndian) { bytes in
            hasher.update(data: Data(bytes))
        }
    }

    private static func isFingerprint(_ value: String) -> Bool {
        value.utf8.count == 64 && value.allSatisfy {
            $0.isHexDigit && !$0.isUppercase
        }
    }

    private static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }
}
