import CryptoKit
import Foundation

public enum VivoCheckpointRuntimeKind: String, Codable, CaseIterable, Sendable {
    case programRuntime
    case exactSSA
    case population
    case physiologicalFlow
    case physiologicalPartition
    case coupledNumiLab
}

public enum VivoCheckpointSectionEncoding: String, Codable, CaseIterable, Sendable {
    case rawBytes
    case float32LittleEndian
    case uint32LittleEndian
    case canonicalJSON
}

public struct VivoCheckpointSection: Codable, Hashable, Sendable {
    public let id: String
    public let encoding: VivoCheckpointSectionEncoding
    public let offset: UInt64
    public let length: UInt64
    public let elementCount: UInt64
    public let elementStride: UInt32
    public let fingerprint: String

    public init(
        id: String,
        encoding: VivoCheckpointSectionEncoding,
        offset: UInt64,
        length: UInt64,
        elementCount: UInt64,
        elementStride: UInt32,
        fingerprint: String
    ) {
        self.id = id
        self.encoding = encoding
        self.offset = offset
        self.length = length
        self.elementCount = elementCount
        self.elementStride = elementStride
        self.fingerprint = fingerprint
    }
}

public struct VivoCheckpointManifest: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let runtime: VivoCheckpointRuntimeKind
    public let artifactFingerprint: String
    public let sourceFingerprint: String?
    public let hostContextFingerprint: String?
    public let experimentFingerprint: String?
    public let couplingFingerprint: String?
    public let parentCheckpointFingerprint: String?
    public let stepIndex: UInt64
    public let logicalTime: Double
    public let randomStreamVersion: UInt32
    public let sections: [VivoCheckpointSection]
    public let metadata: [String: String]
    public let signingKeyID: String?

    public init(
        runtime: VivoCheckpointRuntimeKind,
        artifactFingerprint: String,
        sourceFingerprint: String? = nil,
        hostContextFingerprint: String? = nil,
        experimentFingerprint: String? = nil,
        couplingFingerprint: String? = nil,
        parentCheckpointFingerprint: String? = nil,
        stepIndex: UInt64,
        logicalTime: Double,
        randomStreamVersion: UInt32,
        sections: [VivoCheckpointSection],
        metadata: [String: String] = [:],
        signingKeyID: String? = nil
    ) {
        self.schemaVersion = 1
        self.runtime = runtime
        self.artifactFingerprint = artifactFingerprint
        self.sourceFingerprint = sourceFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.parentCheckpointFingerprint = parentCheckpointFingerprint
        self.stepIndex = stepIndex
        self.logicalTime = logicalTime
        self.randomStreamVersion = randomStreamVersion
        self.sections = sections
        self.metadata = metadata
        self.signingKeyID = signingKeyID
    }
}

public struct VivoCheckpointSectionPayload: Sendable {
    public let id: String
    public let encoding: VivoCheckpointSectionEncoding
    public let elementCount: UInt64
    public let elementStride: UInt32
    public let data: Data

    public init(
        id: String,
        encoding: VivoCheckpointSectionEncoding,
        elementCount: UInt64,
        elementStride: UInt32,
        data: Data
    ) {
        self.id = id
        self.encoding = encoding
        self.elementCount = elementCount
        self.elementStride = elementStride
        self.data = data
    }

    public static func float32(id: String, values: [Float]) -> Self {
        var littleEndian = values.map { $0.bitPattern.littleEndian }
        return littleEndian.withUnsafeBytes { bytes in
            .init(
                id: id,
                encoding: .float32LittleEndian,
                elementCount: UInt64(values.count),
                elementStride: UInt32(MemoryLayout<UInt32>.stride),
                data: Data(bytes)
            )
        }
    }

    public static func uint32(id: String, values: [UInt32]) -> Self {
        var littleEndian = values.map(\.littleEndian)
        return littleEndian.withUnsafeBytes { bytes in
            .init(
                id: id,
                encoding: .uint32LittleEndian,
                elementCount: UInt64(values.count),
                elementStride: UInt32(MemoryLayout<UInt32>.stride),
                data: Data(bytes)
            )
        }
    }

    public static func canonicalJSON<T: Encodable & Sendable>(id: String, value: T) throws -> Self {
        let data = try VivoCheckpointCodec.canonicalEncoder.encode(value)
        return .init(
            id: id,
            encoding: .canonicalJSON,
            elementCount: 1,
            elementStride: UInt32(clamping: data.count),
            data: data
        )
    }
}

public struct VivoCheckpointBuildRequest: Sendable {
    public let runtime: VivoCheckpointRuntimeKind
    public let artifactFingerprint: String
    public let sourceFingerprint: String?
    public let hostContextFingerprint: String?
    public let experimentFingerprint: String?
    public let couplingFingerprint: String?
    public let parentCheckpointFingerprint: String?
    public let stepIndex: UInt64
    public let logicalTime: Double
    public let randomStreamVersion: UInt32
    public let sections: [VivoCheckpointSectionPayload]
    public let metadata: [String: String]

    public init(
        runtime: VivoCheckpointRuntimeKind,
        artifactFingerprint: String,
        sourceFingerprint: String? = nil,
        hostContextFingerprint: String? = nil,
        experimentFingerprint: String? = nil,
        couplingFingerprint: String? = nil,
        parentCheckpointFingerprint: String? = nil,
        stepIndex: UInt64,
        logicalTime: Double,
        randomStreamVersion: UInt32 = 1,
        sections: [VivoCheckpointSectionPayload],
        metadata: [String: String] = [:]
    ) {
        self.runtime = runtime
        self.artifactFingerprint = artifactFingerprint
        self.sourceFingerprint = sourceFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.parentCheckpointFingerprint = parentCheckpointFingerprint
        self.stepIndex = stepIndex
        self.logicalTime = logicalTime
        self.randomStreamVersion = randomStreamVersion
        self.sections = sections
        self.metadata = metadata
    }
}

public struct VivoCheckpointSigner: Sendable {
    public let keyID: String
    private let privateKey: Curve25519.Signing.PrivateKey

    public init(keyID: String, privateKey: Curve25519.Signing.PrivateKey) throws {
        guard !keyID.isEmpty, keyID.utf8.count <= 256 else {
            throw VivoCheckpointError.invalidRequest("signing key identifier is invalid")
        }
        self.keyID = keyID
        self.privateKey = privateKey
    }

    public var publicKey: Curve25519.Signing.PublicKey { privateKey.publicKey }

    func sign(_ digest: Data) throws -> Data {
        try privateKey.signature(for: digest)
    }
}

public struct VivoCheckpointLimits: Sendable, Equatable {
    public var maximumTotalBytes: UInt64
    public var maximumManifestBytes: UInt64
    public var maximumSections: UInt32
    public var maximumMetadataEntries: UInt32
    public var maximumIdentifierBytes: UInt32

    public init(
        maximumTotalBytes: UInt64 = 64 * 1_024 * 1_024 * 1_024,
        maximumManifestBytes: UInt64 = 16 * 1_024 * 1_024,
        maximumSections: UInt32 = 4_096,
        maximumMetadataEntries: UInt32 = 4_096,
        maximumIdentifierBytes: UInt32 = 4_096
    ) {
        self.maximumTotalBytes = maximumTotalBytes
        self.maximumManifestBytes = maximumManifestBytes
        self.maximumSections = maximumSections
        self.maximumMetadataEntries = maximumMetadataEntries
        self.maximumIdentifierBytes = maximumIdentifierBytes
    }
}

public struct VivoDecodedCheckpoint: Sendable {
    public let manifest: VivoCheckpointManifest
    public let packageFingerprint: String
    public let signature: Data?
    private let payload: Data

    init(
        manifest: VivoCheckpointManifest,
        packageFingerprint: String,
        signature: Data?,
        payload: Data
    ) {
        self.manifest = manifest
        self.packageFingerprint = packageFingerprint
        self.signature = signature
        self.payload = payload
    }

    public func section(id: String) throws -> Data {
        guard let section = manifest.sections.first(where: { $0.id == id }) else {
            throw VivoCheckpointError.missingSection(id)
        }
        let start = Int(section.offset)
        let length = Int(section.length)
        guard start >= 0, length >= 0, start <= payload.count, length <= payload.count - start else {
            throw VivoCheckpointError.invalidSection(id, "range is outside the payload")
        }
        return payload.subdata(in: start..<(start + length))
    }

    public func float32Section(id: String) throws -> [Float] {
        guard let descriptor = manifest.sections.first(where: { $0.id == id }),
              descriptor.encoding == .float32LittleEndian,
              descriptor.elementStride == 4 else {
            throw VivoCheckpointError.invalidSection(id, "section is not float32 little-endian")
        }
        let data = try section(id: id)
        guard data.count.isMultiple(of: 4), UInt64(data.count / 4) == descriptor.elementCount else {
            throw VivoCheckpointError.invalidSection(id, "float section size does not match its descriptor")
        }
        return data.withUnsafeBytes { raw in
            let words = raw.bindMemory(to: UInt32.self)
            return words.map { Float(bitPattern: UInt32(littleEndian: $0)) }
        }
    }

    public func uint32Section(id: String) throws -> [UInt32] {
        guard let descriptor = manifest.sections.first(where: { $0.id == id }),
              descriptor.encoding == .uint32LittleEndian,
              descriptor.elementStride == 4 else {
            throw VivoCheckpointError.invalidSection(id, "section is not uint32 little-endian")
        }
        let data = try section(id: id)
        guard data.count.isMultiple(of: 4), UInt64(data.count / 4) == descriptor.elementCount else {
            throw VivoCheckpointError.invalidSection(id, "uint32 section size does not match its descriptor")
        }
        return data.withUnsafeBytes { raw in
            raw.bindMemory(to: UInt32.self).map(UInt32.init(littleEndian:))
        }
    }
}

public enum VivoCheckpointError: Error, LocalizedError, Sendable {
    case invalidRequest(String)
    case invalidHeader(String)
    case unsupportedVersion(UInt16, UInt16)
    case resourceLimit(String)
    case fingerprintMismatch(String)
    case signatureRequired
    case signatureInvalid
    case signingKeyMismatch
    case duplicateSection(String)
    case missingSection(String)
    case invalidSection(String, String)
    case arithmeticOverflow
    case decoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidRequest(let reason): return "invalid checkpoint request: \(reason)"
        case .invalidHeader(let reason): return "invalid checkpoint header: \(reason)"
        case .unsupportedVersion(let major, let minor): return "unsupported checkpoint version \(major).\(minor)"
        case .resourceLimit(let reason): return "checkpoint resource limit exceeded: \(reason)"
        case .fingerprintMismatch(let subject): return "checkpoint fingerprint mismatch: \(subject)"
        case .signatureRequired: return "checkpoint signature is required"
        case .signatureInvalid: return "checkpoint signature is invalid"
        case .signingKeyMismatch: return "checkpoint signing key does not match the manifest"
        case .duplicateSection(let id): return "duplicate checkpoint section \(id)"
        case .missingSection(let id): return "missing checkpoint section \(id)"
        case .invalidSection(let id, let reason): return "invalid checkpoint section \(id): \(reason)"
        case .arithmeticOverflow: return "checkpoint arithmetic overflow"
        case .decoding(let reason): return "checkpoint decoding failed: \(reason)"
        }
    }
}

public enum VivoCheckpointCodec {
    static let canonicalEncoder: JSONEncoder = {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return encoder
    }()

    private static let magic = Data([0x4e, 0x56, 0x43, 0x48, 0x4b, 0x50, 0x54, 0x31]) // NVCHKPT1
    private static let major: UInt16 = 1
    private static let minor: UInt16 = 0
    private static let fixedHeaderBytes: UInt32 = 128
    private static let signedFlag: UInt32 = 1 << 0

    public static func encode(
        _ request: VivoCheckpointBuildRequest,
        signer: VivoCheckpointSigner? = nil,
        limits: VivoCheckpointLimits = .init()
    ) throws -> Data {
        try validate(request: request, limits: limits)

        var payload = Data()
        var descriptors: [VivoCheckpointSection] = []
        descriptors.reserveCapacity(request.sections.count)
        for section in request.sections {
            guard payload.count <= Int.max - section.data.count else {
                throw VivoCheckpointError.arithmeticOverflow
            }
            let offset = UInt64(payload.count)
            payload.append(section.data)
            descriptors.append(.init(
                id: section.id,
                encoding: section.encoding,
                offset: offset,
                length: UInt64(section.data.count),
                elementCount: section.elementCount,
                elementStride: section.elementStride,
                fingerprint: hex(SHA256.hash(data: section.data))
            ))
        }

        let manifest = VivoCheckpointManifest(
            runtime: request.runtime,
            artifactFingerprint: request.artifactFingerprint,
            sourceFingerprint: request.sourceFingerprint,
            hostContextFingerprint: request.hostContextFingerprint,
            experimentFingerprint: request.experimentFingerprint,
            couplingFingerprint: request.couplingFingerprint,
            parentCheckpointFingerprint: request.parentCheckpointFingerprint,
            stepIndex: request.stepIndex,
            logicalTime: request.logicalTime,
            randomStreamVersion: request.randomStreamVersion,
            sections: descriptors,
            metadata: request.metadata,
            signingKeyID: signer?.keyID
        )
        let manifestData = try canonicalEncoder.encode(manifest)
        guard UInt64(manifestData.count) <= limits.maximumManifestBytes else {
            throw VivoCheckpointError.resourceLimit("manifest")
        }

        var authenticated = Data()
        authenticated.append(manifestData)
        authenticated.append(payload)
        let digest = Data(SHA256.hash(data: authenticated))
        let signature = try signer?.sign(digest) ?? Data()

        let total = UInt64(fixedHeaderBytes) + UInt64(manifestData.count) + UInt64(payload.count) + UInt64(signature.count)
        guard total <= limits.maximumTotalBytes, total <= UInt64(Int.max) else {
            throw VivoCheckpointError.resourceLimit("total bytes")
        }

        var header = Data(capacity: Int(fixedHeaderBytes))
        header.append(magic)
        header.appendLittleEndian(major)
        header.appendLittleEndian(minor)
        header.appendLittleEndian(fixedHeaderBytes)
        header.appendLittleEndian(signer == nil ? UInt32(0) : signedFlag)
        header.appendLittleEndian(UInt64(manifestData.count))
        header.appendLittleEndian(UInt64(payload.count))
        header.appendLittleEndian(UInt32(signature.count))
        header.appendLittleEndian(UInt32(descriptors.count))
        header.appendLittleEndian(total)
        header.append(digest)
        header.append(Data(repeating: 0, count: Int(fixedHeaderBytes) - header.count))

        var package = Data(capacity: Int(total))
        package.append(header)
        package.append(manifestData)
        package.append(payload)
        package.append(signature)
        return package
    }

    public static func decode(
        _ data: Data,
        verificationKey: Curve25519.Signing.PublicKey? = nil,
        requireSignature: Bool = false,
        limits: VivoCheckpointLimits = .init()
    ) throws -> VivoDecodedCheckpoint {
        guard UInt64(data.count) <= limits.maximumTotalBytes else {
            throw VivoCheckpointError.resourceLimit("total bytes")
        }
        guard data.count >= Int(fixedHeaderBytes) else {
            throw VivoCheckpointError.invalidHeader("input is shorter than the fixed header")
        }
        var cursor = VivoDataCursor(data: data)
        guard try cursor.read(count: magic.count) == magic else {
            throw VivoCheckpointError.invalidHeader("magic")
        }
        let fileMajor: UInt16 = try cursor.readLittleEndian()
        let fileMinor: UInt16 = try cursor.readLittleEndian()
        guard fileMajor == major else {
            throw VivoCheckpointError.unsupportedVersion(fileMajor, fileMinor)
        }
        let headerBytes: UInt32 = try cursor.readLittleEndian()
        let flags: UInt32 = try cursor.readLittleEndian()
        let manifestBytes: UInt64 = try cursor.readLittleEndian()
        let payloadBytes: UInt64 = try cursor.readLittleEndian()
        let signatureBytes: UInt32 = try cursor.readLittleEndian()
        let sectionCount: UInt32 = try cursor.readLittleEndian()
        let totalBytes: UInt64 = try cursor.readLittleEndian()
        let expectedDigest = try cursor.read(count: 32)

        guard headerBytes == fixedHeaderBytes,
              totalBytes == UInt64(data.count),
              manifestBytes <= limits.maximumManifestBytes,
              sectionCount <= limits.maximumSections else {
            throw VivoCheckpointError.invalidHeader("length or section-count fields")
        }
        let computedTotal = UInt64(headerBytes)
            .addingReportingOverflow(manifestBytes)
        guard !computedTotal.overflow else { throw VivoCheckpointError.arithmeticOverflow }
        let withPayload = computedTotal.partialValue.addingReportingOverflow(payloadBytes)
        guard !withPayload.overflow else { throw VivoCheckpointError.arithmeticOverflow }
        let withSignature = withPayload.partialValue.addingReportingOverflow(UInt64(signatureBytes))
        guard !withSignature.overflow, withSignature.partialValue == totalBytes else {
            throw VivoCheckpointError.invalidHeader("component lengths do not sum to totalBytes")
        }

        let signed = (flags & signedFlag) != 0
        guard signed == (signatureBytes > 0) else {
            throw VivoCheckpointError.invalidHeader("signature flag and length disagree")
        }
        if requireSignature && !signed { throw VivoCheckpointError.signatureRequired }

        let manifestStart = Int(headerBytes)
        let manifestEnd = manifestStart + Int(manifestBytes)
        let payloadEnd = manifestEnd + Int(payloadBytes)
        let signatureEnd = payloadEnd + Int(signatureBytes)
        guard signatureEnd == data.count else {
            throw VivoCheckpointError.invalidHeader("component range")
        }
        let manifestData = data.subdata(in: manifestStart..<manifestEnd)
        let payload = data.subdata(in: manifestEnd..<payloadEnd)
        let signature = signed ? data.subdata(in: payloadEnd..<signatureEnd) : nil

        var authenticated = Data()
        authenticated.append(manifestData)
        authenticated.append(payload)
        let actualDigest = Data(SHA256.hash(data: authenticated))
        guard actualDigest == expectedDigest else {
            throw VivoCheckpointError.fingerprintMismatch("package")
        }
        if let verificationKey {
            guard let signature else { throw VivoCheckpointError.signatureRequired }
            guard verificationKey.isValidSignature(signature, for: actualDigest) else {
                throw VivoCheckpointError.signatureInvalid
            }
        }

        let manifest: VivoCheckpointManifest
        do {
            manifest = try JSONDecoder().decode(VivoCheckpointManifest.self, from: manifestData)
        } catch {
            throw VivoCheckpointError.decoding(error.localizedDescription)
        }
        guard manifest.schemaVersion == 1,
              manifest.sections.count == Int(sectionCount),
              manifest.metadata.count <= Int(limits.maximumMetadataEntries),
              manifest.logicalTime.isFinite else {
            throw VivoCheckpointError.invalidHeader("manifest contract")
        }
        if signed, manifest.signingKeyID == nil { throw VivoCheckpointError.signingKeyMismatch }
        if !signed, manifest.signingKeyID != nil { throw VivoCheckpointError.signingKeyMismatch }

        try validateSections(manifest.sections, payload: payload, limits: limits)
        return .init(
            manifest: manifest,
            packageFingerprint: hex(actualDigest),
            signature: signature,
            payload: payload
        )
    }

    private static func validate(
        request: VivoCheckpointBuildRequest,
        limits: VivoCheckpointLimits
    ) throws {
        guard !request.artifactFingerprint.isEmpty,
              request.artifactFingerprint.utf8.count <= Int(limits.maximumIdentifierBytes),
              request.logicalTime.isFinite,
              request.logicalTime >= 0,
              request.randomStreamVersion > 0,
              request.sections.count <= Int(limits.maximumSections),
              request.metadata.count <= Int(limits.maximumMetadataEntries) else {
            throw VivoCheckpointError.invalidRequest("identity, time, stream version or collection size")
        }
        var identifiers = Set<String>()
        for section in request.sections {
            guard !section.id.isEmpty,
                  section.id.utf8.count <= Int(limits.maximumIdentifierBytes),
                  identifiers.insert(section.id).inserted else {
                throw VivoCheckpointError.duplicateSection(section.id)
            }
            let expected = section.elementCount.multipliedReportingOverflow(by: UInt64(section.elementStride))
            guard !expected.overflow,
                  expected.partialValue == UInt64(section.data.count) || section.encoding == .canonicalJSON else {
                throw VivoCheckpointError.invalidSection(section.id, "element count and stride do not match data length")
            }
        }
        for (key, value) in request.metadata {
            guard key.utf8.count <= Int(limits.maximumIdentifierBytes),
                  value.utf8.count <= Int(limits.maximumIdentifierBytes) else {
                throw VivoCheckpointError.resourceLimit("metadata string")
            }
        }
    }

    private static func validateSections(
        _ sections: [VivoCheckpointSection],
        payload: Data,
        limits: VivoCheckpointLimits
    ) throws {
        var seen = Set<String>()
        var ranges: [(UInt64, UInt64, String)] = []
        for section in sections {
            guard !section.id.isEmpty,
                  section.id.utf8.count <= Int(limits.maximumIdentifierBytes),
                  seen.insert(section.id).inserted else {
                throw VivoCheckpointError.duplicateSection(section.id)
            }
            let end = section.offset.addingReportingOverflow(section.length)
            guard !end.overflow, end.partialValue <= UInt64(payload.count) else {
                throw VivoCheckpointError.invalidSection(section.id, "range exceeds payload")
            }
            if section.encoding != .canonicalJSON {
                let expected = section.elementCount.multipliedReportingOverflow(by: UInt64(section.elementStride))
                guard !expected.overflow, expected.partialValue == section.length else {
                    throw VivoCheckpointError.invalidSection(section.id, "element layout mismatch")
                }
            }
            let bytes = payload.subdata(in: Int(section.offset)..<Int(end.partialValue))
            guard hex(SHA256.hash(data: bytes)) == section.fingerprint else {
                throw VivoCheckpointError.fingerprintMismatch(section.id)
            }
            ranges.append((section.offset, end.partialValue, section.id))
        }
        ranges.sort { $0.0 < $1.0 }
        for index in 1..<ranges.count where ranges[index].0 < ranges[index - 1].1 {
            throw VivoCheckpointError.invalidSection(ranges[index].2, "section overlaps a previous section")
        }
    }

    private static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }
}

private struct VivoDataCursor {
    let data: Data
    private(set) var offset = 0

    mutating func read(count: Int) throws -> Data {
        guard count >= 0, offset <= data.count, count <= data.count - offset else {
            throw VivoCheckpointError.invalidHeader("truncated field")
        }
        defer { offset += count }
        return data.subdata(in: offset..<(offset + count))
    }

    mutating func readLittleEndian<T: FixedWidthInteger>() throws -> T {
        let bytes = try read(count: MemoryLayout<T>.size)
        return bytes.withUnsafeBytes { raw in
            T(littleEndian: raw.loadUnaligned(as: T.self))
        }
    }
}

private extension Data {
    mutating func appendLittleEndian<T: FixedWidthInteger>(_ value: T) {
        var littleEndian = value.littleEndian
        Swift.withUnsafeBytes(of: &littleEndian) { append(contentsOf: $0) }
    }
}
