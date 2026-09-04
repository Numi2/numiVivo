import CryptoKit
import Foundation

public enum VivoGeneratedFingerprintSourceCompiler {
    /// Compiles an editable JSON object into a typed artifact whose canonical
    /// identity is the SHA-256 digest of the same typed object with an empty
    /// `fingerprint` property. This is the identity convention used by the
    /// current Swift model artifacts.
    public static func compile<T: Codable & Sendable>(
        _ type: T.Type,
        from source: Data,
        limits: VivoArtifactLoadLimits = .init(),
        validate: @Sendable (T) throws -> Void
    ) throws -> VivoCompiledSourceArtifact<T> {
        try limits.validate()
        guard source.count <= limits.maximumBytes else {
            throw VivoArtifactLoadError.fileTooLarge(
                actual: source.count,
                maximum: limits.maximumBytes
            )
        }
        try validateStructure(source, limits: limits)

        let root: Any
        do {
            root = try JSONSerialization.jsonObject(with: source)
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }
        guard var sourceObject = root as? [String: Any] else {
            throw VivoArtifactLoadError.invalidJSON("source root must be a JSON object")
        }
        if sourceObject["schemaVersion"] == nil {
            sourceObject["schemaVersion"] = 1
        }
        sourceObject["fingerprint"] = ""

        let unsignedData: Data
        do {
            unsignedData = try JSONSerialization.data(
                withJSONObject: sourceObject,
                options: [.sortedKeys, .withoutEscapingSlashes]
            )
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }

        let unsigned: T
        do {
            unsigned = try JSONDecoder().decode(T.self, from: unsignedData)
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }

        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        let canonicalUnsigned: Data
        do {
            canonicalUnsigned = try encoder.encode(unsigned)
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }
        let artifactFingerprint = hex(SHA256.hash(data: canonicalUnsigned))
        sourceObject["fingerprint"] = artifactFingerprint

        let compiledData: Data
        do {
            compiledData = try JSONSerialization.data(
                withJSONObject: sourceObject,
                options: [.sortedKeys, .withoutEscapingSlashes]
            )
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }
        let compiled: T
        do {
            compiled = try JSONDecoder().decode(T.self, from: compiledData)
            try validate(compiled)
        } catch let error as VivoArtifactLoadError {
            throw error
        } catch {
            throw VivoArtifactLoadError.typeValidation(error.localizedDescription)
        }

        let canonicalCompiled: Data
        do {
            canonicalCompiled = try encoder.encode(compiled)
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }
        guard let finalObject = try JSONSerialization.jsonObject(with: canonicalCompiled) as? [String: Any],
              finalObject["fingerprint"] as? String == artifactFingerprint else {
            throw VivoArtifactLoadError.fingerprintMismatch(
                subject: String(describing: T.self),
                expected: artifactFingerprint,
                actual: "missing-or-rewritten"
            )
        }

        return .init(
            value: compiled,
            sourceFingerprint: hex(SHA256.hash(data: source)),
            artifactFingerprint: artifactFingerprint,
            sourceBytes: source.count
        )
    }

    public static func exactSSAModel(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoExactSSAModel> {
        try compile(VivoExactSSAModel.self, from: source, limits: limits) { model in
            try model.validate()
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
        var nodeCount = 0
        func visit(_ value: Any, depth: Int) throws {
            guard depth <= limits.maximumDepth else {
                throw VivoArtifactLoadError.structureLimit("depth")
            }
            nodeCount += 1
            guard nodeCount <= limits.maximumNodes else {
                throw VivoArtifactLoadError.structureLimit("node count")
            }
            switch value {
            case let value as String:
                guard value.utf8.count <= limits.maximumStringBytes else {
                    throw VivoArtifactLoadError.structureLimit("string bytes")
                }
            case let value as [Any]:
                guard value.count <= limits.maximumArrayElements else {
                    throw VivoArtifactLoadError.structureLimit("array elements")
                }
                for child in value { try visit(child, depth: depth + 1) }
            case let value as [String: Any]:
                guard value.count <= limits.maximumObjectMembers else {
                    throw VivoArtifactLoadError.structureLimit("object members")
                }
                for (key, child) in value {
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

    private static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }
}
