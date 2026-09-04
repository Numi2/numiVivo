import CryptoKit
import Foundation

public struct VivoCompiledSourceArtifact<Value: Sendable>: Sendable {
    public let value: Value
    public let sourceFingerprint: String
    public let artifactFingerprint: String
    public let sourceBytes: Int

    public init(
        value: Value,
        sourceFingerprint: String,
        artifactFingerprint: String,
        sourceBytes: Int
    ) {
        self.value = value
        self.sourceFingerprint = sourceFingerprint
        self.artifactFingerprint = artifactFingerprint
        self.sourceBytes = sourceBytes
    }
}

public enum VivoSourceArtifactCompiler {
    public static func physiologicalPartitionModel(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoPhysiologicalPartitionModel> {
        let decoded: VivoPhysiologicalPartitionModel = try decodeFingerprintSource(
            VivoPhysiologicalPartitionModel.self,
            from: source,
            limits: limits
        )
        let compiled = try VivoPhysiologicalPartitionModel(
            name: decoded.name,
            compartments: decoded.compartments,
            analytes: decoded.analytes,
            edges: decoded.edges
        )
        try compiled.validate()
        return .init(
            value: compiled,
            sourceFingerprint: sourceFingerprint(source),
            artifactFingerprint: compiled.fingerprint,
            sourceBytes: source.count
        )
    }

    public static func populationModel(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoPopulationModel> {
        let decoded: VivoPopulationModel = try decodeFingerprintSource(
            VivoPopulationModel.self,
            from: source,
            limits: limits
        )
        let compiled = try VivoPopulationModel(
            name: decoded.name,
            grid: decoded.grid,
            regulatorFields: decoded.regulatorFields,
            phenotypes: decoded.phenotypes,
            transitions: decoded.transitions,
            interactions: decoded.interactions
        )
        try compiled.validate()
        return .init(
            value: compiled,
            sourceFingerprint: sourceFingerprint(source),
            artifactFingerprint: compiled.fingerprint,
            sourceBytes: source.count
        )
    }

    public static func surrogateContract(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoSurrogateContract> {
        let decoded: VivoSurrogateContract = try decodeFingerprintSource(
            VivoSurrogateContract.self,
            from: source,
            limits: limits
        )
        let compiled = try VivoSurrogateContract(
            id: decoded.id,
            modelFingerprint: decoded.modelFingerprint,
            trainingDataFingerprint: decoded.trainingDataFingerprint,
            mechanismPackFingerprint: decoded.mechanismPackFingerprint,
            hostContextFingerprint: decoded.hostContextFingerprint,
            inputs: decoded.inputs,
            outputs: decoded.outputs,
            uncertaintyKind: decoded.uncertaintyKind,
            maximumNormalizedUncertainty: decoded.maximumNormalizedUncertainty,
            maximumNormalizedExtrapolation: decoded.maximumNormalizedExtrapolation,
            maximumConsecutiveAcceptedSteps: decoded.maximumConsecutiveAcceptedSteps,
            mandatoryAuthorityInterval: decoded.mandatoryAuthorityInterval
        )
        try compiled.validate()
        return .init(
            value: compiled,
            sourceFingerprint: sourceFingerprint(source),
            artifactFingerprint: compiled.fingerprint,
            sourceBytes: source.count
        )
    }

    public static func campaignDefinition(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoCampaignDefinition> {
        let loaded = try VivoValidatedArtifactLoader.decode(
            VivoCampaignDefinition.self,
            from: source,
            limits: limits
        ) { definition in
            try definition.validate()
        }
        let fingerprint = try canonicalFingerprint(loaded.value)
        return .init(
            value: loaded.value,
            sourceFingerprint: loaded.sourceFingerprint,
            artifactFingerprint: fingerprint,
            sourceBytes: loaded.sourceBytes
        )
    }

    public static func compileCampaign(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init()
    ) throws -> VivoCompiledSourceArtifact<VivoCampaignManifest> {
        let definition = try campaignDefinition(from: source, limits: limits)
        let manifest = try VivoCampaignCompiler().compile(definition.value)
        return .init(
            value: manifest,
            sourceFingerprint: definition.sourceFingerprint,
            artifactFingerprint: manifest.fingerprint,
            sourceBytes: definition.sourceBytes
        )
    }

    public static func canonicalJSON<T: Encodable & Sendable>(
        _ artifact: VivoCompiledSourceArtifact<T>
    ) throws -> Data {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        do {
            return try encoder.encode(artifact.value)
        } catch {
            throw VivoArtifactLoadError.decoding(error.localizedDescription)
        }
    }

    private static func decodeFingerprintSource<T: Decodable & Sendable>(
        _ type: T.Type,
        from source: Data,
        limits: VivoArtifactLoadLimits
    ) throws -> T {
        try limits.validate()
        guard source.count <= limits.maximumBytes else {
            throw VivoArtifactLoadError.fileTooLarge(
                actual: source.count,
                maximum: limits.maximumBytes
            )
        }
        let root: Any
        do {
            root = try JSONSerialization.jsonObject(with: source)
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }
        guard var object = root as? [String: Any] else {
            throw VivoArtifactLoadError.invalidJSON("source root must be an object")
        }
        if object["schemaVersion"] == nil {
            object["schemaVersion"] = 1
        }
        if object["fingerprint"] == nil {
            object["fingerprint"] = ""
        }
        let normalized: Data
        do {
            normalized = try JSONSerialization.data(
                withJSONObject: object,
                options: [.sortedKeys, .withoutEscapingSlashes]
            )
        } catch {
            throw VivoArtifactLoadError.invalidJSON(error.localizedDescription)
        }
        let loaded = try VivoValidatedArtifactLoader.decode(
            T.self,
            from: normalized,
            limits: limits
        )
        return loaded.value
    }

    private static func sourceFingerprint(_ source: Data) -> String {
        hex(SHA256.hash(data: source))
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

    private static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }
}
