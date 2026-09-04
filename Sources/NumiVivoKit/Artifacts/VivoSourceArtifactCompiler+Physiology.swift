import CryptoKit
import Foundation

public extension VivoSourceArtifactCompiler {
    static func physiologyModel(
        from source: Data,
        limits: VivoArtifactLoadLimits = .init(),
        compiler: VivoPhysiologyCompiler = .init()
    ) throws -> VivoCompiledSourceArtifact<PreparedVivoPhysiologyModel> {
        let loaded = try VivoValidatedArtifactLoader.decode(
            VivoPhysiologyPack.self,
            from: source,
            limits: limits
        ) { model in
            try model.validateStructure()
        }
        let prepared = try compiler.compile(loaded.value)
        return VivoCompiledSourceArtifact(
            value: prepared,
            sourceFingerprint: loaded.sourceFingerprint,
            artifactFingerprint: prepared.fingerprint,
            sourceBytes: loaded.sourceBytes
        )
    }

    static func molecularPhysiologyCoupling(
        from source: Data,
        programPack: VivoProgramPack,
        physiology: PreparedVivoPhysiologyModel,
        molecularLaneCount: UInt32,
        limits: VivoArtifactLoadLimits = .init(),
        compiler: VivoMolecularPhysiologyCouplingCompiler = .init()
    ) throws -> VivoCompiledSourceArtifact<PreparedVivoMolecularPhysiologyBridge> {
        let loaded = try VivoValidatedArtifactLoader.decode(
            VivoMolecularPhysiologyCouplingPack.self,
            from: source,
            limits: limits
        ) { coupling in
            try coupling.policy.validate()
        }
        let prepared = try compiler.compile(
            loaded.value,
            programPack: programPack,
            physiology: physiology,
            molecularLaneCount: molecularLaneCount
        )
        return VivoCompiledSourceArtifact(
            value: prepared,
            sourceFingerprint: loaded.sourceFingerprint,
            artifactFingerprint: prepared.fingerprint,
            sourceBytes: loaded.sourceBytes
        )
    }
}
