import Foundation

private struct VivoPhysiologyExecutionIdentity: Encodable {
    let numericalVersion: UInt32
    let model: PreparedVivoPhysiologyModel
    let configuration: VivoPhysiologyRuntimeConfiguration
}

/// Preferred standalone resume format. v1 state's source fingerprint alone
/// cannot identify compiled CSR tables, clearance overrides or numerical policy.
public struct VivoPhysiologyResumeCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/physiology-resume-checkpoint/v1"
    public let schema: String
    public let executableFingerprint: VivoFingerprint
    public let state: VivoPhysiologyCheckpoint

    public init(model: PreparedVivoPhysiologyModel,
                configuration: VivoPhysiologyRuntimeConfiguration,
                state: VivoPhysiologyCheckpoint) throws {
        try VivoPhysiologyModelValidator.validate(model)
        try configuration.validate(for: model)
        try state.validate()
        guard state.modelFingerprint == model.fingerprint,
              state.pairCount == model.pairCount,
              state.environmentCount == model.environmentCount else {
            throw VivoRuntimeError.invalidConfiguration("physiology resume shape or source identity mismatch")
        }
        schema = Self.schema
        executableFingerprint = try Self.identity(model: model, configuration: configuration)
        self.state = state
    }
    public func fingerprint() throws -> VivoFingerprint {
        try state.validate()
        guard schema == Self.schema else {
            throw VivoRuntimeError.invalidConfiguration("unknown physiology resume schema")
        }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
    static func identity(model: PreparedVivoPhysiologyModel,
                         configuration: VivoPhysiologyRuntimeConfiguration) throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
            VivoPhysiologyExecutionIdentity(numericalVersion: 2, model: model, configuration: configuration)
        ))
    }
}

extension VivoPhysiologyRuntime {
    public func resumeCheckpoint() async throws -> VivoPhysiologyResumeCheckpoint {
        let captured = try await checkpoint()
        return try .init(model: model, configuration: configuration, state: captured)
    }
    public func restore(_ checkpoint: VivoPhysiologyResumeCheckpoint) throws {
        guard checkpoint.schema == VivoPhysiologyResumeCheckpoint.schema,
              try checkpoint.executableFingerprint == VivoPhysiologyResumeCheckpoint.identity(
                model: model, configuration: configuration
              ) else {
            throw VivoRuntimeError.invalidConfiguration("physiology prepared tables or numerical settings differ from resume identity")
        }
        try restore(checkpoint.state)
    }
}
