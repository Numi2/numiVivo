import Foundation

/// A protocol-stage change is not checkpoint restoration: changing thermostat,
/// pressure, or numerical settings changes the configuration identity. Preserve
/// positions, velocities and physical time explicitly, and start a new runtime
/// whose configuration identity names the new stage.
public struct VivoMDStageTransfer: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/md-stage-transfer/v1"
    public let schema: String
    public let sourceStage: String
    public let destinationStage: String
    public let systemFingerprint: VivoFingerprint
    public let sourceConfigurationFingerprint: VivoFingerprint
    public let destinationConfigurationFingerprint: VivoFingerprint
    public let sourceAcceptedStep: UInt64
    public let timePS: Double
    public let preservesVelocities: Bool

    public init(sourceStage: String, destinationStage: String,
                state: VivoMDStateSnapshot,
                destinationConfiguration: VivoMDConfiguration) throws {
        guard !sourceStage.isEmpty, !destinationStage.isEmpty,
              sourceStage != destinationStage, !state.positionsNM.isEmpty else {
            throw VivoArtifactValidationError.invalid("MD stage transfer requires distinct stage identifiers and nonempty state")
        }
        try state.validate(particleCount: state.positionsNM.count)
        self.schema = Self.schema
        self.sourceStage = sourceStage
        self.destinationStage = destinationStage
        self.systemFingerprint = state.systemFingerprint
        self.sourceConfigurationFingerprint = state.configurationFingerprint
        self.destinationConfigurationFingerprint = try destinationConfiguration.fingerprint()
        self.sourceAcceptedStep = state.stepIndex
        self.timePS = state.timePS
        self.preservesVelocities = true
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}
