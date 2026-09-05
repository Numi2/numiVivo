import Foundation

/// A stage-transfer description and its single preparation entry point. This
/// consolidates the previous same-named struct/enum without dropping either the
/// snapshot-description initializer or the checkpoint-based protocol API.
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
        try state.validate(particleCount: state.positionsNM.count)
        schema = Self.schema
        self.sourceStage = sourceStage
        self.destinationStage = destinationStage
        systemFingerprint = state.systemFingerprint
        sourceConfigurationFingerprint = state.configurationFingerprint
        destinationConfigurationFingerprint = try destinationConfiguration.fingerprint()
        sourceAcceptedStep = state.stepIndex
        timePS = state.timePS
        preservesVelocities = true
        try validate()
    }

    public func validate() throws {
        guard schema == Self.schema, !sourceStage.isEmpty, !destinationStage.isEmpty,
              sourceStage.utf8.count <= 128, destinationStage.utf8.count <= 128,
              sourceStage != destinationStage, timePS.isFinite, timePS >= 0,
              preservesVelocities else {
            throw VivoArtifactValidationError.invalid("invalid velocity-preserving MD stage-transfer description")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    /// Explicit reconfiguration rather than restore with an identity bypass.
    /// Positions, periodic cell, physical time and accepted RNG step survive.
    /// Only the declared velocity initialization may change velocity state.
    public static func prepare(checkpoint: VivoMDCheckpoint,
                               source: VivoMDConfiguration, destination: VivoMDConfiguration,
                               particleCount: Int,
                               velocityInitialization: VivoMDVelocityInitialization = .preserve,
                               thermalizationSeed: UInt64? = nil) throws -> VivoMDStageTransition {
        try checkpoint.validate(particleCount: particleCount)
        let sourceFingerprint = try source.fingerprint()
        let destinationFingerprint = try destination.fingerprint()
        guard checkpoint.configurationFingerprint == sourceFingerprint else {
            throw VivoArtifactValidationError.incompatible("MD stage source configuration does not identify the checkpoint")
        }
        if velocityInitialization == .maxwellBoltzmann {
            guard thermalizationSeed != nil, destination.targetTemperatureK != nil else {
                throw VivoArtifactValidationError.invalid("Maxwell-Boltzmann initialization requires an explicit seed and target temperature")
            }
        } else if thermalizationSeed != nil {
            throw VivoArtifactValidationError.invalid("thermalization seed supplied without Maxwell-Boltzmann initialization")
        }
        var next = checkpoint
        next.configurationFingerprint = destinationFingerprint
        if velocityInitialization != .preserve {
            next.velocitiesNMPerPS = [VivoVector3D](repeating: .zero, count: particleCount)
        }
        try next.validate(particleCount: particleCount)
        return .init(schema: "numivivo.org/md-stage-transition/v1",
                     sourceCheckpoint: try checkpoint.fingerprint(),
                     sourceConfiguration: sourceFingerprint,
                     destinationConfiguration: destinationFingerprint,
                     velocityInitialization: velocityInitialization,
                     thermalizationSeed: thermalizationSeed,
                     destinationCheckpoint: next)
    }
}
