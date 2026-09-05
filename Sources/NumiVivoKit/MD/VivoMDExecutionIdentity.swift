import Foundation

/// Increment whenever executable numerical semantics change. A source/config
/// hash alone does not identify the algorithm used to continue a trajectory.
public enum VivoMDExecutionIdentity {
    public static let current = "numivivo.org/md-metal-numerics/v2"
}

public enum VivoMDVelocityInitialization: String, Codable, Sendable {
    case preserve
    case zero
    case maxwellBoltzmann
}

public struct VivoMDStageTransition: Codable, Sendable, Equatable {
    public let schema: String
    public let sourceCheckpoint: VivoFingerprint
    public let sourceConfiguration: VivoFingerprint
    public let destinationConfiguration: VivoFingerprint
    public let velocityInitialization: VivoMDVelocityInitialization
    public let thermalizationSeed: UInt64?
    public let destinationCheckpoint: VivoMDCheckpoint
}

public enum VivoMDStageTransfer {
    /// Explicit configuration transition, not checkpoint restoration with a
    /// bypassed identity check. Positions, volume, physical time and accepted RNG
    /// step survive; velocity changes are declared and separately certified.
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
        return .init(schema: "numivivo.org/md-stage-transition/v1",
                     sourceCheckpoint: try checkpoint.fingerprint(), sourceConfiguration: sourceFingerprint,
                     destinationConfiguration: destinationFingerprint,
                     velocityInitialization: velocityInitialization, thermalizationSeed: thermalizationSeed,
                     destinationCheckpoint: next)
    }
}
