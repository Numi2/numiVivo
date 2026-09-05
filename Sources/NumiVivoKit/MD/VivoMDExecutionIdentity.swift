import Foundation

/// Increment when executable numerical semantics change. Source/configuration
/// identity alone does not identify the algorithm used to continue a trajectory.
public enum VivoMDExecutionIdentity {
    public static let current = "numivivo.org/md-metal-numerics/v2"
}

public enum VivoMDVelocityInitialization: String, Codable, Sendable {
    case preserve
    case zero
    case maxwellBoltzmann
}

/// Preparation record. For Maxwell-Boltzmann initialization the destination
/// checkpoint deliberately contains zero velocities until the runtime samples
/// and projects them. The protocol records the resulting entry checkpoint too.
public struct VivoMDStageTransition: Codable, Sendable, Equatable {
    public let schema: String
    public let sourceCheckpoint: VivoFingerprint
    public let sourceConfiguration: VivoFingerprint
    public let destinationConfiguration: VivoFingerprint
    public let velocityInitialization: VivoMDVelocityInitialization
    public let thermalizationSeed: UInt64?
    public let destinationCheckpoint: VivoMDCheckpoint
}
