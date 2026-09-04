import Foundation

/// Configuration-bound, accepted-boundary molecular checkpoint. The older v2
/// state payload alone has no spatial-grid descriptor and cannot establish that
/// two equally-sized spatial layouts describe the same simulation.
public struct VivoMolecularResumeCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/molecular-resume-checkpoint/v1"

    public let schema: String
    public let configuration: VivoRuntimeConfiguration
    public let state: VivoMolecularCheckpoint

    public init(configuration: VivoRuntimeConfiguration, state: VivoMolecularCheckpoint) throws {
        self.schema = Self.schema
        self.configuration = configuration
        self.state = state
        try validate()
    }

    public func validate() throws {
        try state.validate()
        if let grid = configuration.spatialGrid { try grid.validate() }
        guard schema == Self.schema,
              configuration.environmentCount > 0,
              configuration.laneCount > 0,
              configuration.laneCount == state.laneCount,
              configuration.environmentCount == state.parameterEnvironmentCount,
              configuration.fidelity == state.fidelity,
              configuration.seed == state.seed,
              configuration.timeStep.isFinite,
              configuration.minimumTimeStep.isFinite,
              configuration.maximumTimeStep.isFinite,
              configuration.minimumTimeStep > 0,
              configuration.timeStep >= configuration.minimumTimeStep,
              configuration.timeStep <= configuration.maximumTimeStep,
              configuration.maximumSubsteps > 0,
              configuration.eventCapacity > 0,
              configuration.privateHeapHeadroom.isFinite,
              configuration.privateHeapHeadroom >= 1 else {
            throw VivoArtifactValidationError.invalid(
                "molecular resume checkpoint configuration is invalid or disagrees with its state"
            )
        }
        if configuration.fidelity.rawValue >= VivoFidelity.spatial.rawValue,
           configuration.spatialGrid == nil {
            throw VivoArtifactValidationError.invalid(
                "spatial molecular resume checkpoints require an explicit grid"
            )
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}
