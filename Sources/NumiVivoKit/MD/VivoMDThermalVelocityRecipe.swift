import Foundation

/// Physical parameters for one Maxwell-Boltzmann draw. Distinct from the
/// preserve/zero/maxwellBoltzmann stage-transition policy enum.
public struct VivoMDThermalVelocityRecipe: Codable, Sendable, Equatable {
    public let temperatureK: Double
    public let seed: UInt64
    public let removeCenterOfMass: Bool
    public init(temperatureK: Double, seed: UInt64, removeCenterOfMass: Bool = true) {
        self.temperatureK=temperatureK;self.seed=seed;self.removeCenterOfMass=removeCenterOfMass
    }
    public func validate() throws {
        guard temperatureK.isFinite,temperatureK>0 else {
            throw VivoArtifactValidationError.invalid("thermal velocity temperature must be finite and positive")
        }
    }
}
