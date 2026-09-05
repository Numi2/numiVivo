import Foundation

public struct VivoMDMinimizationConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: UInt32
    public var forceToleranceKJPerMolNM: Double
    public var initialStepScale: Double
    public var minimumStepScale: Double
    public var maximumStepScale: Double
    public var maximumDisplacementNM: Double
    public var acceptedStepGrowth: Double
    public var rejectedStepShrink: Double

    public init(maximumIterations: UInt32 = 10_000,
                forceToleranceKJPerMolNM: Double = 10,
                initialStepScale: Double = 1e-5,
                minimumStepScale: Double = 1e-12,
                maximumStepScale: Double = 1e-2,
                maximumDisplacementNM: Double = 0.01,
                acceptedStepGrowth: Double = 1.2,
                rejectedStepShrink: Double = 0.5) {
        self.maximumIterations = maximumIterations
        self.forceToleranceKJPerMolNM = forceToleranceKJPerMolNM
        self.initialStepScale = initialStepScale
        self.minimumStepScale = minimumStepScale
        self.maximumStepScale = maximumStepScale
        self.maximumDisplacementNM = maximumDisplacementNM
        self.acceptedStepGrowth = acceptedStepGrowth
        self.rejectedStepShrink = rejectedStepShrink
    }

    public func validate() throws {
        guard maximumIterations > 0,
              forceToleranceKJPerMolNM.isFinite, forceToleranceKJPerMolNM > 0,
              initialStepScale.isFinite, initialStepScale > 0,
              minimumStepScale.isFinite, minimumStepScale > 0,
              maximumStepScale.isFinite, maximumStepScale >= initialStepScale,
              minimumStepScale <= initialStepScale,
              maximumDisplacementNM.isFinite, maximumDisplacementNM > 0,
              acceptedStepGrowth.isFinite, acceptedStepGrowth > 1,
              rejectedStepShrink.isFinite, rejectedStepShrink > 0,
              rejectedStepShrink < 1 else {
            throw VivoArtifactValidationError.invalid("MD minimization configuration is invalid")
        }
    }
}

public struct VivoMDMinimizationCertificate: Codable, Sendable, Equatable {
    public var systemFingerprint: VivoFingerprint
    public var configurationFingerprint: VivoFingerprint
    public var converged: Bool
    public var attemptedIterations: UInt32
    public var acceptedIterations: UInt32
    public var rejectedIterations: UInt32
    public var initialPotentialEnergyKJPerMol: Double
    public var finalPotentialEnergyKJPerMol: Double
    public var finalMaximumForceKJPerMolNM: Double
    public var finalStepScale: Double
}

struct VivoMDMinimizationReduction {
    var potentialEnergyBits: UInt32 = Float(0).bitPattern
    var maximumForceBits: UInt32 = Float(0).bitPattern
    var nonFiniteCount: UInt32 = 0
    var reserved: UInt32 = 0
    mutating func clear() {
        potentialEnergyBits = Float(0).bitPattern
        maximumForceBits = Float(0).bitPattern
        nonFiniteCount = 0
        reserved = 0
    }
    var potentialEnergy: Float { Float(bitPattern: potentialEnergyBits) }
    var maximumForce: Float { Float(bitPattern: maximumForceBits) }
}
