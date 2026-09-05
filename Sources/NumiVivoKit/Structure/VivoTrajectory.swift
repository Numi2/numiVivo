import Foundation

public struct VivoTrajectoryFrame: Codable, Sendable, Equatable {
    public var step: UInt64
    public var timePS: Double
    public var positionsNM: [VivoVector3D]
    public var velocitiesNMPerPS: [VivoVector3D]?
    public var forcesKJPerMolPerNM: [VivoVector3D]?
    public var periodicCell: VivoPeriodicCell?
    public var potentialEnergyKJPerMol: Double?
    public var kineticEnergyKJPerMol: Double?
    public var temperatureK: Double?
    public var pressureBar: Double?

    public init(step: UInt64, timePS: Double, positionsNM: [VivoVector3D],
                velocitiesNMPerPS: [VivoVector3D]? = nil,
                forcesKJPerMolPerNM: [VivoVector3D]? = nil,
                periodicCell: VivoPeriodicCell? = nil,
                potentialEnergyKJPerMol: Double? = nil,
                kineticEnergyKJPerMol: Double? = nil,
                temperatureK: Double? = nil, pressureBar: Double? = nil) {
        self.step = step; self.timePS = timePS; self.positionsNM = positionsNM
        self.velocitiesNMPerPS = velocitiesNMPerPS; self.forcesKJPerMolPerNM = forcesKJPerMolPerNM
        self.periodicCell = periodicCell; self.potentialEnergyKJPerMol = potentialEnergyKJPerMol
        self.kineticEnergyKJPerMol = kineticEnergyKJPerMol
        self.temperatureK = temperatureK; self.pressureBar = pressureBar
    }

    public func validate(atomCount: Int) throws {
        guard timePS.isFinite, timePS >= 0, positionsNM.count == atomCount,
              positionsNM.allSatisfy(\.isFinite) else {
            throw VivoArtifactValidationError.invalid("trajectory frame has invalid time or positions")
        }
        if let velocitiesNMPerPS {
            guard velocitiesNMPerPS.count == atomCount, velocitiesNMPerPS.allSatisfy(\.isFinite) else {
                throw VivoArtifactValidationError.invalid("trajectory velocities must be finite and match atom count")
            }
        }
        if let forcesKJPerMolPerNM {
            guard forcesKJPerMolPerNM.count == atomCount, forcesKJPerMolPerNM.allSatisfy(\.isFinite) else {
                throw VivoArtifactValidationError.invalid("trajectory forces must be finite and match atom count")
            }
        }
        if let periodicCell, !periodicCell.isValid {
            throw VivoArtifactValidationError.invalid("trajectory periodic cell is invalid")
        }
        for (name, value) in [("potential energy", potentialEnergyKJPerMol),
                              ("kinetic energy", kineticEnergyKJPerMol),
                              ("temperature", temperatureK), ("pressure", pressureBar)] {
            if let value, !value.isFinite { throw VivoArtifactValidationError.invalid("trajectory \(name) must be finite") }
        }
        if let temperatureK, temperatureK < 0 {
            throw VivoArtifactValidationError.invalid("trajectory temperature cannot be negative")
        }
    }
}

public struct VivoTrajectory: Codable, Sendable, Equatable {
    public static let schemaVersion: UInt32 = 1
    public var schemaVersion: UInt32
    public var structureFingerprint: VivoFingerprint?
    public var atomCount: UInt32
    public var frames: [VivoTrajectoryFrame]
    public var metadata: [String: String]

    public init(structureFingerprint: VivoFingerprint? = nil, atomCount: UInt32,
                frames: [VivoTrajectoryFrame] = [], metadata: [String: String] = [:]) {
        self.schemaVersion = Self.schemaVersion
        self.structureFingerprint = structureFingerprint
        self.atomCount = atomCount; self.frames = frames; self.metadata = metadata
    }

    public func validate() throws {
        guard schemaVersion == Self.schemaVersion else {
            throw VivoArtifactValidationError.incompatible("unsupported trajectory schema \(schemaVersion)")
        }
        guard atomCount > 0 else { throw VivoArtifactValidationError.invalid("trajectory atom count must be positive") }
        var previousStep: UInt64?
        var previousTime: Double?
        for frame in frames {
            try frame.validate(atomCount: Int(atomCount))
            if let previousStep, frame.step <= previousStep {
                throw VivoArtifactValidationError.invalid("trajectory steps must increase strictly")
            }
            if let previousTime, frame.timePS <= previousTime {
                throw VivoArtifactValidationError.invalid("trajectory times must increase strictly")
            }
            previousStep = frame.step; previousTime = frame.timePS
        }
    }
}

public extension VivoPeriodicCell {
    /// Cartesian displacement mapped to the nearest image in a triclinic cell.
    func minimumImage(_ displacement: VivoVector3D) throws -> VivoVector3D {
        let det = a.dot(b.cross(c))
        guard det.isFinite, abs(det) > 1e-18 else {
            throw VivoArtifactValidationError.invalid("cannot apply minimum image to degenerate cell")
        }
        let fractional = VivoVector3D(displacement.dot(b.cross(c)) / det,
                                      displacement.dot(c.cross(a)) / det,
                                      displacement.dot(a.cross(b)) / det)
        let wrapped = VivoVector3D(fractional.x - fractional.x.rounded(),
                                   fractional.y - fractional.y.rounded(),
                                   fractional.z - fractional.z.rounded())
        return a * wrapped.x + b * wrapped.y + c * wrapped.z
    }
}
