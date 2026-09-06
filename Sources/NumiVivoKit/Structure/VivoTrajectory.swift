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
    /// The same bounded reciprocal-lattice closest-image search used by velocity
    /// preparation. Fractional rounding is NOT sufficient for a skew cell.
    /// This also fixes existing topology, selection and solvent-promotion callers.
    /// It does not relax the Metal runtime's independently enforced cell limits.
    func minimumImage(_ displacement: VivoVector3D) throws -> VivoVector3D {
        try VivoMDPreparationGeometry.minimumImage(displacement, cell: self)
    }

    /// alpha=(b,c), beta=(a,c), gamma=(a,b); lengths in nm, angles in degrees.
    /// Produces a right-handed crystallographic cell with a along +x and b in xy.
    static func crystallographic(lengthsNM lengths: VivoVector3D,
                                anglesDegrees angles: VivoVector3D) throws -> Self {
        guard lengths.isFinite, angles.isFinite,
              lengths.x > 0, lengths.y > 0, lengths.z > 0,
              [angles.x, angles.y, angles.z].allSatisfy({ $0 > 0 && $0 < 180 }) else {
            throw VivoArtifactValidationError.invalid("invalid crystallographic lengths or angles")
        }
        let alpha = angles.x * .pi / 180, beta = angles.y * .pi / 180
        let gamma = angles.z * .pi / 180
        let cg = cos(gamma), sg = sin(gamma), cb = cos(beta)
        guard abs(sg) > 1e-10 else {
            throw VivoArtifactValidationError.invalid("degenerate crystallographic gamma")
        }
        let cy = (cos(alpha) - cb * cg) / sg
        let cz2 = 1 - cb * cb - cy * cy
        guard cz2.isFinite, cz2 > 1e-14 else {
            throw VivoArtifactValidationError.invalid("crystallographic angles do not define a nondegenerate cell")
        }
        let result = Self(a: .init(lengths.x, 0, 0),
                          b: .init(lengths.y * cg, lengths.y * sg, 0),
                          c: .init(lengths.z * cb, lengths.z * cy, lengths.z * sqrt(cz2)))
        guard result.isValid else { throw VivoArtifactValidationError.invalid("crystallographic cell overflow") }
        return result
    }

    /// Perpendicular distances between opposing cell faces, not vector lengths.
    func faceHeightsNM() throws -> VivoVector3D {
        guard isValid else { throw VivoArtifactValidationError.invalid("invalid cell for face heights") }
        let heights = VivoVector3D(volumeNM3 / b.cross(c).norm,
                                   volumeNM3 / c.cross(a).norm,
                                   volumeNM3 / a.cross(b).norm)
        guard heights.isFinite, min(heights.x, min(heights.y, heights.z)) > 0 else {
            throw VivoArtifactValidationError.invalid("cell face height overflow")
        }
        return heights
    }

    /// Sufficient (conservative for some primitive cells) single-image cutoff.
    /// Include the neighbor skin in radiusNM when checking a pair-list radius.
    func validateSingleImageRadius(_ radiusNM: Double) throws {
        let h = try faceHeightsNM()
        guard radiusNM.isFinite, radiusNM > 0,
              2 * radiusNM < min(h.x, min(h.y, h.z)) else {
            throw VivoArtifactValidationError.incompatible("pair radius must be below half every periodic face height")
        }
    }
}
