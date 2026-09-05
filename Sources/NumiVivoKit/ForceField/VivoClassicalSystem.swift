import Foundation

public enum VivoMixingRule: String, Codable, Sendable, CaseIterable {
    case lorentzBerthelot
    case geometric
    /// Pair coefficients are authoritative. Particle sigma/epsilon remain useful
    /// diagnostics, but execution must consult `nonbondedTypePairs`.
    case explicitPairTable
}

public enum VivoParticleRole: String, Codable, Sendable, CaseIterable {
    case atom
    case virtualSite
    case drude
}

public struct VivoClassicalParticle: Codable, Sendable, Equatable {
    public let index: UInt32
    public var atomIndex: UInt32?
    public var typeIdentifier: String
    public var role: VivoParticleRole
    /// Atomic/particle mass in daltons.
    public var massDa: Double
    /// Charge in elementary-charge units.
    public var chargeE: Double
    /// Lennard-Jones sigma in nanometres.
    public var sigmaNM: Double
    /// Lennard-Jones epsilon in kJ mol^-1.
    public var epsilonKJPerMol: Double

    public init(index: UInt32, atomIndex: UInt32?, typeIdentifier: String,
                role: VivoParticleRole = .atom, massDa: Double, chargeE: Double,
                sigmaNM: Double, epsilonKJPerMol: Double) {
        self.index = index; self.atomIndex = atomIndex; self.typeIdentifier = typeIdentifier
        self.role = role; self.massDa = massDa; self.chargeE = chargeE
        self.sigmaNM = sigmaNM; self.epsilonKJPerMol = epsilonKJPerMol
    }
}

public struct VivoHarmonicBond: Codable, Sendable, Equatable, Hashable {
    public var a: UInt32
    public var b: UInt32
    public var lengthNM: Double
    /// Potential is 1/2 k (r-r0)^2. Units kJ mol^-1 nm^-2.
    public var forceConstant: Double
    public init(a: UInt32, b: UInt32, lengthNM: Double, forceConstant: Double) {
        self.a = a; self.b = b; self.lengthNM = lengthNM; self.forceConstant = forceConstant
    }
}

public struct VivoHarmonicAngle: Codable, Sendable, Equatable, Hashable {
    public var a: UInt32
    public var b: UInt32
    public var c: UInt32
    public var angleRadians: Double
    /// Potential is 1/2 k (theta-theta0)^2. Units kJ mol^-1 rad^-2.
    public var forceConstant: Double
    public init(a: UInt32, b: UInt32, c: UInt32, angleRadians: Double, forceConstant: Double) {
        self.a = a; self.b = b; self.c = c; self.angleRadians = angleRadians; self.forceConstant = forceConstant
    }
}

public struct VivoPeriodicTorsion: Codable, Sendable, Equatable, Hashable {
    public var a: UInt32
    public var b: UInt32
    public var c: UInt32
    public var d: UInt32
    public var periodicity: UInt16
    public var phaseRadians: Double
    /// Potential is k * (1 + cos(n*phi-phase)). Units kJ mol^-1.
    public var barrierKJPerMol: Double
    public var improper: Bool

    public init(a: UInt32, b: UInt32, c: UInt32, d: UInt32,
                periodicity: UInt16, phaseRadians: Double,
                barrierKJPerMol: Double, improper: Bool = false) {
        self.a = a; self.b = b; self.c = c; self.d = d
        self.periodicity = periodicity; self.phaseRadians = phaseRadians
        self.barrierKJPerMol = barrierKJPerMol; self.improper = improper
    }
}

public struct VivoDistanceConstraint: Codable, Sendable, Equatable, Hashable {
    public var a: UInt32
    public var b: UInt32
    public var distanceNM: Double
    public init(a: UInt32, b: UInt32, distanceNM: Double) {
        self.a = a; self.b = b; self.distanceNM = distanceNM
    }
}

/// Direct 12-6 Lennard-Jones coefficients for a force-field type pair:
/// E(r) = C12/r^12 - C6/r^6. Units are kJ mol^-1 nm^12 and kJ mol^-1 nm^6.
/// This representation preserves precombined force-field tables such as AMBER's
/// NONBONDED_PARM_INDEX/A/B arrays without inferring a mixing rule.
public struct VivoNonbondedTypePair: Codable, Sendable, Equatable, Hashable {
    public var typeA: String
    public var typeB: String
    public var c12KJNM12PerMol: Double
    public var c6KJNM6PerMol: Double

    public init(typeA: String, typeB: String,
                c12KJNM12PerMol: Double, c6KJNM6PerMol: Double) {
        self.typeA = typeA; self.typeB = typeB
        self.c12KJNM12PerMol = c12KJNM12PerMol
        self.c6KJNM6PerMol = c6KJNM6PerMol
    }
}

/// Explicit pair policy after topology preprocessing. A scale of zero excludes
/// the corresponding interaction; 1 means normal nonbonded interaction.
public struct VivoNonbondedException: Codable, Sendable, Equatable, Hashable {
    public var a: UInt32
    public var b: UInt32
    public var coulombScale: Double
    public var lennardJonesScale: Double
    public var sigmaOverrideNM: Double?
    public var epsilonOverrideKJPerMol: Double?

    public init(a: UInt32, b: UInt32, coulombScale: Double, lennardJonesScale: Double,
                sigmaOverrideNM: Double? = nil, epsilonOverrideKJPerMol: Double? = nil) {
        self.a = a; self.b = b; self.coulombScale = coulombScale
        self.lennardJonesScale = lennardJonesScale
        self.sigmaOverrideNM = sigmaOverrideNM; self.epsilonOverrideKJPerMol = epsilonOverrideKJPerMol
    }
}

public struct VivoClassicalSystem: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/classical-system/v1"
    public var schema: String
    public var identifier: String
    public var structureFingerprint: VivoFingerprint
    public var parameterSourceFingerprints: [VivoFingerprint]
    public var mixingRule: VivoMixingRule
    public var particles: [VivoClassicalParticle]
    public var bonds: [VivoHarmonicBond]
    public var angles: [VivoHarmonicAngle]
    public var torsions: [VivoPeriodicTorsion]
    public var constraints: [VivoDistanceConstraint]
    /// Optional because v1 documents created before explicit-pair preservation
    /// decode with nil. Required when `mixingRule == .explicitPairTable`.
    public var nonbondedTypePairs: [VivoNonbondedTypePair]?
    public var nonbondedExceptions: [VivoNonbondedException]
    public var metadata: [String: String]

    public init(identifier: String, structureFingerprint: VivoFingerprint,
                parameterSourceFingerprints: [VivoFingerprint] = [],
                mixingRule: VivoMixingRule = .lorentzBerthelot,
                particles: [VivoClassicalParticle], bonds: [VivoHarmonicBond] = [],
                angles: [VivoHarmonicAngle] = [], torsions: [VivoPeriodicTorsion] = [],
                constraints: [VivoDistanceConstraint] = [],
                nonbondedTypePairs: [VivoNonbondedTypePair]? = nil,
                nonbondedExceptions: [VivoNonbondedException] = [],
                metadata: [String: String] = [:]) {
        self.schema = Self.schema; self.identifier = identifier
        self.structureFingerprint = structureFingerprint
        self.parameterSourceFingerprints = parameterSourceFingerprints
        self.mixingRule = mixingRule; self.particles = particles; self.bonds = bonds
        self.angles = angles; self.torsions = torsions; self.constraints = constraints
        self.nonbondedTypePairs = nonbondedTypePairs
        self.nonbondedExceptions = nonbondedExceptions; self.metadata = metadata
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoClassicalSystemValidator.validate(self)
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public enum VivoClassicalSystemValidator {
    public static func validate(_ system: VivoClassicalSystem, atomCount: UInt32? = nil) throws {
        guard system.schema == VivoClassicalSystem.schema,
              !system.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              !system.particles.isEmpty, system.particles.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("classical system schema, identifier or particle count is invalid")
        }
        var particleTypes = Set<String>()
        for (index, particle) in system.particles.enumerated() {
            guard particle.index == UInt32(index), !particle.typeIdentifier.isEmpty,
                  particle.massDa.isFinite, particle.massDa >= 0,
                  particle.chargeE.isFinite, particle.sigmaNM.isFinite, particle.sigmaNM >= 0,
                  particle.epsilonKJPerMol.isFinite, particle.epsilonKJPerMol >= 0 else {
                throw VivoArtifactValidationError.invalid("classical particle \(index) has invalid index or parameters")
            }
            particleTypes.insert(particle.typeIdentifier)
            if particle.role == .atom, particle.massDa <= 0 {
                throw VivoArtifactValidationError.invalid("ordinary atom particle \(index) requires positive mass")
            }
            if let atom = particle.atomIndex, let atomCount, atom >= atomCount {
                throw VivoArtifactValidationError.invalid("classical particle \(index) references absent structure atom \(atom)")
            }
        }
        let count = UInt32(system.particles.count)
        var bondedPairs = Set<UInt64>()
        for bond in system.bonds {
            try pair(bond.a, bond.b, count: count, label: "bond")
            guard bond.lengthNM.isFinite, bond.lengthNM > 0,
                  bond.forceConstant.isFinite, bond.forceConstant >= 0,
                  bondedPairs.insert(key(bond.a, bond.b)).inserted else {
                throw VivoArtifactValidationError.invalid("classical bond is invalid or duplicated")
            }
        }
        var angleKeys = Set<String>()
        for angle in system.angles {
            guard distinct([angle.a, angle.b, angle.c]), [angle.a, angle.b, angle.c].allSatisfy({ $0 < count }),
                  angle.angleRadians.isFinite, angle.angleRadians > 0, angle.angleRadians < Double.pi,
                  angle.forceConstant.isFinite, angle.forceConstant >= 0 else {
                throw VivoArtifactValidationError.invalid("classical angle is invalid")
            }
            let forward = "\(angle.a):\(angle.b):\(angle.c)", reverse = "\(angle.c):\(angle.b):\(angle.a)"
            guard angleKeys.insert(min(forward, reverse)).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate classical angle")
            }
        }
        for torsion in system.torsions {
            guard distinct([torsion.a, torsion.b, torsion.c, torsion.d]),
                  [torsion.a,torsion.b,torsion.c,torsion.d].allSatisfy({ $0 < count }),
                  torsion.periodicity > 0, torsion.phaseRadians.isFinite,
                  torsion.barrierKJPerMol.isFinite else {
                throw VivoArtifactValidationError.invalid("classical torsion has invalid indices or parameters")
            }
        }
        var constraintPairs = Set<UInt64>()
        for constraint in system.constraints {
            try pair(constraint.a, constraint.b, count: count, label: "constraint")
            guard constraint.distanceNM.isFinite, constraint.distanceNM > 0,
                  constraintPairs.insert(key(constraint.a, constraint.b)).inserted else {
                throw VivoArtifactValidationError.invalid("distance constraint is invalid or duplicated")
            }
        }
        if system.mixingRule == .explicitPairTable {
            guard let pairs = system.nonbondedTypePairs, !pairs.isEmpty else {
                throw VivoArtifactValidationError.invalid("explicitPairTable requires nonbonded type-pair coefficients")
            }
            var keys = Set<String>()
            for value in pairs {
                guard particleTypes.contains(value.typeA), particleTypes.contains(value.typeB),
                      value.c12KJNM12PerMol.isFinite, value.c12KJNM12PerMol >= 0,
                      value.c6KJNM6PerMol.isFinite, value.c6KJNM6PerMol >= 0 else {
                    throw VivoArtifactValidationError.invalid("nonbonded type-pair coefficient is invalid or references an unused type")
                }
                let pairKey = value.typeA <= value.typeB ? "\(value.typeA)\u{0}\(value.typeB)" : "\(value.typeB)\u{0}\(value.typeA)"
                guard keys.insert(pairKey).inserted else {
                    throw VivoArtifactValidationError.invalid("duplicate nonbonded type-pair coefficient")
                }
            }
            let required = particleTypes.count * (particleTypes.count + 1) / 2
            guard pairs.count == required else {
                throw VivoArtifactValidationError.invalid("explicit nonbonded table is incomplete: expected \(required) type pairs, found \(pairs.count)")
            }
        } else if let pairs = system.nonbondedTypePairs, !pairs.isEmpty {
            throw VivoArtifactValidationError.invalid("nonbonded type-pair coefficients require explicitPairTable mixing rule")
        }
        var exceptionPairs = Set<UInt64>()
        for exception in system.nonbondedExceptions {
            try pair(exception.a, exception.b, count: count, label: "nonbonded exception")
            guard exception.coulombScale.isFinite, exception.coulombScale >= 0,
                  exception.lennardJonesScale.isFinite, exception.lennardJonesScale >= 0,
                  exceptionPairs.insert(key(exception.a, exception.b)).inserted else {
                throw VivoArtifactValidationError.invalid("nonbonded exception is invalid or duplicated")
            }
            if let sigma = exception.sigmaOverrideNM, !sigma.isFinite || sigma < 0 {
                throw VivoArtifactValidationError.invalid("nonbonded sigma override is invalid")
            }
            if let epsilon = exception.epsilonOverrideKJPerMol, !epsilon.isFinite || epsilon < 0 {
                throw VivoArtifactValidationError.invalid("nonbonded epsilon override is invalid")
            }
        }
    }

    private static func pair(_ a: UInt32, _ b: UInt32, count: UInt32, label: String) throws {
        guard a != b, a < count, b < count else {
            throw VivoArtifactValidationError.invalid("\(label) particle pair is invalid")
        }
    }
    private static func distinct(_ values: [UInt32]) -> Bool { Set(values).count == values.count }
    private static func key(_ a: UInt32, _ b: UInt32) -> UInt64 { UInt64(min(a,b)) << 32 | UInt64(max(a,b)) }
}
