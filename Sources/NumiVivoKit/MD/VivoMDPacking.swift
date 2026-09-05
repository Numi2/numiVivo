import Foundation

/// Immutable GPU-facing representation. The hot nonbonded path always consumes
/// a dense type-pair C12/C6 matrix, regardless of the source force-field mixing
/// rule. This keeps scientific semantics in the compiler rather than Metal.
public struct VivoMDPackedSystem: Sendable {
    public struct Bond: Sendable {
        public var atoms: SIMD2<UInt32>
        public var parameters: SIMD2<Float> // r0, k
    }
    public struct Angle: Sendable {
        public var atoms: SIMD4<UInt32> // a,b,c,unused
        public var parameters: SIMD2<Float> // theta0,k
    }
    public struct Torsion: Sendable {
        public var atoms: SIMD4<UInt32>
        public var parameters: SIMD4<Float> // periodicity,phase,k,improperFlag
    }
    public struct Constraint: Sendable {
        public var atoms: SIMD2<UInt32>
        public var distanceNM: Float
    }
    public struct PairException: Sendable {
        public var atoms: SIMD2<UInt32>
        public var scales: SIMD2<Float> // coulomb, LJ
        public var overrideC12C6: SIMD2<Float>
        public var flags: UInt32 // bit0 = explicit LJ override
    }
    public struct Incidence: Sendable {
        public var termIndex: UInt32
        public var localIndex: UInt32
    }

    public let systemFingerprint: VivoFingerprint
    public let particleCount: UInt32
    public let typeCount: UInt32
    /// x mass Da, y inverse mass Da^-1, z charge e, w reserved.
    public let particleDynamics: [SIMD4<Float>]
    public let particleTypeIndices: [UInt32]
    /// Row-major typeCount x typeCount. x=C12, y=C6.
    public let typePairC12C6: [SIMD2<Float>]
    public let bonds: [Bond]
    public let angles: [Angle]
    public let torsions: [Torsion]
    public let constraints: [Constraint]
    public let pairExceptions: [PairException]
    public let exceptionOffsets: [UInt32]
    public let exceptionPartners: [UInt32]
    public let exceptionIndices: [UInt32]
    public let bondOffsets: [UInt32]
    public let bondIncidence: [Incidence]
    public let angleOffsets: [UInt32]
    public let angleIncidence: [Incidence]
    public let torsionOffsets: [UInt32]
    public let torsionIncidence: [Incidence]
}

public enum VivoMDSystemPacker {
    public static func pack(_ system: VivoClassicalSystem) throws -> VivoMDPackedSystem {
        try VivoClassicalSystemValidator.validate(system)
        guard system.particles.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("MD particle count exceeds UInt32 execution ABI")
        }
        let fingerprint = try system.fingerprint()
        let typeNames = Array(Set(system.particles.map(\.typeIdentifier))).sorted()
        guard !typeNames.isEmpty, typeNames.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("MD force-field type table is empty or exceeds UInt32")
        }
        let typeIndex = Dictionary(uniqueKeysWithValues: typeNames.enumerated().map { ($0.element, UInt32($0.offset)) })

        var particleDynamics: [SIMD4<Float>] = []
        var particleTypes: [UInt32] = []
        particleDynamics.reserveCapacity(system.particles.count)
        particleTypes.reserveCapacity(system.particles.count)
        for particle in system.particles {
            let mass = try finiteFloat(particle.massDa, label: "particle mass")
            let inverseMass: Float = particle.massDa > 0 ? try finiteFloat(1 / particle.massDa, label: "inverse particle mass") : 0
            let charge = try finiteFloat(particle.chargeE, label: "particle charge")
            guard let type = typeIndex[particle.typeIdentifier] else {
                throw VivoArtifactValidationError.invalid("internal MD type-index construction failure")
            }
            particleDynamics.append(.init(mass, inverseMass, charge, 0))
            particleTypes.append(type)
        }

        let pairMatrix = try compilePairMatrix(system: system, typeNames: typeNames)
        let bonds = try system.bonds.map {
            VivoMDPackedSystem.Bond(atoms: .init($0.a, $0.b),
                                    parameters: .init(try finiteFloat($0.lengthNM, label: "bond length"),
                                                      try finiteFloat($0.forceConstant, label: "bond force constant")))
        }
        let angles = try system.angles.map {
            VivoMDPackedSystem.Angle(atoms: .init($0.a, $0.b, $0.c, 0),
                                     parameters: .init(try finiteFloat($0.angleRadians, label: "angle"),
                                                       try finiteFloat($0.forceConstant, label: "angle force constant")))
        }
        let torsions = try system.torsions.map {
            VivoMDPackedSystem.Torsion(atoms: .init($0.a, $0.b, $0.c, $0.d),
                                       parameters: .init(Float($0.periodicity),
                                                         try finiteFloat($0.phaseRadians, label: "torsion phase"),
                                                         try finiteFloat($0.barrierKJPerMol, label: "torsion barrier"),
                                                         $0.improper ? 1 : 0))
        }
        let constraints = try system.constraints.map {
            VivoMDPackedSystem.Constraint(atoms: .init($0.a, $0.b),
                                          distanceNM: try finiteFloat($0.distanceNM, label: "constraint distance"))
        }
        let pairExceptions = try system.nonbondedExceptions.map { value -> VivoMDPackedSystem.PairException in
            var override = SIMD2<Float>(repeating: 0)
            var flags: UInt32 = 0
            if value.sigmaOverrideNM != nil || value.epsilonOverrideKJPerMol != nil {
                guard let sigma = value.sigmaOverrideNM, let epsilon = value.epsilonOverrideKJPerMol else {
                    throw VivoArtifactValidationError.invalid("nonbonded exception must provide both sigma and epsilon overrides")
                }
                let coefficients = try ljFromSigmaEpsilon(sigma: sigma, epsilon: epsilon)
                override = coefficients
                flags |= 1
            }
            return .init(atoms: .init(value.a, value.b),
                         scales: .init(try finiteFloat(value.coulombScale, label: "exception Coulomb scale"),
                                       try finiteFloat(value.lennardJonesScale, label: "exception LJ scale")),
                         overrideC12C6: override, flags: flags)
        }.sorted { lhs, rhs in
            pairKey(lhs.atoms.x, lhs.atoms.y) < pairKey(rhs.atoms.x, rhs.atoms.y)
        }

        let exceptionAdjacency = try pairAdjacency(pairExceptions, particleCount: system.particles.count)
        let bondAdjacency = try incidence(termCount: bonds.count, particleCount: system.particles.count) { index in
            let atoms = bonds[index].atoms
            return [(atoms.x, 0), (atoms.y, 1)]
        }
        let angleAdjacency = try incidence(termCount: angles.count, particleCount: system.particles.count) { index in
            let atoms = angles[index].atoms
            return [(atoms.x, 0), (atoms.y, 1), (atoms.z, 2)]
        }
        let torsionAdjacency = try incidence(termCount: torsions.count, particleCount: system.particles.count) { index in
            let atoms = torsions[index].atoms
            return [(atoms.x, 0), (atoms.y, 1), (atoms.z, 2), (atoms.w, 3)]
        }

        return .init(systemFingerprint: fingerprint,
                     particleCount: UInt32(system.particles.count), typeCount: UInt32(typeNames.count),
                     particleDynamics: particleDynamics, particleTypeIndices: particleTypes,
                     typePairC12C6: pairMatrix, bonds: bonds, angles: angles, torsions: torsions,
                     constraints: constraints, pairExceptions: pairExceptions,
                     exceptionOffsets: exceptionAdjacency.offsets,
                     exceptionPartners: exceptionAdjacency.partners,
                     exceptionIndices: exceptionAdjacency.indices,
                     bondOffsets: bondAdjacency.offsets, bondIncidence: bondAdjacency.entries,
                     angleOffsets: angleAdjacency.offsets, angleIncidence: angleAdjacency.entries,
                     torsionOffsets: torsionAdjacency.offsets, torsionIncidence: torsionAdjacency.entries)
    }

    private static func compilePairMatrix(system: VivoClassicalSystem,
                                          typeNames: [String]) throws -> [SIMD2<Float>] {
        var byType: [String: VivoClassicalParticle] = [:]
        for particle in system.particles where byType[particle.typeIdentifier] == nil {
            byType[particle.typeIdentifier] = particle
        }
        let count = typeNames.count
        var output = [SIMD2<Float>](repeating: .zero, count: count * count)
        if system.mixingRule == .explicitPairTable {
            guard let pairs = system.nonbondedTypePairs else {
                throw VivoArtifactValidationError.invalid("explicit pair table disappeared after validation")
            }
            var lookup: [String: VivoNonbondedTypePair] = [:]
            for pair in pairs { lookup[typePairKey(pair.typeA, pair.typeB)] = pair }
            for i in 0..<count {
                for j in i..<count {
                    guard let pair = lookup[typePairKey(typeNames[i], typeNames[j])] else {
                        throw VivoArtifactValidationError.invalid("explicit pair table is incomplete during MD packing")
                    }
                    let value = SIMD2<Float>(try finiteFloat(pair.c12KJNM12PerMol, label: "C12"),
                                             try finiteFloat(pair.c6KJNM6PerMol, label: "C6"))
                    output[i * count + j] = value
                    output[j * count + i] = value
                }
            }
            return output
        }
        for i in 0..<count {
            guard let left = byType[typeNames[i]] else { throw VivoArtifactValidationError.invalid("missing representative particle type") }
            for j in i..<count {
                guard let right = byType[typeNames[j]] else { throw VivoArtifactValidationError.invalid("missing representative particle type") }
                let sigma: Double
                let epsilon: Double
                switch system.mixingRule {
                case .lorentzBerthelot:
                    sigma = 0.5 * (left.sigmaNM + right.sigmaNM)
                    epsilon = sqrt(left.epsilonKJPerMol * right.epsilonKJPerMol)
                case .geometric:
                    sigma = sqrt(left.sigmaNM * right.sigmaNM)
                    epsilon = sqrt(left.epsilonKJPerMol * right.epsilonKJPerMol)
                case .explicitPairTable:
                    preconditionFailure("handled above")
                }
                let value = try ljFromSigmaEpsilon(sigma: sigma, epsilon: epsilon)
                output[i * count + j] = value
                output[j * count + i] = value
            }
        }
        return output
    }

    private static func ljFromSigmaEpsilon(sigma: Double, epsilon: Double) throws -> SIMD2<Float> {
        guard sigma.isFinite, sigma >= 0, epsilon.isFinite, epsilon >= 0 else {
            throw VivoArtifactValidationError.invalid("LJ sigma/epsilon is invalid")
        }
        if sigma == 0 || epsilon == 0 { return .zero }
        let sigma6 = pow(sigma, 6)
        let c6 = 4 * epsilon * sigma6
        let c12 = c6 * sigma6
        return .init(try finiteFloat(c12, label: "derived C12"),
                     try finiteFloat(c6, label: "derived C6"))
    }

    private static func finiteFloat(_ value: Double, label: String) throws -> Float {
        guard value.isFinite, abs(value) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoArtifactValidationError.invalid("\(label) cannot be represented by the Wave B FP32 execution ABI")
        }
        let converted = Float(value)
        if value != 0, converted == 0 {
            throw VivoArtifactValidationError.invalid("\(label) underflows the Wave B FP32 execution ABI")
        }
        return converted
    }

    private static func typePairKey(_ a: String, _ b: String) -> String {
        a <= b ? a + "\u{0}" + b : b + "\u{0}" + a
    }
    private static func pairKey(_ a: UInt32, _ b: UInt32) -> UInt64 {
        UInt64(min(a,b)) << 32 | UInt64(max(a,b))
    }

    private static func incidence(termCount: Int, particleCount: Int,
                                  members: (Int) -> [(UInt32, UInt32)]) throws
    -> (offsets: [UInt32], entries: [VivoMDPackedSystem.Incidence]) {
        var buckets = [[VivoMDPackedSystem.Incidence]](repeating: [], count: particleCount)
        for term in 0..<termCount {
            guard let termIndex = UInt32(exactly: term) else {
                throw VivoArtifactValidationError.invalid("MD term count exceeds UInt32")
            }
            for (particle, local) in members(term) {
                guard particle < UInt32(particleCount) else {
                    throw VivoArtifactValidationError.invalid("MD incidence references absent particle")
                }
                buckets[Int(particle)].append(.init(termIndex: termIndex, localIndex: local))
            }
        }
        var offsets: [UInt32] = [0]
        var entries: [VivoMDPackedSystem.Incidence] = []
        for bucket in buckets {
            entries.append(contentsOf: bucket)
            guard let next = UInt32(exactly: entries.count) else {
                throw VivoArtifactValidationError.invalid("MD incidence table exceeds UInt32")
            }
            offsets.append(next)
        }
        return (offsets, entries)
    }

    private static func pairAdjacency(_ values: [VivoMDPackedSystem.PairException], particleCount: Int) throws
    -> (offsets: [UInt32], partners: [UInt32], indices: [UInt32]) {
        var buckets = [[(UInt32, UInt32)]](repeating: [], count: particleCount)
        for (index, value) in values.enumerated() {
            guard let tableIndex = UInt32(exactly: index), value.atoms.x < UInt32(particleCount), value.atoms.y < UInt32(particleCount) else {
                throw VivoArtifactValidationError.invalid("MD exception adjacency exceeds execution ABI")
            }
            buckets[Int(value.atoms.x)].append((value.atoms.y, tableIndex))
            buckets[Int(value.atoms.y)].append((value.atoms.x, tableIndex))
        }
        var offsets: [UInt32] = [0]
        var partners: [UInt32] = []
        var indices: [UInt32] = []
        for var bucket in buckets {
            bucket.sort { $0.0 < $1.0 }
            for item in bucket { partners.append(item.0); indices.append(item.1) }
            guard let next = UInt32(exactly: partners.count) else {
                throw VivoArtifactValidationError.invalid("MD exception adjacency exceeds UInt32")
            }
            offsets.append(next)
        }
        return (offsets, partners, indices)
    }
}
