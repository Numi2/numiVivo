import Foundation

public struct VivoForceFieldCompilationOptions: Codable, Sendable, Equatable {
    public var requireAllBondParameters: Bool
    public var requireAllAngleParameters: Bool
    public var requireAllProperTorsions: Bool

    public init(requireAllBondParameters: Bool = true,
                requireAllAngleParameters: Bool = true,
                requireAllProperTorsions: Bool = true) {
        self.requireAllBondParameters = requireAllBondParameters
        self.requireAllAngleParameters = requireAllAngleParameters
        self.requireAllProperTorsions = requireAllProperTorsions
    }
}

public struct VivoForceFieldCompilationReport: Codable, Sendable, Equatable {
    public var particleCount: UInt32
    public var bondCount: UInt32
    public var angleCount: UInt32
    public var torsionTermCount: UInt32
    public var exclusionCount: UInt32
    public var scaled14Count: UInt32
}

public struct VivoCompiledForceField: Codable, Sendable, Equatable {
    public var system: VivoClassicalSystem
    public var report: VivoForceFieldCompilationReport
}

public enum VivoForceFieldCompiler {
    public static func compile(structure: VivoMolecularStructure,
                               library: VivoForceFieldLibrary,
                               assignment suppliedAssignment: VivoForceFieldAssignment? = nil,
                               options: VivoForceFieldCompilationOptions = .init()) throws -> VivoCompiledForceField {
        _ = try VivoStructureValidator.validate(structure)
        try library.validate()
        let assignment = try suppliedAssignment ?? VivoResidueTemplateAssigner.assign(structure: structure, library: library)
        guard assignment.atomTypeByAtom.count == structure.atoms.count,
              assignment.chargeByAtomE.count == structure.atoms.count else {
            throw VivoArtifactValidationError.invalid("force-field assignment does not match structure atom count")
        }
        guard assignment.unresolvedAtoms.isEmpty else {
            throw VivoArtifactValidationError.unresolved("force-field assignment has unresolved atoms: \(assignment.unresolvedAtoms.prefix(32).map(String.init).joined(separator: ","))")
        }
        let atomTypes = Dictionary(uniqueKeysWithValues: library.atomTypes.map { ($0.identifier, $0) })
        var particles: [VivoClassicalParticle] = []
        particles.reserveCapacity(structure.atoms.count)
        for atom in structure.atoms {
            let typeID = assignment.atomTypeByAtom[Int(atom.index)]
            guard let type = atomTypes[typeID] else {
                throw VivoArtifactValidationError.unresolved("atom \(atom.index) references unknown force-field atom type '\(typeID)'")
            }
            if let expected = type.elementAtomicNumber, expected != atom.element.atomicNumber {
                throw VivoArtifactValidationError.incompatible("atom \(atom.index) element does not match force-field type '\(typeID)'")
            }
            particles.append(.init(index: atom.index, atomIndex: atom.index,
                                   typeIdentifier: typeID, massDa: type.massDa,
                                   chargeE: assignment.chargeByAtomE[Int(atom.index)],
                                   sigmaNM: type.sigmaNM, epsilonKJPerMol: type.epsilonKJPerMol))
        }

        let topology = try VivoStructureIndex(validated: structure)
        var bonds: [VivoHarmonicBond] = []
        bonds.reserveCapacity(structure.bonds.count)
        for edge in structure.bonds {
            let ta = particles[Int(edge.atomA)].typeIdentifier
            let tb = particles[Int(edge.atomB)].typeIdentifier
            if let parameter = bondParameter(ta, tb, library: library) {
                bonds.append(.init(a: edge.atomA, b: edge.atomB,
                                   lengthNM: parameter.lengthNM, forceConstant: parameter.forceConstant))
            } else if options.requireAllBondParameters {
                throw VivoArtifactValidationError.unresolved("missing bond parameter for \(ta)-\(tb) at atoms \(edge.atomA)-\(edge.atomB)")
            }
        }

        var angles: [VivoHarmonicAngle] = []
        for center in structure.atoms.indices {
            let neighbors = topology.bondedAtoms[center]
            if neighbors.count < 2 { continue }
            for i in neighbors.indices {
                for j in neighbors.indices where j > i {
                    let a = neighbors[i], b = UInt32(center), c = neighbors[j]
                    let types = (particles[Int(a)].typeIdentifier,
                                 particles[Int(b)].typeIdentifier,
                                 particles[Int(c)].typeIdentifier)
                    if let parameter = angleParameter(types.0, types.1, types.2, library: library) {
                        angles.append(.init(a: a, b: b, c: c,
                                            angleRadians: parameter.angleRadians,
                                            forceConstant: parameter.forceConstant))
                    } else if options.requireAllAngleParameters {
                        throw VivoArtifactValidationError.unresolved("missing angle parameter for \(types.0)-\(types.1)-\(types.2) at atoms \(a)-\(b)-\(c)")
                    }
                }
            }
        }

        var torsions: [VivoPeriodicTorsion] = []
        var seenProper = Set<String>()
        for central in structure.bonds {
            let b = central.atomA, c = central.atomB
            for a in topology.bondedAtoms[Int(b)] where a != c {
                for d in topology.bondedAtoms[Int(c)] where d != b && d != a {
                    let forward = "\(a):\(b):\(c):\(d)"
                    let reverse = "\(d):\(c):\(b):\(a)"
                    let canonical = min(forward, reverse)
                    guard seenProper.insert(canonical).inserted else { continue }
                    let types = [particles[Int(a)].typeIdentifier, particles[Int(b)].typeIdentifier,
                                 particles[Int(c)].typeIdentifier, particles[Int(d)].typeIdentifier]
                    let parameters = properTorsions(types, library: library)
                    if parameters.isEmpty, options.requireAllProperTorsions {
                        throw VivoArtifactValidationError.unresolved("missing proper torsion parameters for \(types.joined(separator: "-")) at atoms \(canonical)")
                    }
                    for parameter in parameters {
                        torsions.append(.init(a: a, b: b, c: c, d: d,
                                              periodicity: parameter.periodicity,
                                              phaseRadians: parameter.phaseRadians,
                                              barrierKJPerMol: parameter.barrierKJPerMol))
                    }
                }
            }
        }

        let pairPolicy = pairExceptions(topology: topology, library: library)
        let structureFingerprint = try VivoStructureCodec.fingerprint(structure)
        let libraryFingerprint = try library.fingerprint()
        let system = VivoClassicalSystem(identifier: "\(structure.identifier):\(library.identifier)",
                                         structureFingerprint: structureFingerprint,
                                         parameterSourceFingerprints: [libraryFingerprint],
                                         mixingRule: library.mixingRule,
                                         particles: particles, bonds: bonds, angles: angles,
                                         torsions: torsions, constraints: [],
                                         nonbondedExceptions: pairPolicy.exceptions,
                                         metadata: ["forceField": library.identifier,
                                                    "forceFieldVersion": library.version])
        try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(structure.atoms.count))
        return .init(system: system,
                     report: .init(particleCount: UInt32(particles.count), bondCount: UInt32(bonds.count),
                                   angleCount: UInt32(angles.count), torsionTermCount: UInt32(torsions.count),
                                   exclusionCount: pairPolicy.exclusions, scaled14Count: pairPolicy.scaled14))
    }

    private static func bondParameter(_ a: String, _ b: String,
                                      library: VivoForceFieldLibrary) -> VivoForceFieldBondParameter? {
        library.bondParameters.first { ($0.typeA == a && $0.typeB == b) || ($0.typeA == b && $0.typeB == a) }
    }

    private static func angleParameter(_ a: String, _ b: String, _ c: String,
                                       library: VivoForceFieldLibrary) -> VivoForceFieldAngleParameter? {
        library.angleParameters.first {
            ($0.typeA == a && $0.typeB == b && $0.typeC == c) ||
            ($0.typeA == c && $0.typeB == b && $0.typeC == a)
        }
    }

    private static func properTorsions(_ types: [String],
                                       library: VivoForceFieldLibrary) -> [VivoForceFieldTorsionParameter] {
        let proper = library.torsionParameters.filter { !$0.improper }
        let exact = proper.filter {
            [$0.typeA,$0.typeB,$0.typeC,$0.typeD] == types ||
            [$0.typeD,$0.typeC,$0.typeB,$0.typeA] == types
        }
        if !exact.isEmpty { return exact }
        return proper.filter { parameter in
            let forward = (parameter.typeA == "*" || parameter.typeA == types[0]) &&
                          parameter.typeB == types[1] && parameter.typeC == types[2] &&
                          (parameter.typeD == "*" || parameter.typeD == types[3])
            let reverse = (parameter.typeD == "*" || parameter.typeD == types[0]) &&
                          parameter.typeC == types[1] && parameter.typeB == types[2] &&
                          (parameter.typeA == "*" || parameter.typeA == types[3])
            return forward || reverse
        }
    }

    private static func pairExceptions(topology: VivoStructureIndex,
                                       library: VivoForceFieldLibrary) -> (exceptions: [VivoNonbondedException], exclusions: UInt32, scaled14: UInt32) {
        var shortest: [UInt64: Int] = [:]
        for origin in topology.bondedAtoms.indices {
            var frontier: [(UInt32, Int)] = topology.bondedAtoms[origin].map { ($0, 1) }
            var visited: [UInt32: Int] = [UInt32(origin): 0]
            var cursor = 0
            while cursor < frontier.count {
                let (atom, depth) = frontier[cursor]; cursor += 1
                guard depth <= 3 else { continue }
                if let old = visited[atom], old <= depth { continue }
                visited[atom] = depth
                let a = UInt32(origin), b = atom
                if a != b {
                    let key = pairKey(a,b)
                    shortest[key] = min(shortest[key] ?? Int.max, depth)
                }
                if depth < 3 {
                    for neighbor in topology.bondedAtoms[Int(atom)] { frontier.append((neighbor, depth + 1)) }
                }
            }
        }
        var result: [VivoNonbondedException] = []
        var exclusions: UInt32 = 0, scaled14: UInt32 = 0
        for (encoded, depth) in shortest.sorted(by: { $0.key < $1.key }) {
            let a = UInt32(encoded >> 32), b = UInt32(encoded & 0xffff_ffff)
            if depth <= 2 {
                result.append(.init(a: a, b: b, coulombScale: 0, lennardJonesScale: 0)); exclusions += 1
            } else if depth == 3 {
                result.append(.init(a: a, b: b, coulombScale: library.coulomb14Scale,
                                    lennardJonesScale: library.lennardJones14Scale)); scaled14 += 1
            }
        }
        return (result, exclusions, scaled14)
    }

    private static func pairKey(_ a: UInt32, _ b: UInt32) -> UInt64 {
        UInt64(min(a,b)) << 32 | UInt64(max(a,b))
    }
}
