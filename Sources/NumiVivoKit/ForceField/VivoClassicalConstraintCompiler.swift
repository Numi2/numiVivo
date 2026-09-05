import Foundation

/// Deterministic lowering from chemical identity + equilibrium force-field geometry
/// into the holonomic distance constraints consumed by Wave B MD.
public enum VivoClassicalConstraintPolicy: String, Codable, Sendable, CaseIterable {
    case none
    /// Constrain every physical bond involving hydrogen. For O-H-H water residues,
    /// also constrain H-H from the equilibrium O-H lengths and H-O-H angle.
    case hydrogenBondsAndRigidWater
}

public struct VivoClassicalConstraintReport: Codable, Sendable, Equatable {
    public var policy: VivoClassicalConstraintPolicy
    public var existingConstraintCount: UInt32
    public var addedHydrogenBondCount: UInt32
    public var addedWaterHHCount: UInt32
    public var finalConstraintCount: UInt32
}

public enum VivoClassicalConstraintCompiler {
    public static func compile(system source: VivoClassicalSystem,
                               structure: VivoMolecularStructure,
                               policy: VivoClassicalConstraintPolicy) throws
    -> (system: VivoClassicalSystem, report: VivoClassicalConstraintReport) {
        try VivoStructureValidator.validate(structure)
        try VivoClassicalSystemValidator.validate(source, atomCount: UInt32(structure.atoms.count))
        guard source.structureFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(structure))) else {
            // Structure documents fingerprint the canonical structure payload. Refuse
            // constraint derivation against another molecular identity.
            throw VivoArtifactValidationError.incompatible("constraint compilation structure does not match classical-system structure fingerprint")
        }
        let existingCount = source.constraints.count
        guard policy != .none else {
            return (source, .init(policy: policy, existingConstraintCount: UInt32(existingCount),
                                  addedHydrogenBondCount: 0, addedWaterHHCount: 0,
                                  finalConstraintCount: UInt32(existingCount)))
        }

        let particleCount = source.particles.count
        var atomForParticle = [VivoMolecularAtom?](repeating: nil, count: particleCount)
        for particle in source.particles {
            if let atomIndex = particle.atomIndex {
                guard atomIndex < UInt32(structure.atoms.count) else {
                    throw VivoArtifactValidationError.invalid("particle-to-structure atom mapping is outside structure")
                }
                atomForParticle[Int(particle.index)] = structure.atoms[Int(atomIndex)]
            }
        }

        func key(_ a: UInt32, _ b: UInt32) -> UInt64 {
            UInt64(min(a,b)) << 32 | UInt64(max(a,b))
        }
        var constraints: [UInt64: VivoDistanceConstraint] = [:]
        for value in source.constraints {
            let k = key(value.a, value.b)
            if let previous = constraints[k], abs(previous.distanceNM - value.distanceNM) > 1e-10 {
                throw VivoArtifactValidationError.incompatible("conflicting pre-existing distance constraints")
            }
            constraints[k] = value
        }

        var addedHydrogen = 0
        var bondLength: [UInt64: Double] = [:]
        for bond in source.bonds {
            bondLength[key(bond.a, bond.b)] = bond.lengthNM
            guard let aa = atomForParticle[Int(bond.a)], let bb = atomForParticle[Int(bond.b)] else { continue }
            guard source.particles[Int(bond.a)].role == .atom,
                  source.particles[Int(bond.b)].role == .atom else { continue }
            if aa.element.atomicNumber == 1 || bb.element.atomicNumber == 1 {
                let k = key(bond.a, bond.b)
                if let previous = constraints[k] {
                    guard abs(previous.distanceNM - bond.lengthNM) <= 1e-8 else {
                        throw VivoArtifactValidationError.incompatible("hydrogen-bond constraint disagrees with equilibrium bond length")
                    }
                } else {
                    constraints[k] = .init(a: min(bond.a,bond.b), b: max(bond.a,bond.b), distanceNM: bond.lengthNM)
                    addedHydrogen += 1
                }
            }
        }

        // Rigid-water H-H closure. A qualifying residue contains exactly one O and
        // two H physical particles; extra virtual sites do not change that identity.
        var physicalByResidue: [UInt32: [UInt32]] = [:]
        for particle in source.particles where particle.role == .atom {
            guard let atom = atomForParticle[Int(particle.index)], let residue = atom.residueIndex else { continue }
            physicalByResidue[residue, default: []].append(particle.index)
        }
        var waterResidues = Set<UInt32>()
        for (residue, particles) in physicalByResidue where particles.count == 3 {
            let zs = particles.compactMap { atomForParticle[Int($0)]?.element.atomicNumber }.sorted()
            if zs == [1,1,8] { waterResidues.insert(residue) }
        }

        var addedWaterHH = 0
        for angle in source.angles {
            guard let a = atomForParticle[Int(angle.a)],
                  let b = atomForParticle[Int(angle.b)],
                  let c = atomForParticle[Int(angle.c)],
                  a.element.atomicNumber == 1,
                  b.element.atomicNumber == 8,
                  c.element.atomicNumber == 1,
                  let ra = a.residueIndex, let rb = b.residueIndex, let rc = c.residueIndex,
                  ra == rb, rb == rc, waterResidues.contains(rb) else { continue }
            guard let r1 = bondLength[key(angle.a, angle.b)],
                  let r2 = bondLength[key(angle.b, angle.c)] else {
                throw VivoArtifactValidationError.unresolved("rigid water angle is missing one or both O-H equilibrium bonds")
            }
            let hh2 = r1*r1 + r2*r2 - 2*r1*r2*cos(angle.angleRadians)
            guard hh2.isFinite, hh2 > 0 else {
                throw VivoArtifactValidationError.invalid("rigid water equilibrium geometry is degenerate")
            }
            let hh = sqrt(hh2), k = key(angle.a, angle.c)
            if let previous = constraints[k] {
                guard abs(previous.distanceNM - hh) <= 1e-8 else {
                    throw VivoArtifactValidationError.incompatible("water H-H constraint disagrees with equilibrium geometry")
                }
            } else {
                constraints[k] = .init(a: min(angle.a,angle.c), b: max(angle.a,angle.c), distanceNM: hh)
                addedWaterHH += 1
            }
        }

        var output = source
        output.constraints = constraints.values.sorted { lhs, rhs in key(lhs.a,lhs.b) < key(rhs.a,rhs.b) }
        output.metadata["constraintPolicy"] = policy.rawValue
        output.metadata["constraintCount"] = String(output.constraints.count)
        try VivoClassicalSystemValidator.validate(output, atomCount: UInt32(structure.atoms.count))
        guard output.constraints.count <= Int(UInt32.max), existingCount <= Int(UInt32.max),
              addedHydrogen <= Int(UInt32.max), addedWaterHH <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("constraint report exceeds UInt32")
        }
        return (output, .init(policy: policy,
                              existingConstraintCount: UInt32(existingCount),
                              addedHydrogenBondCount: UInt32(addedHydrogen),
                              addedWaterHHCount: UInt32(addedWaterHH),
                              finalConstraintCount: UInt32(output.constraints.count)))
    }
}
