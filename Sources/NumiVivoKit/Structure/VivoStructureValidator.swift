import Foundation

public struct VivoStructureValidationReport: Codable, Sendable, Equatable {
    public var atomCount: UInt32
    public var bondCount: UInt32
    public var residueCount: UInt32
    public var chainCount: UInt32
    public var conformerCount: UInt32
    public var connectedComponents: UInt32
    public var hasPeriodicCell: Bool

    public init(atomCount: UInt32, bondCount: UInt32, residueCount: UInt32,
                chainCount: UInt32, conformerCount: UInt32,
                connectedComponents: UInt32, hasPeriodicCell: Bool) {
        self.atomCount = atomCount; self.bondCount = bondCount
        self.residueCount = residueCount; self.chainCount = chainCount
        self.conformerCount = conformerCount; self.connectedComponents = connectedComponents
        self.hasPeriodicCell = hasPeriodicCell
    }
}

public enum VivoStructureValidator {
    @discardableResult
    public static func validate(_ structure: VivoMolecularStructure) throws -> VivoStructureValidationReport {
        guard structure.schemaVersion == VivoMolecularStructure.schemaVersion else {
            throw VivoArtifactValidationError.incompatible("unsupported molecular structure schema \(structure.schemaVersion)")
        }
        guard !structure.identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoArtifactValidationError.invalid("molecular structure identifier must not be empty")
        }
        guard structure.atoms.count <= Int(UInt32.max), structure.bonds.count <= Int(UInt32.max),
              structure.residues.count <= Int(UInt32.max), structure.chains.count <= Int(UInt32.max),
              structure.conformers.count <= Int(UInt32.max) else {
            throw VivoArtifactValidationError.invalid("molecular structure exceeds UInt32 indexing capacity")
        }

        for (index, atom) in structure.atoms.enumerated() {
            guard atom.index == UInt32(index) else {
                throw VivoArtifactValidationError.invalid("atom indices must be dense and equal to array order")
            }
            guard atom.element.atomicNumber > 0, atom.element.atomicNumber <= 118,
                  !atom.name.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
                throw VivoArtifactValidationError.invalid("atom \(index) has invalid element or empty name")
            }
            if let residue = atom.residueIndex, Int(residue) >= structure.residues.count {
                throw VivoArtifactValidationError.invalid("atom \(index) references absent residue \(residue)")
            }
            if let occupancy = atom.occupancy, !occupancy.isFinite || occupancy < 0 || occupancy > 1 {
                throw VivoArtifactValidationError.invalid("atom \(index) occupancy must be finite in [0,1]")
            }
            if let bFactor = atom.bFactor, !bFactor.isFinite || bFactor < 0 {
                throw VivoArtifactValidationError.invalid("atom \(index) B-factor must be finite and nonnegative")
            }
        }

        var seenBonds = Set<UInt64>()
        var adjacency = [[UInt32]](repeating: [], count: structure.atoms.count)
        for (index, bond) in structure.bonds.enumerated() {
            guard bond.atomA != bond.atomB,
                  Int(bond.atomA) < structure.atoms.count,
                  Int(bond.atomB) < structure.atoms.count else {
                throw VivoArtifactValidationError.invalid("bond \(index) has invalid atom endpoints")
            }
            let low = min(bond.atomA, bond.atomB), high = max(bond.atomA, bond.atomB)
            let key = UInt64(low) << 32 | UInt64(high)
            guard seenBonds.insert(key).inserted else {
                throw VivoArtifactValidationError.invalid("duplicate bond between atoms \(low) and \(high)")
            }
            adjacency[Int(low)].append(high)
            adjacency[Int(high)].append(low)
        }

        var residueMembership = [Int](repeating: 0, count: structure.atoms.count)
        for (index, residue) in structure.residues.enumerated() {
            guard residue.index == UInt32(index) else {
                throw VivoArtifactValidationError.invalid("residue indices must be dense and equal to array order")
            }
            if let chain = residue.chainIndex, Int(chain) >= structure.chains.count {
                throw VivoArtifactValidationError.invalid("residue \(index) references absent chain \(chain)")
            }
            var local = Set<UInt32>()
            for atom in residue.atomIndices {
                guard Int(atom) < structure.atoms.count, local.insert(atom).inserted else {
                    throw VivoArtifactValidationError.invalid("residue \(index) has invalid or duplicate atom membership")
                }
                residueMembership[Int(atom)] += 1
                guard structure.atoms[Int(atom)].residueIndex == residue.index else {
                    throw VivoArtifactValidationError.invalid("residue membership disagrees with atom \(atom) residueIndex")
                }
            }
        }
        guard residueMembership.allSatisfy({ $0 <= 1 }) else {
            throw VivoArtifactValidationError.invalid("an atom belongs to multiple residues")
        }
        for atom in structure.atoms where atom.residueIndex != nil {
            guard residueMembership[Int(atom.index)] == 1 else {
                throw VivoArtifactValidationError.invalid("atom \(atom.index) names a residue but is absent from its atom list")
            }
        }

        var residueChainMembership = [Int](repeating: 0, count: structure.residues.count)
        var chainIdentifiers = Set<String>()
        for (index, chain) in structure.chains.enumerated() {
            guard chain.index == UInt32(index), chainIdentifiers.insert(chain.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("chain indices must be dense and chain identifiers unique")
            }
            var local = Set<UInt32>()
            for residue in chain.residueIndices {
                guard Int(residue) < structure.residues.count, local.insert(residue).inserted else {
                    throw VivoArtifactValidationError.invalid("chain \(index) has invalid or duplicate residue membership")
                }
                residueChainMembership[Int(residue)] += 1
                guard structure.residues[Int(residue)].chainIndex == chain.index else {
                    throw VivoArtifactValidationError.invalid("chain membership disagrees with residue \(residue) chainIndex")
                }
            }
        }
        guard residueChainMembership.allSatisfy({ $0 <= 1 }) else {
            throw VivoArtifactValidationError.invalid("a residue belongs to multiple chains")
        }
        for residue in structure.residues where residue.chainIndex != nil {
            guard residueChainMembership[Int(residue.index)] == 1 else {
                throw VivoArtifactValidationError.invalid("residue \(residue.index) names a chain but is absent from its residue list")
            }
        }

        var conformerIDs = Set<String>()
        for conformer in structure.conformers {
            guard conformerIDs.insert(conformer.identifier).inserted,
                  conformer.positionsNM.count == structure.atoms.count,
                  conformer.positionsNM.allSatisfy(\.isFinite) else {
                throw VivoArtifactValidationError.invalid("conformer identifiers must be unique and each conformer must contain finite coordinates for every atom")
            }
            if let energy = conformer.potentialEnergyKJPerMol, !energy.isFinite {
                throw VivoArtifactValidationError.invalid("conformer energy must be finite when present")
            }
        }
        if let cell = structure.periodicCell, !cell.isValid {
            throw VivoArtifactValidationError.invalid("periodic cell must be finite and have positive nondegenerate volume")
        }

        var visited = [Bool](repeating: false, count: structure.atoms.count)
        var components: UInt32 = 0
        for seed in structure.atoms.indices where !visited[seed] {
            components += 1
            var stack = [seed]
            visited[seed] = true
            while let node = stack.popLast() {
                for neighbor in adjacency[node] {
                    let n = Int(neighbor)
                    if !visited[n] { visited[n] = true; stack.append(n) }
                }
            }
        }
        return .init(atomCount: UInt32(structure.atoms.count), bondCount: UInt32(structure.bonds.count),
                     residueCount: UInt32(structure.residues.count), chainCount: UInt32(structure.chains.count),
                     conformerCount: UInt32(structure.conformers.count), connectedComponents: components,
                     hasPeriodicCell: structure.periodicCell != nil)
    }
}

/// Immutable indexing facade. It is intentionally derived rather than serialized
/// so duplicated caches can never disagree with the authoritative structure.
public struct VivoStructureIndex: Sendable {
    public let bondedAtoms: [[UInt32]]
    public let atomsByElement: [UInt16: [UInt32]]
    public let atomsByResidue: [[UInt32]]
    public let residuesByChain: [[UInt32]]

    public init(validated structure: VivoMolecularStructure) throws {
        _ = try VivoStructureValidator.validate(structure)
        var adjacency = [[UInt32]](repeating: [], count: structure.atoms.count)
        for bond in structure.bonds {
            adjacency[Int(bond.atomA)].append(bond.atomB)
            adjacency[Int(bond.atomB)].append(bond.atomA)
        }
        for index in adjacency.indices { adjacency[index].sort() }
        bondedAtoms = adjacency
        atomsByElement = Dictionary(grouping: structure.atoms, by: { $0.element.atomicNumber })
            .mapValues { $0.map(\.index).sorted() }
        atomsByResidue = structure.residues.map { $0.atomIndices.sorted() }
        residuesByChain = structure.chains.map { $0.residueIndices.sorted() }
    }
}
