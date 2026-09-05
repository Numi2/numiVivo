import Foundation

public enum VivoAlternateLocationPolicy: Codable, Sendable, Equatable {
    case highestOccupancy
    case preferred(String)
}

public struct VivoStructureRewrite: Codable, Sendable, Equatable {
    public var structure: VivoMolecularStructure
    /// Original atom index -> rewritten atom index; nil means removed.
    public var oldToNew: [UInt32?]
    /// Rewritten atom index -> original atom index.
    public var newToOld: [UInt32]

    public func validate(originalAtomCount: Int) throws {
        _ = try VivoStructureValidator.validate(structure)
        guard oldToNew.count == originalAtomCount, newToOld.count == structure.atoms.count else {
            throw VivoArtifactValidationError.invalid("structure rewrite mapping shape is inconsistent")
        }
        for (new, old) in newToOld.enumerated() {
            guard Int(old) < oldToNew.count, oldToNew[Int(old)] == UInt32(new) else {
                throw VivoArtifactValidationError.invalid("structure rewrite maps are not mutually consistent")
            }
        }
    }
}

public enum VivoAlternateLocationResolver {
    /// Resolves alternate PDB/mmCIF atom locations by semantic site. Selection is
    /// deterministic and applied identically to every conformer.
    public static func resolve(_ structure: VivoMolecularStructure,
                               policy: VivoAlternateLocationPolicy = .highestOccupancy) throws -> VivoStructureRewrite {
        _ = try VivoStructureValidator.validate(structure)
        let keys = try VivoAtomMapper.semanticKeys(structure)
        struct Site: Hashable {
            var chain: String?
            var sequence: Int32?
            var insertion: String?
            var residue: String?
            var atomName: String
            var element: UInt16
        }
        let sites = keys.map { Site(chain: $0.chainIdentifier, sequence: $0.residueSequenceNumber,
                                    insertion: $0.residueInsertionCode, residue: $0.residueName,
                                    atomName: $0.atomName, element: $0.elementAtomicNumber) }
        var groups: [Site: [UInt32]] = [:]
        for (index, site) in sites.enumerated() { groups[site, default: []].append(UInt32(index)) }

        var retained = Set<UInt32>()
        for group in groups.values {
            if group.count == 1 { retained.insert(group[0]); continue }
            let selected: UInt32
            switch policy {
            case .highestOccupancy:
                selected = best(group, structure: structure)
            case .preferred(let label):
                let exact = group.filter { structure.atoms[Int($0)].alternateLocation == label }
                if exact.count > 1 {
                    throw VivoArtifactValidationError.invalid("alternate-location label '\(label)' is duplicated within one atom site")
                }
                selected = exact.first ?? best(group, structure: structure)
            }
            retained.insert(selected)
        }
        return try VivoStructureSlicer.slice(structure, atomIndices: retained.sorted())
    }

    private static func best(_ group: [UInt32], structure: VivoMolecularStructure) -> UInt32 {
        group.sorted { lhs, rhs in
            let a = structure.atoms[Int(lhs)], b = structure.atoms[Int(rhs)]
            let ao = a.occupancy ?? 0, bo = b.occupancy ?? 0
            if ao != bo { return ao > bo }
            let ar = rank(a.alternateLocation), br = rank(b.alternateLocation)
            if ar != br { return ar < br }
            return lhs < rhs
        }[0]
    }

    private static func rank(_ label: String?) -> Int {
        guard let label else { return 0 }
        if label == "A" { return 1 }
        return 2 + Int(label.utf8.first ?? 255)
    }
}

public enum VivoStructureSlicer {
    /// Extracts an atom subset and rewrites all residue, chain, bond and conformer
    /// references. No hidden coordinates or stale indices are retained.
    public static func slice(_ source: VivoMolecularStructure,
                             atomIndices requested: [UInt32],
                             identifier: String? = nil) throws -> VivoStructureRewrite {
        _ = try VivoStructureValidator.validate(source)
        let unique = Array(Set(requested)).sorted()
        guard unique.allSatisfy({ Int($0) < source.atoms.count }) else {
            throw VivoArtifactValidationError.invalid("structure slice contains an out-of-range atom")
        }
        guard !unique.isEmpty else { throw VivoArtifactValidationError.invalid("structure slice cannot be empty") }
        var oldToNew = [UInt32?](repeating: nil, count: source.atoms.count)
        for (new, old) in unique.enumerated() { oldToNew[Int(old)] = UInt32(new) }

        let selectedSet = Set(unique)
        var keptResidues: [UInt32] = []
        for residue in source.residues where residue.atomIndices.contains(where: { selectedSet.contains($0) }) {
            keptResidues.append(residue.index)
        }
        let keptResidueSet = Set(keptResidues)
        var oldResidueToNew: [UInt32: UInt32] = [:]
        for (new, old) in keptResidues.enumerated() { oldResidueToNew[old] = UInt32(new) }

        var keptChains: [UInt32] = []
        for chain in source.chains where chain.residueIndices.contains(where: { keptResidueSet.contains($0) }) {
            keptChains.append(chain.index)
        }
        var oldChainToNew: [UInt32: UInt32] = [:]
        for (new, old) in keptChains.enumerated() { oldChainToNew[old] = UInt32(new) }

        var atoms: [VivoMolecularAtom] = []
        atoms.reserveCapacity(unique.count)
        for (new, old) in unique.enumerated() {
            var atom = source.atoms[Int(old)]
            atom = VivoMolecularAtom(index: UInt32(new), name: atom.name, element: atom.element,
                                     isotopeMassNumber: atom.isotopeMassNumber, formalCharge: atom.formalCharge,
                                     residueIndex: atom.residueIndex.flatMap { oldResidueToNew[$0] },
                                     sourceSerial: atom.sourceSerial, alternateLocation: atom.alternateLocation,
                                     occupancy: atom.occupancy, bFactor: atom.bFactor, isHetero: atom.isHetero)
            atoms.append(atom)
        }

        var residues: [VivoMolecularResidue] = []
        residues.reserveCapacity(keptResidues.count)
        for old in keptResidues {
            let sourceResidue = source.residues[Int(old)]
            let newIndex = oldResidueToNew[old]!
            residues.append(.init(index: newIndex, name: sourceResidue.name,
                                  chainIndex: sourceResidue.chainIndex.flatMap { oldChainToNew[$0] },
                                  sequenceNumber: sourceResidue.sequenceNumber,
                                  insertionCode: sourceResidue.insertionCode,
                                  atomIndices: sourceResidue.atomIndices.compactMap { oldToNew[Int($0)] }))
        }
        var chains: [VivoMolecularChain] = []
        chains.reserveCapacity(keptChains.count)
        for old in keptChains {
            let sourceChain = source.chains[Int(old)]
            chains.append(.init(index: oldChainToNew[old]!, identifier: sourceChain.identifier,
                                residueIndices: sourceChain.residueIndices.compactMap { oldResidueToNew[$0] }))
        }
        let bonds = source.bonds.compactMap { bond -> VivoMolecularBond? in
            guard let a = oldToNew[Int(bond.atomA)], let b = oldToNew[Int(bond.atomB)] else { return nil }
            return .init(atomA: a, atomB: b, order: bond.order, stereo: bond.stereo)
        }
        let conformers = source.conformers.map { conformer in
            VivoMolecularConformer(identifier: conformer.identifier,
                                   positionsNM: unique.map { conformer.positionsNM[Int($0)] },
                                   potentialEnergyKJPerMol: conformer.potentialEnergyKJPerMol)
        }
        let output = VivoMolecularStructure(identifier: identifier ?? source.identifier,
                                            atoms: atoms, bonds: bonds, residues: residues, chains: chains,
                                            conformers: conformers, periodicCell: source.periodicCell,
                                            metadata: source.metadata)
        var rewrite = VivoStructureRewrite(structure: output, oldToNew: oldToNew, newToOld: unique)
        try rewrite.validate(originalAtomCount: source.atoms.count)
        return rewrite
    }

    public static func connectedComponents(_ source: VivoMolecularStructure) throws -> [[UInt32]] {
        let index = try VivoStructureIndex(validated: source)
        var visited = [Bool](repeating: false, count: source.atoms.count)
        var result: [[UInt32]] = []
        for seed in source.atoms.indices where !visited[seed] {
            var component: [UInt32] = []
            var stack = [UInt32(seed)]
            visited[seed] = true
            while let atom = stack.popLast() {
                component.append(atom)
                for neighbor in index.bondedAtoms[Int(atom)] where !visited[Int(neighbor)] {
                    visited[Int(neighbor)] = true; stack.append(neighbor)
                }
            }
            result.append(component.sorted())
        }
        return result.sorted { lhs, rhs in lhs.count == rhs.count ? (lhs.first ?? 0) < (rhs.first ?? 0) : lhs.count > rhs.count }
    }
}
