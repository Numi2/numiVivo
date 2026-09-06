import Foundation

public enum VivoSolventPromotionDistanceAtoms: String, Codable, Sendable {
    case allAtoms
    case heavyAtoms
}

public struct VivoQMMMSolventPromotionPolicy: Codable, Sendable, Equatable {
    public var radiusNM: Double
    public var solventResidueNames: [String]
    public var distanceAtoms: VivoSolventPromotionDistanceAtoms
    public var periodicMinimumImage: Bool
    public var maximumPromotedResidues: Int
    public var maximumPromotedAtoms: Int
    public init(radiusNM: Double,
                solventResidueNames: [String] = ["HOH","OPC","SOL","TIP3","TIP3P","WAT"],
                distanceAtoms: VivoSolventPromotionDistanceAtoms = .heavyAtoms,
                periodicMinimumImage: Bool = true,
                maximumPromotedResidues: Int = 256,
                maximumPromotedAtoms: Int = 2048) {
        self.radiusNM = radiusNM; self.solventResidueNames = solventResidueNames
        self.distanceAtoms = distanceAtoms; self.periodicMinimumImage = periodicMinimumImage
        self.maximumPromotedResidues = maximumPromotedResidues; self.maximumPromotedAtoms = maximumPromotedAtoms
    }
    public func validate() throws {
        let normalized = solventResidueNames.map { $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() }
        guard radiusNM.isFinite, radiusNM > 0, radiusNM <= 10, !normalized.isEmpty,
              normalized.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 64 }), Set(normalized).count == normalized.count,
              (1...100000).contains(maximumPromotedResidues), (1...1000000).contains(maximumPromotedAtoms) else {
            throw VivoChemistryError.invalid("solvent-shell radius, residue identities or capacity")
        }
    }
}

public struct VivoQMMMPromotedSolventResidue: Codable, Sendable, Equatable {
    public let residueIndex: UInt32
    public let residueName: String
    public let atomIndices: [UInt32]
    public let minimumDistanceNM: Double
    public let addedElectrons: Int
}

public struct VivoQMMMSolventPromotionResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-solvent-promotion/v1"
    public let schema: String
    public let policy: VivoQMMMSolventPromotionPolicy
    public let originalQMAtomIndices: [UInt32]
    public let promotedResidues: [VivoQMMMPromotedSolventResidue]
    public let qmAtomIndices: [UInt32]
    public let addedAlphaElectrons: Int
    public let addedBetaElectrons: Int
    public let periodicCellUsed: Bool
    public let requiresUnwrappedFiniteClusterForQMMM: Bool
    public let request: VivoQMMMRegionRequest
    public let interpretation: String
}

/// Standalone selection retains its explicit no-unwrapping contract. The MD
/// interface instead calls this selector on its reconstructed finite coordinates.
public enum VivoQMMMSolventPromotion {
    public static let interpretation = "whole-molecule explicit-solvent QM promotion with residue provenance; closed-shell electron assignment is a declared preparation assumption, not a spin-state determination; periodic selection does not itself unwrap coordinates; no solvent free-energy convergence implied"

    public static func promote(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                               particlePositionsNM positions: [VivoVector3D], base: VivoQMMMRegionRequest,
                               policy: VivoQMMMSolventPromotionPolicy) throws -> VivoQMMMSolventPromotionResult {
        try select(document: document, system: system, positions: positions, base: base, policy: policy,
                   distanceCell: policy.periodicMinimumImage ? document.structure.periodicCell : nil,
                   requiresUnwrapping: document.structure.periodicCell != nil)
    }

    static func select(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                       positions: [VivoVector3D], base: VivoQMMMRegionRequest,
                       policy: VivoQMMMSolventPromotionPolicy, distanceCell: VivoPeriodicCell?,
                       requiresUnwrapping: Bool, maximumDistancePairs: Int = 50_000_000) throws -> VivoQMMMSolventPromotionResult {
        try policy.validate()
        let topology = try VivoQMMMTopology(document: document, system: system)
        let structure = document.structure, n = structure.atoms.count
        guard positions.count == system.particles.count, positions.allSatisfy(\.isFinite), distanceCell?.isValid != false,
              !base.qmAtomIndices.isEmpty, Set(base.qmAtomIndices).count == base.qmAtomIndices.count,
              base.qmAtomIndices.allSatisfy({ Int($0) < n }),
              (0...100000).contains(base.alphaElectrons), (0...100000).contains(base.betaElectrons), maximumDistancePairs > 0 else {
            throw VivoChemistryError.invalid("solvent promotion frame or base QM region")
        }
        let original = Set(base.qmAtomIndices)
        let allowed = Set(policy.solventResidueNames.map { $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased() })
        let solventResidues = structure.residues.filter {
            allowed.contains($0.name.trimmingCharacters(in: .whitespacesAndNewlines).uppercased())
        }
        var residueGroups: [Int: [VivoMolecularResidue]] = [:]
        for residue in solventResidues {
            let groups = Set(residue.atomIndices.map { topology.atomMolecule[Int($0)] })
            guard groups.count == 1, let group = groups.first else {
                throw VivoChemistryError.invalid("solvent residue \(residue.index) lacks connected molecular topology")
            }
            residueGroups[group, default: []].append(residue)
        }
        var work = 0
        func minimumDistance(_ atoms: [UInt32]) throws -> Double {
            let distanceAtoms = atoms.filter {
                policy.distanceAtoms == .allAtoms || structure.atoms[Int($0)].element.atomicNumber > 1
            }
            guard !distanceAtoms.isEmpty else { throw VivoChemistryError.unsupported("heavy-atom solvent selection has no heavy atom") }
            var minimum = Double.infinity
            for atom in distanceAtoms { for q in base.qmAtomIndices.sorted() {
                guard work < maximumDistancePairs else { throw VivoChemistryError.resourceLimit("solvent selection pair-work budget") }
                work += 1
                let d = positions[Int(topology.atomToParticle[Int(atom)])] - positions[Int(topology.atomToParticle[Int(q)])]
                let nearest = try VivoMDPreparationGeometry.minimumImage(d, cell: distanceCell)
                minimum = min(minimum, nearest.norm)
            } }
            guard minimum.isFinite else { throw VivoChemistryError.invalid("solvent distance overflow") }
            return minimum
        }
        func electronCount(_ atoms: [UInt32]) -> Int {
            atoms.reduce(0) { $0 + Int(structure.atoms[Int($1)].element.atomicNumber) - Int(structure.atoms[Int($1)].formalCharge) }
        }
        struct Candidate {
            let group: Int
            let residues: [VivoMolecularResidue]
            let distance: Double
            let electrons: Int
        }
        var candidates: [Candidate] = []
        for group in residueGroups.keys.sorted() {
            let atoms = topology.molecules[group], residues = residueGroups[group]!
            let overlap = atoms.filter { original.contains($0) }
            guard overlap.isEmpty || overlap.count == atoms.count else {
                throw VivoChemistryError.invalid("base QM region partially selects a solvent molecule; complete molecular ownership is required")
            }
            if overlap.count == atoms.count { continue }
            guard Set(residues.flatMap(\.atomIndices)) == Set(atoms) else {
                throw VivoChemistryError.unsupported("solvent-labeled residue is connected to non-solvent atoms; select its complete molecule explicitly")
            }
            let minimum = try minimumDistance(atoms)
            if minimum <= policy.radiusNM {
                let electrons = electronCount(atoms)
                guard electrons >= 0, electrons % 2 == 0 else {
                    throw VivoChemistryError.unsupported("automatic solvent promotion requires a declared closed-shell molecule")
                }
                candidates.append(.init(group: group, residues: residues.sorted { $0.index < $1.index },
                                        distance: minimum, electrons: electrons))
            }
        }
        candidates.sort { $0.distance == $1.distance ? $0.group < $1.group : $0.distance < $1.distance }
        let atomCount = candidates.reduce(0) { $0 + topology.molecules[$1.group].count }
        let residueCount = candidates.reduce(0) { $0 + $1.residues.count }
        guard residueCount <= policy.maximumPromotedResidues, atomCount <= policy.maximumPromotedAtoms else {
            throw VivoChemistryError.resourceLimit("whole solvent shell exceeds explicit capacity; truncation is not permitted")
        }
        var selected = original, added = 0, records: [VivoQMMMPromotedSolventResidue] = []
        for item in candidates {
            selected.formUnion(topology.molecules[item.group]); added += item.electrons / 2
            for residue in item.residues {
                records.append(.init(residueIndex: residue.index, residueName: residue.name,
                    atomIndices: residue.atomIndices.sorted(), minimumDistanceNM: item.distance,
                    addedElectrons: electronCount(residue.atomIndices)))
            }
        }
        guard added <= 100000 - base.alphaElectrons, added <= 100000 - base.betaElectrons else {
            throw VivoChemistryError.resourceLimit("promoted electronic sector exceeds supported electron count")
        }
        var request = base
        request.qmAtomIndices = selected.sorted(); request.alphaElectrons += added; request.betaElectrons += added
        return .init(schema: VivoQMMMSolventPromotionResult.schema, policy: policy,
            originalQMAtomIndices: base.qmAtomIndices.sorted(), promotedResidues: records, qmAtomIndices: request.qmAtomIndices,
            addedAlphaElectrons: added, addedBetaElectrons: added, periodicCellUsed: distanceCell != nil,
            requiresUnwrappedFiniteClusterForQMMM: requiresUnwrapping, request: request, interpretation: interpretation)
    }
}

public extension VivoQMMMCompiler {
    static func promoteSolventRegion(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                     particlePositionsNM positions: [VivoVector3D], request: VivoQMMMRegionRequest,
                                     policy: VivoQMMMSolventPromotionPolicy) throws -> VivoQMMMSolventPromotionResult {
        try VivoQMMMSolventPromotion.promote(document: document, system: system, particlePositionsNM: positions, base: request, policy: policy)
    }
}
