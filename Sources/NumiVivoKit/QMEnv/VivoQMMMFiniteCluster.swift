import Foundation

/// Preparation-time limits. These do not alter the MD integrator or its state.
public struct VivoQMMMClusterConfiguration: Codable, Sendable, Equatable {
    public var cycleToleranceNM: Double
    public var virtualSiteToleranceNM: Double
    public var maximumImagePairs: Int
    public init(cycleToleranceNM: Double = 1e-6, virtualSiteToleranceNM: Double = 1e-5,
                maximumImagePairs: Int = 50_000_000) {
        self.cycleToleranceNM = cycleToleranceNM
        self.virtualSiteToleranceNM = virtualSiteToleranceNM
        self.maximumImagePairs = maximumImagePairs
    }
    public func validate() throws {
        guard cycleToleranceNM.isFinite, cycleToleranceNM > 0, cycleToleranceNM < 0.01,
              virtualSiteToleranceNM.isFinite, virtualSiteToleranceNM > 0, virtualSiteToleranceNM < 0.01,
              maximumImagePairs > 0 else {
            throw VivoChemistryError.invalid("QM/MM cluster tolerances or image-work capacity")
        }
    }
}

public struct VivoQMMMParticleImage: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let structureAtomIndex: UInt32?
    public let moleculeIndex: Int
    /// Integer coefficients of the lattice translation ADDED to a physical
    /// particle's source position. Virtual sites are reconstructed, not shifted.
    public let latticeImage: [Int64]?
}

public struct VivoQMMMFiniteCluster: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-finite-cluster/v1"
    public let schema: String
    public let structureFingerprint: VivoFingerprint
    public let classicalSystemFingerprint: VivoFingerprint
    public let sourceFrameFingerprint: VivoFingerprint
    /// The sampled cell, including NPT changes; never substituted by the PDB cell.
    public let sourcePeriodicCell: VivoPeriodicCell?
    public let configuration: VivoQMMMClusterConfiguration
    public let anchorAtomIndex: UInt32
    public let atomToParticle: [UInt32]
    /// All source particle slots are retained, including virtual sites and Z1 atoms.
    public let particlePositionsNM: [VivoVector3D]
    public let particleImages: [VivoQMMMParticleImage]
    public let moleculeAtomIndices: [[UInt32]]
    public let linearVirtualSites: [VivoLinearVirtualSite]
    public let virtualSiteDefinitions: [VivoDependentSite]?
    public let imagePairEvaluations: Int

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

/// A single topology authority for molecular closure and periodic reconstruction.
/// A residue label is NOT a molecular connectivity definition.
struct VivoQMMMTopology {
    let atomToParticle: [UInt32]
    let particleToAtom: [UInt32?]
    let adjacency: [[Int]]
    let molecules: [[UInt32]]
    let atomMolecule: [Int]
    let sites: [VivoLinearVirtualSite]
    let siteGraph: VivoVirtualSiteGraph

    init(document: VivoMolecularStructureDocument, system: VivoClassicalSystem) throws {
        try document.validate()
        let n = document.structure.atoms.count
        guard n > 0, n <= Int(UInt32.max) else { throw VivoChemistryError.invalid("QM/MM topology atom count") }
        try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(n))
        guard document.structureFingerprint == system.structureFingerprint else {
            throw VivoChemistryError.invalid("QM/MM structure and classical system identities differ")
        }
        var map = [Int](repeating: -1, count: n)
        var reverse = [UInt32?](repeating: nil, count: system.particles.count)
        for p in system.particles {
            guard p.role != .drude else { throw VivoChemistryError.unsupported("Drude QM/MM reconstruction") }
            if p.role == .atom {
                guard let a = p.atomIndex, Int(a) < n, map[Int(a)] == -1 else {
                    throw VivoChemistryError.invalid("QM/MM requires a physical atom/particle bijection")
                }
                map[Int(a)] = Int(p.index); reverse[Int(p.index)] = a
            } else if p.atomIndex != nil {
                throw VivoChemistryError.invalid("a virtual site cannot own a chemical atom")
            }
        }
        guard map.allSatisfy({ $0 >= 0 }) else { throw VivoChemistryError.invalid("QM/MM unmapped physical atom") }
        let rules = (system.linearVirtualSites ?? []).sorted { $0.siteParticle < $1.siteParticle }
        let siteGraph = try system.resolvedVirtualSiteGraph()
        var edges = [Set<Int>](repeating: [], count: n)
        func joinAtoms(_ a: Int, _ b: Int) {
            edges[a].insert(b); edges[b].insert(a)
        }
        func joinParticles(_ a: UInt32, _ b: UInt32) throws {
            guard let aa = reverse[Int(a)], let bb = reverse[Int(b)] else {
                throw VivoChemistryError.unsupported("QM/MM connectivity must join physical particles")
            }
            joinAtoms(Int(aa), Int(bb))
        }
        for bond in document.structure.bonds { joinAtoms(Int(bond.atomA), Int(bond.atomB)) }
        for bond in system.bonds { try joinParticles(bond.a, bond.b) }
        for constraint in system.constraints { try joinParticles(constraint.a, constraint.b) }
        // Construction dependencies keep all parents and their site in one image.
        // They do not invent chemical bonds or authorize a link-atom cut.
        for parents in siteGraph.physicalAncestors {
            for parent in parents.dropFirst() { try joinParticles(parents[0], parent) }
        }
        let graph = edges.map { $0.sorted() }
        var membership = [Int](repeating: -1, count: n), groups: [[UInt32]] = []
        for root in 0..<n where membership[root] == -1 {
            let group = groups.count
            var queue = [root], cursor = 0
            membership[root] = group
            while cursor < queue.count {
                let atom = queue[cursor]; cursor += 1
                for neighbor in graph[atom] where membership[neighbor] == -1 {
                    membership[neighbor] = group; queue.append(neighbor)
                }
            }
            groups.append(queue.sorted().map(UInt32.init))
        }
        atomToParticle = map.map(UInt32.init); particleToAtom = reverse
        adjacency = graph; molecules = groups; atomMolecule = membership; sites = rules; self.siteGraph = siteGraph
    }

    func rebuildVirtualSites(_ positions: inout [VivoVector3D]) throws {
        positions = try siteGraph.construct(positionsNM: positions).positionsNM
    }
}

public enum VivoQMMMClusterBuilder {
    /// Reconstruct all physical molecules before placing whole images around the
    /// base QM region. No coordinates are cropped, and no source index is renumbered.
    public static func reconstruct(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                   particlePositionsNM source: [VivoVector3D], periodicCell cell: VivoPeriodicCell?,
                                   sourceFrameFingerprint: VivoFingerprint, qmAtomIndices: [UInt32],
                                   solventPolicy: VivoQMMMSolventPromotionPolicy? = nil,
                                   configuration cfg: VivoQMMMClusterConfiguration = .init()) throws -> VivoQMMMFiniteCluster {
        try cfg.validate(); try solventPolicy?.validate()
        let topology = try VivoQMMMTopology(document: document, system: system)
        guard source.count == system.particles.count, source.allSatisfy(\.isFinite), cell?.isValid != false,
              !qmAtomIndices.isEmpty, Set(qmAtomIndices).count == qmAtomIndices.count,
              qmAtomIndices.allSatisfy({ Int($0) < topology.atomToParticle.count }) else {
            throw VivoChemistryError.invalid("QM/MM source frame, cell or base selection")
        }
        if cell != nil, solventPolicy?.periodicMinimumImage == false {
            throw VivoChemistryError.invalid("periodic MD extraction requires periodic solvent imaging; nonperiodic selection is available on a finite frame")
        }
        let seeds = qmAtomIndices.sorted(), anchor = seeds[0]
        let anchorGroup = topology.atomMolecule[Int(anchor)]
        let seedGroups = Set(seeds.map { topology.atomMolecule[Int($0)] })
        var positions = source
        var visited = [Bool](repeating: false, count: topology.atomToParticle.count)
        var work = 0
        func image(_ d: VivoVector3D) throws -> VivoVector3D {
            guard work < cfg.maximumImagePairs else { throw VivoChemistryError.resourceLimit("QM/MM periodic image-work budget") }
            work += 1
            return try VivoMDPreparationGeometry.minimumImage(d, cell: cell)
        }
        func p(_ atom: UInt32) -> Int { Int(topology.atomToParticle[Int(atom)]) }
        // Bond/constraint/dependency traversal is performed in chemical atom order,
        // not particle order; a permuted MD particle table remains unambiguous.
        for (group, atoms) in topology.molecules.enumerated() {
            let root = group == anchorGroup ? anchor : atoms[0]
            var queue = [Int(root)], cursor = 0
            visited[Int(root)] = true
            while cursor < queue.count {
                let a = queue[cursor]; cursor += 1
                for b in topology.adjacency[a] {
                    let displacement = try image(source[p(UInt32(b))] - source[p(UInt32(a))])
                    let candidate = positions[p(UInt32(a))] + displacement
                    if visited[b] {
                        guard (candidate - positions[p(UInt32(b))]).norm <= cfg.cycleToleranceNM else {
                            throw VivoChemistryError.unsupported("molecule has inconsistent periodic cycle/winding; no unique finite molecular image")
                        }
                    } else {
                        positions[p(UInt32(b))] = candidate
                        visited[b] = true; queue.append(b)
                    }
                }
            }
        }
        // Validate sampled virtual sites modulo the current cell, using unwrapped
        // parents. Wrapped-parent averaging is never used as a construction rule.
        let siteState = try topology.siteGraph.construct(positionsNM: positions)
        for site in topology.siteGraph.sites {
            let residual = try image(source[Int(site.siteParticle)] - siteState.positionsNM[Int(site.siteParticle)])
            guard residual.norm <= cfg.virtualSiteToleranceNM else {
                throw VivoChemistryError.invalid("stale or inconsistent sampled virtual site \(site.siteParticle)")
            }
        }
        func translate(_ group: Int, by shift: VivoVector3D) {
            for atom in topology.molecules[group] { positions[p(atom)] = positions[p(atom)] + shift }
        }
        // Find a single molecule translation, never independent atom translations.
        func nearestShift(atoms: [UInt32], references: [UInt32]) throws -> VivoVector3D {
            var best = Double.infinity, shift = VivoVector3D.zero
            for a in atoms.sorted() { for q in references.sorted() {
                let d = positions[p(a)] - positions[p(q)]
                let nearest = try image(d), distance = nearest.squaredNorm
                guard distance.isFinite else { throw VivoChemistryError.invalid("QM/MM image distance overflow") }
                // Sorted identity order provides deterministic equal-distance ties.
                if distance < best { best = distance; shift = nearest - d }
            } }
            guard best.isFinite else { throw VivoChemistryError.invalid("empty molecule imaging selection") }
            return shift
        }
        if cell != nil {
            var placed = seeds.filter { topology.atomMolecule[Int($0)] == anchorGroup }
            for group in seedGroups.sorted() where group != anchorGroup {
                let localSeeds = seeds.filter { topology.atomMolecule[Int($0)] == group }
                translate(group, by: try nearestShift(atoms: localSeeds, references: placed))
                placed += localSeeds
            }
            let names = Set((solventPolicy?.solventResidueNames ?? []).map {
                $0.trimmingCharacters(in: .whitespacesAndNewlines).uppercased()
            })
            for (group, atoms) in topology.molecules.enumerated() where !seedGroups.contains(group) {
                let solvent = atoms.allSatisfy { atom in
                    guard let r = document.structure.atoms[Int(atom)].residueIndex else { return false }
                    return names.contains(document.structure.residues[Int(r)].name.trimmingCharacters(in: .whitespacesAndNewlines).uppercased())
                }
                let distanceAtoms = solvent && solventPolicy?.distanceAtoms == .heavyAtoms
                    ? atoms.filter { document.structure.atoms[Int($0)].element.atomicNumber > 1 } : atoms
                translate(group, by: try nearestShift(atoms: distanceAtoms, references: seeds))
            }
        }
        try topology.rebuildVirtualSites(&positions)
        guard positions.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("nonfinite finite-cluster coordinates") }
        func latticeImage(_ shift: VivoVector3D) throws -> [Int64] {
            guard let cell else {
                guard shift.norm <= cfg.cycleToleranceNM else { throw VivoChemistryError.invalid("nonperiodic physical coordinates changed") }
                return [0, 0, 0]
            }
            let det = cell.a.dot(cell.b.cross(cell.c))
            let values = [cell.b.cross(cell.c).dot(shift) / det,
                          cell.c.cross(cell.a).dot(shift) / det,
                          cell.a.cross(cell.b).dot(shift) / det]
            guard values.allSatisfy({ $0.isFinite && abs($0) < 0x1p50 }) else {
                throw VivoChemistryError.resourceLimit("QM/MM lattice image integer range")
            }
            let integers = values.map { Int64($0.rounded()) }
            let rebuilt = cell.a * Double(integers[0]) + cell.b * Double(integers[1]) + cell.c * Double(integers[2])
            guard (rebuilt - shift).norm <= cfg.cycleToleranceNM else {
                throw VivoChemistryError.invalid("finite-cluster image is not an integer lattice translation")
            }
            return integers
        }
        let siteParents = Dictionary(uniqueKeysWithValues: zip(topology.siteGraph.sites,topology.siteGraph.physicalAncestors).map { ($0.0.siteParticle,$0.1[0]) })
        let images = try system.particles.map { particle -> VivoQMMMParticleImage in
            if let atom = topology.particleToAtom[Int(particle.index)] {
                return .init(particleIndex: particle.index, structureAtomIndex: atom,
                    moleculeIndex: topology.atomMolecule[Int(atom)],
                    latticeImage: try latticeImage(positions[Int(particle.index)] - source[Int(particle.index)]))
            }
            guard let ancestor = siteParents[particle.index], let parent = topology.particleToAtom[Int(ancestor)] else {
                throw VivoChemistryError.invalid("unowned virtual site in QM/MM image mapping")
            }
            return .init(particleIndex: particle.index, structureAtomIndex: nil,
                moleculeIndex: topology.atomMolecule[Int(parent)], latticeImage: nil)
        }
        return .init(schema: VivoQMMMFiniteCluster.schema, structureFingerprint: document.structureFingerprint,
            classicalSystemFingerprint: try system.fingerprint(), sourceFrameFingerprint: sourceFrameFingerprint,
            sourcePeriodicCell: cell, configuration: cfg, anchorAtomIndex: anchor,
            atomToParticle: topology.atomToParticle, particlePositionsNM: positions, particleImages: images,
            moleculeAtomIndices: topology.molecules, linearVirtualSites: topology.sites, virtualSiteDefinitions: system.virtualSiteDefinitions, imagePairEvaluations: work)
    }
}
