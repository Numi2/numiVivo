import Foundation

public enum VivoQMMMTermOwner: String, Codable, Sendable { case retainedClassical, electronic }
public enum VivoQMMMConstraintPolicy: String, Codable, Sendable {
    case rejectQMConstraints
    /// The supplied constraints remain the same holonomic manifold in every partition.
    case preserveExplicitManifold
}
public enum VivoQMMMEmbeddingBoundary: String, Codable, Sendable {
    case finiteCluster
    /// Localized electronic model with periodic electrostatics, not band-structure QM.
    case periodicElectrostatic
}
public struct VivoQMMMChargeRecipient: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let weight: Double
    public init(particleIndex: UInt32, weight: Double) { self.particleIndex = particleIndex; self.weight = weight }
}
public struct VivoQMMMChargeTransfer: Codable, Sendable, Equatable {
    public let donorParticle: UInt32
    public let recipients: [VivoQMMMChargeRecipient]
    public init(donorParticle: UInt32, recipients: [VivoQMMMChargeRecipient]) {
        self.donorParticle = donorParticle; self.recipients = recipients
    }
}
public struct VivoQMMMHamiltonianConfiguration: Codable, Sendable, Equatable {
    public var region: VivoQMMMRegionRequest
    public var boundary: VivoQMMMEmbeddingBoundary
    public var constraintPolicy: VivoQMMMConstraintPolicy
    /// Fixed, topology-defined coefficients. These affect embedding, not MM self energy.
    public var boundaryChargeTransfers: [VivoQMMMChargeTransfer]
    public var totalCellChargeE: Double
    public var maximumGeneratedExceptions: Int
    public init(region: VivoQMMMRegionRequest, boundary: VivoQMMMEmbeddingBoundary,
                constraintPolicy: VivoQMMMConstraintPolicy = .rejectQMConstraints,
                boundaryChargeTransfers: [VivoQMMMChargeTransfer] = [], totalCellChargeE: Double = 0,
                maximumGeneratedExceptions: Int = 1_000_000) {
        self.region = region; self.boundary = boundary; self.constraintPolicy = constraintPolicy
        self.boundaryChargeTransfers = boundaryChargeTransfers; self.totalCellChargeE = totalCellChargeE
        self.maximumGeneratedExceptions = maximumGeneratedExceptions
    }
}

/// Compiles term ownership into the EXISTING classical representation. The same
/// packer, bonded, neighbor, LJ, PME and constraint kernels evaluate the retained
/// system. This object does not introduce another classical force-field engine.
public struct VivoQMMMHamiltonianPlan: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/qmmm-hamiltonian-plan/v1"
    public static let bondedConvention = "single-link-count-policy-v1: remove bonds with 2 QM, angles with >=2 QM, proper/improper periodic torsions with >=3 QM; preserve declared constraint manifold"
    public let schema: String
    public let sourceSystemFingerprint: VivoFingerprint
    public let structureFingerprint: VivoFingerprint
    public let configuration: VivoQMMMHamiltonianConfiguration
    public let atomToParticle: [UInt32]
    public let qmParticles: [UInt32]
    public let boundaryMMParticles: [UInt32]
    public let bondOwners: [VivoQMMMTermOwner]
    public let angleOwners: [VivoQMMMTermOwner]
    public let torsionOwners: [VivoQMMMTermOwner]
    public let embeddingChargesE: [Double]
    public let retainedSystem: VivoClassicalSystem
    public let expectedQMChargeE: Double
    public let totalMMChargeE: Double
    public let bondedConvention: String

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
    public func validate(document: VivoMolecularStructureDocument, source: VivoClassicalSystem,
                         budget: VivoChemistryBudget = .init()) throws {
        guard try Self.compile(document: document, system: source, configuration: configuration, budget: budget) == self else {
            throw VivoChemistryError.invalid("QM/MM Hamiltonian ownership plan does not reconstruct")
        }
    }
    public static func compile(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                               configuration input: VivoQMMMHamiltonianConfiguration,
                               budget: VivoChemistryBudget = .init()) throws -> Self {
        try budget.validate()
        let topology = try VivoQMMMTopology(document: document, system: system)
        _ = try budget.elements([system.particles.count, 8])
        var cfg = input
        cfg.region.qmAtomIndices.sort()
        cfg.boundaryChargeTransfers = cfg.boundaryChargeTransfers.map {
            .init(donorParticle: $0.donorParticle, recipients: $0.recipients.sorted { $0.particleIndex < $1.particleIndex })
        }.sorted { $0.donorParticle < $1.donorParticle }
        let selected = Set(cfg.region.qmAtomIndices)
        guard selected.count == cfg.region.qmAtomIndices.count,
              selected.allSatisfy({ Int($0) < topology.atomToParticle.count }),
              (0...100000).contains(cfg.region.alphaElectrons), (0...100000).contains(cfg.region.betaElectrons),
              cfg.totalCellChargeE.isFinite, cfg.maximumGeneratedExceptions >= 0,
              cfg.region.hydrogenLinkDistancesNM.values.allSatisfy({ $0.isFinite && $0 > 0 && $0 < 1 }) else {
            throw VivoChemistryError.invalid("QM/MM Hamiltonian selection, electron sector or resource capacity")
        }
        if selected.isEmpty, cfg.region.alphaElectrons != 0 || cfg.region.betaElectrons != 0 {
            throw VivoChemistryError.invalid("all-MM Hamiltonian has no electronic occupations")
        }
        if cfg.boundary == .periodicElectrostatic, abs(cfg.totalCellChargeE) > 1e-8 {
            throw VivoChemistryError.unsupported("periodic QM/MM v1 requires an explicitly neutral prepared cell")
        }
        var qm = Set(selected.map { topology.atomToParticle[Int($0)] })
        for site in topology.siteGraph.sites {
            let owned = site.parentParticles.filter { qm.contains($0) }.count
            guard owned == 0 || owned == site.parentParticles.count else {
                throw VivoChemistryError.unsupported("dependent site crosses the QM/MM ownership boundary")
            }
            if owned > 0 { qm.insert(site.siteParticle) }
        }
        var boundary = Set<UInt32>(), linkCount = 0
        let chemicalPairs = Set(document.structure.bonds.map { pairKey($0.atomA, $0.atomB) })
        for bond in document.structure.bonds where selected.contains(bond.atomA) != selected.contains(bond.atomB) {
            let q = selected.contains(bond.atomA) ? bond.atomA : bond.atomB
            let m = q == bond.atomA ? bond.atomB : bond.atomA
            guard bond.order == .single, cfg.region.hydrogenLinkDistancesNM[document.structure.atoms[Int(q)].element.symbol] != nil else {
                throw VivoChemistryError.unsupported("unparameterized or nonsingle QM/MM boundary bond")
            }
            boundary.insert(topology.atomToParticle[Int(m)]); linkCount += 1
        }
        for bond in system.bonds where qm.contains(bond.a) != qm.contains(bond.b) {
            guard let a = topology.particleToAtom[Int(bond.a)], let b = topology.particleToAtom[Int(bond.b)],
                  chemicalPairs.contains(pairKey(a,b)) else {
                throw VivoChemistryError.unsupported("classical cut bond lacks chemical bond-order identity")
            }
        }
        if cfg.constraintPolicy == .rejectQMConstraints,
           system.constraints.contains(where: { qm.contains($0.a) || qm.contains($0.b) }) {
            throw VivoChemistryError.unsupported("QM constraints require explicit preserveExplicitManifold policy")
        }
        for ancestors in topology.siteGraph.physicalAncestors where ancestors.contains(where: { boundary.contains($0) }) {
            throw VivoChemistryError.unsupported("boundary atom with a dependent charge site needs a parameterized boundary-site model")
        }
        func owner(_ particles: [UInt32], threshold: Int) -> VivoQMMMTermOwner {
            particles.filter { qm.contains($0) }.count >= threshold ? .electronic : .retainedClassical
        }
        let bonds = system.bonds.map { owner([$0.a,$0.b], threshold: 2) }
        let angles = system.angles.map { owner([$0.a,$0.b,$0.c], threshold: 2) }
        let torsions = system.torsions.map { owner([$0.a,$0.b,$0.c,$0.d], threshold: 3) }
        var retained = system
        retained.bonds = zip(system.bonds,bonds).filter { $0.1 == .retainedClassical }.map(\.0)
        retained.angles = zip(system.angles,angles).filter { $0.1 == .retainedClassical }.map(\.0)
        retained.torsions = zip(system.torsions,torsions).filter { $0.1 == .retainedClassical }.map(\.0)
        for index in qm { retained.particles[Int(index)].chargeE = 0 }
        // Preserve original physical charges on the MM side. Embedding transforms
        // are a distinct channel and must not change retained MM self interactions.
        var embedding = retained.particles.map(\.chargeE), donors = Set<UInt32>()
        for transfer in cfg.boundaryChargeTransfers {
            guard boundary.contains(transfer.donorParticle), donors.insert(transfer.donorParticle).inserted,
                  !transfer.recipients.isEmpty, Set(transfer.recipients.map(\.particleIndex)).count == transfer.recipients.count,
                  transfer.recipients.allSatisfy({ Int($0.particleIndex) < system.particles.count && $0.weight.isFinite && $0.weight >= 0 }),
                  abs(transfer.recipients.reduce(0) { $0 + $1.weight } - 1) < 1e-12,
                  let donorAtom = topology.particleToAtom[Int(transfer.donorParticle)] else {
                throw VivoChemistryError.invalid("boundary-charge transfer donor, recipients or weights")
            }
            for recipient in transfer.recipients {
                guard !qm.contains(recipient.particleIndex), !boundary.contains(recipient.particleIndex),
                      let atom = topology.particleToAtom[Int(recipient.particleIndex)],
                      topology.atomMolecule[Int(atom)] == topology.atomMolecule[Int(donorAtom)] else {
                    throw VivoChemistryError.invalid("boundary charge recipient must be a nonboundary MM atom in the same molecule")
                }
                embedding[Int(recipient.particleIndex)] += system.particles[Int(transfer.donorParticle)].chargeE * recipient.weight
            }
            embedding[Int(transfer.donorParticle)] = 0
        }
        guard boundary.allSatisfy({ system.particles[Int($0)].chargeE == 0 || donors.contains($0) }),
              embedding.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("every charged boundary endpoint needs an explicit charge-conserving embedding transfer")
        }
        let mmCharge = retained.particles.reduce(0) { $0 + $1.chargeE }
        let qmCharge = selected.reduce(Double(linkCount)) { $0 + Double(document.structure.atoms[Int($1)].element.atomicNumber) }
            - Double(cfg.region.alphaElectrons + cfg.region.betaElectrons)
        guard abs(embedding.reduce(0,+) - mmCharge) < 1e-8,
              abs(qmCharge + mmCharge - cfg.totalCellChargeE) < 1e-6 else {
            throw VivoChemistryError.invalid("electronic sector, caps, MM charges and declared cell charge disagree")
        }
        let ordered = qm.sorted()
        let product = ordered.count.multipliedReportingOverflow(by: max(0,ordered.count-1))
        guard !product.overflow, product.partialValue / 2 <= cfg.maximumGeneratedExceptions,
              product.partialValue / 2 <= budget.maximumOperatorApplications else {
            throw VivoChemistryError.resourceLimit("QM-QM retained-term exception capacity")
        }
        var exceptions = Dictionary(uniqueKeysWithValues: system.nonbondedExceptions.map { (pairKey($0.a,$0.b),$0) })
        for i in ordered.indices { for j in 0..<i {
            let a = ordered[j], b = ordered[i]
            exceptions[pairKey(a,b)] = .init(a: a, b: b, coulombScale: 0, lennardJonesScale: 0)
        } }
        retained.nonbondedExceptions = exceptions.keys.sorted().map { exceptions[$0]! }
        let sourceID = try system.fingerprint()
        let configID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(cfg))
        retained.metadata["numivivo.qmmm.sourceSystem"] = sourceID.hex
        retained.metadata["numivivo.qmmm.configuration"] = configID.hex
        retained.metadata["numivivo.qmmm.accounting"] = "retained-MM-self+retained-bonded+MM-MM-LJ+QM-MM-LJ; electronic contribution absent"
        try VivoClassicalSystemValidator.validate(retained)
        return .init(schema: Self.schema, sourceSystemFingerprint: sourceID, structureFingerprint: document.structureFingerprint,
            configuration: cfg, atomToParticle: topology.atomToParticle, qmParticles: ordered, boundaryMMParticles: boundary.sorted(),
            bondOwners: bonds, angleOwners: angles, torsionOwners: torsions, embeddingChargesE: embedding, retainedSystem: retained,
            expectedQMChargeE: qmCharge, totalMMChargeE: mmCharge, bondedConvention: Self.bondedConvention)
    }
    private static func pairKey(_ a: UInt32, _ b: UInt32) -> UInt64 { UInt64(min(a,b)) << 32 | UInt64(max(a,b)) }
}
