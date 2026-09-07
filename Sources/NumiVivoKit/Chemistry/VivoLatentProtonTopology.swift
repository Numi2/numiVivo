import Foundation

public enum VivoLatentAtomIdentity: Codable, Sendable, Equatable, Hashable {
    case source(UInt32)
    case generatedHydrogen(String)

    fileprivate var stableKey: String {
        switch self {
        case .source(let index): return "S:\(index)"
        case .generatedHydrogen(let identifier): return "H:\(identifier)"
        }
    }
}

public struct VivoLatentProtonSlot: Codable, Sendable, Equatable {
    public let identity: VivoLatentAtomIdentity
    public let unionParticleIndex: UInt32
    public let presentStateIdentifiers: [String]
    public let anchorBondCount: Int
    public let anchorAngleCount: Int
    public let anchorTorsionCount: Int
    public let constrained: Bool
}

public struct VivoLatentProtonStateSystem: Codable, Sendable, Equatable {
    public let stateIdentifier: String
    public let boundProtonOffset: Int
    public let ghostParticleIndices: [UInt32]
    public let system: VivoClassicalSystem
}

public struct VivoLatentProtonTopologyResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/latent-proton-topology/v1"
    public let schema: String
    public let catalogEvidenceFingerprint: VivoFingerprint
    public let unionStructure: VivoMolecularStructure
    public let unionStructureFingerprint: VivoFingerprint
    public let atomIdentities: [VivoLatentAtomIdentity]
    public let latentProtons: [VivoLatentProtonSlot]
    public let states: [VivoLatentProtonStateSystem]
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

/// Compiles prepared fixed-topology microstates into an always-present proton
/// manifold. Every possible proton remains a positive-mass atom in every state.
/// When absent chemically, its charge/LJ interactions are zeroed while a common
/// set of valence anchors is retained in every state. Those anchor partition
/// factors are therefore state-independent bookkeeping degrees of freedom and
/// must be covered by the explicit reference/reservoir free-energy calibration.
///
/// v1 is deliberately strict: native compiled systems must be atom-only, use a
/// Lorentz-Berthelot or geometric mixing rule, and carry no polarization/virtual
/// sites. This prevents silently changing the measure while creating the union.
public enum VivoLatentProtonTopologyCompiler {
    public static let interpretation = "Common positive-mass latent-proton topology for discrete constant-pH and alchemical sampling. All source atoms plus every explicitly generated proton are retained. Chemically absent protons are nonbonded ghosts with zero charge and zero Lennard-Jones energy, while proton-involving valence anchors are identical across all states. The construction removes atom creation/deletion from state moves; absolute proton chemistry still requires explicit calibrated semigrand reference free energies."

    private struct StateView {
        let state: VivoPreparedChemicalState
        let structure: VivoMolecularStructure
        let system: VivoClassicalSystem
        let identityByPreparedAtom: [VivoLatentAtomIdentity]
        let unionByPreparedAtom: [UInt32]
        let present: Set<VivoLatentAtomIdentity>
    }

    private struct BondKey: Hashable { let a: UInt32; let b: UInt32 }
    private struct AngleKey: Hashable { let a: UInt32; let b: UInt32; let c: UInt32 }
    private struct TorsionKey: Hashable { let a: UInt32; let b: UInt32; let c: UInt32; let d: UInt32; let periodicity: UInt16; let improper: Bool }
    private struct ConstraintKey: Hashable { let a: UInt32; let b: UInt32 }
    private struct ExceptionKey: Hashable { let a: UInt32; let b: UInt32 }

    public static func compile(catalogRequest: VivoChemicalStateCatalogRequest,
                               catalog: VivoChemicalStateCatalogResult) throws -> VivoLatentProtonTopologyResult {
        try VivoChemicalStateCatalog.validate(catalog, request: catalogRequest)
        guard catalog.states.count >= 2,
              catalog.states.count == catalogRequest.entries.count else {
            throw VivoChemistryError.invalid("latent-proton catalog state count")
        }
        let source = catalogRequest.entries[0].preparation.structure
        _ = try VivoStructureValidator.validate(source)

        var placementByIdentifier: [String: VivoHydrogenPlacement] = [:]
        for entry in catalogRequest.entries {
            for placement in entry.preparation.addHydrogens {
                if let existing = placementByIdentifier[placement.identifier] {
                    guard existing == placement else {
                        throw VivoChemistryError.invalid("generated proton identity has inconsistent placement across chemical states")
                    }
                } else {
                    placementByIdentifier[placement.identifier] = placement
                }
            }
        }

        var identitySet = Set<VivoLatentAtomIdentity>()
        for atom in source.atoms { identitySet.insert(.source(atom.index)) }
        for identifier in placementByIdentifier.keys { identitySet.insert(.generatedHydrogen(identifier)) }
        let identities = identitySet.sorted { $0.stableKey < $1.stableKey }
        guard identities.count <= Int(UInt32.max) else {
            throw VivoChemistryError.resourceLimit("latent-proton union atom count")
        }
        let unionIndex = Dictionary(uniqueKeysWithValues: identities.enumerated().map { ($0.element, UInt32($0.offset)) })

        var views: [StateView] = []
        views.reserveCapacity(catalog.states.count)
        for stateIndex in catalog.states.indices {
            let state = catalog.states[stateIndex]
            let result = state.preparationResult
            guard let compiled = result.compiledForceField else {
                throw VivoChemistryError.unsupported("latent-proton topology requires native compiled force fields for every state")
            }
            let system = compiled.system
            try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(result.structure.atoms.count))
            guard system.particles.count == result.structure.atoms.count,
                  system.particles.enumerated().allSatisfy({ index, particle in
                      particle.role == .atom && particle.atomIndex == UInt32(index) && particle.massDa > 0
                  }),
                  system.linearVirtualSites?.isEmpty != false,
                  system.virtualSiteDefinitions?.isEmpty != false,
                  system.polarization == nil,
                  system.mixingRule != .explicitPairTable,
                  system.nonbondedTypePairs?.isEmpty != false else {
                throw VivoChemistryError.unsupported("latent-proton v1 requires atom-only nonpolarizable native systems without explicit pair tables")
            }
            var generatedByPrepared: [UInt32: String] = [:]
            for (identifier, prepared) in result.hydrogenIndices {
                guard generatedByPrepared[prepared] == nil else {
                    throw VivoChemistryError.invalid("multiple generated proton identities map to one prepared atom")
                }
                generatedByPrepared[prepared] = identifier
            }
            var atomIDs: [VivoLatentAtomIdentity] = []
            atomIDs.reserveCapacity(result.structure.atoms.count)
            for prepared in result.structure.atoms.indices {
                let p = UInt32(prepared)
                let identity: VivoLatentAtomIdentity
                if prepared < result.preparedToSource.count, let sourceIndex = result.preparedToSource[prepared] {
                    identity = .source(sourceIndex)
                } else if let identifier = generatedByPrepared[p] {
                    identity = .generatedHydrogen(identifier)
                } else {
                    throw VivoChemistryError.invalid("prepared atom lacks persistent source or generated-proton identity")
                }
                guard unionIndex[identity] != nil else {
                    throw VivoChemistryError.invalid("prepared atom identity missing from latent union")
                }
                atomIDs.append(identity)
            }
            guard Set(atomIDs).count == atomIDs.count else {
                throw VivoChemistryError.invalid("chemical state maps multiple prepared atoms to one latent identity")
            }
            let mapped = try atomIDs.map { identity -> UInt32 in
                guard let value = unionIndex[identity] else { throw VivoChemistryError.invalid("latent union mapping") }
                return value
            }
            views.append(.init(state: state, structure: result.structure, system: system,
                               identityByPreparedAtom: atomIDs, unionByPreparedAtom: mapped,
                               present: Set(atomIDs)))
        }

        let presentInAll = identities.filter { identity in views.allSatisfy { $0.present.contains(identity) } }
        let latent = Set(identities).subtracting(presentInAll)
        guard !latent.isEmpty else {
            throw VivoChemistryError.invalid("latent-proton compiler received no state-varying atom identities")
        }
        for identity in latent {
            let atomicNumbers = views.compactMap { view -> UInt16? in
                guard let prepared = view.identityByPreparedAtom.firstIndex(of: identity) else { return nil }
                return view.structure.atoms[prepared].element.atomicNumber
            }
            guard !atomicNumbers.isEmpty, atomicNumbers.allSatisfy({ $0 == 1 }) else {
                throw VivoChemistryError.invalid("only hydrogen atoms may vary across latent-proton states")
            }
        }

        let unionStructure = try makeUnionStructure(source: source, identities: identities,
                                                    unionIndex: unionIndex,
                                                    placements: placementByIdentifier,
                                                    views: views)
        let unionStructureID = try VivoStructureCodec.fingerprint(unionStructure)

        var templateParticle: [VivoLatentAtomIdentity: VivoClassicalParticle] = [:]
        for identity in identities {
            for view in views {
                if let prepared = view.identityByPreparedAtom.firstIndex(of: identity) {
                    let particle = view.system.particles[prepared]
                    if let previous = templateParticle[identity] {
                        guard abs(previous.massDa - particle.massDa) <= 1e-8 else {
                            throw VivoChemistryError.invalid("persistent atom mass differs across chemical states")
                        }
                    } else { templateParticle[identity] = particle }
                }
            }
        }
        guard templateParticle.count == identities.count else {
            throw VivoChemistryError.invalid("latent union particle templates incomplete")
        }

        var anchorBonds: [BondKey: VivoHarmonicBond] = [:]
        var anchorAngles: [AngleKey: VivoHarmonicAngle] = [:]
        var anchorTorsions: [TorsionKey: VivoPeriodicTorsion] = [:]
        var anchorConstraints: [ConstraintKey: VivoDistanceConstraint] = [:]
        for view in views {
            for bond in view.system.bonds {
                let mapped = VivoHarmonicBond(a:view.unionByPreparedAtom[Int(bond.a)], b:view.unionByPreparedAtom[Int(bond.b)],
                                              lengthNM:bond.lengthNM, forceConstant:bond.forceConstant)
                if latent.contains(identities[Int(mapped.a)]) || latent.contains(identities[Int(mapped.b)]) {
                    try mergeBond(mapped, into:&anchorBonds)
                }
            }
            for angle in view.system.angles {
                let mapped = VivoHarmonicAngle(a:view.unionByPreparedAtom[Int(angle.a)], b:view.unionByPreparedAtom[Int(angle.b)],
                                               c:view.unionByPreparedAtom[Int(angle.c)], angleRadians:angle.angleRadians,
                                               forceConstant:angle.forceConstant)
                if [mapped.a,mapped.b,mapped.c].contains(where:{ latent.contains(identities[Int($0)]) }) {
                    try mergeAngle(mapped, into:&anchorAngles)
                }
            }
            for torsion in view.system.torsions {
                let mapped = VivoPeriodicTorsion(a:view.unionByPreparedAtom[Int(torsion.a)], b:view.unionByPreparedAtom[Int(torsion.b)],
                                                 c:view.unionByPreparedAtom[Int(torsion.c)], d:view.unionByPreparedAtom[Int(torsion.d)],
                                                 periodicity:torsion.periodicity, phaseRadians:torsion.phaseRadians,
                                                 barrierKJPerMol:torsion.barrierKJPerMol, improper:torsion.improper)
                if [mapped.a,mapped.b,mapped.c,mapped.d].contains(where:{ latent.contains(identities[Int($0)]) }) {
                    try mergeTorsion(mapped, into:&anchorTorsions)
                }
            }
            for constraint in view.system.constraints {
                let mapped = VivoDistanceConstraint(a:view.unionByPreparedAtom[Int(constraint.a)],
                                                    b:view.unionByPreparedAtom[Int(constraint.b)], distanceNM:constraint.distanceNM)
                if latent.contains(identities[Int(mapped.a)]) || latent.contains(identities[Int(mapped.b)]) {
                    try mergeConstraint(mapped, into:&anchorConstraints)
                }
            }
        }
        guard latent.allSatisfy({ identity in
            guard let index = unionIndex[identity] else { return false }
            return anchorBonds.values.contains(where:{ $0.a == index || $0.b == index }) ||
                   anchorConstraints.values.contains(where:{ $0.a == index || $0.b == index })
        }) else {
            throw VivoChemistryError.unsupported("every latent proton requires a common bond or distance-constraint anchor")
        }

        var commonNonlatentConstraints: [ConstraintKey: VivoDistanceConstraint]? = nil
        for view in views {
            var mapped: [ConstraintKey: VivoDistanceConstraint] = [:]
            for constraint in view.system.constraints {
                let value = VivoDistanceConstraint(a:view.unionByPreparedAtom[Int(constraint.a)],
                                                   b:view.unionByPreparedAtom[Int(constraint.b)], distanceNM:constraint.distanceNM)
                if !latent.contains(identities[Int(value.a)]), !latent.contains(identities[Int(value.b)]) {
                    let key = constraintKey(value)
                    guard mapped[key] == nil else { throw VivoChemistryError.invalid("duplicate mapped nonlatent constraint") }
                    mapped[key] = value
                }
            }
            if let common = commonNonlatentConstraints {
                guard common == mapped else {
                    throw VivoChemistryError.unsupported("nonlatent distance constraints differ across chemical states")
                }
            } else { commonNonlatentConstraints = mapped }
        }
        let commonConstraints = Array((commonNonlatentConstraints ?? [:]).values) + Array(anchorConstraints.values)
        let sortedConstraints = commonConstraints.sorted(by: constraintOrdering)

        var outputStates: [VivoLatentProtonStateSystem] = []
        outputStates.reserveCapacity(views.count)
        for view in views {
            var particles: [VivoClassicalParticle] = []
            particles.reserveCapacity(identities.count)
            var ghosts: [UInt32] = []
            for (i, identity) in identities.enumerated() {
                let u = UInt32(i)
                if let prepared = view.identityByPreparedAtom.firstIndex(of: identity) {
                    var particle = view.system.particles[prepared]
                    particle = .init(index:u, atomIndex:u, typeIdentifier:particle.typeIdentifier,
                                     role:.atom, massDa:particle.massDa, chargeE:particle.chargeE,
                                     sigmaNM:particle.sigmaNM, epsilonKJPerMol:particle.epsilonKJPerMol)
                    particles.append(particle)
                } else {
                    guard latent.contains(identity), let template = templateParticle[identity] else {
                        throw VivoChemistryError.invalid("only latent protons may be absent from a state")
                    }
                    particles.append(.init(index:u, atomIndex:u,
                        typeIdentifier:"LATENT:\(identity.stableKey)", role:.atom,
                        massDa:template.massDa, chargeE:0, sigmaNM:0, epsilonKJPerMol:0))
                    ghosts.append(u)
                }
            }

            var bonds: [VivoHarmonicBond] = []
            for bond in view.system.bonds {
                let value = VivoHarmonicBond(a:view.unionByPreparedAtom[Int(bond.a)], b:view.unionByPreparedAtom[Int(bond.b)],
                                             lengthNM:bond.lengthNM, forceConstant:bond.forceConstant)
                if ![value.a,value.b].contains(where:{ latent.contains(identities[Int($0)]) }) { bonds.append(value) }
            }
            bonds.append(contentsOf:anchorBonds.values)
            bonds = uniqueBonds(bonds)

            var angles: [VivoHarmonicAngle] = []
            for angle in view.system.angles {
                let value = VivoHarmonicAngle(a:view.unionByPreparedAtom[Int(angle.a)], b:view.unionByPreparedAtom[Int(angle.b)],
                                              c:view.unionByPreparedAtom[Int(angle.c)], angleRadians:angle.angleRadians,
                                              forceConstant:angle.forceConstant)
                if ![value.a,value.b,value.c].contains(where:{ latent.contains(identities[Int($0)]) }) { angles.append(value) }
            }
            angles.append(contentsOf:anchorAngles.values)
            angles = uniqueAngles(angles)

            var torsions: [VivoPeriodicTorsion] = []
            for torsion in view.system.torsions {
                let value = VivoPeriodicTorsion(a:view.unionByPreparedAtom[Int(torsion.a)], b:view.unionByPreparedAtom[Int(torsion.b)],
                                                c:view.unionByPreparedAtom[Int(torsion.c)], d:view.unionByPreparedAtom[Int(torsion.d)],
                                                periodicity:torsion.periodicity, phaseRadians:torsion.phaseRadians,
                                                barrierKJPerMol:torsion.barrierKJPerMol, improper:torsion.improper)
                if ![value.a,value.b,value.c,value.d].contains(where:{ latent.contains(identities[Int($0)]) }) { torsions.append(value) }
            }
            torsions.append(contentsOf:anchorTorsions.values)
            torsions = uniqueTorsions(torsions)

            var exceptions: [VivoNonbondedException] = []
            for exception in view.system.nonbondedExceptions {
                exceptions.append(.init(a:view.unionByPreparedAtom[Int(exception.a)],
                                        b:view.unionByPreparedAtom[Int(exception.b)],
                                        coulombScale:exception.coulombScale,
                                        lennardJonesScale:exception.lennardJonesScale,
                                        sigmaOverrideNM:exception.sigmaOverrideNM,
                                        epsilonOverrideKJPerMol:exception.epsilonOverrideKJPerMol))
            }
            exceptions = uniqueExceptions(exceptions)
            var metadata = view.system.metadata
            metadata["latentProtonTopology"] = "numivivo.org/latent-proton-topology/v1"
            metadata["latentChemicalState"] = view.state.identifier
            metadata["latentGhostParticles"] = ghosts.map(String.init).joined(separator:",")
            metadata["chemicalStateCatalogEvidence"] = catalog.evidenceFingerprint.hex
            let system = VivoClassicalSystem(identifier:"latent-\(view.state.identifier)",
                structureFingerprint:unionStructureID,
                parameterSourceFingerprints:view.system.parameterSourceFingerprints,
                mixingRule:view.system.mixingRule, particles:particles,
                bonds:bonds.sorted(by:bondOrdering), angles:angles.sorted(by:angleOrdering),
                torsions:torsions.sorted(by:torsionOrdering), constraints:sortedConstraints,
                linearVirtualSites:nil, virtualSiteDefinitions:nil, nonbondedTypePairs:nil,
                nonbondedExceptions:exceptions.sorted(by:exceptionOrdering), metadata:metadata,
                polarization:nil)
            try VivoClassicalSystemValidator.validate(system, atomCount:UInt32(unionStructure.atoms.count))
            outputStates.append(.init(stateIdentifier:view.state.identifier,
                                      boundProtonOffset:view.state.boundProtonOffset,
                                      ghostParticleIndices:ghosts, system:system))
        }
        let manifolds = try outputStates.map { try physicalManifoldFingerprint($0.system) }
        guard Set(manifolds).count == 1 else {
            throw VivoChemistryError.invalid("latent-proton compiler failed to produce one physical manifold")
        }

        var slots: [VivoLatentProtonSlot] = []
        for identity in latent.sorted(by:{ $0.stableKey < $1.stableKey }) {
            let index = unionIndex[identity]!
            let present = views.filter { $0.present.contains(identity) }.map { $0.state.identifier }
            slots.append(.init(identity:identity, unionParticleIndex:index,
                presentStateIdentifiers:present,
                anchorBondCount:anchorBonds.values.filter({ $0.a == index || $0.b == index }).count,
                anchorAngleCount:anchorAngles.values.filter({ [$0.a,$0.b,$0.c].contains(index) }).count,
                anchorTorsionCount:anchorTorsions.values.filter({ [$0.a,$0.b,$0.c,$0.d].contains(index) }).count,
                constrained:anchorConstraints.values.contains(where:{ $0.a == index || $0.b == index })))
        }
        struct Evidence: Codable {
            let schema: String
            let catalog: VivoFingerprint
            let unionStructure: VivoFingerprint
            let identities: [VivoLatentAtomIdentity]
            let slots: [VivoLatentProtonSlot]
            let states: [VivoLatentProtonStateSystem]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema:"numivivo.org/latent-proton-topology-evidence/v1",
            catalog:catalog.evidenceFingerprint, unionStructure:unionStructureID,
            identities:identities, slots:slots, states:outputStates)))
        return .init(schema:VivoLatentProtonTopologyResult.schema,
                     catalogEvidenceFingerprint:catalog.evidenceFingerprint,
                     unionStructure:unionStructure,
                     unionStructureFingerprint:unionStructureID,
                     atomIdentities:identities, latentProtons:slots, states:outputStates,
                     interpretation:interpretation, evidenceFingerprint:evidenceID)
    }

    private static func makeUnionStructure(source: VivoMolecularStructure,
                                           identities: [VivoLatentAtomIdentity],
                                           unionIndex: [VivoLatentAtomIdentity: UInt32],
                                           placements: [String: VivoHydrogenPlacement],
                                           views: [StateView]) throws -> VivoMolecularStructure {
        var atoms: [VivoAtom] = []
        atoms.reserveCapacity(identities.count)
        var residues = source.residues
        for r in residues.indices { residues[r].atomIndices = [] }
        var positions = source.conformers.map { [VivoVector3D](repeating:.zero, count:identities.count) }
        var bonds: [VivoBond] = []
        var sourceToUnion: [UInt32: UInt32] = [:]
        for atom in source.atoms {
            guard let u = unionIndex[.source(atom.index)] else { throw VivoChemistryError.invalid("source atom missing from union") }
            sourceToUnion[atom.index] = u
            var copy = atom; copy.index = u
            atoms.append(copy)
            if let residue = copy.residueIndex { residues[Int(residue)].atomIndices.append(u) }
            for c in source.conformers.indices { positions[c][Int(u)] = source.conformers[c].positionsNM[Int(atom.index)] }
        }
        for bond in source.bonds {
            guard let a = sourceToUnion[bond.atomA], let b = sourceToUnion[bond.atomB] else {
                throw VivoChemistryError.invalid("source bond missing from union mapping")
            }
            bonds.append(.init(atomA:a, atomB:b, order:bond.order))
        }
        var occupiedNames = Dictionary(grouping:atoms.compactMap { atom -> (UInt32,String)? in
            atom.residueIndex.map { ($0,atom.name) }
        }, by:{ $0.0 }).mapValues { Set($0.map(\.1)) }
        for identity in identities {
            guard case .generatedHydrogen(let identifier) = identity,
                  let u = unionIndex[identity], let placement = placements[identifier],
                  let parent = sourceToUnion[placement.parent] else { continue }
            let parentAtom = atoms[Int(parent)]
            let residue = parentAtom.residueIndex
            var name = placement.name
            if let residue, occupiedNames[residue, default:[]].contains(name) {
                name = "LH\(u)"
            }
            if let residue { occupiedNames[residue, default:[]].insert(name) }
            atoms.append(.init(index:u, name:name, element:.init(atomicNumber:1, symbol:"H"),
                               residueIndex:residue, isHetero:parentAtom.isHetero))
            if let residue { residues[Int(residue)].atomIndices.append(u) }
            bonds.append(.init(atomA:parent, atomB:u))
            guard let template = views.first(where:{ $0.present.contains(identity) }),
                  let prepared = template.identityByPreparedAtom.firstIndex(of:identity),
                  template.structure.conformers.count == source.conformers.count else {
                throw VivoChemistryError.invalid("generated proton lacks coordinate template")
            }
            for c in source.conformers.indices {
                positions[c][Int(u)] = template.structure.conformers[c].positionsNM[prepared]
            }
        }
        atoms.sort { $0.index < $1.index }
        for i in atoms.indices { guard atoms[i].index == UInt32(i) else { throw VivoChemistryError.invalid("union atom ordering") } }
        for r in residues.indices { residues[r].atomIndices.sort() }
        var conformers = source.conformers
        for c in conformers.indices {
            conformers[c].positionsNM = positions[c]
            conformers[c].potentialEnergyKJPerMol = nil
        }
        var metadata = source.metadata
        metadata["latentProtonUnion"] = "true"
        metadata["latentProtonUnionSchema"] = "numivivo.org/latent-proton-topology/v1"
        let union = VivoMolecularStructure(identifier:source.identifier + "-latent-proton-union",
            atoms:atoms, residues:residues, bonds:bonds, conformers:conformers,
            periodicCell:source.periodicCell, metadata:metadata)
        _ = try VivoStructureValidator.validate(union)
        return union
    }

    private static func physicalManifoldFingerprint(_ system: VivoClassicalSystem) throws -> VivoFingerprint {
        struct Particle: Codable { let index:UInt32; let atomIndex:UInt32?; let role:VivoParticleRole; let massDa:Double }
        struct Signature: Codable {
            let schema:String; let particles:[Particle]; let constraints:[VivoDistanceConstraint]
            let linear:[VivoLinearVirtualSite]; let dependent:[VivoDependentSite]
        }
        let signature = Signature(schema:"numivivo.org/constant-ph-physical-manifold/v1",
            particles:system.particles.map { .init(index:$0.index, atomIndex:$0.atomIndex, role:$0.role, massDa:$0.massDa) },
            constraints:system.constraints, linear:system.linearVirtualSites ?? [], dependent:system.virtualSiteDefinitions ?? [])
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(signature))
    }

    private static func bondKey(_ v:VivoHarmonicBond)->BondKey { .init(a:min(v.a,v.b),b:max(v.a,v.b)) }
    private static func angleKey(_ v:VivoHarmonicAngle)->AngleKey {
        v.a <= v.c ? .init(a:v.a,b:v.b,c:v.c) : .init(a:v.c,b:v.b,c:v.a)
    }
    private static func torsionKey(_ v:VivoPeriodicTorsion)->TorsionKey {
        let forward=(v.a,v.b,v.c,v.d), reverse=(v.d,v.c,v.b,v.a)
        let useForward = [forward.0,forward.1,forward.2,forward.3].lexicographicallyPrecedes([reverse.0,reverse.1,reverse.2,reverse.3])
        let q = useForward ? forward : reverse
        return .init(a:q.0,b:q.1,c:q.2,d:q.3,periodicity:v.periodicity,improper:v.improper)
    }
    private static func constraintKey(_ v:VivoDistanceConstraint)->ConstraintKey { .init(a:min(v.a,v.b),b:max(v.a,v.b)) }
    private static func exceptionKey(_ v:VivoNonbondedException)->ExceptionKey { .init(a:min(v.a,v.b),b:max(v.a,v.b)) }

    private static func mergeBond(_ v:VivoHarmonicBond,into map:inout [BondKey:VivoHarmonicBond]) throws {
        let k=bondKey(v); if let x=map[k],x != v && !(x.a==v.b && x.b==v.a && x.lengthNM==v.lengthNM && x.forceConstant==v.forceConstant) { throw VivoChemistryError.unsupported("latent proton bond anchor parameters differ across states") }; map[k]=canonical(v)
    }
    private static func mergeAngle(_ v:VivoHarmonicAngle,into map:inout [AngleKey:VivoHarmonicAngle]) throws {
        let k=angleKey(v); let c=canonical(v); if let x=map[k],x != c { throw VivoChemistryError.unsupported("latent proton angle anchor parameters differ across states") }; map[k]=c
    }
    private static func mergeTorsion(_ v:VivoPeriodicTorsion,into map:inout [TorsionKey:VivoPeriodicTorsion]) throws {
        let k=torsionKey(v); let c=canonical(v); if let x=map[k],x != c { throw VivoChemistryError.unsupported("latent proton torsion anchor parameters differ across states") }; map[k]=c
    }
    private static func mergeConstraint(_ v:VivoDistanceConstraint,into map:inout [ConstraintKey:VivoDistanceConstraint]) throws {
        let k=constraintKey(v); let c=canonical(v); if let x=map[k],x != c { throw VivoChemistryError.unsupported("latent proton constraint anchor parameters differ across states") }; map[k]=c
    }

    private static func canonical(_ v:VivoHarmonicBond)->VivoHarmonicBond { v.a<=v.b ? v:.init(a:v.b,b:v.a,lengthNM:v.lengthNM,forceConstant:v.forceConstant) }
    private static func canonical(_ v:VivoHarmonicAngle)->VivoHarmonicAngle { v.a<=v.c ? v:.init(a:v.c,b:v.b,c:v.a,angleRadians:v.angleRadians,forceConstant:v.forceConstant) }
    private static func canonical(_ v:VivoPeriodicTorsion)->VivoPeriodicTorsion {
        let f=[v.a,v.b,v.c,v.d],r=[v.d,v.c,v.b,v.a]
        return f.lexicographicallyPrecedes(r) || f==r ? v:.init(a:v.d,b:v.c,c:v.b,d:v.a,periodicity:v.periodicity,phaseRadians:v.phaseRadians,barrierKJPerMol:v.barrierKJPerMol,improper:v.improper)
    }
    private static func canonical(_ v:VivoDistanceConstraint)->VivoDistanceConstraint { v.a<=v.b ? v:.init(a:v.b,b:v.a,distanceNM:v.distanceNM) }
    private static func canonical(_ v:VivoNonbondedException)->VivoNonbondedException { v.a<=v.b ? v:.init(a:v.b,b:v.a,coulombScale:v.coulombScale,lennardJonesScale:v.lennardJonesScale,sigmaOverrideNM:v.sigmaOverrideNM,epsilonOverrideKJPerMol:v.epsilonOverrideKJPerMol) }

    private static func uniqueBonds(_ values:[VivoHarmonicBond])->[VivoHarmonicBond] { Array(Dictionary(uniqueKeysWithValues:values.map { (bondKey($0),canonical($0)) }).values) }
    private static func uniqueAngles(_ values:[VivoHarmonicAngle])->[VivoHarmonicAngle] { Array(Dictionary(uniqueKeysWithValues:values.map { (angleKey($0),canonical($0)) }).values) }
    private static func uniqueTorsions(_ values:[VivoPeriodicTorsion])->[VivoPeriodicTorsion] { Array(Dictionary(uniqueKeysWithValues:values.map { (torsionKey($0),canonical($0)) }).values) }
    private static func uniqueExceptions(_ values:[VivoNonbondedException])->[VivoNonbondedException] { Array(Dictionary(uniqueKeysWithValues:values.map { (exceptionKey($0),canonical($0)) }).values) }

    private static func bondOrdering(_ a:VivoHarmonicBond,_ b:VivoHarmonicBond)->Bool { let x=bondKey(a),y=bondKey(b); return (x.a,x.b)<(y.a,y.b) }
    private static func angleOrdering(_ a:VivoHarmonicAngle,_ b:VivoHarmonicAngle)->Bool { let x=angleKey(a),y=angleKey(b); return (x.a,x.b,x.c)<(y.a,y.b,y.c) }
    private static func torsionOrdering(_ a:VivoPeriodicTorsion,_ b:VivoPeriodicTorsion)->Bool { let x=torsionKey(a),y=torsionKey(b); if x.a != y.a{return x.a<y.a};if x.b != y.b{return x.b<y.b};if x.c != y.c{return x.c<y.c};if x.d != y.d{return x.d<y.d};if x.periodicity != y.periodicity{return x.periodicity<y.periodicity};return x.improper==false && y.improper==true }
    private static func constraintOrdering(_ a:VivoDistanceConstraint,_ b:VivoDistanceConstraint)->Bool { let x=constraintKey(a),y=constraintKey(b); return (x.a,x.b)<(y.a,y.b) }
    private static func exceptionOrdering(_ a:VivoNonbondedException,_ b:VivoNonbondedException)->Bool { let x=exceptionKey(a),y=exceptionKey(b); return (x.a,x.b)<(y.a,y.b) }
}
