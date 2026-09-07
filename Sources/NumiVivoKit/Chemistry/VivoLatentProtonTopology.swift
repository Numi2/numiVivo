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

/// Builds one positive-mass particle manifold for an explicit prepared chemical-
/// state catalog. All source atoms remain first and retain their source indices;
/// explicitly generated hydrogens are appended in stable identifier order.
/// Chemically absent protons keep their mass and common valence anchors but have
/// zero charge/LJ interactions. The anchor contribution is state independent and
/// therefore belongs in the explicit state/reference calibration, not an implicit
/// atom-insertion Jacobian.
public enum VivoLatentProtonTopologyCompiler {
    public static let interpretation = "Common positive-mass latent-proton topology for discrete constant-pH and alchemical sampling. Source atoms retain source ordering; every explicitly generated proton is appended once. Chemically absent protons are zero-charge/zero-Lennard-Jones ghosts with common valence anchors. The construction removes atom creation/deletion from state moves, while absolute proton chemistry remains conditional on explicit calibrated semigrand reference free energies."

    private struct StateView {
        let state: VivoPreparedChemicalState
        let structure: VivoMolecularStructure
        let system: VivoClassicalSystem
        let identities: [VivoLatentAtomIdentity]
        let unionIndexByPreparedAtom: [UInt32]
        let present: Set<VivoLatentAtomIdentity>
    }
    private struct BondKey: Hashable { let a:UInt32; let b:UInt32 }
    private struct AngleKey: Hashable { let a:UInt32; let b:UInt32; let c:UInt32 }
    private struct TorsionKey: Hashable { let a:UInt32; let b:UInt32; let c:UInt32; let d:UInt32; let periodicity:UInt16; let improper:Bool }
    private struct ConstraintKey: Hashable { let a:UInt32; let b:UInt32 }
    private struct ExceptionKey: Hashable { let a:UInt32; let b:UInt32 }

    public static func compile(catalogRequest: VivoChemicalStateCatalogRequest,
                               catalog: VivoChemicalStateCatalogResult) throws -> VivoLatentProtonTopologyResult {
        try VivoChemicalStateCatalog.validate(catalog, request: catalogRequest)
        guard catalog.states.count >= 2, catalog.states.count == catalogRequest.entries.count,
              let firstEntry = catalogRequest.entries.first else {
            throw VivoChemistryError.invalid("latent-proton catalog state count")
        }
        let source = firstEntry.preparation.structure
        _ = try VivoStructureValidator.validate(source)

        var placements:[String:VivoHydrogenPlacement] = [:]
        for entry in catalogRequest.entries {
            for placement in entry.preparation.addHydrogens {
                if let old = placements[placement.identifier] {
                    guard old == placement else {
                        throw VivoChemistryError.invalid("generated proton identity has inconsistent placement across states")
                    }
                } else { placements[placement.identifier] = placement }
            }
        }
        var atomIdentities = source.atoms.map { VivoLatentAtomIdentity.source($0.index) }
        atomIdentities.append(contentsOf: placements.keys.sorted().map { .generatedHydrogen($0) })
        guard atomIdentities.count <= Int(UInt32.max), Set(atomIdentities).count == atomIdentities.count else {
            throw VivoChemistryError.resourceLimit("latent-proton union atom identity capacity")
        }
        let unionByIdentity = Dictionary(uniqueKeysWithValues: atomIdentities.enumerated().map { ($0.element, UInt32($0.offset)) })

        var views:[StateView] = []
        for state in catalog.states {
            let prepared = state.preparationResult
            guard let compiled = prepared.compiledForceField else {
                throw VivoChemistryError.unsupported("latent-proton topology requires a compiled native force field for every state")
            }
            let system = compiled.system
            try VivoClassicalSystemValidator.validate(system, atomCount: UInt32(prepared.structure.atoms.count))
            guard system.particles.count == prepared.structure.atoms.count,
                  system.particles.enumerated().allSatisfy({ pair in
                      pair.element.role == .atom && pair.element.atomIndex == UInt32(pair.offset) && pair.element.massDa > 0
                  }),
                  system.linearVirtualSites?.isEmpty != false,
                  system.virtualSiteDefinitions?.isEmpty != false,
                  system.polarization == nil,
                  system.mixingRule != .explicitPairTable,
                  system.nonbondedTypePairs?.isEmpty != false else {
                throw VivoChemistryError.unsupported("latent-proton v1 requires atom-only nonpolarizable systems without explicit pair tables")
            }
            var generated:[UInt32:String] = [:]
            for (identifier,index) in prepared.hydrogenIndices {
                guard generated[index] == nil else {
                    throw VivoChemistryError.invalid("multiple generated proton identities map to one prepared atom")
                }
                generated[index] = identifier
            }
            var identities:[VivoLatentAtomIdentity] = []
            for p in prepared.structure.atoms.indices {
                let identity:VivoLatentAtomIdentity
                if p < prepared.preparedToSource.count, let sourceIndex = prepared.preparedToSource[p] {
                    identity = .source(sourceIndex)
                } else if let identifier = generated[UInt32(p)] {
                    identity = .generatedHydrogen(identifier)
                } else {
                    throw VivoChemistryError.invalid("prepared atom lacks persistent source/generated-proton identity")
                }
                guard unionByIdentity[identity] != nil else { throw VivoChemistryError.invalid("prepared atom absent from latent union") }
                identities.append(identity)
            }
            guard Set(identities).count == identities.count else {
                throw VivoChemistryError.invalid("chemical state duplicates a latent atom identity")
            }
            let mapped = try identities.map { identity -> UInt32 in
                guard let index = unionByIdentity[identity] else { throw VivoChemistryError.invalid("latent union atom mapping") }
                return index
            }
            views.append(.init(state:state, structure:prepared.structure, system:system,
                               identities:identities, unionIndexByPreparedAtom:mapped,
                               present:Set(identities)))
        }

        let persistent = Set(atomIdentities.filter { identity in views.allSatisfy { $0.present.contains(identity) } })
        let latent = Set(atomIdentities).subtracting(persistent)
        guard !latent.isEmpty else { throw VivoChemistryError.invalid("chemical-state catalog has no varying proton topology") }
        for identity in latent {
            let elements = views.compactMap { view -> UInt16? in
                guard let p = view.identities.firstIndex(of:identity) else { return nil }
                return view.structure.atoms[p].element.atomicNumber
            }
            guard !elements.isEmpty, elements.allSatisfy({ $0 == 1 }) else {
                throw VivoChemistryError.invalid("only hydrogens may vary across latent-proton states")
            }
        }

        let unionStructure = try makeUnionStructure(source:source, atomIdentities:atomIdentities,
                                                    unionByIdentity:unionByIdentity, placements:placements,
                                                    views:views)
        let unionStructureID = try VivoStructureCodec.fingerprint(unionStructure)

        var templateParticle:[VivoLatentAtomIdentity:VivoClassicalParticle] = [:]
        for identity in atomIdentities {
            for view in views {
                guard let p = view.identities.firstIndex(of:identity) else { continue }
                let particle = view.system.particles[p]
                if let old = templateParticle[identity] {
                    guard abs(old.massDa-particle.massDa) <= 1e-8 else {
                        throw VivoChemistryError.invalid("persistent particle mass differs across states")
                    }
                } else { templateParticle[identity] = particle }
            }
        }
        guard templateParticle.count == atomIdentities.count else { throw VivoChemistryError.invalid("latent particle templates incomplete") }

        var anchorBonds:[BondKey:VivoHarmonicBond] = [:]
        var anchorAngles:[AngleKey:VivoHarmonicAngle] = [:]
        var anchorTorsions:[TorsionKey:VivoPeriodicTorsion] = [:]
        var anchorConstraints:[ConstraintKey:VivoDistanceConstraint] = [:]
        for view in views {
            for term in view.system.bonds {
                let v = VivoHarmonicBond(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], lengthNM:term.lengthNM, forceConstant:term.forceConstant)
                if touchesLatent([v.a,v.b], latent:latent, identities:atomIdentities) { try mergeBond(v,&anchorBonds) }
            }
            for term in view.system.angles {
                let v = VivoHarmonicAngle(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], c:view.unionIndexByPreparedAtom[Int(term.c)], angleRadians:term.angleRadians, forceConstant:term.forceConstant)
                if touchesLatent([v.a,v.b,v.c], latent:latent, identities:atomIdentities) { try mergeAngle(v,&anchorAngles) }
            }
            for term in view.system.torsions {
                let v = VivoPeriodicTorsion(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], c:view.unionIndexByPreparedAtom[Int(term.c)], d:view.unionIndexByPreparedAtom[Int(term.d)], periodicity:term.periodicity, phaseRadians:term.phaseRadians, barrierKJPerMol:term.barrierKJPerMol, improper:term.improper)
                if touchesLatent([v.a,v.b,v.c,v.d], latent:latent, identities:atomIdentities) { try mergeTorsion(v,&anchorTorsions) }
            }
            for term in view.system.constraints {
                let v = VivoDistanceConstraint(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], distanceNM:term.distanceNM)
                if touchesLatent([v.a,v.b], latent:latent, identities:atomIdentities) { try mergeConstraint(v,&anchorConstraints) }
            }
        }
        for identity in latent {
            guard let u=unionByIdentity[identity],
                  anchorBonds.values.contains(where:{ $0.a==u || $0.b==u }) || anchorConstraints.values.contains(where:{ $0.a==u || $0.b==u }) else {
                throw VivoChemistryError.unsupported("every latent proton requires a common bond or distance anchor")
            }
        }

        var commonNonlatentConstraints:[ConstraintKey:VivoDistanceConstraint]? = nil
        for view in views {
            var mapped:[ConstraintKey:VivoDistanceConstraint] = [:]
            for term in view.system.constraints {
                let v=canonical(VivoDistanceConstraint(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], distanceNM:term.distanceNM))
                if !touchesLatent([v.a,v.b], latent:latent, identities:atomIdentities) {
                    let key=constraintKey(v); guard mapped[key]==nil else { throw VivoChemistryError.invalid("duplicate mapped nonlatent constraint") }; mapped[key]=v
                }
            }
            if let old=commonNonlatentConstraints {
                guard old==mapped else { throw VivoChemistryError.unsupported("nonlatent constraints differ across states") }
            } else { commonNonlatentConstraints=mapped }
        }
        var constraints=Array((commonNonlatentConstraints ?? [:]).values)+Array(anchorConstraints.values)
        constraints.sort(by:constraintLess)

        var output:[VivoLatentProtonStateSystem] = []
        for view in views {
            var particles:[VivoClassicalParticle] = []
            var ghosts:[UInt32] = []
            for (i,identity) in atomIdentities.enumerated() {
                let u=UInt32(i)
                if let p=view.identities.firstIndex(of:identity) {
                    let sourceParticle=view.system.particles[p]
                    particles.append(.init(index:u, atomIndex:u, typeIdentifier:sourceParticle.typeIdentifier,
                                           role:.atom, massDa:sourceParticle.massDa, chargeE:sourceParticle.chargeE,
                                           sigmaNM:sourceParticle.sigmaNM, epsilonKJPerMol:sourceParticle.epsilonKJPerMol))
                } else {
                    guard latent.contains(identity),let template=templateParticle[identity] else { throw VivoChemistryError.invalid("nonlatent atom absent from chemical state") }
                    particles.append(.init(index:u, atomIndex:u, typeIdentifier:"LATENT-\(u)", role:.atom,
                                           massDa:template.massDa, chargeE:0, sigmaNM:0, epsilonKJPerMol:0))
                    ghosts.append(u)
                }
            }
            var bonds:[VivoHarmonicBond]=[]
            for term in view.system.bonds {
                let v=canonical(VivoHarmonicBond(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], lengthNM:term.lengthNM, forceConstant:term.forceConstant))
                if !touchesLatent([v.a,v.b],latent:latent,identities:atomIdentities){bonds.append(v)}
            }
            bonds.append(contentsOf:anchorBonds.values); bonds=try uniqueBonds(bonds); bonds.sort(by:bondLess)
            var angles:[VivoHarmonicAngle]=[]
            for term in view.system.angles {
                let v=canonical(VivoHarmonicAngle(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], c:view.unionIndexByPreparedAtom[Int(term.c)], angleRadians:term.angleRadians, forceConstant:term.forceConstant))
                if !touchesLatent([v.a,v.b,v.c],latent:latent,identities:atomIdentities){angles.append(v)}
            }
            angles.append(contentsOf:anchorAngles.values); angles=try uniqueAngles(angles); angles.sort(by:angleLess)
            var torsions:[VivoPeriodicTorsion]=[]
            for term in view.system.torsions {
                let v=canonical(VivoPeriodicTorsion(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], c:view.unionIndexByPreparedAtom[Int(term.c)], d:view.unionIndexByPreparedAtom[Int(term.d)], periodicity:term.periodicity, phaseRadians:term.phaseRadians, barrierKJPerMol:term.barrierKJPerMol, improper:term.improper))
                if !touchesLatent([v.a,v.b,v.c,v.d],latent:latent,identities:atomIdentities){torsions.append(v)}
            }
            torsions.append(contentsOf:anchorTorsions.values); torsions=try uniqueTorsions(torsions); torsions.sort(by:torsionLess)
            var exceptions:[VivoNonbondedException]=view.system.nonbondedExceptions.map { term in
                canonical(.init(a:view.unionIndexByPreparedAtom[Int(term.a)], b:view.unionIndexByPreparedAtom[Int(term.b)], coulombScale:term.coulombScale, lennardJonesScale:term.lennardJonesScale, sigmaOverrideNM:term.sigmaOverrideNM, epsilonOverrideKJPerMol:term.epsilonOverrideKJPerMol))
            }
            exceptions=try uniqueExceptions(exceptions); exceptions.sort(by:exceptionLess)
            var metadata=view.system.metadata
            metadata["latentProtonTopology"]="numivivo.org/latent-proton-topology/v1"
            metadata["latentChemicalState"]=view.state.identifier
            metadata["latentGhostParticles"]=ghosts.map(String.init).joined(separator:",")
            metadata["chemicalStateCatalogEvidence"]=catalog.evidenceFingerprint.hex
            let system=VivoClassicalSystem(identifier:"latent-\(view.state.identifier)", structureFingerprint:unionStructureID,
                parameterSourceFingerprints:view.system.parameterSourceFingerprints, mixingRule:view.system.mixingRule,
                particles:particles,bonds:bonds,angles:angles,torsions:torsions,constraints:constraints,
                linearVirtualSites:nil,virtualSiteDefinitions:nil,nonbondedTypePairs:nil,
                nonbondedExceptions:exceptions,metadata:metadata,polarization:nil)
            try VivoClassicalSystemValidator.validate(system, atomCount:UInt32(unionStructure.atoms.count))
            output.append(.init(stateIdentifier:view.state.identifier,boundProtonOffset:view.state.boundProtonOffset,
                                ghostParticleIndices:ghosts,system:system))
        }
        let manifolds=try output.map { try physicalManifoldFingerprint($0.system) }
        guard Set(manifolds).count==1 else { throw VivoChemistryError.invalid("latent-proton states do not share one physical manifold") }

        let slots:[VivoLatentProtonSlot]=latent.sorted(by:{$0.stableKey<$1.stableKey}).map { identity in
            let u=unionByIdentity[identity]!
            return .init(identity:identity,unionParticleIndex:u,
                         presentStateIdentifiers:views.filter{$0.present.contains(identity)}.map{$0.state.identifier},
                         anchorBondCount:anchorBonds.values.filter{$0.a==u||$0.b==u}.count,
                         anchorAngleCount:anchorAngles.values.filter{[$0.a,$0.b,$0.c].contains(u)}.count,
                         anchorTorsionCount:anchorTorsions.values.filter{[$0.a,$0.b,$0.c,$0.d].contains(u)}.count,
                         constrained:anchorConstraints.values.contains(where:{$0.a==u||$0.b==u}))
        }
        struct Evidence:Codable{let schema:String;let catalog:VivoFingerprint;let union:VivoFingerprint;let identities:[VivoLatentAtomIdentity];let slots:[VivoLatentProtonSlot];let states:[VivoLatentProtonStateSystem]}
        let evidenceID=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(schema:"numivivo.org/latent-proton-topology-evidence/v1",catalog:catalog.evidenceFingerprint,union:unionStructureID,identities:atomIdentities,slots:slots,states:output)))
        return .init(schema:VivoLatentProtonTopologyResult.schema,catalogEvidenceFingerprint:catalog.evidenceFingerprint,
                     unionStructure:unionStructure,unionStructureFingerprint:unionStructureID,atomIdentities:atomIdentities,
                     latentProtons:slots,states:output,interpretation:interpretation,evidenceFingerprint:evidenceID)
    }

    private static func makeUnionStructure(source:VivoMolecularStructure,atomIdentities:[VivoLatentAtomIdentity],
                                           unionByIdentity:[VivoLatentAtomIdentity:UInt32],placements:[String:VivoHydrogenPlacement],
                                           views:[StateView]) throws -> VivoMolecularStructure {
        var atoms:[VivoMolecularAtom]=[]
        var residues=source.residues; for i in residues.indices{residues[i].atomIndices=[]}
        var conformerPositions=source.conformers.map{_ in [VivoVector3D](repeating:.zero,count:atomIdentities.count)}
        var sourceToUnion:[UInt32:UInt32]=[:]
        for atom in source.atoms {
            guard let u=unionByIdentity[.source(atom.index)] else{throw VivoChemistryError.invalid("source atom missing from latent union")}
            sourceToUnion[atom.index]=u
            atoms.append(.init(index:u,name:atom.name,element:atom.element,isotopeMassNumber:atom.isotopeMassNumber,
                               formalCharge:atom.formalCharge,residueIndex:atom.residueIndex,sourceSerial:atom.sourceSerial,
                               alternateLocation:atom.alternateLocation,occupancy:atom.occupancy,bFactor:atom.bFactor,isHetero:atom.isHetero))
            if let r=atom.residueIndex{residues[Int(r)].atomIndices.append(u)}
            for c in source.conformers.indices{conformerPositions[c][Int(u)]=source.conformers[c].positionsNM[Int(atom.index)]}
        }
        var bonds:[VivoMolecularBond]=try source.bonds.map{bond in
            guard let a=sourceToUnion[bond.atomA],let b=sourceToUnion[bond.atomB] else{throw VivoChemistryError.invalid("source bond missing from latent union")}
            return .init(atomA:a,atomB:b,order:bond.order,stereo:bond.stereo)
        }
        var occupied:[UInt32:Set<String>]=[:]
        for atom in atoms{if let r=atom.residueIndex{occupied[r,default:[]].insert(atom.name)}}
        for identifier in placements.keys.sorted() {
            let identity=VivoLatentAtomIdentity.generatedHydrogen(identifier)
            guard let u=unionByIdentity[identity],let placement=placements[identifier],let parent=sourceToUnion[placement.parent] else{throw VivoChemistryError.invalid("generated proton union placement")}
            let parentAtom=atoms[Int(parent)],residue=parentAtom.residueIndex
            var name=placement.name;if let r=residue,occupied[r,default:[]].contains(name){name="LH\(u)"};if let r=residue{occupied[r,default:[]].insert(name)}
            atoms.append(.init(index:u,name:name,element:.init(atomicNumber:1,symbol:"H"),residueIndex:residue,isHetero:parentAtom.isHetero))
            if let r=residue{residues[Int(r)].atomIndices.append(u)}
            bonds.append(.init(atomA:parent,atomB:u))
            guard let view=views.first(where:{$0.present.contains(identity)}),let p=view.identities.firstIndex(of:identity),view.structure.conformers.count==source.conformers.count else{throw VivoChemistryError.invalid("generated proton coordinate template")}
            for c in source.conformers.indices{conformerPositions[c][Int(u)]=view.structure.conformers[c].positionsNM[p]}
        }
        guard atoms.enumerated().allSatisfy({$0.element.index==UInt32($0.offset)}) else{throw VivoChemistryError.invalid("latent union atom ordering")}
        for i in residues.indices{residues[i].atomIndices.sort()}
        var conformers=source.conformers;for c in conformers.indices{conformers[c].positionsNM=conformerPositions[c];conformers[c].potentialEnergyKJPerMol=nil}
        var metadata=source.metadata;metadata["latentProtonUnion"]="true";metadata["latentProtonUnionSchema"]="numivivo.org/latent-proton-topology/v1"
        let result=VivoMolecularStructure(identifier:source.identifier+"-latent-proton-union",atoms:atoms,bonds:bonds,
                                          residues:residues,chains:source.chains,conformers:conformers,periodicCell:source.periodicCell,metadata:metadata)
        _=try VivoStructureValidator.validate(result);return result
    }

    private static func touchesLatent(_ indices:[UInt32],latent:Set<VivoLatentAtomIdentity>,identities:[VivoLatentAtomIdentity])->Bool{indices.contains{latent.contains(identities[Int($0)])}}
    private static func bondKey(_ v:VivoHarmonicBond)->BondKey{.init(a:min(v.a,v.b),b:max(v.a,v.b))}
    private static func angleKey(_ v:VivoHarmonicAngle)->AngleKey{v.a<=v.c ? .init(a:v.a,b:v.b,c:v.c):.init(a:v.c,b:v.b,c:v.a)}
    private static func torsionKey(_ v:VivoPeriodicTorsion)->TorsionKey{let f=[v.a,v.b,v.c,v.d],r=[v.d,v.c,v.b,v.a],q=f.lexicographicallyPrecedes(r) ? f:r;return .init(a:q[0],b:q[1],c:q[2],d:q[3],periodicity:v.periodicity,improper:v.improper)}
    private static func constraintKey(_ v:VivoDistanceConstraint)->ConstraintKey{.init(a:min(v.a,v.b),b:max(v.a,v.b))}
    private static func exceptionKey(_ v:VivoNonbondedException)->ExceptionKey{.init(a:min(v.a,v.b),b:max(v.a,v.b))}
    private static func canonical(_ v:VivoHarmonicBond)->VivoHarmonicBond{v.a<=v.b ? v:.init(a:v.b,b:v.a,lengthNM:v.lengthNM,forceConstant:v.forceConstant)}
    private static func canonical(_ v:VivoHarmonicAngle)->VivoHarmonicAngle{v.a<=v.c ? v:.init(a:v.c,b:v.b,c:v.a,angleRadians:v.angleRadians,forceConstant:v.forceConstant)}
    private static func canonical(_ v:VivoPeriodicTorsion)->VivoPeriodicTorsion{let f=[v.a,v.b,v.c,v.d],r=[v.d,v.c,v.b,v.a];return f.lexicographicallyPrecedes(r) ? v:.init(a:v.d,b:v.c,c:v.b,d:v.a,periodicity:v.periodicity,phaseRadians:v.phaseRadians,barrierKJPerMol:v.barrierKJPerMol,improper:v.improper)}
    private static func canonical(_ v:VivoDistanceConstraint)->VivoDistanceConstraint{v.a<=v.b ? v:.init(a:v.b,b:v.a,distanceNM:v.distanceNM)}
    private static func canonical(_ v:VivoNonbondedException)->VivoNonbondedException{v.a<=v.b ? v:.init(a:v.b,b:v.a,coulombScale:v.coulombScale,lennardJonesScale:v.lennardJonesScale,sigmaOverrideNM:v.sigmaOverrideNM,epsilonOverrideKJPerMol:v.epsilonOverrideKJPerMol)}
    private static func mergeBond(_ value:VivoHarmonicBond,_ map:inout[BondKey:VivoHarmonicBond])throws{let v=canonical(value),k=bondKey(v);if let x=map[k],x != v{throw VivoChemistryError.unsupported("latent proton bond-anchor parameters differ across states")};map[k]=v}
    private static func mergeAngle(_ value:VivoHarmonicAngle,_ map:inout[AngleKey:VivoHarmonicAngle])throws{let v=canonical(value),k=angleKey(v);if let x=map[k],x != v{throw VivoChemistryError.unsupported("latent proton angle-anchor parameters differ across states")};map[k]=v}
    private static func mergeTorsion(_ value:VivoPeriodicTorsion,_ map:inout[TorsionKey:VivoPeriodicTorsion])throws{let v=canonical(value),k=torsionKey(v);if let x=map[k],x != v{throw VivoChemistryError.unsupported("latent proton torsion-anchor parameters differ across states")};map[k]=v}
    private static func mergeConstraint(_ value:VivoDistanceConstraint,_ map:inout[ConstraintKey:VivoDistanceConstraint])throws{let v=canonical(value),k=constraintKey(v);if let x=map[k],x != v{throw VivoChemistryError.unsupported("latent proton constraint-anchor parameters differ across states")};map[k]=v}
    private static func uniqueBonds(_ values:[VivoHarmonicBond])throws->[VivoHarmonicBond]{var m:[BondKey:VivoHarmonicBond]=[:];for v in values{try mergeBond(v,&m)};return Array(m.values)}
    private static func uniqueAngles(_ values:[VivoHarmonicAngle])throws->[VivoHarmonicAngle]{var m:[AngleKey:VivoHarmonicAngle]=[:];for v in values{try mergeAngle(v,&m)};return Array(m.values)}
    private static func uniqueTorsions(_ values:[VivoPeriodicTorsion])throws->[VivoPeriodicTorsion]{var m:[TorsionKey:VivoPeriodicTorsion]=[:];for v in values{try mergeTorsion(v,&m)};return Array(m.values)}
    private static func uniqueExceptions(_ values:[VivoNonbondedException])throws->[VivoNonbondedException]{var m:[ExceptionKey:VivoNonbondedException]=[:];for value in values{let v=canonical(value),k=exceptionKey(v);if let x=m[k],x != v{throw VivoChemistryError.unsupported("nonbonded exception differs after latent mapping")};m[k]=v};return Array(m.values)}
    private static func pairLess(_ a:(UInt32,UInt32),_ b:(UInt32,UInt32))->Bool{a.0 != b.0 ? a.0<b.0:a.1<b.1}
    private static func bondLess(_ a:VivoHarmonicBond,_ b:VivoHarmonicBond)->Bool{let x=bondKey(a),y=bondKey(b);return pairLess((x.a,x.b),(y.a,y.b))}
    private static func angleLess(_ a:VivoHarmonicAngle,_ b:VivoHarmonicAngle)->Bool{let x=angleKey(a),y=angleKey(b);if x.a != y.a{return x.a<y.a};if x.b != y.b{return x.b<y.b};return x.c<y.c}
    private static func torsionLess(_ a:VivoPeriodicTorsion,_ b:VivoPeriodicTorsion)->Bool{let x=torsionKey(a),y=torsionKey(b);if x.a != y.a{return x.a<y.a};if x.b != y.b{return x.b<y.b};if x.c != y.c{return x.c<y.c};if x.d != y.d{return x.d<y.d};if x.periodicity != y.periodicity{return x.periodicity<y.periodicity};return !x.improper && y.improper}
    private static func constraintLess(_ a:VivoDistanceConstraint,_ b:VivoDistanceConstraint)->Bool{let x=constraintKey(a),y=constraintKey(b);return pairLess((x.a,x.b),(y.a,y.b))}
    private static func exceptionLess(_ a:VivoNonbondedException,_ b:VivoNonbondedException)->Bool{let x=exceptionKey(a),y=exceptionKey(b);return pairLess((x.a,x.b),(y.a,y.b))}

    private static func physicalManifoldFingerprint(_ system:VivoClassicalSystem)throws->VivoFingerprint{
        struct P:Codable{let index:UInt32;let atomIndex:UInt32?;let role:VivoParticleRole;let massDa:Double}
        struct S:Codable{let schema:String;let particles:[P];let constraints:[VivoDistanceConstraint];let linear:[VivoLinearVirtualSite];let dependent:[VivoDependentSite]}
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(S(schema:"numivivo.org/constant-ph-physical-manifold/v1",particles:system.particles.map{.init(index:$0.index,atomIndex:$0.atomIndex,role:$0.role,massDa:$0.massDa)},constraints:system.constraints,linear:system.linearVirtualSites ?? [],dependent:system.virtualSiteDefinitions ?? [])))
    }
}
