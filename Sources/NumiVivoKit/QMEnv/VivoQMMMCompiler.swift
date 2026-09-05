import Foundation

public struct VivoQMMMLinkAtom: Codable, Sendable, Equatable {
    public let qmAtomIndex: UInt32
    public let mmAtomIndex: UInt32
    public let qmNucleusIndex: Int
    public let linkNucleusIndex: Int
    public let distanceBohr: Double
}
public struct VivoQMMMPreparedRegion: Codable, Sendable, Equatable {
    public let structureFingerprint: VivoFingerprint
    public let classicalSystemFingerprint: VivoFingerprint
    public let snapshotFingerprint: VivoFingerprint
    public let frozenEnvironmentFingerprint: VivoFingerprint
    public let regionPolicyFingerprint: VivoFingerprint
    public let electronicSystem: VivoElectronicSystem
    public let qmAtomIndices: [UInt32]
    public let excludedQMParticles: [UInt32]
    public let z1ZeroedMMParticles: [UInt32]
    public let links: [VivoQMMMLinkAtom]
    public let energyConvention: String
}
public struct VivoQMMMRegionRequest: Codable, Sendable, Equatable {
    public var qmAtomIndices: [UInt32]
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public var coordinatesAreUnwrappedFiniteCluster: Bool
    public var hydrogenLinkDistancesNM: [String:Double]
    public init(qmAtomIndices: [UInt32], alphaElectrons: Int, betaElectrons: Int,
                coordinatesAreUnwrappedFiniteCluster: Bool = false,
                hydrogenLinkDistancesNM: [String:Double] = ["C":0.109,"N":0.102,"O":0.096,"S":0.134,"P":0.142]) {
        self.qmAtomIndices=qmAtomIndices;self.alphaElectrons=alphaElectrons;self.betaElectrons=betaElectrons
        self.coordinatesAreUnwrappedFiniteCluster=coordinatesAreUnwrappedFiniteCluster
        self.hydrogenLinkDistancesNM=hydrogenLinkDistancesNM
    }
}
private struct VivoQMMMEnvironmentPosition: Codable, Sendable { let particleIndex: UInt32; let positionNM: VivoVector3D }
private struct VivoQMMMEnvironmentIdentity: Codable, Sendable {
    let classicalSystem: VivoFingerprint; let regionPolicy: VivoFingerprint
    let z1ZeroedParticles: [UInt32]; let mmPositions: [VivoQMMMEnvironmentPosition]
}
private struct VivoQMMMSnapshotIdentity: Codable, Sendable {
    let structure: VivoFingerprint; let classicalSystem: VivoFingerprint; let particlePositionsNM: [VivoVector3D]
}
public struct VivoQMMMLennardJonesResult: Codable, Sendable, Equatable {
    public let energyHartree: Double
    public let qmForcesHartreePerBohr: [VivoVector3D]
}
public enum VivoQMMMCompiler {
    private static func vector(_ p: VivoVector3D) -> SIMD3<Double> { .init(p.x,p.y,p.z) }
    private static func key(_ a: UInt32,_ b: UInt32) -> UInt64 { UInt64(min(a,b))<<32 | UInt64(max(a,b)) }
    private static func mapping(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                positions: [VivoVector3D], request: VivoQMMMRegionRequest) throws -> [Int] {
        try document.validate()
        guard document.structure.atoms.count<=Int(UInt32.max) else { throw VivoChemistryError.resourceLimit("QM/MM atom index width") }
        try VivoClassicalSystemValidator.validate(system,atomCount:UInt32(document.structure.atoms.count))
        let n=document.structure.atoms.count
        guard system.structureFingerprint==document.structureFingerprint,positions.count==system.particles.count,
              positions.allSatisfy(\.isFinite),!request.qmAtomIndices.isEmpty,
              Set(request.qmAtomIndices).count==request.qmAtomIndices.count,request.qmAtomIndices.allSatisfy({Int($0)<n}),
              request.hydrogenLinkDistancesNM.values.allSatisfy({$0.isFinite && $0>0 && $0<1}) else {
            throw VivoChemistryError.invalid("QM/MM source identity, coordinates or region request")
        }
        if document.structure.periodicCell != nil && !request.coordinatesAreUnwrappedFiniteCluster {
            throw VivoChemistryError.unsupported("periodic QM/MM requires an explicitly unwrapped finite cluster; Ewald QM/MM is not implemented")
        }
        var map=[Int](repeating:-1,count:n)
        for (index,particle) in system.particles.enumerated() {
            if particle.role == .drude { throw VivoChemistryError.unsupported("polarizable Drude QM/MM") }
            if particle.role == .atom {
                guard let atom=particle.atomIndex,Int(atom)<n,map[Int(atom)] == -1 else { throw VivoChemistryError.invalid("QM/MM requires a physical-atom/particle bijection") }
                map[Int(atom)]=index
            } else if particle.atomIndex != nil { throw VivoChemistryError.invalid("virtual site cannot own a chemical atom in QM/MM") }
        }
        guard map.allSatisfy({$0>=0}) else { throw VivoChemistryError.invalid("unmapped physical atom in classical system") }
        let rules=system.linearVirtualSites ?? []
        guard Set(rules.map(\.siteParticle))==Set(system.particles.filter{$0.role == .virtualSite}.map(\.index)) else {
            throw VivoChemistryError.invalid("all QM/MM virtual-site positions require construction rules")
        }
        for rule in rules {
            var p=VivoVector3D.zero
            for (parent,weight) in zip(rule.parentParticles,rule.weights) { p=p+positions[Int(parent)]*weight }
            guard (p-positions[Int(rule.siteParticle)]).norm<1e-8 else { throw VivoChemistryError.invalid("stale virtual-site coordinates") }
        }
        return map
    }
    public static func prepare(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                               particlePositionsNM positions: [VivoVector3D], request: VivoQMMMRegionRequest) throws -> VivoQMMMPreparedRegion {
        let map=try mapping(document:document,system:system,positions:positions,request:request)
        let selected=Set(request.qmAtomIndices),atoms=request.qmAtomIndices.sorted()
        var nuclei=atoms.map { atom in
            VivoQMNucleus(atomicNumber:Int(document.structure.atoms[Int(atom)].element.atomicNumber),
                          positionBohr:vector(positions[map[Int(atom)]])/VivoAtomicUnits.bohrInNM,structureAtomIndex:atom)
        }
        let nucleusIndex=Dictionary(uniqueKeysWithValues:atoms.enumerated().map{($0.element,$0.offset)})
        var excluded=Set(atoms.map{UInt32(map[Int($0)])}),zeroed=Set<UInt32>(),links:[VivoQMMMLinkAtom]=[]
        let crossing=document.structure.bonds.filter{selected.contains($0.atomA) != selected.contains($0.atomB)}.sorted {
            key($0.atomA,$0.atomB)<key($1.atomA,$1.atomB)
        }
        for bond in crossing {
            guard bond.order == .single else { throw VivoChemistryError.unsupported("QM/MM boundary cuts must be explicit single bonds; promote the conjugated/metal region") }
            let q=selected.contains(bond.atomA) ? bond.atomA : bond.atomB,m=q==bond.atomA ? bond.atomB : bond.atomA
            let symbol=document.structure.atoms[Int(q)].element.symbol
            guard let distance=request.hydrogenLinkDistancesNM[symbol] else { throw VivoChemistryError.unsupported("no explicit hydrogen-link distance for \(symbol)") }
            let rq=positions[map[Int(q)]],rm=positions[map[Int(m)]],direction=rm-rq,length=direction.norm
            guard length.isFinite,length>distance else { throw VivoChemistryError.invalid("QM/MM cut bond is shorter than its link distance") }
            let position=rq+direction*(distance/length),linkIndex=nuclei.count
            nuclei.append(.init(atomicNumber:1,positionBohr:vector(position)/VivoAtomicUnits.bohrInNM))
            links.append(.init(qmAtomIndex:q,mmAtomIndex:m,qmNucleusIndex:nucleusIndex[q]!,linkNucleusIndex:linkIndex,
                               distanceBohr:distance/VivoAtomicUnits.bohrInNM))
            zeroed.insert(UInt32(map[Int(m)]))
        }
        for rule in system.linearVirtualSites ?? [] {
            let parentQM=rule.parentParticles.filter{excluded.contains($0)}.count
            if parentQM==rule.parentParticles.count { excluded.insert(rule.siteParticle) }
            else if parentQM != 0 { throw VivoChemistryError.unsupported("virtual site spans QM/MM partition; promote its complete molecule") }
            else if rule.parentParticles.contains(where:{zeroed.contains($0)}) {
                throw VivoChemistryError.unsupported("Z1 on a boundary atom owning an MM virtual site requires an explicit charge model")
            }
        }
        let charges=system.particles.filter{!excluded.contains($0.index) && !zeroed.contains($0.index) && $0.chargeE != 0}.map {
            VivoQMPointCharge(chargeE:$0.chargeE,positionBohr:vector(positions[Int($0.index)])/VivoAtomicUnits.bohrInNM,
                              classicalParticleIndex:$0.index)
        }
        let electronic=VivoElectronicSystem(nuclei:nuclei,alphaElectrons:request.alphaElectrons,betaElectrons:request.betaElectrons,pointCharges:charges)
        try electronic.validate()
        let classicalFingerprint=try system.fingerprint()
        let snapshot=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(VivoQMMMSnapshotIdentity(
            structure:document.structureFingerprint,classicalSystem:classicalFingerprint,particlePositionsNM:positions)))
        var canonicalRequest=request;canonicalRequest.qmAtomIndices=atoms
        let policyFingerprint=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(canonicalRequest))
        let environment=VivoQMMMEnvironmentIdentity(classicalSystem:classicalFingerprint,regionPolicy:policyFingerprint,
            z1ZeroedParticles:zeroed.sorted(),mmPositions:system.particles.filter{!excluded.contains($0.index)}.map{
                .init(particleIndex:$0.index,positionNM:positions[Int($0.index)])
            })
        let environmentFingerprint=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(environment))
        return .init(structureFingerprint:document.structureFingerprint,classicalSystemFingerprint:classicalFingerprint,
                     snapshotFingerprint:snapshot,frozenEnvironmentFingerprint:environmentFingerprint,regionPolicyFingerprint:policyFingerprint,electronicSystem:electronic,qmAtomIndices:atoms,
                     excludedQMParticles:excluded.sorted(),z1ZeroedMMParticles:zeroed.sorted(),links:links,
                     energyConvention:"finite-cluster-electrostatic-Z1; E_QM+nucleus-MM+QM-MM-LJ; frozen-MM-self omitted; link atoms have no LJ")
    }
    public static func lennardJones(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                    particlePositionsNM positions: [VivoVector3D], request: VivoQMMMRegionRequest,
                                    budget: VivoChemistryBudget = .init()) throws -> VivoQMMMLennardJonesResult {
        let map=try mapping(document:document,system:system,positions:positions,request:request)
        let prepared=try prepare(document:document,system:system,particlePositionsNM:positions,request:request)
        let excluded=Set(prepared.excludedQMParticles),mm=system.particles.filter{!excluded.contains($0.index)}
        let (work,overflow)=request.qmAtomIndices.count.multipliedReportingOverflow(by:mm.count)
        guard !overflow,work<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("QM/MM pair-work budget") }
        let exceptions=Dictionary(uniqueKeysWithValues:system.nonbondedExceptions.map{(key($0.a,$0.b),$0)})
        func typeKey(_ a: String,_ b: String) -> String { a<=b ? a+"\u{0}"+b : b+"\u{0}"+a }
        let table=Dictionary(uniqueKeysWithValues:(system.nonbondedTypePairs ?? []).map{(typeKey($0.typeA,$0.typeB),$0)})
        var energy=0.0,forces=[VivoVector3D](repeating:.zero,count:request.qmAtomIndices.count)
        for (slot,atom) in request.qmAtomIndices.enumerated() {
            let q=system.particles[map[Int(atom)]]
            for m in mm {
                let exception=exceptions[key(q.index,m.index)],scale=exception?.lennardJonesScale ?? 1
                if scale==0 { continue }
                var c12=0.0,c6=0.0
                if exception?.sigmaOverrideNM != nil || exception?.epsilonOverrideKJPerMol != nil {
                    guard let sigma=exception?.sigmaOverrideNM,let epsilon=exception?.epsilonOverrideKJPerMol else {
                        throw VivoChemistryError.unsupported("QM/MM requires paired sigma/epsilon exception overrides")
                    }
                    c6=4*epsilon*pow(sigma,6);c12=4*epsilon*pow(sigma,12)
                } else if system.mixingRule == .explicitPairTable {
                    guard let pair=table[typeKey(q.typeIdentifier,m.typeIdentifier)] else { throw VivoChemistryError.invalid("missing QM/MM type pair") }
                    c6=pair.c6KJNM6PerMol;c12=pair.c12KJNM12PerMol
                } else {
                    let sigma=system.mixingRule == .geometric ? sqrt(q.sigmaNM*m.sigmaNM) : 0.5*(q.sigmaNM+m.sigmaNM)
                    let epsilon=sqrt(q.epsilonKJPerMol*m.epsilonKJPerMol)
                    c6=4*epsilon*pow(sigma,6);c12=4*epsilon*pow(sigma,12)
                }
                if c6==0 && c12==0 { continue }
                let delta=positions[Int(q.index)]-positions[Int(m.index)],r2=delta.squaredNorm
                guard r2.isFinite,r2>1e-20 else { throw VivoChemistryError.invalid("overlapping QM/MM LJ centers") }
                let inverse2=1/r2,inverse6=inverse2*inverse2*inverse2,inverse12=inverse6*inverse6
                energy += scale*(c12*inverse12-c6*inverse6)
                let coefficient=scale*(12*c12*inverse12-6*c6*inverse6)*inverse2
                forces[slot]=forces[slot]+delta*coefficient
            }
        }
        guard energy.isFinite,forces.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite QM/MM LJ result") }
        return .init(energyHartree:energy/VivoAtomicUnits.hartreeInKJPerMol,
                     qmForcesHartreePerBohr:forces.map{$0*(VivoAtomicUnits.bohrInNM/VivoAtomicUnits.hartreeInKJPerMol)})
    }
    public static func promoteSolventResidues(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                              particlePositionsNM positions: [VivoVector3D], request: VivoQMMMRegionRequest,
                                              radiusNM: Double, solventResidueNames: Set<String> = ["HOH","WAT","SOL"]) throws -> [UInt32] {
        let map=try mapping(document:document,system:system,positions:positions,request:request)
        guard radiusNM.isFinite,radiusNM>0 else { throw VivoChemistryError.invalid("solvent promotion radius") }
        var selected=Set(request.qmAtomIndices)
        for residue in document.structure.residues where solventResidueNames.contains(residue.name) {
            let close=residue.atomIndices.contains { atom in
                request.qmAtomIndices.contains { q in (positions[map[Int(atom)]]-positions[map[Int(q)]]).squaredNorm<=radiusNM*radiusNM }
            }
            if close { selected.formUnion(residue.atomIndices) }
        }
        return selected.sorted()
    }
}
