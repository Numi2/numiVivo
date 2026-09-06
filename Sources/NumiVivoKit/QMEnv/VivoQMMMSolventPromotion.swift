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
        self.radiusNM=radiusNM;self.solventResidueNames=solventResidueNames
        self.distanceAtoms=distanceAtoms;self.periodicMinimumImage=periodicMinimumImage
        self.maximumPromotedResidues=maximumPromotedResidues;self.maximumPromotedAtoms=maximumPromotedAtoms
    }
    public func validate() throws {
        let normalized=solventResidueNames.map { $0.trimmingCharacters(in:.whitespacesAndNewlines).uppercased() }
        guard radiusNM.isFinite,radiusNM>0,radiusNM<=10,!normalized.isEmpty,
              normalized.allSatisfy({!$0.isEmpty && $0.utf8.count<=64}),Set(normalized).count==normalized.count,
              (1...100000).contains(maximumPromotedResidues),(1...1000000).contains(maximumPromotedAtoms) else {
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
    public static let schema="numivivo.org/qmmm-solvent-promotion/v1"
    public let schema:String
    public let policy:VivoQMMMSolventPromotionPolicy
    public let originalQMAtomIndices:[UInt32]
    public let promotedResidues:[VivoQMMMPromotedSolventResidue]
    public let qmAtomIndices:[UInt32]
    public let addedAlphaElectrons:Int
    public let addedBetaElectrons:Int
    public let periodicCellUsed:Bool
    public let requiresUnwrappedFiniteClusterForQMMM:Bool
    public let request:VivoQMMMRegionRequest
    public let interpretation:String
}

/// Deterministic whole-residue explicit-solvent promotion. Selection can use
/// periodic minimum-image distances, but the resulting finite QM/MM calculation
/// still requires explicitly unwrapped coordinates when the source is periodic.
/// Electron assignment is automatic only for newly promoted closed-shell solvent
/// residues; odd-electron/ambiguous solvent must be selected with an explicit
/// user-authored QM region instead of an invented spin state.
public enum VivoQMMMSolventPromotion {
    public static let interpretation="whole-residue explicit-solvent QM promotion by declared geometric shell; closed-shell electron count from atomic numbers/formal charges; periodic selection does not itself unwrap coordinates; no solvent free-energy convergence or transition-state qualification implied"

    private static func atomParticleMap(document:VivoMolecularStructureDocument,system:VivoClassicalSystem,
                                        positions:[VivoVector3D]) throws -> [Int] {
        try document.validate()
        try VivoClassicalSystemValidator.validate(system,atomCount:UInt32(document.structure.atoms.count))
        guard system.structureFingerprint==document.structureFingerprint,positions.count==system.particles.count,
              positions.allSatisfy(\.isFinite) else {throw VivoChemistryError.invalid("solvent promotion source identity or coordinates")}
        var map=[Int](repeating:-1,count:document.structure.atoms.count)
        for (i,p) in system.particles.enumerated() where p.role == .atom {
            guard let atom=p.atomIndex,Int(atom)<map.count,map[Int(atom)] == -1 else {
                throw VivoChemistryError.invalid("solvent promotion requires atom/particle bijection")
            }
            map[Int(atom)]=i
        }
        guard map.allSatisfy({$0>=0}) else {throw VivoChemistryError.invalid("solvent promotion has unmapped physical atom")}
        return map
    }

    public static func promote(document:VivoMolecularStructureDocument,system:VivoClassicalSystem,
                               particlePositionsNM positions:[VivoVector3D],base:VivoQMMMRegionRequest,
                               policy:VivoQMMMSolventPromotionPolicy) throws -> VivoQMMMSolventPromotionResult {
        try policy.validate()
        let map=try atomParticleMap(document:document,system:system,positions:positions)
        let structure=document.structure,n=structure.atoms.count
        guard !base.qmAtomIndices.isEmpty,Set(base.qmAtomIndices).count==base.qmAtomIndices.count,
              base.qmAtomIndices.allSatisfy({Int($0)<n}),base.alphaElectrons>=0,base.betaElectrons>=0 else {
            throw VivoChemistryError.invalid("solvent promotion base QM region")
        }
        let original=Set(base.qmAtomIndices),allowed=Set(policy.solventResidueNames.map{$0.trimmingCharacters(in:.whitespacesAndNewlines).uppercased()})
        // A preselected solvent residue must already be whole. Otherwise the
        // missing spin/electron ownership of the partial molecule is ambiguous.
        for residue in structure.residues where allowed.contains(residue.name.uppercased()) {
            let overlap=residue.atomIndices.filter{original.contains($0)}
            guard overlap.isEmpty || overlap.count==residue.atomIndices.count else {
                throw VivoChemistryError.invalid("base QM region partially selects solvent residue \(residue.index); whole-molecule ownership is required")
            }
        }
        let seedAtoms=base.qmAtomIndices
        let cell=policy.periodicMinimumImage ? structure.periodicCell:nil
        if let cell {guard cell.isValid else {throw VivoChemistryError.invalid("invalid periodic cell for solvent promotion")}}
        func distance(_ a:UInt32,_ b:UInt32) throws -> Double {
            var d=positions[map[Int(a)]]-positions[map[Int(b)]]
            if let cell {d=try cell.minimumImage(d)}
            return d.norm
        }
        var candidates:[(residue:VivoMolecularResidue,distance:Double,electrons:Int)]=[]
        for residue in structure.residues where allowed.contains(residue.name.uppercased()) {
            guard !residue.atomIndices.isEmpty else {throw VivoChemistryError.invalid("empty solvent residue \(residue.index)")}
            if residue.atomIndices.allSatisfy({original.contains($0)}) {continue}
            let distanceAtoms=residue.atomIndices.filter { atom in
                policy.distanceAtoms == .allAtoms || structure.atoms[Int(atom)].element.atomicNumber>1
            }
            guard !distanceAtoms.isEmpty else {throw VivoChemistryError.unsupported("heavy-atom solvent distance requested for residue without heavy atoms")}
            var minimum=Double.infinity
            for a in distanceAtoms {for q in seedAtoms {minimum=min(minimum,try distance(a,q))}}
            if minimum<=policy.radiusNM {
                let electrons=residue.atomIndices.reduce(0) { partial,atom in
                    partial+Int(structure.atoms[Int(atom)].element.atomicNumber)-Int(structure.atoms[Int(atom)].formalCharge)
                }
                guard electrons>=0,electrons%2==0 else {
                    throw VivoChemistryError.unsupported("automatic solvent promotion requires a closed-shell residue; residue \(residue.index) has \(electrons) electrons")
                }
                candidates.append((residue,minimum,electrons))
            }
        }
        candidates.sort { a,b in a.distance == b.distance ? a.residue.index<b.residue.index:a.distance<b.distance }
        let atomCount=candidates.reduce(0){$0+$1.residue.atomIndices.count}
        guard candidates.count<=policy.maximumPromotedResidues,atomCount<=policy.maximumPromotedAtoms else {
            throw VivoChemistryError.resourceLimit("solvent shell contains \(candidates.count) residues/\(atomCount) atoms; increase explicit promotion capacity rather than silently truncating")
        }
        var selected=original,alpha=0,beta=0,records:[VivoQMMMPromotedSolventResidue]=[]
        for item in candidates {
            selected.formUnion(item.residue.atomIndices)
            alpha+=item.electrons/2;beta+=item.electrons/2
            records.append(.init(residueIndex:item.residue.index,residueName:item.residue.name,
                atomIndices:item.residue.atomIndices.sorted(),minimumDistanceNM:item.distance,addedElectrons:item.electrons))
        }
        var request=base
        request.qmAtomIndices=selected.sorted()
        request.alphaElectrons += alpha;request.betaElectrons += beta
        return .init(schema:VivoQMMMSolventPromotionResult.schema,policy:policy,
            originalQMAtomIndices:base.qmAtomIndices.sorted(),promotedResidues:records,qmAtomIndices:request.qmAtomIndices,
            addedAlphaElectrons:alpha,addedBetaElectrons:beta,periodicCellUsed:cell != nil,
            requiresUnwrappedFiniteClusterForQMMM:structure.periodicCell != nil,request:request,interpretation:interpretation)
    }
}

public extension VivoQMMMCompiler {
    /// New strict promotion API. The legacy `promoteSolventResidues` remains for
    /// source compatibility but does not carry electron/provenance contracts.
    static func promoteSolventRegion(document:VivoMolecularStructureDocument,system:VivoClassicalSystem,
                                     particlePositionsNM positions:[VivoVector3D],request:VivoQMMMRegionRequest,
                                     policy:VivoQMMMSolventPromotionPolicy) throws -> VivoQMMMSolventPromotionResult {
        try VivoQMMMSolventPromotion.promote(document:document,system:system,particlePositionsNM:positions,base:request,policy:policy)
    }
}
