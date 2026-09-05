import Foundation

public struct VivoQMRegionConfiguration: Codable, Sendable {
    public var atomIndices: [UInt32]
    /// Charge of the capped QM region, NOT the sum of fractional MM charges.
    public var molecularCharge: Int
    public var spinProjectionTwice: Int
    public var promoteWaterWithinNM: Double?
    public var waterResidueNames: [String]
    public var imageAnchorNM: VivoVector3D?
    public init(atomIndices: [UInt32], molecularCharge: Int, spinProjectionTwice: Int = 0,
                promoteWaterWithinNM: Double? = nil, waterResidueNames: [String] = ["HOH","WAT","OPC"],
                imageAnchorNM: VivoVector3D? = nil) {
        self.atomIndices = atomIndices; self.molecularCharge = molecularCharge
        self.spinProjectionTwice = spinProjectionTwice; self.promoteWaterWithinNM = promoteWaterWithinNM
        self.waterResidueNames = waterResidueNames; self.imageAnchorNM = imageAnchorNM
    }
}

public struct VivoQMLinkAtom: Codable, Sendable {
    public let qmAtomIndex: UInt32
    public let mmAtomIndex: UInt32
    public let nucleusIndex: Int
    public let lengthNM: Double
}

public struct VivoQMRegion: Codable, Sendable {
    public let structureFingerprint: VivoFingerprint
    public let systemFingerprint: VivoFingerprint
    public let checkpointFingerprint: VivoFingerprint
    public let frozenEnvironmentFingerprint: VivoFingerprint
    public let configuration: VivoQMRegionConfiguration
    public let qmAtomIndices: [UInt32]
    public let nuclei: [VivoQMNucleus]
    public let pointCharges: [VivoQMPointCharge]
    public let linkAtoms: [VivoQMLinkAtom]
    public let zeroedBoundaryParticleIndices: [UInt32]
    public let promotedWaterResidueIndices: [UInt32]
    public let qmMMVanDerWaalsKJPerMol: Double
    public let electrostaticConvention: String

    public func problem(basisIdentity: String, basis: [VivoCartesianGaussian]) throws -> VivoGaussianProblem {
        let result = VivoGaussianProblem(basisIdentity:basisIdentity,nuclei:nuclei,pointCharges:pointCharges,
                                         basis:basis,molecularCharge:configuration.molecularCharge,
                                         spinProjectionTwice:configuration.spinProjectionTwice)
        try result.validate()
        return result
    }
    /// Adds only the cross-region LJ term. QM/MM nuclear electrostatics is
    /// already included by the integral engine. Frozen MM-MM energy is omitted;
    /// comparison requires the same frozen environment and energy convention.
    public func integrals(basisIdentity: String, basis: [VivoCartesianGaussian],
                          schwarzThreshold: Double = 1e-12) throws -> VivoGaussianIntegralSet {
        guard qmMMVanDerWaalsKJPerMol.isFinite else { throw VivoQMError.invalid("QM/MM LJ constant") }
        var result = try VivoGaussianIntegralEngine.build(problem(basisIdentity:basisIdentity,basis:basis),
                                                           schwarzThreshold:schwarzThreshold)
        result.constantEnergyHartree += qmMMVanDerWaalsKJPerMol/VivoQMUnits.hartreeInKJPerMol
        try result.validate()
        return result
    }
}

public enum VivoQMEnvironmentBuilder {
    private static func validateConfiguration(_ config: VivoQMRegionConfiguration, atomCount: Int) throws {
        guard !config.atomIndices.isEmpty, Set(config.atomIndices).count == config.atomIndices.count,
              config.atomIndices.allSatisfy({ Int($0) < atomCount }),
              (-1000...1000).contains(config.molecularCharge),
              (-96...96).contains(config.spinProjectionTwice), config.waterResidueNames.count <= 32,
              config.imageAnchorNM?.isFinite != false else {
            throw VivoQMError.invalid("QM region selection, charge, or spin sector")
        }
        if let radius = config.promoteWaterWithinNM {
            guard radius.isFinite, (0.01...2).contains(radius) else { throw VivoQMError.invalid("QM water-promotion radius") }
        }
    }
    /// This implementation deliberately supports finite nonperiodic regions or
    /// orthogonal-box nearest-image clusters. It does NOT silently approximate
    /// general triclinic unwrapping or periodic Ewald QM electrostatics.
    private static func displacement(_ value: VivoVector3D, cell: VivoPeriodicCell?) throws -> VivoVector3D {
        guard let cell else { return value }
        let axes = [cell.a,cell.b,cell.c]
        guard cell.isValid, abs(cell.a.dot(cell.b)) < 1e-10*cell.a.norm*cell.b.norm,
              abs(cell.a.dot(cell.c)) < 1e-10*cell.a.norm*cell.c.norm,
              abs(cell.b.dot(cell.c)) < 1e-10*cell.b.norm*cell.c.norm else {
            throw VivoQMError.unsupported("QM/MM finite-cluster extraction requires an orthogonal cell or an explicitly unwrapped nonperiodic snapshot")
        }
        var result = value
        for axis in axes { result = result-axis*(value.dot(axis)/axis.squaredNorm).rounded() }
        return result
    }

    public static func build(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                             checkpoint: VivoMDCheckpoint, configuration config: VivoQMRegionConfiguration) throws -> VivoQMRegion {
        try document.validate()
        let structure = document.structure
        try VivoClassicalSystemValidator.validate(system,atomCount:UInt32(structure.atoms.count))
        try checkpoint.validate(particleCount:system.particles.count)
        let systemID = try system.fingerprint()
        guard system.structureFingerprint == document.structureFingerprint,
              checkpoint.systemFingerprint == systemID else { throw VivoQMError.invalid("QM region structure/system/checkpoint identity mismatch") }
        try validateConfiguration(config,atomCount:structure.atoms.count)
        guard !system.particles.contains(where: { $0.role == .drude }) else {
            throw VivoQMError.unsupported("polarizable Drude QM/MM environment")
        }
        var particleForAtom: [UInt32:UInt32] = [:]
        for particle in system.particles where particle.role == .atom {
            guard let atom = particle.atomIndex, particleForAtom[atom] == nil else {
                throw VivoQMError.invalid("QM/MM requires a bijection between chemical atoms and physical particles")
            }
            particleForAtom[atom] = particle.index
        }
        guard particleForAtom.count == structure.atoms.count else { throw VivoQMError.invalid("QM/MM atom/particle map is incomplete") }
        let position: (UInt32) -> VivoVector3D = { checkpoint.positionsNM[Int(particleForAtom[$0]!)] }
        _ = try displacement(.zero,cell:checkpoint.periodicCell)
        var selected = Set(config.atomIndices), promoted: [UInt32] = []
        if let radius = config.promoteWaterWithinNM {
            for residue in structure.residues where config.waterResidueNames.contains(residue.name) {
                var near = false
                for atom in residue.atomIndices {
                    for seed in config.atomIndices {
                        if try displacement(position(atom)-position(seed),cell:checkpoint.periodicCell).norm <= radius { near = true; break }
                    }
                    if near { break }
                }
                if !near { continue }
                let numbers = residue.atomIndices.map { Int(structure.atoms[Int($0)].element.atomicNumber) }.sorted()
                guard numbers == [1,1,8] else { throw VivoQMError.invalid("water promotion requires a complete H2O residue, not a name-only guess") }
                selected.formUnion(residue.atomIndices); promoted.append(residue.index)
            }
        }
        let atomIndices = selected.sorted(), center = config.imageAnchorNM ?? position(config.atomIndices[0])
        var qmPositions: [UInt32:VivoVector3D] = [:]
        for atom in atomIndices {
            let delta = try displacement(position(atom)-center,cell:checkpoint.periodicCell)
            if let cell = checkpoint.periodicCell {
                guard delta.norm < 0.45*min(cell.a.norm,cell.b.norm,cell.c.norm) else {
                    throw VivoQMError.invalid("QM region is too large for the finite nearest-image cluster convention")
                }
            }
            qmPositions[atom] = center+delta
        }
        // Check that imaging did not tear a selected covalent bond.
        for bond in structure.bonds where selected.contains(bond.atomA) && selected.contains(bond.atomB) {
            let local = qmPositions[bond.atomB]!-qmPositions[bond.atomA]!
            let minimum = try displacement(position(bond.atomB)-position(bond.atomA),cell:checkpoint.periodicCell)
            guard (local-minimum).norm < 1e-7 else { throw VivoQMError.invalid("QM selection crosses inconsistent periodic images; unwrap the molecule first") }
        }
        func bohr(_ p: VivoVector3D) -> [Double] {
            // Translation does not change finite-cluster physics. Centering also
            // avoids loss of significance for large absolute MD coordinates.
            let d = (p-center)/VivoQMUnits.bohrInNM
            return [d.x,d.y,d.z]
        }
        var nuclei = atomIndices.map { atom in
            VivoQMNucleus(sourceAtomID:"atom-\(atom)",nuclearCharge:Int(structure.atoms[Int(atom)].element.atomicNumber),
                          positionBohr:bohr(qmPositions[atom]!))
        }
        let linkLengthsNM: [UInt16:Double] = [6:0.109,7:0.102,8:0.096,9:0.092,15:0.142,16:0.134,17:0.128,35:0.141,53:0.161]
        var boundary = Set<UInt32>(), links: [VivoQMLinkAtom] = []
        let sortedBonds = structure.bonds.sorted {
            let a = $0.canonicalPair, b = $1.canonicalPair
            return a.0 == b.0 ? a.1 < b.1 : a.0 < b.0
        }
        for bond in sortedBonds where selected.contains(bond.atomA) != selected.contains(bond.atomB) {
            guard bond.order == .single else { throw VivoQMError.unsupported("QM/MM hydrogen links may cut only an explicitly single bond") }
            let qm = selected.contains(bond.atomA) ? bond.atomA : bond.atomB
            let mm = qm == bond.atomA ? bond.atomB : bond.atomA
            guard let length = linkLengthsNM[structure.atoms[Int(qm)].element.atomicNumber] else {
                throw VivoQMError.unsupported("no hydrogen-link length for this boundary element")
            }
            let direction = try displacement(position(mm)-position(qm),cell:checkpoint.periodicCell)
            guard direction.norm > length+1e-5 else { throw VivoQMError.invalid("QM/MM cut bond is shorter than its link bond") }
            let link = qmPositions[qm]!+direction*(length/direction.norm)
            links.append(.init(qmAtomIndex:qm,mmAtomIndex:mm,nucleusIndex:nuclei.count,lengthNM:length))
            nuclei.append(.init(sourceAtomID:"link-\(qm)-\(mm)",nuclearCharge:1,positionBohr:bohr(link)))
            boundary.insert(particleForAtom[mm]!)
        }
        let qmParticles = Set(atomIndices.map { particleForAtom[$0]! })
        var omitted = qmParticles
        let virtualDefinitions = Dictionary(uniqueKeysWithValues:(system.linearVirtualSites ?? []).map { ($0.siteParticle,$0) })
        for particle in system.particles where particle.role == .virtualSite {
            guard let rule = virtualDefinitions[particle.index] else { throw VivoQMError.unsupported("QM/MM requires resolved virtual-site ownership") }
            let inside = rule.parentParticles.filter { qmParticles.contains($0) }.count
            if inside == rule.parentParticles.count { omitted.insert(particle.index) }
            else if inside != 0 { throw VivoQMError.invalid("QM boundary splits a virtual-site parent group") }
            if rule.parentParticles.contains(where: { boundary.contains($0) }) { boundary.insert(particle.index) }
        }
        let mmParticles = system.particles.filter { !omitted.contains($0.index) }
        let charges = try mmParticles.filter { !boundary.contains($0.index) && $0.chargeE != 0 }.map { particle in
            VivoQMPointCharge(sourceParticleID:"particle-\(particle.index)",chargeE:particle.chargeE,
                              positionBohr:bohr(center + (try displacement(checkpoint.positionsNM[Int(particle.index)]-center,cell:checkpoint.periodicCell))))
        }
        func pairKey(_ a: UInt32,_ b: UInt32) -> UInt64 { UInt64(min(a,b))<<32 | UInt64(max(a,b)) }
        let exceptions = Dictionary(uniqueKeysWithValues:system.nonbondedExceptions.map { (pairKey($0.a,$0.b),$0) })
        func typeKey(_ a: String,_ b: String) -> String { a <= b ? a+"\u{0}"+b : b+"\u{0}"+a }
        let table = Dictionary(uniqueKeysWithValues:(system.nonbondedTypePairs ?? []).map { (typeKey($0.typeA,$0.typeB),$0) })
        var lj = 0.0
        for q in qmParticles.sorted() { for m in mmParticles {
            let a = system.particles[Int(q)], exception = exceptions[pairKey(q,m.index)]
            let scale = exception?.lennardJonesScale ?? 1
            if scale == 0 { continue }
            var c6: Double, c12: Double
            if system.mixingRule == .explicitPairTable, exception?.sigmaOverrideNM == nil, exception?.epsilonOverrideKJPerMol == nil {
                guard let pair = table[typeKey(a.typeIdentifier,m.typeIdentifier)] else { throw VivoQMError.invalid("QM/MM missing explicit LJ pair") }
                c6 = pair.c6KJNM6PerMol; c12 = pair.c12KJNM12PerMol
            } else {
                if system.mixingRule == .explicitPairTable && (exception?.sigmaOverrideNM == nil || exception?.epsilonOverrideKJPerMol == nil) {
                    throw VivoQMError.unsupported("partial LJ override on an explicit type-pair table")
                }
                let sigma = exception?.sigmaOverrideNM ?? (system.mixingRule == .geometric ? sqrt(a.sigmaNM*m.sigmaNM) : 0.5*(a.sigmaNM+m.sigmaNM))
                let epsilon = exception?.epsilonOverrideKJPerMol ?? sqrt(a.epsilonKJPerMol*m.epsilonKJPerMol)
                c6 = 4*epsilon*pow(sigma,6); c12 = 4*epsilon*pow(sigma,12)
            }
            if c6 == 0 && c12 == 0 { continue }
            let mmPosition = center + (try displacement(checkpoint.positionsNM[Int(m.index)]-center,cell:checkpoint.periodicCell))
            let distance = (qmPositions[a.atomIndex!]!-mmPosition).norm
            guard distance > 1e-8 else { throw VivoQMError.invalid("overlapping QM/MM LJ sites") }
            lj += scale*(c12/pow(distance,12)-c6/pow(distance,6))
        } }
        guard lj.isFinite else { throw VivoQMError.convergence("non-finite QM/MM LJ energy") }
        struct FrozenEnvironment: Encodable {
            let schema: String
            let systemFingerprint: VivoFingerprint
            let qmAtomIndices: [UInt32]
            let charges: [VivoQMPointCharge]
            let mmParticleIndices: [UInt32]
            let mmPositionsBohr: [[Double]]
            let zeroedBoundary: [UInt32]
        }
        let frozen = FrozenEnvironment(schema:"numivivo.org/frozen-qm-environment/v1",systemFingerprint:systemID,
            qmAtomIndices:atomIndices,charges:charges,mmParticleIndices:mmParticles.map(\.index),
            mmPositionsBohr:try mmParticles.map { bohr(center + (try displacement(checkpoint.positionsNM[Int($0.index)]-center,cell:checkpoint.periodicCell))) },
            zeroedBoundary:boundary.sorted())
        let frozenID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(frozen))
        return .init(structureFingerprint:document.structureFingerprint,systemFingerprint:systemID,
                     checkpointFingerprint:try checkpoint.fingerprint(),frozenEnvironmentFingerprint:frozenID,configuration:config,qmAtomIndices:atomIndices,
                     nuclei:nuclei,pointCharges:charges,linkAtoms:links,zeroedBoundaryParticleIndices:boundary.sorted(),
                     promotedWaterResidueIndices:promoted.sorted(),qmMMVanDerWaalsKJPerMol:lj,
                     electrostaticConvention:"finite-cluster electrostatic embedding; Z1 boundary charges; frozen MM-MM constant omitted; no periodic QM Ewald")
    }

    /// Union water promotion across snapshots, then use the returned configuration
    /// unchanged to rebuild every point. This fixes chemical membership, not the
    /// individual geometry or the requirement to validate each path Hamiltonian.
    public static func commonPathSelection(document: VivoMolecularStructureDocument, system: VivoClassicalSystem,
                                            checkpoints: [VivoMDCheckpoint], configuration: VivoQMRegionConfiguration) throws -> VivoQMRegionConfiguration {
        guard (1...64).contains(checkpoints.count) else { throw VivoQMError.capacity("QM region path snapshot count") }
        var selected = Set(configuration.atomIndices)
        for checkpoint in checkpoints {
            let region = try build(document:document,system:system,checkpoint:checkpoint,configuration:configuration)
            selected.formUnion(region.qmAtomIndices)
        }
        var result = configuration
        result.atomIndices = selected.sorted(); result.promoteWaterWithinNM = nil
        if result.imageAnchorNM == nil {
            guard let seed = configuration.atomIndices.first,
                  let particle = system.particles.first(where: { $0.role == .atom && $0.atomIndex == seed }) else {
                throw VivoQMError.invalid("path image-anchor atom mapping")
            }
            result.imageAnchorNM = checkpoints[0].positionsNM[Int(particle.index)]
        }
        return result
    }
}
