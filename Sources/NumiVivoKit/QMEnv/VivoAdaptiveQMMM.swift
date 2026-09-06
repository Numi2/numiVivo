import Foundation

public struct VivoAdaptiveSolventMolecule: Codable, Sendable, Equatable {
    public let identifier: String
    public let atomIndices: [UInt32]
    /// Geometry-independent energy reference for promoting this whole molecule.
    /// It cannot be inferred from a single instantaneous QM energy.
    public let promotionReferenceKJPerMol: Double
    /// U_comp(s)=sum coefficients[j]*s^(j+1); optional, explicitly parameterized.
    public let compensationCoefficientsKJPerMol: [Double]
    public let parameterProvenance: String
    public init(identifier: String,atomIndices: [UInt32],promotionReferenceKJPerMol: Double,
                compensationCoefficientsKJPerMol: [Double] = [],parameterProvenance: String) {
        self.identifier=identifier;self.atomIndices=atomIndices;self.promotionReferenceKJPerMol=promotionReferenceKJPerMol
        self.compensationCoefficientsKJPerMol=compensationCoefficientsKJPerMol;self.parameterProvenance=parameterProvenance
    }
}
public struct VivoAdaptiveQMMMConfiguration: Codable, Sendable, Equatable {
    public var fixedCore: VivoQMMMRegionRequest
    public var centerAtomIndices: [UInt32]
    public var molecules: [VivoAdaptiveSolventMolecule]
    public var innerRadiusNM: Double
    public var outerRadiusNM: Double
    public var maximumSwitchingMolecules: Int
    public var maximumActiveQMAtoms: Int
    public var maximumResidentPartitions: Int
    public var maximumResidentParticleSlots: Int
    public var boundaryChargeTransfers: [VivoQMMMChargeTransfer]
    public var constraintPolicy: VivoQMMMConstraintPolicy
    public init(fixedCore: VivoQMMMRegionRequest,centerAtomIndices: [UInt32],molecules: [VivoAdaptiveSolventMolecule],
                innerRadiusNM: Double,outerRadiusNM: Double,maximumSwitchingMolecules: Int = 4,
                maximumActiveQMAtoms: Int = 128,maximumResidentPartitions: Int = 2,
                maximumResidentParticleSlots: Int = 250_000,boundaryChargeTransfers: [VivoQMMMChargeTransfer] = [],
                constraintPolicy: VivoQMMMConstraintPolicy = .rejectQMConstraints) {
        self.fixedCore=fixedCore;self.centerAtomIndices=centerAtomIndices;self.molecules=molecules
        self.innerRadiusNM=innerRadiusNM;self.outerRadiusNM=outerRadiusNM
        self.maximumSwitchingMolecules=maximumSwitchingMolecules;self.maximumActiveQMAtoms=maximumActiveQMAtoms
        self.maximumResidentPartitions=maximumResidentPartitions;self.maximumResidentParticleSlots=maximumResidentParticleSlots
        self.boundaryChargeTransfers=boundaryChargeTransfers;self.constraintPolicy=constraintPolicy
    }
    public func validate() throws {
        guard !fixedCore.qmAtomIndices.isEmpty,Set(fixedCore.qmAtomIndices).count==fixedCore.qmAtomIndices.count,
              !centerAtomIndices.isEmpty,Set(centerAtomIndices).count==centerAtomIndices.count,
              Set(centerAtomIndices).isSubset(of:Set(fixedCore.qmAtomIndices)),!molecules.isEmpty,
              Set(molecules.map(\.identifier)).count==molecules.count,innerRadiusNM.isFinite,innerRadiusNM>0,
              outerRadiusNM.isFinite,outerRadiusNM>innerRadiusNM,outerRadiusNM-innerRadiusNM>1e-6,
              (0...8).contains(maximumSwitchingMolecules),maximumActiveQMAtoms>0,
              (1...16).contains(maximumResidentPartitions),maximumResidentParticleSlots>0 else {
            throw VivoChemistryError.invalid("adaptive region geometry, identity or execution capacity")
        }
        var selected=Set(fixedCore.qmAtomIndices)
        for molecule in molecules {
            guard !molecule.identifier.isEmpty,!molecule.atomIndices.isEmpty,!molecule.parameterProvenance.isEmpty,
                  molecule.promotionReferenceKJPerMol.isFinite,molecule.compensationCoefficientsKJPerMol.count<=6,
                  molecule.compensationCoefficientsKJPerMol.allSatisfy(\.isFinite),
                  Set(molecule.atomIndices).count==molecule.atomIndices.count,
                  selected.isDisjoint(with:Set(molecule.atomIndices)) else {
                throw VivoChemistryError.invalid("adaptive molecules must be distinct, complete and have explicit energy references")
            }
            selected.formUnion(molecule.atomIndices)
        }
    }
}
public struct VivoAdaptiveMassWeight: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let weight: Double
}
public struct VivoAdaptiveParticleGradient: Codable, Sendable, Equatable {
    public let particleIndex: UInt32
    public let gradientPerNM: VivoVector3D
}
public struct VivoAdaptivePartitionWeight: Codable, Sendable, Equatable {
    public let selectedMoleculeIndices: [Int]
    public let weight: Double
    public let gradients: [VivoAdaptiveParticleGradient]
    public let affineStrainDerivative: VivoQMMatrix
    public let energyReferenceKJPerMol: Double
}
public struct VivoAdaptiveWeightState: Codable, Sendable, Equatable {
    public let geometryFingerprint: VivoFingerprint
    public let partitions: [VivoAdaptivePartitionWeight]
    public let switchingMoleculeIndices: [Int]
    public let compensationEnergyKJPerMol: Double
    public let compensationForcesKJPerMolNM: [VivoVector3D]
    public let compensationAffineStrainDerivativeKJPerMol: VivoQMMatrix
}

/// Exact bounded multipartition weights. No nearest-N membership, hysteresis,
/// omission of nonzero partitions or deletion of transition forces is permitted.
public struct VivoAdaptiveQMMMPlan: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/adaptive-qmmm-plan/v1"
    public let schema: String
    public let configuration: VivoAdaptiveQMMMConfiguration
    public let structureFingerprint: VivoFingerprint
    public let systemFingerprint: VivoFingerprint
    public let moleculeElectronCounts: [Int]
    public let atomToParticle: [UInt32]
    public let centerParticleWeights: [VivoAdaptiveMassWeight]
    public let moleculeParticleWeights: [[VivoAdaptiveMassWeight]]

    public static func compile(document: VivoMolecularStructureDocument,system: VivoClassicalSystem,
                               configuration cfg: VivoAdaptiveQMMMConfiguration) throws -> Self {
        try cfg.validate()
        let topology=try VivoQMMMTopology(document:document,system:system)
        let n=document.structure.atoms.count
        guard (cfg.fixedCore.qmAtomIndices+cfg.molecules.flatMap(\.atomIndices)).allSatisfy({Int($0)<n}),
              Set(cfg.centerAtomIndices.map { topology.atomMolecule[Int($0)] }).count==1 else {
            throw VivoChemistryError.invalid("adaptive core center must lie in one connected molecular component")
        }
        let resident=system.particles.count.multipliedReportingOverflow(by:cfg.maximumResidentPartitions+1)
        guard !resident.overflow,resident.partialValue<=cfg.maximumResidentParticleSlots else {
            throw VivoChemistryError.resourceLimit("adaptive retained-kernel residency exceeds declared particle-slot capacity")
        }
        func weights(_ atoms: [UInt32]) throws -> [VivoAdaptiveMassWeight] {
            let particles=atoms.sorted().map { topology.atomToParticle[Int($0)] }
            let mass=particles.reduce(0) { $0+system.particles[Int($1)].massDa }
            guard mass.isFinite,mass>0 else { throw VivoChemistryError.invalid("adaptive molecule mass") }
            return particles.map { .init(particleIndex:$0,weight:system.particles[Int($0)].massDa/mass) }
        }
        var electrons: [Int]=[],particleWeights: [[VivoAdaptiveMassWeight]]=[]
        for molecule in cfg.molecules {
            let atoms=Set(molecule.atomIndices),group=topology.atomMolecule[Int(molecule.atomIndices[0])]
            guard atoms==Set(topology.molecules[group]) else { throw VivoChemistryError.invalid("adaptive solvent entry is not a whole connected molecule") }
            let charge=molecule.atomIndices.reduce(0) { $0+Int(document.structure.atoms[Int($1)].formalCharge) }
            let count=molecule.atomIndices.reduce(0) { $0+Int(document.structure.atoms[Int($1)].element.atomicNumber) }
            var members=Set(molecule.atomIndices.map { topology.atomToParticle[Int($0)] })
            for (site,parents) in zip(topology.siteGraph.sites,topology.siteGraph.physicalAncestors) where parents.allSatisfy({members.contains($0)}) {
                members.insert(site.siteParticle)
            }
            let mmCharge=members.reduce(0) { $0+system.particles[Int($1)].chargeE }
            guard charge==0,count%2==0,abs(mmCharge)<1e-8 else {
                throw VivoChemistryError.unsupported("adaptive solvent v1 requires neutral closed-shell molecules in both charge representations")
            }
            if system.constraints.contains(where: { members.contains($0.a) || members.contains($0.b) }) {
                throw VivoChemistryError.unsupported("adaptive solvent must use a fixed flexible model; changing solvent constraint manifolds is not supported")
            }
            for site in system.polarization?.sites ?? [] {
                let ids=[site.particleIndex]+[site.zParticle,site.xParticle].compactMap({$0})
                let included=ids.filter { members.contains($0) }.count
                guard included==0 || included==ids.count else { throw VivoChemistryError.unsupported("adaptive solvent intersects an anisotropic polarization frame") }
            }
            electrons.append(count);particleWeights.append(try weights(molecule.atomIndices))
        }
        return .init(schema:Self.schema,configuration:cfg,structureFingerprint:document.structureFingerprint,
            systemFingerprint:try system.fingerprint(),moleculeElectronCounts:electrons,atomToParticle:topology.atomToParticle,
            centerParticleWeights:try weights(cfg.centerAtomIndices),moleculeParticleWeights:particleWeights)
    }
    public func weights(document: VivoMolecularStructureDocument,system: VivoClassicalSystem,
                        geometry: VivoMDCandidateGeometry,budget: VivoChemistryBudget = .init()) throws -> VivoAdaptiveWeightState {
        guard try Self.compile(document:document,system:system,configuration:configuration)==self else {
            throw VivoChemistryError.invalid("adaptive plan source or decoded membership differs")
        }
        try budget.validate()
        let cfg=configuration
        guard let cell=geometry.periodicCell else { throw VivoChemistryError.unsupported("adaptive QM/MM v1 requires a periodic source cell") }
        try cell.validateSingleImageRadius(cfg.outerRadiusNM)
        let cluster=try VivoQMMMClusterBuilder.reconstruct(document:document,system:system,
            particlePositionsNM:geometry.particlePositionsNM,periodicCell:cell,sourceFrameFingerprint:geometry.fingerprint(),
            qmAtomIndices:cfg.fixedCore.qmAtomIndices,configuration:.init(maximumImagePairs:budget.maximumOperatorApplications))
        let positions=cluster.particlePositionsNM
        func center(_ weights: [VivoAdaptiveMassWeight]) -> VivoVector3D {
            weights.reduce(.zero) { $0+positions[Int($1.particleIndex)]*$1.weight }
        }
        let core=center(centerParticleWeights)
        struct Switch {
            let value: Double;let derivative: Double;let displacement: VivoVector3D
            let gradients: [VivoAdaptiveParticleGradient]
            let strain: VivoQMMatrix
        }
        var switches: [Switch]=[],mixed: [Int]=[],inside: [Int]=[]
        var compensation=0.0,compForce=[VivoVector3D](repeating:.zero,count:system.particles.count),compStrain=VivoQMMatrix(3,3)
        for (i,weights) in moleculeParticleWeights.enumerated() {
            let d=try cell.minimumImage(center(weights)-core),r=d.norm
            var value=0.0,derivative=0.0
            if r<=cfg.innerRadiusNM { value=1;inside.append(i) }
            else if r<cfg.outerRadiusNM {
                let width=cfg.outerRadiusNM-cfg.innerRadiusNM,t=(r-cfg.innerRadiusNM)/width
                value=pow(1-t,3)*(1+3*t+6*t*t)
                derivative = -30*t*t*(1-t)*(1-t)/width;mixed.append(i)
            }
            let radial=r>0 ? d*(derivative/r):.zero
            var gradients: [VivoAdaptiveParticleGradient]=[]
            for p in weights { gradients.append(.init(particleIndex:p.particleIndex,gradientPerNM:radial*p.weight)) }
            for p in centerParticleWeights { gradients.append(.init(particleIndex:p.particleIndex,gradientPerNM:radial*(-p.weight))) }
            var strain=VivoQMMatrix(3,3)
            let g=[radial.x,radial.y,radial.z],displacement=[d.x,d.y,d.z]
            for a in 0..<3 { for b in 0..<3 { strain[a,b]=g[a]*displacement[b] } }
            switches.append(.init(value:value,derivative:derivative,displacement:d,gradients:gradients,strain:strain))
            var potential=0.0,slope=0.0
            for (j,c) in cfg.molecules[i].compensationCoefficientsKJPerMol.enumerated() {
                potential+=c*pow(value,Double(j+1));slope+=Double(j+1)*c*pow(value,Double(j))
            }
            compensation+=potential
            for term in gradients { compForce[Int(term.particleIndex)]=compForce[Int(term.particleIndex)]-term.gradientPerNM*slope }
            compStrain=try compStrain.adding(strain,scale:slope)
        }
        guard mixed.count<=cfg.maximumSwitchingMolecules else {
            throw VivoChemistryError.resourceLimit("too many switching molecules for the declared exact adaptive Hamiltonian; no partition was truncated")
        }
        let possible=inside+mixed
        let activeCount=cfg.fixedCore.qmAtomIndices.count+possible.reduce(0) { $0+cfg.molecules[$1].atomIndices.count }
        guard activeCount<=cfg.maximumActiveQMAtoms else { throw VivoChemistryError.resourceLimit("adaptive QM atom capacity") }
        let partitions=1<<mixed.count
        _ = try budget.elements([partitions,max(1,activeCount),6])
        var output: [VivoAdaptivePartitionWeight]=[]
        for mask in 0..<partitions {
            var selected=inside,factors: [Double]=[]
            for (bit,i) in mixed.enumerated() {
                let on=(mask&(1<<bit)) != 0
                if on { selected.append(i) }
                factors.append(on ? switches[i].value:1-switches[i].value)
            }
            let weight=factors.reduce(1,*)
            var gradient: [UInt32:VivoVector3D]=[:],strain=VivoQMMatrix(3,3)
            for (bit,i) in mixed.enumerated() {
                let coefficient=factors.enumerated().filter { $0.offset != bit }.reduce(1.0) { $0*$1.element }*((mask&(1<<bit)) != 0 ? 1:-1)
                for term in switches[i].gradients { gradient[term.particleIndex,default:.zero]=gradient[term.particleIndex,default:.zero]+term.gradientPerNM*coefficient }
                strain=try strain.adding(switches[i].strain,scale:coefficient)
            }
            let reference=selected.reduce(0) { $0+cfg.molecules[$1].promotionReferenceKJPerMol }
            guard weight>=0,weight.isFinite,reference.isFinite else { throw VivoChemistryError.convergence("nonfinite adaptive weight or energy reference") }
            output.append(.init(selectedMoleculeIndices:selected.sorted(),weight:weight,
                gradients:gradient.keys.sorted().map { .init(particleIndex:$0,gradientPerNM:gradient[$0]!) },
                affineStrainDerivative:strain,energyReferenceKJPerMol:reference))
        }
        guard abs(output.reduce(0) { $0+$1.weight }-1)<1e-12,compensation.isFinite,compForce.allSatisfy(\.isFinite),compStrain.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("adaptive weights failed partition of unity or compensation is nonfinite")
        }
        return .init(geometryFingerprint:try geometry.fingerprint(),partitions:output,switchingMoleculeIndices:mixed,
            compensationEnergyKJPerMol:compensation,compensationForcesKJPerMolNM:compForce,
            compensationAffineStrainDerivativeKJPerMol:compStrain)
    }
    public func region(selectedMoleculeIndices indices: [Int]) throws -> VivoQMMMRegionRequest {
        guard Set(indices).count==indices.count,indices.allSatisfy({$0>=0 && $0<configuration.molecules.count}) else {
            throw VivoChemistryError.invalid("adaptive partition indices")
        }
        var region=configuration.fixedCore
        for i in indices { region.qmAtomIndices+=configuration.molecules[i].atomIndices;region.alphaElectrons+=moleculeElectronCounts[i]/2;region.betaElectrons+=moleculeElectronCounts[i]/2 }
        region.qmAtomIndices.sort();return region
    }
}

public struct VivoAdaptiveHamiltonianResult: Codable, Sendable, Equatable {
    public let geometryFingerprint: VivoFingerprint
    public let energyKJPerMol: Double
    public let physicalParticleForcesKJPerMolNM: [VivoVector3D]
    public let partitionEnergiesKJPerMol: [Double]
    public let maximumNormalizedResidual: Double
    public let weightState: VivoAdaptiveWeightState
}
public enum VivoAdaptiveHamiltonian {
    /// Energies must be COMPLETE partition Hamiltonians, not QM-only corrections.
    /// Each nonzero partition is evaluated; normalized force mixing is not used.
    public static func combine(weights: VivoAdaptiveWeightState,particleCount: Int,
                               evaluations: [VivoMDHamiltonianEvaluation]) throws -> VivoAdaptiveHamiltonianResult {
        guard !weights.partitions.isEmpty,weights.partitions.count==evaluations.count,
              weights.compensationForcesKJPerMolNM.count==particleCount else { throw VivoChemistryError.invalid("adaptive complete-energy input shape") }
        var energy=weights.compensationEnergyKJPerMol,forces=weights.compensationForcesKJPerMolNM
        for (partition,value) in zip(weights.partitions,evaluations) {
            guard value.energyKJPerMol.isFinite,value.physicalParticleForcesKJPerMolNM.count==particleCount,
                  value.physicalParticleForcesKJPerMolNM.allSatisfy(\.isFinite),
                  value.normalizedForceResidual?.isFinite == true,(value.normalizedForceResidual ?? 2)<=1,(value.normalizedForceResidual ?? -1)>=0,partition.weight.isFinite,partition.weight>=0,
                  partition.energyReferenceKJPerMol.isFinite,partition.gradients.allSatisfy({Int($0.particleIndex)<particleCount && $0.gradientPerNM.isFinite}),
                  try value.evaluatedGeometry.fingerprint()==weights.geometryFingerprint else {
                throw VivoChemistryError.invalid("adaptive partition energy, force or evaluated geometry differs")
            }
            let shifted=value.energyKJPerMol-partition.energyReferenceKJPerMol
            energy+=partition.weight*shifted
            for i in forces.indices { forces[i]=forces[i]+value.physicalParticleForcesKJPerMolNM[i]*partition.weight }
            for term in partition.gradients { forces[Int(term.particleIndex)]=forces[Int(term.particleIndex)]-term.gradientPerNM*shifted }
        }
        guard abs(weights.partitions.reduce(0) { $0+$1.weight }-1)<1e-12,energy.isFinite,forces.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("adaptive energy/force assembly failed")
        }
        return .init(geometryFingerprint:weights.geometryFingerprint,energyKJPerMol:energy,physicalParticleForcesKJPerMolNM:forces,
            partitionEnergiesKJPerMol:evaluations.map(\.energyKJPerMol),maximumNormalizedResidual:evaluations.compactMap(\.normalizedForceResidual).max() ?? 0,weightState:weights)
    }
}
