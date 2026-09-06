import Foundation
@preconcurrency import Metal

public struct VivoAdaptiveBasisSpecification: Codable, Sendable, Equatable {
    /// Shell nucleusIndex refers to original structure atomIndex, including all
    /// atoms eligible for QM. No atom-index-to-element parameter guessing occurs.
    public let atomBasis: VivoGaussianBasis
    /// Hydrogen shell shapes; nucleusIndex is rebound to each artificial cap.
    public let linkHydrogenShells: [VivoGaussianShell]
    public init(atomBasis: VivoGaussianBasis,linkHydrogenShells: [VivoGaussianShell] = []) {
        self.atomBasis=atomBasis;self.linkHydrogenShells=linkHydrogenShells
    }
    func basis(document: VivoMolecularStructureDocument,region: VivoQMMMRegionRequest) throws -> VivoGaussianBasis {
        try atomBasis.validate(nucleusCount:document.structure.atoms.count)
        let atoms=region.qmAtomIndices.sorted(),selected=Set(atoms)
        let map=Dictionary(uniqueKeysWithValues:atoms.enumerated().map { (Int($0.element),$0.offset) })
        var shells=atomBasis.shells.compactMap { shell -> VivoGaussianShell? in
            guard let slot=map[shell.nucleusIndex] else { return nil }
            return .init(nucleusIndex:slot,angularMomentum:shell.angularMomentum,primitives:shell.primitives)
        }
        guard Set(shells.map(\.nucleusIndex))==Set(0..<atoms.count) else {
            throw VivoChemistryError.invalid("adaptive basis is missing an eligible physical QM atom")
        }
        let cuts=document.structure.bonds.filter { selected.contains($0.atomA) != selected.contains($0.atomB) }
        if !cuts.isEmpty,linkHydrogenShells.isEmpty { throw VivoChemistryError.invalid("adaptive core boundary requires explicit hydrogen-link basis shells") }
        for i in cuts.indices {
            for shell in linkHydrogenShells { shells.append(.init(nucleusIndex:atoms.count+i,angularMomentum:shell.angularMomentum,primitives:shell.primitives)) }
        }
        let basis=VivoGaussianBasis(identifier:atomBasis.identifier+"-adaptive-mapped",shells:shells,source:atomBasis.source)
        try basis.validate(nucleusCount:atoms.count+cuts.count);return basis
    }
}

public struct VivoAdaptiveDynamicsSetup: Sendable {
    public let plan: VivoAdaptiveQMMMPlan
    /// Permanent classical baseline. Polarization is owned by the full partition
    /// evaluation and is intentionally not independently added to this baseline.
    public let retainedBaselineSystem: VivoClassicalSystem
    public let provider: VivoMDCandidateForceProvider
}

private actor VivoAdaptivePartitionExecutor {
    let document: VivoMolecularStructureDocument
    let source: VivoClassicalSystem
    let plan: VivoAdaptiveQMMMPlan
    let dynamics: VivoMDConfiguration
    let basis: VivoAdaptiveBasisSpecification
    let electronic: VivoQMMMDynamicsElectronicConfiguration
    let device: MTLDevice
    let baseline: VivoMDMetalRuntime
    var cache: [String:VivoMDMetalRuntime]=[:]
    var recency: [String]=[]
    var busy=false
    init(document: VivoMolecularStructureDocument,source: VivoClassicalSystem,plan: VivoAdaptiveQMMMPlan,
         dynamics: VivoMDConfiguration,basis: VivoAdaptiveBasisSpecification,electronic: VivoQMMMDynamicsElectronicConfiguration,
         device: MTLDevice,baseline: VivoMDMetalRuntime) {
        self.document=document;self.source=source;self.plan=plan;self.dynamics=dynamics;self.basis=basis
        self.electronic=electronic;self.device=device;self.baseline=baseline
    }
    private func partition(_ indices: [Int],geometry: VivoMDCandidateGeometry) async throws -> VivoMDMetalRuntime {
        let key=indices.map(String.init).joined(separator:",")
        if let runtime=cache[key] { recency.removeAll { $0==key };recency.append(key);return runtime }
        let region=try plan.region(selectedMoleculeIndices:indices)
        let ownership=try VivoQMMMHamiltonianPlan.compile(document:document,system:source,
            configuration:.init(region:region,boundary:.periodicElectrostatic,constraintPolicy:plan.configuration.constraintPolicy,
                boundaryChargeTransfers:plan.configuration.boundaryChargeTransfers),budget:electronic.budget)
        var cfg=electronic
        cfg.basis=try basis.basis(document:document,region:region)
        cfg.periodic.polarization=ownership.retainedSystem.polarization
        let provider: VivoMDCandidateForceProvider
        if cfg.periodic.reciprocalMesh == nil {
            provider=try .hartreeFock(document:document,sourceSystem:source,plan:ownership,configuration:cfg)
        } else {
            provider=try await .hartreeFockPME(document:document,sourceSystem:source,plan:ownership,configuration:cfg,device:device)
        }
        // Evict before allocation. Cached contexts contain no physical history;
        // eviction changes resource reuse, never partition membership or energy.
        while cache.count>=plan.configuration.maximumResidentPartitions,let oldest=recency.first {
            cache.removeValue(forKey:oldest);recency.removeFirst()
        }
        let initial=VivoClassicalInitialState(systemFingerprint:try ownership.retainedSystem.fingerprint(),
            positionsNM:geometry.particlePositionsNM,periodicCell:geometry.periodicCell)
        let runtime=try await VivoMDMetalRuntime.make(system:ownership.retainedSystem,initialState:initial,
            configuration:dynamics,device:device,forceProvider:provider)
        cache[key]=runtime;recency.append(key);return runtime
    }
    func evaluate(_ geometry: VivoMDCandidateGeometry) async throws -> (VivoAdaptiveHamiltonianResult,VivoMDHamiltonianEvaluation) {
        guard !busy else { throw VivoChemistryError.invalid("adaptive executor already has an in-flight geometry") }
        busy=true;defer { busy=false }
        let weights=try plan.weights(document:document,system:source,geometry:geometry,budget:electronic.budget)
        var values: [VivoMDHamiltonianEvaluation]=[]
        for entry in weights.partitions {
            try Task.checkCancellation()
            let runtime=try await partition(entry.selectedMoleculeIndices,geometry:geometry)
            values.append(try await runtime.evaluateHamiltonian(at:geometry))
        }
        let result=try VivoAdaptiveHamiltonian.combine(weights:weights,particleCount:source.particles.count,evaluations:values)
        return (result,try await baseline.evaluateHamiltonian(at:geometry))
    }
}

public enum VivoAdaptiveQMMMDynamics {
    public static func make(document: VivoMolecularStructureDocument,system: VivoClassicalSystem,
                            initialGeometry: VivoMDCandidateGeometry,dynamics: VivoMDConfiguration,
                            adaptive: VivoAdaptiveQMMMConfiguration,basis: VivoAdaptiveBasisSpecification,
                            electronic: VivoQMMMDynamicsElectronicConfiguration,device requested: MTLDevice? = nil) async throws -> VivoAdaptiveDynamicsSetup {
        try dynamics.validate();try electronic.budget.validate();try electronic.periodic.validate()
        guard dynamics.electrostatics == .pme,dynamics.relativeDielectric==1,initialGeometry.periodicCell != nil,
              dynamics.pmeGridDimensions != nil else {
            throw VivoChemistryError.invalid("adaptive periodic Hamiltonian requires vacuum electrostatics and a fixed retained PME mesh")
        }
        let plan=try VivoAdaptiveQMMMPlan.compile(document:document,system:system,configuration:adaptive)
        _ = try plan.weights(document:document,system:system,geometry:initialGeometry,budget:electronic.budget)
        let device=try requested ?? VivoMetalDeviceSelector.productionDevice()
        var base=system;base.polarization=nil
        base.metadata["numivivo.adaptive.sourceSystem"]=(try system.fingerprint()).hex
        base.metadata["numivivo.adaptive.role"]="permanent-classical-baseline; complete partition correction added by provider"
        let baseID=try base.fingerprint()
        let initial=VivoClassicalInitialState(systemFingerprint:baseID,positionsNM:initialGeometry.particlePositionsNM,periodicCell:initialGeometry.periodicCell)
        let baseline=try await VivoMDMetalRuntime.make(system:base,initialState:initial,configuration:dynamics,device:device)
        let executor=VivoAdaptivePartitionExecutor(document:document,source:system,plan:plan,dynamics:dynamics,
            basis:basis,electronic:electronic,device:device,baseline:baseline)
        struct Identity: Encodable {
            let method: String;let plan: VivoAdaptiveQMMMPlan;let dynamics: VivoMDConfiguration
            let basis: VivoAdaptiveBasisSpecification;let electronic: VivoQMMMDynamicsElectronicConfiguration
        }
        let id=try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(method:"complete-partition-adaptive-hf-md/v1",
            plan:plan,dynamics:dynamics,basis:basis,electronic:electronic)))
        let provider=try VivoMDCandidateForceProvider(fingerprint:id,retainedSystemFingerprint:baseID,
            boundary:.periodicElectrostatic,supportsCellMoves:true,maximumAcceptedResidual:1,
            molecularConnectivitySystem:system,evaluate:{ geometry in
                let (result,baseline)=try await executor.evaluate(geometry)
                let force=zip(result.physicalParticleForcesKJPerMolNM,baseline.physicalParticleForcesKJPerMolNM).map(-)
                return try .init(providerFingerprint:id,geometry:geometry,
                    additionalEnergyKJPerMol:result.energyKJPerMol-baseline.energyKJPerMol,
                    physicalParticleForcesKJPerMolNM:force,
                    derivativeMethod:"exact-bounded-multipartition-complete-Hamiltonians; C2-energy-weights-and-transition-forces; fixed-energy-references; stationary-electronic-and-polarization-checks-per-partition",
                    convergenceResidual:result.maximumNormalizedResidual,requiredResidual:1)
            })
        try provider.validate(system:base,configuration:dynamics,cell:initialGeometry.periodicCell)
        return .init(plan:plan,retainedBaselineSystem:base,provider:provider)
    }
}
