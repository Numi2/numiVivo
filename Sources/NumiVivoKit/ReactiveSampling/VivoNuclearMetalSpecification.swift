import Foundation

/// Serializable preparation only: executable closures are constructed by the
/// existing fixed/adaptive QM/MM factory, never decoded from an artifact.
public enum VivoNuclearElectronicSpecification: Codable, Sendable, Equatable {
    case classical
    case fixed(document: VivoMolecularStructureDocument, plan: VivoQMMMHamiltonianPlan,
               electronic: VivoQMMMDynamicsElectronicConfiguration)
    case adaptive(document: VivoMolecularStructureDocument, configuration: VivoAdaptiveQMMMConfiguration,
                  basis: VivoAdaptiveBasisSpecification, electronic: VivoQMMMDynamicsElectronicConfiguration)
}
public struct VivoNuclearMetalSpecification: Codable, Sendable, Equatable {
    public let system: VivoClassicalSystem
    public let initialState: VivoClassicalInitialState
    public let dynamics: VivoMDConfiguration
    public let electronic: VivoNuclearElectronicSpecification
    public let maximumMeshPoints: Int
    public init(system: VivoClassicalSystem, initialState: VivoClassicalInitialState, dynamics: VivoMDConfiguration,
                electronic: VivoNuclearElectronicSpecification = .classical, maximumMeshPoints: Int = 1_048_576) {
        self.system = system; self.initialState = initialState; self.dynamics = dynamics
        self.electronic = electronic; self.maximumMeshPoints = maximumMeshPoints
    }
    public func definition() throws -> VivoNuclearPotentialDefinition {
        try VivoClassicalSystemValidator.validate(system); try initialState.validate(particleCount: system.particles.count)
        try dynamics.validate()
        guard initialState.systemFingerprint == (try system.fingerprint()), system.constraints.isEmpty,
              !system.particles.contains(where: { $0.role == .drude }), dynamics.ensemble != .npt,
              (1...16_777_216).contains(maximumMeshPoints) else {
            throw VivoChemistryError.invalid("nuclear specification source, unconstrained manifold or fixed-cell contract")
        }
        let atoms = system.particles.filter { $0.role == .atom }
        guard atoms.allSatisfy({ $0.atomIndex != nil }) else { throw VivoChemistryError.invalid("nuclear atom mapping") }
        struct Identity: Encodable {
            let schema: String; let numericalContract: String; let specification: VivoNuclearMetalSpecification
        }
        let id = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            schema: "numivivo.org/nuclear-metal-complete-potential/v1;explicit-FP32-coordinate-projection",
            numericalContract: VivoMDExecutionIdentity.current,specification: self)))
        let document: VivoMolecularStructureDocument?
        switch electronic { case .classical: document = nil; case .fixed(let d,_,_), .adaptive(let d,_,_,_): document = d }
        let groups = initialState.periodicCell == nil ? [] : try VivoNuclearMolecularGroups.build(system: system,document: document)
        return try .init(hamiltonianFingerprint: id,atomIndices: atoms.map { $0.atomIndex! },particleIndices: atoms.map(\.index),
                         massesDa: atoms.map(\.massDa),periodicCell: initialState.periodicCell,periodicMoleculeGroups: groups,coordinateEvaluation: .projectedFP32)
    }
    /// Conservative reservation in addition to the kernels' own hard limits.
    /// Electronic scratch and the classical mesh are admitted together.
    public func reservationBytes(budget: VivoChemistryBudget) throws -> Int {
        _ = try definition(); try budget.validate()
        let n = system.particles.count, types = Set(system.particles.map(\.typeIdentifier)).count
        let neighbors = min(max(n-1,1),Int(dynamics.resolvedMaximumNeighborsPerParticle))
        var bytes = Double(n)*Double(2048+neighbors*16)+Double(types)*Double(types)*32
        bytes += Double(system.bonds.count+system.angles.count+system.torsions.count+system.nonbondedExceptions.count)*512
        if dynamics.electrostatics == .pme, let cell = initialState.periodicCell {
            if dynamics.pmeGridDimensions == nil {
                for v in [cell.a,cell.b,cell.c] {
                    guard v.norm/dynamics.resolvedPMEGridSpacingNM <= 65_536 else { throw VivoChemistryError.resourceLimit("nuclear PME axis budget") }
                }
            }
            // Admit the same explicit-or-spacing-derived mesh that the runtime
            // allocates, for both the point cap and combined byte reservation.
            let mesh = try VivoPMEPlan.make(cell: cell,cutoffNM: dynamics.cutoffNM,
                tolerance: dynamics.resolvedPMETolerance,targetGridSpacingNM: dynamics.resolvedPMEGridSpacingNM,
                fixedGridDimensions: dynamics.pmeGridDimensions)
            guard mesh.gridPointCount <= UInt64(maximumMeshPoints) else { throw VivoChemistryError.resourceLimit("nuclear PME mesh budget") }
            bytes += Double(mesh.gridPointCount)*256
        }
        switch electronic {
        case .classical: break
        case .fixed(let document,let plan,let config):
            try plan.validate(document: document,source: system,budget: config.budget)
            try config.periodic.validate(); try config.cluster.validate()
        case .adaptive(let document,let adaptive,_,let config):
            _ = try VivoAdaptiveQMMMPlan.compile(document: document,system: system,configuration: adaptive)
            try config.periodic.validate(); try config.cluster.validate()
            bytes *= Double(adaptive.maximumResidentPartitions+2)
        }
        switch electronic {
        case .classical: break
        case .fixed(_,_,let config), .adaptive(_,_,_,let config):
            try config.budget.validate()
            guard config.budget.maximumOperatorApplications <= budget.maximumOperatorApplications,
                  config.budget.maximumBasisFunctions <= budget.maximumBasisFunctions,
                  config.budget.maximumDeterminants <= budget.maximumDeterminants else {
                throw VivoChemistryError.resourceLimit("nuclear electronic work exceeds parent budget")
            }
            bytes += Double(config.budget.maximumBytes)
        }
        guard bytes.isFinite, bytes < Double(Int.max), bytes <= Double(budget.maximumBytes) else {
            throw VivoChemistryError.resourceLimit("combined nuclear classical/electronic allocation budget")
        }
        return Int(bytes.rounded(.up))
    }
    public func make(budget: VivoChemistryBudget = .init()) async throws -> VivoNuclearPotential {
        _ = try reservationBytes(budget: budget)
        let bound = try definition()
        let actual: VivoNuclearPotential
        switch electronic {
        case .classical:
            actual = try await .metal(system: system,initialState: initialState,configuration: dynamics)
        case .fixed(let document,let plan,let config):
            let geometry = try VivoMDCandidateGeometry(particlePositionsNM: initialState.positionsNM,periodicCell: initialState.periodicCell)
            let setup = try await VivoQMMMFreeEnergyForceFactory.make(document: document,sourceSystem: system,
                initialGeometry: geometry,dynamics: dynamics,mode: .fixed(plan: plan,configuration: config))
            let initial = VivoClassicalInitialState(systemFingerprint: try setup.retainedSystem.fingerprint(),
                positionsNM: initialState.positionsNM,periodicCell: initialState.periodicCell,sourceTimePS: initialState.sourceTimePS)
            actual = try await .metal(system: setup.retainedSystem,initialState: initial,configuration: dynamics,provider: setup.provider)
        case .adaptive(let document,let config,let basis,let electronic):
            let geometry = try VivoMDCandidateGeometry(particlePositionsNM: initialState.positionsNM,periodicCell: initialState.periodicCell)
            let setup = try await VivoQMMMFreeEnergyForceFactory.make(document: document,sourceSystem: system,
                initialGeometry: geometry,dynamics: dynamics,mode: .adaptive(configuration: config,basis: basis,electronic: electronic))
            let initial = VivoClassicalInitialState(systemFingerprint: try setup.retainedSystem.fingerprint(),
                positionsNM: initialState.positionsNM,periodicCell: initialState.periodicCell,sourceTimePS: initialState.sourceTimePS)
            actual = try await .metal(system: setup.retainedSystem,initialState: initial,configuration: dynamics,provider: setup.provider)
        }
        guard bound.atomIndices == actual.definition.atomIndices, bound.particleIndices == actual.definition.particleIndices,
              bound.massesDa == actual.definition.massesDa, bound.periodicCell == actual.definition.periodicCell,
              bound.coordinateEvaluation == actual.definition.coordinateEvaluation else {
            throw VivoChemistryError.invalid("nuclear force factory changed the declared manifold")
        }
        return try .init(definition: bound) { q in
            let result = try await actual.checked(q)
            return try .init(definition: bound,positionsNM: q,energyKJPerMol: result.energyKJPerMol,
                forcesKJPerMolNM: result.forcesKJPerMolNM,normalizedConvergenceResidual: result.normalizedConvergenceResidual)
        }
    }
}
