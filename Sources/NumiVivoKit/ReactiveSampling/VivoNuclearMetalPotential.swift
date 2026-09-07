import Foundation

public extension VivoNuclearPotential {
    /// The existing complete classical+BO force evaluator is the only authority.
    /// All independent atoms are included. Constrained quantum coordinates and
    /// massive Drude oscillators require different measures and are rejected.
    static func metal(system: VivoClassicalSystem, initialState: VivoClassicalInitialState,
                      configuration: VivoMDConfiguration,
                      provider: VivoMDCandidateForceProvider? = nil) async throws -> Self {
        try VivoClassicalSystemValidator.validate(system); try configuration.validate()
        guard system.constraints.isEmpty, !system.particles.contains(where: { $0.role == .drude }),
              configuration.ensemble != .npt else {
            throw VivoChemistryError.unsupported("nuclear/surrogate equilibrium adapter requires unconstrained atoms and a fixed cell; no massive Drude phase space")
        }
        let atoms = system.particles.filter { $0.role == .atom }
        guard atoms.allSatisfy({ $0.atomIndex != nil }), !atoms.isEmpty else {
            throw VivoChemistryError.invalid("nuclear potential requires physical chemical atom identities")
        }
        var probeConfig = configuration
        probeConfig.ensemble = .nve; probeConfig.thermostat = .none
        probeConfig.targetTemperatureK = nil; probeConfig.frictionPerPS = nil
        probeConfig.barostat = .none; probeConfig.targetPressureBar = nil
        let runtime = try await VivoMDMetalRuntime.make(system: system, initialState: initialState,
            configuration: probeConfig, forceProvider: provider)
        struct Identity: Encodable { let system: VivoFingerprint; let execution: VivoFingerprint; let projection: String }
        let id = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(system: system.fingerprint(),
            execution: VivoMDCandidateForceProvider.executionFingerprint(configuration: probeConfig, provider: provider),
            projection: "explicit-nearest-FP32-coordinate-projection;complete-Metal-plus-BO-potential/v1")))
        let definition = try VivoNuclearPotentialDefinition(hamiltonianFingerprint: id,
            atomIndices: atoms.map { $0.atomIndex! }, particleIndices: atoms.map(\.index),
            massesDa: atoms.map(\.massDa), periodicCell: initialState.periodicCell,
            periodicMoleculeGroups: initialState.periodicCell == nil ? [] : try VivoNuclearMolecularGroups.build(system: provider?.molecularConnectivitySystem ?? system), coordinateEvaluation: .projectedFP32)
        let owner = VivoNuclearMetalProbe(runtime: runtime, definition: definition,
            initialPositions: initialState.positionsNM)
        return try .init(definition: definition, evaluate: { try await owner.evaluate($0) })
    }
}
private actor VivoNuclearMetalProbe {
    let runtime: VivoMDMetalRuntime
    let definition: VivoNuclearPotentialDefinition
    let initialPositions: [VivoVector3D]
    var inFlight = false
    init(runtime: VivoMDMetalRuntime, definition: VivoNuclearPotentialDefinition, initialPositions: [VivoVector3D]) {
        self.runtime = runtime; self.definition = definition; self.initialPositions = initialPositions
    }
    func evaluate(_ positions: [VivoVector3D]) async throws -> VivoNuclearPotentialEvaluation {
        guard !inFlight else { throw VivoChemistryError.invalid("nuclear Metal probe is already in flight") }
        inFlight = true; defer { inFlight = false }
        guard positions.count == definition.atomIndices.count, positions.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("nuclear Metal coordinate layout")
        }
        var all = initialPositions
        for (i, slot) in definition.particleIndices.enumerated() {
            let p = positions[i]
            all[Int(slot)] = .init(Double(Float(p.x)),Double(Float(p.y)),Double(Float(p.z)))
        }
        guard all.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("nuclear coordinate exceeds FP32 range") }
        let geometry = try VivoMDCandidateGeometry(particlePositionsNM: all, periodicCell: definition.periodicCell)
        let evaluation = try await runtime.evaluateHamiltonian(at: geometry)
        guard let residual = evaluation.normalizedForceResidual,
              evaluation.evaluatedGeometry.periodicCell == definition.periodicCell else {
            throw VivoChemistryError.convergence("nuclear probe has no convergence/cell contract")
        }
        for slot in definition.particleIndices {
            guard evaluation.evaluatedGeometry.particlePositionsNM[Int(slot)] == all[Int(slot)] else {
                throw VivoChemistryError.invalid("nuclear Metal probe changed an independent atom image")
            }
        }
        return try .init(definition: definition, positionsNM: positions, energyKJPerMol: evaluation.energyKJPerMol,
            forcesKJPerMolNM: definition.particleIndices.map { evaluation.physicalParticleForcesKJPerMolNM[Int($0)] },
            normalizedConvergenceResidual: residual)
    }
}
