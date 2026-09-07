import Foundation

/// A fixed-cell alchemical proposal, not a new physical force field. Intermediate
/// U(lambda) = (1-lambda) U_A + lambda U_B blends COMPLETE endpoint Hamiltonians.
/// Electronic/polarization response is converged independently at both endpoints.
public struct VivoNCMCMetalConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/ncmc-metal-configuration/v1"
    public var schema: String = Self.schema
    public var schedule: VivoNCMCLambdaSchedule
    public var timeStepPS: Double
    public var maximumEndpointEvaluations: Int
    public var maximumAbsoluteWorkKJPerMol: Double
    public var accountingAbsoluteToleranceKJPerMol: Double
    public init(schedule: VivoNCMCLambdaSchedule, timeStepPS: Double = 0.0005,
                maximumEndpointEvaluations: Int = 100_000,
                maximumAbsoluteWorkKJPerMol: Double = 1e8,
                accountingAbsoluteToleranceKJPerMol: Double = 1e-7) {
        self.schedule = schedule; self.timeStepPS = timeStepPS
        self.maximumEndpointEvaluations = maximumEndpointEvaluations
        self.maximumAbsoluteWorkKJPerMol = maximumAbsoluteWorkKJPerMol
        self.accountingAbsoluteToleranceKJPerMol = accountingAbsoluteToleranceKJPerMol
    }
    public func validate() throws {
        try schedule.validate()
        guard schema == Self.schema, timeStepPS.isFinite, timeStepPS > 0, timeStepPS <= 0.01,
              Float(timeStepPS).isFinite, Float(timeStepPS) > 0,
              (1...10_000_000).contains(maximumEndpointEvaluations),
              maximumAbsoluteWorkKJPerMol.isFinite, maximumAbsoluteWorkKJPerMol > 0,
              maximumAbsoluteWorkKJPerMol <= 1e12,
              accountingAbsoluteToleranceKJPerMol.isFinite,
              accountingAbsoluteToleranceKJPerMol > 0, accountingAbsoluteToleranceKJPerMol <= 1e-3 else {
            throw VivoChemistryError.invalid("NCMC Metal configuration")
        }
        // No reliance on cache hits for admission. Two endpoint evaluations per
        // probe; each midpoint MD segment requires two force probes per step.
        let intervals = schedule.lambdas.count - 1
        let perInterval = 4 * Int(schedule.propagationStepsPerLambda) + 16
        let (required, overflow) = intervals.multipliedReportingOverflow(by: perInterval)
        guard !overflow, required <= maximumEndpointEvaluations else {
            throw VivoChemistryError.resourceLimit("NCMC endpoint-evaluation admission exceeds declared budget")
        }
    }
}

public struct VivoConstantPHNCMCMetalSetup: Sendable {
    public let executableStates: [VivoConstantPHExecutableState]
    public let switchEngine: VivoConstantPHNCMCSwitchEngine
    public let configuration: VivoNCMCMetalConfiguration
}

/// Native switching uses the EXISTING velocity-Verlet/RATTLE Metal runtime. It
/// owns no replacement bonded, nonbonded, QM, polarization or site-force code.
public enum VivoConstantPHNCMCMetalFactory {
    public static let interpretation = "Complete-endpoint linear-Hamiltonian NCMC using native Metal NVE propagation and symmetric midpoint switching. Fixed-lambda work is measured as delta(U+K), including integration error; no stochastic heat/path ratio is fabricated. Fixed masses, constraints, dependent sites and cell are required. Float32 trajectory arithmetic and constrained solver tolerances still require numerical reversibility/convergence qualification. Linear interpolation is not a soft-core alchemical model or pKa calibration."

    public static func make(specifications: [VivoConstantPHMetalStateSpecification],
                            samplingTemperatureK: Double,
                            configuration: VivoNCMCMetalConfiguration) throws -> VivoConstantPHNCMCMetalSetup {
        try configuration.validate()
        let states = try VivoConstantPHMetalStateFactory.make(specifications: specifications,
                                                              samplingTemperatureK: samplingTemperatureK)
        guard states.count >= 2, let first = specifications.first, let manifold = states.first?.physicalManifoldFingerprint else {
            throw VivoChemistryError.invalid("NCMC requires at least two endpoint Hamiltonians")
        }
        guard specifications.allSatisfy({ !$0.system.particles.contains(where: { $0.role == .drude }) }) else {
            throw VivoChemistryError.unsupported("NCMC massive Drude variables require their own phase-space contract")
        }
        let byID = Dictionary(uniqueKeysWithValues: specifications.map { ($0.identifier, $0) })
        let stateIDs = Dictionary(uniqueKeysWithValues: states.map { ($0.identifier, $0.hamiltonianFingerprint) })
        var switching = first.configuration
        switching.timeStepPS = configuration.timeStepPS
        switching.ensemble = .nve; switching.thermostat = .none; switching.targetTemperatureK = nil
        switching.frictionPerPS = nil; switching.barostat = .none; switching.targetPressureBar = nil
        try switching.validate()
        let switchingConfiguration = switching
        let carriers = first.system.particles.map { p in
            VivoClassicalParticle(index: p.index, atomIndex: p.atomIndex, typeIdentifier: "NCMC-carrier",
                                  role: p.role, massDa: p.massDa, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
        }
        // This carrier owns only the invariant mass/constraint/site manifold.
        // ALL physical potential terms are supplied by the two endpoint evaluators.
        let carrier = VivoClassicalSystem(identifier: "ncmc-common-manifold", structureFingerprint: first.system.structureFingerprint,
            particles: carriers, constraints: first.system.constraints,
            linearVirtualSites: first.system.linearVirtualSites, virtualSiteDefinitions: first.system.virtualSiteDefinitions)
        try VivoClassicalSystemValidator.validate(carrier)
        struct Identity: Codable {
            let schema: String; let endpoints: [String: VivoFingerprint]; let manifold: VivoFingerprint
            let settings: VivoNCMCMetalConfiguration; let switching: VivoMDConfiguration
        }
        let engineID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Identity(
            schema: "numivivo.org/ncmc-metal-engine/v1", endpoints: stateIDs, manifold: manifold,
            settings: configuration, switching: switchingConfiguration)))
        let engine = VivoConstantPHNCMCSwitchEngine(fingerprint: engineID, physicalManifoldFingerprint: manifold) { initial, from, to in
            guard from.identifier != to.identifier,
                  from.physicalManifoldFingerprint == manifold, to.physicalManifoldFingerprint == manifold,
                  stateIDs[from.identifier] == from.hamiltonianFingerprint,
                  stateIDs[to.identifier] == to.hamiltonianFingerprint,
                  let a = byID[from.identifier], let b = byID[to.identifier] else {
                throw VivoChemistryError.invalid("NCMC endpoint transplant or physical-manifold mismatch")
            }
            try validateState(initial, system: carrier)
            let left = try await restore(a.system, switchingConfiguration, a.forceProvider, initial)
            let right = try await restore(b.system, switchingConfiguration, b.forceProvider, initial)
            let endpoints = VivoNCMCMetalEndpoints(left: left, right: right,
                maximumEvaluations: configuration.maximumEndpointEvaluations)
            struct PairIdentity: Codable { let engine: VivoFingerprint; let from: String; let to: String }
            let pairID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(PairIdentity(
                engine: engineID, from: from.identifier, to: to.identifier)))
            let lambdaHamiltonian = VivoNCMCLambdaHamiltonian(fingerprint: pairID,
                physicalManifoldFingerprint: manifold, fromStateIdentifier: from.identifier, toStateIdentifier: to.identifier,
                completePotentialEnergyKJPerMol: { state, lambda in
                    try validateState(state, system: carrier)
                    return try await endpoints.evaluate(geometry(state), lambda: lambda).energyKJPerMol
                }, propagate: { state, lambda, steps in
                    try validateState(state, system: carrier)
                    struct StageIdentity: Codable { let pair: VivoFingerprint; let lambda: Double }
                    let stageID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(StageIdentity(pair: pairID, lambda: lambda)))
                    let carrierID = try carrier.fingerprint()
                    let provider = try VivoMDCandidateForceProvider(fingerprint: stageID,
                        retainedSystemFingerprint: carrierID,
                        boundary: state.periodicCell == nil ? .finiteCluster : .periodicElectrostatic,
                        supportsCellMoves: false, maximumAcceptedResidual: 1, molecularConnectivitySystem: carrier) { candidate in
                            let result = try await endpoints.evaluate(candidate, lambda: lambda)
                            return try .init(providerFingerprint: stageID, geometry: candidate,
                                additionalEnergyKJPerMol: result.energyKJPerMol,
                                physicalParticleForcesKJPerMolNM: result.physicalParticleForcesKJPerMolNM,
                                derivativeMethod: "complete-endpoint linear Hamiltonian; analytic endpoint forces",
                                convergenceResidual: result.normalizedForceResidual ?? 0, requiredResidual: 1)
                        }
                    let runtime = try await restore(carrier, switchingConfiguration, provider, state)
                    let before = try await physical(runtime)
                    try samePhysicalCoordinatesAndMomenta(before, state, system: carrier)
                    let u0 = try await endpoints.evaluate(geometry(before), lambda: lambda).energyKJPerMol
                    let k0 = try kineticEnergy(before, system: carrier)
                    for _ in 0..<steps {
                        try Task.checkCancellation()
                        let certificate = try await runtime.step()
                        guard certificate.committed else { throw VivoMDRuntimeError.candidateRejected(certificate.statusFlags) }
                    }
                    let after = try await physical(runtime)
                    let u1 = try await endpoints.evaluate(geometry(after), lambda: lambda).energyKJPerMol
                    let k1 = try kineticEnergy(after, system: carrier)
                    return .init(finalPhysicalState: after, shadowWorkKJPerMol: (u1 - u0) + (k1 - k0),
                        steps: steps, method: "native Metal NVE; full fixed-lambda delta potential plus kinetic energy")
                })
            let schedule = try from.identifier < to.identifier ? configuration.schedule : configuration.schedule.reversed()
            let u0 = try await endpoints.evaluate(geometry(initial), lambda: 0).energyKJPerMol
            let k0 = try kineticEnergy(initial, system: carrier)
            let result = try await VivoNCMCProtocolKernel.run(initial: initial, hamiltonian: lambdaHamiltonian,
                schedule: schedule, maximumAbsoluteWorkKJPerMol: configuration.maximumAbsoluteWorkKJPerMol)
            let final = result.finalPhysicalState
            let u1 = try await endpoints.evaluate(geometry(final), lambda: 1).energyKJPerMol
            let k1 = try kineticEnergy(final, system: carrier)
            let deltaH = (u1 - u0) + (k1 - k0)
            // Sum of perturbation plus deterministic shadow work must telescope.
            let tolerance = configuration.accountingAbsoluteToleranceKJPerMol
                + 128 * Double.ulpOfOne * Double(result.steps.count) * max(1, max(abs(u0), abs(u1)))
            guard deltaH.isFinite, abs(deltaH - result.totalNonequilibriumWorkKJPerMol) <= tolerance else {
                throw VivoChemistryError.convergence("NCMC complete-Hamiltonian work does not telescope")
            }
            return try .init(switchEngineFingerprint: engineID, physicalManifoldFingerprint: manifold,
                fromStateIdentifier: from.identifier, toStateIdentifier: to.identifier,
                initialPhysicalState: initial, finalPhysicalState: final,
                protocolWorkKJPerMol: result.perturbationWorkKJPerMol, shadowWorkKJPerMol: result.shadowWorkKJPerMol,
                switchingSteps: UInt64(schedule.lambdas.count - 1) * schedule.propagationStepsPerLambda,
                interpretation: interpretation)
        }
        return .init(executableStates: states, switchEngine: engine, configuration: configuration)
    }

    private static func geometry(_ state: VivoConstantPHPhysicalState) throws -> VivoMDCandidateGeometry {
        try .init(particlePositionsNM: state.positionsNM, periodicCell: state.periodicCell)
    }
    private static func physical(_ runtime: VivoMDMetalRuntime) async throws -> VivoConstantPHPhysicalState {
        let snapshot = try await runtime.snapshot()
        return try .init(stepIndex: snapshot.stepIndex, timePS: snapshot.timePS, positionsNM: snapshot.positionsNM,
                         velocitiesNMPerPS: snapshot.velocitiesNMPerPS, periodicCell: snapshot.periodicCell)
    }
    private static func restore(_ system: VivoClassicalSystem, _ configuration: VivoMDConfiguration,
                                _ provider: VivoMDCandidateForceProvider?, _ state: VivoConstantPHPhysicalState) async throws -> VivoMDMetalRuntime {
        let checkpoint = VivoMDCheckpoint(systemFingerprint: try system.fingerprint(),
            configurationFingerprint: try VivoMDCandidateForceProvider.executionFingerprint(configuration: configuration, provider: provider),
            acceptedStep: state.stepIndex, timePS: state.timePS, positionsNM: state.positionsNM,
            velocitiesNMPerPS: state.velocitiesNMPerPS, periodicCell: state.periodicCell)
        return try await VivoMDMetalRuntime.restore(system: system, configuration: configuration, checkpoint: checkpoint, forceProvider: provider)
    }
    private static func validateState(_ state: VivoConstantPHPhysicalState, system: VivoClassicalSystem) throws {
        try state.validate()
        guard state.positionsNM.count == system.particles.count else { throw VivoChemistryError.invalid("NCMC particle shape") }
        for particle in system.particles where particle.role != .virtualSite {
            let i = Int(particle.index), p = state.positionsNM[i], v = state.velocitiesNMPerPS[i]
            guard [p.x, p.y, p.z, v.x, v.y, v.z].allSatisfy({ Double(Float($0)) == $0 }) else {
                throw VivoChemistryError.invalid("NCMC Metal requires explicit FP32 physical coordinates and velocities; unaccounted rounding is not a proposal")
            }
        }
    }
    private static func samePhysicalCoordinatesAndMomenta(_ a: VivoConstantPHPhysicalState, _ b: VivoConstantPHPhysicalState,
                                                         system: VivoClassicalSystem) throws {
        for p in system.particles where p.role != .virtualSite {
            let i = Int(p.index)
            guard a.positionsNM[i] == b.positionsNM[i], a.velocitiesNMPerPS[i] == b.velocitiesNMPerPS[i] else {
                throw VivoChemistryError.invalid("NCMC restore changed the independent phase-space coordinates")
            }
        }
        guard a.periodicCell == b.periodicCell else { throw VivoChemistryError.invalid("NCMC restore changed the cell") }
    }
    private static func kineticEnergy(_ state: VivoConstantPHPhysicalState, system: VivoClassicalSystem) throws -> Double {
        // Da * (nm/ps)^2 is kJ/mol. The target uses the declared invariant masses;
        // numerical force/mass rounding is included in measured proposal work.
        var value = 0.0
        for p in system.particles where p.role != .virtualSite {
            value += 0.5 * p.massDa * state.velocitiesNMPerPS[Int(p.index)].squaredNorm
        }
        guard value.isFinite else { throw VivoChemistryError.convergence("NCMC kinetic energy overflow") }
        return value
    }
}

/// At most two endpoint runtimes and one active midpoint runtime are resident.
/// Cache only an exact geometry; changing lambda reuses converged endpoint
/// energies/forces, never a stale density at a different conformation.
private actor VivoNCMCMetalEndpoints {
    let left: VivoMDMetalRuntime
    let right: VivoMDMetalRuntime
    let maximumEvaluations: Int
    private var evaluations = 0
    private var cachedGeometry: VivoMDCandidateGeometry?
    private var cached: (VivoMDHamiltonianEvaluation, VivoMDHamiltonianEvaluation)?
    private var inFlight = false
    init(left: VivoMDMetalRuntime, right: VivoMDMetalRuntime, maximumEvaluations: Int) {
        self.left = left; self.right = right; self.maximumEvaluations = maximumEvaluations
    }
    func evaluate(_ geometry: VivoMDCandidateGeometry, lambda: Double) async throws -> VivoMDHamiltonianEvaluation {
        guard !inFlight else { throw VivoChemistryError.invalid("NCMC endpoint evaluator already in flight") }
        inFlight = true; defer { inFlight = false }
        guard lambda.isFinite, (0...1).contains(lambda) else { throw VivoChemistryError.invalid("NCMC lambda outside [0,1]") }
        let a: VivoMDHamiltonianEvaluation, b: VivoMDHamiltonianEvaluation
        if cachedGeometry == geometry, let cached {
            (a, b) = cached
        } else {
            guard evaluations <= maximumEvaluations - 2 else { throw VivoChemistryError.resourceLimit("NCMC endpoint-evaluation budget") }
            evaluations += 2
            a = try await left.evaluateHamiltonian(at: geometry)
            b = try await right.evaluateHamiltonian(at: geometry)
            guard a.evaluatedGeometry.periodicCell == geometry.periodicCell, b.evaluatedGeometry.periodicCell == geometry.periodicCell else {
                throw VivoChemistryError.invalid("NCMC endpoint probe changed periodic geometry")
            }
            for p in left.system.particles where p.role != .virtualSite {
                let i = Int(p.index)
                guard a.evaluatedGeometry.particlePositionsNM[i] == geometry.particlePositionsNM[i],
                      b.evaluatedGeometry.particlePositionsNM[i] == geometry.particlePositionsNM[i] else {
                    throw VivoChemistryError.invalid("NCMC endpoint probe changed independent coordinates")
                }
            }
            cachedGeometry = geometry; cached = (a, b)
        }
        if lambda == 0 { return a }
        if lambda == 1 { return b }
        let energy = a.energyKJPerMol + lambda * (b.energyKJPerMol - a.energyKJPerMol)
        let forces = zip(a.physicalParticleForcesKJPerMolNM, b.physicalParticleForcesKJPerMolNM).map { $0.0 + ($0.1 - $0.0) * lambda }
        guard energy.isFinite, forces.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("NCMC Hamiltonian interpolation overflow") }
        return .init(systemFingerprint: a.systemFingerprint, configurationFingerprint: a.configurationFingerprint,
            evaluatedGeometry: geometry, energyKJPerMol: energy, physicalParticleForcesKJPerMolNM: forces,
            normalizedForceResidual: max(a.normalizedForceResidual ?? 0, b.normalizedForceResidual ?? 0))
    }
}
