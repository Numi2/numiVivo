import Foundation

public struct VivoConstantPHNCMCConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-ncmc/v2"
    public var schema: String
    public var base: VivoConstantPHConfiguration
    /// Maximum absolute reported protocol/shadow work accepted as numerically sane.
    public var maximumAbsoluteWorkKJPerMol: Double
    public init(base: VivoConstantPHConfiguration,
                maximumAbsoluteWorkKJPerMol: Double = 1.0e8) {
        schema = Self.schema
        self.base = base
        self.maximumAbsoluteWorkKJPerMol = maximumAbsoluteWorkKJPerMol
    }
    public func validate() throws {
        guard schema == Self.schema,
              maximumAbsoluteWorkKJPerMol.isFinite,
              maximumAbsoluteWorkKJPerMol > 0,
              maximumAbsoluteWorkKJPerMol <= 1.0e12 else {
            throw VivoChemistryError.invalid("constant-pH NCMC configuration")
        }
        try base.validate()
    }
    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

/// One complete nonequilibrium candidate. Work is reported in the physical
/// forward direction from `fromStateIdentifier` to `toStateIdentifier`.
/// `logReverseOverForwardPathProbability` covers any nonsymmetric stochastic
/// protocol generation not already represented by the discrete neighbor proposal.
public struct VivoConstantPHNCMCSwitchResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-ncmc-switch-result/v2"
    public let schema: String
    public let switchEngineFingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let fromStateIdentifier: String
    public let toStateIdentifier: String
    public let initialPhysicalStateFingerprint: VivoFingerprint
    public let finalPhysicalState: VivoConstantPHPhysicalState
    public let finalPhysicalStateFingerprint: VivoFingerprint
    public let protocolWorkKJPerMol: Double
    public let shadowWorkKJPerMol: Double
    public let logReverseOverForwardPathProbability: Double
    public let switchingSteps: UInt64
    public let interpretation: String

    public init(switchEngineFingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                fromStateIdentifier: String,
                toStateIdentifier: String,
                initialPhysicalState: VivoConstantPHPhysicalState,
                finalPhysicalState: VivoConstantPHPhysicalState,
                protocolWorkKJPerMol: Double,
                shadowWorkKJPerMol: Double,
                logReverseOverForwardPathProbability: Double = 0,
                switchingSteps: UInt64,
                interpretation: String) throws {
        schema = Self.schema
        self.switchEngineFingerprint = switchEngineFingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.fromStateIdentifier = fromStateIdentifier
        self.toStateIdentifier = toStateIdentifier
        self.initialPhysicalStateFingerprint = try initialPhysicalState.fingerprint()
        self.finalPhysicalState = finalPhysicalState
        self.finalPhysicalStateFingerprint = try finalPhysicalState.fingerprint()
        self.protocolWorkKJPerMol = protocolWorkKJPerMol
        self.shadowWorkKJPerMol = shadowWorkKJPerMol
        self.logReverseOverForwardPathProbability = logReverseOverForwardPathProbability
        self.switchingSteps = switchingSteps
        self.interpretation = interpretation
    }

    public var totalNonequilibriumWorkKJPerMol: Double {
        protocolWorkKJPerMol + shadowWorkKJPerMol
    }

    public func validate(engine: VivoConstantPHNCMCSwitchEngine,
                         initial: VivoConstantPHPhysicalState,
                         from: VivoConstantPHExecutableState,
                         to: VivoConstantPHExecutableState,
                         maximumAbsoluteWorkKJPerMol: Double) throws {
        try initial.validate(); try finalPhysicalState.validate()
        guard schema == Self.schema,
              switchEngineFingerprint == engine.fingerprint,
              physicalManifoldFingerprint == engine.physicalManifoldFingerprint,
              physicalManifoldFingerprint == from.physicalManifoldFingerprint,
              physicalManifoldFingerprint == to.physicalManifoldFingerprint,
              fromStateIdentifier == from.identifier,
              toStateIdentifier == to.identifier,
              initialPhysicalStateFingerprint == (try initial.fingerprint()),
              finalPhysicalStateFingerprint == (try finalPhysicalState.fingerprint()),
              finalPhysicalState.positionsNM.count == initial.positionsNM.count,
              finalPhysicalState.velocitiesNMPerPS.count == initial.velocitiesNMPerPS.count,
              switchingSteps > 0,
              initial.stepIndex <= UInt64.max - switchingSteps,
              finalPhysicalState.stepIndex == initial.stepIndex + switchingSteps,
              finalPhysicalState.timePS > initial.timePS,
              finalPhysicalState.periodicCell == initial.periodicCell,
              protocolWorkKJPerMol.isFinite,
              shadowWorkKJPerMol.isFinite,
              logReverseOverForwardPathProbability.isFinite,
              abs(protocolWorkKJPerMol) <= maximumAbsoluteWorkKJPerMol,
              abs(shadowWorkKJPerMol) <= maximumAbsoluteWorkKJPerMol,
              totalNonequilibriumWorkKJPerMol.isFinite,
              !interpretation.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
            throw VivoChemistryError.invalid("constant-pH NCMC switch identity, work or physical manifold")
        }
    }
}

/// Executable NCMC switching authority. The closure owns interpolation schedules,
/// force evaluation, stochastic propagation, protocol work and integrator/shadow
/// work. Endpoint potential differences are insufficient and are not accepted as
/// a substitute for the exact reported nonequilibrium work.
public struct VivoConstantPHNCMCSwitchEngine: Sendable {
    public let fingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let switchCandidate: @Sendable (
        _ initial: VivoConstantPHPhysicalState,
        _ from: VivoConstantPHExecutableState,
        _ to: VivoConstantPHExecutableState
    ) async throws -> VivoConstantPHNCMCSwitchResult

    public init(fingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                switchCandidate: @escaping @Sendable (
                    VivoConstantPHPhysicalState,
                    VivoConstantPHExecutableState,
                    VivoConstantPHExecutableState
                ) async throws -> VivoConstantPHNCMCSwitchResult) {
        self.fingerprint = fingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.switchCandidate = switchCandidate
    }
}

public struct VivoConstantPHNCMCAttempt: Codable, Sendable, Equatable {
    public let attemptIndex: Int
    public let fromStateIdentifier: String
    public let proposedStateIdentifier: String
    public let preSwitchPhysicalStateFingerprint: VivoFingerprint
    public let proposedPhysicalStateFingerprint: VivoFingerprint
    public let fromEndpointPotentialEnergyKJPerMol: Double
    public let proposedEndpointPotentialEnergyKJPerMol: Double
    public let endpointPotentialEnergyDifferenceKJPerMol: Double
    public let protocolWorkKJPerMol: Double
    public let shadowWorkKJPerMol: Double
    public let totalNonequilibriumWorkKJPerMol: Double
    public let semigrandReservoirDifferenceKJPerMol: Double
    public let logProposalRatio: Double
    public let logPathProbabilityRatio: Double
    public let logAcceptanceProbability: Double
    public let logUniform: Double
    public let switchingSteps: UInt64
    public let accepted: Bool
}

public struct VivoConstantPHNCMCCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-ncmc-checkpoint/v2"
    public var schema: String
    public let configurationFingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let executableHamiltonianFingerprints: [String: VivoFingerprint]
    public let switchEngineFingerprint: VivoFingerprint
    public let completedAttempts: Int
    public let currentStateIdentifier: String
    public let physicalState: VivoConstantPHPhysicalState
    public let randomState: UInt64
    public let visitCounts: [String: Int]
    public let acceptedMoves: Int
    public let attempts: [VivoConstantPHNCMCAttempt]

    public init(configurationFingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                executableHamiltonianFingerprints: [String: VivoFingerprint],
                switchEngineFingerprint: VivoFingerprint,
                completedAttempts: Int,
                currentStateIdentifier: String,
                physicalState: VivoConstantPHPhysicalState,
                randomState: UInt64,
                visitCounts: [String: Int],
                acceptedMoves: Int,
                attempts: [VivoConstantPHNCMCAttempt]) {
        schema = Self.schema
        self.configurationFingerprint = configurationFingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.executableHamiltonianFingerprints = executableHamiltonianFingerprints
        self.switchEngineFingerprint = switchEngineFingerprint
        self.completedAttempts = completedAttempts
        self.currentStateIdentifier = currentStateIdentifier
        self.physicalState = physicalState
        self.randomState = randomState
        self.visitCounts = visitCounts
        self.acceptedMoves = acceptedMoves
        self.attempts = attempts
    }
}

public struct VivoConstantPHNCMCResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-ncmc-result/v2"
    public let schema: String
    public let configurationFingerprint: VivoFingerprint
    public let finalCheckpoint: VivoConstantPHNCMCCheckpoint
    public let populations: [String: Double]
    public let acceptanceFraction: Double
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

public typealias VivoConstantPHNCMCProgressSink = @Sendable (VivoConstantPHNCMCCheckpoint) async throws -> Void

/// Alternating fixed-state MD and exact-work NCMC state proposals on one common
/// particle/mass manifold. Rejected candidates discard the complete NCMC endpoint
/// and retain the pre-switch propagated state; accepted candidates commit the
/// NCMC endpoint atomically with the proposed chemical state.
public actor VivoConstantPHNCMC {
    private static let gasConstantKJ = 0.00831446261815324
    public static let interpretation = "Discrete constant-pH NCMC over a fixed particle/mass manifold. Each state proposal uses a complete switching-engine nonequilibrium work (protocol plus shadow/integration work) and any explicit reverse/forward path-probability correction. Endpoint potential-energy differences are retained only as diagnostics. Rejected switching trajectories restore pre-switch coordinates and reverse momenta; failed evaluations do not commit an attempt."

    public let configuration: VivoConstantPHNCMCConfiguration
    private let runtime: [String: VivoConstantPHExecutableState]
    private let definitions: [String: VivoConstantPHStateDefinition]
    private let manifold: VivoFingerprint
    private let hamiltonians: [String: VivoFingerprint]
    private let switchEngine: VivoConstantPHNCMCSwitchEngine
    private var committed: VivoConstantPHNCMCCheckpoint
    private var inFlight = false

    public init(configuration: VivoConstantPHNCMCConfiguration,
                initialPhysicalState: VivoConstantPHPhysicalState,
                executableStates: [VivoConstantPHExecutableState],
                switchEngine: VivoConstantPHNCMCSwitchEngine,
                checkpoint: VivoConstantPHNCMCCheckpoint? = nil) throws {
        try configuration.validate(); try initialPhysicalState.validate()
        let base = configuration.base
        guard executableStates.count == base.states.count,
              Set(executableStates.map(\.identifier)) == Set(base.states.map(\.identifier)),
              Set(executableStates.map(\.physicalManifoldFingerprint)).count == 1,
              let manifold = executableStates.first?.physicalManifoldFingerprint,
              switchEngine.physicalManifoldFingerprint == manifold else {
            throw VivoChemistryError.invalid("constant-pH NCMC executable state set, switch engine or physical manifold")
        }
        let runtime = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0) })
        let hamiltonians = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0.hamiltonianFingerprint) })
        let definitions = Dictionary(uniqueKeysWithValues: base.states.map { ($0.identifier, $0) })
        let configurationID = try configuration.fingerprint()
        if let checkpoint {
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians,
                              switchEngineFingerprint: switchEngine.fingerprint)
            committed = checkpoint
        } else {
            var counts = Dictionary(uniqueKeysWithValues: base.states.map { ($0.identifier, 0) })
            counts[base.initialStateIdentifier] = 1
            committed = .init(configurationFingerprint: configurationID,
                physicalManifoldFingerprint: manifold,
                executableHamiltonianFingerprints: hamiltonians,
                switchEngineFingerprint: switchEngine.fingerprint,
                completedAttempts: 0,
                currentStateIdentifier: base.initialStateIdentifier,
                physicalState: initialPhysicalState,
                randomState: base.seed,
                visitCounts: counts,
                acceptedMoves: 0,
                attempts: [])
        }
        self.configuration = configuration
        self.runtime = runtime
        self.definitions = definitions
        self.manifold = manifold
        self.hamiltonians = hamiltonians
        self.switchEngine = switchEngine
    }

    public func checkpoint() -> VivoConstantPHNCMCCheckpoint { committed }

    public func run(progress: VivoConstantPHNCMCProgressSink? = nil) async throws -> VivoConstantPHNCMCResult {
        guard !inFlight else { throw VivoChemistryError.invalid("constant-pH NCMC run already in flight") }
        inFlight = true
        defer { inFlight = false }
        let base = configuration.base
        var checkpoint = committed
        var random = VivoSplitMix64(state: checkpoint.randomState)
        let rt = Self.gasConstantKJ * base.temperatureK
        let protonSlope = rt * log(10.0) * (base.targetPH - base.referencePH)

        while checkpoint.completedAttempts < base.attemptCount {
            try Task.checkCancellation()
            guard let currentDefinition = definitions[checkpoint.currentStateIdentifier],
                  let currentRuntime = runtime[checkpoint.currentStateIdentifier] else {
                throw VivoChemistryError.invalid("constant-pH NCMC current state disappeared")
            }
            let propagated = try await currentRuntime.propagate(checkpoint.physicalState, base.mdStepsPerAttempt)
            try propagated.validate()
            guard propagated.positionsNM.count == checkpoint.physicalState.positionsNM.count,
                  checkpoint.physicalState.stepIndex <= UInt64.max - base.mdStepsPerAttempt,
                  propagated.stepIndex == checkpoint.physicalState.stepIndex + base.mdStepsPerAttempt,
                  propagated.timePS > checkpoint.physicalState.timePS,
                  propagated.periodicCell == checkpoint.physicalState.periodicCell else {
                throw VivoChemistryError.invalid("constant-pH NCMC propagator changed physical manifold or reversed the MD clock")
            }
            let neighbors = currentDefinition.neighbors
            let selector = min(Int(random.unitInterval() * Double(neighbors.count)), neighbors.count - 1)
            let proposedIdentifier = neighbors[selector]
            guard let proposedDefinition = definitions[proposedIdentifier],
                  let proposedRuntime = runtime[proposedIdentifier] else {
                throw VivoChemistryError.invalid("constant-pH NCMC proposal state disappeared")
            }
            let switched = try await switchEngine.switchCandidate(propagated, currentRuntime, proposedRuntime)
            try switched.validate(engine: switchEngine, initial: propagated,
                                  from: currentRuntime, to: proposedRuntime,
                                  maximumAbsoluteWorkKJPerMol: configuration.maximumAbsoluteWorkKJPerMol)

            async let fromEndpointValue = currentRuntime.potentialEnergyKJPerMol(propagated)
            async let proposedEndpointValue = proposedRuntime.potentialEnergyKJPerMol(switched.finalPhysicalState)
            let fromEndpoint = try await fromEndpointValue
            let proposedEndpoint = try await proposedEndpointValue
            guard fromEndpoint.isFinite, proposedEndpoint.isFinite else {
                throw VivoChemistryError.convergence("nonfinite constant-pH NCMC endpoint energy")
            }

            let fromReservoir = currentDefinition.referenceSemigrandBiasKJPerMol
                + Double(currentDefinition.boundProtonOffset) * protonSlope
            let proposedReservoir = proposedDefinition.referenceSemigrandBiasKJPerMol
                + Double(proposedDefinition.boundProtonOffset) * protonSlope
            let reservoirDelta = proposedReservoir - fromReservoir
            let logProposalRatio = log(Double(currentDefinition.neighbors.count) / Double(proposedDefinition.neighbors.count))
            let logPathRatio = switched.logReverseOverForwardPathProbability
            let exponent = -(switched.totalNonequilibriumWorkKJPerMol + reservoirDelta) / rt
                + logProposalRatio + logPathRatio
            guard exponent.isFinite else {
                throw VivoChemistryError.convergence("constant-pH NCMC acceptance exponent overflow")
            }
            let logAcceptance = min(0.0, exponent)
            let logUniform = log(max(random.unitInterval(), Double.leastNonzeroMagnitude))
            let accepted = logUniform < logAcceptance
            let nextState = accepted ? proposedIdentifier : checkpoint.currentStateIdentifier
            // Rejection restores the pre-switch coordinates and reverses momentum.
            // Returning the unchanged forward momenta breaks the NCMC invariant
            // measure for deterministic propagation without a fresh Maxwell draw.
            var nextPhysical = accepted ? switched.finalPhysicalState : propagated
            if !accepted { nextPhysical.velocitiesNMPerPS = propagated.velocitiesNMPerPS.map { $0 * -1 } }
            var counts = checkpoint.visitCounts
            counts[nextState, default: 0] += 1
            let attempt = VivoConstantPHNCMCAttempt(
                attemptIndex: checkpoint.completedAttempts,
                fromStateIdentifier: checkpoint.currentStateIdentifier,
                proposedStateIdentifier: proposedIdentifier,
                preSwitchPhysicalStateFingerprint: try propagated.fingerprint(),
                proposedPhysicalStateFingerprint: try switched.finalPhysicalState.fingerprint(),
                fromEndpointPotentialEnergyKJPerMol: fromEndpoint,
                proposedEndpointPotentialEnergyKJPerMol: proposedEndpoint,
                endpointPotentialEnergyDifferenceKJPerMol: proposedEndpoint - fromEndpoint,
                protocolWorkKJPerMol: switched.protocolWorkKJPerMol,
                shadowWorkKJPerMol: switched.shadowWorkKJPerMol,
                totalNonequilibriumWorkKJPerMol: switched.totalNonequilibriumWorkKJPerMol,
                semigrandReservoirDifferenceKJPerMol: reservoirDelta,
                logProposalRatio: logProposalRatio,
                logPathProbabilityRatio: logPathRatio,
                logAcceptanceProbability: logAcceptance,
                logUniform: logUniform,
                switchingSteps: switched.switchingSteps,
                accepted: accepted)
            var attempts = checkpoint.attempts
            attempts.append(attempt)
            checkpoint = .init(configurationFingerprint: checkpoint.configurationFingerprint,
                physicalManifoldFingerprint: manifold,
                executableHamiltonianFingerprints: hamiltonians,
                switchEngineFingerprint: switchEngine.fingerprint,
                completedAttempts: checkpoint.completedAttempts + 1,
                currentStateIdentifier: nextState,
                physicalState: nextPhysical,
                randomState: random.state,
                visitCounts: counts,
                acceptedMoves: checkpoint.acceptedMoves + (accepted ? 1 : 0),
                attempts: attempts)
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians,
                              switchEngineFingerprint: switchEngine.fingerprint)
            try Task.checkCancellation()
            committed = checkpoint
            if let progress { try await progress(checkpoint) }
        }

        let totalVisits = checkpoint.visitCounts.values.reduce(0, +)
        let populations = Dictionary(uniqueKeysWithValues: base.states.map { state in
            (state.identifier, totalVisits > 0 ? Double(checkpoint.visitCounts[state.identifier, default: 0]) / Double(totalVisits) : 0)
        })
        let result = VivoConstantPHNCMCResult(schema: VivoConstantPHNCMCResult.schema,
            configurationFingerprint: checkpoint.configurationFingerprint,
            finalCheckpoint: checkpoint,
            populations: populations,
            acceptanceFraction: checkpoint.completedAttempts > 0 ? Double(checkpoint.acceptedMoves) / Double(checkpoint.completedAttempts) : 0,
            interpretation: Self.interpretation,
            evidenceFingerprint: try evidenceFingerprint(checkpoint: checkpoint, populations: populations))
        return result
    }

    private func evidenceFingerprint(checkpoint: VivoConstantPHNCMCCheckpoint,
                                     populations: [String: Double]) throws -> VivoFingerprint {
        struct Evidence: Codable {
            let schema: String
            let configuration: VivoConstantPHNCMCConfiguration
            let checkpoint: VivoConstantPHNCMCCheckpoint
            let populations: [String: Double]
        }
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/constant-ph-ncmc-evidence/v2",
            configuration: configuration,
            checkpoint: checkpoint,
            populations: populations)))
    }

    public static func validate(checkpoint: VivoConstantPHNCMCCheckpoint,
                                configuration: VivoConstantPHNCMCConfiguration,
                                manifold: VivoFingerprint,
                                hamiltonians: [String: VivoFingerprint],
                                switchEngineFingerprint: VivoFingerprint) throws {
        try configuration.validate(); try checkpoint.physicalState.validate()
        let base = configuration.base
        guard checkpoint.schema == VivoConstantPHNCMCCheckpoint.schema,
              checkpoint.configurationFingerprint == (try configuration.fingerprint()),
              checkpoint.physicalManifoldFingerprint == manifold,
              checkpoint.executableHamiltonianFingerprints == hamiltonians,
              checkpoint.switchEngineFingerprint == switchEngineFingerprint,
              (0...base.attemptCount).contains(checkpoint.completedAttempts),
              base.states.contains(where: { $0.identifier == checkpoint.currentStateIdentifier }),
              checkpoint.visitCounts.keys.sorted() == base.states.map(\.identifier).sorted(),
              checkpoint.visitCounts.values.allSatisfy({ $0 >= 0 && $0 <= base.attemptCount + 1 }),
              checkpoint.visitCounts.values.reduce(0, +) == checkpoint.completedAttempts + 1,
              checkpoint.acceptedMoves >= 0,
              checkpoint.acceptedMoves <= checkpoint.completedAttempts,
              checkpoint.attempts.count == checkpoint.completedAttempts else {
            throw VivoChemistryError.invalid("constant-pH NCMC checkpoint identity or counters")
        }
        let definitions = Dictionary(uniqueKeysWithValues: base.states.map { ($0.identifier, $0) })
        var random = VivoSplitMix64(state: base.seed)
        var current = base.initialStateIdentifier
        var counts = Dictionary(uniqueKeysWithValues: base.states.map { ($0.identifier, 0) })
        counts[current] = 1
        var acceptedCount = 0
        let rt = Self.gasConstantKJ * base.temperatureK
        let protonSlope = rt * log(10) * (base.targetPH - base.referencePH)
        func agrees(_ a: Double, _ b: Double) -> Bool {
            a.isFinite && b.isFinite && abs(a - b) <= 1e-10 * max(1, abs(a), abs(b))
        }
        for (index, attempt) in checkpoint.attempts.enumerated() {
            guard let from = definitions[current], attempt.fromStateIdentifier == current else {
                throw VivoChemistryError.invalid("NCMC checkpoint state history")
            }
            let selected = min(Int(random.unitInterval() * Double(from.neighbors.count)), from.neighbors.count - 1)
            guard let to = definitions[from.neighbors[selected]], attempt.proposedStateIdentifier == to.identifier else {
                throw VivoChemistryError.invalid("NCMC checkpoint proposal/RNG history")
            }
            let logUniform = log(max(random.unitInterval(), Double.leastNonzeroMagnitude))
            let reservoir = (to.referenceSemigrandBiasKJPerMol + Double(to.boundProtonOffset) * protonSlope)
                - (from.referenceSemigrandBiasKJPerMol + Double(from.boundProtonOffset) * protonSlope)
            let proposal = log(Double(from.neighbors.count) / Double(to.neighbors.count))
            let acceptance = min(0, -(attempt.totalNonequilibriumWorkKJPerMol + reservoir) / rt
                + proposal + attempt.logPathProbabilityRatio)
            guard agrees(attempt.logUniform, logUniform), agrees(attempt.semigrandReservoirDifferenceKJPerMol, reservoir),
                  agrees(attempt.logProposalRatio, proposal), agrees(attempt.logAcceptanceProbability, acceptance),
                  attempt.accepted == (logUniform < acceptance),
                  agrees(attempt.endpointPotentialEnergyDifferenceKJPerMol,
                         attempt.proposedEndpointPotentialEnergyKJPerMol - attempt.fromEndpointPotentialEnergyKJPerMol),
                  abs(attempt.protocolWorkKJPerMol) <= configuration.maximumAbsoluteWorkKJPerMol,
                  abs(attempt.shadowWorkKJPerMol) <= configuration.maximumAbsoluteWorkKJPerMol else {
                throw VivoChemistryError.invalid("NCMC checkpoint acceptance/work history")
            }
            if attempt.accepted { current = to.identifier; acceptedCount += 1 }
            counts[current, default: 0] += 1
            guard attempt.attemptIndex == index,
                  attempt.fromEndpointPotentialEnergyKJPerMol.isFinite,
                  attempt.proposedEndpointPotentialEnergyKJPerMol.isFinite,
                  attempt.endpointPotentialEnergyDifferenceKJPerMol.isFinite,
                  attempt.protocolWorkKJPerMol.isFinite,
                  attempt.shadowWorkKJPerMol.isFinite,
                  attempt.totalNonequilibriumWorkKJPerMol.isFinite,
                  abs(attempt.totalNonequilibriumWorkKJPerMol - (attempt.protocolWorkKJPerMol + attempt.shadowWorkKJPerMol)) <= 1e-10,
                  attempt.semigrandReservoirDifferenceKJPerMol.isFinite,
                  attempt.logProposalRatio.isFinite,
                  attempt.logPathProbabilityRatio.isFinite,
                  attempt.logAcceptanceProbability.isFinite,
                  attempt.logAcceptanceProbability <= 0,
                  attempt.logUniform.isFinite,
                  attempt.logUniform <= 0,
                  attempt.switchingSteps > 0 else {
                throw VivoChemistryError.invalid("constant-pH NCMC attempt trace")
            }
        }
        guard current == checkpoint.currentStateIdentifier, counts == checkpoint.visitCounts,
              acceptedCount == checkpoint.acceptedMoves, random.state == checkpoint.randomState else {
            throw VivoChemistryError.invalid("NCMC checkpoint final state/counters/RNG do not reconstruct")
        }
    }
}
