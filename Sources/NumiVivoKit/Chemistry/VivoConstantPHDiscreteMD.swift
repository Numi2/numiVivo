import Foundation

public struct VivoConstantPHStateDefinition: Codable, Sendable, Equatable {
    public var identifier: String
    /// Relative proton count. Only differences enter the semigrand weight.
    public var boundProtonOffset: Int
    /// State correction at referencePH. The executable state's energy callback
    /// returns the complete physical potential and must NOT include this term.
    public var referenceSemigrandBiasKJPerMol: Double
    public var origin: VivoKineticOrigin
    public var evidence: VivoKineticEvidence
    /// Explicit reciprocal proposal graph. Proposals are uniform over neighbors.
    public var neighbors: [String]

    public init(identifier: String, boundProtonOffset: Int,
                referenceSemigrandBiasKJPerMol: Double,
                origin: VivoKineticOrigin, evidence: VivoKineticEvidence,
                neighbors: [String]) {
        self.identifier = identifier
        self.boundProtonOffset = boundProtonOffset
        self.referenceSemigrandBiasKJPerMol = referenceSemigrandBiasKJPerMol
        self.origin = origin
        self.evidence = evidence
        self.neighbors = neighbors
    }

    public func validate() throws {
        guard !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              (-128...128).contains(boundProtonOffset),
              referenceSemigrandBiasKJPerMol.isFinite,
              !neighbors.isEmpty, neighbors.count <= 4096,
              Set(neighbors).count == neighbors.count,
              !neighbors.contains(identifier) else {
            throw VivoChemistryError.invalid("constant-pH state definition")
        }
        try evidence.validate(origin: origin)
    }
}

public struct VivoConstantPHConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-discrete-md/v1"
    public var schema: String
    public var identifier: String
    public var temperatureK: Double
    public var referencePH: Double
    public var targetPH: Double
    public var states: [VivoConstantPHStateDefinition]
    public var initialStateIdentifier: String
    public var mdStepsPerAttempt: UInt64
    public var attemptCount: Int
    public var seed: UInt64

    public init(identifier: String, temperatureK: Double, referencePH: Double,
                targetPH: Double, states: [VivoConstantPHStateDefinition],
                initialStateIdentifier: String, mdStepsPerAttempt: UInt64,
                attemptCount: Int, seed: UInt64) {
        schema = Self.schema
        self.identifier = identifier
        self.temperatureK = temperatureK
        self.referencePH = referencePH
        self.targetPH = targetPH
        self.states = states
        self.initialStateIdentifier = initialStateIdentifier
        self.mdStepsPerAttempt = mdStepsPerAttempt
        self.attemptCount = attemptCount
        self.seed = seed
    }

    public func validate() throws {
        guard schema == Self.schema,
              !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              temperatureK.isFinite, temperatureK > 0,
              referencePH.isFinite, targetPH.isFinite,
              (-10...30).contains(referencePH), (-10...30).contains(targetPH),
              !states.isEmpty, states.count <= 4096,
              Set(states.map(\.identifier)).count == states.count,
              states.contains(where: { $0.identifier == initialStateIdentifier }),
              mdStepsPerAttempt > 0, mdStepsPerAttempt <= 10_000_000,
              (1...1_000_000).contains(attemptCount) else {
            throw VivoChemistryError.invalid("constant-pH configuration")
        }
        for state in states { try state.validate() }
        let identifiers = Set(states.map(\.identifier))
        for state in states {
            guard state.neighbors.allSatisfy(identifiers.contains) else {
                throw VivoChemistryError.invalid("constant-pH proposal graph references an unknown state")
            }
            for neighbor in state.neighbors {
                guard states.first(where: { $0.identifier == neighbor })?.neighbors.contains(state.identifier) == true else {
                    throw VivoChemistryError.invalid("constant-pH proposal graph must be reciprocal")
                }
            }
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

/// State independent coordinates/momenta. All executable protonation states in
/// v1 must share this physical particle/mass manifold. This makes kinetic energy
/// cancel in a state switch and permits velocities to remain unchanged.
public struct VivoConstantPHPhysicalState: Codable, Sendable, Equatable {
    public var stepIndex: UInt64
    public var timePS: Double
    public var positionsNM: [VivoVector3D]
    public var velocitiesNMPerPS: [VivoVector3D]
    public var periodicCell: VivoPeriodicCell?

    public init(stepIndex: UInt64, timePS: Double,
                positionsNM: [VivoVector3D], velocitiesNMPerPS: [VivoVector3D],
                periodicCell: VivoPeriodicCell?) throws {
        self.stepIndex = stepIndex
        self.timePS = timePS
        self.positionsNM = positionsNM
        self.velocitiesNMPerPS = velocitiesNMPerPS
        self.periodicCell = periodicCell
        try validate()
    }

    public func validate() throws {
        guard !positionsNM.isEmpty, positionsNM.count == velocitiesNMPerPS.count,
              positionsNM.allSatisfy(\.isFinite), velocitiesNMPerPS.allSatisfy(\.isFinite),
              timePS.isFinite, timePS >= 0, periodicCell?.isValid != false else {
            throw VivoChemistryError.invalid("constant-pH physical state")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoConstantPHExecutableState: Sendable {
    public let identifier: String
    public let physicalManifoldFingerprint: VivoFingerprint
    public let hamiltonianFingerprint: VivoFingerprint
    /// Complete physical potential energy for this state at the supplied common
    /// conformation. It excludes the reference bias/proton-reservoir term.
    public let potentialEnergyKJPerMol: @Sendable (VivoConstantPHPhysicalState) async throws -> Double
    /// Propagate exactly mdSteps under this state's complete Hamiltonian.
    public let propagate: @Sendable (VivoConstantPHPhysicalState, UInt64) async throws -> VivoConstantPHPhysicalState

    public init(identifier: String, physicalManifoldFingerprint: VivoFingerprint,
                hamiltonianFingerprint: VivoFingerprint,
                potentialEnergyKJPerMol: @escaping @Sendable (VivoConstantPHPhysicalState) async throws -> Double,
                propagate: @escaping @Sendable (VivoConstantPHPhysicalState, UInt64) async throws -> VivoConstantPHPhysicalState) {
        self.identifier = identifier
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.hamiltonianFingerprint = hamiltonianFingerprint
        self.potentialEnergyKJPerMol = potentialEnergyKJPerMol
        self.propagate = propagate
    }
}

public struct VivoConstantPHAttempt: Codable, Sendable, Equatable {
    public let attemptIndex: Int
    public let physicalStateFingerprint: VivoFingerprint
    public let fromStateIdentifier: String
    public let proposedStateIdentifier: String
    public let fromPotentialEnergyKJPerMol: Double
    public let proposedPotentialEnergyKJPerMol: Double
    public let semigrandEnergyDifferenceKJPerMol: Double
    public let logProposalRatio: Double
    public let logAcceptanceProbability: Double
    public let logUniform: Double
    public let accepted: Bool
}

public struct VivoConstantPHCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-discrete-md-checkpoint/v1"
    public var schema: String
    public let configurationFingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let executableHamiltonianFingerprints: [String: VivoFingerprint]
    public let completedAttempts: Int
    public let currentStateIdentifier: String
    public let physicalState: VivoConstantPHPhysicalState
    public let randomState: UInt64
    public let visitCounts: [String: Int]
    public let acceptedMoves: Int
    public let attempts: [VivoConstantPHAttempt]

    public init(configurationFingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                executableHamiltonianFingerprints: [String: VivoFingerprint],
                completedAttempts: Int, currentStateIdentifier: String,
                physicalState: VivoConstantPHPhysicalState, randomState: UInt64,
                visitCounts: [String: Int], acceptedMoves: Int,
                attempts: [VivoConstantPHAttempt]) {
        schema = Self.schema
        self.configurationFingerprint = configurationFingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.executableHamiltonianFingerprints = executableHamiltonianFingerprints
        self.completedAttempts = completedAttempts
        self.currentStateIdentifier = currentStateIdentifier
        self.physicalState = physicalState
        self.randomState = randomState
        self.visitCounts = visitCounts
        self.acceptedMoves = acceptedMoves
        self.attempts = attempts
    }
}

public struct VivoConstantPHResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-discrete-md-result/v1"
    public let schema: String
    public let configurationFingerprint: VivoFingerprint
    public let finalCheckpoint: VivoConstantPHCheckpoint
    public let populations: [String: Double]
    public let acceptanceFraction: Double
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

public typealias VivoConstantPHProgressSink = @Sendable (VivoConstantPHCheckpoint) async throws -> Void

/// Discrete-state semigrand MD/MC. The propagation and energy functions are
/// executable runtime resources and are intentionally not serialized. This v1
/// contract requires a common particle/mass manifold across states; adding or
/// deleting physical atoms during a proposal is unsupported.
public actor VivoConstantPHDiscreteMD {
    private static let gasConstantKJ = 0.00831446261815324
    public static let interpretation = "Alternating fixed-state MD and discrete reciprocal-graph Metropolis-Hastings protonation moves in the semigrand ensemble. Current and proposed complete Hamiltonians are evaluated at the exact same conformation. State reference corrections and the proton chemical potential are explicit. State discovery, pKa prediction, atom-manifold changes and free-energy accuracy remain separate requirements."

    public let configuration: VivoConstantPHConfiguration
    private let runtime: [String: VivoConstantPHExecutableState]
    private let manifold: VivoFingerprint
    private let hamiltonians: [String: VivoFingerprint]
    private var committed: VivoConstantPHCheckpoint
    private var inFlight = false

    public init(configuration: VivoConstantPHConfiguration,
                initialPhysicalState: VivoConstantPHPhysicalState,
                executableStates: [VivoConstantPHExecutableState],
                checkpoint: VivoConstantPHCheckpoint? = nil) throws {
        try configuration.validate(); try initialPhysicalState.validate()
        guard executableStates.count == configuration.states.count,
              Set(executableStates.map(\.identifier)) == Set(configuration.states.map(\.identifier)),
              Set(executableStates.map(\.physicalManifoldFingerprint)).count == 1,
              let manifold = executableStates.first?.physicalManifoldFingerprint else {
            throw VivoChemistryError.invalid("constant-pH executable state set or physical manifold")
        }
        let runtime = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0) })
        let hamiltonians = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0.hamiltonianFingerprint) })
        let configurationID = try configuration.fingerprint()
        if let checkpoint {
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians)
            self.committed = checkpoint
        } else {
            var counts = Dictionary(uniqueKeysWithValues: configuration.states.map { ($0.identifier, 0) })
            counts[configuration.initialStateIdentifier] = 1
            self.committed = .init(configurationFingerprint: configurationID,
                physicalManifoldFingerprint: manifold,
                executableHamiltonianFingerprints: hamiltonians,
                completedAttempts: 0,
                currentStateIdentifier: configuration.initialStateIdentifier,
                physicalState: initialPhysicalState,
                randomState: configuration.seed,
                visitCounts: counts, acceptedMoves: 0, attempts: [])
        }
        self.configuration = configuration
        self.runtime = runtime
        self.manifold = manifold
        self.hamiltonians = hamiltonians
    }

    public func checkpoint() -> VivoConstantPHCheckpoint { committed }

    public func run(progress: VivoConstantPHProgressSink? = nil) async throws -> VivoConstantPHResult {
        guard !inFlight else { throw VivoChemistryError.invalid("constant-pH run already in flight") }
        inFlight = true
        defer { inFlight = false }
        let states = Dictionary(uniqueKeysWithValues: configuration.states.map { ($0.identifier, $0) })
        var checkpoint = committed
        var random = VivoSplitMix64(state: checkpoint.randomState)
        let rt = Self.gasConstantKJ * configuration.temperatureK
        let protonSlope = rt * log(10.0) * (configuration.targetPH - configuration.referencePH)

        while checkpoint.completedAttempts < configuration.attemptCount {
            try Task.checkCancellation()
            guard let currentDefinition = states[checkpoint.currentStateIdentifier],
                  let currentRuntime = runtime[checkpoint.currentStateIdentifier] else {
                throw VivoChemistryError.invalid("constant-pH current state disappeared")
            }
            let propagated = try await currentRuntime.propagate(checkpoint.physicalState, configuration.mdStepsPerAttempt)
            try propagated.validate()
            guard propagated.positionsNM.count == checkpoint.physicalState.positionsNM.count,
                  propagated.stepIndex >= checkpoint.physicalState.stepIndex,
                  propagated.timePS >= checkpoint.physicalState.timePS else {
                throw VivoChemistryError.invalid("constant-pH propagator changed physical manifold or reversed the MD clock")
            }
            let neighbors = currentDefinition.neighbors
            let selector = min(Int(random.unitInterval() * Double(neighbors.count)), neighbors.count - 1)
            let proposedIdentifier = neighbors[selector]
            guard let proposedDefinition = states[proposedIdentifier],
                  let proposedRuntime = runtime[proposedIdentifier] else {
                throw VivoChemistryError.invalid("constant-pH proposal state disappeared")
            }
            async let fromEnergyValue = currentRuntime.potentialEnergyKJPerMol(propagated)
            async let proposedEnergyValue = proposedRuntime.potentialEnergyKJPerMol(propagated)
            let fromEnergy = try await fromEnergyValue
            let proposedEnergy = try await proposedEnergyValue
            guard fromEnergy.isFinite, proposedEnergy.isFinite else {
                throw VivoChemistryError.convergence("nonfinite constant-pH state energy")
            }
            let fromReservoir = currentDefinition.referenceSemigrandBiasKJPerMol
                + Double(currentDefinition.boundProtonOffset) * protonSlope
            let proposedReservoir = proposedDefinition.referenceSemigrandBiasKJPerMol
                + Double(proposedDefinition.boundProtonOffset) * protonSlope
            let delta = (proposedEnergy + proposedReservoir) - (fromEnergy + fromReservoir)
            let logProposalRatio = log(Double(currentDefinition.neighbors.count) / Double(proposedDefinition.neighbors.count))
            let logAcceptance = min(0.0, -delta / rt + logProposalRatio)
            let logUniform = log(max(random.unitInterval(), Double.leastNonzeroMagnitude))
            let accepted = logUniform < logAcceptance
            let nextState = accepted ? proposedIdentifier : checkpoint.currentStateIdentifier
            var counts = checkpoint.visitCounts
            counts[nextState, default: 0] += 1
            let attempt = VivoConstantPHAttempt(attemptIndex: checkpoint.completedAttempts,
                physicalStateFingerprint: try propagated.fingerprint(),
                fromStateIdentifier: checkpoint.currentStateIdentifier,
                proposedStateIdentifier: proposedIdentifier,
                fromPotentialEnergyKJPerMol: fromEnergy,
                proposedPotentialEnergyKJPerMol: proposedEnergy,
                semigrandEnergyDifferenceKJPerMol: delta,
                logProposalRatio: logProposalRatio,
                logAcceptanceProbability: logAcceptance,
                logUniform: logUniform,
                accepted: accepted)
            var trace = checkpoint.attempts
            trace.append(attempt)
            checkpoint = .init(configurationFingerprint: checkpoint.configurationFingerprint,
                physicalManifoldFingerprint: manifold,
                executableHamiltonianFingerprints: hamiltonians,
                completedAttempts: checkpoint.completedAttempts + 1,
                currentStateIdentifier: nextState,
                physicalState: propagated,
                randomState: random.state,
                visitCounts: counts,
                acceptedMoves: checkpoint.acceptedMoves + (accepted ? 1 : 0),
                attempts: trace)
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians)
            committed = checkpoint
            if let progress { try await progress(checkpoint) }
        }
        let totalVisits = checkpoint.visitCounts.values.reduce(0, +)
        guard totalVisits == configuration.attemptCount + 1 else {
            throw VivoChemistryError.invalid("constant-pH visit accounting")
        }
        let populations = checkpoint.visitCounts.mapValues { Double($0) / Double(totalVisits) }
        struct Evidence: Codable {
            let schema: String
            let configuration: VivoConstantPHConfiguration
            let checkpoint: VivoConstantPHCheckpoint
            let populations: [String: Double]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/constant-ph-discrete-md-evidence/v1",
            configuration: configuration, checkpoint: checkpoint, populations: populations)))
        return .init(schema: VivoConstantPHResult.schema,
            configurationFingerprint: try configuration.fingerprint(),
            finalCheckpoint: checkpoint, populations: populations,
            acceptanceFraction: Double(checkpoint.acceptedMoves) / Double(configuration.attemptCount),
            interpretation: Self.interpretation, evidenceFingerprint: evidenceID)
    }

    private static func validate(checkpoint: VivoConstantPHCheckpoint,
                                 configuration: VivoConstantPHConfiguration,
                                 manifold: VivoFingerprint,
                                 hamiltonians: [String: VivoFingerprint]) throws {
        try checkpoint.physicalState.validate()
        guard checkpoint.schema == VivoConstantPHCheckpoint.schema,
              checkpoint.configurationFingerprint == (try configuration.fingerprint()),
              checkpoint.physicalManifoldFingerprint == manifold,
              checkpoint.executableHamiltonianFingerprints == hamiltonians,
              checkpoint.completedAttempts >= 0,
              checkpoint.completedAttempts <= configuration.attemptCount,
              configuration.states.contains(where: { $0.identifier == checkpoint.currentStateIdentifier }),
              checkpoint.attempts.count == checkpoint.completedAttempts,
              checkpoint.acceptedMoves >= 0, checkpoint.acceptedMoves <= checkpoint.completedAttempts,
              Set(checkpoint.visitCounts.keys) == Set(configuration.states.map(\.identifier)),
              checkpoint.visitCounts.values.allSatisfy({ $0 >= 0 }),
              checkpoint.visitCounts.values.reduce(0, +) == checkpoint.completedAttempts + 1 else {
            throw VivoChemistryError.invalid("constant-pH checkpoint identity or accounting")
        }
    }
}