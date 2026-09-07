import Foundation

public struct VivoConstantPHReplicaExchangeConfiguration: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-replica-exchange/v1"
    public var schema: String
    public var identifier: String
    public var temperatureK: Double
    public var referencePH: Double
    public var pHLadder: [Double]
    public var states: [VivoConstantPHStateDefinition]
    public var initialStateIdentifier: String
    public var mdStepsPerCycle: UInt64
    public var cycleCount: Int
    public var seed: UInt64

    public init(identifier: String, temperatureK: Double, referencePH: Double,
                pHLadder: [Double], states: [VivoConstantPHStateDefinition],
                initialStateIdentifier: String, mdStepsPerCycle: UInt64,
                cycleCount: Int, seed: UInt64) {
        schema = Self.schema; self.identifier = identifier; self.temperatureK = temperatureK
        self.referencePH = referencePH; self.pHLadder = pHLadder; self.states = states
        self.initialStateIdentifier = initialStateIdentifier; self.mdStepsPerCycle = mdStepsPerCycle
        self.cycleCount = cycleCount; self.seed = seed
    }

    public func validate() throws {
        guard schema == Self.schema,
              !identifier.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              identifier.utf8.count <= 512,
              temperatureK.isFinite, temperatureK > 0,
              referencePH.isFinite, (-10...30).contains(referencePH),
              pHLadder.count >= 2, pHLadder.count <= 128,
              pHLadder.allSatisfy({ $0.isFinite && (-10...30).contains($0) }),
              zip(pHLadder, pHLadder.dropFirst()).allSatisfy({ $0.0 < $0.1 }),
              !states.isEmpty, states.count <= 4096,
              Set(states.map(\.identifier)).count == states.count,
              states.contains(where: { $0.identifier == initialStateIdentifier }),
              mdStepsPerCycle > 0, mdStepsPerCycle <= 10_000_000,
              (1...1_000_000).contains(cycleCount) else {
            throw VivoChemistryError.invalid("constant-pH replica-exchange configuration")
        }
        // Reuse the discrete sampler's graph validation without creating a runtime.
        let probe = VivoConstantPHConfiguration(identifier: identifier + ".graph-validation",
            temperatureK: temperatureK, referencePH: referencePH, targetPH: pHLadder[0],
            states: states, initialStateIdentifier: initialStateIdentifier,
            mdStepsPerAttempt: mdStepsPerCycle, attemptCount: 1, seed: seed)
        try probe.validate()
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoConstantPHReplicaLane: Codable, Sendable, Equatable {
    public let laneIndex: Int
    public let targetPH: Double
    public let walkerIdentifier: Int
    public let chemicalStateIdentifier: String
    public let physicalState: VivoConstantPHPhysicalState
}

public struct VivoConstantPHReplicaSwapAttempt: Codable, Sendable, Equatable {
    public let cycleIndex: Int
    public let lowerLaneIndex: Int
    public let upperLaneIndex: Int
    public let lowerPH: Double
    public let upperPH: Double
    public let lowerStateIdentifier: String
    public let upperStateIdentifier: String
    public let logAcceptanceProbability: Double
    public let logUniform: Double
    public let accepted: Bool
}

public struct VivoConstantPHReplicaPopulation: Codable, Sendable, Equatable {
    public let targetPH: Double
    public let stateIdentifier: String
    public let visits: Int
    public let population: Double
}

public struct VivoConstantPHReplicaExchangeCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-replica-exchange-checkpoint/v1"
    public var schema: String
    public let configurationFingerprint: VivoFingerprint
    public let physicalManifoldFingerprint: VivoFingerprint
    public let executableHamiltonianFingerprints: [String: VivoFingerprint]
    public let completedCycles: Int
    public let lanes: [VivoConstantPHReplicaLane]
    public let randomState: UInt64
    public let chemicalMoveAttempts: Int
    public let chemicalMoveAccepts: Int
    public let swapAttempts: [VivoConstantPHReplicaSwapAttempt]
    /// key format is stable and internal: laneIndex + U+001F + state identifier.
    public let populationCounts: [String: Int]

    public init(configurationFingerprint: VivoFingerprint,
                physicalManifoldFingerprint: VivoFingerprint,
                executableHamiltonianFingerprints: [String: VivoFingerprint],
                completedCycles: Int, lanes: [VivoConstantPHReplicaLane], randomState: UInt64,
                chemicalMoveAttempts: Int, chemicalMoveAccepts: Int,
                swapAttempts: [VivoConstantPHReplicaSwapAttempt], populationCounts: [String: Int]) {
        schema = Self.schema; self.configurationFingerprint = configurationFingerprint
        self.physicalManifoldFingerprint = physicalManifoldFingerprint
        self.executableHamiltonianFingerprints = executableHamiltonianFingerprints
        self.completedCycles = completedCycles; self.lanes = lanes; self.randomState = randomState
        self.chemicalMoveAttempts = chemicalMoveAttempts; self.chemicalMoveAccepts = chemicalMoveAccepts
        self.swapAttempts = swapAttempts; self.populationCounts = populationCounts
    }
}

public struct VivoConstantPHReplicaExchangeResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/constant-ph-replica-exchange-result/v1"
    public let schema: String
    public let configurationFingerprint: VivoFingerprint
    public let finalCheckpoint: VivoConstantPHReplicaExchangeCheckpoint
    public let populations: [VivoConstantPHReplicaPopulation]
    public let chemicalMoveAcceptanceFraction: Double
    public let swapAcceptanceFraction: Double
    public let interpretation: String
    public let evidenceFingerprint: VivoFingerprint
}

public typealias VivoConstantPHReplicaExchangeProgressSink = @Sendable (VivoConstantPHReplicaExchangeCheckpoint) async throws -> Void

/// Fixed-temperature pH replica exchange. Physical configurations and chemical
/// state identities swap between fixed pH lanes. The Hamiltonian energy cancels
/// in the pH-label exchange criterion, leaving only the proton-reservoir term.
public actor VivoConstantPHReplicaExchange {
    private static let gasConstantKJ = 0.00831446261815324
    private static let separator = "\u{1F}"
    public static let interpretation = "Fixed-temperature pH replica exchange combined with discrete-state semigrand MD/MC. Each pH lane performs a complete-Hamiltonian chemical-state move, then adjacent lanes exchange full physical/chemical configurations using the exact proton-reservoir Metropolis criterion. State rules, reference corrections and pH ladder remain explicit inputs."

    public let configuration: VivoConstantPHReplicaExchangeConfiguration
    private let runtime: [String: VivoConstantPHExecutableState]
    private let manifold: VivoFingerprint
    private let hamiltonians: [String: VivoFingerprint]
    private var committed: VivoConstantPHReplicaExchangeCheckpoint
    private var inFlight = false

    public init(configuration: VivoConstantPHReplicaExchangeConfiguration,
                initialPhysicalStates: [VivoConstantPHPhysicalState],
                executableStates: [VivoConstantPHExecutableState],
                checkpoint: VivoConstantPHReplicaExchangeCheckpoint? = nil) throws {
        try configuration.validate()
        guard initialPhysicalStates.count == configuration.pHLadder.count else {
            throw VivoChemistryError.invalid("constant-pH replica initial-state count")
        }
        for state in initialPhysicalStates { try state.validate() }
        guard executableStates.count == configuration.states.count,
              Set(executableStates.map(\.identifier)) == Set(configuration.states.map(\.identifier)),
              Set(executableStates.map(\.physicalManifoldFingerprint)).count == 1,
              let manifold = executableStates.first?.physicalManifoldFingerprint else {
            throw VivoChemistryError.invalid("constant-pH replica executable state set")
        }
        let runtime = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0) })
        let hamiltonians = Dictionary(uniqueKeysWithValues: executableStates.map { ($0.identifier, $0.hamiltonianFingerprint) })
        let configID = try configuration.fingerprint()
        if let checkpoint {
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians)
            self.committed = checkpoint
        } else {
            let lanes = configuration.pHLadder.indices.map { index in
                VivoConstantPHReplicaLane(laneIndex: index, targetPH: configuration.pHLadder[index],
                    walkerIdentifier: index, chemicalStateIdentifier: configuration.initialStateIdentifier,
                    physicalState: initialPhysicalStates[index])
            }
            var counts: [String: Int] = [:]
            for lane in lanes { counts[Self.populationKey(lane.laneIndex, lane.chemicalStateIdentifier)] = 1 }
            self.committed = .init(configurationFingerprint: configID,
                physicalManifoldFingerprint: manifold, executableHamiltonianFingerprints: hamiltonians,
                completedCycles: 0, lanes: lanes, randomState: configuration.seed,
                chemicalMoveAttempts: 0, chemicalMoveAccepts: 0, swapAttempts: [], populationCounts: counts)
        }
        self.configuration = configuration; self.runtime = runtime
        self.manifold = manifold; self.hamiltonians = hamiltonians
    }

    public func checkpoint() -> VivoConstantPHReplicaExchangeCheckpoint { committed }

    public func run(progress: VivoConstantPHReplicaExchangeProgressSink? = nil) async throws -> VivoConstantPHReplicaExchangeResult {
        guard !inFlight else { throw VivoChemistryError.invalid("constant-pH replica exchange already in flight") }
        inFlight = true; defer { inFlight = false }
        let definitions = Dictionary(uniqueKeysWithValues: configuration.states.map { ($0.identifier, $0) })
        let rt = Self.gasConstantKJ * configuration.temperatureK
        var checkpoint = committed
        var random = VivoSplitMix64(state: checkpoint.randomState)
        while checkpoint.completedCycles < configuration.cycleCount {
            try Task.checkCancellation()
            var lanes = checkpoint.lanes
            var moveAccepts = checkpoint.chemicalMoveAccepts
            var moveAttempts = checkpoint.chemicalMoveAttempts
            // One fixed-pH MD/MC move per lane. Lane order is deterministic.
            for index in lanes.indices {
                let lane = lanes[index]
                guard let currentDefinition = definitions[lane.chemicalStateIdentifier],
                      let currentRuntime = runtime[lane.chemicalStateIdentifier] else {
                    throw VivoChemistryError.invalid("constant-pH replica current state disappeared")
                }
                let propagated = try await currentRuntime.propagate(lane.physicalState, configuration.mdStepsPerCycle)
                try propagated.validate()
                guard propagated.positionsNM.count == lane.physicalState.positionsNM.count,
                      propagated.stepIndex >= lane.physicalState.stepIndex,
                      propagated.timePS >= lane.physicalState.timePS else {
                    throw VivoChemistryError.invalid("constant-pH replica propagator changed manifold or reversed clock")
                }
                let selector = min(Int(random.unitInterval() * Double(currentDefinition.neighbors.count)), currentDefinition.neighbors.count - 1)
                let proposedIdentifier = currentDefinition.neighbors[selector]
                guard let proposedDefinition = definitions[proposedIdentifier], let proposedRuntime = runtime[proposedIdentifier] else {
                    throw VivoChemistryError.invalid("constant-pH replica proposal state disappeared")
                }
                async let fromValue = currentRuntime.potentialEnergyKJPerMol(propagated)
                async let proposedValue = proposedRuntime.potentialEnergyKJPerMol(propagated)
                let fromEnergy = try await fromValue, proposedEnergy = try await proposedValue
                guard fromEnergy.isFinite, proposedEnergy.isFinite else { throw VivoChemistryError.convergence("nonfinite constant-pH replica energy") }
                let slope = rt * log(10.0) * (lane.targetPH - configuration.referencePH)
                let delta = proposedEnergy + proposedDefinition.referenceSemigrandBiasKJPerMol + Double(proposedDefinition.boundProtonOffset) * slope
                    - fromEnergy - currentDefinition.referenceSemigrandBiasKJPerMol - Double(currentDefinition.boundProtonOffset) * slope
                let logProposalRatio = log(Double(currentDefinition.neighbors.count) / Double(proposedDefinition.neighbors.count))
                let logAcceptance = min(0.0, -delta / rt + logProposalRatio)
                let accepted = log(max(random.unitInterval(), Double.leastNonzeroMagnitude)) < logAcceptance
                moveAttempts += 1; if accepted { moveAccepts += 1 }
                lanes[index] = .init(laneIndex: lane.laneIndex, targetPH: lane.targetPH,
                    walkerIdentifier: lane.walkerIdentifier,
                    chemicalStateIdentifier: accepted ? proposedIdentifier : lane.chemicalStateIdentifier,
                    physicalState: propagated)
            }
            // Alternating even/odd neighbor pairing gives every adjacent edge a chance.
            var swaps = checkpoint.swapAttempts
            let parity = checkpoint.completedCycles & 1
            var lower = parity
            while lower + 1 < lanes.count {
                let upper = lower + 1
                let a = lanes[lower], b = lanes[upper]
                guard let da = definitions[a.chemicalStateIdentifier], let db = definitions[b.chemicalStateIdentifier] else {
                    throw VivoChemistryError.invalid("constant-pH replica swap state disappeared")
                }
                let logRatio = log(10.0) * Double(da.boundProtonOffset - db.boundProtonOffset) * (a.targetPH - b.targetPH)
                let logAcceptance = min(0.0, logRatio)
                let logUniform = log(max(random.unitInterval(), Double.leastNonzeroMagnitude))
                let accepted = logUniform < logAcceptance
                swaps.append(.init(cycleIndex: checkpoint.completedCycles,
                    lowerLaneIndex: lower, upperLaneIndex: upper,
                    lowerPH: a.targetPH, upperPH: b.targetPH,
                    lowerStateIdentifier: a.chemicalStateIdentifier, upperStateIdentifier: b.chemicalStateIdentifier,
                    logAcceptanceProbability: logAcceptance, logUniform: logUniform, accepted: accepted))
                if accepted {
                    lanes[lower] = .init(laneIndex: lower, targetPH: a.targetPH,
                        walkerIdentifier: b.walkerIdentifier, chemicalStateIdentifier: b.chemicalStateIdentifier,
                        physicalState: b.physicalState)
                    lanes[upper] = .init(laneIndex: upper, targetPH: b.targetPH,
                        walkerIdentifier: a.walkerIdentifier, chemicalStateIdentifier: a.chemicalStateIdentifier,
                        physicalState: a.physicalState)
                }
                lower += 2
            }
            var counts = checkpoint.populationCounts
            for lane in lanes { counts[Self.populationKey(lane.laneIndex, lane.chemicalStateIdentifier), default: 0] += 1 }
            checkpoint = .init(configurationFingerprint: checkpoint.configurationFingerprint,
                physicalManifoldFingerprint: manifold, executableHamiltonianFingerprints: hamiltonians,
                completedCycles: checkpoint.completedCycles + 1, lanes: lanes, randomState: random.state,
                chemicalMoveAttempts: moveAttempts, chemicalMoveAccepts: moveAccepts,
                swapAttempts: swaps, populationCounts: counts)
            try Self.validate(checkpoint: checkpoint, configuration: configuration,
                              manifold: manifold, hamiltonians: hamiltonians)
            committed = checkpoint
            if let progress { try await progress(checkpoint) }
        }
        let visitsPerLane = configuration.cycleCount + 1
        var populations: [VivoConstantPHReplicaPopulation] = []
        for laneIndex in configuration.pHLadder.indices {
            for state in configuration.states {
                let visits = checkpoint.populationCounts[Self.populationKey(laneIndex, state.identifier), default: 0]
                populations.append(.init(targetPH: configuration.pHLadder[laneIndex], stateIdentifier: state.identifier,
                    visits: visits, population: Double(visits) / Double(visitsPerLane)))
            }
        }
        let swapAccepted = checkpoint.swapAttempts.filter(\.accepted).count
        let swapFraction = checkpoint.swapAttempts.isEmpty ? 0 : Double(swapAccepted) / Double(checkpoint.swapAttempts.count)
        struct Evidence: Codable {
            let schema: String; let configuration: VivoConstantPHReplicaExchangeConfiguration
            let checkpoint: VivoConstantPHReplicaExchangeCheckpoint; let populations: [VivoConstantPHReplicaPopulation]
        }
        let evidenceID = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(Evidence(
            schema: "numivivo.org/constant-ph-replica-exchange-evidence/v1",
            configuration: configuration, checkpoint: checkpoint, populations: populations)))
        return .init(schema: VivoConstantPHReplicaExchangeResult.schema,
            configurationFingerprint: try configuration.fingerprint(), finalCheckpoint: checkpoint,
            populations: populations,
            chemicalMoveAcceptanceFraction: Double(checkpoint.chemicalMoveAccepts) / Double(checkpoint.chemicalMoveAttempts),
            swapAcceptanceFraction: swapFraction,
            interpretation: Self.interpretation, evidenceFingerprint: evidenceID)
    }

    private static func populationKey(_ lane: Int, _ state: String) -> String { "\(lane)\(separator)\(state)" }

    private static func validate(checkpoint: VivoConstantPHReplicaExchangeCheckpoint,
                                 configuration: VivoConstantPHReplicaExchangeConfiguration,
                                 manifold: VivoFingerprint,
                                 hamiltonians: [String: VivoFingerprint]) throws {
        guard checkpoint.schema == VivoConstantPHReplicaExchangeCheckpoint.schema,
              checkpoint.configurationFingerprint == (try configuration.fingerprint()),
              checkpoint.physicalManifoldFingerprint == manifold,
              checkpoint.executableHamiltonianFingerprints == hamiltonians,
              checkpoint.completedCycles >= 0, checkpoint.completedCycles <= configuration.cycleCount,
              checkpoint.lanes.count == configuration.pHLadder.count,
              checkpoint.chemicalMoveAttempts == checkpoint.completedCycles * configuration.pHLadder.count,
              checkpoint.chemicalMoveAccepts >= 0, checkpoint.chemicalMoveAccepts <= checkpoint.chemicalMoveAttempts,
              Set(checkpoint.lanes.map(\.walkerIdentifier)) == Set(configuration.pHLadder.indices),
              checkpoint.populationCounts.values.allSatisfy({ $0 >= 0 }) else {
            throw VivoChemistryError.invalid("constant-pH replica checkpoint identity or accounting")
        }
        let stateIDs = Set(configuration.states.map(\.identifier))
        for lane in checkpoint.lanes {
            try lane.physicalState.validate()
            guard lane.laneIndex >= 0, lane.laneIndex < configuration.pHLadder.count,
                  lane.targetPH == configuration.pHLadder[lane.laneIndex],
                  stateIDs.contains(lane.chemicalStateIdentifier) else {
                throw VivoChemistryError.invalid("constant-pH replica lane identity")
            }
        }
        let expectedVisits = (checkpoint.completedCycles + 1) * configuration.pHLadder.count
        guard checkpoint.populationCounts.values.reduce(0, +) == expectedVisits else {
            throw VivoChemistryError.invalid("constant-pH replica population accounting")
        }
    }
}