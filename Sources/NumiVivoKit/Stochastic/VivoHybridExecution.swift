import CryptoKit
import Foundation

/// Numerical authorities used by the executable hybrid Metal ABI (version 1).
/// A species belongs to exactly one authority, including propensity-only inputs.
public enum VivoHybridGPUAuthority: UInt32, Codable, Sendable {
    case constant = 0
    case exactSSA = 1
    case tauLeap = 2
    case deterministicRK2 = 3
}

public struct VivoHybridRuntimeConfiguration: Codable, Equatable, Sendable {
    public var laneCount: UInt32
    public var timeStep: Float
    public var minimumTimeStep: Float
    public var maximumTimeStep: Float
    public var seed: UInt64
    public var exactEventsPerDispatch: UInt32
    public var maximumExactDispatches: UInt32
    public var workingSetFraction: Double
    public var maximumPublications: UInt32

    public init(laneCount: UInt32, timeStep: Float = 0.01,
                minimumTimeStep: Float = 1e-7, maximumTimeStep: Float = 60,
                seed: UInt64 = 0x4e554d495649564f,
                exactEventsPerDispatch: UInt32 = 256,
                maximumExactDispatches: UInt32 = 64,
                workingSetFraction: Double = 0.75,
                maximumPublications: UInt32 = 65_536) {
        self.laneCount = laneCount
        self.timeStep = timeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.seed = seed
        self.exactEventsPerDispatch = exactEventsPerDispatch
        self.maximumExactDispatches = maximumExactDispatches
        self.workingSetFraction = workingSetFraction
        self.maximumPublications = maximumPublications
    }

    public func validate() throws {
        guard laneCount > 0, timeStep.isFinite, minimumTimeStep.isFinite,
              maximumTimeStep.isFinite, minimumTimeStep > 0,
              timeStep >= minimumTimeStep, timeStep <= maximumTimeStep,
              exactEventsPerDispatch > 0, maximumExactDispatches > 0,
              UInt64(exactEventsPerDispatch) * UInt64(maximumExactDispatches) <= UInt64(UInt32.max),
              workingSetFraction.isFinite, workingSetFraction > 0,
              workingSetFraction <= 0.95, maximumPublications > 0 else {
            throw VivoHybridExecutionError.invalidConfiguration
        }
    }
}

public enum VivoHybridExecutionError: Error, LocalizedError, Sendable {
    case invalidConfiguration
    case invalidPlan(String)
    case invalidState(String)
    case incompatibleCheckpoint(String)
    case transactionConflict
    case missingTransaction
    case resourceLimit(String)
    case metal(String)

    public var errorDescription: String? {
        switch self {
        case .invalidConfiguration: return "invalid hybrid runtime configuration"
        case .invalidPlan(let reason): return "invalid executable hybrid plan: \(reason)"
        case .invalidState(let reason): return "invalid hybrid state: \(reason)"
        case .incompatibleCheckpoint(let reason): return "incompatible hybrid checkpoint: \(reason)"
        case .transactionConflict: return "hybrid runtime has an in-flight operation or prepared transaction"
        case .missingTransaction: return "hybrid transaction identifier does not match the prepared candidate"
        case .resourceLimit(let reason): return "hybrid resource limit: \(reason)"
        case .metal(let reason): return "hybrid Metal failure: \(reason)"
        }
    }
}

public struct VivoHybridGPUStatus: Codable, Equatable, Sendable {
    public let flags: UInt32
    public let firstInvalidLane: UInt32?
    public let firstInvalidReaction: UInt32?
    public let unfinishedExactLanes: UInt32
    public let maximumExactEvents: UInt32

    public var valid: Bool { flags == 0 }
    public var complete: Bool { valid && unfinishedExactLanes == 0 }
    // flags: 1 nonfinite/invalid propensity; 2 count underflow; 4 count overflow;
    // 8 continuous-state bounds; 16 Poisson sampler budget; 32 stalled event time.
}

public struct VivoHybridPublicationRequest: Codable, Equatable, Sendable {
    public let speciesIndex: UInt32
    public let laneIndex: UInt32
    public init(speciesIndex: UInt32, laneIndex: UInt32) {
        self.speciesIndex = speciesIndex
        self.laneIndex = laneIndex
    }
}

public struct VivoHybridPublication: Codable, Equatable, Sendable {
    public let speciesIndex: UInt32
    public let laneIndex: UInt32
    public let authority: VivoHybridGPUAuthority
    public let value: Double
    /// Preserves all 32 count bits. Never reconstruct this from an FP32 value.
    public let exactCount: UInt32?
}

public enum VivoHybridStepDisposition: String, Codable, Sendable {
    case prepared
    case committed
    case rejected
    case exactWorkBudgetExceeded
}

public struct VivoHybridStepCertificate: Codable, Equatable, Sendable {
    public let transactionID: UUID
    public let disposition: VivoHybridStepDisposition
    public let modelFingerprint: String
    public let planFingerprint: String
    public let stepIndex: UInt64
    public let timeBefore: Double
    public let timeAfter: Double
    public let requestedTimeStep: Float
    public let exactDispatches: UInt32
    public let status: VivoHybridGPUStatus
    public let publications: [VivoHybridPublication]
    public var canCommit: Bool { disposition == .prepared && status.complete }
}

public struct VivoHybridStateSnapshot: Codable, Equatable, Sendable {
    public let modelFingerprint: String
    public let planFingerprint: String
    public let stepIndex: UInt64
    public let timeSeconds: Double
    public let laneCount: UInt32
    public let authorities: [VivoHybridGPUAuthority]
    public let continuousValues: [Float]
    public let counts: [UInt32]

    public func value(species: UInt32, lane: UInt32) -> Double? {
        guard Int(species) < authorities.count, lane < laneCount else { return nil }
        let index = UInt64(species) * UInt64(laneCount) + UInt64(lane)
        guard index < UInt64(counts.count), index < UInt64(continuousValues.count) else { return nil }
        return authorities[Int(species)] == .deterministicRK2
            ? Double(continuousValues[Int(index)]) : Double(counts[Int(index)])
    }
}

/// Complete accepted-boundary state. The random namespace is (seed, accepted
/// step, lane, reaction/cohort, draw). No uncommitted SSA continuation is exported.
public struct VivoHybridCheckpoint: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let numericalABIVersion: UInt32
    public let modelFingerprint: String
    public let planFingerprint: String
    public let seed: UInt64
    public let stepIndex: UInt64
    public let timeSeconds: Double
    public let laneCount: UInt32
    public let speciesCount: UInt32
    public let continuousFP32LE: Data
    public let countsUInt32LE: Data

    public init(modelFingerprint: String, planFingerprint: String, seed: UInt64,
                snapshot: VivoHybridStateSnapshot) {
        schemaVersion = 1
        numericalABIVersion = 1
        self.modelFingerprint = modelFingerprint
        self.planFingerprint = planFingerprint
        self.seed = seed
        stepIndex = snapshot.stepIndex
        timeSeconds = snapshot.timeSeconds
        laneCount = snapshot.laneCount
        speciesCount = UInt32(snapshot.authorities.count)
        continuousFP32LE = VivoLittleEndianFP32.encode(snapshot.continuousValues)
        countsUInt32LE = Self.encodeCounts(snapshot.counts)
    }

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    static func encodeCounts(_ counts: [UInt32]) -> Data {
        var result = Data(capacity: counts.count * 4)
        for count in counts {
            var value = count.littleEndian
            withUnsafeBytes(of: &value) { result.append(contentsOf: $0) }
        }
        return result
    }

    static func decodeCounts(_ data: Data) throws -> [UInt32] {
        guard data.count.isMultiple(of: 4) else {
            throw VivoHybridExecutionError.incompatibleCheckpoint("count payload is truncated")
        }
        return data.withUnsafeBytes { bytes in
            stride(from: 0, to: bytes.count, by: 4).map {
                UInt32(littleEndian: bytes.loadUnaligned(fromByteOffset: $0, as: UInt32.self))
            }
        }
    }
}

// These records are mirrored verbatim in NumiVivoHybridExecution.metal.
struct VivoHybridGPUReaction: Sendable {
    var law: UInt32
    var reactantA: UInt32
    var reactantB: UInt32
    var changeOffset: UInt32
    var changeCount: UInt32
    var reserved: UInt32 = 0
    var reserved1: UInt32 = 0
    var reserved2: UInt32 = 0
    var rate: Float
    var parameter1: Float
    var parameter2: Float
    var parameter3: Float
}
struct VivoHybridGPUChange: Sendable { var index: UInt32; var delta: Int32 }
struct VivoHybridGPUCohort: Sendable { var reactionOffset: UInt32; var reactionCount: UInt32 }
struct VivoHybridGPUCommand: Sendable {
    var laneCount: UInt32
    var speciesCount: UInt32
    var reactionCount: UInt32
    var exactCohortCount: UInt32
    var stepLow: UInt32
    var stepHigh: UInt32
    var stage: UInt32
    var eventsPerDispatch: UInt32
    var dt: Float
    var reserved: UInt32 = 0
    var seed: UInt64
}

public struct VivoCompiledHybridExecution: Sendable {
    public let modelFingerprint: String
    public let planFingerprint: String
    public let species: [String]
    public let authorities: [VivoHybridGPUAuthority]
    public let reactionAuthorities: [VivoHybridGPUAuthority]
    public let recommendedMaximumStep: Double
    let reactions: [VivoHybridGPUReaction]
    let changes: [VivoHybridGPUChange]
    let incidenceOffsets: [UInt32]
    let incidence: [VivoHybridGPUChange]
    let exactCohorts: [VivoHybridGPUCohort]
    let exactReactionIndices: [UInt32]
}

public enum VivoHybridExecutionCompiler {
    /// Lowers a validated source model and an authority plan to executable sparse
    /// tables. Plan membership is checked against actual propensity dependencies,
    /// not trusted from the caller's reaction descriptors.
    public static func compile(model: VivoExactSSAModel,
                               plan: VivoHybridStochasticPlan) throws -> VivoCompiledHybridExecution {
        try model.validate()
        guard plan.schemaVersion == 1, plan.speciesCount == UInt32(model.species.count) else {
            throw VivoHybridExecutionError.invalidPlan("version or species count mismatch")
        }
        let verifiedModel = try VivoExactSSAModel(species: model.species, reactions: model.reactions)
        guard verifiedModel.fingerprint == model.fingerprint else {
            throw VivoHybridExecutionError.invalidPlan("model fingerprint mismatch; compile source before execution")
        }
        var speciesOwners = [UInt32?](repeating: nil, count: model.species.count)
        var reactionOwners = [UInt32?](repeating: nil, count: model.reactions.count)
        var authorities = [VivoHybridGPUAuthority](repeating: .constant, count: model.species.count)
        var reactionAuthorities = [VivoHybridGPUAuthority](repeating: .constant, count: model.reactions.count)
        var seenCohorts = Set<UInt32>()
        var exactCohorts: [VivoHybridGPUCohort] = []
        var exactIndices: [UInt32] = []
        let ordered = plan.cohorts.sorted { $0.cohortID < $1.cohortID }
        for cohort in ordered {
            guard seenCohorts.insert(cohort.cohortID).inserted,
                  !cohort.speciesIndices.isEmpty, !cohort.reactionIndices.isEmpty,
                  cohort.recommendedMaximumStep.isFinite, cohort.recommendedMaximumStep > 0 else {
                throw VivoHybridExecutionError.invalidPlan("duplicate, empty, or invalid cohort")
            }
            let authority: VivoHybridGPUAuthority
            switch cohort.mode {
            case .exactSSA: authority = .exactSSA
            case .tauLeap: authority = .tauLeap
            case .deterministicRK2: authority = .deterministicRK2
            case .spatialSplit:
                throw VivoHybridExecutionError.invalidPlan("spatial cohorts require a transport participant; they cannot execute as well-mixed reactions")
            }
            for species in cohort.speciesIndices {
                guard Int(species) < speciesOwners.count, speciesOwners[Int(species)] == nil else {
                    throw VivoHybridExecutionError.invalidPlan("overlapping or out-of-range species ownership")
                }
                speciesOwners[Int(species)] = cohort.cohortID
                authorities[Int(species)] = authority
            }
            for reaction in cohort.reactionIndices {
                guard Int(reaction) < reactionOwners.count, reactionOwners[Int(reaction)] == nil else {
                    throw VivoHybridExecutionError.invalidPlan("overlapping or out-of-range reaction ownership")
                }
                reactionOwners[Int(reaction)] = cohort.cohortID
                reactionAuthorities[Int(reaction)] = authority
            }
            if authority == .exactSSA {
                exactCohorts.append(.init(reactionOffset: UInt32(exactIndices.count),
                                          reactionCount: UInt32(cohort.reactionIndices.count)))
                exactIndices.append(contentsOf: cohort.reactionIndices.sorted())
            }
        }
        guard reactionOwners.allSatisfy({ $0 != nil }) else {
            throw VivoHybridExecutionError.invalidPlan("every reaction must have one authority")
        }
        let uncoupled = Set(plan.uncoupledSpecies)
        guard uncoupled.count == plan.uncoupledSpecies.count,
              uncoupled.allSatisfy({ Int($0) < speciesOwners.count }),
              Set(speciesOwners.indices.filter { speciesOwners[$0] == nil }.map(UInt32.init)) == uncoupled else {
            throw VivoHybridExecutionError.invalidPlan("uncoupled species do not equal the unowned species set")
        }
        var records: [VivoHybridGPUReaction] = []
        var changes: [VivoHybridGPUChange] = []
        var incidence = [[VivoHybridGPUChange]](repeating: [], count: model.species.count)
        for (index, reaction) in model.reactions.enumerated() {
            let dependencies = Set(reaction.changes.map(\.speciesIndex) +
                                   [reaction.reactantA, reaction.reactantB].compactMap { $0 })
            guard dependencies.allSatisfy({ speciesOwners[Int($0)] == reactionOwners[index] }) else {
                throw VivoHybridExecutionError.invalidPlan("reaction \(reaction.id) reads or writes another cohort's species")
            }
            guard changes.count <= Int(UInt32.max) - reaction.changes.count else {
                throw VivoHybridExecutionError.resourceLimit("stoichiometry exceeds UInt32")
            }
            records.append(.init(law: reaction.law.rawValue,
                                 reactantA: reaction.reactantA ?? .max,
                                 reactantB: reaction.reactantB ?? .max,
                                 changeOffset: UInt32(changes.count),
                                 changeCount: UInt32(reaction.changes.count),
                                 rate: reaction.rate, parameter1: reaction.parameter1,
                                 parameter2: reaction.parameter2, parameter3: reaction.parameter3))
            for change in reaction.changes.sorted(by: { $0.speciesIndex < $1.speciesIndex }) {
                changes.append(.init(index: change.speciesIndex, delta: change.delta))
                incidence[Int(change.speciesIndex)].append(.init(index: UInt32(index), delta: change.delta))
            }
        }
        var offsets: [UInt32] = [0]
        var flatIncidence: [VivoHybridGPUChange] = []
        for row in incidence {
            flatIncidence.append(contentsOf: row)
            offsets.append(UInt32(flatIncidence.count))
        }
        struct Identity: Encodable {
            let abi: UInt32
            let model: String
            let speciesAuthorities: [VivoHybridGPUAuthority]
            let reactionAuthorities: [VivoHybridGPUAuthority]
            let exactReactionGroups: [[UInt32]]
        }
        let identity = Identity(abi: 1, model: verifiedModel.fingerprint,
                                speciesAuthorities: authorities, reactionAuthorities: reactionAuthorities,
                                exactReactionGroups: ordered.filter { $0.mode == .exactSSA }.map { $0.reactionIndices.sorted() })
        let digest = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity)).hex
        return .init(modelFingerprint: verifiedModel.fingerprint, planFingerprint: digest,
                     species: model.species, authorities: authorities,
                     reactionAuthorities: reactionAuthorities,
                     recommendedMaximumStep: ordered.map(\.recommendedMaximumStep).min() ?? 60,
                     reactions: records, changes: changes, incidenceOffsets: offsets,
                     incidence: flatIncidence, exactCohorts: exactCohorts,
                     exactReactionIndices: exactIndices)
    }
}
