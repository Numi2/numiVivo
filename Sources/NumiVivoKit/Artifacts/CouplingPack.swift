import Foundation

public enum VivoSubsystem: String, Codable, Sendable, CaseIterable {
    case numiVivo
    case numiTissue
    case numanX
    case numiBrain
    case experiment
    case external
}

public enum VivoCouplingInterpolation: String, Codable, Sendable {
    case hold
    case linear
    case conservativeAverage
    case nearest
    case barycentric
}

public enum VivoCouplingAggregation: String, Codable, Sendable {
    case direct
    case mean
    case sum
    case minimum
    case maximum
    case volumeWeightedMean
}

public enum VivoCouplingAuthority: String, Codable, Sendable {
    case sourceAuthoritative
    case destinationAuthoritative
    case negotiated
}

public enum VivoLaneSelection: Codable, Sendable, Equatable {
    case all
    case lane(UInt32)
    case environment(UInt32)
    case voxel(UInt32)
    case explicit([UInt32])

    private enum CodingKeys: String, CodingKey {
        case kind
        case index
        case indices
    }

    private enum Kind: String, Codable {
        case all
        case lane
        case environment
        case voxel
        case explicit
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .all:
            self = .all
        case .lane:
            self = .lane(try container.decode(UInt32.self, forKey: .index))
        case .environment:
            self = .environment(try container.decode(UInt32.self, forKey: .index))
        case .voxel:
            self = .voxel(try container.decode(UInt32.self, forKey: .index))
        case .explicit:
            self = .explicit(try container.decode([UInt32].self, forKey: .indices))
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .all:
            try container.encode(Kind.all, forKey: .kind)
        case .lane(let index):
            try container.encode(Kind.lane, forKey: .kind)
            try container.encode(index, forKey: .index)
        case .environment(let index):
            try container.encode(Kind.environment, forKey: .kind)
            try container.encode(index, forKey: .index)
        case .voxel(let index):
            try container.encode(Kind.voxel, forKey: .kind)
            try container.encode(index, forKey: .index)
        case .explicit(let indices):
            try container.encode(Kind.explicit, forKey: .kind)
            try container.encode(indices, forKey: .indices)
        }
    }

    public func resolve(configuration: VivoRuntimeConfiguration) throws -> [UInt32] {
        let laneCount = configuration.laneCount
        let voxelCount = configuration.voxelCount
        switch self {
        case .all:
            return Array(0..<laneCount)
        case .lane(let lane):
            guard lane < laneCount else {
                throw VivoArtifactValidationError.invalid("coupling lane \(lane) exceeds laneCount \(laneCount)")
            }
            return [lane]
        case .environment(let environment):
            guard environment < configuration.environmentCount else {
                throw VivoArtifactValidationError.invalid("coupling environment \(environment) is out of bounds")
            }
            let base = environment.multipliedReportingOverflow(by: voxelCount)
            guard !base.overflow else {
                throw VivoArtifactValidationError.invalid("coupling environment lane offset overflow")
            }
            return (0..<voxelCount).map { base.partialValue + $0 }
        case .voxel(let voxel):
            guard voxel < voxelCount else {
                throw VivoArtifactValidationError.invalid("coupling voxel \(voxel) is out of bounds")
            }
            return (0..<configuration.environmentCount).map { environment in
                environment * voxelCount + voxel
            }
        case .explicit(let lanes):
            guard !lanes.isEmpty else {
                throw VivoArtifactValidationError.invalid("explicit coupling lane selection cannot be empty")
            }
            guard Set(lanes).count == lanes.count else {
                throw VivoArtifactValidationError.invalid("explicit coupling lane selection contains duplicates")
            }
            guard lanes.allSatisfy({ $0 < laneCount }) else {
                throw VivoArtifactValidationError.invalid("explicit coupling lane selection contains an out-of-bounds lane")
            }
            return lanes.sorted()
        }
    }
}

public struct VivoCouplingEndpoint: Codable, Sendable, Equatable {
    public var subsystem: VivoSubsystem
    public var channel: String
    public var unit: String
    public var lanes: VivoLaneSelection

    public init(
        subsystem: VivoSubsystem,
        channel: String,
        unit: String,
        lanes: VivoLaneSelection = .all
    ) {
        self.subsystem = subsystem
        self.channel = channel
        self.unit = unit
        self.lanes = lanes
    }
}

public struct VivoCouplingChannel: Codable, Sendable, Equatable {
    public var identifier: String
    public var source: VivoCouplingEndpoint
    public var destination: VivoCouplingEndpoint
    public var updateMode: VivoCouplingMode
    public var interpolation: VivoCouplingInterpolation
    public var aggregation: VivoCouplingAggregation
    public var authority: VivoCouplingAuthority
    public var cadence: VivoQuantity
    public var latency: VivoQuantity
    public var destinationBounds: VivoInterval?
    public var required: Bool
    public var transactional: Bool
    public var evidence: [VivoEvidenceReference]

    public init(
        identifier: String,
        source: VivoCouplingEndpoint,
        destination: VivoCouplingEndpoint,
        updateMode: VivoCouplingMode = .replace,
        interpolation: VivoCouplingInterpolation = .hold,
        aggregation: VivoCouplingAggregation = .direct,
        authority: VivoCouplingAuthority = .sourceAuthoritative,
        cadence: VivoQuantity,
        latency: VivoQuantity = .init(value: 0, unit: "s"),
        destinationBounds: VivoInterval? = nil,
        required: Bool = true,
        transactional: Bool = true,
        evidence: [VivoEvidenceReference] = []
    ) {
        self.identifier = identifier
        self.source = source
        self.destination = destination
        self.updateMode = updateMode
        self.interpolation = interpolation
        self.aggregation = aggregation
        self.authority = authority
        self.cadence = cadence
        self.latency = latency
        self.destinationBounds = destinationBounds
        self.required = required
        self.transactional = transactional
        self.evidence = evidence
    }
}

public struct VivoCouplingPolicy: Codable, Sendable, Equatable {
    public enum FailureAction: String, Codable, Sendable {
        case rejectGlobalStep
        case substep
        case holdLastValue
        case isolateSubsystem
        case reversibleShutdown
        case permanentShutdown
    }

    public var failureAction: FailureAction
    public var maximumStaleness: VivoQuantity
    public var requireAllParticipants: Bool
    public var requireMatchingStepIndex: Bool
    public var requireMatchingProgramFingerprint: Bool
    public var maximumSubsteps: UInt32

    public init(
        failureAction: FailureAction = .rejectGlobalStep,
        maximumStaleness: VivoQuantity,
        requireAllParticipants: Bool = true,
        requireMatchingStepIndex: Bool = true,
        requireMatchingProgramFingerprint: Bool = true,
        maximumSubsteps: UInt32 = 16
    ) {
        self.failureAction = failureAction
        self.maximumStaleness = maximumStaleness
        self.requireAllParticipants = requireAllParticipants
        self.requireMatchingStepIndex = requireMatchingStepIndex
        self.requireMatchingProgramFingerprint = requireMatchingProgramFingerprint
        self.maximumSubsteps = maximumSubsteps
    }
}

public struct VivoCouplingPack: Codable, Sendable, Equatable {
    public var programFingerprint: VivoFingerprint
    public var hostContextFingerprint: VivoFingerprint?
    public var channels: [VivoCouplingChannel]
    public var policy: VivoCouplingPolicy

    public init(
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint? = nil,
        channels: [VivoCouplingChannel],
        policy: VivoCouplingPolicy
    ) {
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.channels = channels
        self.policy = policy
    }

    public func validate(
        programPack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        units: VivoUnitSystem = .standard
    ) throws {
        guard programFingerprint == programPack.header.contentFingerprint else {
            throw VivoArtifactValidationError.incompatible("coupling pack references a different ProgramPack")
        }
        guard !channels.isEmpty else {
            throw VivoArtifactValidationError.invalid("coupling pack must contain at least one channel")
        }
        var identifiers = Set<String>()
        for channel in channels {
            guard !channel.identifier.isEmpty, identifiers.insert(channel.identifier).inserted else {
                throw VivoArtifactValidationError.invalid("coupling channel identifiers must be non-empty and unique")
            }
            guard channel.source.subsystem != channel.destination.subsystem ||
                  channel.source.channel != channel.destination.channel else {
                throw VivoArtifactValidationError.invalid("coupling channel \(channel.identifier) connects an endpoint to itself")
            }
            guard !channel.source.channel.isEmpty, !channel.destination.channel.isEmpty else {
                throw VivoArtifactValidationError.invalid("coupling endpoints require channel identifiers")
            }
            guard units.definition(for: channel.source.unit) != nil,
                  units.definition(for: channel.destination.unit) != nil,
                  units.areCompatible(channel.source.unit, channel.destination.unit) else {
                throw VivoArtifactValidationError.incompatible(
                    "coupling channel \(channel.identifier) has incompatible units \(channel.source.unit) and \(channel.destination.unit)"
                )
            }
            try channel.cadence.validate(label: "channels.\(channel.identifier).cadence", nonnegative: true)
            try channel.latency.validate(label: "channels.\(channel.identifier).latency", nonnegative: true)
            let cadenceSeconds = try units.convert(channel.cadence, to: "s")
            guard cadenceSeconds > 0 else {
                throw VivoArtifactValidationError.invalid("coupling cadence must be positive")
            }
            _ = try channel.source.lanes.resolve(configuration: configuration)
            _ = try channel.destination.lanes.resolve(configuration: configuration)
            try channel.destinationBounds?.validate(label: "channels.\(channel.identifier).destinationBounds")
        }
        try policy.maximumStaleness.validate(label: "coupling.policy.maximumStaleness", nonnegative: true)
        guard policy.maximumSubsteps > 0 else {
            throw VivoArtifactValidationError.invalid("coupling policy maximumSubsteps must be positive")
        }
    }
}

public struct VivoUnitTransform: Codable, Sendable, Equatable {
    public let scale: Double
    public let offset: Double

    public init(scale: Double, offset: Double) {
        self.scale = scale
        self.offset = offset
    }

    public func apply(_ value: Double) throws -> Double {
        let converted = value * scale + offset
        guard converted.isFinite else {
            throw VivoArtifactValidationError.invalid("coupling unit transform produced a non-finite value")
        }
        return converted
    }
}

public struct PreparedVivoCouplingChannel: Codable, Sendable, Equatable {
    public enum RuntimeDirection: String, Codable, Sendable {
        case inbound
        case outbound
    }

    public let identifier: String
    public let direction: RuntimeDirection
    public let peerSubsystem: VivoSubsystem
    public let peerChannel: String
    public let speciesIndex: UInt32
    public let lanes: [UInt32]
    public let updateMode: VivoCouplingMode
    public let aggregation: VivoCouplingAggregation
    public let interpolation: VivoCouplingInterpolation
    public let authority: VivoCouplingAuthority
    public let cadenceSeconds: Double
    public let latencySeconds: Double
    public let transform: VivoUnitTransform
    public let destinationBounds: VivoInterval?
    public let required: Bool
    public let transactional: Bool
}

public struct PreparedVivoCouplingPlan: Codable, Sendable, Equatable {
    public let fingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint?
    public let inbound: [PreparedVivoCouplingChannel]
    public let outbound: [PreparedVivoCouplingChannel]
    public let policy: VivoCouplingPolicy
}

public struct VivoCouplingCompiler: Sendable {
    private let units: VivoUnitSystem

    public init(units: VivoUnitSystem = .standard) {
        self.units = units
    }

    public func compile(
        _ coupling: VivoCouplingPack,
        for pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration
    ) throws -> PreparedVivoCouplingPlan {
        try coupling.validate(programPack: pack, configuration: configuration, units: units)
        let species = try pack.speciesMetadata()
        var indices: [String: Int] = [:]
        for (index, value) in species.enumerated() {
            guard indices.updateValue(index, forKey: value.identifier) == nil else {
                throw VivoArtifactValidationError.invalid("duplicate ProgramPack species \(value.identifier)")
            }
        }

        var inbound: [PreparedVivoCouplingChannel] = []
        var outbound: [PreparedVivoCouplingChannel] = []
        for channel in coupling.channels.sorted(by: { $0.identifier < $1.identifier }) {
            if channel.destination.subsystem == .numiVivo {
                guard channel.source.subsystem != .numiVivo else {
                    throw VivoArtifactValidationError.invalid(
                        "coupling channel \(channel.identifier) cannot be internal to NumiVivo"
                    )
                }
                guard let index = indices[channel.destination.channel] else {
                    throw VivoArtifactValidationError.unresolved(
                        "inbound coupling destination \(channel.destination.channel) is absent from ProgramPack species"
                    )
                }
                guard species[index].isExternallyOwned || species[index].isInput else {
                    throw VivoArtifactValidationError.incompatible(
                        "inbound coupling destination \(channel.destination.channel) is internally owned"
                    )
                }
                let transform = try unitTransform(from: channel.source.unit, to: species[index].unit)
                let lanes = try channel.destination.lanes.resolve(configuration: configuration)
                let bounds = try convertBounds(
                    channel.destinationBounds,
                    from: channel.destination.unit,
                    to: species[index].unit
                )
                inbound.append(.init(
                    identifier: channel.identifier,
                    direction: .inbound,
                    peerSubsystem: channel.source.subsystem,
                    peerChannel: channel.source.channel,
                    speciesIndex: UInt32(index),
                    lanes: lanes,
                    updateMode: channel.updateMode,
                    aggregation: channel.aggregation,
                    interpolation: channel.interpolation,
                    authority: channel.authority,
                    cadenceSeconds: try units.convert(channel.cadence, to: "s"),
                    latencySeconds: try units.convert(channel.latency, to: "s"),
                    transform: transform,
                    destinationBounds: bounds,
                    required: channel.required,
                    transactional: channel.transactional
                ))
            } else if channel.source.subsystem == .numiVivo {
                guard let index = indices[channel.source.channel] else {
                    throw VivoArtifactValidationError.unresolved(
                        "outbound coupling source \(channel.source.channel) is absent from ProgramPack species"
                    )
                }
                let transform = try unitTransform(from: species[index].unit, to: channel.destination.unit)
                let lanes = try channel.source.lanes.resolve(configuration: configuration)
                outbound.append(.init(
                    identifier: channel.identifier,
                    direction: .outbound,
                    peerSubsystem: channel.destination.subsystem,
                    peerChannel: channel.destination.channel,
                    speciesIndex: UInt32(index),
                    lanes: lanes,
                    updateMode: channel.updateMode,
                    aggregation: channel.aggregation,
                    interpolation: channel.interpolation,
                    authority: channel.authority,
                    cadenceSeconds: try units.convert(channel.cadence, to: "s"),
                    latencySeconds: try units.convert(channel.latency, to: "s"),
                    transform: transform,
                    destinationBounds: channel.destinationBounds,
                    required: channel.required,
                    transactional: channel.transactional
                ))
            } else {
                throw VivoArtifactValidationError.incompatible(
                    "coupling channel \(channel.identifier) does not have NumiVivo as source or destination"
                )
            }
        }

        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(coupling))
        return PreparedVivoCouplingPlan(
            fingerprint: fingerprint,
            programFingerprint: coupling.programFingerprint,
            hostContextFingerprint: coupling.hostContextFingerprint,
            inbound: inbound,
            outbound: outbound,
            policy: coupling.policy
        )
    }

    private func unitTransform(from source: String, to destination: String) throws -> VivoUnitTransform {
        let zero = try units.convert(0, from: source, to: destination)
        let one = try units.convert(1, from: source, to: destination)
        return VivoUnitTransform(scale: one - zero, offset: zero)
    }

    private func convertBounds(
        _ bounds: VivoInterval?,
        from sourceUnit: String,
        to destinationUnit: String
    ) throws -> VivoInterval? {
        guard let bounds else { return nil }
        let lower = try units.convert(bounds.minimum, from: sourceUnit, to: destinationUnit)
        let upper = try units.convert(bounds.maximum, from: sourceUnit, to: destinationUnit)
        return VivoInterval(minimum: min(lower, upper), maximum: max(lower, upper))
    }
}

public struct VivoCouplingValue: Codable, Sendable, Equatable {
    public let channelIdentifier: String
    public let sourceStepIndex: UInt32
    public let sourceTime: Double
    public let values: [Double]

    public init(
        channelIdentifier: String,
        sourceStepIndex: UInt32,
        sourceTime: Double,
        values: [Double]
    ) {
        self.channelIdentifier = channelIdentifier
        self.sourceStepIndex = sourceStepIndex
        self.sourceTime = sourceTime
        self.values = values
    }
}

public struct VivoCouplingFrame: Codable, Sendable, Equatable {
    public let transactionIdentifier: UUID
    public let planFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let proposedStepIndex: UInt32
    public let timeBefore: Double
    public let proposedTimeAfter: Double
    public let values: [VivoCouplingValue]

    public init(
        transactionIdentifier: UUID,
        planFingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint,
        proposedStepIndex: UInt32,
        timeBefore: Double,
        proposedTimeAfter: Double,
        values: [VivoCouplingValue]
    ) {
        self.transactionIdentifier = transactionIdentifier
        self.planFingerprint = planFingerprint
        self.programFingerprint = programFingerprint
        self.proposedStepIndex = proposedStepIndex
        self.timeBefore = timeBefore
        self.proposedTimeAfter = proposedTimeAfter
        self.values = values
    }
}

public enum VivoCouplingVote: Codable, Sendable, Equatable {
    case accept
    case reject(reason: String)
    case requestSubstep(count: UInt32, reason: String)
    case reversibleShutdown(reason: String)
    case permanentShutdown(reason: String)

    private enum CodingKeys: String, CodingKey {
        case kind
        case reason
        case count
    }

    private enum Kind: String, Codable {
        case accept
        case reject
        case requestSubstep
        case reversibleShutdown
        case permanentShutdown
    }

    public init(from decoder: Decoder) throws {
        let container = try decoder.container(keyedBy: CodingKeys.self)
        switch try container.decode(Kind.self, forKey: .kind) {
        case .accept:
            self = .accept
        case .reject:
            self = .reject(reason: try container.decode(String.self, forKey: .reason))
        case .requestSubstep:
            self = .requestSubstep(
                count: try container.decode(UInt32.self, forKey: .count),
                reason: try container.decode(String.self, forKey: .reason)
            )
        case .reversibleShutdown:
            self = .reversibleShutdown(reason: try container.decode(String.self, forKey: .reason))
        case .permanentShutdown:
            self = .permanentShutdown(reason: try container.decode(String.self, forKey: .reason))
        }
    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)
        switch self {
        case .accept:
            try container.encode(Kind.accept, forKey: .kind)
        case .reject(let reason):
            try container.encode(Kind.reject, forKey: .kind)
            try container.encode(reason, forKey: .reason)
        case .requestSubstep(let count, let reason):
            try container.encode(Kind.requestSubstep, forKey: .kind)
            try container.encode(count, forKey: .count)
            try container.encode(reason, forKey: .reason)
        case .reversibleShutdown(let reason):
            try container.encode(Kind.reversibleShutdown, forKey: .kind)
            try container.encode(reason, forKey: .reason)
        case .permanentShutdown(let reason):
            try container.encode(Kind.permanentShutdown, forKey: .kind)
            try container.encode(reason, forKey: .reason)
        }
    }
}

public protocol VivoTransactionalParticipant: Actor {
    var subsystem: VivoSubsystem { get }
    func prepare(frame: VivoCouplingFrame) async throws -> VivoCouplingVote
    func commit(transactionIdentifier: UUID) async throws
    func rollback(transactionIdentifier: UUID) async
}
