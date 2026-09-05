import Foundation
@preconcurrency import Metal

public struct VivoSpatialGrid: Sendable, Codable, Equatable {
    public enum Boundary: UInt32, Sendable, Codable {
        case noFlux = 0
        case periodic = 1
        case absorbing = 2
    }

    public var width: UInt32
    public var height: UInt32
    public var depth: UInt32
    public var spacingX: Float
    public var spacingY: Float
    public var spacingZ: Float
    public var boundary: Boundary

    public init(
        width: UInt32,
        height: UInt32,
        depth: UInt32,
        spacingX: Float,
        spacingY: Float,
        spacingZ: Float,
        boundary: Boundary = .noFlux
    ) {
        self.width = width
        self.height = height
        self.depth = depth
        self.spacingX = spacingX
        self.spacingY = spacingY
        self.spacingZ = spacingZ
        self.boundary = boundary
    }

    public var voxelCount: UInt32 {
        let xy = width.multipliedReportingOverflow(by: height)
        guard !xy.overflow else { return 0 }
        let xyz = xy.partialValue.multipliedReportingOverflow(by: depth)
        return xyz.overflow ? 0 : xyz.partialValue
    }

    public func validate() throws {
        guard width > 0, height > 0, depth > 0, voxelCount > 0 else {
            throw VivoRuntimeError.invalidConfiguration("spatial dimensions must be positive and must not overflow UInt32")
        }
        guard spacingX.isFinite, spacingY.isFinite, spacingZ.isFinite,
              spacingX > 0, spacingY > 0, spacingZ > 0 else {
            throw VivoRuntimeError.invalidConfiguration("spatial spacing must be finite and positive")
        }
    }
}

public struct VivoRuntimeSeed: Sendable, Codable, Equatable, Hashable {
    public var low: UInt64
    public var high: UInt64

    public init(low: UInt64, high: UInt64) {
        self.low = low
        self.high = high
    }

    public static let deterministic = VivoRuntimeSeed(
        low: 0x4e756d695669766f,
        high: 0x4d6f6c6563756c65
    )
}

public struct VivoRuntimeConfiguration: Sendable, Codable, Equatable {
    public var fidelity: VivoFidelity
    public var environmentCount: UInt32
    public var spatialGrid: VivoSpatialGrid?
    public var timeStep: Float
    public var minimumTimeStep: Float
    public var maximumTimeStep: Float
    public var maximumSubsteps: UInt32
    public var eventCapacity: UInt32
    public var seed: VivoRuntimeSeed
    public var privateHeapHeadroom: Double

    public init(
        fidelity: VivoFidelity,
        environmentCount: UInt32 = 1,
        spatialGrid: VivoSpatialGrid? = nil,
        timeStep: Float = 0.01,
        minimumTimeStep: Float = 0.000_001,
        maximumTimeStep: Float = 1,
        maximumSubsteps: UInt32 = 16,
        eventCapacity: UInt32 = 16_384,
        seed: VivoRuntimeSeed = .deterministic,
        privateHeapHeadroom: Double = 1.15
    ) {
        self.fidelity = fidelity
        self.environmentCount = environmentCount
        self.spatialGrid = spatialGrid
        self.timeStep = timeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.maximumSubsteps = maximumSubsteps
        self.eventCapacity = eventCapacity
        self.seed = seed
        self.privateHeapHeadroom = privateHeapHeadroom
    }

    public var voxelCount: UInt32 { spatialGrid?.voxelCount ?? 1 }

    public var laneCount: UInt32 {
        let product = environmentCount.multipliedReportingOverflow(by: voxelCount)
        return product.overflow ? 0 : product.partialValue
    }

    public func validate(for pack: VivoProgramPack) throws {
        guard environmentCount > 0, laneCount > 0 else {
            throw VivoRuntimeError.invalidConfiguration("environmentCount × voxelCount must be positive and must not overflow UInt32")
        }
        guard timeStep.isFinite, minimumTimeStep.isFinite, maximumTimeStep.isFinite,
              minimumTimeStep > 0, maximumTimeStep >= minimumTimeStep,
              timeStep >= minimumTimeStep, timeStep <= maximumTimeStep else {
            throw VivoRuntimeError.invalidConfiguration("time-step bounds are inconsistent or non-finite")
        }
        guard maximumSubsteps > 0 else {
            throw VivoRuntimeError.invalidConfiguration("maximumSubsteps must be positive")
        }
        guard eventCapacity > 0 else {
            throw VivoRuntimeError.invalidConfiguration("eventCapacity must be positive")
        }
        guard privateHeapHeadroom.isFinite, privateHeapHeadroom >= 1 else {
            throw VivoRuntimeError.invalidConfiguration("privateHeapHeadroom must be finite and at least 1.0")
        }
        if let spatialGrid { try spatialGrid.validate() }
        if fidelity.rawValue < pack.header.fidelity.rawValue {
            throw VivoRuntimeError.invalidConfiguration(
                "runtime fidelity \(fidelity.label) is below ProgramPack fidelity \(pack.header.fidelity.label)"
            )
        }
        if fidelity.rawValue >= VivoFidelity.spatial.rawValue, spatialGrid == nil {
            throw VivoRuntimeError.invalidConfiguration("F3 and F4 runtimes require a spatial grid")
        }
    }
}

public enum VivoRuntimeError: Error, Sendable, CustomStringConvertible {
    case metalUnavailable
    case commandQueueUnavailable
    case invalidConfiguration(String)
    case incompatibleDevice(String)
    case allocationFailed(String)
    case pipelineEncodingFailed(String)
    case commandFailed(String)
    case packError(String)
    case runtimeStopped(String)

    public var description: String {
        switch self {
        case .metalUnavailable: "No Metal device is available."
        case .commandQueueUnavailable: "The Metal device could not create a command queue."
        case .invalidConfiguration(let message): "Invalid NumiVivo runtime configuration: \(message)"
        case .incompatibleDevice(let message): "Metal device is incompatible: \(message)"
        case .allocationFailed(let message): "NumiVivo allocation failed: \(message)"
        case .pipelineEncodingFailed(let message): "NumiVivo pipeline encoding failed: \(message)"
        case .commandFailed(let message): "NumiVivo command failed: \(message)"
        case .packError(let message): "NumiVivo ProgramPack error: \(message)"
        case .runtimeStopped(let message): "NumiVivo runtime is stopped: \(message)"
        }
    }
}

@frozen
public struct VivoRuntimeCommandABI: Sendable {
    public static let abiVersion: UInt32 = 1

    public var abiVersion: UInt32
    public var fidelity: UInt32
    public var laneCount: UInt32
    public var speciesCount: UInt32
    public var parameterCount: UInt32
    public var parameterStride: UInt32
    public var reactionCount: UInt32
    public var expressionCount: UInt32
    public var ruleCount: UInt32
    public var monitorCount: UInt32
    public var temporalStateCount: UInt32
    public var stepIndex: UInt32
    public var substepIndex: UInt32
    public var flags: UInt32
    public var voxelCount: UInt32
    public var eventCapacity: UInt32

    public var dt: Float
    public var absoluteTime: Float
    public var minimumDt: Float
    public var maximumDt: Float

    public var gridX: UInt32
    public var gridY: UInt32
    public var gridZ: UInt32
    public var boundaryMode: UInt32

    public var spacingX: Float
    public var spacingY: Float
    public var spacingZ: Float
    public var reservedFloat: Float

    public var seedLow: UInt64
    public var seedHigh: UInt64

    public init(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        stepIndex: UInt32,
        substepIndex: UInt32,
        dt: Float,
        absoluteTime: Float
    ) throws {
        let grid = configuration.spatialGrid
        self.abiVersion = Self.abiVersion
        self.fidelity = configuration.fidelity.rawValue
        self.laneCount = configuration.laneCount
        self.speciesCount = pack.runtimeContract.speciesCount
        self.parameterCount = pack.runtimeContract.parameterCount
        self.parameterStride = configuration.environmentCount
        self.reactionCount = pack.runtimeContract.reactionCount
        self.expressionCount = try pack.section(.expressions).count
        self.ruleCount = pack.runtimeContract.ruleCount
        self.monitorCount = pack.runtimeContract.monitorCount
        self.temporalStateCount = pack.runtimeContract.temporalStateCount
        self.stepIndex = stepIndex
        self.substepIndex = substepIndex
        self.flags = 0
        self.voxelCount = configuration.voxelCount
        self.eventCapacity = configuration.eventCapacity
        self.dt = dt
        self.absoluteTime = absoluteTime
        self.minimumDt = configuration.minimumTimeStep
        self.maximumDt = configuration.maximumTimeStep
        self.gridX = grid?.width ?? 1
        self.gridY = grid?.height ?? 1
        self.gridZ = grid?.depth ?? 1
        self.boundaryMode = grid?.boundary.rawValue ?? VivoSpatialGrid.Boundary.noFlux.rawValue
        self.spacingX = grid?.spacingX ?? 1
        self.spacingY = grid?.spacingY ?? 1
        self.spacingZ = grid?.spacingZ ?? 1
        self.reservedFloat = 0
        self.seedLow = configuration.seed.low
        self.seedHigh = configuration.seed.high
    }

    public static func validateMemoryLayout() throws {
        guard MemoryLayout<Self>.size == 128, MemoryLayout<Self>.stride == 128 else {
            throw VivoRuntimeError.incompatibleDevice(
                "Swift runtime command layout is \(MemoryLayout<Self>.size)/\(MemoryLayout<Self>.stride), expected 128/128"
            )
        }
    }
}

@frozen
public struct VivoRuntimeStatusABI: Sendable {
    public var flags: UInt32 = 0
    public var violationCount: UInt32 = 0
    public var firstLane: UInt32 = .max
    public var firstSpecies: UInt32 = .max
    public var firstMonitor: UInt32 = .max
    public var requestedResponse: UInt32 = 0
    public var requiredSubsteps: UInt32 = 1
    public var eventCount: UInt32 = 0
    public var eventDropped: UInt32 = 0
    public var maximumViolationBits: UInt32 = 0
    public var maximumRateBits: UInt32 = 0
    public var reservedAtomic: UInt32 = 0
    public var reserved0: UInt32 = 0
    public var reserved1: UInt32 = 0
    public var reserved2: UInt32 = 0
    public var reserved3: UInt32 = 0

    public init() {}

    public static func validateMemoryLayout() throws {
        guard MemoryLayout<Self>.size == 64, MemoryLayout<Self>.stride == 64 else {
            throw VivoRuntimeError.incompatibleDevice(
                "Swift runtime status layout is \(MemoryLayout<Self>.size)/\(MemoryLayout<Self>.stride), expected 64/64"
            )
        }
    }
}

public struct VivoRuntimeStatus: Sendable, Codable, Equatable {
    public struct Flags: OptionSet, Sendable, Codable, Hashable {
        public let rawValue: UInt32
        public init(rawValue: UInt32) { self.rawValue = rawValue }

        public static let nonFinite = Flags(rawValue: 1 << 0)
        public static let bounds = Flags(rawValue: 1 << 1)
        public static let monitor = Flags(rawValue: 1 << 2)
        public static let reject = Flags(rawValue: 1 << 3)
        public static let substep = Flags(rawValue: 1 << 4)
        public static let reversibleShutdown = Flags(rawValue: 1 << 5)
        public static let permanentShutdown = Flags(rawValue: 1 << 6)
        public static let randomSaturated = Flags(rawValue: 1 << 7)
        public static let eventOverflow = Flags(rawValue: 1 << 8)
        public static let invalidRate = Flags(rawValue: 1 << 9)
        public static let abiMismatch = Flags(rawValue: 1 << 10)
    }

    public let flags: Flags
    public let violationCount: UInt32
    public let firstLane: UInt32?
    public let firstSpecies: UInt32?
    public let firstMonitor: UInt32?
    public let requestedResponse: UInt32
    public let requiredSubsteps: UInt32
    public let eventCount: UInt32
    public let eventDropped: UInt32
    public let maximumViolation: Float
    public let maximumRate: Float

    init(raw: VivoRuntimeStatusABI) {
        self.flags = Flags(rawValue: raw.flags)
        self.violationCount = raw.violationCount
        self.firstLane = raw.firstLane == .max ? nil : raw.firstLane
        self.firstSpecies = raw.firstSpecies == .max ? nil : raw.firstSpecies
        self.firstMonitor = raw.firstMonitor == .max ? nil : raw.firstMonitor
        self.requestedResponse = raw.requestedResponse
        self.requiredSubsteps = max(raw.requiredSubsteps, 1)
        self.eventCount = raw.eventCount
        self.eventDropped = raw.eventDropped
        self.maximumViolation = Float(bitPattern: raw.maximumViolationBits)
        self.maximumRate = Float(bitPattern: raw.maximumRateBits)
    }

    public var blocksCommit: Bool {
        !flags.intersection([.nonFinite, .bounds, .reject, .reversibleShutdown, .permanentShutdown, .abiMismatch]).isEmpty
    }
}

@frozen
public struct VivoSpeciesTransportABI: Sendable {
    public var diffusion: Float
    public var membranePermeability: Float
    public var decayRate: Float
    public var flags: UInt32

    public init(diffusion: Float = 0, membranePermeability: Float = 0, decayRate: Float = 0, flags: UInt32 = 0) {
        self.diffusion = diffusion
        self.membranePermeability = membranePermeability
        self.decayRate = decayRate
        self.flags = flags
    }
}

@frozen
public struct VivoCouplingUpdateABI: Sendable {
    public var speciesIndex: UInt32
    public var laneIndex: UInt32
    public var mode: UInt32
    public var value: Float

    public init(speciesIndex: UInt32, laneIndex: UInt32, mode: UInt32 = 0, value: Float) {
        self.speciesIndex = speciesIndex
        self.laneIndex = laneIndex
        self.mode = mode
        self.value = value
    }
}

@frozen
public struct VivoPublicationRequestABI: Sendable {
    public var speciesIndex: UInt32
    public var laneIndex: UInt32
    public var outputIndex: UInt32
    public var flags: UInt32

    public init(speciesIndex: UInt32, laneIndex: UInt32, outputIndex: UInt32, flags: UInt32 = 0) {
        self.speciesIndex = speciesIndex
        self.laneIndex = laneIndex
        self.outputIndex = outputIndex
        self.flags = flags
    }
}

@frozen
public struct VivoEventABI: Sendable {
    public var laneIndex: UInt32
    public var stepIndex: UInt32
    public var actionKind: UInt32
    public var subject: UInt32
    public var value: Float
    public var absoluteTime: Float
    public var flags: UInt32
    public var reserved: UInt32

    public init(
        laneIndex: UInt32 = 0,
        stepIndex: UInt32 = 0,
        actionKind: UInt32 = 0,
        subject: UInt32 = 0,
        value: Float = 0,
        absoluteTime: Float = 0,
        flags: UInt32 = 0,
        reserved: UInt32 = 0
    ) {
        self.laneIndex = laneIndex
        self.stepIndex = stepIndex
        self.actionKind = actionKind
        self.subject = subject
        self.value = value
        self.absoluteTime = absoluteTime
        self.flags = flags
        self.reserved = reserved
    }
}

public struct VivoEvent: Sendable, Codable, Equatable {
    public let laneIndex: UInt32
    public let stepIndex: UInt32
    public let actionKind: UInt32
    public let subject: UInt32
    public let value: Float
    public let absoluteTime: Float
    public let flags: UInt32

    init(_ raw: VivoEventABI) {
        laneIndex = raw.laneIndex
        stepIndex = raw.stepIndex
        actionKind = raw.actionKind
        subject = raw.subject
        value = raw.value
        absoluteTime = raw.absoluteTime
        flags = raw.flags
    }
}
