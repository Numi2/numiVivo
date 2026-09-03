import Foundation

public enum VivoFidelity: UInt32, Codable, Sendable, CaseIterable {
    case f0Logic = 0
    case f1Deterministic = 1
    case f2Stochastic = 2
    case f3Spatial = 3
    case f4Tissue = 4
}

public enum VivoStepMode: UInt32, Codable, Sendable {
    case deterministic = 1
    case stochastic = 2
}

public enum VivoStepStatus: UInt32, Codable, Sendable {
    case committed = 0
    case rejected = 1
    case substepRequired = 2
    case reversibleShutdown = 3
    case permanentShutdown = 4
}

public enum VivoInputUpdateMode: UInt32, Codable, Sendable {
    case set = 0
    case add = 1
}

public struct VivoDiagnosticFlags: OptionSet, Codable, Sendable, Hashable {
    public let rawValue: UInt32

    public init(rawValue: UInt32) {
        self.rawValue = rawValue
    }

    public static let nonFinite = Self(rawValue: 1 << 0)
    public static let boundViolation = Self(rawValue: 1 << 1)
    public static let monitor = Self(rawValue: 1 << 2)
    public static let reject = Self(rawValue: 1 << 3)
    public static let substep = Self(rawValue: 1 << 4)
    public static let reversibleShutdown = Self(rawValue: 1 << 5)
    public static let permanentShutdown = Self(rawValue: 1 << 6)
    public static let expressionFault = Self(rawValue: 1 << 7)
    public static let stochasticFallback = Self(rawValue: 1 << 8)
    public static let fluxTruncation = Self(rawValue: 1 << 9)
    public static let eventOverflow = Self(rawValue: 1 << 10)
}

public struct VivoInputUpdate: Sendable, Hashable {
    public var signal: String
    public var cellIndex: UInt32
    public var value: Float
    public var mode: VivoInputUpdateMode

    public init(
        signal: String,
        cellIndex: UInt32,
        value: Float,
        mode: VivoInputUpdateMode = .set
    ) {
        self.signal = signal
        self.cellIndex = cellIndex
        self.value = value
        self.mode = mode
    }
}

public struct VivoRuntimeConfiguration: Sendable, Hashable {
    public var cellCapacity: Int
    public var activeCellCount: Int
    public var mode: VivoStepMode
    public var seed: UInt64
    public var maximumSubstep: Float
    public var minimumSubstep: Float
    public var maximumInternalSubsteps: Int
    public var eventCapacity: Int
    public var maximumPrivateBytes: Int
    public var maximumDelayBytes: Int
    public var label: String

    public init(
        cellCapacity: Int,
        activeCellCount: Int? = nil,
        mode: VivoStepMode = .stochastic,
        seed: UInt64 = 0,
        maximumSubstep: Float = 0.05,
        minimumSubstep: Float = 0.001,
        maximumInternalSubsteps: Int = 4_096,
        eventCapacity: Int = 65_536,
        maximumPrivateBytes: Int = 8 * 1_024 * 1_024 * 1_024,
        maximumDelayBytes: Int = 2 * 1_024 * 1_024 * 1_024,
        label: String = "NumiVivoRuntime"
    ) {
        self.cellCapacity = cellCapacity
        self.activeCellCount = activeCellCount ?? cellCapacity
        self.mode = mode
        self.seed = seed
        self.maximumSubstep = maximumSubstep
        self.minimumSubstep = minimumSubstep
        self.maximumInternalSubsteps = maximumInternalSubsteps
        self.eventCapacity = eventCapacity
        self.maximumPrivateBytes = maximumPrivateBytes
        self.maximumDelayBytes = maximumDelayBytes
        self.label = label
    }
}

public struct VivoRuntimeMemoryReport: Sendable, Codable, Hashable {
    public var privateBytes: Int
    public var sharedBytes: Int
    public var stateBytes: Int
    public var temporalStateBytes: Int
    public var reactionFluxBytes: Int
    public var delayBytes: Int
    public var refractoryBytes: Int
    public var staticTableBytes: Int
    public var delaySlotCount: Int

    public var totalBytes: Int { privateBytes + sharedBytes }
}

public struct VivoStepDiagnostics: Sendable, Codable, Hashable {
    public var flags: UInt32
    public var nonFiniteCount: UInt32
    public var boundViolationCount: UInt32
    public var monitorViolationCount: UInt32
    public var shutdownCount: UInt32
    public var expressionFaultCount: UInt32
    public var stochasticFallbackCount: UInt32
    public var fluxTruncationCount: UInt32
    public var eventOverflowCount: UInt32
    public var firstCell: UInt32?
    public var firstSubject: UInt32?
    public var requestedResponse: UInt32
    public var maximumSeverity: UInt32

    public var categories: VivoDiagnosticFlags { .init(rawValue: flags) }
    public var isClean: Bool { flags == 0 }
}

public struct VivoEvent: Sendable, Codable, Hashable {
    public var cellIndex: UInt32
    public var kind: UInt32
    public var subject: UInt32
    public var logicalStep: UInt64
    public var value0: Float
    public var value1: Float
    public var flags: UInt32
}

public struct VivoStepResult: Sendable, Codable, Hashable {
    public var status: VivoStepStatus
    public var requestedDeltaTime: Float
    public var executedSubsteps: Int
    public var committedLogicalStep: UInt64
    public var committedTime: Double
    public var stateVersion: UInt32
    public var diagnostics: VivoStepDiagnostics
    public var events: [VivoEvent]
}

public struct VivoStateSlice: Sendable, Hashable {
    public var species: String
    public var values: [Float]
}

public enum VivoRuntimeError: Error, Sendable, CustomStringConvertible {
    case metalUnavailable
    case invalidConfiguration(String)
    case invalidProgramPack(String)
    case compiler(status: Int32, report: String)
    case missingShaderResource(String)
    case metalLibrary(String)
    case missingKernel(String)
    case pipeline(String)
    case allocation(String)
    case resourceLimit(String)
    case unknownSignal(String)
    case inputNotExternallyOwned(String)
    case invalidCellIndex(UInt32)
    case commandEncoding(String)
    case commandExecution(String)
    case runtimeShutdown(VivoStepStatus)

    public var description: String {
        switch self {
        case .metalUnavailable:
            return "No compatible Metal device is available."
        case .invalidConfiguration(let message),
             .invalidProgramPack(let message),
             .metalLibrary(let message),
             .pipeline(let message),
             .allocation(let message),
             .resourceLimit(let message),
             .commandEncoding(let message),
             .commandExecution(let message):
            return message
        case .compiler(let status, let report):
            return "NumiVivo compiler failed with status \(status): \(report)"
        case .missingShaderResource(let name):
            return "Missing Metal shader resource '\(name)'."
        case .missingKernel(let name):
            return "Missing Metal kernel '\(name)'."
        case .unknownSignal(let signal):
            return "Program does not declare signal '\(signal)'."
        case .inputNotExternallyOwned(let signal):
            return "Signal '\(signal)' is not externally owned and cannot be staged as an input."
        case .invalidCellIndex(let index):
            return "Cell index \(index) is outside the active cell range."
        case .runtimeShutdown(let status):
            return "Runtime is in terminal status \(status)."
        }
    }
}

// Exact host mirrors of the Metal ABI. VivoMetalProgram validates their byte
// widths before creating pipelines or argument buffers.
struct NVivoStepUniformsSwift: Sendable {
    var activeCellCount: UInt32
    var cellCapacity: UInt32
    var speciesCount: UInt32
    var parameterCount: UInt32
    var reactionCount: UInt32
    var ruleCount: UInt32
    var monitorCount: UInt32
    var temporalStateCount: UInt32
    var deltaTime: Float
    var absoluteTime: Float
    var logicalStepLow: UInt32
    var logicalStepHigh: UInt32
    var seedLow: UInt32
    var seedHigh: UInt32
    var mode: UInt32
    var substepIndex: UInt32
    var eventCapacity: UInt32
    var featureFlags: UInt32
    var delaySlotCount: UInt32
    var delayWriteSlot: UInt32
}

struct NVivoCohortUniformsSwift: Sendable {
    var reactionOffset: UInt32
    var reactionCount: UInt32
    var dispatchCellCount: UInt32
    var reserved: UInt32 = 0
}

struct NVivoInputUpdateSwift: Sendable {
    var speciesIndex: UInt32
    var cellIndex: UInt32
    var value: Float
    var mode: UInt32
}

struct NVivoDiagnosticsSwift: Sendable {
    var flags: UInt32
    var nonFiniteCount: UInt32
    var boundViolationCount: UInt32
    var monitorViolationCount: UInt32
    var shutdownCount: UInt32
    var expressionFaultCount: UInt32
    var stochasticFallbackCount: UInt32
    var fluxTruncationCount: UInt32
    var eventOverflowCount: UInt32
    var firstCell: UInt32
    var firstSubject: UInt32
    var requestedResponse: UInt32
    var maximumSeverity: UInt32
    var reserved0: UInt32
    var reserved1: UInt32
    var reserved2: UInt32

    var publicValue: VivoStepDiagnostics {
        VivoStepDiagnostics(
            flags: flags,
            nonFiniteCount: nonFiniteCount,
            boundViolationCount: boundViolationCount,
            monitorViolationCount: monitorViolationCount,
            shutdownCount: shutdownCount,
            expressionFaultCount: expressionFaultCount,
            stochasticFallbackCount: stochasticFallbackCount,
            fluxTruncationCount: fluxTruncationCount,
            eventOverflowCount: eventOverflowCount,
            firstCell: firstCell == UInt32.max ? nil : firstCell,
            firstSubject: firstSubject == UInt32.max ? nil : firstSubject,
            requestedResponse: requestedResponse,
            maximumSeverity: maximumSeverity
        )
    }
}

struct NVivoPublicationSwift: Sendable {
    var committedStepLow: UInt32
    var committedStepHigh: UInt32
    var stateVersion: UInt32
    var status: UInt32
    var diagnosticFlags: UInt32
    var shutdownState: UInt32
    var activeCellCount: UInt32
    var eventCount: UInt32
}

struct NVivoEventSwift: Sendable {
    var cellIndex: UInt32
    var kind: UInt32
    var subject: UInt32
    var logicalStepLow: UInt32
    var value0: Float
    var value1: Float
    var flags: UInt32
    var logicalStepHigh: UInt32
}

extension UInt64 {
    var nvivoLow: UInt32 { UInt32(truncatingIfNeeded: self) }
    var nvivoHigh: UInt32 { UInt32(truncatingIfNeeded: self >> 32) }

    init(nvivoLow low: UInt32, high: UInt32) {
        self = UInt64(low) | (UInt64(high) << 32)
    }
}
