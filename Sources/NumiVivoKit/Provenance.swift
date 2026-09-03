import CryptoKit
import Foundation

public struct VivoExecutionConfigurationSnapshot: Codable, Sendable, Hashable {
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

    public init(_ configuration: VivoRuntimeConfiguration) {
        self.cellCapacity = configuration.cellCapacity
        self.activeCellCount = configuration.activeCellCount
        self.mode = configuration.mode
        self.seed = configuration.seed
        self.maximumSubstep = configuration.maximumSubstep
        self.minimumSubstep = configuration.minimumSubstep
        self.maximumInternalSubsteps = configuration.maximumInternalSubsteps
        self.eventCapacity = configuration.eventCapacity
        self.maximumPrivateBytes = configuration.maximumPrivateBytes
        self.maximumDelayBytes = configuration.maximumDelayBytes
    }
}

public struct VivoCumulativeDiagnostics: Codable, Sendable, Hashable {
    public var flags: UInt32 = 0
    public var rejectedSteps: UInt64 = 0
    public var substepRequests: UInt64 = 0
    public var nonFiniteCount: UInt64 = 0
    public var boundViolationCount: UInt64 = 0
    public var monitorViolationCount: UInt64 = 0
    public var shutdownCount: UInt64 = 0
    public var expressionFaultCount: UInt64 = 0
    public var stochasticFallbackCount: UInt64 = 0
    public var fluxTruncationCount: UInt64 = 0
    public var eventOverflowCount: UInt64 = 0

    mutating func incorporate(_ step: VivoStepResult) {
        flags |= step.diagnostics.flags
        if step.status == .rejected { rejectedSteps &+= 1 }
        if step.status == .substepRequired { substepRequests &+= 1 }
        nonFiniteCount &+= UInt64(step.diagnostics.nonFiniteCount)
        boundViolationCount &+= UInt64(step.diagnostics.boundViolationCount)
        monitorViolationCount &+= UInt64(step.diagnostics.monitorViolationCount)
        shutdownCount &+= UInt64(step.diagnostics.shutdownCount)
        expressionFaultCount &+= UInt64(step.diagnostics.expressionFaultCount)
        stochasticFallbackCount &+= UInt64(step.diagnostics.stochasticFallbackCount)
        fluxTruncationCount &+= UInt64(step.diagnostics.fluxTruncationCount)
        eventOverflowCount &+= UInt64(step.diagnostics.eventOverflowCount)
    }
}

public struct VivoExecutionCertificate: Codable, Sendable, Hashable {
    public static let schema = "numivivo.org/execution-certificate/v1"

    public var schema: String
    public var runID: UUID
    public var parentRunID: UUID?
    public var programSourceFingerprint: String
    public var programContentFingerprint: String
    public var programPackMajor: UInt16
    public var programPackMinor: UInt16
    public var compilerABI: UInt32
    public var fidelity: VivoFidelity
    public var runtimeVersion: String
    public var runtimeConfiguration: VivoExecutionConfigurationSnapshot
    public var device: VivoMetalDeviceInfo
    public var startedAt: Date
    public var finishedAt: Date
    public var requestedSteps: UInt64
    public var committedSteps: UInt64
    public var committedLogicalStep: UInt64
    public var committedSimulationTime: Double
    public var stateVersion: UInt32
    public var finalStatus: VivoStepStatus
    public var cumulativeDiagnostics: VivoCumulativeDiagnostics
    public var eventCount: UInt64
    public var eventStreamSHA256: String
    public var observationRecordCount: UInt64
    public var observationStreamSHA256: String
    public var metadata: [String: String]
}

public struct VivoCertifiedExecution: Codable, Sendable, Hashable {
    public var certificate: VivoExecutionCertificate
    public var certificateSHA256: String

    public func canonicalJSON() throws -> Data {
        let encoder = JSONEncoder()
        encoder.dateEncodingStrategy = .iso8601
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return try encoder.encode(self)
    }
}

public actor VivoExecutionRecorder {
    private let runID: UUID
    private let parentRunID: UUID?
    private let pack: VivoProgramPack
    private let configuration: VivoExecutionConfigurationSnapshot
    private let device: VivoMetalDeviceInfo
    private let runtimeVersion: String
    private let startedAt: Date
    private let metadata: [String: String]

    private var requestedSteps: UInt64 = 0
    private var committedSteps: UInt64 = 0
    private var committedLogicalStep: UInt64 = 0
    private var committedSimulationTime: Double = 0
    private var stateVersion: UInt32 = 0
    private var finalStatus: VivoStepStatus = .committed
    private var cumulativeDiagnostics = VivoCumulativeDiagnostics()
    private var eventCount: UInt64 = 0
    private var observationRecordCount: UInt64 = 0
    private var eventHasher = SHA256()
    private var observationHasher = SHA256()
    private var finalized: VivoCertifiedExecution?

    public init(
        programPack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device: VivoMetalDeviceInfo,
        runtimeVersion: String = "0.1.0",
        runID: UUID = UUID(),
        parentRunID: UUID? = nil,
        startedAt: Date = Date(),
        metadata: [String: String] = [:]
    ) {
        self.runID = runID
        self.parentRunID = parentRunID
        self.pack = programPack
        self.configuration = VivoExecutionConfigurationSnapshot(configuration)
        self.device = device
        self.runtimeVersion = runtimeVersion
        self.startedAt = startedAt
        self.metadata = metadata
    }

    public func record(step: VivoStepResult) throws {
        guard finalized == nil else {
            throw VivoRuntimeError.invalidConfiguration("Execution recorder is already finalized.")
        }
        requestedSteps &+= 1
        if step.status == .committed { committedSteps &+= 1 }
        committedLogicalStep = step.committedLogicalStep
        committedSimulationTime = step.committedTime
        stateVersion = step.stateVersion
        finalStatus = step.status
        cumulativeDiagnostics.incorporate(step)

        for event in step.events {
            eventHasher.update(data: Self.canonicalEvent(event))
            eventCount &+= 1
        }
    }

    public func record(observations: [VivoStateSlice], logicalStep: UInt64) throws {
        guard finalized == nil else {
            throw VivoRuntimeError.invalidConfiguration("Execution recorder is already finalized.")
        }
        for observation in observations.sorted(by: { $0.species < $1.species }) {
            observationHasher.update(
                data: Self.canonicalObservation(observation, logicalStep: logicalStep)
            )
            observationRecordCount &+= 1
        }
    }

    public func finalize(at finishedAt: Date = Date()) throws -> VivoCertifiedExecution {
        if let finalized { return finalized }

        let eventDigest = Self.hex(eventHasher.finalize())
        let observationDigest = Self.hex(observationHasher.finalize())
        let certificate = VivoExecutionCertificate(
            schema: VivoExecutionCertificate.schema,
            runID: runID,
            parentRunID: parentRunID,
            programSourceFingerprint: pack.header.sourceFingerprint,
            programContentFingerprint: pack.header.contentFingerprint,
            programPackMajor: pack.header.major,
            programPackMinor: pack.header.minor,
            compilerABI: pack.header.compilerABI,
            fidelity: pack.header.fidelity,
            runtimeVersion: runtimeVersion,
            runtimeConfiguration: configuration,
            device: device,
            startedAt: startedAt,
            finishedAt: finishedAt,
            requestedSteps: requestedSteps,
            committedSteps: committedSteps,
            committedLogicalStep: committedLogicalStep,
            committedSimulationTime: committedSimulationTime,
            stateVersion: stateVersion,
            finalStatus: finalStatus,
            cumulativeDiagnostics: cumulativeDiagnostics,
            eventCount: eventCount,
            eventStreamSHA256: eventDigest,
            observationRecordCount: observationRecordCount,
            observationStreamSHA256: observationDigest,
            metadata: metadata
        )
        let encoder = JSONEncoder()
        encoder.dateEncodingStrategy = .iso8601
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        let certificateBytes = try encoder.encode(certificate)
        let output = VivoCertifiedExecution(
            certificate: certificate,
            certificateSHA256: Self.hex(SHA256.hash(data: certificateBytes))
        )
        finalized = output
        return output
    }

    private static func canonicalEvent(_ event: VivoEvent) -> Data {
        var data = Data()
        data.appendLittleEndian(event.cellIndex)
        data.appendLittleEndian(event.kind)
        data.appendLittleEndian(event.subject)
        data.appendLittleEndian(event.logicalStep)
        data.appendLittleEndian(event.value0.bitPattern)
        data.appendLittleEndian(event.value1.bitPattern)
        data.appendLittleEndian(event.flags)
        return data
    }

    private static func canonicalObservation(
        _ observation: VivoStateSlice,
        logicalStep: UInt64
    ) -> Data {
        var data = Data()
        data.appendLittleEndian(logicalStep)
        let name = Data(observation.species.utf8)
        data.appendLittleEndian(UInt32(clamping: name.count))
        data.append(name)
        data.appendLittleEndian(UInt64(observation.values.count))
        for value in observation.values {
            data.appendLittleEndian(value.bitPattern)
        }
        return data
    }

    private static func hex<D: Sequence>(_ digest: D) -> String where D.Element == UInt8 {
        let alphabet = Array("0123456789abcdef".utf8)
        var output = [UInt8]()
        output.reserveCapacity(64)
        for byte in digest {
            output.append(alphabet[Int(byte >> 4)])
            output.append(alphabet[Int(byte & 0x0f)])
        }
        return String(decoding: output, as: UTF8.self)
    }
}

private extension Data {
    mutating func appendLittleEndian<T: FixedWidthInteger>(_ value: T) {
        var little = value.littleEndian
        Swift.withUnsafeBytes(of: &little) { append(contentsOf: $0) }
    }
}
