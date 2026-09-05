import Foundation
import Metal

public enum VivoCouplingMode: UInt32, Sendable, Codable {
    case replace = 0
    case add = 1
    case relaxHalfway = 2
}
public struct VivoCouplingUpdate: Sendable, Codable, Equatable {
    public let speciesIndex: UInt32
    public let laneIndex: UInt32
    public let mode: VivoCouplingMode
    public let value: Float
    public init(speciesIndex: UInt32, laneIndex: UInt32, mode: VivoCouplingMode = .replace, value: Float) {
        self.speciesIndex = speciesIndex
        self.laneIndex = laneIndex
        self.mode = mode
        self.value = value
    }
    var abi: VivoCouplingUpdateABI {
        .init(speciesIndex: speciesIndex, laneIndex: laneIndex, mode: mode.rawValue, value: value)
    }
}
public struct VivoPublicationRequest: Sendable, Codable, Equatable {
    public let speciesIndex: UInt32
    public let laneIndex: UInt32
    public let flags: UInt32
    public init(speciesIndex: UInt32, laneIndex: UInt32, flags: UInt32 = 0) {
        self.speciesIndex = speciesIndex
        self.laneIndex = laneIndex
        self.flags = flags
    }
}
public struct VivoStepRequest: Sendable, Codable, Equatable {
    public var timeStep: Float?
    public var coupling: [VivoCouplingUpdate]
    public var publications: [VivoPublicationRequest]
    public var permitAdaptiveReduction: Bool
    public init(timeStep: Float? = nil, coupling: [VivoCouplingUpdate] = [],
                publications: [VivoPublicationRequest] = [], permitAdaptiveReduction: Bool = true) {
        self.timeStep = timeStep
        self.coupling = coupling
        self.publications = publications
        self.permitAdaptiveReduction = permitAdaptiveReduction
    }
}
public enum VivoStepDisposition: String, Sendable, Codable {
    case committed
    case committedWithReducedStep
    case rejected
    case reversibleShutdown
    case permanentShutdown
}
public struct VivoStepCertificate: Sendable, Codable, Equatable {
    public let disposition: VivoStepDisposition
    public let sourceFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let fidelity: VivoFidelity
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let stepIndex: UInt32
    public let attemptedStep: Float
    public let acceptedStep: Float?
    public let attemptCount: UInt32
    public let timeBefore: Float
    public let timeAfter: Float
    public let status: VivoRuntimeStatus
    public var committed: Bool { disposition == .committed || disposition == .committedWithReducedStep }
}
public struct VivoStepResult: Sendable, Codable, Equatable {
    public let certificate: VivoStepCertificate
    public let events: [VivoEvent]
    public let publications: [Float]
}
public struct VivoStateSnapshot: Sendable, Codable, Equatable {
    public let sourceFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let absoluteTime: Float
    public let speciesCount: UInt32
    public let laneCount: UInt32
    public let values: [Float]
    public func value(species: UInt32, lane: UInt32) -> Float? {
        guard species < speciesCount, lane < laneCount else { return nil }
        let index = UInt64(species) * UInt64(laneCount) + UInt64(lane)
        guard index < UInt64(values.count) else { return nil }
        return values[Int(index)]
    }
}

/// Stable public facade over the single authoritative ProgramPack transaction
/// engine. It does not maintain another arena, clock, stochastic stream or set
/// of Metal bindings that can diverge from the coupling runtime.
public actor VivoRuntime {
    public typealias Lifecycle = VivoTransactionalMolecularRuntime.Lifecycle
    public nonisolated let pack: VivoProgramPack
    public nonisolated let configuration: VivoRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities
    private let engine: VivoTransactionalMolecularRuntime

    public static func make(pack: VivoProgramPack, configuration: VivoRuntimeConfiguration,
                            device: MTLDevice? = nil) async throws -> VivoRuntime {
        let engine = try await VivoTransactionalMolecularRuntime.make(pack: pack, configuration: configuration, device: device)
        return VivoRuntime(pack: pack, configuration: configuration, engine: engine)
    }
    private init(pack: VivoProgramPack, configuration: VivoRuntimeConfiguration, engine: VivoTransactionalMolecularRuntime) {
        self.pack = pack
        self.configuration = configuration
        self.engine = engine
        capabilities = engine.capabilities
    }
    public func lifecycle() async -> Lifecycle { await engine.lifecycle() }
    public func time() async -> Float { await engine.time() }
    public func timeSeconds() async -> Double { await engine.timeSeconds() }
    public func stepIndex() async -> UInt32 { await engine.stepIndex() }
    public func hasPendingTransaction() async -> Bool { await engine.hasPendingTransaction() }
    public func resume() async throws { try await engine.resume() }
    public func stopPermanently(reason: String) async { await engine.stopPermanently(reason: reason) }
    public func setTransport(_ values: [VivoSpeciesTransportABI]) async throws { try await engine.setTransport(values) }
    public func setVelocity(_ values: [SIMD4<Float>]) async throws { try await engine.setVelocity(values) }
    public func setVolumeFractions(_ values: [Float]) async throws { try await engine.setVolumeFractions(values) }
    public func apply(_ context: PreparedVivoHostContext) async throws { try await engine.apply(context) }
    public func apply(intervention: PreparedVivoIntervention.Operation) async throws { try await engine.apply(intervention: intervention) }
    /// Initial coupling is applied atomically by apply; no uncommitted facade queue exists.
    public func primePendingInputs() async throws { try await engine.requireAcceptedBoundary() }
    public func snapshot() async throws -> VivoStateSnapshot { try await engine.snapshot() }
    public func checkpoint() async throws -> VivoMolecularCheckpoint { try await engine.checkpoint() }
    public func resumeCheckpoint() async throws -> VivoMolecularResumeCheckpoint { try await engine.resumeCheckpoint() }
    public func restore(_ checkpoint: VivoMolecularCheckpoint) async throws { try await engine.restore(checkpoint) }
    public func restore(_ checkpoint: VivoMolecularResumeCheckpoint) async throws { try await engine.restore(checkpoint) }

    public func step(_ request: VivoStepRequest = .init()) async throws -> VivoStepResult {
        let result = try await engine.step(request)
        let certificate = result.certificate
        return .init(certificate: .init(disposition: certificate.disposition,
                                        sourceFingerprint: certificate.sourceFingerprint,
                                        programFingerprint: certificate.programFingerprint,
                                        fidelity: certificate.fidelity, deviceName: certificate.deviceName,
                                        deviceRegistryID: certificate.deviceRegistryID, stepIndex: certificate.stepIndex,
                                        attemptedStep: certificate.requestedTimeStep, acceptedStep: certificate.acceptedTimeStep,
                                        attemptCount: certificate.attemptCount, timeBefore: certificate.timeBefore,
                                        timeAfter: certificate.timeAfter, status: certificate.status),
                     events: result.events, publications: result.publications)
    }
}
