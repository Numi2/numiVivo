import Foundation
@preconcurrency import Metal

public struct VivoMolecularPhysiologyStepRequest: Codable, Sendable, Equatable {
    public var timeStepSeconds: Double?
    public var molecularCoupling: [VivoCouplingUpdate]
    public var molecularPublications: [VivoPublicationRequest]
    public var physiologyPreUpdates: [VivoPhysiologyStateUpdate]
    public var physiologyPostUpdates: [VivoPhysiologyStateUpdate]
    public var physiologyPublications: [VivoPhysiologyPublicationRequest]
    public var permitAdaptiveReduction: Bool
    public init(timeStepSeconds: Double? = nil, molecularCoupling: [VivoCouplingUpdate] = [],
                molecularPublications: [VivoPublicationRequest] = [], physiologyPreUpdates: [VivoPhysiologyStateUpdate] = [],
                physiologyPostUpdates: [VivoPhysiologyStateUpdate] = [], physiologyPublications: [VivoPhysiologyPublicationRequest] = [],
                permitAdaptiveReduction: Bool = true) {
        self.timeStepSeconds = timeStepSeconds; self.molecularCoupling = molecularCoupling
        self.molecularPublications = molecularPublications; self.physiologyPreUpdates = physiologyPreUpdates
        self.physiologyPostUpdates = physiologyPostUpdates; self.physiologyPublications = physiologyPublications
        self.permitAdaptiveReduction = permitAdaptiveReduction
    }
}
public enum VivoMolecularPhysiologyDisposition: String, Codable, Sendable {
    case committed, committedWithReducedStep, committedWithoutConvergence, rejected, reversibleShutdown, permanentShutdown
}
public struct VivoMolecularPhysiologyCertificate: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoMolecularPhysiologyDisposition
    public let bridgeFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let physiologyFingerprint: VivoFingerprint
    /// The full prepared model/bridge/configuration identity, not only the source
    /// fingerprints embedded by their compilers. Nil only for legacy records.
    public let executableFingerprint: VivoFingerprint?
    public let stepIndex: UInt32
    public let timeBeforeSeconds: Double
    public let timeAfterSeconds: Double
    public let requestedTimeStepSeconds: Double
    public let acceptedTimeStepSeconds: Double?
    public let stepNegotiationCount: UInt32
    public let fixedPointIterationCount: UInt32
    public let finalResidual: Double
    public let converged: Bool
    public let molecularStatus: VivoRuntimeStatus?
    public let physiologyStatus: VivoPhysiologyRuntimeStatus?
    public let message: String
    public var committed: Bool {
        disposition == .committed || disposition == .committedWithReducedStep || disposition == .committedWithoutConvergence
    }
}
public struct VivoMolecularPhysiologyStepResult: Codable, Sendable, Equatable {
    public let certificate: VivoMolecularPhysiologyCertificate
    public let molecularCertificate: VivoMolecularTransactionCertificate?
    public let physiologyCertificate: VivoPhysiologyStepCertificate?
    public let molecularEvents: [VivoEvent]
    public let bridgeExposureValues: [Float]
    public let bridgeFeedbackValues: [Float]
    public let molecularPublications: [Float]
    public let physiologyPublications: [Float]
}
public struct VivoMolecularPhysiologyCheckpoint: Codable, Sendable, Equatable {
    public let schemaVersion: UInt32
    public let executableFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let timeSeconds: Double
    public let previousFeedbackValues: [Float]
    public let hasAcceptedFeedback: Bool
    public let molecular: VivoMolecularResumeCheckpoint
    public let physiology: VivoPhysiologyCheckpoint
    public func fingerprint() throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self)) }
}

private struct CouplingExecutionIdentity: Encodable {
    let numericalVersion: UInt32
    let program: VivoFingerprint
    let physiology: PreparedVivoPhysiologyModel
    let bridge: PreparedVivoMolecularPhysiologyBridge
    let molecularConfiguration: VivoRuntimeConfiguration
    let physiologyConfiguration: VivoPhysiologyRuntimeConfiguration
}

/// Owns its participants exclusively. Public observations are blocked throughout
/// preparation AND release. Release is logically atomic to this API, not a
/// crash-atomic distributed transaction. Unexpected release failure poisons the
/// coordinator rather than publishing a mixed state or claiming rollback.
public actor VivoMolecularPhysiologyCoordinator {
    public nonisolated let bridge: PreparedVivoMolecularPhysiologyBridge
    public nonisolated let programPack: VivoProgramPack
    public nonisolated let physiologyModel: PreparedVivoPhysiologyModel
    public nonisolated let molecularConfiguration: VivoRuntimeConfiguration
    public nonisolated let physiologyConfiguration: VivoPhysiologyRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities
    public nonisolated let executableFingerprint: VivoFingerprint
    private let device: MTLDevice
    private var molecular: VivoTransactionalMolecularRuntime
    private var physiology: VivoPhysiologyRuntime
    private var acceptedStepIndex: UInt32 = 0
    private var acceptedTime: Double = 0
    private var previousFeedbackValues: [Float]
    private var hasAcceptedFeedback = false
    private var inFlight = false
    private var failure: String?

    public static func make(programPack: VivoProgramPack, molecularConfiguration: VivoRuntimeConfiguration,
                            physiologyModel: PreparedVivoPhysiologyModel,
                            physiologyConfiguration: VivoPhysiologyRuntimeConfiguration = .init(),
                            bridge: PreparedVivoMolecularPhysiologyBridge, device requestedDevice: MTLDevice? = nil) async throws -> VivoMolecularPhysiologyCoordinator {
        try VivoProgramExecutionContract.validate(pack: programPack, configuration: molecularConfiguration)
        try VivoPhysiologyModelValidator.validate(physiologyModel)
        try physiologyConfiguration.validate(for: physiologyModel)
        try bridge.validate(program: programPack, physiology: physiologyModel, molecularLaneCount: molecularConfiguration.laneCount)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        // Sequential allocation lets the second allocator observe the first
        // participant's footprint. The overall working-set budget is not free.
        let molecular = try await VivoTransactionalMolecularRuntime.make(pack: programPack, configuration: molecularConfiguration, device: device)
        let physiology = try await VivoPhysiologyRuntime.make(model: physiologyModel, configuration: physiologyConfiguration, device: device)
        let identity = CouplingExecutionIdentity(numericalVersion: 2, program: programPack.header.contentFingerprint,
                                                  physiology: physiologyModel, bridge: bridge,
                                                  molecularConfiguration: molecularConfiguration, physiologyConfiguration: physiologyConfiguration)
        return try .init(program: programPack, molecularConfiguration: molecularConfiguration, physiologyModel: physiologyModel,
                          physiologyConfiguration: physiologyConfiguration, bridge: bridge, device: device,
                          fingerprint: VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(identity)), molecular: molecular, physiology: physiology)
    }
    private init(program: VivoProgramPack, molecularConfiguration: VivoRuntimeConfiguration,
                 physiologyModel: PreparedVivoPhysiologyModel, physiologyConfiguration: VivoPhysiologyRuntimeConfiguration,
                 bridge: PreparedVivoMolecularPhysiologyBridge, device: MTLDevice, fingerprint: VivoFingerprint,
                 molecular: VivoTransactionalMolecularRuntime, physiology: VivoPhysiologyRuntime) {
        programPack = program; self.molecularConfiguration = molecularConfiguration; self.physiologyModel = physiologyModel
        self.physiologyConfiguration = physiologyConfiguration; self.bridge = bridge; self.device = device
        self.molecular = molecular; self.physiology = physiology; executableFingerprint = fingerprint
        capabilities = VivoMetalCapabilities(device: device)
        previousFeedbackValues = Array(repeating: 0, count: bridge.molecularToPhysiology.count)
    }
    private func requireIdle() throws {
        if let failure { throw VivoRuntimeError.runtimeStopped(failure) }
        guard !inFlight else { throw VivoRuntimeError.invalidConfiguration("coupling operation is already in flight") }
    }
    public func stepIndex() -> UInt32 { acceptedStepIndex }
    public func timeSeconds() throws -> Double { try requireIdle(); return acceptedTime }
    public func molecularSnapshot() async throws -> VivoStateSnapshot {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        return try await molecular.snapshot()
    }
    public func physiologySnapshot() async throws -> VivoPhysiologySnapshot {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        return try await physiology.snapshot()
    }
    public func physiologyCheckpoint() async throws -> VivoPhysiologyCheckpoint {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        return try await physiology.checkpoint()
    }
    public func checkpoint() async throws -> VivoMolecularPhysiologyCheckpoint {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        let molecular = try await molecular.resumeCheckpoint()
        let physiology = try await physiology.checkpoint()
        guard molecular.state.absoluteTimeSeconds == acceptedTime, physiology.absoluteTimeSeconds == acceptedTime,
              molecular.state.stepIndex == acceptedStepIndex, physiology.stepIndex == acceptedStepIndex else {
            throw VivoRuntimeError.invalidConfiguration("coupling checkpoint participants disagree with the accepted generation")
        }
        return .init(schemaVersion: 1, executableFingerprint: executableFingerprint, stepIndex: acceptedStepIndex,
                      timeSeconds: acceptedTime, previousFeedbackValues: previousFeedbackValues, hasAcceptedFeedback: hasAcceptedFeedback,
                      molecular: molecular, physiology: physiology)
    }
    public func restore(_ checkpoint: VivoMolecularPhysiologyCheckpoint) async throws {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        guard checkpoint.schemaVersion == 1, checkpoint.executableFingerprint == executableFingerprint,
              checkpoint.timeSeconds.isFinite, checkpoint.timeSeconds >= 0,
              checkpoint.molecular.state.absoluteTimeSeconds == checkpoint.timeSeconds,
              checkpoint.physiology.absoluteTimeSeconds == checkpoint.timeSeconds,
              checkpoint.molecular.state.stepIndex == checkpoint.stepIndex, checkpoint.physiology.stepIndex == checkpoint.stepIndex,
              checkpoint.previousFeedbackValues.count == bridge.molecularToPhysiology.count,
              checkpoint.previousFeedbackValues.allSatisfy(\.isFinite) else {
            throw VivoRuntimeError.invalidConfiguration("coupling checkpoint execution identity, clock, state or feedback mismatch")
        }
        try Task.checkCancellation()
        let nextMolecular = try await VivoTransactionalMolecularRuntime.make(pack: programPack, configuration: molecularConfiguration, device: device)
        try await nextMolecular.restore(checkpoint.molecular)
        let nextPhysiology = try await VivoPhysiologyRuntime.make(model: physiologyModel, configuration: physiologyConfiguration, device: device)
        try await nextPhysiology.restore(checkpoint.physiology)
        try Task.checkCancellation()
        // All fallible construction/upload precedes this non-suspending owner swap.
        molecular = nextMolecular; physiology = nextPhysiology
        acceptedStepIndex = checkpoint.stepIndex; acceptedTime = checkpoint.timeSeconds
        previousFeedbackValues = checkpoint.previousFeedbackValues; hasAcceptedFeedback = checkpoint.hasAcceptedFeedback
    }

    public func step(_ request: VivoMolecularPhysiologyStepRequest = .init()) async throws -> VivoMolecularPhysiologyStepResult {
        try requireIdle()
        inFlight = true
        defer { inFlight = false }
        try Task.checkCancellation()
        guard acceptedStepIndex < UInt32.max else { throw VivoRuntimeError.invalidConfiguration("coupling step index exhausted") }
        let id = UUID(), before = acceptedTime
        let minimum = max(Double(molecularConfiguration.minimumTimeStep), physiologyModel.minimumTimeStepSeconds)
        let maximum = min(Double(molecularConfiguration.maximumTimeStep), physiologyModel.maximumTimeStepSeconds)
        let requested = request.timeStepSeconds ?? min(Double(molecularConfiguration.timeStep), physiologyModel.preferredTimeStepSeconds)
        guard requested.isFinite, requested >= minimum, requested <= maximum else { throw VivoRuntimeError.invalidConfiguration("invalid common timestep") }
        guard request.molecularCoupling.count <= 65_536, request.molecularPublications.count <= 65_536,
              request.physiologyPreUpdates.count <= physiologyConfiguration.maximumTransformsPerStep,
              request.physiologyPostUpdates.count <= physiologyConfiguration.maximumTransformsPerStep else {
            throw VivoRuntimeError.invalidConfiguration("coupling request exceeds boundary capacities")
        }
        let moleculeTime = await molecular.timeSeconds(), physiologyTime = await physiology.time()
        guard moleculeTime == before, physiologyTime == before else {
            failure = "participant clocks diverged from the accepted generation"
            throw VivoRuntimeError.runtimeStopped(failure!)
        }
        // A realized stochastic failure must not select a shorter random horizon.
        // The molecular participant also receives permitAdaptiveReduction=false.
        let mayReduce = request.permitAdaptiveReduction && molecularConfiguration.fidelity != .stochastic
        var dt = try canonicalStep(requested, minimum: minimum)
        var lastMolecular: VivoPreparedMolecularStep?, lastPhysiology: VivoPreparedPhysiologyStep?
        var lastExposure: [Float] = [], lastFeedback: [Float] = []
        var lastResidual = Double.greatestFiniteMagnitude
        var finalIterations: UInt32 = 0
        var negotiations: UInt32 = 0
        var molecularPending = false, physiologyPending = false, releasing = false
        do {
            negotiation: for negotiation in 0..<bridge.policy.maximumStepNegotiations {
                negotiations = negotiation + 1
                var guess = previousFeedbackValues
                var haveGuess = hasAcceptedFeedback || guess.isEmpty
                var iterationController = try VivoCouplingIterationController(relaxation: bridge.policy.relaxation)
                for iteration in 0..<bridge.policy.maximumFixedPointIterations {
                    finalIterations = iteration + 1
                    try Task.checkCancellation()
                    let feedbackUpdates: [VivoPhysiologyStateUpdate] = haveGuess ? zip(bridge.molecularToPhysiology, guess).map {
                        .init(pairIndex: $0.0.pairIndex, environmentIndex: $0.0.environmentIndex, mode: $0.0.physiologyMode, value: $0.1)
                    } : []
                    let physiologyPublications = bridge.physiologyToMolecular.map {
                        VivoPhysiologyPublicationRequest(pairIndex: $0.pairIndex, environmentIndex: $0.environmentIndex)
                    } + request.physiologyPublications
                    let p = try await physiology.prepareStep(.init(timeStepSeconds: dt, preUpdates: request.physiologyPreUpdates,
                                                                    postUpdates: request.physiologyPostUpdates + feedbackUpdates,
                                                                    publications: physiologyPublications, permitAdaptiveReduction: false), transactionID: id)
                    lastPhysiology = p; physiologyPending = p.canCommit
                    // Compare actual horizons exactly, not a tolerance scaled by 1s.
                    if !p.canCommit || p.candidateTimeStep != dt {
                        if physiologyPending { try await physiology.discardPreparedStep(transactionID: id); physiologyPending = false }
                        guard mayReduce else { break negotiation }
                        let candidate = p.canCommit ? p.candidateTimeStep : dt * 0.5
                        let smaller = try canonicalStep(min(candidate, dt), minimum: minimum)
                        guard smaller < dt else { break negotiation }
                        dt = smaller
                        continue negotiation
                    }
                    guard p.publications.count == physiologyPublications.count else { throw VivoRuntimeError.commandFailed("missing physiology publications") }
                    let exposureCount = bridge.physiologyToMolecular.count
                    let exposure = try zip(bridge.physiologyToMolecular, p.publications.prefix(exposureCount)).map {
                        try $0.0.transfer.apply($0.1, subject: $0.0.identifier)
                    }
                    lastExposure = exposure
                    let coupling = zip(bridge.physiologyToMolecular, exposure).map {
                        VivoCouplingUpdate(speciesIndex: $0.0.molecularSpeciesIndex, laneIndex: $0.0.molecularLaneIndex,
                                           mode: $0.0.molecularMode, value: $0.1)
                    } + request.molecularCoupling
                    let publications = bridge.molecularToPhysiology.map {
                        VivoPublicationRequest(speciesIndex: $0.molecularSpeciesIndex, laneIndex: $0.molecularLaneIndex)
                    } + request.molecularPublications
                    let m = try await molecular.prepareStep(.init(timeStep: Float(dt), coupling: coupling, publications: publications,
                                                                  permitAdaptiveReduction: false), transactionID: id)
                    lastMolecular = m; molecularPending = m.canCommit
                    if !m.canCommit || Double(m.candidateTimeStep) != dt {
                        if molecularPending { try await molecular.discardPreparedStep(transactionID: id); molecularPending = false }
                        try await physiology.discardPreparedStep(transactionID: id); physiologyPending = false
                        if m.disposition == .permanentShutdown || m.disposition == .reversibleShutdown { break negotiation }
                        guard mayReduce else { break negotiation }
                        dt = try canonicalStep(dt * 0.5, minimum: minimum)
                        continue negotiation
                    }
                    guard m.publications.count == publications.count else { throw VivoRuntimeError.commandFailed("missing molecular publications") }
                    let feedbackCount = bridge.molecularToPhysiology.count
                    let feedback = try zip(bridge.molecularToPhysiology, m.publications.prefix(feedbackCount)).map {
                        try $0.0.transfer.apply($0.1, subject: $0.0.identifier)
                    }
                    lastFeedback = feedback
                    let residual = haveGuess ? try VivoCouplingIterationController.residual(candidate: feedback, reference: guess) : Double.greatestFiniteMagnitude
                    lastResidual = residual
                    let converged = feedback.isEmpty || (haveGuess && residual <= bridge.policy.convergenceTolerance)
                    let last = iteration + 1 == bridge.policy.maximumFixedPointIterations
                    if converged || (last && !bridge.policy.requireConvergence) {
                        try Task.checkCancellation()
                        releasing = true
                        // Private participants cannot be mutated by another client;
                        // both candidates are validated before any release. Do not
                        // observe cancellation between these non-I/O actor commits.
                        let mc = try await molecular.commitPreparedStep(transactionID: id)
                        molecularPending = false
                        let pc = try await physiology.commitPreparedStep(transactionID: id)
                        physiologyPending = false
                        guard pc.timeAfter == before + dt else { throw VivoRuntimeError.commandFailed("release clock invariant failed") }
                        acceptedTime = before + dt
                        let index = acceptedStepIndex
                        acceptedStepIndex += 1
                        previousFeedbackValues = feedback; hasAcceptedFeedback = true
                        releasing = false
                        let disposition: VivoMolecularPhysiologyDisposition = !converged ? .committedWithoutConvergence :
                            (dt == requested ? .committed : .committedWithReducedStep)
                        return .init(certificate: .init(transactionID: id, disposition: disposition, bridgeFingerprint: bridge.fingerprint,
                                                        programFingerprint: programPack.header.contentFingerprint, physiologyFingerprint: physiologyModel.fingerprint,
                                                        executableFingerprint: executableFingerprint, stepIndex: index, timeBeforeSeconds: before,
                                                        timeAfterSeconds: acceptedTime, requestedTimeStepSeconds: requested, acceptedTimeStepSeconds: dt,
                                                        stepNegotiationCount: negotiations, fixedPointIterationCount: finalIterations, finalResidual: residual,
                                                        converged: converged, molecularStatus: m.status, physiologyStatus: p.status,
                                                        message: converged ? "joint candidate released" : "nonconverged candidate released under explicit policy"),
                                     molecularCertificate: mc, physiologyCertificate: pc, molecularEvents: m.events,
                                     bridgeExposureValues: exposure, bridgeFeedbackValues: feedback,
                                     molecularPublications: Array(m.publications.dropFirst(feedbackCount)),
                                     physiologyPublications: Array(p.publications.dropFirst(exposureCount)))
                    }
                    try await molecular.discardPreparedStep(transactionID: id); molecularPending = false
                    try await physiology.discardPreparedStep(transactionID: id); physiologyPending = false
                    guess = haveGuess ? try iterationController.next(current: guess, proposed: feedback) : feedback
                    haveGuess = true
                }
                break // Nonconvergence is not evidence that a smaller random horizon is unbiased.
            }
        } catch {
            if releasing {
                failure = "joint release failed; coordinator quarantined, no rollback claimed: \(error)"
            } else {
                do {
                    if molecularPending { try await molecular.discardPreparedStep(transactionID: id) }
                    if physiologyPending { try await physiology.discardPreparedStep(transactionID: id) }
                } catch { failure = "candidate cleanup failed: \(error)" }
            }
            throw error
        }
        let disposition: VivoMolecularPhysiologyDisposition = lastMolecular?.disposition == .permanentShutdown ? .permanentShutdown :
            (lastMolecular?.disposition == .reversibleShutdown ? .reversibleShutdown : .rejected)
        return .init(certificate: .init(transactionID: id, disposition: disposition, bridgeFingerprint: bridge.fingerprint,
                                        programFingerprint: programPack.header.contentFingerprint, physiologyFingerprint: physiologyModel.fingerprint,
                                        executableFingerprint: executableFingerprint, stepIndex: acceptedStepIndex,
                                        timeBeforeSeconds: before, timeAfterSeconds: before, requestedTimeStepSeconds: requested,
                                        acceptedTimeStepSeconds: nil, stepNegotiationCount: negotiations, fixedPointIterationCount: finalIterations,
                                        finalResidual: lastResidual, converged: false, molecularStatus: lastMolecular?.status,
                                        physiologyStatus: lastPhysiology?.status, message: "joint candidate not released"),
                     molecularCertificate: nil, physiologyCertificate: nil, molecularEvents: [],
                     bridgeExposureValues: lastExposure, bridgeFeedbackValues: lastFeedback,
                     molecularPublications: [], physiologyPublications: [])
    }
    private func canonicalStep(_ value: Double, minimum: Double) throws -> Double {
        guard value.isFinite, value > 0, value <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoRuntimeError.invalidConfiguration("common timestep is not FP32 representable")
        }
        var rounded = Float(value)
        if Double(rounded) > value { rounded = rounded.nextDown }
        guard rounded > 0, Double(rounded) >= minimum else {
            throw VivoRuntimeError.invalidConfiguration("no common FP32 timestep fits before this boundary; reduce the minimum timestep")
        }
        return Double(rounded)
    }
}
