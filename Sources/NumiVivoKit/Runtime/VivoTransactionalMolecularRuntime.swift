import Foundation
@preconcurrency import Metal
import NumiVivoShaders

/// ProgramPack transaction owner. Every operation that suspends while using GPU
/// resources holds an explicit reservation; committed readers never see a candidate.
public actor VivoTransactionalMolecularRuntime {
    public enum Lifecycle: Sendable, Codable, Equatable {
        case ready
        case reversiblyStopped(reason: String)
        case permanentlyStopped(reason: String)
        case failed(reason: String)
    }
    private struct AttemptResult {
        let status: VivoRuntimeStatus
        let events: [VivoEvent]
        let publications: [Float]
    }
    public nonisolated let pack: VivoProgramPack
    public nonisolated let configuration: VivoRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities
    private let commandQueue: MTLCommandQueue
    private var arena: VivoMetalArena
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let species: [VivoProgramPack.SpeciesMetadata]
    private let containsCountValuedSpecies: Bool
    private var lifecycleState: Lifecycle = .ready
    private var absoluteTimeSeconds: Double = 0
    private var nextStepIndex: UInt32 = 0
    private var pending: VivoPreparedMolecularStep?
    private var inFlight = false

    public static func make(pack: VivoProgramPack, configuration: VivoRuntimeConfiguration,
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoTransactionalMolecularRuntime {
        try VivoProgramExecutionContract.validate(pack: pack, configuration: configuration)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard let queue = device.makeCommandQueue() else { throw VivoRuntimeError.commandQueueUnavailable }
        queue.label = "NumiVivo.TransactionalMolecularQueue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in [NumiVivoKernel.clearStatus, .prepareTransaction, .applyCouplingUpdates,
                       .f1HeunPredict, .f1HeunCorrect, .f2SampleReactions, .f2ApplyReactions,
                       .f3Transport, .executeRules, .evaluateMonitors, .validateShadow, .publish] {
            pipelines[kernel] = try await catalog.pipeline(kernel)
        }
        let arena = try VivoMetalArena(device: device, commandQueue: queue, pack: pack, configuration: configuration)
        return try VivoTransactionalMolecularRuntime(pack: pack, configuration: configuration,
                                                      queue: queue, arena: arena, pipelines: pipelines)
    }
    private init(pack: VivoProgramPack, configuration: VivoRuntimeConfiguration,
                 queue: MTLCommandQueue, arena: VivoMetalArena,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline]) throws {
        self.pack = pack
        self.configuration = configuration
        commandQueue = queue
        self.arena = arena
        self.pipelines = pipelines
        species = try pack.speciesMetadata()
        containsCountValuedSpecies = species.contains { $0.isCountValued }
        capabilities = arena.capabilities
    }
    public func lifecycle() -> Lifecycle { lifecycleState }
    public func time() -> Float { Float(absoluteTimeSeconds) }
    public func timeSeconds() -> Double { absoluteTimeSeconds }
    public func stepIndex() -> UInt32 { nextStepIndex }
    public func hasPendingTransaction() -> Bool { inFlight || pending != nil }

    public func resume() throws {
        try requireUnreserved()
        switch lifecycleState {
        case .reversiblyStopped: lifecycleState = .ready
        case .ready: return
        case .permanentlyStopped(let reason), .failed(let reason): throw VivoRuntimeError.runtimeStopped(reason)
        }
    }
    public func requireAcceptedBoundary() throws { try requireUnreserved() }
    public func stopReversibly(reason: String) throws {
        try requireReadyAndUnreserved()
        guard !reason.isEmpty else { throw VivoRuntimeError.invalidConfiguration("empty reversible shutdown reason") }
        lifecycleState = .reversiblyStopped(reason: reason)
    }
    public func stopPermanently(reason: String) {
        // Already-submitted GPU work may finish, but cannot become publishable.
        pending = nil
        lifecycleState = .permanentlyStopped(reason: reason)
    }
    public func setTransport(_ values: [VivoSpeciesTransportABI]) throws {
        try requireReadyAndUnreserved()
        try validateTransport(values)
        try arena.replaceTransport(values, commandQueue: commandQueue)
    }
    private func validateTransport(_ values: [VivoSpeciesTransportABI]) throws {
        guard values.count == species.count else { throw VivoRuntimeError.invalidConfiguration("transport/species count mismatch") }
        for (index, value) in values.enumerated() {
            guard value.diffusion.isFinite, value.diffusion >= 0,
                  value.membranePermeability == 0,
                  value.decayRate.isFinite, value.decayRate >= 0 else {
                throw VivoRuntimeError.invalidConfiguration("diffusion/decay must be finite and nonnegative; nonzero membrane permeability requires a membrane-interface backend")
            }
            if species[index].isExternallyOwned, value.diffusion != 0 || value.decayRate != 0 {
                throw VivoRuntimeError.invalidConfiguration("transport cannot mutate externally owned species")
            }
        }
    }
    public func setVelocity(_ values: [SIMD4<Float>]) throws {
        try requireReadyAndUnreserved()
        guard values.allSatisfy({ $0.x.isFinite && $0.y.isFinite && $0.z.isFinite && $0.w.isFinite }) else {
            throw VivoRuntimeError.invalidConfiguration("velocity field contains non-finite values")
        }
        try arena.replaceVelocity(values, commandQueue: commandQueue)
    }
    public func setVolumeFractions(_ values: [Float]) throws {
        try requireReadyAndUnreserved()
        try arena.replaceVolumeFraction(values, commandQueue: commandQueue)
    }

    public func prepareStep(_ request: VivoStepRequest = .init(), transactionID: UUID = UUID()) async throws -> VivoPreparedMolecularStep {
        try requireReadyAndUnreserved()
        try Task.checkCancellation()
        let requested = request.timeStep ?? configuration.timeStep
        guard requested.isFinite, requested >= configuration.minimumTimeStep,
              requested <= configuration.maximumTimeStep, nextStepIndex < UInt32.max,
              absoluteTimeSeconds + Double(requested) > absoluteTimeSeconds,
              Float(absoluteTimeSeconds + Double(requested)).isFinite else {
            throw VivoRuntimeError.invalidConfiguration("molecular time step, clock, or step index cannot advance")
        }
        try validate(request)
        inFlight = true
        defer { inFlight = false }
        try arena.write(couplingUpdates: request.coupling.map(\.abi))
        try arena.write(publicationRequests: request.publications.enumerated().map { index, value in
            .init(speciesIndex: value.speciesIndex, laneIndex: value.laneIndex, outputIndex: UInt32(index), flags: value.flags)
        })
        var dt = requested
        for attempt in 0..<configuration.maximumSubsteps {
            guard Float(absoluteTimeSeconds + Double(dt)) > Float(absoluteTimeSeconds) else {
                throw VivoRuntimeError.invalidConfiguration("time step is below the ProgramPack FP32 absolute-clock resolution")
            }
            let output = try await executeAttempt(dt: dt, attempt: attempt,
                                                  couplingCount: request.coupling.count,
                                                  publicationCount: request.publications.count)
            try ensureReadyForMutation()
            try Task.checkCancellation()
            let status = output.status
            let disposition: VivoPreparedMolecularDisposition
            if status.flags.contains(.permanentShutdown) || status.requestedResponse == 5 {
                lifecycleState = .permanentlyStopped(reason: shutdownReason(status))
                disposition = .permanentShutdown
            } else if status.flags.contains(.reversibleShutdown) || status.requestedResponse == 4 {
                lifecycleState = .reversiblyStopped(reason: shutdownReason(status))
                disposition = .reversibleShutdown
            } else if !status.blocksCommit && !status.flags.contains(.invalidRate) &&
                        !status.flags.contains(.eventOverflow) && !status.flags.contains(.randomSaturated) &&
                        !status.flags.contains(.substep) && status.requestedResponse != 1 && status.requestedResponse != 3 {
                disposition = .prepared
            } else {
                let reduced = dt * 0.5
                if request.permitAdaptiveReduction, attempt + 1 < configuration.maximumSubsteps,
                   reduced >= configuration.minimumTimeStep, reduced < dt {
                    dt = reduced
                    continue
                }
                disposition = .requiresSmallerStep
            }
            let value = VivoPreparedMolecularStep(transactionID: transactionID, disposition: disposition,
                                                   sourceFingerprint: pack.header.sourceFingerprint,
                                                   programFingerprint: pack.header.contentFingerprint,
                                                   fidelity: configuration.fidelity, stepIndex: nextStepIndex,
                                                   timeBefore: Float(absoluteTimeSeconds), requestedTimeStep: requested,
                                                   candidateTimeStep: dt, attemptCount: attempt + 1, status: status,
                                                   events: output.events,
                                                   publications: disposition == .prepared ? output.publications : [])
            if value.canCommit { pending = value }
            return value
        }
        throw VivoRuntimeError.commandFailed("molecular attempt budget exhausted")
    }
    public func commitPreparedStep(transactionID: UUID) throws -> VivoMolecularTransactionCertificate {
        try ensureReadyForMutation()
        guard !inFlight, let candidate = pending, candidate.transactionID == transactionID, candidate.canCommit else {
            throw VivoRuntimeError.invalidConfiguration("molecular transaction does not match an eligible candidate")
        }
        arena.commit()
        absoluteTimeSeconds += Double(candidate.candidateTimeStep)
        nextStepIndex += 1
        pending = nil
        return certificate(candidate, disposition: candidate.candidateTimeStep == candidate.requestedTimeStep ? .committed : .committedWithReducedStep,
                           accepted: candidate.candidateTimeStep, timeAfter: Float(absoluteTimeSeconds))
    }
    public func discardPreparedStep(transactionID: UUID) throws {
        guard !inFlight, pending?.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration("molecular transaction does not match the prepared candidate")
        }
        pending = nil
    }
    public func step(_ request: VivoStepRequest = .init()) async throws -> VivoMolecularTransactionResult {
        let candidate = try await prepareStep(request)
        if candidate.canCommit {
            return .init(certificate: try commitPreparedStep(transactionID: candidate.transactionID),
                          events: candidate.events, publications: candidate.publications)
        }
        let disposition: VivoStepDisposition = candidate.disposition == .permanentShutdown ? .permanentShutdown :
            (candidate.disposition == .reversibleShutdown ? .reversibleShutdown : .rejected)
        return .init(certificate: certificate(candidate, disposition: disposition, accepted: nil, timeAfter: candidate.timeBefore),
                      events: candidate.events, publications: [])
    }
    public func snapshot() async throws -> VivoStateSnapshot {
        try requireUnreserved()
        inFlight = true
        defer { inFlight = false }
        guard let command = commandQueue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        blit.copy(from: arena.currentState, sourceOffset: 0, to: arena.stateReadbackBuffer, destinationOffset: 0,
                  size: arena.capacities.stateElements * 4)
        blit.endEncoding()
        try await complete(command)
        let values = Array(UnsafeBufferPointer(start: arena.stateReadbackBuffer.contents().assumingMemoryBound(to: Float.self), count: arena.capacities.stateElements))
        return .init(sourceFingerprint: pack.header.sourceFingerprint, programFingerprint: pack.header.contentFingerprint,
                      stepIndex: nextStepIndex, absoluteTime: Float(absoluteTimeSeconds),
                      speciesCount: pack.runtimeContract.speciesCount, laneCount: configuration.laneCount, values: values)
    }

    public func checkpoint() throws -> VivoMolecularCheckpoint {
        try requireReadyAndUnreserved()
        let data = try arena.captureCheckpoint(commandQueue: commandQueue)
        let parameterScalars = Int(pack.runtimeContract.parameterCount) * Int(configuration.environmentCount)
        let value = VivoMolecularCheckpoint(programFingerprint: pack.header.contentFingerprint,
                                            sourceProgramFingerprint: pack.header.sourceFingerprint,
                                            fidelity: configuration.fidelity, seed: configuration.seed,
                                            stepIndex: nextStepIndex, absoluteTimeSeconds: absoluteTimeSeconds,
                                            speciesCount: pack.runtimeContract.speciesCount, laneCount: configuration.laneCount,
                                            parameterCount: pack.runtimeContract.parameterCount,
                                            parameterEnvironmentCount: configuration.environmentCount,
                                            temporalStateCount: pack.runtimeContract.temporalStateCount,
                                            stateFP32LE: VivoLittleEndianFP32.encode(data.state),
                                            parametersFP32LE: VivoLittleEndianFP32.encode(Array(data.parameters.prefix(parameterScalars))),
                                            temporalStateFP32LE: VivoLittleEndianFP32.encode(data.temporalState),
                                            transportRecordLE: VivoTransportRecordLE.encode(data.transport),
                                            velocityFP32LE: VivoLittleEndianFP32.encode(data.velocity),
                                            volumeFractionFP32LE: VivoLittleEndianFP32.encode(data.volumeFractions))
        try value.validate()
        return value
    }
    public func resumeCheckpoint() throws -> VivoMolecularResumeCheckpoint {
        try VivoMolecularResumeCheckpoint(configuration: configuration, state: checkpoint())
    }
    public func restore(_ checkpoint: VivoMolecularCheckpoint) throws {
        guard configuration.spatialGrid == nil else {
            throw VivoRuntimeError.invalidConfiguration("raw v2 checkpoints omit spatial layout; restore the configuration-bound resume checkpoint")
        }
        try restoreState(checkpoint)
    }
    public func restore(_ checkpoint: VivoMolecularResumeCheckpoint) throws {
        try checkpoint.validate()
        guard checkpoint.configuration == configuration else {
            throw VivoRuntimeError.invalidConfiguration("resume configuration differs, including spatial grid, seed or numerical settings")
        }
        try restoreState(checkpoint.state)
    }
    private func restoreState(_ checkpoint: VivoMolecularCheckpoint) throws {
        try requireReadyAndUnreserved()
        try checkpoint.validate()
        guard checkpoint.programFingerprint == pack.header.contentFingerprint,
              checkpoint.sourceProgramFingerprint == pack.header.sourceFingerprint,
              checkpoint.fidelity == configuration.fidelity, checkpoint.seed == configuration.seed,
              checkpoint.speciesCount == pack.runtimeContract.speciesCount,
              checkpoint.laneCount == configuration.laneCount,
              checkpoint.parameterCount == pack.runtimeContract.parameterCount,
              checkpoint.parameterEnvironmentCount == configuration.environmentCount,
              checkpoint.temporalStateCount == pack.runtimeContract.temporalStateCount,
              checkpoint.absoluteTimeSeconds <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoRuntimeError.invalidConfiguration("checkpoint identity, layout, seed, fidelity or clock mismatch")
        }
        let values = try VivoLittleEndianFP32.decode(checkpoint.stateFP32LE)
        for (index, value) in values.enumerated() {
            let metadata = species[index / Int(configuration.laneCount)]
            guard value >= metadata.minimum, value <= metadata.maximum else {
                throw VivoRuntimeError.invalidConfiguration("checkpoint violates species bounds")
            }
            if metadata.isCountValued, (value < 0 || value.rounded() != value || value > 16_777_216) {
                throw VivoRuntimeError.invalidConfiguration("count state is not an exactly representable FP32 integer")
            }
        }
        let parameters = try VivoLittleEndianFP32.decode(checkpoint.parametersFP32LE)
        let metadata = try pack.parameterMetadata()
        for (index, value) in parameters.enumerated() {
            let bounds = metadata[index / Int(configuration.environmentCount)]
            guard Double(value) >= bounds.minimum, Double(value) <= bounds.maximum else {
                throw VivoRuntimeError.invalidConfiguration("checkpoint violates parameter bounds")
            }
        }
        let transport = try VivoTransportRecordLE.decode(checkpoint.transportRecordLE)
        try validateTransport(transport)
        let decoded = VivoMolecularArenaCheckpointData(state: values,
                                                       parameters: parameters.isEmpty ? [0] : parameters,
                                                       temporalState: try VivoLittleEndianFP32.decode(checkpoint.temporalStateFP32LE),
                                                       transport: transport,
                                                       velocity: try VivoLittleEndianFP32.decode(checkpoint.velocityFP32LE),
                                                       volumeFractions: try VivoLittleEndianFP32.decode(checkpoint.volumeFractionFP32LE))
        let replacement = try VivoMetalArena(device: arena.device, commandQueue: commandQueue,
                                              pack: pack, configuration: configuration)
        try replacement.restoreCheckpoint(decoded, commandQueue: commandQueue)
        arena = replacement
        absoluteTimeSeconds = checkpoint.absoluteTimeSeconds
        nextStepIndex = checkpoint.stepIndex
    }

    private func executeAttempt(dt: Float, attempt: UInt32, couplingCount: Int, publicationCount: Int) async throws -> AttemptResult {
        guard let command = commandQueue.makeCommandBuffer() else { throw VivoRuntimeError.commandQueueUnavailable }
        command.label = "NumiVivo.Molecular.\(nextStepIndex).\(attempt)"
        let uniforms = try VivoRuntimeCommandABI(pack: pack, configuration: configuration, stepIndex: nextStepIndex,
                                                 substepIndex: attempt, dt: dt, absoluteTime: Float(absoluteTimeSeconds + Double(dt)))
        try encode(.clearStatus, count: 1, command: command, buffers: [arena.statusBuffer])
        try encode(.prepareTransaction, count: max(arena.capacities.stateElements, arena.capacities.temporalElements), command: command,
                   buffers: [arena.currentState, arena.baseState, arena.candidateState, arena.temporalCurrent, arena.temporalCandidate, arena.statusBuffer], uniforms: uniforms)
        if couplingCount > 0 {
            try encode(.applyCouplingUpdates, count: couplingCount, command: command,
                       buffers: [arena.baseState, arena.couplingBuffer, arena.statusBuffer], uniforms: uniforms, extraCount: UInt32(couplingCount))
        }
        switch configuration.fidelity {
        case .logic: try copy(command, from: arena.baseState, to: arena.candidateState)
        case .deterministic:
            try deterministic(command, uniforms: uniforms, source: arena.baseState, predictor: arena.stageState, destination: arena.candidateState)
        case .stochastic:
            try stochastic(command, uniforms: uniforms, source: arena.baseState, destination: arena.candidateState)
        case .spatial, .tissue:
            let doubled = attempt.multipliedReportingOverflow(by: 2)
            guard !doubled.overflow, doubled.partialValue < UInt32.max else { throw VivoRuntimeError.invalidConfiguration("split-step random namespace overflow") }
            var first = uniforms
            first.dt *= 0.5
            first.absoluteTime = Float(absoluteTimeSeconds + Double(first.dt))
            first.substepIndex = doubled.partialValue
            if containsCountValuedSpecies {
                try stochastic(command, uniforms: first, source: arena.baseState, destination: arena.candidateState)
            } else {
                try deterministic(command, uniforms: first, source: arena.baseState, predictor: arena.stageState, destination: arena.candidateState)
            }
            try encode(.f3Transport, count: arena.capacities.stateElements, command: command,
                       buffers: [arena.candidateState, arena.stageState, arena.transport, arena.velocity, arena.volumeFraction, arena.statusBuffer], uniforms: uniforms)
            var second = first
            second.absoluteTime = uniforms.absoluteTime
            second.substepIndex += 1
            if containsCountValuedSpecies {
                try stochastic(command, uniforms: second, source: arena.stageState, destination: arena.candidateState)
            } else {
                try deterministic(command, uniforms: second, source: arena.stageState, predictor: arena.baseState, destination: arena.candidateState)
            }
        }
        let controls = [arena.candidateState, arena.temporalCandidate, arena.parameterBuffer, arena.programBuffer, arena.eventBuffer, arena.statusBuffer]
        if pack.runtimeContract.ruleCount > 0 {
            try encode(.executeRules, count: Int(configuration.laneCount), command: command, buffers: controls, uniforms: uniforms)
        }
        if pack.runtimeContract.monitorCount > 0 {
            try encode(.evaluateMonitors, count: Int(configuration.laneCount), command: command, buffers: controls, uniforms: uniforms)
        }
        try encode(.validateShadow, count: arena.capacities.stateElements, command: command,
                   buffers: [arena.candidateState, arena.programBuffer, arena.statusBuffer], uniforms: uniforms)
        if publicationCount > 0 {
            try encode(.publish, count: publicationCount, command: command,
                       buffers: [arena.candidateState, arena.publicationRequestBuffer, arena.publicationOutputBuffer, arena.statusBuffer],
                       uniforms: uniforms, extraCount: UInt32(publicationCount))
        }
        try await complete(command)
        let status = arena.status()
        return .init(status: status, events: arena.events(status: status), publications: try arena.publicationValues(count: publicationCount))
    }
    private func deterministic(_ command: MTLCommandBuffer, uniforms: VivoRuntimeCommandABI,
                               source: MTLBuffer, predictor: MTLBuffer, destination: MTLBuffer) throws {
        try encode(.f1HeunPredict, count: arena.capacities.stateElements, command: command,
                   buffers: [source, predictor, arena.derivativeK1, arena.parameterBuffer, arena.programBuffer, arena.statusBuffer], uniforms: uniforms)
        try encode(.f1HeunCorrect, count: arena.capacities.stateElements, command: command,
                   buffers: [source, predictor, arena.derivativeK1, destination, arena.parameterBuffer, arena.programBuffer, arena.statusBuffer], uniforms: uniforms)
    }
    private func stochastic(_ command: MTLCommandBuffer, uniforms: VivoRuntimeCommandABI,
                            source: MTLBuffer, destination: MTLBuffer) throws {
        try encode(.f2SampleReactions, count: try work(pack.runtimeContract.reactionCount), command: command,
                   buffers: [source, arena.parameterBuffer, arena.reactionEvents, arena.programBuffer, arena.statusBuffer], uniforms: uniforms)
        try encode(.f2ApplyReactions, count: arena.capacities.stateElements, command: command,
                   buffers: [source, destination, arena.reactionEvents, arena.programBuffer, arena.statusBuffer], uniforms: uniforms)
    }
    private func encode(_ kernel: NumiVivoKernel, count: Int, command: MTLCommandBuffer,
                        buffers: [MTLBuffer], uniforms: VivoRuntimeCommandABI? = nil, extraCount: UInt32? = nil) throws {
        guard count > 0 else { return }
        guard let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.pipelineEncodingFailed(kernel.rawValue)
        }
        defer { encoder.endEncoding() }
        encoder.label = kernel.rawValue
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        if var uniforms {
            encoder.setBytes(&uniforms, length: MemoryLayout<VivoRuntimeCommandABI>.stride, index: buffers.count)
        }
        if var extraCount { encoder.setBytes(&extraCount, length: 4, index: buffers.count + 1) }
        if kernel == .f3Transport {
            // Species ownership is consulted by the finite-volume pass as well
            // as reactions. Externally-owned fields cannot be advected locally.
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 8)
        }
        encoder.dispatchThreads(pipeline.gridSize(for: count), threadsPerThreadgroup: pipeline.threadgroupSize(for: count, preferred: capabilities.recommendedThreadsPerThreadgroup))
    }
    private func copy(_ command: MTLCommandBuffer, from source: MTLBuffer, to destination: MTLBuffer) throws {
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoRuntimeError.commandQueueUnavailable }
        blit.copy(from: source, sourceOffset: 0, to: destination, destinationOffset: 0, size: arena.capacities.stateElements * 4)
        blit.endEncoding()
    }
    private func complete(_ command: MTLCommandBuffer) async throws {
        do {
            try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
                command.addCompletedHandler { completed in
                    if let error = completed.error { continuation.resume(throwing: VivoRuntimeError.commandFailed(String(describing: error))) }
                    else { continuation.resume(returning: ()) }
                }
                command.commit()
            }
        } catch {
            if case .permanentlyStopped = lifecycleState {} else { lifecycleState = .failed(reason: String(describing: error)) }
            throw error
        }
    }
    private func validate(_ request: VivoStepRequest) throws {
        guard request.coupling.count <= arena.capacities.couplingUpdates,
              request.publications.count <= arena.capacities.publicationRequests else { throw VivoRuntimeError.invalidConfiguration("transaction boundary exceeds capacity") }
        var destinations = Set<UInt64>()
        for update in request.coupling {
            guard Int(update.speciesIndex) < species.count, update.laneIndex < configuration.laneCount, update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("invalid coupling update")
            }
            let metadata = species[Int(update.speciesIndex)]
            guard metadata.isExternallyOwned || metadata.isInput,
                  destinations.insert(UInt64(update.speciesIndex) << 32 | UInt64(update.laneIndex)).inserted else {
                throw VivoRuntimeError.invalidConfiguration("internally owned or multiply-written coupling destination")
            }
        }
        for value in request.publications {
            guard Int(value.speciesIndex) < species.count, value.laneIndex < configuration.laneCount else {
                throw VivoRuntimeError.invalidConfiguration("invalid publication index")
            }
        }
    }
    private func requireUnreserved() throws {
        guard !inFlight, pending == nil else { throw VivoRuntimeError.invalidConfiguration("molecular operation or prepared transaction is already active") }
    }
    private func requireReadyAndUnreserved() throws { try ensureReadyForMutation(); try requireUnreserved() }
    private func ensureReadyForMutation() throws {
        switch lifecycleState {
        case .ready: return
        case .reversiblyStopped(let reason), .permanentlyStopped(let reason), .failed(let reason): throw VivoRuntimeError.runtimeStopped(reason)
        }
    }
    private func work(_ count: UInt32) throws -> Int {
        let value = UInt64(count) * UInt64(configuration.laneCount)
        guard value <= UInt64(UInt32.max), value <= UInt64(Int.max) else { throw VivoRuntimeError.invalidConfiguration("dispatch grid exceeds UInt32") }
        return Int(value)
    }
    private func shutdownReason(_ status: VivoRuntimeStatus) -> String {
        "monitor \(status.firstMonitor.map(String.init) ?? "unknown") requested response \(status.requestedResponse)"
    }
    private func certificate(_ value: VivoPreparedMolecularStep, disposition: VivoStepDisposition,
                             accepted: Float?, timeAfter: Float) -> VivoMolecularTransactionCertificate {
        .init(transactionID: value.transactionID, disposition: disposition, sourceFingerprint: value.sourceFingerprint,
              programFingerprint: value.programFingerprint, fidelity: value.fidelity,
              deviceName: capabilities.deviceName, deviceRegistryID: capabilities.registryID,
              stepIndex: value.stepIndex, timeBefore: value.timeBefore, timeAfter: timeAfter,
              requestedTimeStep: value.requestedTimeStep, acceptedTimeStep: accepted,
              attemptCount: value.attemptCount, status: value.status)
    }
}
