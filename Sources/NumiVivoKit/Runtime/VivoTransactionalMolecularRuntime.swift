import Foundation
import Metal
import NumiVivoShaders

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

    private struct PendingStep {
        let prepared: VivoPreparedMolecularStep
        let requestedTimeStep: Float
    }

    public nonisolated let pack: VivoProgramPack
    public nonisolated let configuration: VivoRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities

    private let commandQueue: MTLCommandQueue
    private let arena: VivoMetalArena
    private let pipelines: VivoTransactionalRuntimePipelines
    private let containsCountValuedSpecies: Bool

    private var lifecycleState: Lifecycle = .ready
    private var absoluteTime: Float = 0
    private var nextStepIndex: UInt32 = 0
    private var pending: PendingStep?

    public static func make(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoTransactionalMolecularRuntime {
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard let queue = device.makeCommandQueue() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.TransactionalMolecularQueue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let pipelines = try await VivoTransactionalRuntimePipelines.load(from: catalog)
        let arena = try VivoMetalArena(
            device: device,
            commandQueue: queue,
            pack: pack,
            configuration: configuration
        )
        let countValued = try pack.speciesMetadata().contains(where: { $0.isCountValued })
        return VivoTransactionalMolecularRuntime(
            pack: pack,
            configuration: configuration,
            commandQueue: queue,
            arena: arena,
            pipelines: pipelines,
            containsCountValuedSpecies: countValued
        )
    }

    private init(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        commandQueue: MTLCommandQueue,
        arena: VivoMetalArena,
        pipelines: VivoTransactionalRuntimePipelines,
        containsCountValuedSpecies: Bool
    ) {
        self.pack = pack
        self.configuration = configuration
        self.commandQueue = commandQueue
        self.arena = arena
        self.pipelines = pipelines
        self.capabilities = arena.capabilities
        self.containsCountValuedSpecies = containsCountValuedSpecies
    }

    public func lifecycle() -> Lifecycle { lifecycleState }
    public func time() -> Float { absoluteTime }
    public func stepIndex() -> UInt32 { nextStepIndex }
    public func hasPendingTransaction() -> Bool { pending != nil }

    public func resume() throws {
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "cannot resume while a molecular transaction is pending"
            )
        }
        switch lifecycleState {
        case .reversiblyStopped:
            lifecycleState = .ready
        case .ready:
            return
        case .permanentlyStopped(let reason):
            throw VivoRuntimeError.runtimeStopped("permanent shutdown: \(reason)")
        case .failed(let reason):
            throw VivoRuntimeError.runtimeStopped("runtime failure: \(reason)")
        }
    }

    public func stopPermanently(reason: String) {
        pending = nil
        lifecycleState = .permanentlyStopped(reason: reason)
    }

    public func setTransport(_ values: [VivoSpeciesTransportABI]) throws {
        try ensureReadyForMutation()
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "cannot change transport while a molecular transaction is pending"
            )
        }
        try arena.replaceTransport(values, commandQueue: commandQueue)
    }

    public func setVelocity(_ values: [SIMD4<Float>]) throws {
        try ensureReadyForMutation()
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "cannot change velocity while a molecular transaction is pending"
            )
        }
        try arena.replaceVelocity(values, commandQueue: commandQueue)
    }

    public func setVolumeFractions(_ values: [Float]) throws {
        try ensureReadyForMutation()
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "cannot change volume fractions while a molecular transaction is pending"
            )
        }
        try arena.replaceVolumeFraction(values, commandQueue: commandQueue)
    }

    public func prepareStep(
        _ request: VivoStepRequest = .init(),
        transactionID: UUID = UUID()
    ) async throws -> VivoPreparedMolecularStep {
        try ensureReadyForMutation()
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular runtime already has a prepared transaction"
            )
        }

        let requestedStep = request.timeStep ?? configuration.timeStep
        guard requestedStep.isFinite,
              requestedStep >= configuration.minimumTimeStep,
              requestedStep <= configuration.maximumTimeStep else {
            throw VivoRuntimeError.invalidConfiguration(
                "requested molecular time step is outside configured finite bounds"
            )
        }
        try validate(request: request)

        let coupling = request.coupling.map { $0.abi }
        let publications = request.publications.enumerated().map { index, request in
            VivoPublicationRequestABI(
                speciesIndex: request.speciesIndex,
                laneIndex: request.laneIndex,
                outputIndex: UInt32(index),
                flags: request.flags
            )
        }
        try arena.write(couplingUpdates: coupling)
        try arena.write(publicationRequests: publications)

        var dt = requestedStep
        var attempt: UInt32 = 0

        while attempt < configuration.maximumSubsteps {
            let result = try await executeAttempt(
                dt: dt,
                attemptIndex: attempt,
                couplingCount: coupling.count,
                publicationCount: publications.count
            )
            let status = result.status

            if status.flags.contains(.permanentShutdown) || status.requestedResponse == 5 {
                lifecycleState = .permanentlyStopped(reason: shutdownReason(status: status))
                return terminalPrepared(
                    transactionID: transactionID,
                    disposition: .permanentShutdown,
                    requestedStep: requestedStep,
                    candidateStep: dt,
                    attemptCount: attempt + 1,
                    result: result
                )
            }
            if status.flags.contains(.reversibleShutdown) || status.requestedResponse == 4 {
                lifecycleState = .reversiblyStopped(reason: shutdownReason(status: status))
                return terminalPrepared(
                    transactionID: transactionID,
                    disposition: .reversibleShutdown,
                    requestedStep: requestedStep,
                    candidateStep: dt,
                    attemptCount: attempt + 1,
                    result: result
                )
            }

            let unsupportedClamp = status.requestedResponse == 1
            let hardFailure = status.blocksCommit ||
                status.flags.contains(.invalidRate) ||
                status.flags.contains(.eventOverflow) ||
                unsupportedClamp
            let asksForSubdivision = status.flags.contains(.substep) ||
                status.requestedResponse == 3

            if !hardFailure && !asksForSubdivision {
                let prepared = VivoPreparedMolecularStep(
                    transactionID: transactionID,
                    disposition: .prepared,
                    sourceFingerprint: pack.header.sourceFingerprint,
                    programFingerprint: pack.header.contentFingerprint,
                    fidelity: configuration.fidelity,
                    stepIndex: nextStepIndex,
                    timeBefore: absoluteTime,
                    requestedTimeStep: requestedStep,
                    candidateTimeStep: dt,
                    attemptCount: attempt + 1,
                    status: status,
                    events: result.events,
                    publications: result.publications
                )
                pending = PendingStep(
                    prepared: prepared,
                    requestedTimeStep: requestedStep
                )
                return prepared
            }

            let reduced = dt * 0.5
            let mayRetry = request.permitAdaptiveReduction &&
                attempt + 1 < configuration.maximumSubsteps &&
                reduced >= configuration.minimumTimeStep &&
                reduced < dt
            guard mayRetry else {
                return VivoPreparedMolecularStep(
                    transactionID: transactionID,
                    disposition: reduced < dt ? .requiresSmallerStep : .rejected,
                    sourceFingerprint: pack.header.sourceFingerprint,
                    programFingerprint: pack.header.contentFingerprint,
                    fidelity: configuration.fidelity,
                    stepIndex: nextStepIndex,
                    timeBefore: absoluteTime,
                    requestedTimeStep: requestedStep,
                    candidateTimeStep: dt,
                    attemptCount: attempt + 1,
                    status: status,
                    events: result.events,
                    publications: []
                )
            }
            dt = reduced
            attempt &+= 1
        }

        throw VivoRuntimeError.commandFailed(
            "molecular prepare loop exhausted without producing a candidate"
        )
    }

    public func commitPreparedStep(
        transactionID: UUID
    ) throws -> VivoMolecularTransactionCertificate {
        guard let pending,
              pending.prepared.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular transaction identifier does not match the prepared candidate"
            )
        }
        guard pending.prepared.canCommit else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular candidate is not eligible for commit"
            )
        }

        arena.commit()
        let before = absoluteTime
        absoluteTime += pending.prepared.candidateTimeStep
        nextStepIndex &+= 1
        self.pending = nil

        return VivoMolecularTransactionCertificate(
            transactionID: transactionID,
            disposition: pending.prepared.candidateTimeStep == pending.requestedTimeStep
                ? .committed
                : .committedWithReducedStep,
            sourceFingerprint: pack.header.sourceFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            fidelity: configuration.fidelity,
            deviceName: capabilities.deviceName,
            deviceRegistryID: capabilities.registryID,
            stepIndex: pending.prepared.stepIndex,
            timeBefore: before,
            timeAfter: absoluteTime,
            requestedTimeStep: pending.requestedTimeStep,
            acceptedTimeStep: pending.prepared.candidateTimeStep,
            attemptCount: pending.prepared.attemptCount,
            status: pending.prepared.status
        )
    }

    public func discardPreparedStep(transactionID: UUID) throws {
        guard let pending,
              pending.prepared.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration(
                "molecular transaction identifier does not match the prepared candidate"
            )
        }
        self.pending = nil
    }

    public func step(
        _ request: VivoStepRequest = .init()
    ) async throws -> VivoMolecularTransactionResult {
        let transactionID = UUID()
        let prepared = try await prepareStep(request, transactionID: transactionID)
        guard prepared.canCommit else {
            return VivoMolecularTransactionResult(
                certificate: rejectedCertificate(for: prepared),
                events: prepared.events,
                publications: []
            )
        }
        let certificate = try commitPreparedStep(transactionID: transactionID)
        return VivoMolecularTransactionResult(
            certificate: certificate,
            events: prepared.events,
            publications: prepared.publications
        )
    }

    public func snapshot() async throws -> VivoStateSnapshot {
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration(
                "cannot snapshot molecular state while a candidate transaction is pending"
            )
        }
        guard let command = commandQueue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.TransactionalMolecularSnapshot"
        blit.copy(
            from: arena.currentState,
            sourceOffset: 0,
            to: arena.stateReadbackBuffer,
            destinationOffset: 0,
            size: arena.stateReadbackBuffer.length
        )
        blit.endEncoding()
        try await complete(command)

        let pointer = arena.stateReadbackBuffer.contents().assumingMemoryBound(to: Float.self)
        let values = Array(
            UnsafeBufferPointer(start: pointer, count: arena.capacities.stateElements)
        )
        return VivoStateSnapshot(
            sourceFingerprint: pack.header.sourceFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            stepIndex: nextStepIndex,
            absoluteTime: absoluteTime,
            speciesCount: pack.runtimeContract.speciesCount,
            laneCount: configuration.laneCount,
            values: values
        )
    }

    private func executeAttempt(
        dt: Float,
        attemptIndex: UInt32,
        couplingCount: Int,
        publicationCount: Int
    ) async throws -> AttemptResult {
        guard let command = commandQueue.makeCommandBuffer() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.TransactionalMolecular.Step.\(nextStepIndex).Attempt.\(attemptIndex)"

        let runtimeCommand = try VivoRuntimeCommandABI(
            pack: pack,
            configuration: configuration,
            stepIndex: nextStepIndex,
            substepIndex: attemptIndex,
            dt: dt,
            absoluteTime: absoluteTime + dt
        )

        try encodeClearStatus(command)
        try encodePrepare(command, runtimeCommand: runtimeCommand)
        if couplingCount > 0 {
            try encodeCoupling(
                command,
                runtimeCommand: runtimeCommand,
                count: couplingCount
            )
        }

        switch configuration.fidelity {
        case .logic:
            try encodeCopy(
                command,
                source: arena.baseState,
                destination: arena.candidateState
            )
        case .deterministic:
            try encodeDeterministicReaction(
                command,
                runtimeCommand: runtimeCommand,
                source: arena.baseState,
                predictor: arena.stageState,
                destination: arena.candidateState
            )
        case .stochastic:
            try encodeStochasticReaction(
                command,
                runtimeCommand: runtimeCommand,
                source: arena.baseState,
                destination: arena.candidateState
            )
        case .spatial, .tissue:
            try encodeSpatialSplit(command, runtimeCommand: runtimeCommand)
        }

        try encodeRules(command, runtimeCommand: runtimeCommand)
        try encodeMonitors(command, runtimeCommand: runtimeCommand)
        try encodeValidation(command, runtimeCommand: runtimeCommand)
        if publicationCount > 0 {
            try encodePublications(
                command,
                runtimeCommand: runtimeCommand,
                count: publicationCount
            )
        }

        try await complete(command)
        let status = arena.status()
        return AttemptResult(
            status: status,
            events: arena.events(status: status),
            publications: try arena.publicationValues(count: publicationCount)
        )
    }

    private func encodeSpatialSplit(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        var firstHalf = runtimeCommand
        firstHalf.dt = runtimeCommand.dt * 0.5
        firstHalf.absoluteTime = absoluteTime + firstHalf.dt

        if containsCountValuedSpecies {
            try encodeStochasticReaction(
                commandBuffer,
                runtimeCommand: firstHalf,
                source: arena.baseState,
                destination: arena.candidateState
            )
        } else {
            try encodeDeterministicReaction(
                commandBuffer,
                runtimeCommand: firstHalf,
                source: arena.baseState,
                predictor: arena.stageState,
                destination: arena.candidateState
            )
        }

        try encodeTransport(
            commandBuffer,
            runtimeCommand: runtimeCommand,
            source: arena.candidateState,
            destination: arena.stageState
        )

        var secondHalf = firstHalf
        secondHalf.absoluteTime = absoluteTime + runtimeCommand.dt
        if containsCountValuedSpecies {
            try encodeStochasticReaction(
                commandBuffer,
                runtimeCommand: secondHalf,
                source: arena.stageState,
                destination: arena.candidateState
            )
        } else {
            try encodeDeterministicReaction(
                commandBuffer,
                runtimeCommand: secondHalf,
                source: arena.stageState,
                predictor: arena.baseState,
                destination: arena.candidateState
            )
        }
    }

    private func encodeClearStatus(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(
            pipelines.clearStatus,
            count: 1,
            label: "NumiVivo.Transactional.ClearStatus",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 0)
        }
    }

    private func encodePrepare(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.prepare,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.Prepare",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.currentState, offset: 0, index: 0)
            encoder.setBuffer(arena.baseState, offset: 0, index: 1)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 2)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 3)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 6)
        }
    }

    private func encodeCoupling(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        count: Int
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.applyCoupling,
            count: count,
            label: "NumiVivo.Transactional.Coupling",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.baseState, offset: 0, index: 0)
            encoder.setBuffer(arena.couplingBuffer, offset: 0, index: 1)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 3)
            var count32 = UInt32(count)
            encoder.setBytes(&count32, length: MemoryLayout<UInt32>.stride, index: 4)
        }
    }

    private func encodeDeterministicReaction(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        predictor: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.heunPredict,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.HeunPredict",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(predictor, offset: 0, index: 1)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 2)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 6)
        }
        try encode(
            pipelines.heunCorrect,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.HeunCorrect",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(predictor, offset: 0, index: 1)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 2)
            encoder.setBuffer(destination, offset: 0, index: 3)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 6)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 7)
        }
    }

    private func encodeStochasticReaction(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        arena.write(command: runtimeCommand)
        let reactionWork = try checkedWorkItems(
            UInt64(pack.runtimeContract.reactionCount),
            UInt64(configuration.laneCount),
            label: "stochastic reaction"
        )
        try encode(
            pipelines.sampleReactions,
            count: reactionWork,
            label: "NumiVivo.Transactional.SampleReactions",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 1)
            encoder.setBuffer(arena.reactionEvents, offset: 0, index: 2)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 5)
        }
        try encode(
            pipelines.applyReactions,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.ApplyReactions",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(destination, offset: 0, index: 1)
            encoder.setBuffer(arena.reactionEvents, offset: 0, index: 2)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 5)
        }
    }

    private func encodeTransport(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.transport,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.Transport",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(destination, offset: 0, index: 1)
            encoder.setBuffer(arena.transport, offset: 0, index: 2)
            encoder.setBuffer(arena.velocity, offset: 0, index: 3)
            encoder.setBuffer(arena.volumeFraction, offset: 0, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 6)
        }
    }

    private func encodeRules(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        guard pack.runtimeContract.ruleCount > 0 else { return }
        arena.write(command: runtimeCommand)
        let work = try checkedWorkItems(
            UInt64(pack.runtimeContract.ruleCount),
            UInt64(configuration.laneCount),
            label: "rule"
        )
        try encode(
            pipelines.rules,
            count: work,
            label: "NumiVivo.Transactional.Rules",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 1)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.eventBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 6)
        }
    }

    private func encodeMonitors(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        guard pack.runtimeContract.monitorCount > 0 else { return }
        arena.write(command: runtimeCommand)
        let work = try checkedWorkItems(
            UInt64(pack.runtimeContract.monitorCount),
            UInt64(configuration.laneCount),
            label: "monitor"
        )
        try encode(
            pipelines.monitors,
            count: work,
            label: "NumiVivo.Transactional.Monitors",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 1)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.eventBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 6)
        }
    }

    private func encodeValidation(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.validate,
            count: arena.capacities.stateElements,
            label: "NumiVivo.Transactional.Validate",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.programBuffer, offset: 0, index: 1)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 3)
        }
    }

    private func encodePublications(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        count: Int
    ) throws {
        arena.write(command: runtimeCommand)
        try encode(
            pipelines.publish,
            count: count,
            label: "NumiVivo.Transactional.Publish",
            commandBuffer: commandBuffer
        ) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.publicationRequestBuffer, offset: 0, index: 1)
            encoder.setBuffer(arena.publicationOutputBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 3)
            encoder.setBuffer(arena.commandBuffer, offset: 0, index: 4)
            var count32 = UInt32(count)
            encoder.setBytes(&count32, length: MemoryLayout<UInt32>.stride, index: 5)
        }
    }

    private func encodeCopy(
        _ commandBuffer: MTLCommandBuffer,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        guard let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandFailed("could not create blit encoder")
        }
        blit.label = "NumiVivo.Transactional.CopyState"
        blit.copy(
            from: source,
            sourceOffset: 0,
            to: destination,
            destinationOffset: 0,
            size: arena.capacities.stateElements * MemoryLayout<Float>.stride
        )
        blit.endEncoding()
    }

    private func encode(
        _ pipeline: NumiVivoPipeline,
        count: Int,
        label: String,
        commandBuffer: MTLCommandBuffer,
        bindings: (MTLComputeCommandEncoder) throws -> Void
    ) throws {
        guard count > 0 else { return }
        guard let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandFailed(
                "could not create compute encoder for \(label)"
            )
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline.state)
        try bindings(encoder)
        encoder.dispatchThreads(
            pipeline.gridSize(for: count),
            threadsPerThreadgroup: pipeline.threadgroupSize(
                for: count,
                preferred: capabilities.recommendedThreadsPerThreadgroup
            )
        )
        encoder.endEncoding()
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if let error = completed.error {
                    continuation.resume(
                        throwing: VivoRuntimeError.commandFailed(String(describing: error))
                    )
                } else {
                    continuation.resume(returning: ())
                }
            }
            command.commit()
        }
    }

    private func validate(request: VivoStepRequest) throws {
        guard request.coupling.count <= arena.capacities.couplingUpdates else {
            throw VivoRuntimeError.invalidConfiguration(
                "coupling request exceeds the allocated capacity"
            )
        }
        guard request.publications.count <= arena.capacities.publicationRequests else {
            throw VivoRuntimeError.invalidConfiguration(
                "publication request exceeds the allocated capacity"
            )
        }
        for update in request.coupling {
            guard update.speciesIndex < pack.runtimeContract.speciesCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "coupling species index is outside the ProgramPack"
                )
            }
            guard update.laneIndex < configuration.laneCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "coupling lane index is outside the runtime"
                )
            }
            guard update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration(
                    "coupling value must be finite"
                )
            }
        }
        for publication in request.publications {
            guard publication.speciesIndex < pack.runtimeContract.speciesCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "publication species index is outside the ProgramPack"
                )
            }
            guard publication.laneIndex < configuration.laneCount else {
                throw VivoRuntimeError.invalidConfiguration(
                    "publication lane index is outside the runtime"
                )
            }
        }
    }

    private func ensureReadyForMutation() throws {
        switch lifecycleState {
        case .ready:
            return
        case .reversiblyStopped(let reason):
            throw VivoRuntimeError.runtimeStopped("reversible shutdown: \(reason)")
        case .permanentlyStopped(let reason):
            throw VivoRuntimeError.runtimeStopped("permanent shutdown: \(reason)")
        case .failed(let reason):
            throw VivoRuntimeError.runtimeStopped("runtime failure: \(reason)")
        }
    }

    private func terminalPrepared(
        transactionID: UUID,
        disposition: VivoPreparedMolecularDisposition,
        requestedStep: Float,
        candidateStep: Float,
        attemptCount: UInt32,
        result: AttemptResult
    ) -> VivoPreparedMolecularStep {
        VivoPreparedMolecularStep(
            transactionID: transactionID,
            disposition: disposition,
            sourceFingerprint: pack.header.sourceFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            fidelity: configuration.fidelity,
            stepIndex: nextStepIndex,
            timeBefore: absoluteTime,
            requestedTimeStep: requestedStep,
            candidateTimeStep: candidateStep,
            attemptCount: attemptCount,
            status: result.status,
            events: result.events,
            publications: []
        )
    }

    private func rejectedCertificate(
        for prepared: VivoPreparedMolecularStep
    ) -> VivoMolecularTransactionCertificate {
        let disposition: VivoStepDisposition
        switch prepared.disposition {
        case .permanentShutdown:
            disposition = .permanentShutdown
        case .reversibleShutdown:
            disposition = .reversibleShutdown
        case .prepared, .requiresSmallerStep, .rejected:
            disposition = .rejected
        }
        return VivoMolecularTransactionCertificate(
            transactionID: prepared.transactionID,
            disposition: disposition,
            sourceFingerprint: pack.header.sourceFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            fidelity: configuration.fidelity,
            deviceName: capabilities.deviceName,
            deviceRegistryID: capabilities.registryID,
            stepIndex: prepared.stepIndex,
            timeBefore: prepared.timeBefore,
            timeAfter: prepared.timeBefore,
            requestedTimeStep: prepared.requestedTimeStep,
            acceptedTimeStep: nil,
            attemptCount: prepared.attemptCount,
            status: prepared.status
        )
    }

    private func shutdownReason(status: VivoRuntimeStatus) -> String {
        if let monitor = status.firstFailingMonitor {
            return "monitor \(monitor) requested response \(status.requestedResponse)"
        }
        return "runtime requested response \(status.requestedResponse)"
    }

    private func checkedWorkItems(
        _ left: UInt64,
        _ right: UInt64,
        label: String
    ) throws -> Int {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow,
              result.partialValue <= UInt64(Int.max) else {
            throw VivoRuntimeError.invalidConfiguration(
                "\(label) work-item count overflow"
            )
        }
        return Int(result.partialValue)
    }
}

private struct VivoTransactionalRuntimePipelines: Sendable {
    let clearStatus: NumiVivoPipeline
    let prepare: NumiVivoPipeline
    let applyCoupling: NumiVivoPipeline
    let heunPredict: NumiVivoPipeline
    let heunCorrect: NumiVivoPipeline
    let sampleReactions: NumiVivoPipeline
    let applyReactions: NumiVivoPipeline
    let transport: NumiVivoPipeline
    let rules: NumiVivoPipeline
    let monitors: NumiVivoPipeline
    let validate: NumiVivoPipeline
    let publish: NumiVivoPipeline

    static func load(
        from catalog: NumiVivoPipelineCatalog
    ) async throws -> VivoTransactionalRuntimePipelines {
        async let clearStatus = catalog.pipeline(.clearStatus)
        async let prepare = catalog.pipeline(.prepareTransaction)
        async let applyCoupling = catalog.pipeline(.applyCouplingUpdates)
        async let heunPredict = catalog.pipeline(.f1HeunPredict)
        async let heunCorrect = catalog.pipeline(.f1HeunCorrect)
        async let sampleReactions = catalog.pipeline(.f2SampleReactions)
        async let applyReactions = catalog.pipeline(.f2ApplyReactions)
        async let transport = catalog.pipeline(.f3Transport)
        async let rules = catalog.pipeline(.executeRules)
        async let monitors = catalog.pipeline(.evaluateMonitors)
        async let validate = catalog.pipeline(.validateShadow)
        async let publish = catalog.pipeline(.publish)
        return try await VivoTransactionalRuntimePipelines(
            clearStatus: clearStatus,
            prepare: prepare,
            applyCoupling: applyCoupling,
            heunPredict: heunPredict,
            heunCorrect: heunCorrect,
            sampleReactions: sampleReactions,
            applyReactions: applyReactions,
            transport: transport,
            rules: rules,
            monitors: monitors,
            validate: validate,
            publish: publish
        )
    }
}
