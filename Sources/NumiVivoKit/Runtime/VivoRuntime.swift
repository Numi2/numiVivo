import Foundation
import Metal
import NumiVivoShaders

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

    public init(
        speciesIndex: UInt32,
        laneIndex: UInt32,
        mode: VivoCouplingMode = .replace,
        value: Float
    ) {
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

    public init(
        timeStep: Float? = nil,
        coupling: [VivoCouplingUpdate] = [],
        publications: [VivoPublicationRequest] = [],
        permitAdaptiveReduction: Bool = true
    ) {
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

    public var committed: Bool {
        disposition == .committed || disposition == .committedWithReducedStep
    }
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
        return values[Int(species * laneCount + lane)]
    }
}

public actor VivoRuntime {
    public enum Lifecycle: Sendable, Codable, Equatable {
        case ready
        case reversiblyStopped(reason: String)
        case permanentlyStopped(reason: String)
        case failed(reason: String)
    }

    public nonisolated let pack: VivoProgramPack
    public nonisolated let configuration: VivoRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities

    private let device: MTLDevice
    private let commandQueue: MTLCommandQueue
    private let arena: VivoMetalArena
    private let pipelines: RuntimePipelines
    private let containsCountValuedSpecies: Bool

    private var lifecycleState: Lifecycle = .ready
    private var absoluteTime: Float = 0
    private var nextStepIndex: UInt32 = 0

    public static func make(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoRuntime {
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard let queue = device.makeCommandQueue() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.RuntimeQueue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let pipelines = try await RuntimePipelines.load(from: catalog)
        let arena = try VivoMetalArena(
            device: device,
            commandQueue: queue,
            pack: pack,
            configuration: configuration
        )
        let countValued = try pack.speciesMetadata().contains(where: \ .isCountValued)
        return VivoRuntime(
            pack: pack,
            configuration: configuration,
            device: device,
            queue: queue,
            arena: arena,
            pipelines: pipelines,
            containsCountValuedSpecies: countValued
        )
    }

    private init(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device: MTLDevice,
        queue: MTLCommandQueue,
        arena: VivoMetalArena,
        pipelines: RuntimePipelines,
        containsCountValuedSpecies: Bool
    ) {
        self.pack = pack
        self.configuration = configuration
        self.device = device
        self.commandQueue = queue
        self.arena = arena
        self.pipelines = pipelines
        self.capabilities = arena.capabilities
        self.containsCountValuedSpecies = containsCountValuedSpecies
    }

    public func lifecycle() -> Lifecycle { lifecycleState }
    public func time() -> Float { absoluteTime }
    public func stepIndex() -> UInt32 { nextStepIndex }

    public func resume() throws {
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
        lifecycleState = .permanentlyStopped(reason: reason)
    }

    public func setTransport(_ values: [VivoSpeciesTransportABI]) throws {
        try ensureReadyForMutation()
        try arena.replaceTransport(values, commandQueue: commandQueue)
    }

    public func setVelocity(_ values: [SIMD4<Float>]) throws {
        try ensureReadyForMutation()
        try arena.replaceVelocity(values, commandQueue: commandQueue)
    }

    public func setVolumeFractions(_ values: [Float]) throws {
        try ensureReadyForMutation()
        try arena.replaceVolumeFraction(values, commandQueue: commandQueue)
    }

    public func step(_ request: VivoStepRequest = .init()) async throws -> VivoStepResult {
        try ensureReadyForMutation()
        let requestedStep = request.timeStep ?? configuration.timeStep
        guard requestedStep.isFinite,
              requestedStep >= configuration.minimumTimeStep,
              requestedStep <= configuration.maximumTimeStep else {
            throw VivoRuntimeError.invalidConfiguration("requested time step is outside configured finite bounds")
        }
        try validate(request: request)

        let coupling = request.coupling.map(\ .abi)
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
        var lastAttempt: AttemptResult?

        while attempt < configuration.maximumSubsteps {
            let result = try await executeAttempt(
                dt: dt,
                attemptIndex: attempt,
                couplingCount: coupling.count,
                publicationCount: publications.count
            )
            lastAttempt = result
            let status = result.status

            if status.flags.contains(.permanentShutdown) || status.requestedResponse == 5 {
                lifecycleState = .permanentlyStopped(reason: shutdownReason(status: status))
                return makeResult(
                    disposition: .permanentShutdown,
                    requestedStep: requestedStep,
                    acceptedStep: nil,
                    attemptCount: attempt + 1,
                    status: status,
                    events: result.events,
                    publications: []
                )
            }
            if status.flags.contains(.reversibleShutdown) || status.requestedResponse == 4 {
                lifecycleState = .reversiblyStopped(reason: shutdownReason(status: status))
                return makeResult(
                    disposition: .reversibleShutdown,
                    requestedStep: requestedStep,
                    acceptedStep: nil,
                    attemptCount: attempt + 1,
                    status: status,
                    events: result.events,
                    publications: []
                )
            }

            let unsupportedClamp = status.requestedResponse == 1
            let hardFailure = status.blocksCommit ||
                status.flags.contains(.invalidRate) ||
                status.flags.contains(.eventOverflow) ||
                unsupportedClamp
            let asksForSubdivision = status.flags.contains(.substep) || status.requestedResponse == 3
            if !hardFailure && !asksForSubdivision {
                arena.commit()
                let before = absoluteTime
                absoluteTime += dt
                nextStepIndex &+= 1
                let disposition: VivoStepDisposition = dt == requestedStep ? .committed : .committedWithReducedStep
                let certificate = VivoStepCertificate(
                    disposition: disposition,
                    sourceFingerprint: pack.header.sourceFingerprint,
                    programFingerprint: pack.header.contentFingerprint,
                    fidelity: configuration.fidelity,
                    deviceName: capabilities.deviceName,
                    deviceRegistryID: capabilities.registryID,
                    stepIndex: nextStepIndex - 1,
                    attemptedStep: requestedStep,
                    acceptedStep: dt,
                    attemptCount: attempt + 1,
                    timeBefore: before,
                    timeAfter: absoluteTime,
                    status: status
                )
                return VivoStepResult(
                    certificate: certificate,
                    events: result.events,
                    publications: result.publications
                )
            }

            let reduced = dt * 0.5
            let mayRetry = request.permitAdaptiveReduction &&
                attempt + 1 < configuration.maximumSubsteps &&
                reduced >= configuration.minimumTimeStep &&
                reduced < dt
            guard mayRetry else {
                return makeResult(
                    disposition: .rejected,
                    requestedStep: requestedStep,
                    acceptedStep: nil,
                    attemptCount: attempt + 1,
                    status: status,
                    events: result.events,
                    publications: []
                )
            }
            dt = reduced
            attempt &+= 1
        }

        let status = lastAttempt?.status ?? VivoRuntimeStatus(raw: .init())
        return makeResult(
            disposition: .rejected,
            requestedStep: requestedStep,
            acceptedStep: nil,
            attemptCount: attempt,
            status: status,
            events: lastAttempt?.events ?? [],
            publications: []
        )
    }

    public func snapshot() async throws -> VivoStateSnapshot {
        guard let command = commandQueue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.StateSnapshot"
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
        let values = Array(UnsafeBufferPointer(start: pointer, count: arena.capacities.stateElements))
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
        command.label = "NumiVivo.Step.\(nextStepIndex).Attempt.\(attemptIndex)"

        let finalCommand = try VivoRuntimeCommandABI(
            pack: pack,
            configuration: configuration,
            stepIndex: nextStepIndex,
            substepIndex: attemptIndex,
            dt: dt,
            absoluteTime: absoluteTime + dt
        )

        try encodeClearStatus(command)
        try encodePrepare(command, runtimeCommand: finalCommand)
        if couplingCount > 0 {
            try encodeCoupling(command, runtimeCommand: finalCommand, count: couplingCount)
        }

        switch configuration.fidelity {
        case .logic:
            try encodeCopy(command, source: arena.baseState, destination: arena.candidateState)
        case .deterministic:
            try encodeDeterministicReaction(
                command,
                runtimeCommand: finalCommand,
                source: arena.baseState,
                predictor: arena.stageState,
                destination: arena.candidateState
            )
        case .stochastic:
            try encodeStochasticReaction(
                command,
                runtimeCommand: finalCommand,
                source: arena.baseState,
                destination: arena.candidateState
            )
        case .spatial, .tissue:
            try encodeSpatialSplit(command, runtimeCommand: finalCommand)
        }

        try encodeRules(command, runtimeCommand: finalCommand)
        try encodeMonitors(command, runtimeCommand: finalCommand)
        try encodeValidation(command, runtimeCommand: finalCommand)
        if publicationCount > 0 {
            try encodePublications(command, runtimeCommand: finalCommand, count: publicationCount)
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
        let stochastic = containsCountValuedSpecies

        if stochastic {
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
        if stochastic {
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
        try encode(pipelines.clearStatus, count: 1, label: "NumiVivo.ClearStatus", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 0)
        }
    }

    private func encodePrepare(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        let count = max(arena.capacities.stateElements, arena.capacities.temporalElements)
        try encode(pipelines.prepareTransaction, count: count, label: "NumiVivo.Prepare", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.currentState, offset: 0, index: 0)
            encoder.setBuffer(arena.baseState, offset: 0, index: 1)
            encoder.setBuffer(arena.stageState, offset: 0, index: 2)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 3)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 4)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 5)
            set(runtimeCommand, on: encoder, index: 6)
        }
    }

    private func encodeCoupling(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        count: Int
    ) throws {
        var count32 = UInt32(count)
        try encode(pipelines.applyCouplingUpdates, count: count, label: "NumiVivo.Coupling", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.couplingBuffer, offset: 0, index: 0)
            encoder.setBytes(&count32, length: MemoryLayout<UInt32>.stride, index: 1)
            encoder.setBuffer(arena.baseState, offset: 0, index: 2)
            try arena.bindProgramSection(.species, to: encoder, index: 3)
            set(runtimeCommand, on: encoder, index: 4)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 5)
        }
    }

    private func encodeDeterministicReaction(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        predictor: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        let count = arena.capacities.stateElements
        try encode(pipelines.f1HeunPredict, count: count, label: "NumiVivo.F1.Predict", commandBuffer: commandBuffer) { encoder in
            try bindReactionProgram(to: encoder)
            encoder.setBuffer(source, offset: 0, index: 7)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 8)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 9)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 10)
            encoder.setBuffer(predictor, offset: 0, index: 11)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 12)
            set(runtimeCommand, on: encoder, index: 13)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 14)
        }
        try encode(pipelines.f1HeunCorrect, count: count, label: "NumiVivo.F1.Correct", commandBuffer: commandBuffer) { encoder in
            try bindReactionProgram(to: encoder)
            encoder.setBuffer(source, offset: 0, index: 7)
            encoder.setBuffer(predictor, offset: 0, index: 8)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 9)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 10)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 11)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 12)
            encoder.setBuffer(destination, offset: 0, index: 13)
            set(runtimeCommand, on: encoder, index: 14)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 15)
        }
    }

    private func encodeStochasticReaction(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        let eventCount = arena.capacities.reactionEventElements
        try encode(pipelines.f2SampleReactions, count: eventCount, label: "NumiVivo.F2.Sample", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.reactions, to: encoder, index: 0)
            try arena.bindProgramSection(.stoichiometry, to: encoder, index: 1)
            try arena.bindProgramSection(.reactionParameterIndices, to: encoder, index: 2)
            try arena.bindProgramSection(.expressions, to: encoder, index: 3)
            encoder.setBuffer(source, offset: 0, index: 4)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 5)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 6)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 7)
            encoder.setBuffer(arena.reactionEvents, offset: 0, index: 8)
            set(runtimeCommand, on: encoder, index: 9)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 10)
        }
        try encode(pipelines.f2ApplyReactions, count: arena.capacities.stateElements, label: "NumiVivo.F2.Apply", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.species, to: encoder, index: 0)
            try arena.bindProgramSection(.speciesIncidenceOffsets, to: encoder, index: 1)
            try arena.bindProgramSection(.speciesIncidence, to: encoder, index: 2)
            encoder.setBuffer(arena.reactionEvents, offset: 0, index: 3)
            encoder.setBuffer(source, offset: 0, index: 4)
            encoder.setBuffer(destination, offset: 0, index: 5)
            set(runtimeCommand, on: encoder, index: 6)
        }
    }

    private func encodeTransport(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        try encode(pipelines.f3Transport, count: arena.capacities.stateElements, label: "NumiVivo.F3.Transport", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.species, to: encoder, index: 0)
            encoder.setBuffer(arena.transport, offset: 0, index: 1)
            encoder.setBuffer(source, offset: 0, index: 2)
            encoder.setBuffer(destination, offset: 0, index: 3)
            encoder.setBuffer(arena.velocity, offset: 0, index: 4)
            encoder.setBuffer(arena.volumeFraction, offset: 0, index: 5)
            set(runtimeCommand, on: encoder, index: 6)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 7)
        }
    }

    private func encodeRules(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        guard pack.runtimeContract.ruleCount > 0 else { return }
        try encode(pipelines.executeRules, count: Int(configuration.laneCount), label: "NumiVivo.Rules", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.species, to: encoder, index: 0)
            try arena.bindProgramSection(.rules, to: encoder, index: 1)
            try arena.bindProgramSection(.actions, to: encoder, index: 2)
            try arena.bindProgramSection(.expressions, to: encoder, index: 3)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 4)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 5)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 6)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 7)
            encoder.setBuffer(arena.eventBuffer, offset: 0, index: 8)
            set(runtimeCommand, on: encoder, index: 9)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 10)
        }
    }

    private func encodeMonitors(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        guard pack.runtimeContract.monitorCount > 0 else { return }
        try encode(pipelines.evaluateMonitors, count: Int(configuration.laneCount), label: "NumiVivo.Monitors", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.monitors, to: encoder, index: 0)
            try arena.bindProgramSection(.expressions, to: encoder, index: 1)
            encoder.setBuffer(arena.parameterBuffer, offset: 0, index: 2)
            encoder.setBuffer(arena.temporalCurrent, offset: 0, index: 3)
            encoder.setBuffer(arena.temporalCandidate, offset: 0, index: 4)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 5)
            set(runtimeCommand, on: encoder, index: 6)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 7)
        }
    }

    private func encodeValidation(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI
    ) throws {
        try encode(pipelines.validateShadow, count: arena.capacities.stateElements, label: "NumiVivo.Validate", commandBuffer: commandBuffer) { encoder in
            try arena.bindProgramSection(.species, to: encoder, index: 0)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 1)
            set(runtimeCommand, on: encoder, index: 2)
            encoder.setBuffer(arena.statusBuffer, offset: 0, index: 3)
        }
    }

    private func encodePublications(
        _ commandBuffer: MTLCommandBuffer,
        runtimeCommand: VivoRuntimeCommandABI,
        count: Int
    ) throws {
        var count32 = UInt32(count)
        try encode(pipelines.publish, count: count, label: "NumiVivo.Publish", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.publicationRequestBuffer, offset: 0, index: 0)
            encoder.setBytes(&count32, length: MemoryLayout<UInt32>.stride, index: 1)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 2)
            encoder.setBuffer(arena.publicationOutputBuffer, offset: 0, index: 3)
            set(runtimeCommand, on: encoder, index: 4)
        }
    }

    private func bindReactionProgram(to encoder: MTLComputeCommandEncoder) throws {
        try arena.bindProgramSection(.species, to: encoder, index: 0)
        try arena.bindProgramSection(.reactions, to: encoder, index: 1)
        try arena.bindProgramSection(.stoichiometry, to: encoder, index: 2)
        try arena.bindProgramSection(.reactionParameterIndices, to: encoder, index: 3)
        try arena.bindProgramSection(.expressions, to: encoder, index: 4)
        try arena.bindProgramSection(.speciesIncidenceOffsets, to: encoder, index: 5)
        try arena.bindProgramSection(.speciesIncidence, to: encoder, index: 6)
    }

    private func encodeCopy(
        _ commandBuffer: MTLCommandBuffer,
        source: MTLBuffer,
        destination: MTLBuffer
    ) throws {
        guard let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.pipelineEncodingFailed("could not create state-copy blit encoder")
        }
        blit.label = "NumiVivo.CopyCandidate"
        blit.copy(from: source, sourceOffset: 0, to: destination, destinationOffset: 0, size: arena.stateReadbackBuffer.length)
        blit.endEncoding()
    }

    private func encode(
        _ pipeline: NumiVivoPipeline,
        count: Int,
        label: String,
        commandBuffer: MTLCommandBuffer,
        body: (MTLComputeCommandEncoder) throws -> Void
    ) throws {
        guard count > 0 else { return }
        guard let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.pipelineEncodingFailed("could not create encoder for \(label)")
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline.state)
        do {
            try body(encoder)
            let preferred = capabilities.threadgroupSize(
                elementCount: count,
                executionWidth: pipeline.executionWidth,
                pipelineMaximum: pipeline.maximumThreadsPerThreadgroup
            )
            encoder.dispatchThreads(
                MTLSize(width: count, height: 1, depth: 1),
                threadsPerThreadgroup: preferred
            )
            encoder.endEncoding()
        } catch {
            encoder.endEncoding()
            throw error
        }
    }

    private func set(
        _ command: VivoRuntimeCommandABI,
        on encoder: MTLComputeCommandEncoder,
        index: Int
    ) {
        var command = command
        encoder.setBytes(&command, length: MemoryLayout<VivoRuntimeCommandABI>.stride, index: index)
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if let error = completed.error {
                    continuation.resume(throwing: VivoRuntimeError.commandFailed(String(describing: error)))
                } else {
                    continuation.resume(returning: ())
                }
            }
            command.commit()
        }
    }

    private func validate(request: VivoStepRequest) throws {
        for update in request.coupling {
            guard update.speciesIndex < pack.runtimeContract.speciesCount,
                  update.laneIndex < configuration.laneCount,
                  update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("coupling update references an invalid species/lane or non-finite value")
            }
        }
        guard request.publications.count <= Int(UInt32.max) else {
            throw VivoRuntimeError.invalidConfiguration("publication request count exceeds UInt32")
        }
        for publication in request.publications {
            guard publication.speciesIndex < pack.runtimeContract.speciesCount,
                  publication.laneIndex < configuration.laneCount else {
                throw VivoRuntimeError.invalidConfiguration("publication request references an invalid species or lane")
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

    private func shutdownReason(status: VivoRuntimeStatus) -> String {
        if let monitor = status.firstMonitor {
            return "monitor \(monitor) requested response \(status.requestedResponse)"
        }
        return "runtime status flags 0x\(String(status.flags.rawValue, radix: 16))"
    }

    private func makeResult(
        disposition: VivoStepDisposition,
        requestedStep: Float,
        acceptedStep: Float?,
        attemptCount: UInt32,
        status: VivoRuntimeStatus,
        events: [VivoEvent],
        publications: [Float]
    ) -> VivoStepResult {
        let certificate = VivoStepCertificate(
            disposition: disposition,
            sourceFingerprint: pack.header.sourceFingerprint,
            programFingerprint: pack.header.contentFingerprint,
            fidelity: configuration.fidelity,
            deviceName: capabilities.deviceName,
            deviceRegistryID: capabilities.registryID,
            stepIndex: nextStepIndex,
            attemptedStep: requestedStep,
            acceptedStep: acceptedStep,
            attemptCount: attemptCount,
            timeBefore: absoluteTime,
            timeAfter: acceptedStep.map { absoluteTime + $0 } ?? absoluteTime,
            status: status
        )
        return VivoStepResult(certificate: certificate, events: events, publications: publications)
    }
}

private struct AttemptResult {
    let status: VivoRuntimeStatus
    let events: [VivoEvent]
    let publications: [Float]
}

private struct RuntimePipelines: Sendable {
    let clearStatus: NumiVivoPipeline
    let prepareTransaction: NumiVivoPipeline
    let applyCouplingUpdates: NumiVivoPipeline
    let f1HeunPredict: NumiVivoPipeline
    let f1HeunCorrect: NumiVivoPipeline
    let f2SampleReactions: NumiVivoPipeline
    let f2ApplyReactions: NumiVivoPipeline
    let f3Transport: NumiVivoPipeline
    let executeRules: NumiVivoPipeline
    let evaluateMonitors: NumiVivoPipeline
    let validateShadow: NumiVivoPipeline
    let publish: NumiVivoPipeline

    static func load(from catalog: NumiVivoPipelineCatalog) async throws -> RuntimePipelines {
        async let clearStatus = catalog.pipeline(.clearStatus)
        async let prepareTransaction = catalog.pipeline(.prepareTransaction)
        async let applyCouplingUpdates = catalog.pipeline(.applyCouplingUpdates)
        async let f1HeunPredict = catalog.pipeline(.f1HeunPredict)
        async let f1HeunCorrect = catalog.pipeline(.f1HeunCorrect)
        async let f2SampleReactions = catalog.pipeline(.f2SampleReactions)
        async let f2ApplyReactions = catalog.pipeline(.f2ApplyReactions)
        async let f3Transport = catalog.pipeline(.f3Transport)
        async let executeRules = catalog.pipeline(.executeRules)
        async let evaluateMonitors = catalog.pipeline(.evaluateMonitors)
        async let validateShadow = catalog.pipeline(.validateShadow)
        async let publish = catalog.pipeline(.publish)
        return try await RuntimePipelines(
            clearStatus: clearStatus,
            prepareTransaction: prepareTransaction,
            applyCouplingUpdates: applyCouplingUpdates,
            f1HeunPredict: f1HeunPredict,
            f1HeunCorrect: f1HeunCorrect,
            f2SampleReactions: f2SampleReactions,
            f2ApplyReactions: f2ApplyReactions,
            f3Transport: f3Transport,
            executeRules: executeRules,
            evaluateMonitors: evaluateMonitors,
            validateShadow: validateShadow,
            publish: publish
        )
    }
}
