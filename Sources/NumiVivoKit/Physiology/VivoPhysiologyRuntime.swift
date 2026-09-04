import Foundation
import Metal
import NumiVivoShaders

public actor VivoPhysiologyRuntime {
    private struct Pending {
        let prepared: VivoPreparedPhysiologyStep
        let nextDoseCursor: Int
    }
    private struct Output {
        let status: VivoPhysiologyRuntimeStatus
        let publications: [Float]
    }
    private struct Accumulator {
        var replacement: Float?
        var delta: Double = 0
        var minimum: Float?
        var maximum: Float?
    }
    private struct DosePlan {
        let pre: [VivoPhysiologyStateUpdate]
        let post: [VivoPhysiologyStateUpdate]
        let identifiers: [String]
        let nextCursor: Int
        let nextBoundary: Double?
    }
    public nonisolated let model: PreparedVivoPhysiologyModel
    public nonisolated let configuration: VivoPhysiologyRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities
    private let device: MTLDevice
    private let commandQueue: MTLCommandQueue
    private var arena: VivoPhysiologyMetalArena
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private var absoluteTimeSeconds: Double = 0
    private var nextStepIndex: UInt32 = 0
    private var doseCursor = 0
    private var pending: Pending?
    private var inFlight = false
    private var failure: String?

    public static func make(model: PreparedVivoPhysiologyModel,
                            configuration: VivoPhysiologyRuntimeConfiguration = .init(),
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoPhysiologyRuntime {
        try VivoPhysiologyModelValidator.validate(model)
        try configuration.validate(for: model)
        guard configuration.maximumTransformsPerStep <= Int(UInt32.max), configuration.maximumPublicationsPerStep <= Int(UInt32.max) else {
            throw VivoRuntimeError.invalidConfiguration("physiology boundary capacity exceeds UInt32")
        }
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard let queue = device.makeCommandQueue() else { throw VivoRuntimeError.commandQueueUnavailable }
        let catalog = try NumiVivoPipelineCatalog(device: device)
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in [NumiVivoKernel.physiologyClearStatus, .physiologyPrepareTransaction,
                       .physiologyApplyTransforms, .physiologyHeunPredict, .physiologyHeunCorrect,
                       .physiologyValidateCandidate, .physiologyPublish] {
            pipelines[kernel] = try await catalog.pipeline(kernel)
        }
        let arena = try VivoPhysiologyMetalArena(device: device, commandQueue: queue, model: model, configuration: configuration)
        return .init(model: model, configuration: configuration, device: device, queue: queue, arena: arena, pipelines: pipelines)
    }
    private init(model: PreparedVivoPhysiologyModel, configuration: VivoPhysiologyRuntimeConfiguration,
                 device: MTLDevice, queue: MTLCommandQueue, arena: VivoPhysiologyMetalArena,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline]) {
        self.model = model; self.configuration = configuration; self.device = device
        commandQueue = queue; self.arena = arena; self.pipelines = pipelines
        capabilities = arena.capabilities
    }
    public func time() -> Double { absoluteTimeSeconds }
    public func stepIndex() -> UInt32 { nextStepIndex }
    public func hasPendingTransaction() -> Bool { inFlight || pending != nil }
    private func requireIdle() throws {
        if let failure { throw VivoRuntimeError.runtimeStopped(failure) }
        guard !inFlight, pending == nil else { throw VivoRuntimeError.invalidConfiguration("physiology operation or candidate already active") }
    }
    public func prepareStep(_ request: VivoPhysiologyPrepareRequest = .init(), transactionID: UUID = UUID()) async throws -> VivoPreparedPhysiologyStep {
        try requireIdle()
        try Task.checkCancellation()
        try validate(request)
        let requested = request.timeStepSeconds ?? model.preferredTimeStepSeconds
        guard requested.isFinite, requested >= model.minimumTimeStepSeconds, requested <= model.maximumTimeStepSeconds,
              nextStepIndex < UInt32.max else { throw VivoRuntimeError.invalidConfiguration("invalid physiology timestep or exhausted step index") }
        inFlight = true
        defer { inFlight = false }
        var dt = try alignedTimeStep(requested)
        for attempt in 0..<configuration.maximumSubsteps {
            let doses = try dosePlan(dt: dt)
            let pre = try normalize(doses.pre + request.preUpdates, dt: dt)
            let post = try normalize(request.postUpdates + doses.post, dt: dt)
            let publications = request.publications.enumerated().map { index, value in
                VivoPhysiologyPublicationRequestABI(pairIndex: value.pairIndex, environmentIndex: value.environmentIndex,
                                                     outputIndex: UInt32(index), flags: value.flags)
            }
            let result = try await execute(transactionID: transactionID, dt: dt, pre: pre, post: post,
                                           publications: publications, appliedDoseCount: doses.identifiers.count)
            try Task.checkCancellation()
            let retry = result.status.blocksCommit || result.status.flags.contains(.requiresSubstep) || result.status.flags.contains(.excessiveDerivative)
            if retry {
                let suggested = Double(result.status.suggestedTimeStep)
                let reduced = min(dt * 0.5, suggested.isFinite && suggested > 0 ? suggested : dt * 0.5)
                if request.permitAdaptiveReduction, attempt + 1 < configuration.maximumSubsteps,
                   reduced >= model.minimumTimeStepSeconds, reduced < dt {
                    dt = try alignedTimeStep(reduced)
                    continue
                }
            }
            let value = VivoPreparedPhysiologyStep(transactionID: transactionID, disposition: retry ? .requiresSmallerStep : .prepared,
                                                    modelFingerprint: model.fingerprint, stepIndex: nextStepIndex,
                                                    timeBefore: absoluteTimeSeconds, requestedTimeStep: requested,
                                                    candidateTimeStep: dt, attemptCount: attempt + 1,
                                                    status: result.status, publications: retry ? [] : result.publications,
                                                    appliedDoseIdentifiers: doses.identifiers, nextBoundarySeconds: doses.nextBoundary)
            if value.canCommit { pending = .init(prepared: value, nextDoseCursor: doses.nextCursor) }
            return value
        }
        throw VivoRuntimeError.commandFailed("physiology attempt budget exhausted")
    }
    public func commitPreparedStep(transactionID: UUID) throws -> VivoPhysiologyStepCertificate {
        guard !inFlight, failure == nil, let candidate = pending,
              candidate.prepared.transactionID == transactionID, candidate.prepared.canCommit else {
            throw VivoRuntimeError.invalidConfiguration("physiology transaction is not eligible for commit")
        }
        arena.commit()
        absoluteTimeSeconds += candidate.prepared.candidateTimeStep
        doseCursor = candidate.nextDoseCursor
        nextStepIndex += 1
        pending = nil
        return certificate(candidate.prepared, committed: true)
    }
    public func discardPreparedStep(transactionID: UUID) throws {
        guard !inFlight, pending?.prepared.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration("physiology transaction identifier mismatch")
        }
        pending = nil
    }
    public func step(_ request: VivoPhysiologyPrepareRequest = .init()) async throws -> VivoPhysiologyStepResult {
        let prepared = try await prepareStep(request)
        if !prepared.canCommit { return .init(certificate: certificate(prepared, committed: false), publications: []) }
        return .init(certificate: try commitPreparedStep(transactionID: prepared.transactionID), publications: prepared.publications)
    }
    private func certificate(_ value: VivoPreparedPhysiologyStep, committed: Bool) -> VivoPhysiologyStepCertificate {
        .init(transactionID: value.transactionID,
              disposition: !committed ? .rejected : (value.candidateTimeStep == value.requestedTimeStep ? .committed : .committedWithReducedStep),
              modelFingerprint: model.fingerprint, deviceName: capabilities.deviceName, deviceRegistryID: capabilities.registryID,
              stepIndex: value.stepIndex, timeBefore: value.timeBefore,
              timeAfter: committed ? absoluteTimeSeconds : value.timeBefore,
              requestedTimeStep: value.requestedTimeStep, acceptedTimeStep: committed ? value.candidateTimeStep : nil,
              attemptCount: value.attemptCount, status: value.status, appliedDoseIdentifiers: value.appliedDoseIdentifiers)
    }

    public func snapshot() async throws -> VivoPhysiologySnapshot {
        try requireIdle()
        inFlight = true
        defer { inFlight = false }
        let step = nextStepIndex, time = absoluteTimeSeconds
        let values = try await readCurrentState()
        return .init(modelFingerprint: model.fingerprint, stepIndex: step, absoluteTimeSeconds: time,
                      pairCount: model.pairCount, environmentCount: model.environmentCount, values: values)
    }
    public func checkpoint() async throws -> VivoPhysiologyCheckpoint {
        // Hold one reservation across readback AND metadata capture. Calling the
        // public snapshot() then reading a later cursor would mix generations.
        try requireIdle()
        inFlight = true
        defer { inFlight = false }
        let step = nextStepIndex, time = absoluteTimeSeconds, cursor = doseCursor
        let values = try await readCurrentState()
        return .init(modelFingerprint: model.fingerprint, stepIndex: step, absoluteTimeSeconds: time,
                      doseCursor: cursor, pairCount: model.pairCount, environmentCount: model.environmentCount,
                      stateFP32LE: VivoLittleEndianFP32.encode(values))
    }
    public func restore(_ checkpoint: VivoPhysiologyCheckpoint) throws {
        try requireIdle()
        try checkpoint.validate()
        guard checkpoint.modelFingerprint == model.fingerprint, checkpoint.pairCount == model.pairCount,
              checkpoint.environmentCount == model.environmentCount, checkpoint.doseCursor <= model.doses.count else {
            throw VivoRuntimeError.invalidConfiguration("physiology checkpoint identity or shape mismatch")
        }
        // At an accepted boundary every dose strictly before time has been
        // consumed, and no future dose may already have advanced the cursor.
        let epsilon = boundaryTolerance(at: checkpoint.absoluteTimeSeconds)
        for (index, dose) in model.doses.enumerated() {
            if index < checkpoint.doseCursor, dose.timeSeconds > checkpoint.absoluteTimeSeconds + epsilon {
                throw VivoRuntimeError.invalidConfiguration("checkpoint cursor has consumed a future dose")
            }
            if index >= checkpoint.doseCursor, dose.timeSeconds < checkpoint.absoluteTimeSeconds - epsilon {
                throw VivoRuntimeError.invalidConfiguration("checkpoint cursor would replay a past dose")
            }
        }
        let values = try VivoLittleEndianFP32.decode(checkpoint.stateFP32LE)
        try validateRestored(values)
        let replacement = try VivoPhysiologyMetalArena(device: device, commandQueue: commandQueue,
                                                        model: model, configuration: configuration)
        try replacement.uploadCurrentState(values, commandQueue: commandQueue, label: "NumiVivo.Physiology.Restore")
        arena = replacement
        absoluteTimeSeconds = checkpoint.absoluteTimeSeconds
        nextStepIndex = checkpoint.stepIndex
        doseCursor = checkpoint.doseCursor
    }

    private func execute(transactionID: UUID, dt: Double, pre: [VivoPhysiologyStateTransformABI],
                         post: [VivoPhysiologyStateTransformABI], publications: [VivoPhysiologyPublicationRequestABI],
                         appliedDoseCount: Int) async throws -> Output {
        guard Float(dt).isFinite, Float(dt) > 0, absoluteTimeSeconds + dt > absoluteTimeSeconds,
              Float(absoluteTimeSeconds + dt).isFinite else {
            throw VivoRuntimeError.invalidConfiguration("physiology interval is not representable")
        }
        var uniforms = try VivoPhysiologyRuntimeCommandABI(model: model, stepIndex: nextStepIndex, dt: Float(dt),
                                                           absoluteTime: Float(absoluteTimeSeconds + dt), transactionID: transactionID,
                                                           preTransformCount: pre.count, postTransformCount: post.count,
                                                           publicationCount: publications.count, appliedDoseCount: appliedDoseCount,
                                                           maximumDerivative: configuration.maximumAbsoluteDerivative)
        uniforms.reservedWord0 = configuration.boundTolerance.bitPattern
        try arena.write(pre: pre); try arena.write(post: post); try arena.write(publications: publications)
        guard let command = commandQueue.makeCommandBuffer() else { throw VivoRuntimeError.commandQueueUnavailable }
        command.label = "NumiVivo.Physiology.\(nextStepIndex)"
        try encode(.physiologyClearStatus, count: 1, command: command, buffers: [(1, arena.status)], uniforms: uniforms, uniformIndex: 0)
        try encode(.physiologyPrepareTransaction, count: arena.capacities.stateElements, command: command,
                   buffers: [(0,arena.currentState),(1,arena.baseState),(2,arena.candidateState),(3,arena.stageState),(4,arena.derivativeK1)],
                   uniforms: uniforms, uniformIndex: 5)
        if !pre.isEmpty { try transforms(command, source: arena.preTransforms, destination: arena.baseState, count: pre.count, uniforms: uniforms) }
        try encode(.physiologyHeunPredict, count: arena.capacities.stateElements, command: command,
                   buffers: [(0,arena.baseState),(1,arena.stageState),(2,arena.derivativeK1),(3,arena.incidenceOffsets),(4,arena.incidence),(5,arena.clearances),(7,arena.status)],
                   uniforms: uniforms, uniformIndex: 6)
        try encode(.physiologyHeunCorrect, count: arena.capacities.stateElements, command: command,
                   buffers: [(0,arena.baseState),(1,arena.stageState),(2,arena.derivativeK1),(3,arena.candidateState),(4,arena.incidenceOffsets),(5,arena.incidence),(6,arena.clearances),(8,arena.status)],
                   uniforms: uniforms, uniformIndex: 7)
        if !post.isEmpty { try transforms(command, source: arena.postTransforms, destination: arena.candidateState, count: post.count, uniforms: uniforms) }
        try encode(.physiologyValidateCandidate, count: arena.capacities.stateElements, command: command,
                   buffers: [(0,arena.candidateState),(1,arena.bounds),(3,arena.status)], uniforms: uniforms, uniformIndex: 2)
        if !publications.isEmpty {
            try encode(.physiologyPublish, count: publications.count, command: command,
                       buffers: [(0,arena.candidateState),(1,arena.publicationRequests),(2,arena.publicationOutput),(4,arena.status)],
                       uniforms: uniforms, uniformIndex: 3)
        }
        try await complete(command)
        return .init(status: arena.runtimeStatus(), publications: try arena.publicationValues(count: publications.count))
    }
    private func transforms(_ command: MTLCommandBuffer, source: MTLBuffer, destination: MTLBuffer, count: Int,
                            uniforms: VivoPhysiologyRuntimeCommandABI) throws {
        try encode(.physiologyApplyTransforms, count: count, command: command,
                   buffers: [(0,source),(1,destination),(3,arena.status)], uniforms: uniforms, uniformIndex: 2, countIndex: 4)
    }
    private func encode(_ kernel: NumiVivoKernel, count: Int, command: MTLCommandBuffer,
                        buffers: [(Int, MTLBuffer)], uniforms: VivoPhysiologyRuntimeCommandABI,
                        uniformIndex: Int, countIndex: Int? = nil) throws {
        guard count > 0 else { return }
        guard count <= Int(UInt32.max), let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.pipelineEncodingFailed(kernel.rawValue)
        }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers { encoder.setBuffer(buffer, offset: 0, index: index) }
        var copy = uniforms
        encoder.setBytes(&copy, length: MemoryLayout<VivoPhysiologyRuntimeCommandABI>.stride, index: uniformIndex)
        if let countIndex {
            var value = UInt32(count)
            encoder.setBytes(&value, length: 4, index: countIndex)
        }
        encoder.dispatchThreads(pipeline.gridSize(for: count), threadsPerThreadgroup: pipeline.threadgroupSize(for: count, preferred: capabilities.recommendedThreadsPerThreadgroup))
    }
    private func complete(_ command: MTLCommandBuffer) async throws {
        do {
            try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
                command.addCompletedHandler { completed in
                    guard completed.status == .completed else {
                        continuation.resume(throwing: VivoRuntimeError.commandFailed(String(describing: completed.error))); return
                    }
                    continuation.resume(returning: ())
                }
                command.commit()
            }
        } catch { failure = String(describing: error); throw error }
    }
    private func readCurrentState() async throws -> [Float] {
        guard let command = commandQueue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        blit.copy(from: arena.currentState, sourceOffset: 0, to: arena.stateReadback, destinationOffset: 0, size: arena.stateReadback.length)
        blit.endEncoding()
        try await complete(command)
        return arena.readbackValues()
    }
    private func validate(_ request: VivoPhysiologyPrepareRequest) throws {
        guard request.preUpdates.count <= configuration.maximumTransformsPerStep,
              request.postUpdates.count <= configuration.maximumTransformsPerStep,
              request.publications.count <= configuration.maximumPublicationsPerStep else {
            throw VivoRuntimeError.invalidConfiguration("physiology request exceeds capacity")
        }
        for update in request.preUpdates + request.postUpdates {
            guard update.pairIndex < model.pairCount, update.environmentIndex < model.environmentCount, update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("invalid physiology update")
            }
        }
        for value in request.publications {
            guard value.pairIndex < model.pairCount, value.environmentIndex < model.environmentCount, value.flags == 0 else {
                throw VivoRuntimeError.invalidConfiguration("invalid physiology publication")
            }
        }
    }
    private func normalize(_ updates: [VivoPhysiologyStateUpdate], dt: Double) throws -> [VivoPhysiologyStateTransformABI] {
        var values: [UInt64: Accumulator] = [:]
        for update in updates {
            guard update.pairIndex < model.pairCount, update.environmentIndex < model.environmentCount, update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("invalid physiology transform")
            }
            let key = UInt64(update.pairIndex) << 32 | UInt64(update.environmentIndex)
            var value = values[key] ?? Accumulator()
            switch update.mode {
            case .replace:
                guard value.replacement == nil else { throw VivoRuntimeError.invalidConfiguration("multiple replacements for one physiology element") }
                value.replacement = update.value
            case .add: value.delta += Double(update.value)
            case .rate: value.delta += Double(update.value) * dt
            case .minimum: value.minimum = max(value.minimum ?? -.greatestFiniteMagnitude, update.value)
            case .maximum: value.maximum = min(value.maximum ?? .greatestFiniteMagnitude, update.value)
            }
            guard value.delta.isFinite, abs(value.delta) <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoRuntimeError.invalidConfiguration("physiology transform overflow")
            }
            values[key] = value
            guard values.count <= configuration.maximumTransformsPerStep else {
                throw VivoRuntimeError.invalidConfiguration("normalized physiology transform capacity exceeded")
            }
        }
        return try values.sorted { $0.key < $1.key }.map { key, value in
            let minimum = value.minimum ?? -Float.greatestFiniteMagnitude, maximum = value.maximum ?? Float.greatestFiniteMagnitude
            guard minimum <= maximum else { throw VivoRuntimeError.invalidConfiguration("conflicting physiology transform bounds") }
            var flags: UInt32 = 0
            if value.replacement != nil { flags |= VivoPhysiologyStateTransformABI.replaceFlag }
            if value.minimum != nil { flags |= VivoPhysiologyStateTransformABI.minimumFlag }
            if value.maximum != nil { flags |= VivoPhysiologyStateTransformABI.maximumFlag }
            return .init(pairIndex: UInt32(key >> 32), environmentIndex: UInt32(truncatingIfNeeded: key), flags: flags,
                          replacement: value.replacement ?? 0, additiveDelta: Float(value.delta), minimum: minimum, maximum: maximum)
        }
    }
    private func alignedTimeStep(_ requested: Double) throws -> Double {
        var value = requested
        let epsilon = boundaryTolerance(at: absoluteTimeSeconds)
        for dose in model.doses {
            for boundary in [dose.timeSeconds, dose.endTimeSeconds] where boundary > absoluteTimeSeconds + epsilon {
                value = min(value, boundary - absoluteTimeSeconds)
            }
        }
        // Do not round a truncated interval back UP across a scheduled event.
        guard value >= model.minimumTimeStepSeconds else {
            throw VivoRuntimeError.invalidConfiguration("a physiology event boundary is closer than the minimum timestep; reduce the configured minimum")
        }
        guard value > 0, absoluteTimeSeconds + value > absoluteTimeSeconds else {
            throw VivoRuntimeError.invalidConfiguration("physiology clock cannot advance")
        }
        return value
    }
    private func dosePlan(dt: Double) throws -> DosePlan {
        let epsilon = boundaryTolerance(at: absoluteTimeSeconds)
        var pre: [VivoPhysiologyStateUpdate] = [], post: [VivoPhysiologyStateUpdate] = []
        var identifiers: [String] = []
        var cursor = doseCursor
        while cursor < model.doses.count, model.doses[cursor].timeSeconds <= absoluteTimeSeconds + epsilon {
            let dose = model.doses[cursor]
            if dose.kind == .concentrationDelta {
                for environment in dose.environments {
                    guard pre.count < configuration.maximumTransformsPerStep else {
                        throw VivoRuntimeError.invalidConfiguration("expanded bolus capacity exceeded")
                    }
                    pre.append(.init(pairIndex: dose.pairIndex, environmentIndex: environment, mode: .add, value: dose.value))
                }
                identifiers.append(dose.identifier)
            }
            cursor += 1
        }
        for dose in model.doses where dose.kind == .concentrationInfusion &&
            dose.timeSeconds <= absoluteTimeSeconds + epsilon && dose.endTimeSeconds > absoluteTimeSeconds + epsilon {
            // Symmetric source kicks retain second-order forcing integration.
            // Applying the entire infusion before clearance is only first order.
            let half = dose.value * 0.5
            guard dose.value == 0 || half > 0 else { throw VivoRuntimeError.invalidConfiguration("infusion half-rate underflows FP32") }
            for environment in dose.environments {
                let update = VivoPhysiologyStateUpdate(pairIndex: dose.pairIndex, environmentIndex: environment, mode: .rate, value: half)
                guard pre.count < configuration.maximumTransformsPerStep,
                      post.count < configuration.maximumTransformsPerStep else {
                    throw VivoRuntimeError.invalidConfiguration("expanded infusion capacity exceeded")
                }
                pre.append(update); post.append(update)
            }
            identifiers.append(dose.identifier)
        }
        guard pre.count <= configuration.maximumTransformsPerStep, post.count <= configuration.maximumTransformsPerStep else {
            throw VivoRuntimeError.invalidConfiguration("expanded dose updates exceed configured transform capacity")
        }
        let end = absoluteTimeSeconds + dt
        var boundary: Double?
        for dose in model.doses {
            for candidate in [dose.timeSeconds, dose.endTimeSeconds]
                where candidate > absoluteTimeSeconds + epsilon && candidate <= end + epsilon {
                boundary = min(boundary ?? candidate, candidate)
            }
        }
        return .init(pre: pre, post: post, identifiers: Array(Set(identifiers)).sorted(), nextCursor: cursor, nextBoundary: boundary)
    }
    private func boundaryTolerance(at time: Double) -> Double { max(1e-12, abs(time) * 8 * Double.ulpOfOne) }
    private func validateRestored(_ values: [Float]) throws {
        guard values.count == model.initialState.count else { throw VivoRuntimeError.invalidConfiguration("physiology restore shape mismatch") }
        for (index, value) in values.enumerated() {
            let pair = index / Int(model.environmentCount)
            let analyte = model.analytes[pair / model.compartments.count]
            let low = Float(analyte.minimum), high = Float(min(analyte.maximum, Double(Float.greatestFiniteMagnitude)))
            guard value.isFinite, value >= low - configuration.boundTolerance, value <= high + configuration.boundTolerance else {
                throw VivoRuntimeError.invalidConfiguration("physiology restore violates analyte bounds")
            }
        }
    }
}
