import Foundation
import Metal
import NumiVivoShaders

private struct VivoPhysiologyPipelines: Sendable {
    let clearStatus: NumiVivoPipeline
    let prepare: NumiVivoPipeline
    let applyTransforms: NumiVivoPipeline
    let predict: NumiVivoPipeline
    let correct: NumiVivoPipeline
    let validate: NumiVivoPipeline
    let publish: NumiVivoPipeline

    static func load(from catalog: NumiVivoPipelineCatalog) async throws -> Self {
        try await Self(
            clearStatus: catalog.pipeline(.physiologyClearStatus),
            prepare: catalog.pipeline(.physiologyPrepareTransaction),
            applyTransforms: catalog.pipeline(.physiologyApplyTransforms),
            predict: catalog.pipeline(.physiologyHeunPredict),
            correct: catalog.pipeline(.physiologyHeunCorrect),
            validate: catalog.pipeline(.physiologyValidateCandidate),
            publish: catalog.pipeline(.physiologyPublish)
        )
    }
}

public actor VivoPhysiologyRuntime {
    private struct PendingStep {
        let prepared: VivoPreparedPhysiologyStep
        let requestedTimeStep: Double
        let nextDoseCursor: Int
    }

    private struct AttemptOutput {
        let status: VivoPhysiologyRuntimeStatus
        let publications: [Float]
    }

    private struct TransformAccumulator {
        var replacement: Float?
        var additiveDelta: Double = 0
        var minimum: Float?
        var maximum: Float?
    }

    public nonisolated let model: PreparedVivoPhysiologyModel
    public nonisolated let configuration: VivoPhysiologyRuntimeConfiguration
    public nonisolated let capabilities: VivoMetalCapabilities

    private let device: MTLDevice
    private let commandQueue: MTLCommandQueue
    private let arena: VivoPhysiologyMetalArena
    private let pipelines: VivoPhysiologyPipelines

    private var absoluteTimeSeconds: Double = 0
    private var nextStepIndex: UInt32 = 0
    private var doseCursor: Int = 0
    private var pending: PendingStep?

    public static func make(
        model: PreparedVivoPhysiologyModel,
        configuration: VivoPhysiologyRuntimeConfiguration = .init(),
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoPhysiologyRuntime {
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard let queue = device.makeCommandQueue() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.Physiology.CommandQueue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let pipelines = try await VivoPhysiologyPipelines.load(from: catalog)
        let arena = try VivoPhysiologyMetalArena(
            device: device,
            commandQueue: queue,
            model: model,
            configuration: configuration
        )
        return VivoPhysiologyRuntime(
            model: model,
            configuration: configuration,
            device: device,
            commandQueue: queue,
            arena: arena,
            pipelines: pipelines
        )
    }

    private init(
        model: PreparedVivoPhysiologyModel,
        configuration: VivoPhysiologyRuntimeConfiguration,
        device: MTLDevice,
        commandQueue: MTLCommandQueue,
        arena: VivoPhysiologyMetalArena,
        pipelines: VivoPhysiologyPipelines
    ) {
        self.model = model
        self.configuration = configuration
        self.device = device
        self.commandQueue = commandQueue
        self.arena = arena
        self.pipelines = pipelines
        self.capabilities = arena.capabilities
    }

    public func time() -> Double { absoluteTimeSeconds }
    public func stepIndex() -> UInt32 { nextStepIndex }
    public func hasPendingTransaction() -> Bool { pending != nil }

    public func prepareStep(
        _ request: VivoPhysiologyPrepareRequest = .init(),
        transactionID: UUID = UUID()
    ) async throws -> VivoPreparedPhysiologyStep {
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration("physiology runtime already has a prepared transaction")
        }
        try validate(request)

        let requested = request.timeStepSeconds ?? model.preferredTimeStepSeconds
        guard requested.isFinite,
              requested >= model.minimumTimeStepSeconds,
              requested <= model.maximumTimeStepSeconds else {
            throw VivoRuntimeError.invalidConfiguration("requested physiology time step is outside prepared-model bounds")
        }
        var dt = alignedTimeStep(requested)
        var attempt: UInt32 = 0
        var lastStatus = VivoPhysiologyRuntimeStatus(raw: .init())
        var lastBoundary: Double?

        while attempt < configuration.maximumSubsteps {
            let dosePlan = try dosePlan(time: absoluteTimeSeconds, dt: dt)
            let pre = try normalize(
                dosePlan.updates + request.preUpdates,
                dt: dt,
                capacity: configuration.maximumTransformsPerStep,
                label: "pre"
            )
            let post = try normalize(
                request.postUpdates,
                dt: dt,
                capacity: configuration.maximumTransformsPerStep,
                label: "post"
            )
            let publications = try preparePublications(request.publications)
            let output = try await executeAttempt(
                transactionID: transactionID,
                dt: dt,
                pre: pre,
                post: post,
                publications: publications,
                appliedDoseCount: dosePlan.identifiers.count
            )
            lastStatus = output.status
            lastBoundary = dosePlan.nextBoundary

            let needsReduction = output.status.flags.contains(.requiresSubstep) ||
                output.status.flags.contains(.excessiveDerivative) ||
                output.status.blocksCommit
            if !needsReduction {
                let prepared = VivoPreparedPhysiologyStep(
                    transactionID: transactionID,
                    disposition: .prepared,
                    modelFingerprint: model.fingerprint,
                    stepIndex: nextStepIndex,
                    timeBefore: absoluteTimeSeconds,
                    requestedTimeStep: requested,
                    candidateTimeStep: dt,
                    attemptCount: attempt + 1,
                    status: output.status,
                    publications: output.publications,
                    appliedDoseIdentifiers: dosePlan.identifiers,
                    nextBoundarySeconds: dosePlan.nextBoundary
                )
                pending = PendingStep(
                    prepared: prepared,
                    requestedTimeStep: requested,
                    nextDoseCursor: dosePlan.nextCursor
                )
                return prepared
            }

            let suggested = Double(output.status.suggestedTimeStep)
            let halved = dt * 0.5
            let candidate = min(
                halved,
                suggested.isFinite && suggested > 0 ? suggested : halved
            )
            let canRetry = request.permitAdaptiveReduction &&
                attempt + 1 < configuration.maximumSubsteps &&
                candidate >= model.minimumTimeStepSeconds &&
                candidate < dt
            guard canRetry else {
                return VivoPreparedPhysiologyStep(
                    transactionID: transactionID,
                    disposition: candidate < dt ? .requiresSmallerStep : .rejected,
                    modelFingerprint: model.fingerprint,
                    stepIndex: nextStepIndex,
                    timeBefore: absoluteTimeSeconds,
                    requestedTimeStep: requested,
                    candidateTimeStep: dt,
                    attemptCount: attempt + 1,
                    status: output.status,
                    publications: [],
                    appliedDoseIdentifiers: dosePlan.identifiers,
                    nextBoundarySeconds: dosePlan.nextBoundary
                )
            }
            dt = alignedTimeStep(candidate)
            attempt &+= 1
        }

        return VivoPreparedPhysiologyStep(
            transactionID: transactionID,
            disposition: .rejected,
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            timeBefore: absoluteTimeSeconds,
            requestedTimeStep: requested,
            candidateTimeStep: dt,
            attemptCount: attempt,
            status: lastStatus,
            publications: [],
            appliedDoseIdentifiers: [],
            nextBoundarySeconds: lastBoundary
        )
    }

    public func commitPreparedStep(transactionID: UUID) throws -> VivoPhysiologyStepCertificate {
        guard let pending, pending.prepared.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration("physiology transaction identifier does not match the prepared candidate")
        }
        guard pending.prepared.canCommit else {
            throw VivoRuntimeError.invalidConfiguration("physiology candidate is not eligible for commit")
        }
        arena.commit()
        let before = absoluteTimeSeconds
        absoluteTimeSeconds += pending.prepared.candidateTimeStep
        doseCursor = pending.nextDoseCursor
        nextStepIndex &+= 1
        self.pending = nil
        return VivoPhysiologyStepCertificate(
            transactionID: transactionID,
            disposition: pending.prepared.candidateTimeStep == pending.requestedTimeStep
                ? .committed
                : .committedWithReducedStep,
            modelFingerprint: model.fingerprint,
            deviceName: capabilities.deviceName,
            deviceRegistryID: capabilities.registryID,
            stepIndex: pending.prepared.stepIndex,
            timeBefore: before,
            timeAfter: absoluteTimeSeconds,
            requestedTimeStep: pending.requestedTimeStep,
            acceptedTimeStep: pending.prepared.candidateTimeStep,
            attemptCount: pending.prepared.attemptCount,
            status: pending.prepared.status,
            appliedDoseIdentifiers: pending.prepared.appliedDoseIdentifiers
        )
    }

    public func discardPreparedStep(transactionID: UUID) throws {
        guard let pending, pending.prepared.transactionID == transactionID else {
            throw VivoRuntimeError.invalidConfiguration("physiology transaction identifier does not match the prepared candidate")
        }
        self.pending = nil
    }

    public func step(_ request: VivoPhysiologyPrepareRequest = .init()) async throws -> VivoPhysiologyStepResult {
        let transactionID = UUID()
        let prepared = try await prepareStep(request, transactionID: transactionID)
        guard prepared.canCommit else {
            let certificate = VivoPhysiologyStepCertificate(
                transactionID: transactionID,
                disposition: .rejected,
                modelFingerprint: model.fingerprint,
                deviceName: capabilities.deviceName,
                deviceRegistryID: capabilities.registryID,
                stepIndex: prepared.stepIndex,
                timeBefore: prepared.timeBefore,
                timeAfter: prepared.timeBefore,
                requestedTimeStep: prepared.requestedTimeStep,
                acceptedTimeStep: nil,
                attemptCount: prepared.attemptCount,
                status: prepared.status,
                appliedDoseIdentifiers: prepared.appliedDoseIdentifiers
            )
            return VivoPhysiologyStepResult(certificate: certificate, publications: [])
        }
        let certificate = try commitPreparedStep(transactionID: transactionID)
        return VivoPhysiologyStepResult(certificate: certificate, publications: prepared.publications)
    }

    public func snapshot() async throws -> VivoPhysiologySnapshot {
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration("cannot snapshot physiology while a candidate transaction is pending")
        }
        let values = try await readCurrentState()
        return VivoPhysiologySnapshot(
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            absoluteTimeSeconds: absoluteTimeSeconds,
            pairCount: model.pairCount,
            environmentCount: model.environmentCount,
            values: values
        )
    }

    public func checkpoint() async throws -> VivoPhysiologyCheckpoint {
        let snapshot = try await snapshot()
        var data = Data(capacity: snapshot.values.count * MemoryLayout<UInt32>.stride)
        for value in snapshot.values {
            var bits = value.bitPattern.littleEndian
            withUnsafeBytes(of: &bits) { data.append(contentsOf: $0) }
        }
        return VivoPhysiologyCheckpoint(
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            absoluteTimeSeconds: absoluteTimeSeconds,
            doseCursor: doseCursor,
            pairCount: model.pairCount,
            environmentCount: model.environmentCount,
            stateFP32LE: data
        )
    }

    public func restore(_ checkpoint: VivoPhysiologyCheckpoint) throws {
        guard pending == nil else {
            throw VivoRuntimeError.invalidConfiguration("cannot restore physiology while a candidate transaction is pending")
        }
        try checkpoint.validate()
        guard checkpoint.modelFingerprint == model.fingerprint,
              checkpoint.pairCount == model.pairCount,
              checkpoint.environmentCount == model.environmentCount,
              checkpoint.doseCursor <= model.doses.count else {
            throw VivoRuntimeError.invalidConfiguration("physiology checkpoint is incompatible with the runtime model")
        }
        var values: [Float] = []
        values.reserveCapacity(model.initialState.count)
        checkpoint.stateFP32LE.withUnsafeBytes { raw in
            for offset in stride(from: 0, to: raw.count, by: 4) {
                let stored = raw.loadUnaligned(fromByteOffset: offset, as: UInt32.self)
                values.append(Float(bitPattern: UInt32(littleEndian: stored)))
            }
        }
        try validateRestored(values)
        try arena.uploadCurrentState(
            values,
            commandQueue: commandQueue,
            label: "NumiVivo.Physiology.Restore"
        )
        absoluteTimeSeconds = checkpoint.absoluteTimeSeconds
        nextStepIndex = checkpoint.stepIndex
        doseCursor = checkpoint.doseCursor
    }

    private func executeAttempt(
        transactionID: UUID,
        dt: Double,
        pre: [VivoPhysiologyStateTransformABI],
        post: [VivoPhysiologyStateTransformABI],
        publications: [VivoPhysiologyPublicationRequestABI],
        appliedDoseCount: Int
    ) async throws -> AttemptOutput {
        guard dt <= Double(Float.greatestFiniteMagnitude),
              absoluteTimeSeconds + dt <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoRuntimeError.invalidConfiguration("physiology time is not representable by the FP32 Metal ABI")
        }
        var command = try VivoPhysiologyRuntimeCommandABI(
            model: model,
            stepIndex: nextStepIndex,
            dt: Float(dt),
            absoluteTime: Float(absoluteTimeSeconds + dt),
            transactionID: transactionID,
            preTransformCount: pre.count,
            postTransformCount: post.count,
            publicationCount: publications.count,
            appliedDoseCount: appliedDoseCount,
            maximumDerivative: configuration.maximumAbsoluteDerivative
        )
        command.reservedWord0 = configuration.boundTolerance.bitPattern
        arena.write(command: command)
        try arena.write(pre: pre)
        try arena.write(post: post)
        try arena.write(publications: publications)

        guard let commandBuffer = commandQueue.makeCommandBuffer() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        commandBuffer.label = "NumiVivo.Physiology.Step.\(nextStepIndex)"
        try encodeClearStatus(commandBuffer)
        try encodePrepare(commandBuffer)
        if !pre.isEmpty {
            try encodeTransforms(commandBuffer, source: arena.preTransforms, destination: arena.baseState, count: pre.count, label: "PreTransforms")
        }
        try encodePredict(commandBuffer)
        try encodeCorrect(commandBuffer)
        if !post.isEmpty {
            try encodeTransforms(commandBuffer, source: arena.postTransforms, destination: arena.candidateState, count: post.count, label: "PostTransforms")
        }
        try encodeValidation(commandBuffer)
        if !publications.isEmpty {
            try encodePublications(commandBuffer, count: publications.count)
        }
        try await complete(commandBuffer)
        return AttemptOutput(
            status: arena.runtimeStatus(),
            publications: try arena.publicationValues(count: publications.count)
        )
    }

    private func encodeClearStatus(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(pipelines.clearStatus, count: 1, label: "NumiVivo.Physiology.ClearStatus", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.command, offset: 0, index: 0)
            encoder.setBuffer(arena.status, offset: 0, index: 1)
        }
    }

    private func encodePrepare(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(pipelines.prepare, count: arena.capacities.stateElements, label: "NumiVivo.Physiology.Prepare", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.currentState, offset: 0, index: 0)
            encoder.setBuffer(arena.baseState, offset: 0, index: 1)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 2)
            encoder.setBuffer(arena.stageState, offset: 0, index: 3)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 4)
            encoder.setBuffer(arena.command, offset: 0, index: 5)
        }
    }

    private func encodeTransforms(
        _ commandBuffer: MTLCommandBuffer,
        source: MTLBuffer,
        destination: MTLBuffer,
        count: Int,
        label: String
    ) throws {
        try encode(pipelines.applyTransforms, count: count, label: "NumiVivo.Physiology.\(label)", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(source, offset: 0, index: 0)
            encoder.setBuffer(destination, offset: 0, index: 1)
            encoder.setBuffer(arena.command, offset: 0, index: 2)
            encoder.setBuffer(arena.status, offset: 0, index: 3)
            var value = UInt32(count)
            encoder.setBytes(&value, length: MemoryLayout<UInt32>.stride, index: 4)
        }
    }

    private func encodePredict(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(pipelines.predict, count: arena.capacities.stateElements, label: "NumiVivo.Physiology.HeunPredict", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.baseState, offset: 0, index: 0)
            encoder.setBuffer(arena.stageState, offset: 0, index: 1)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 2)
            encoder.setBuffer(arena.incidenceOffsets, offset: 0, index: 3)
            encoder.setBuffer(arena.incidence, offset: 0, index: 4)
            encoder.setBuffer(arena.clearances, offset: 0, index: 5)
            encoder.setBuffer(arena.command, offset: 0, index: 6)
            encoder.setBuffer(arena.status, offset: 0, index: 7)
        }
    }

    private func encodeCorrect(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(pipelines.correct, count: arena.capacities.stateElements, label: "NumiVivo.Physiology.HeunCorrect", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.baseState, offset: 0, index: 0)
            encoder.setBuffer(arena.stageState, offset: 0, index: 1)
            encoder.setBuffer(arena.derivativeK1, offset: 0, index: 2)
            encoder.setBuffer(arena.candidateState, offset: 0, index: 3)
            encoder.setBuffer(arena.incidenceOffsets, offset: 0, index: 4)
            encoder.setBuffer(arena.incidence, offset: 0, index: 5)
            encoder.setBuffer(arena.clearances, offset: 0, index: 6)
            encoder.setBuffer(arena.command, offset: 0, index: 7)
            encoder.setBuffer(arena.status, offset: 0, index: 8)
        }
    }

    private func encodeValidation(_ commandBuffer: MTLCommandBuffer) throws {
        try encode(pipelines.validate, count: arena.capacities.stateElements, label: "NumiVivo.Physiology.Validate", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.bounds, offset: 0, index: 1)
            encoder.setBuffer(arena.command, offset: 0, index: 2)
            encoder.setBuffer(arena.status, offset: 0, index: 3)
        }
    }

    private func encodePublications(_ commandBuffer: MTLCommandBuffer, count: Int) throws {
        try encode(pipelines.publish, count: count, label: "NumiVivo.Physiology.Publish", commandBuffer: commandBuffer) { encoder in
            encoder.setBuffer(arena.candidateState, offset: 0, index: 0)
            encoder.setBuffer(arena.publicationRequests, offset: 0, index: 1)
            encoder.setBuffer(arena.publicationOutput, offset: 0, index: 2)
            encoder.setBuffer(arena.command, offset: 0, index: 3)
            encoder.setBuffer(arena.status, offset: 0, index: 4)
        }
    }

    private func encode(
        _ pipeline: NumiVivoPipeline,
        count: Int,
        label: String,
        commandBuffer: MTLCommandBuffer,
        bindings: (MTLComputeCommandEncoder) -> Void
    ) throws {
        guard let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandFailed("could not create compute encoder for \(label)")
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline.state)
        bindings(encoder)
        encoder.dispatchThreads(
            pipeline.gridSize(for: count),
            threadsPerThreadgroup: pipeline.threadgroupSize(for: count, preferred: capabilities.recommendedThreadsPerThreadgroup)
        )
        encoder.endEncoding()
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

    private func readCurrentState() async throws -> [Float] {
        guard let command = commandQueue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Physiology.Readback"
        blit.copy(
            from: arena.currentState,
            sourceOffset: 0,
            to: arena.stateReadback,
            destinationOffset: 0,
            size: arena.stateReadback.length
        )
        blit.endEncoding()
        try await complete(command)
        return arena.readbackValues()
    }

    private func validate(_ request: VivoPhysiologyPrepareRequest) throws {
        guard request.preUpdates.count <= configuration.maximumTransformsPerStep,
              request.postUpdates.count <= configuration.maximumTransformsPerStep,
              request.publications.count <= configuration.maximumPublicationsPerStep else {
            throw VivoRuntimeError.invalidConfiguration("physiology request exceeds configured transform or publication capacity")
        }
        for update in request.preUpdates + request.postUpdates {
            guard update.pairIndex < model.pairCount,
                  update.environmentIndex < model.environmentCount,
                  update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("physiology state update is out of bounds or non-finite")
            }
        }
        for request in request.publications {
            guard request.pairIndex < model.pairCount,
                  request.environmentIndex < model.environmentCount else {
                throw VivoRuntimeError.invalidConfiguration("physiology publication request is out of bounds")
            }
        }
    }

    private func preparePublications(
        _ requests: [VivoPhysiologyPublicationRequest]
    ) throws -> [VivoPhysiologyPublicationRequestABI] {
        try requests.enumerated().map { index, request in
            guard index <= Int(UInt32.max) else {
                throw VivoRuntimeError.invalidConfiguration("physiology publication index exceeds UInt32")
            }
            return VivoPhysiologyPublicationRequestABI(
                pairIndex: request.pairIndex,
                environmentIndex: request.environmentIndex,
                outputIndex: UInt32(index),
                flags: request.flags
            )
        }
    }

    private func normalize(
        _ updates: [VivoPhysiologyStateUpdate],
        dt: Double,
        capacity: Int,
        label: String
    ) throws -> [VivoPhysiologyStateTransformABI] {
        var accumulated: [UInt64: TransformAccumulator] = [:]
        accumulated.reserveCapacity(updates.count)
        for update in updates {
            guard update.pairIndex < model.pairCount,
                  update.environmentIndex < model.environmentCount,
                  update.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration("physiology \(label) update is invalid")
            }
            let key = UInt64(update.pairIndex) << 32 | UInt64(update.environmentIndex)
            var value = accumulated[key] ?? TransformAccumulator()
            switch update.mode {
            case .replace:
                guard value.replacement == nil else {
                    throw VivoRuntimeError.invalidConfiguration("physiology \(label) updates contain multiple replacements for one state element")
                }
                value.replacement = update.value
            case .add:
                value.additiveDelta += Double(update.value)
            case .rate:
                value.additiveDelta += Double(update.value) * dt
            case .minimum:
                value.minimum = max(value.minimum ?? -Float.greatestFiniteMagnitude, update.value)
            case .maximum:
                value.maximum = min(value.maximum ?? Float.greatestFiniteMagnitude, update.value)
            }
            guard value.additiveDelta.isFinite,
                  abs(value.additiveDelta) <= Double(Float.greatestFiniteMagnitude) else {
                throw VivoRuntimeError.invalidConfiguration("physiology \(label) additive transform overflow")
            }
            accumulated[key] = value
        }
        guard accumulated.count <= capacity else {
            throw VivoRuntimeError.invalidConfiguration("normalized physiology \(label) transform count exceeds capacity")
        }

        return try accumulated.sorted(by: { $0.key < $1.key }).map { key, value in
            let pairIndex = UInt32(key >> 32)
            let environmentIndex = UInt32(key & 0xffff_ffff)
            let minimum = value.minimum ?? -Float.greatestFiniteMagnitude
            let maximum = value.maximum ?? Float.greatestFiniteMagnitude
            guard minimum <= maximum else {
                throw VivoRuntimeError.invalidConfiguration("physiology \(label) minimum exceeds maximum")
            }
            var flags: UInt32 = 0
            if value.replacement != nil { flags |= VivoPhysiologyStateTransformABI.replaceFlag }
            if value.minimum != nil { flags |= VivoPhysiologyStateTransformABI.minimumFlag }
            if value.maximum != nil { flags |= VivoPhysiologyStateTransformABI.maximumFlag }
            return VivoPhysiologyStateTransformABI(
                pairIndex: pairIndex,
                environmentIndex: environmentIndex,
                flags: flags,
                replacement: value.replacement ?? 0,
                additiveDelta: Float(value.additiveDelta),
                minimum: minimum,
                maximum: maximum
            )
        }
    }

    private func alignedTimeStep(_ requested: Double) -> Double {
        var value = min(max(requested, model.minimumTimeStepSeconds), model.maximumTimeStepSeconds)
        let epsilon = boundaryTolerance(at: absoluteTimeSeconds)
        let end = absoluteTimeSeconds + value
        for dose in model.doses {
            for boundary in [dose.timeSeconds, dose.endTimeSeconds] {
                if boundary > absoluteTimeSeconds + epsilon,
                   boundary < end - epsilon {
                    value = min(value, boundary - absoluteTimeSeconds)
                }
            }
        }
        return max(value, model.minimumTimeStepSeconds)
    }

    private func dosePlan(
        time: Double,
        dt: Double
    ) throws -> (updates: [VivoPhysiologyStateUpdate], identifiers: [String], nextCursor: Int, nextBoundary: Double?) {
        let epsilon = boundaryTolerance(at: time)
        var updates: [VivoPhysiologyStateUpdate] = []
        var identifiers: [String] = []
        var nextCursor = doseCursor

        while nextCursor < model.doses.count,
              model.doses[nextCursor].timeSeconds <= time + epsilon {
            let dose = model.doses[nextCursor]
            if dose.kind == .concentrationDelta {
                for environment in dose.environments {
                    updates.append(.init(
                        pairIndex: dose.pairIndex,
                        environmentIndex: environment,
                        mode: .add,
                        value: dose.value
                    ))
                }
                identifiers.append(dose.identifier)
            }
            nextCursor += 1
        }

        for dose in model.doses where dose.kind == .concentrationInfusion {
            if dose.timeSeconds <= time + epsilon,
               dose.endTimeSeconds > time + epsilon {
                for environment in dose.environments {
                    updates.append(.init(
                        pairIndex: dose.pairIndex,
                        environmentIndex: environment,
                        mode: .rate,
                        value: dose.value
                    ))
                }
                identifiers.append(dose.identifier)
            }
        }

        let end = time + dt
        var nextBoundary: Double?
        for dose in model.doses {
            for boundary in [dose.timeSeconds, dose.endTimeSeconds] {
                if boundary > time + epsilon, boundary <= end + epsilon {
                    nextBoundary = min(nextBoundary ?? boundary, boundary)
                }
            }
        }
        return (updates, Array(Set(identifiers)).sorted(), nextCursor, nextBoundary)
    }

    private func boundaryTolerance(at time: Double) -> Double {
        max(1e-12, abs(time) * 8 * Double.ulpOfOne)
    }

    private func validateRestored(_ values: [Float]) throws {
        guard values.count == model.initialState.count,
              values.allSatisfy(\.isFinite) else {
            throw VivoRuntimeError.invalidConfiguration("restored physiology values are invalid")
        }
        let environments = Int(model.environmentCount)
        for pair in 0..<Int(model.pairCount) {
            let analyteIndex = pair / model.compartments.count
            let analyte = model.analytes[analyteIndex]
            for environment in 0..<environments {
                let value = Double(values[pair * environments + environment])
                guard value >= analyte.minimum - Double(configuration.boundTolerance),
                      value <= analyte.maximum + Double(configuration.boundTolerance) else {
                    throw VivoRuntimeError.invalidConfiguration("restored physiology state violates analyte bounds")
                }
            }
        }
    }
}
