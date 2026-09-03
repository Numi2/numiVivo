@preconcurrency import Metal
import CryptoKit
import Foundation

public actor NumiVivoTransactionalRuntime {
    private struct StepPlan: Sendable {
        let substeps: Int
        let substepDeltaTime: Float
        let stableLimit: Float
    }

    private struct HashSegment: @unchecked Sendable {
        let buffer: MTLBuffer
        let byteCount: Int
        let chunkCount: Int
        let outputWordOffset: Int
        let domain: UInt32
    }

    private struct HashPlan: @unchecked Sendable {
        let buffer: MTLBuffer
        let segments: [HashSegment]
        let wordCount: Int
        let chunkBytes: Int
    }

    private struct PendingTransaction: Sendable {
        let epoch: VivoTransactionEpoch
        let receipt: VivoPrepareReceipt
        let finalUniforms: NVivoStepUniformsSwift
        let requestedDeltaTime: Float
        let executedSubsteps: Int
        let projectedTime: Double
        let projectedDelayWriteSlot: Int
        let currentStateVersion: UInt32
        let diagnostics: VivoStepDiagnostics
        let events: [VivoEvent]
    }

    private struct CommittedTransaction: Sendable {
        let epochID: UUID
        let receipt: VivoPrepareReceipt
        let result: VivoStepResult
    }

    private let program: VivoMetalProgram
    private let participantIdentifier: String
    private let hashChunkBytes: Int
    private var logicalStep: UInt64 = 0
    private var absoluteTime: Double = 0
    private var delayWriteSlot: Int = 0
    private var terminalStatus: VivoStepStatus?
    private var pending: PendingTransaction?
    private var lastCommitted: CommittedTransaction?

    public nonisolated var participantID: String { participantIdentifier }
    public nonisolated var programPack: VivoProgramPack { program.pack }
    public nonisolated var deviceInfo: VivoMetalDeviceInfo { program.deviceInfo }
    public nonisolated var memoryReport: VivoRuntimeMemoryReport { program.memoryReport }
    public nonisolated var configuration: VivoRuntimeConfiguration { program.configuration }

    public init(
        programPack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        participantID: String = "NumiVivo",
        hashChunkBytes: Int = 4_096,
        device: MTLDevice? = nil
    ) throws {
        guard !participantID.isEmpty else {
            throw VivoRuntimeError.invalidConfiguration("Transactional participantID cannot be empty.")
        }
        guard hashChunkBytes >= 64,
              hashChunkBytes <= Int(UInt32.max) else {
            throw VivoRuntimeError.invalidConfiguration(
                "hashChunkBytes must be in 64...UInt32.max."
            )
        }
        guard MemoryLayout<NVivoHashUniformsSwift>.stride == 32 else {
            throw VivoRuntimeError.invalidConfiguration(
                "Swift compiler layout does not match NumiVivo state-attestation ABI v1."
            )
        }
        self.program = try VivoMetalProgram(
            pack: programPack,
            configuration: configuration,
            device: device
        )
        self.participantIdentifier = participantID
        self.hashChunkBytes = hashChunkBytes
    }

    public init(
        sourceJSON: Data,
        compilerConfiguration: VivoCompilerConfiguration = .init(),
        runtimeConfiguration: VivoRuntimeConfiguration,
        participantID: String = "NumiVivo",
        hashChunkBytes: Int = 4_096,
        device: MTLDevice? = nil
    ) throws {
        let output = try VivoCompiler.compile(
            json: sourceJSON,
            configuration: compilerConfiguration
        )
        try self.init(
            programPack: output.programPack,
            configuration: runtimeConfiguration,
            participantID: participantID,
            hashChunkBytes: hashChunkBytes,
            device: device
        )
    }

    public func snapshot() -> VivoTransactionalRuntimeSnapshot {
        let publication = program.readPublication()
        return VivoTransactionalRuntimeSnapshot(
            committedLogicalStep: logicalStep,
            committedTime: absoluteTime,
            stateVersion: publication.stateVersion,
            pendingEpochID: pending?.epoch.id,
            terminalStatus: terminalStatus
        )
    }

    public func prepare(
        epoch: VivoTransactionEpoch,
        inputs: [VivoInputUpdate] = []
    ) async throws -> VivoPreparedStep {
        if let pending {
            throw VivoTransactionalRuntimeError.alreadyPrepared(pending.epoch.id)
        }
        try validate(epoch: epoch)
        if let terminalStatus {
            return try terminalPreparedStep(epoch: epoch, status: terminalStatus)
        }

        guard epoch.deltaTime <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoTransactionalRuntimeError.invalidEpoch(
                "deltaTime exceeds the authoritative FP32 execution range"
            )
        }
        let requestedDeltaTime = Float(epoch.deltaTime)
        let plan = try makeStepPlan(deltaTime: requestedDeltaTime)
        guard plan.substeps <= Int(UInt32.max) else {
            throw VivoRuntimeError.resourceLimit(
                "Internal substep count exceeds the Metal transaction ABI."
            )
        }
        let rawInputs = try resolveInputs(inputs)
        let hashPlan = try makeHashPlan(logicalStep: epoch.sequence)

        guard let commandBuffer = program.commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding(
                "Unable to create a transactional prepare command encoder."
            )
        }
        commandBuffer.label = "\(program.configuration.label).prepare.\(epoch.sequence)"
        encoder.label = "\(program.configuration.label).prepare.compute"
        program.declareArgumentResources(on: encoder)
        encoder.useResource(hashPlan.buffer, usage: [.write])

        var initialUniforms = program.makeUniforms(
            deltaTime: plan.substepDeltaTime,
            absoluteTime: Float(epoch.simulationTime),
            logicalStep: epoch.sequence,
            substepIndex: 0,
            delayWriteSlot: UInt32(delayWriteSlot)
        )
        encode(
            kernel: .prepareStep,
            grid: program.maximumPreparationGrid,
            uniforms: &initialUniforms,
            encoder: encoder
        )
        encoder.memoryBarrier(scope: .buffers)

        var retainedInputBuffer: MTLBuffer?
        if !rawInputs.isEmpty {
            retainedInputBuffer = rawInputs.withUnsafeBufferPointer { pointer in
                guard let baseAddress = pointer.baseAddress else { return nil }
                return program.device.makeBuffer(
                    bytes: baseAddress,
                    length: pointer.count * MemoryLayout<NVivoInputUpdateSwift>.stride,
                    options: [.storageModeShared]
                )
            }
            guard let retainedInputBuffer else {
                encoder.endEncoding()
                throw VivoRuntimeError.allocation(
                    "Unable to allocate transactional input updates."
                )
            }
            retainedInputBuffer.label = "\(program.configuration.label).prepare.inputs"
            let pipeline = program.kernels[.stageInputUpdates]
            encoder.setComputePipelineState(pipeline.state)
            encoder.setBuffer(program.buffers.argumentBuffer, offset: 0, index: 0)
            encoder.setBuffer(retainedInputBuffer, offset: 0, index: 1)
            var count = UInt32(rawInputs.count)
            encoder.setBytes(&count, length: MemoryLayout<UInt32>.stride, index: 2)
            encoder.setBytes(
                &initialUniforms,
                length: MemoryLayout<NVivoStepUniformsSwift>.stride,
                index: 3
            )
            encoder.dispatchThreads(
                MTLSize(width: rawInputs.count, height: 1, depth: 1),
                threadsPerThreadgroup: pipeline.threadsPerThreadgroup()
            )
            encoder.memoryBarrier(scope: .buffers)
        }

        var finalUniforms = initialUniforms
        var projectedDelayWriteSlot = delayWriteSlot
        for substep in 0..<plan.substeps {
            if program.delaySlotCount > 0 {
                projectedDelayWriteSlot = (delayWriteSlot + substep) % program.delaySlotCount
            } else {
                projectedDelayWriteSlot = 0
            }
            let substepEndTime = epoch.simulationTime +
                Double(plan.substepDeltaTime) * Double(substep + 1)
            guard substepEndTime.isFinite,
                  substepEndTime <= Double(Float.greatestFiniteMagnitude) else {
                encoder.endEncoding()
                throw VivoTransactionalRuntimeError.invalidEpoch(
                    "projected simulation time exceeds the authoritative FP32 execution range"
                )
            }
            var uniforms = program.makeUniforms(
                deltaTime: plan.substepDeltaTime,
                absoluteTime: Float(substepEndTime),
                logicalStep: epoch.sequence,
                substepIndex: UInt32(substep),
                delayWriteSlot: UInt32(projectedDelayWriteSlot)
            )

            for cohort in program.pack.cohorts where cohort.reactionCount > 0 {
                var cohortUniforms = NVivoCohortUniformsSwift(
                    reactionOffset: cohort.reactionOffset,
                    reactionCount: cohort.reactionCount,
                    dispatchCellCount: UInt32(program.configuration.activeCellCount)
                )
                let grid = try checkedGrid(
                    Int(cohort.reactionCount),
                    program.configuration.activeCellCount,
                    label: "transactional reaction cohort"
                )
                let pipeline = program.kernels[.evaluateReactionCohort]
                encoder.setComputePipelineState(pipeline.state)
                encoder.setBuffer(program.buffers.argumentBuffer, offset: 0, index: 0)
                encoder.setBytes(
                    &uniforms,
                    length: MemoryLayout<NVivoStepUniformsSwift>.stride,
                    index: 1
                )
                encoder.setBytes(
                    &cohortUniforms,
                    length: MemoryLayout<NVivoCohortUniformsSwift>.stride,
                    index: 2
                )
                encoder.dispatchThreads(
                    MTLSize(width: grid, height: 1, depth: 1),
                    threadsPerThreadgroup: pipeline.threadsPerThreadgroup(
                        preferred: Int(cohort.preferredThreads)
                    )
                )
            }
            if !program.pack.cohorts.isEmpty {
                encoder.memoryBarrier(scope: .buffers)
            }

            let incidenceGrid = try checkedGrid(
                program.pack.runtimeContract.speciesCount,
                program.configuration.activeCellCount,
                label: "transactional incidence"
            )
            encode(
                kernel: .applyIncidence,
                grid: incidenceGrid,
                uniforms: &uniforms,
                encoder: encoder
            )
            encoder.memoryBarrier(scope: .buffers)

            if program.pack.runtimeContract.ruleCount > 0 {
                encode(
                    kernel: .applyRules,
                    grid: program.configuration.activeCellCount,
                    uniforms: &uniforms,
                    encoder: encoder
                )
                encoder.memoryBarrier(scope: .buffers)
            }

            encode(
                kernel: .validateState,
                grid: incidenceGrid,
                uniforms: &uniforms,
                encoder: encoder
            )
            encoder.memoryBarrier(scope: .buffers)

            if program.pack.runtimeContract.monitorCount > 0 {
                encode(
                    kernel: .evaluateMonitors,
                    grid: program.configuration.activeCellCount,
                    uniforms: &uniforms,
                    encoder: encoder
                )
                encoder.memoryBarrier(scope: .buffers)
            }
            finalUniforms = uniforms
        }

        encode(
            kernel: .finalizePrepare,
            grid: 1,
            uniforms: &finalUniforms,
            encoder: encoder
        )
        encoder.memoryBarrier(scope: .buffers)
        encodeHashPlan(
            hashPlan,
            logicalStep: epoch.sequence,
            encoder: encoder
        )
        encoder.endEncoding()

        try await awaitCompletion(commandBuffer)
        _ = retainedInputBuffer

        let publication = program.readPublication()
        let diagnostics = program.readDiagnostics().publicValue
        let events = normalize(events: program.readEvents(count: Int(publication.eventCount)))
        let vote = try prepareVote(publicationStatus: publication.status)
        let preparedFingerprint = try attest(
            hashPlan: hashPlan,
            epoch: epoch,
            publication: publication,
            diagnostics: diagnostics,
            events: events,
            vote: vote
        )
        let token = vote == .prepared ? randomToken(byteCount: 32) : Data()
        let requestedMaximumDeltaTime = vote == .substepRequired
            ? Double(plan.stableLimit)
            : nil
        let receipt = VivoPrepareReceipt(
            participantID: participantIdentifier,
            epochID: epoch.id,
            vote: vote,
            token: token,
            preparedStateSHA256: preparedFingerprint,
            diagnosticFlags: diagnostics.flags,
            requestedMaximumDeltaTime: requestedMaximumDeltaTime,
            metadata: [
                "programContentFingerprint": program.pack.header.contentFingerprint,
                "projectedLogicalStep": String(epoch.sequence),
                "projectedTime": String(epoch.simulationTime + epoch.deltaTime),
                "executedSubsteps": String(plan.substeps),
                "stateVersion": String(publication.stateVersion),
                "deviceRegistryID": String(program.deviceInfo.registryID)
            ]
        )

        let projectedTime = epoch.simulationTime + epoch.deltaTime
        let preparedStep = VivoPreparedStep(
            receipt: receipt,
            requestedDeltaTime: requestedDeltaTime,
            executedSubsteps: plan.substeps,
            projectedLogicalStep: epoch.sequence,
            projectedTime: projectedTime,
            currentStateVersion: publication.stateVersion,
            diagnostics: diagnostics,
            events: vote == .prepared ? events : []
        )

        if vote == .prepared {
            let nextDelaySlot = program.delaySlotCount > 0
                ? (delayWriteSlot + plan.substeps) % program.delaySlotCount
                : 0
            pending = PendingTransaction(
                epoch: epoch,
                receipt: receipt,
                finalUniforms: finalUniforms,
                requestedDeltaTime: requestedDeltaTime,
                executedSubsteps: plan.substeps,
                projectedTime: projectedTime,
                projectedDelayWriteSlot: nextDelaySlot,
                currentStateVersion: publication.stateVersion,
                diagnostics: diagnostics,
                events: events
            )
        } else if vote == .reversibleShutdown {
            terminalStatus = .reversibleShutdown
        } else if vote == .permanentShutdown {
            terminalStatus = .permanentShutdown
        }
        return preparedStep
    }

    @discardableResult
    public func commit(
        epoch: VivoTransactionEpoch,
        receipt: VivoPrepareReceipt
    ) async throws -> VivoStepResult {
        if let committed = lastCommitted,
           committed.epochID == epoch.id,
           constantTimeEqual(committed.receipt.token, receipt.token),
           committed.receipt.preparedStateSHA256 == receipt.preparedStateSHA256 {
            return committed.result
        }
        guard let pending else {
            throw VivoTransactionalRuntimeError.noPreparedTransaction
        }
        try validate(receipt: receipt, epoch: epoch, pending: pending)

        guard let commandBuffer = program.commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding(
                "Unable to create transactional commit command encoder."
            )
        }
        commandBuffer.label = "\(program.configuration.label).commit.\(epoch.sequence)"
        encoder.label = "\(program.configuration.label).commit.compute"
        program.declareArgumentResources(on: encoder)
        var uniforms = pending.finalUniforms
        encode(
            kernel: .commitPrepared,
            grid: program.maximumCommitGrid,
            uniforms: &uniforms,
            encoder: encoder
        )
        encoder.endEncoding()
        try await awaitCompletion(commandBuffer)

        let publication = program.readPublication()
        guard publication.status == VivoStepStatus.committed.rawValue else {
            throw VivoTransactionalRuntimeError.invalidPreparedPublication(
                publication.status
            )
        }
        let publishedStep = UInt64(
            nvivoLow: publication.committedStepLow,
            high: publication.committedStepHigh
        )
        guard publishedStep == epoch.sequence else {
            throw VivoTransactionalRuntimeError.sequenceMismatch(
                expected: epoch.sequence,
                actual: publishedStep
            )
        }

        logicalStep = epoch.sequence
        absoluteTime = pending.projectedTime
        delayWriteSlot = pending.projectedDelayWriteSlot
        let result = VivoStepResult(
            status: .committed,
            requestedDeltaTime: pending.requestedDeltaTime,
            executedSubsteps: pending.executedSubsteps,
            committedLogicalStep: publishedStep,
            committedTime: absoluteTime,
            stateVersion: publication.stateVersion,
            diagnostics: pending.diagnostics,
            events: pending.events
        )
        let committedRecord = CommittedTransaction(
            epochID: epoch.id,
            receipt: pending.receipt,
            result: result
        )
        self.pending = nil
        self.lastCommitted = committedRecord
        return result
    }

    public func rollback(
        epoch: VivoTransactionEpoch,
        receipt: VivoPrepareReceipt
    ) async throws {
        if let committed = lastCommitted,
           committed.epochID == epoch.id,
           constantTimeEqual(committed.receipt.token, receipt.token) {
            return
        }
        guard let pending else {
            throw VivoTransactionalRuntimeError.noPreparedTransaction
        }
        try validate(receipt: receipt, epoch: epoch, pending: pending)

        guard let commandBuffer = program.commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding(
                "Unable to create transactional rollback command encoder."
            )
        }
        commandBuffer.label = "\(program.configuration.label).rollback.\(epoch.sequence)"
        encoder.label = "\(program.configuration.label).rollback.compute"
        program.declareArgumentResources(on: encoder)
        let pipeline = program.kernels[.rollbackPrepared]
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(program.buffers.argumentBuffer, offset: 0, index: 0)
        encoder.dispatchThreads(
            MTLSize(width: 1, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: 1, height: 1, depth: 1)
        )
        encoder.endEncoding()
        try await awaitCompletion(commandBuffer)
        self.pending = nil
    }

    public nonisolated func participant(
        inputProvider: @escaping @Sendable (VivoTransactionEpoch) async throws -> [VivoInputUpdate] = { _ in [] }
    ) -> VivoTransactionalParticipant {
        let runtime = self
        let identifier = participantIdentifier
        return VivoTransactionalParticipant(
            id: identifier,
            prepare: { epoch in
                let inputs = try await inputProvider(epoch)
                return try await runtime.prepare(epoch: epoch, inputs: inputs).receipt
            },
            commit: { epoch, receipt in
                _ = try await runtime.commit(epoch: epoch, receipt: receipt)
            },
            rollback: { epoch, receipt in
                try? await runtime.rollback(epoch: epoch, receipt: receipt)
            }
        )
    }

    public func readState(
        species identifiers: [String],
        cells requestedCells: Range<Int>? = nil
    ) async throws -> [VivoStateSlice] {
        guard pending == nil else {
            throw VivoTransactionalRuntimeError.alreadyPrepared(pending!.epoch.id)
        }
        let speciesByID = program.pack.speciesByID
        let cells = requestedCells ?? 0..<program.configuration.activeCellCount
        guard cells.lowerBound >= 0,
              cells.upperBound <= program.configuration.activeCellCount,
              !cells.isEmpty else {
            throw VivoRuntimeError.invalidConfiguration(
                "Requested cell range is outside active cells."
            )
        }
        let selected = try identifiers.map { identifier -> VivoProgramPack.Species in
            guard let species = speciesByID[identifier] else {
                throw VivoRuntimeError.unknownSignal(identifier)
            }
            return species
        }
        guard !selected.isEmpty else { return [] }

        let valuesPerSpecies = cells.count
        let bytesPerSpecies = try checkedProduct(
            valuesPerSpecies,
            MemoryLayout<Float>.stride,
            label: "transactional state readback"
        )
        let totalBytes = try checkedProduct(
            bytesPerSpecies,
            selected.count,
            label: "transactional state readback"
        )
        guard totalBytes <= program.device.maxBufferLength else {
            throw VivoRuntimeError.resourceLimit(
                "Transactional state readback exceeds the device maximum buffer length."
            )
        }
        guard let readback = program.device.makeBuffer(
            length: totalBytes,
            options: [.storageModeShared]
        ), let commandBuffer = program.commandQueue.makeCommandBuffer(),
           let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.allocation(
                "Unable to allocate transactional state readback resources."
            )
        }
        readback.label = "\(program.configuration.label).transactionalStateReadback"
        commandBuffer.label = "\(program.configuration.label).transactionalReadState"

        for (outputIndex, species) in selected.enumerated() {
            let sourceElement = try checkedProduct(
                Int(species.index),
                program.configuration.cellCapacity,
                label: "transactional state source"
            ) + cells.lowerBound
            let sourceOffset = try checkedProduct(
                sourceElement,
                MemoryLayout<Float>.stride,
                label: "transactional state source bytes"
            )
            blit.copy(
                from: program.buffers.committedState,
                sourceOffset: sourceOffset,
                to: readback,
                destinationOffset: outputIndex * bytesPerSpecies,
                size: bytesPerSpecies
            )
        }
        blit.endEncoding()
        try await awaitCompletion(commandBuffer)

        let pointer = readback.contents().bindMemory(
            to: Float.self,
            capacity: selected.count * valuesPerSpecies
        )
        return selected.enumerated().map { index, species in
            let start = index * valuesPerSpecies
            return VivoStateSlice(
                species: species.id,
                values: Array(
                    UnsafeBufferPointer(
                        start: pointer + start,
                        count: valuesPerSpecies
                    )
                )
            )
        }
    }

    public func resumeAfterReversibleShutdown() throws {
        guard pending == nil else {
            throw VivoTransactionalRuntimeError.alreadyPrepared(pending!.epoch.id)
        }
        guard terminalStatus == .reversibleShutdown else {
            if terminalStatus == .permanentShutdown {
                throw VivoRuntimeError.runtimeShutdown(.permanentShutdown)
            }
            return
        }
        terminalStatus = nil
        let pointer = program.buffers.publication.contents().bindMemory(
            to: NVivoPublicationSwift.self,
            capacity: 1
        )
        pointer.pointee.status = VivoStepStatus.committed.rawValue
        pointer.pointee.shutdownState = 0
        pointer.pointee.diagnosticFlags = 0
    }

    private func validate(epoch: VivoTransactionEpoch) throws {
        guard epoch.simulationTime.isFinite,
              epoch.deltaTime.isFinite,
              epoch.deltaTime > 0 else {
            throw VivoTransactionalRuntimeError.invalidEpoch(
                "simulationTime and deltaTime must be finite and deltaTime must be positive"
            )
        }
        guard epoch.programContentFingerprint.lowercased() ==
              program.pack.header.contentFingerprint.lowercased() else {
            throw VivoTransactionalRuntimeError.programFingerprintMismatch(
                expected: program.pack.header.contentFingerprint,
                actual: epoch.programContentFingerprint
            )
        }
        let next = logicalStep.addingReportingOverflow(1)
        guard !next.overflow else {
            throw VivoRuntimeError.resourceLimit("Logical transaction sequence overflow.")
        }
        guard epoch.sequence == next.partialValue else {
            throw VivoTransactionalRuntimeError.sequenceMismatch(
                expected: next.partialValue,
                actual: epoch.sequence
            )
        }
        let tolerance = max(1.0e-12, abs(absoluteTime) * 1.0e-9)
        guard abs(epoch.simulationTime - absoluteTime) <= tolerance else {
            throw VivoTransactionalRuntimeError.timeMismatch(
                expected: absoluteTime,
                actual: epoch.simulationTime
            )
        }
        guard (epoch.simulationTime + epoch.deltaTime).isFinite else {
            throw VivoTransactionalRuntimeError.invalidEpoch(
                "simulationTime + deltaTime overflows"
            )
        }
    }

    private func validate(
        receipt: VivoPrepareReceipt,
        epoch: VivoTransactionEpoch,
        pending: PendingTransaction
    ) throws {
        guard pending.epoch.id == epoch.id else {
            throw VivoTransactionalRuntimeError.epochMismatch(
                expected: pending.epoch.id,
                actual: epoch.id
            )
        }
        guard receipt.epochID == epoch.id else {
            throw VivoTransactionalRuntimeError.epochMismatch(
                expected: epoch.id,
                actual: receipt.epochID
            )
        }
        guard receipt.participantID == participantIdentifier else {
            throw VivoTransactionalRuntimeError.receiptParticipantMismatch(
                expected: participantIdentifier,
                actual: receipt.participantID
            )
        }
        guard receipt.vote == .prepared else {
            throw VivoTransactionalRuntimeError.invalidPreparedPublication(
                program.readPublication().status
            )
        }
        guard constantTimeEqual(pending.receipt.token, receipt.token) else {
            throw VivoTransactionalRuntimeError.receiptTokenMismatch
        }
        guard pending.receipt.preparedStateSHA256 == receipt.preparedStateSHA256 else {
            throw VivoTransactionalRuntimeError.receiptFingerprintMismatch
        }
    }

    private func terminalPreparedStep(
        epoch: VivoTransactionEpoch,
        status: VivoStepStatus
    ) throws -> VivoPreparedStep {
        let vote: VivoPrepareVote
        switch status {
        case .reversibleShutdown:
            vote = .reversibleShutdown
        case .permanentShutdown:
            vote = .permanentShutdown
        default:
            throw VivoRuntimeError.runtimeShutdown(status)
        }
        let publication = program.readPublication()
        var hasher = SHA256()
        hasher.update(data: Data("NumiVivoTerminalVote/v1\0".utf8))
        hasher.update(data: Data(program.pack.header.contentFingerprint.utf8))
        hasher.update(data: try canonicalEpoch(epoch))
        var terminalRaw = status.rawValue.littleEndian
        Swift.withUnsafeBytes(of: &terminalRaw) {
            hasher.update(data: Data($0))
        }
        let fingerprint = hex(hasher.finalize())
        let diagnostics = VivoStepDiagnostics(
            flags: status == .permanentShutdown
                ? VivoDiagnosticFlags.permanentShutdown.rawValue
                : VivoDiagnosticFlags.reversibleShutdown.rawValue,
            nonFiniteCount: 0,
            boundViolationCount: 0,
            monitorViolationCount: 0,
            shutdownCount: 1,
            expressionFaultCount: 0,
            stochasticFallbackCount: 0,
            fluxTruncationCount: 0,
            eventOverflowCount: 0,
            firstCell: nil,
            firstSubject: nil,
            requestedResponse: status.rawValue,
            maximumSeverity: 0
        )
        let receipt = VivoPrepareReceipt(
            participantID: participantIdentifier,
            epochID: epoch.id,
            vote: vote,
            preparedStateSHA256: fingerprint,
            diagnosticFlags: diagnostics.flags,
            metadata: [
                "terminal": "true",
                "programContentFingerprint": program.pack.header.contentFingerprint
            ]
        )
        return VivoPreparedStep(
            receipt: receipt,
            requestedDeltaTime: Float(epoch.deltaTime),
            executedSubsteps: 0,
            projectedLogicalStep: logicalStep,
            projectedTime: absoluteTime,
            currentStateVersion: publication.stateVersion,
            diagnostics: diagnostics,
            events: []
        )
    }

    private func makeStepPlan(deltaTime: Float) throws -> StepPlan {
        let cohortLimit = program.pack.cohorts.lazy
            .filter {
                $0.reactionCount > 0 &&
                $0.maximumStableStep.isFinite &&
                $0.maximumStableStep > 0
            }
            .map(\.maximumStableStep)
            .min() ?? program.configuration.maximumSubstep
        let stableLimit = min(program.configuration.maximumSubstep, cohortLimit)
        guard stableLimit >= program.configuration.minimumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "Program requires a substep of \(stableLimit) s, below configured minimumSubstep \(program.configuration.minimumSubstep) s."
            )
        }
        let required = max(1, Int(ceil(Double(deltaTime / stableLimit))))
        guard required <= program.configuration.maximumInternalSubsteps else {
            throw VivoRuntimeError.resourceLimit(
                "Transaction requires \(required) internal substeps, exceeding maximumInternalSubsteps \(program.configuration.maximumInternalSubsteps)."
            )
        }
        let substep = deltaTime / Float(required)
        guard substep >= program.configuration.minimumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "Planned transaction substep \(substep) s is below configured minimumSubstep."
            )
        }
        return StepPlan(
            substeps: required,
            substepDeltaTime: substep,
            stableLimit: stableLimit
        )
    }

    private func resolveInputs(_ inputs: [VivoInputUpdate]) throws -> [NVivoInputUpdateSwift] {
        struct Key: Hashable { let species: UInt32; let cell: UInt32 }
        struct Accumulator { var value: Float; var mode: VivoInputUpdateMode }
        var resolved: [Key: Accumulator] = [:]
        resolved.reserveCapacity(inputs.count)

        for input in inputs {
            guard input.value.isFinite else {
                throw VivoRuntimeError.invalidConfiguration(
                    "Input '\(input.signal)' contains a non-finite value."
                )
            }
            guard input.cellIndex < UInt32(program.configuration.activeCellCount) else {
                throw VivoRuntimeError.invalidCellIndex(input.cellIndex)
            }
            guard let species = program.inputSpeciesByID[input.signal] else {
                if program.pack.speciesByID[input.signal] != nil {
                    throw VivoRuntimeError.inputNotExternallyOwned(input.signal)
                }
                throw VivoRuntimeError.unknownSignal(input.signal)
            }
            let key = Key(species: species.index, cell: input.cellIndex)
            switch input.mode {
            case .set:
                resolved[key] = Accumulator(value: input.value, mode: .set)
            case .add:
                if var existing = resolved[key] {
                    existing.value += input.value
                    guard existing.value.isFinite else {
                        throw VivoRuntimeError.invalidConfiguration(
                            "Accumulated input '\(input.signal)' is non-finite."
                        )
                    }
                    resolved[key] = existing
                } else {
                    resolved[key] = Accumulator(value: input.value, mode: .add)
                }
            }
        }

        return resolved
            .sorted { left, right in
                if left.key.species != right.key.species {
                    return left.key.species < right.key.species
                }
                return left.key.cell < right.key.cell
            }
            .map { key, value in
                NVivoInputUpdateSwift(
                    speciesIndex: key.species,
                    cellIndex: key.cell,
                    value: value.value,
                    mode: value.mode.rawValue
                )
            }
    }

    private func makeHashPlan(logicalStep: UInt64) throws -> HashPlan {
        let scalarBytes = MemoryLayout<Float>.stride
        let capacity = program.configuration.cellCapacity
        let contract = program.pack.runtimeContract
        let stateBytes = try checkedProduct(
            try checkedProduct(contract.speciesCount, capacity, label: "attested state elements"),
            scalarBytes,
            label: "attested state bytes"
        )
        let temporalBytes = try checkedProduct(
            try checkedProduct(contract.temporalStateCount, capacity, label: "attested temporal elements"),
            scalarBytes,
            label: "attested temporal bytes"
        )
        let delayBytes = try checkedProduct(
            try checkedProduct(
                try checkedProduct(contract.reactionCount, capacity, label: "attested reaction elements"),
                program.delaySlotCount,
                label: "attested delay elements"
            ),
            scalarBytes,
            label: "attested delay bytes"
        )
        let refractoryBytes = try checkedProduct(
            try checkedProduct(contract.ruleCount, capacity, label: "attested refractory elements"),
            scalarBytes,
            label: "attested refractory bytes"
        )

        let sources: [(MTLBuffer, Int, UInt32)] = [
            (program.buffers.shadowState, stateBytes, 1),
            (program.buffers.shadowTemporalState, temporalBytes, 2),
            (program.buffers.shadowDelayedFlux, delayBytes, 3),
            (program.buffers.shadowRuleRefractory, refractoryBytes, 4)
        ]
        var segments: [HashSegment] = []
        var wordOffset = 0
        for (buffer, byteCount, domain) in sources where byteCount > 0 {
            let chunkCount = (byteCount + hashChunkBytes - 1) / hashChunkBytes
            guard chunkCount <= Int(UInt32.max) else {
                throw VivoRuntimeError.resourceLimit(
                    "Attestation chunk count exceeds the Metal ABI."
                )
            }
            let words = try checkedProduct(chunkCount, 8, label: "attestation digest words")
            guard wordOffset <= Int(UInt32.max),
                  words <= Int(UInt32.max) - wordOffset else {
                throw VivoRuntimeError.resourceLimit(
                    "Attestation digest offsets exceed the Metal ABI."
                )
            }
            segments.append(
                HashSegment(
                    buffer: buffer,
                    byteCount: byteCount,
                    chunkCount: chunkCount,
                    outputWordOffset: wordOffset,
                    domain: domain
                )
            )
            wordOffset += words
        }
        let digestBytes = max(
            MemoryLayout<UInt32>.stride,
            try checkedProduct(wordOffset, MemoryLayout<UInt32>.stride, label: "attestation digest bytes")
        )
        guard digestBytes <= program.device.maxBufferLength,
              let buffer = program.device.makeBuffer(
                length: digestBytes,
                options: [.storageModeShared]
              ) else {
            throw VivoRuntimeError.allocation(
                "Unable to allocate prepared-state attestation output."
            )
        }
        buffer.label = "\(program.configuration.label).attestation.\(logicalStep)"
        buffer.contents().assumingMemoryBound(to: UInt8.self).initialize(
            repeating: 0,
            count: buffer.length
        )
        return HashPlan(
            buffer: buffer,
            segments: segments,
            wordCount: wordOffset,
            chunkBytes: hashChunkBytes
        )
    }

    private func encodeHashPlan(
        _ plan: HashPlan,
        logicalStep: UInt64,
        encoder: MTLComputeCommandEncoder
    ) {
        let pipeline = program.kernels[.sha256StateChunks]
        for segment in plan.segments {
            var uniforms = NVivoHashUniformsSwift(
                byteCountLow: UInt64(segment.byteCount).nvivoLow,
                byteCountHigh: UInt64(segment.byteCount).nvivoHigh,
                chunkBytes: UInt32(plan.chunkBytes),
                chunkCount: UInt32(segment.chunkCount),
                outputWordOffset: UInt32(segment.outputWordOffset),
                domain: segment.domain,
                logicalStepLow: logicalStep.nvivoLow,
                logicalStepHigh: logicalStep.nvivoHigh
            )
            encoder.setComputePipelineState(pipeline.state)
            encoder.setBuffer(segment.buffer, offset: 0, index: 0)
            encoder.setBuffer(plan.buffer, offset: 0, index: 1)
            encoder.setBytes(
                &uniforms,
                length: MemoryLayout<NVivoHashUniformsSwift>.stride,
                index: 2
            )
            encoder.dispatchThreads(
                MTLSize(width: segment.chunkCount, height: 1, depth: 1),
                threadsPerThreadgroup: pipeline.threadsPerThreadgroup()
            )
        }
    }

    private func attest(
        hashPlan: HashPlan,
        epoch: VivoTransactionEpoch,
        publication: NVivoPublicationSwift,
        diagnostics: VivoStepDiagnostics,
        events: [VivoEvent],
        vote: VivoPrepareVote
    ) throws -> String {
        var hasher = SHA256()
        hasher.update(data: Data("NumiVivoPreparedState/v1\0".utf8))
        update(&hasher, string: participantIdentifier)
        update(&hasher, string: program.pack.header.contentFingerprint)
        update(&hasher, data: try canonicalEpoch(epoch))
        update(&hasher, uint32: program.configuration.mode.rawValue)
        update(&hasher, uint64: program.configuration.seed)
        update(&hasher, uint64: UInt64(program.configuration.cellCapacity))
        update(&hasher, uint64: UInt64(program.configuration.activeCellCount))
        update(&hasher, uint32: publication.stateVersion)
        update(&hasher, string: vote.rawValue)

        let words = hashPlan.buffer.contents().bindMemory(
            to: UInt32.self,
            capacity: max(hashPlan.wordCount, 1)
        )
        for segment in hashPlan.segments.sorted(by: { $0.domain < $1.domain }) {
            update(&hasher, uint32: segment.domain)
            update(&hasher, uint64: UInt64(segment.byteCount))
            update(&hasher, uint32: UInt32(segment.chunkCount))
            for chunk in 0..<segment.chunkCount {
                let start = segment.outputWordOffset + chunk * 8
                var digest = Data()
                digest.reserveCapacity(32)
                for wordIndex in 0..<8 {
                    var word = words[start + wordIndex].bigEndian
                    Swift.withUnsafeBytes(of: &word) {
                        digest.append(contentsOf: $0)
                    }
                }
                hasher.update(data: digest)
            }
        }

        update(&hasher, uint32: diagnostics.flags)
        update(&hasher, uint32: diagnostics.nonFiniteCount)
        update(&hasher, uint32: diagnostics.boundViolationCount)
        update(&hasher, uint32: diagnostics.monitorViolationCount)
        update(&hasher, uint32: diagnostics.shutdownCount)
        update(&hasher, uint32: diagnostics.expressionFaultCount)
        update(&hasher, uint32: diagnostics.stochasticFallbackCount)
        update(&hasher, uint32: diagnostics.fluxTruncationCount)
        update(&hasher, uint32: diagnostics.eventOverflowCount)
        update(&hasher, uint32: diagnostics.requestedResponse)
        update(&hasher, uint32: diagnostics.maximumSeverity)
        update(&hasher, uint64: UInt64(events.count))
        for event in events {
            update(&hasher, uint32: event.cellIndex)
            update(&hasher, uint32: event.kind)
            update(&hasher, uint32: event.subject)
            update(&hasher, uint64: event.logicalStep)
            update(&hasher, uint32: event.value0.bitPattern)
            update(&hasher, uint32: event.value1.bitPattern)
            update(&hasher, uint32: event.flags)
        }
        return hex(hasher.finalize())
    }

    private func prepareVote(publicationStatus: UInt32) throws -> VivoPrepareVote {
        switch publicationStatus {
        case 5: return .prepared
        case VivoStepStatus.rejected.rawValue: return .reject
        case VivoStepStatus.substepRequired.rawValue: return .substepRequired
        case VivoStepStatus.reversibleShutdown.rawValue: return .reversibleShutdown
        case VivoStepStatus.permanentShutdown.rawValue: return .permanentShutdown
        default:
            throw VivoTransactionalRuntimeError.invalidPreparedPublication(
                publicationStatus
            )
        }
    }

    private func normalize(events: [VivoEvent]) -> [VivoEvent] {
        events.sorted { left, right in
            if left.logicalStep != right.logicalStep { return left.logicalStep < right.logicalStep }
            if left.cellIndex != right.cellIndex { return left.cellIndex < right.cellIndex }
            if left.kind != right.kind { return left.kind < right.kind }
            if left.subject != right.subject { return left.subject < right.subject }
            if left.value0.bitPattern != right.value0.bitPattern {
                return left.value0.bitPattern < right.value0.bitPattern
            }
            if left.value1.bitPattern != right.value1.bitPattern {
                return left.value1.bitPattern < right.value1.bitPattern
            }
            return left.flags < right.flags
        }
    }

    private func canonicalEpoch(_ epoch: VivoTransactionEpoch) throws -> Data {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return try encoder.encode(epoch)
    }

    private func randomToken(byteCount: Int) -> Data {
        var generator = SystemRandomNumberGenerator()
        var bytes = [UInt8](repeating: 0, count: byteCount)
        for index in bytes.indices {
            bytes[index] = UInt8.random(in: UInt8.min...UInt8.max, using: &generator)
        }
        return Data(bytes)
    }

    private func constantTimeEqual(_ left: Data, _ right: Data) -> Bool {
        guard left.count == right.count else { return false }
        var difference: UInt8 = 0
        for index in left.indices {
            difference |= left[index] ^ right[index]
        }
        return difference == 0
    }

    private func update(_ hasher: inout SHA256, string: String) {
        update(&hasher, data: Data(string.utf8))
    }

    private func update(_ hasher: inout SHA256, data: Data) {
        update(&hasher, uint64: UInt64(data.count))
        hasher.update(data: data)
    }

    private func update(_ hasher: inout SHA256, uint32 value: UInt32) {
        var little = value.littleEndian
        Swift.withUnsafeBytes(of: &little) {
            hasher.update(data: Data($0))
        }
    }

    private func update(_ hasher: inout SHA256, uint64 value: UInt64) {
        var little = value.littleEndian
        Swift.withUnsafeBytes(of: &little) {
            hasher.update(data: Data($0))
        }
    }

    private func hex<D: Sequence>(_ digest: D) -> String where D.Element == UInt8 {
        let alphabet = Array("0123456789abcdef".utf8)
        var output = [UInt8]()
        output.reserveCapacity(64)
        for byte in digest {
            output.append(alphabet[Int(byte >> 4)])
            output.append(alphabet[Int(byte & 0x0f)])
        }
        return String(decoding: output, as: UTF8.self)
    }

    private func encode(
        kernel: VivoMetalKernelLibrary.Kernel,
        grid: Int,
        uniforms: inout NVivoStepUniformsSwift,
        encoder: MTLComputeCommandEncoder
    ) {
        let pipeline = program.kernels[kernel]
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(program.buffers.argumentBuffer, offset: 0, index: 0)
        encoder.setBytes(
            &uniforms,
            length: MemoryLayout<NVivoStepUniformsSwift>.stride,
            index: 1
        )
        encoder.dispatchThreads(
            MTLSize(width: grid, height: 1, depth: 1),
            threadsPerThreadgroup: pipeline.threadsPerThreadgroup()
        )
    }

    private func checkedGrid(_ left: Int, _ right: Int, label: String) throws -> Int {
        let value = try checkedProduct(left, right, label: label)
        guard value <= Int(UInt32.max) else {
            throw VivoRuntimeError.resourceLimit("\(label) grid exceeds the Metal ABI.")
        }
        return value
    }

    private func checkedProduct(_ left: Int, _ right: Int, label: String) throws -> Int {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow, result.partialValue >= 0 else {
            throw VivoRuntimeError.resourceLimit(
                "Integer overflow while computing \(label)."
            )
        }
        return result.partialValue
    }

    private func awaitCompletion(_ commandBuffer: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation {
            (continuation: CheckedContinuation<Void, Error>) in
            commandBuffer.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(
                        throwing: VivoRuntimeError.commandExecution(
                            completed.error?.localizedDescription ??
                            "Metal command buffer failed."
                        )
                    )
                }
            }
            commandBuffer.commit()
        }
    }
}
