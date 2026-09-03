@preconcurrency import Metal
import Foundation

public actor NumiVivoRuntime {
    private let program: VivoMetalProgram
    private var logicalStep: UInt64 = 0
    private var absoluteTime: Double = 0
    private var delayWriteSlot: Int = 0
    private var terminalStatus: VivoStepStatus?

    public nonisolated var programPack: VivoProgramPack { program.pack }
    public nonisolated var deviceInfo: VivoMetalDeviceInfo { program.deviceInfo }
    public nonisolated var memoryReport: VivoRuntimeMemoryReport { program.memoryReport }
    public nonisolated var configuration: VivoRuntimeConfiguration { program.configuration }

    public init(
        programPack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device: MTLDevice? = nil
    ) throws {
        self.program = try VivoMetalProgram(
            pack: programPack,
            configuration: configuration,
            device: device
        )
    }

    public init(
        sourceJSON: Data,
        compilerConfiguration: VivoCompilerConfiguration = .init(),
        runtimeConfiguration: VivoRuntimeConfiguration,
        device: MTLDevice? = nil
    ) throws {
        let output = try VivoCompiler.compile(
            json: sourceJSON,
            configuration: compilerConfiguration
        )
        self.program = try VivoMetalProgram(
            pack: output.programPack,
            configuration: runtimeConfiguration,
            device: device
        )
    }

    public func stateClock() -> (logicalStep: UInt64, time: Double) {
        (logicalStep, absoluteTime)
    }

    public func step(
        deltaTime requestedDeltaTime: Float,
        inputs: [VivoInputUpdate] = []
    ) async throws -> VivoStepResult {
        if let terminalStatus {
            throw VivoRuntimeError.runtimeShutdown(terminalStatus)
        }
        guard requestedDeltaTime.isFinite, requestedDeltaTime > 0 else {
            throw VivoRuntimeError.invalidConfiguration("Step deltaTime must be finite and positive.")
        }
        guard requestedDeltaTime >= program.configuration.minimumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "Step deltaTime \(requestedDeltaTime) is below configured minimumSubstep \(program.configuration.minimumSubstep)."
            )
        }

        let plan = try makeStepPlan(deltaTime: requestedDeltaTime)
        let rawInputs = try resolveInputs(inputs)
        let nextLogicalStep = logicalStep.addingReportingOverflow(1)
        guard !nextLogicalStep.overflow else {
            throw VivoRuntimeError.resourceLimit("Logical step counter overflow.")
        }

        guard let commandBuffer = program.commandQueue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding("Unable to create a Metal step command encoder.")
        }
        commandBuffer.label = "\(program.configuration.label).step.\(nextLogicalStep.partialValue)"
        encoder.label = "\(program.configuration.label).step.compute"
        program.declareArgumentResources(on: encoder)

        var initialUniforms = program.makeUniforms(
            deltaTime: plan.substepDeltaTime,
            absoluteTime: Float(absoluteTime),
            logicalStep: nextLogicalStep.partialValue,
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
                throw VivoRuntimeError.allocation("Unable to allocate staged input updates.")
            }
            retainedInputBuffer.label = "\(program.configuration.label).step.inputs"
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
        var projectedDelaySlot = delayWriteSlot
        for substep in 0..<plan.substeps {
            if program.delaySlotCount > 0 {
                projectedDelaySlot = (delayWriteSlot + substep) % program.delaySlotCount
            } else {
                projectedDelaySlot = 0
            }
            let substepEndTime = absoluteTime + Double(plan.substepDeltaTime) * Double(substep + 1)
            var uniforms = program.makeUniforms(
                deltaTime: plan.substepDeltaTime,
                absoluteTime: Float(substepEndTime),
                logicalStep: nextLogicalStep.partialValue,
                substepIndex: UInt32(substep),
                delayWriteSlot: UInt32(projectedDelaySlot)
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
                    label: "reaction cohort"
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
            if !program.pack.cohorts.isEmpty { encoder.memoryBarrier(scope: .buffers) }

            let incidenceGrid = try checkedGrid(
                program.pack.runtimeContract.speciesCount,
                program.configuration.activeCellCount,
                label: "incidence"
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
            kernel: .commitIfValid,
            grid: program.maximumCommitGrid,
            uniforms: &finalUniforms,
            encoder: encoder
        )
        encoder.endEncoding()

        try await awaitCompletion(commandBuffer)
        _ = retainedInputBuffer

        let publication = program.readPublication()
        let diagnostics = program.readDiagnostics().publicValue
        let status = VivoStepStatus(rawValue: publication.status) ?? .rejected
        let committedStep = UInt64(
            nvivoLow: publication.committedStepLow,
            high: publication.committedStepHigh
        )

        if status == .committed {
            logicalStep = nextLogicalStep.partialValue
            absoluteTime += Double(requestedDeltaTime)
            if program.delaySlotCount > 0 {
                delayWriteSlot = (delayWriteSlot + plan.substeps) % program.delaySlotCount
            }
        } else if status == .reversibleShutdown || status == .permanentShutdown {
            terminalStatus = status
        }

        return VivoStepResult(
            status: status,
            requestedDeltaTime: requestedDeltaTime,
            executedSubsteps: plan.substeps,
            committedLogicalStep: committedStep,
            committedTime: absoluteTime,
            stateVersion: publication.stateVersion,
            diagnostics: diagnostics,
            events: program.readEvents(count: Int(publication.eventCount))
        )
    }

    public func readState(
        species identifiers: [String],
        cells requestedCells: Range<Int>? = nil
    ) async throws -> [VivoStateSlice] {
        let speciesByID = program.pack.speciesByID
        let cells = requestedCells ?? 0..<program.configuration.activeCellCount
        guard cells.lowerBound >= 0,
              cells.upperBound <= program.configuration.activeCellCount,
              !cells.isEmpty else {
            throw VivoRuntimeError.invalidConfiguration("Requested cell range is outside active cells.")
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
            label: "state readback"
        )
        let totalBytes = try checkedProduct(bytesPerSpecies, selected.count, label: "state readback")
        guard totalBytes <= program.device.maxBufferLength else {
            throw VivoRuntimeError.resourceLimit("State readback exceeds the device maximum buffer length.")
        }
        guard let readback = program.device.makeBuffer(
            length: totalBytes,
            options: [.storageModeShared]
        ), let commandBuffer = program.commandQueue.makeCommandBuffer(),
           let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.allocation("Unable to allocate state readback resources.")
        }
        readback.label = "\(program.configuration.label).stateReadback"
        commandBuffer.label = "\(program.configuration.label).readState"

        for (outputIndex, species) in selected.enumerated() {
            let sourceElement = try checkedProduct(
                Int(species.index),
                program.configuration.cellCapacity,
                label: "state readback source"
            ) + cells.lowerBound
            let sourceOffset = try checkedProduct(
                sourceElement,
                MemoryLayout<Float>.stride,
                label: "state readback source bytes"
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
                values: Array(UnsafeBufferPointer(start: pointer + start, count: valuesPerSpecies))
            )
        }
    }

    public func resumeAfterReversibleShutdown() throws {
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

    private struct StepPlan {
        let substeps: Int
        let substepDeltaTime: Float
    }

    private func makeStepPlan(deltaTime: Float) throws -> StepPlan {
        let cohortLimit = program.pack.cohorts.lazy
            .filter { $0.reactionCount > 0 && $0.maximumStableStep.isFinite && $0.maximumStableStep > 0 }
            .map(\.maximumStableStep)
            .min() ?? program.configuration.maximumSubstep
        let stableStep = min(program.configuration.maximumSubstep, cohortLimit)
        guard stableStep >= program.configuration.minimumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "Program requires a substep of \(stableStep) s, below configured minimumSubstep \(program.configuration.minimumSubstep) s."
            )
        }
        let required = max(1, Int(ceil(Double(deltaTime / stableStep))))
        guard required <= program.configuration.maximumInternalSubsteps else {
            throw VivoRuntimeError.resourceLimit(
                "Step requires \(required) internal substeps, exceeding maximumInternalSubsteps \(program.configuration.maximumInternalSubsteps)."
            )
        }
        let substep = deltaTime / Float(required)
        guard substep >= program.configuration.minimumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "Planned substep \(substep) s is below configured minimumSubstep."
            )
        }
        return StepPlan(substeps: required, substepDeltaTime: substep)
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
            throw VivoRuntimeError.resourceLimit("Integer overflow while computing \(label).")
        }
        return result.partialValue
    }

    private func awaitCompletion(_ commandBuffer: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
            commandBuffer.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(
                        throwing: VivoRuntimeError.commandExecution(
                            completed.error?.localizedDescription ?? "Metal command buffer failed."
                        )
                    )
                }
            }
            commandBuffer.commit()
        }
    }
}
