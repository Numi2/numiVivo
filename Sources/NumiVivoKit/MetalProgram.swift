@preconcurrency import Metal
import Foundation

final class VivoMetalProgram: @unchecked Sendable {
    struct Buffers: @unchecked Sendable {
        let committedState: MTLBuffer
        let shadowState: MTLBuffer
        let committedTemporalState: MTLBuffer
        let shadowTemporalState: MTLBuffer
        let parameters: MTLBuffer
        let species: MTLBuffer
        let reactionParameterIndices: MTLBuffer
        let stoichiometry: MTLBuffer
        let reactions: MTLBuffer
        let expressions: MTLBuffer
        let actions: MTLBuffer
        let rules: MTLBuffer
        let monitors: MTLBuffer
        let cohorts: MTLBuffer
        let speciesIncidenceOffsets: MTLBuffer
        let speciesIncidence: MTLBuffer
        let reactionFlux: MTLBuffer
        let diagnostics: MTLBuffer
        let publication: MTLBuffer
        let events: MTLBuffer
        let eventCount: MTLBuffer
        let cellActiveMask: MTLBuffer
        let committedDelayedFlux: MTLBuffer
        let shadowDelayedFlux: MTLBuffer
        let committedRuleRefractory: MTLBuffer
        let shadowRuleRefractory: MTLBuffer
        let argumentBuffer: MTLBuffer
        let sentinel: MTLBuffer
    }

    struct BufferSpec: Sendable {
        let name: String
        let length: Int
    }

    let pack: VivoProgramPack
    let configuration: VivoRuntimeConfiguration
    let device: MTLDevice
    let deviceInfo: VivoMetalDeviceInfo
    let commandQueue: MTLCommandQueue
    let kernels: VivoMetalKernelLibrary
    let privateHeap: MTLHeap
    let buffers: Buffers
    let memoryReport: VivoRuntimeMemoryReport
    let delaySlotCount: Int
    let maximumPreparationGrid: Int
    let maximumCommitGrid: Int
    let inputSpeciesByID: [String: VivoProgramPack.Species]

    init(
        pack: VivoProgramPack,
        configuration: VivoRuntimeConfiguration,
        device requestedDevice: MTLDevice? = nil
    ) throws {
        try Self.validate(configuration: configuration, pack: pack)
        guard let selectedDevice = requestedDevice ?? MTLCreateSystemDefaultDevice() else {
            throw VivoRuntimeError.metalUnavailable
        }
        guard Self.deviceHasUnifiedMemory(selectedDevice) else {
            throw VivoRuntimeError.invalidConfiguration(
                "NumiVivo production execution requires an Apple-silicon unified-memory Metal device."
            )
        }
        guard let queue = selectedDevice.makeCommandQueue() else {
            throw VivoRuntimeError.allocation("Unable to create a Metal command queue.")
        }
        queue.label = "\(configuration.label).queue"
        let kernelLibrary = try VivoMetalKernelLibrary(device: selectedDevice)

        let contract = pack.runtimeContract
        let capacity = configuration.cellCapacity
        let stateElements = try Self.checkedProduct(contract.speciesCount, capacity, label: "state elements")
        let temporalElements = try Self.checkedProduct(contract.temporalStateCount, capacity, label: "temporal-state elements")
        let reactionElements = try Self.checkedProduct(contract.reactionCount, capacity, label: "reaction-flux elements")
        let refractoryElements = try Self.checkedProduct(contract.ruleCount, capacity, label: "rule-refractory elements")

        let delaySlots: Int
        if pack.maximumReactionDelay > 0 {
            let computed = ceil(Double(pack.maximumReactionDelay) / Double(configuration.minimumSubstep)) + 2
            guard computed.isFinite, computed <= Double(UInt32.max) else {
                throw VivoRuntimeError.resourceLimit("Reaction delay horizon exceeds the Metal delay-wheel ABI.")
            }
            delaySlots = max(2, Int(computed))
        } else {
            delaySlots = 0
        }
        let delayElements = try Self.checkedProduct(reactionElements, delaySlots, label: "delayed-flux elements")

        let scalarBytes = MemoryLayout<Float>.stride
        let stateBytes = try Self.checkedProduct(stateElements, scalarBytes, label: "state bytes")
        let temporalBytes = try Self.checkedProduct(temporalElements, scalarBytes, label: "temporal-state bytes")
        let fluxBytes = try Self.checkedProduct(reactionElements, scalarBytes, label: "reaction-flux bytes")
        let delayBytesPerCopy = try Self.checkedProduct(delayElements, scalarBytes, label: "delay bytes")
        let refractoryBytes = try Self.checkedProduct(refractoryElements, scalarBytes, label: "refractory bytes")
        let totalDelayBytes = try Self.checkedProduct(delayBytesPerCopy, 2, label: "transactional delay bytes")
        guard totalDelayBytes <= configuration.maximumDelayBytes else {
            throw VivoRuntimeError.resourceLimit(
                "Exact delayed-reaction history requires \(totalDelayBytes) bytes, exceeding maximumDelayBytes \(configuration.maximumDelayBytes). Increase the explicit budget, increase minimumSubstep, reduce cell capacity, or compile a declared distributed-delay approximation."
            )
        }

        let staticData: [(name: String, data: Data)] = [
            ("parameters", try pack.materializeGPUParameters()),
            ("species", try pack.data(for: .species)),
            ("reactionParameterIndices", try pack.data(for: .reactionParameterIndices)),
            ("stoichiometry", try pack.data(for: .stoichiometry)),
            ("reactions", try pack.data(for: .reactions)),
            ("expressions", try pack.data(for: .expressions)),
            ("actions", try pack.data(for: .actions)),
            ("rules", try pack.data(for: .rules)),
            ("monitors", try pack.data(for: .monitors)),
            ("cohorts", try pack.data(for: .cohorts)),
            ("speciesIncidenceOffsets", try pack.data(for: .speciesIncidenceOffsets)),
            ("speciesIncidence", try pack.data(for: .speciesIncidence))
        ]
        var staticBytes = 0
        for item in staticData {
            staticBytes = try Self.checkedSum(staticBytes, max(item.data.count, 16), label: "static table bytes")
        }

        var specs: [BufferSpec] = [
            .init(name: "sentinel", length: 256),
            .init(name: "committedState", length: max(stateBytes, 16)),
            .init(name: "shadowState", length: max(stateBytes, 16)),
            .init(name: "committedTemporalState", length: max(temporalBytes, 16)),
            .init(name: "shadowTemporalState", length: max(temporalBytes, 16)),
            .init(name: "reactionFlux", length: max(fluxBytes, 16)),
            .init(name: "committedDelayedFlux", length: max(delayBytesPerCopy, 16)),
            .init(name: "shadowDelayedFlux", length: max(delayBytesPerCopy, 16)),
            .init(name: "committedRuleRefractory", length: max(refractoryBytes, 16)),
            .init(name: "shadowRuleRefractory", length: max(refractoryBytes, 16)),
            .init(name: "cellActiveMask", length: max(try Self.checkedProduct(capacity, MemoryLayout<UInt32>.stride, label: "cell mask bytes"), 16))
        ]
        specs.append(contentsOf: staticData.map { .init(name: $0.name, length: max($0.data.count, 16)) })

        let heapPlan = try Self.makeHeap(device: selectedDevice, specs: specs)
        guard heapPlan.heap.size <= configuration.maximumPrivateBytes else {
            throw VivoRuntimeError.resourceLimit(
                "Private Metal heap requires \(heapPlan.heap.size) bytes, exceeding maximumPrivateBytes \(configuration.maximumPrivateBytes)."
            )
        }
        let recommendedWorkingSet = Self.recommendedWorkingSet(device: selectedDevice)
        if recommendedWorkingSet > 0, heapPlan.heap.size > recommendedWorkingSet {
            throw VivoRuntimeError.resourceLimit(
                "Private Metal heap requires \(heapPlan.heap.size) bytes, above the device recommended working set \(recommendedWorkingSet)."
            )
        }

        func privateBuffer(_ name: String) throws -> MTLBuffer {
            guard let buffer = heapPlan.buffers[name] else {
                throw VivoRuntimeError.allocation("Missing planned private buffer '\(name)'.")
            }
            buffer.label = "\(configuration.label).\(name)"
            return buffer
        }

        let eventBytes = try Self.checkedProduct(
            configuration.eventCapacity,
            MemoryLayout<NVivoEventSwift>.stride,
            label: "event buffer bytes"
        )
        guard let diagnostics = selectedDevice.makeBuffer(
            length: MemoryLayout<NVivoDiagnosticsSwift>.stride,
            options: [.storageModeShared]
        ), let publication = selectedDevice.makeBuffer(
            length: MemoryLayout<NVivoPublicationSwift>.stride,
            options: [.storageModeShared]
        ), let events = selectedDevice.makeBuffer(
            length: eventBytes,
            options: [.storageModeShared]
        ), let eventCount = selectedDevice.makeBuffer(
            length: MemoryLayout<UInt32>.stride,
            options: [.storageModeShared]
        ), let argumentBuffer = selectedDevice.makeBuffer(
            length: kernelLibrary.argumentEncoder.encodedLength,
            options: [.storageModeShared]
        ) else {
            throw VivoRuntimeError.allocation("Unable to allocate shared runtime boundary buffers.")
        }
        diagnostics.label = "\(configuration.label).diagnostics"
        publication.label = "\(configuration.label).publication"
        events.label = "\(configuration.label).events"
        eventCount.label = "\(configuration.label).eventCount"
        argumentBuffer.label = "\(configuration.label).arguments"
        Self.zero(diagnostics)
        Self.zero(publication)
        Self.zero(events)
        Self.zero(eventCount)
        Self.zero(argumentBuffer)

        let constructedBuffers = Buffers(
            committedState: try privateBuffer("committedState"),
            shadowState: try privateBuffer("shadowState"),
            committedTemporalState: try privateBuffer("committedTemporalState"),
            shadowTemporalState: try privateBuffer("shadowTemporalState"),
            parameters: try privateBuffer("parameters"),
            species: try privateBuffer("species"),
            reactionParameterIndices: try privateBuffer("reactionParameterIndices"),
            stoichiometry: try privateBuffer("stoichiometry"),
            reactions: try privateBuffer("reactions"),
            expressions: try privateBuffer("expressions"),
            actions: try privateBuffer("actions"),
            rules: try privateBuffer("rules"),
            monitors: try privateBuffer("monitors"),
            cohorts: try privateBuffer("cohorts"),
            speciesIncidenceOffsets: try privateBuffer("speciesIncidenceOffsets"),
            speciesIncidence: try privateBuffer("speciesIncidence"),
            reactionFlux: try privateBuffer("reactionFlux"),
            diagnostics: diagnostics,
            publication: publication,
            events: events,
            eventCount: eventCount,
            cellActiveMask: try privateBuffer("cellActiveMask"),
            committedDelayedFlux: try privateBuffer("committedDelayedFlux"),
            shadowDelayedFlux: try privateBuffer("shadowDelayedFlux"),
            committedRuleRefractory: try privateBuffer("committedRuleRefractory"),
            shadowRuleRefractory: try privateBuffer("shadowRuleRefractory"),
            argumentBuffer: argumentBuffer,
            sentinel: try privateBuffer("sentinel")
        )

        kernelLibrary.argumentEncoder.setArgumentBuffer(argumentBuffer, offset: 0)
        let argumentBuffers: [MTLBuffer] = [
            constructedBuffers.committedState,
            constructedBuffers.shadowState,
            constructedBuffers.committedTemporalState,
            constructedBuffers.shadowTemporalState,
            constructedBuffers.parameters,
            constructedBuffers.species,
            constructedBuffers.reactionParameterIndices,
            constructedBuffers.stoichiometry,
            constructedBuffers.reactions,
            constructedBuffers.expressions,
            constructedBuffers.actions,
            constructedBuffers.rules,
            constructedBuffers.monitors,
            constructedBuffers.cohorts,
            constructedBuffers.speciesIncidenceOffsets,
            constructedBuffers.speciesIncidence,
            constructedBuffers.reactionFlux,
            constructedBuffers.diagnostics,
            constructedBuffers.publication,
            constructedBuffers.events,
            constructedBuffers.eventCount,
            constructedBuffers.cellActiveMask,
            constructedBuffers.committedDelayedFlux,
            constructedBuffers.shadowDelayedFlux,
            constructedBuffers.committedRuleRefractory,
            constructedBuffers.shadowRuleRefractory
        ]
        for (index, buffer) in argumentBuffers.enumerated() {
            kernelLibrary.argumentEncoder.setBuffer(buffer, offset: 0, index: index)
        }

        let sharedBytes = diagnostics.length + publication.length + events.length +
                          eventCount.length + argumentBuffer.length
        let report = VivoRuntimeMemoryReport(
            privateBytes: heapPlan.heap.size,
            sharedBytes: sharedBytes,
            stateBytes: stateBytes * 2,
            temporalStateBytes: temporalBytes * 2,
            reactionFluxBytes: fluxBytes,
            delayBytes: totalDelayBytes,
            refractoryBytes: refractoryBytes * 2,
            staticTableBytes: staticBytes,
            delaySlotCount: delaySlots
        )
        let preparationGrid = max(1, stateElements, temporalElements, reactionElements, delayElements, refractoryElements)
        let commitGrid = max(1, stateElements, temporalElements, delayElements, refractoryElements)
        guard preparationGrid <= Int(UInt32.max), commitGrid <= Int(UInt32.max) else {
            throw VivoRuntimeError.resourceLimit(
                "A one-dimensional Metal dispatch exceeds UInt32.max elements. Reduce program or cell capacity."
            )
        }
        let inputMap = Dictionary(
            uniqueKeysWithValues: pack.species.filter(\.isExternallyOwned).map { ($0.id, $0) }
        )

        self.pack = pack
        self.configuration = configuration
        self.device = selectedDevice
        self.deviceInfo = selectedDevice.nvivoInfo
        self.commandQueue = queue
        self.kernels = kernelLibrary
        self.privateHeap = heapPlan.heap
        self.buffers = constructedBuffers
        self.memoryReport = report
        self.delaySlotCount = delaySlots
        self.maximumPreparationGrid = preparationGrid
        self.maximumCommitGrid = commitGrid
        self.inputSpeciesByID = inputMap

        try uploadAndInitialize(staticData: staticData)
    }

    func declareArgumentResources(on encoder: MTLComputeCommandEncoder) {
        encoder.useHeap(privateHeap)
        encoder.useResource(buffers.diagnostics, usage: [.read, .write])
        encoder.useResource(buffers.publication, usage: [.read, .write])
        encoder.useResource(buffers.events, usage: [.read, .write])
        encoder.useResource(buffers.eventCount, usage: [.read, .write])
    }

    private func uploadAndInitialize(staticData: [(name: String, data: Data)]) throws {
        guard let commandBuffer = commandQueue.makeCommandBuffer(),
              let blit = commandBuffer.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding("Unable to create initialization command encoders.")
        }
        commandBuffer.label = "\(configuration.label).initialize"
        var stagingBuffers: [MTLBuffer] = []
        stagingBuffers.reserveCapacity(staticData.count + 1)

        let destinations: [String: MTLBuffer] = [
            "parameters": buffers.parameters,
            "species": buffers.species,
            "reactionParameterIndices": buffers.reactionParameterIndices,
            "stoichiometry": buffers.stoichiometry,
            "reactions": buffers.reactions,
            "expressions": buffers.expressions,
            "actions": buffers.actions,
            "rules": buffers.rules,
            "monitors": buffers.monitors,
            "cohorts": buffers.cohorts,
            "speciesIncidenceOffsets": buffers.speciesIncidenceOffsets,
            "speciesIncidence": buffers.speciesIncidence
        ]

        for item in staticData where !item.data.isEmpty {
            guard let destination = destinations[item.name] else {
                throw VivoRuntimeError.commandEncoding("No destination for static table '\(item.name)'.")
            }
            let staging = item.data.withUnsafeBytes { rawBuffer -> MTLBuffer? in
                guard let baseAddress = rawBuffer.baseAddress else { return nil }
                return device.makeBuffer(
                    bytes: baseAddress,
                    length: rawBuffer.count,
                    options: [.storageModeShared]
                )
            }
            guard let staging else {
                throw VivoRuntimeError.allocation("Unable to allocate staging buffer for '\(item.name)'.")
            }
            staging.label = "\(configuration.label).staging.\(item.name)"
            stagingBuffers.append(staging)
            blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: item.data.count)
        }

        var activeMask = [UInt32](repeating: 0, count: configuration.cellCapacity)
        for index in 0..<configuration.activeCellCount { activeMask[index] = 1 }
        let activeStaging = activeMask.withUnsafeBufferPointer { pointer -> MTLBuffer? in
            guard let baseAddress = pointer.baseAddress else { return nil }
            return device.makeBuffer(
                bytes: baseAddress,
                length: pointer.count * MemoryLayout<UInt32>.stride,
                options: [.storageModeShared]
            )
        }
        guard let activeStaging else {
            throw VivoRuntimeError.allocation("Unable to allocate active-cell mask staging buffer.")
        }
        activeStaging.label = "\(configuration.label).staging.cellActiveMask"
        stagingBuffers.append(activeStaging)
        blit.copy(
            from: activeStaging,
            sourceOffset: 0,
            to: buffers.cellActiveMask,
            destinationOffset: 0,
            size: activeMask.count * MemoryLayout<UInt32>.stride
        )
        blit.endEncoding()

        guard let compute = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.commandEncoding("Unable to create initialization compute encoder.")
        }
        compute.label = "\(configuration.label).initialize.compute"
        declareArgumentResources(on: compute)
        let pipeline = kernels[.initializeProgram]
        compute.setComputePipelineState(pipeline.state)
        compute.setBuffer(buffers.argumentBuffer, offset: 0, index: 0)
        var uniforms = makeUniforms(
            deltaTime: configuration.minimumSubstep,
            absoluteTime: 0,
            logicalStep: 0,
            substepIndex: 0,
            delayWriteSlot: 0
        )
        compute.setBytes(&uniforms, length: MemoryLayout<NVivoStepUniformsSwift>.stride, index: 1)
        compute.dispatchThreads(
            MTLSize(width: maximumPreparationGrid, height: 1, depth: 1),
            threadsPerThreadgroup: pipeline.threadsPerThreadgroup()
        )
        compute.endEncoding()

        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        guard commandBuffer.status == .completed else {
            throw VivoRuntimeError.commandExecution(
                "Metal initialization failed: \(commandBuffer.error?.localizedDescription ?? "unknown error")"
            )
        }
        _ = stagingBuffers
    }

    func makeUniforms(
        deltaTime: Float,
        absoluteTime: Float,
        logicalStep: UInt64,
        substepIndex: UInt32,
        delayWriteSlot: UInt32
    ) -> NVivoStepUniformsSwift {
        NVivoStepUniformsSwift(
            activeCellCount: UInt32(configuration.activeCellCount),
            cellCapacity: UInt32(configuration.cellCapacity),
            speciesCount: UInt32(pack.runtimeContract.speciesCount),
            parameterCount: UInt32(pack.runtimeContract.parameterCount),
            reactionCount: UInt32(pack.runtimeContract.reactionCount),
            ruleCount: UInt32(pack.runtimeContract.ruleCount),
            monitorCount: UInt32(pack.runtimeContract.monitorCount),
            temporalStateCount: UInt32(pack.runtimeContract.temporalStateCount),
            deltaTime: deltaTime,
            absoluteTime: absoluteTime,
            logicalStepLow: logicalStep.nvivoLow,
            logicalStepHigh: logicalStep.nvivoHigh,
            seedLow: configuration.seed.nvivoLow,
            seedHigh: configuration.seed.nvivoHigh,
            mode: configuration.mode.rawValue,
            substepIndex: substepIndex,
            eventCapacity: UInt32(configuration.eventCapacity),
            featureFlags: pack.runtimeContract.featureFlags,
            delaySlotCount: UInt32(delaySlotCount),
            delayWriteSlot: delayWriteSlot
        )
    }

    func readPublication() -> NVivoPublicationSwift {
        buffers.publication.contents().load(as: NVivoPublicationSwift.self)
    }

    func readDiagnostics() -> NVivoDiagnosticsSwift {
        buffers.diagnostics.contents().load(as: NVivoDiagnosticsSwift.self)
    }

    func readEvents(count: Int) -> [VivoEvent] {
        let boundedCount = max(0, min(count, configuration.eventCapacity))
        guard boundedCount > 0 else { return [] }
        let pointer = buffers.events.contents().bindMemory(to: NVivoEventSwift.self, capacity: configuration.eventCapacity)
        return (0..<boundedCount).map { index in
            let event = pointer[index]
            return VivoEvent(
                cellIndex: event.cellIndex,
                kind: event.kind,
                subject: event.subject,
                logicalStep: UInt64(nvivoLow: event.logicalStepLow, high: event.logicalStepHigh),
                value0: event.value0,
                value1: event.value1,
                flags: event.flags
            )
        }
    }

    private static func validate(
        configuration: VivoRuntimeConfiguration,
        pack: VivoProgramPack
    ) throws {
        guard configuration.cellCapacity > 0 else {
            throw VivoRuntimeError.invalidConfiguration("cellCapacity must be positive.")
        }
        guard configuration.activeCellCount > 0,
              configuration.activeCellCount <= configuration.cellCapacity else {
            throw VivoRuntimeError.invalidConfiguration("activeCellCount must be in 1...cellCapacity.")
        }
        guard configuration.maximumSubstep.isFinite,
              configuration.maximumSubstep > 0,
              configuration.minimumSubstep.isFinite,
              configuration.minimumSubstep > 0,
              configuration.minimumSubstep <= configuration.maximumSubstep else {
            throw VivoRuntimeError.invalidConfiguration(
                "minimumSubstep and maximumSubstep must be finite, positive, and ordered."
            )
        }
        guard configuration.maximumInternalSubsteps > 0,
              configuration.eventCapacity > 0,
              configuration.maximumPrivateBytes > 0,
              configuration.maximumDelayBytes >= 0 else {
            throw VivoRuntimeError.invalidConfiguration("Runtime capacities must be positive.")
        }
        guard configuration.cellCapacity <= Int(UInt32.max),
              configuration.activeCellCount <= Int(UInt32.max),
              configuration.eventCapacity <= Int(UInt32.max) else {
            throw VivoRuntimeError.resourceLimit("Runtime capacities exceed the Metal ABI's 32-bit fields.")
        }
        guard pack.header.fidelity == .f1Deterministic || pack.header.fidelity == .f2Stochastic else {
            throw VivoRuntimeError.invalidConfiguration(
                "The well-mixed runtime accepts F1 or F2 ProgramPacks. Use the spatial or tissue runtime for F3/F4 packs."
            )
        }
        if configuration.mode == .stochastic, pack.header.fidelity != .f2Stochastic {
            throw VivoRuntimeError.invalidConfiguration(
                "Stochastic execution requires an F2 ProgramPack."
            )
        }
        guard MemoryLayout<NVivoStepUniformsSwift>.stride == 80,
              MemoryLayout<NVivoCohortUniformsSwift>.stride == 16,
              MemoryLayout<NVivoInputUpdateSwift>.stride == 16,
              MemoryLayout<NVivoDiagnosticsSwift>.stride == 64,
              MemoryLayout<NVivoPublicationSwift>.stride == 32,
              MemoryLayout<NVivoEventSwift>.stride == 32 else {
            throw VivoRuntimeError.invalidConfiguration(
                "Swift compiler layout does not match NumiVivo Metal ABI v1."
            )
        }
    }

    private static func checkedProduct(_ left: Int, _ right: Int, label: String) throws -> Int {
        let result = left.multipliedReportingOverflow(by: right)
        guard !result.overflow, result.partialValue >= 0 else {
            throw VivoRuntimeError.resourceLimit("Integer overflow while computing \(label).")
        }
        return result.partialValue
    }

    private static func checkedSum(_ left: Int, _ right: Int, label: String) throws -> Int {
        let result = left.addingReportingOverflow(right)
        guard !result.overflow, result.partialValue >= 0 else {
            throw VivoRuntimeError.resourceLimit("Integer overflow while computing \(label).")
        }
        return result.partialValue
    }

    private static func makeHeap(
        device: MTLDevice,
        specs: [BufferSpec]
    ) throws -> (heap: MTLHeap, buffers: [String: MTLBuffer]) {
        var required = 0
        for spec in specs {
            guard spec.length > 0, spec.length <= device.maxBufferLength else {
                throw VivoRuntimeError.resourceLimit(
                    "Buffer '\(spec.name)' length \(spec.length) is unsupported by the selected Metal device."
                )
            }
            let sizeAndAlign = device.heapBufferSizeAndAlign(length: spec.length, options: [.storageModePrivate])
            required = try align(required, to: sizeAndAlign.align)
            required = try checkedSum(required, sizeAndAlign.size, label: "Metal heap size")
        }

        let descriptor = MTLHeapDescriptor()
        descriptor.size = required
        descriptor.storageMode = .private
        descriptor.cpuCacheMode = .defaultCache
        descriptor.hazardTrackingMode = .tracked
        guard let heap = device.makeHeap(descriptor: descriptor) else {
            throw VivoRuntimeError.allocation("Unable to allocate a \(required)-byte private Metal heap.")
        }
        heap.label = "NumiVivo.ProgramHeap"

        var buffers: [String: MTLBuffer] = [:]
        buffers.reserveCapacity(specs.count)
        for spec in specs {
            guard let buffer = heap.makeBuffer(length: spec.length, options: [.storageModePrivate]) else {
                throw VivoRuntimeError.allocation(
                    "Unable to place private buffer '\(spec.name)' in the program heap."
                )
            }
            buffers[spec.name] = buffer
        }
        return (heap, buffers)
    }

    private static func align(_ value: Int, to alignment: Int) throws -> Int {
        guard alignment > 0, alignment & (alignment - 1) == 0 else {
            throw VivoRuntimeError.allocation("Metal reported an invalid heap alignment.")
        }
        let mask = alignment - 1
        let sum = value.addingReportingOverflow(mask)
        guard !sum.overflow else { throw VivoRuntimeError.resourceLimit("Heap alignment overflow.") }
        return sum.partialValue & ~mask
    }

    private static func zero(_ buffer: MTLBuffer) {
        buffer.contents().assumingMemoryBound(to: UInt8.self).initialize(repeating: 0, count: buffer.length)
    }

    private static func recommendedWorkingSet(device: MTLDevice) -> Int {
        #if os(macOS)
        return Int(clamping: device.recommendedMaxWorkingSetSize)
        #else
        return 0
        #endif
    }

    private static func deviceHasUnifiedMemory(_ device: MTLDevice) -> Bool {
        #if os(macOS)
        return device.hasUnifiedMemory
        #else
        return true
        #endif
    }
}
