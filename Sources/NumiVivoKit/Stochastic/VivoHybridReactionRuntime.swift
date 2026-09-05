import Foundation
@preconcurrency import Metal
import NumiVivoShaders

/// Executes disjoint exact-SSA, tau-leap and RK2 cohorts on one Metal timeline.
/// The two count buffers are UInt32, not Float reinterpretations. CPU readback
/// during a step is limited to a 32-byte status and explicitly requested outputs.
public actor VivoHybridReactionRuntime {
    public nonisolated let execution: VivoCompiledHybridExecution
    public nonisolated let configuration: VivoHybridRuntimeConfiguration
    public nonisolated let deviceName: String
    public nonisolated let registryID: UInt64
    private let device: MTLDevice
    private let queue: MTLCommandQueue
    private let pipelines: NumiVivoHybridPipelineSet
    private var arena: VivoHybridGPUArena
    private var nextStep: UInt64 = 0
    private var absoluteTime: Double = 0
    private var busy = false
    private var failure: String?
    private var pending: VivoHybridStepCertificate?

    public static func make(model: VivoExactSSAModel, plan: VivoHybridStochasticPlan,
                            configuration: VivoHybridRuntimeConfiguration,
                            initialCounts: [UInt32], device requestedDevice: MTLDevice? = nil) async throws -> VivoHybridReactionRuntime {
        try configuration.validate()
        let execution = try VivoHybridExecutionCompiler.compile(model: model, plan: plan)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoHybridExecutionError.metal("an Apple-silicon unified-memory device and command queue are required")
        }
        let count = UInt64(execution.species.count) * UInt64(configuration.laneCount)
        guard count <= UInt64(UInt32.max), count == UInt64(initialCounts.count) else {
            throw VivoHybridExecutionError.invalidState("initial counts must be species-major with one element per species and lane")
        }
        var continuous = [Float](repeating: 0, count: initialCounts.count)
        var counts = initialCounts
        for species in execution.authorities.indices where execution.authorities[species] == .deterministicRK2 {
            let start = species * Int(configuration.laneCount)
            for i in start..<(start + Int(configuration.laneCount)) {
                continuous[i] = Float(counts[i])
                counts[i] = 0
            }
        }
        let pipelines = try await NumiVivoHybridPipelineCache.shared.pipelines(device: device)
        let arena = try VivoHybridGPUArena(device: device, queue: queue, execution: execution,
                                           configuration: configuration, continuous: continuous, counts: counts)
        return VivoHybridReactionRuntime(execution: execution, configuration: configuration,
                                          device: device, queue: queue, pipelines: pipelines, arena: arena)
    }

    private init(execution: VivoCompiledHybridExecution, configuration: VivoHybridRuntimeConfiguration,
                 device: MTLDevice, queue: MTLCommandQueue,
                 pipelines: NumiVivoHybridPipelineSet, arena: VivoHybridGPUArena) {
        self.execution = execution
        self.configuration = configuration
        self.device = device
        self.queue = queue
        self.pipelines = pipelines
        self.arena = arena
        deviceName = device.name
        registryID = device.registryID
    }

    public func timeSeconds() -> Double { absoluteTime }
    public func stepIndex() -> UInt64 { nextStep }
    public func hasPendingTransaction() -> Bool { busy || pending != nil }

    public func prepareStep(deltaTime: Float? = nil,
                            publications: [VivoHybridPublicationRequest] = [],
                            transactionID: UUID = UUID()) async throws -> VivoHybridStepCertificate {
        try requireIdle()
        try Task.checkCancellation()
        let dt = deltaTime ?? configuration.timeStep
        guard dt.isFinite, dt >= configuration.minimumTimeStep, dt <= configuration.maximumTimeStep,
              nextStep < UInt64.max, absoluteTime + Double(dt) > absoluteTime,
              (absoluteTime + Double(dt)).isFinite else {
            throw VivoHybridExecutionError.invalidConfiguration
        }
        guard publications.count <= Int(configuration.maximumPublications) else {
            throw VivoHybridExecutionError.resourceLimit("publication capacity exceeded")
        }
        var indices: [UInt32] = []
        for publication in publications {
            guard Int(publication.speciesIndex) < execution.species.count,
                  publication.laneIndex < configuration.laneCount else {
                throw VivoHybridExecutionError.invalidState("publication index is out of bounds")
            }
            indices.append(UInt32(UInt64(publication.speciesIndex) * UInt64(configuration.laneCount) + UInt64(publication.laneIndex)))
        }
        busy = true
        defer { busy = false }
        try arena.writePublicationIndices(indices)
        var uniforms = command(dt: dt)
        let first = try makeCommand("Prepare")
        guard let blit = first.makeBlitCommandEncoder() else { throw VivoHybridExecutionError.metal("blit encoder unavailable") }
        blit.copy(from: arena.buffer(0), sourceOffset: 0, to: arena.buffer(2), destinationOffset: 0, size: arena.stateBytes)
        blit.copy(from: arena.buffer(0), sourceOffset: 0, to: arena.buffer(1), destinationOffset: 0, size: arena.stateBytes)
        blit.copy(from: arena.buffer(4), sourceOffset: 0, to: arena.buffer(5), destinationOffset: 0, size: arena.stateBytes)
        blit.fill(buffer: arena.buffer(16), range: 0..<arena.buffer(16).length, value: 0)
        blit.endEncoding()
        try encode("clear", count: 1, into: first, uniforms: uniforms)
        if execution.reactionAuthorities.contains(.deterministicRK2) {
            try encode("rk_rates", count: arena.reactionWork, into: first, uniforms: uniforms)
            try encode("rk_predict", count: arena.stateElements, into: first, uniforms: uniforms)
            uniforms.stage = 1
            try encode("rk_rates", count: arena.reactionWork, into: first, uniforms: uniforms)
            try encode("rk_correct", count: arena.stateElements, into: first, uniforms: uniforms)
            uniforms.stage = 0
        }
        if execution.reactionAuthorities.contains(.tauLeap) {
            try encode("tau_sample", count: arena.reactionWork, into: first, uniforms: uniforms)
            try encode("tau_apply", count: arena.stateElements, into: first, uniforms: uniforms)
        }
        var exactDispatches: UInt32 = 0
        if arena.exactWork > 0 {
            try encode("exact_advance", count: arena.exactWork, into: first, uniforms: uniforms)
            exactDispatches = 1
        }
        try await complete(first)
        var status = arena.status()
        // Continue the SAME candidate and SAME per-lane random cursor. Reducing
        // a horizon based on a realized event count can bias a supposedly exact
        // process, so work exhaustion never silently becomes timestep reduction.
        while status.valid && status.unfinishedExactLanes > 0 && exactDispatches < configuration.maximumExactDispatches {
            try Task.checkCancellation()
            let continuation = try makeCommand("ExactContinuation")
            try encode("reset_continuation", count: 1, into: continuation, uniforms: uniforms)
            try encode("exact_advance", count: arena.exactWork, into: continuation, uniforms: uniforms)
            try await complete(continuation)
            exactDispatches += 1
            status = arena.status()
        }
        try Task.checkCancellation()
        guard status.complete else {
            return certificate(transactionID: transactionID,
                               disposition: status.valid ? .exactWorkBudgetExceeded : .rejected,
                               dt: dt, exactDispatches: exactDispatches, status: status, publications: [])
        }
        let final = try makeCommand("ValidatePublish")
        try encode("validate", count: arena.stateElements, into: final, uniforms: uniforms)
        if !indices.isEmpty {
            try encode("publish", count: indices.count, into: final, uniforms: uniforms, publicationCount: UInt32(indices.count))
        }
        try await complete(final)
        try Task.checkCancellation()
        status = arena.status()
        let outputs = status.complete ? arena.publications(requests: publications) : []
        let result = certificate(transactionID: transactionID,
                                 disposition: status.complete ? .prepared : .rejected,
                                 dt: dt, exactDispatches: exactDispatches, status: status, publications: outputs)
        if result.canCommit { pending = result }
        return result
    }

    public func commitPreparedStep(transactionID: UUID) throws -> VivoHybridStepCertificate {
        guard !busy, let candidate = pending, candidate.transactionID == transactionID, candidate.canCommit else {
            throw VivoHybridExecutionError.missingTransaction
        }
        arena.commit()
        absoluteTime += Double(candidate.requestedTimeStep)
        nextStep += 1
        pending = nil
        return .init(transactionID: transactionID, disposition: .committed,
                     modelFingerprint: candidate.modelFingerprint, planFingerprint: candidate.planFingerprint,
                     stepIndex: candidate.stepIndex, timeBefore: candidate.timeBefore, timeAfter: absoluteTime,
                     requestedTimeStep: candidate.requestedTimeStep, exactDispatches: candidate.exactDispatches,
                     status: candidate.status, publications: candidate.publications)
    }

    public func discardPreparedStep(transactionID: UUID) throws {
        guard !busy, pending?.transactionID == transactionID else { throw VivoHybridExecutionError.missingTransaction }
        pending = nil
    }

    public func step(deltaTime: Float? = nil,
                     publications: [VivoHybridPublicationRequest] = []) async throws -> VivoHybridStepCertificate {
        let candidate = try await prepareStep(deltaTime: deltaTime, publications: publications)
        guard candidate.canCommit else { return candidate }
        // No suspension between returning from preparation and flipping owners.
        return try commitPreparedStep(transactionID: candidate.transactionID)
    }

    public func snapshot() async throws -> VivoHybridStateSnapshot {
        try requireIdle()
        busy = true
        defer { busy = false }
        let bytes = arena.stateBytes * 2
        try arena.checkAdditionalAllocation(bytes)
        guard let readback = device.makeBuffer(length: bytes, options: .storageModeShared),
              let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoHybridExecutionError.metal("snapshot readback allocation failed")
        }
        blit.copy(from: arena.buffer(0), sourceOffset: 0, to: readback, destinationOffset: 0, size: arena.stateBytes)
        blit.copy(from: arena.buffer(4), sourceOffset: 0, to: readback, destinationOffset: arena.stateBytes, size: arena.stateBytes)
        blit.endEncoding()
        try await complete(command)
        let continuous = Array(UnsafeBufferPointer(start: readback.contents().assumingMemoryBound(to: Float.self), count: arena.stateElements))
        let counts = Array(UnsafeBufferPointer(start: readback.contents().advanced(by: arena.stateBytes).assumingMemoryBound(to: UInt32.self), count: arena.stateElements))
        return .init(modelFingerprint: execution.modelFingerprint, planFingerprint: execution.planFingerprint,
                     stepIndex: nextStep, timeSeconds: absoluteTime, laneCount: configuration.laneCount,
                     authorities: execution.authorities, continuousValues: continuous, counts: counts)
    }

    public func checkpoint() async throws -> VivoHybridCheckpoint {
        let snapshot = try await snapshot()
        return .init(modelFingerprint: execution.modelFingerprint, planFingerprint: execution.planFingerprint,
                     seed: configuration.seed, snapshot: snapshot)
    }

    /// Restores into a replacement arena, then changes ownership only after all
    /// uploads complete. A failed allocation/upload leaves the old arena intact.
    public func restore(_ checkpoint: VivoHybridCheckpoint) throws {
        try requireIdle()
        guard checkpoint.schemaVersion == 1, checkpoint.numericalABIVersion == 1,
              checkpoint.modelFingerprint == execution.modelFingerprint,
              checkpoint.planFingerprint == execution.planFingerprint,
              checkpoint.seed == configuration.seed,
              checkpoint.laneCount == configuration.laneCount,
              checkpoint.speciesCount == UInt32(execution.species.count),
              checkpoint.timeSeconds.isFinite, checkpoint.timeSeconds >= 0,
              checkpoint.continuousFP32LE.count == arena.stateBytes,
              checkpoint.countsUInt32LE.count == arena.stateBytes else {
            throw VivoHybridExecutionError.incompatibleCheckpoint("identity, version, seed, clock, or shape mismatch")
        }
        let continuous = try VivoLittleEndianFP32.decode(checkpoint.continuousFP32LE)
        let counts = try VivoHybridCheckpoint.decodeCounts(checkpoint.countsUInt32LE)
        for i in continuous.indices {
            let continuousOwner = execution.authorities[i / Int(configuration.laneCount)] == .deterministicRK2
            guard continuous[i].isFinite, continuous[i] >= 0,
                  continuousOwner ? counts[i] == 0 : continuous[i] == 0 else {
                throw VivoHybridExecutionError.incompatibleCheckpoint("invalid state or noncanonical inactive representation")
            }
        }
        let replacement = try VivoHybridGPUArena(device: device, queue: queue, execution: execution,
                                                  configuration: configuration, continuous: continuous, counts: counts)
        arena = replacement
        nextStep = checkpoint.stepIndex
        absoluteTime = checkpoint.timeSeconds
    }

    private func requireIdle() throws {
        if let failure { throw VivoHybridExecutionError.metal(failure) }
        guard !busy, pending == nil else { throw VivoHybridExecutionError.transactionConflict }
    }
    private func command(dt: Float) -> VivoHybridGPUCommand {
        .init(laneCount: configuration.laneCount, speciesCount: UInt32(execution.species.count),
              reactionCount: UInt32(execution.reactions.count), exactCohortCount: UInt32(execution.exactCohorts.count),
              stepLow: UInt32(truncatingIfNeeded: nextStep), stepHigh: UInt32(truncatingIfNeeded: nextStep >> 32),
              stage: 0, eventsPerDispatch: configuration.exactEventsPerDispatch, dt: dt, seed: configuration.seed)
    }
    private func certificate(transactionID: UUID, disposition: VivoHybridStepDisposition, dt: Float,
                             exactDispatches: UInt32, status: VivoHybridGPUStatus,
                             publications: [VivoHybridPublication]) -> VivoHybridStepCertificate {
        .init(transactionID: transactionID, disposition: disposition,
              modelFingerprint: execution.modelFingerprint, planFingerprint: execution.planFingerprint,
              stepIndex: nextStep, timeBefore: absoluteTime, timeAfter: absoluteTime,
              requestedTimeStep: dt, exactDispatches: exactDispatches, status: status, publications: publications)
    }
    private func makeCommand(_ label: String) throws -> MTLCommandBuffer {
        guard let command = queue.makeCommandBuffer() else { throw VivoHybridExecutionError.metal("command buffer unavailable") }
        command.label = "NumiVivo.Hybrid.\(nextStep).\(label)"
        return command
    }
    private func encode(_ name: String, count: Int, into command: MTLCommandBuffer,
                        uniforms: VivoHybridGPUCommand, publicationCount: UInt32 = 0) throws {
        guard count > 0 else { return }
        let pipeline = try pipelines.pipeline(name)
        guard let encoder = command.makeComputeCommandEncoder() else { throw VivoHybridExecutionError.metal("compute encoder unavailable") }
        defer { encoder.endEncoding() }
        encoder.label = "NumiVivo.Hybrid.\(name)"
        encoder.setComputePipelineState(pipeline)
        for index in 0...17 { encoder.setBuffer(arena.buffer(index), offset: 0, index: index) }
        encoder.setBuffer(arena.buffer(19), offset: 0, index: 19)
        encoder.setBuffer(arena.buffer(20), offset: 0, index: 20)
        var copy = uniforms
        encoder.setBytes(&copy, length: MemoryLayout<VivoHybridGPUCommand>.stride, index: 18)
        var publicationCount = publicationCount
        encoder.setBytes(&publicationCount, length: 4, index: 21)
        let width = pipeline.threadExecutionWidth
        let threads = max(1, min(pipeline.maxTotalThreadsPerThreadgroup, width * 4))
        encoder.dispatchThreads(MTLSize(width: count, height: 1, depth: 1),
                                threadsPerThreadgroup: MTLSize(width: threads, height: 1, depth: 1))
    }
    private func complete(_ command: MTLCommandBuffer) async throws {
        do {
            try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
                command.addCompletedHandler { completed in
                    if let error = completed.error {
                        continuation.resume(throwing: VivoHybridExecutionError.metal(String(describing: error)))
                    } else { continuation.resume(returning: ()) }
                }
                command.commit()
            }
        } catch {
            failure = String(describing: error)
            throw error
        }
    }
}

private final class VivoHybridGPUArena: @unchecked Sendable {
    let device: MTLDevice
    let heap: MTLHeap
    let stateElements: Int
    let stateBytes: Int
    let reactionWork: Int
    let exactWork: Int
    private var buffers: [Int: MTLBuffer]
    private let workingSetFraction: Double

    init(device: MTLDevice, queue: MTLCommandQueue, execution: VivoCompiledHybridExecution,
         configuration: VivoHybridRuntimeConfiguration, continuous: [Float], counts: [UInt32]) throws {
        guard MemoryLayout<VivoHybridGPUReaction>.stride == 48,
              MemoryLayout<VivoHybridGPUChange>.stride == 8,
              MemoryLayout<VivoHybridGPUCohort>.stride == 8,
              MemoryLayout<VivoHybridGPUCommand>.stride == 48 else {
            throw VivoHybridExecutionError.invalidPlan("Swift/Metal hybrid ABI layout mismatch")
        }
        let lanes = UInt64(configuration.laneCount)
        func work(_ count: Int) throws -> Int {
            let result = UInt64(count).multipliedReportingOverflow(by: lanes)
            guard !result.overflow, result.partialValue <= UInt64(UInt32.max) else {
                throw VivoHybridExecutionError.resourceLimit("dispatch grid exceeds UInt32")
            }
            return Int(result.partialValue)
        }
        stateElements = try work(execution.species.count)
        stateBytes = stateElements * 4
        reactionWork = try work(execution.reactions.count)
        exactWork = try work(execution.exactCohorts.count)
        guard continuous.count == stateElements, counts.count == stateElements else {
            throw VivoHybridExecutionError.invalidState("arena input shape mismatch")
        }
        func bytes<T>(_ values: [T]) -> Data { values.withUnsafeBytes { Data($0) } }
        let payloads: [Int: Data] = [
            0: bytes(continuous), 4: bytes(counts), 6: bytes(execution.reactions),
            7: bytes(execution.changes), 8: bytes(execution.incidenceOffsets),
            9: bytes(execution.incidence), 10: bytes(execution.authorities.map(\.rawValue)),
            11: bytes(execution.reactionAuthorities.map(\.rawValue)),
            14: bytes(execution.exactCohorts), 15: bytes(execution.exactReactionIndices)
        ]
        var sizes: [Int: Int] = [0: stateBytes, 1: stateBytes, 2: stateBytes, 3: stateBytes,
                                 4: stateBytes, 5: stateBytes, 12: max(4, reactionWork * 4),
                                 13: max(4, reactionWork * 4), 16: max(16, exactWork * 16)]
        for (index, payload) in payloads { sizes[index] = max(payload.count, 16) }
        var heapBytes = 0
        for index in sizes.keys.sorted() {
            let length = sizes[index]!
            guard length <= device.maxBufferLength else { throw VivoHybridExecutionError.resourceLimit("buffer \(index) exceeds maxBufferLength") }
            let requirement = device.heapBufferSizeAndAlign(length: length, options: .storageModePrivate)
            let alignment = requirement.align
            guard alignment > 0 else { throw VivoHybridExecutionError.resourceLimit("invalid Metal heap alignment") }
            let padding = (alignment - heapBytes % alignment) % alignment
            let sum = heapBytes.addingReportingOverflow(padding)
            let next = sum.partialValue.addingReportingOverflow(requirement.size)
            guard !sum.overflow, !next.overflow else { throw VivoHybridExecutionError.resourceLimit("heap size overflow") }
            heapBytes = next.partialValue
        }
        // A small allocation margin covers heap placement variation; the global
        // working-set fraction remains a hard pre-allocation bound.
        let expanded = Double(heapBytes) * 1.05
        guard expanded < Double(Int.max) else { throw VivoHybridExecutionError.resourceLimit("heap headroom overflow") }
        heapBytes = Int(expanded.rounded(.up))
        let publicationBytes = Int(configuration.maximumPublications) * 16
        let requestBytes = Int(configuration.maximumPublications) * 4
        let stagingBytes = payloads.values.reduce(0) { $0 + $1.count }
        let additional = UInt64(heapBytes) + UInt64(publicationBytes) + UInt64(requestBytes) + 32 + UInt64(stagingBytes)
        let budget = UInt64((Double(device.recommendedMaxWorkingSetSize) * configuration.workingSetFraction).rounded(.down))
        guard publicationBytes <= device.maxBufferLength, requestBytes <= device.maxBufferLength,
              additional <= budget, UInt64(device.currentAllocatedSize) <= budget - additional else {
            throw VivoHybridExecutionError.resourceLimit("heap, shared output and initialization staging exceed the working-set budget")
        }
        let descriptor = MTLHeapDescriptor()
        descriptor.size = heapBytes
        descriptor.storageMode = .private
        descriptor.hazardTrackingMode = .tracked
        descriptor.label = "NumiVivo.Hybrid.PrivateHeap"
        guard let heap = device.makeHeap(descriptor: descriptor) else { throw VivoHybridExecutionError.metal("private heap allocation failed") }
        var buffers: [Int: MTLBuffer] = [:]
        for index in sizes.keys.sorted() {
            guard let buffer = heap.makeBuffer(length: sizes[index]!, options: .storageModePrivate) else {
                throw VivoHybridExecutionError.metal("private buffer \(index) allocation failed")
            }
            buffer.label = "NumiVivo.Hybrid.Buffer.\(index)"
            buffers[index] = buffer
        }
        for (index, length) in [(17, 32), (19, requestBytes), (20, publicationBytes)] {
            guard let buffer = device.makeBuffer(length: length, options: .storageModeShared) else { throw VivoHybridExecutionError.metal("shared buffer allocation failed") }
            buffers[index] = buffer
        }
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else { throw VivoHybridExecutionError.metal("initialization command unavailable") }
        for index in buffers.keys.sorted() { blit.fill(buffer: buffers[index]!, range: 0..<buffers[index]!.length, value: 0) }
        var staging: [MTLBuffer] = []
        for index in payloads.keys.sorted() {
            let data = payloads[index]!
            if data.isEmpty { continue }
            let buffer = data.withUnsafeBytes { raw in
                device.makeBuffer(bytes: raw.baseAddress!, length: raw.count, options: .storageModeShared)
            }
            guard let buffer else { blit.endEncoding(); throw VivoHybridExecutionError.metal("initialization staging allocation failed") }
            staging.append(buffer)
            blit.copy(from: buffer, sourceOffset: 0, to: buffers[index]!, destinationOffset: 0, size: data.count)
        }
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        withExtendedLifetime(staging) {}
        if let error = command.error { throw VivoHybridExecutionError.metal(String(describing: error)) }
        self.device = device
        self.heap = heap
        self.buffers = buffers
        workingSetFraction = configuration.workingSetFraction
    }

    func buffer(_ index: Int) -> MTLBuffer { buffers[index]! }
    func commit() {
        let oldFloat = buffers[0]!
        buffers[0] = buffers[2]
        buffers[2] = oldFloat
        let oldCount = buffers[4]!
        buffers[4] = buffers[5]
        buffers[5] = oldCount
    }
    func checkAdditionalAllocation(_ bytes: Int) throws {
        let budget = UInt64((Double(device.recommendedMaxWorkingSetSize) * workingSetFraction).rounded(.down))
        guard bytes > 0, bytes <= device.maxBufferLength, UInt64(bytes) <= budget,
              UInt64(device.currentAllocatedSize) <= budget - UInt64(bytes) else {
            throw VivoHybridExecutionError.resourceLimit("readback exceeds remaining Metal working-set budget")
        }
    }
    func writePublicationIndices(_ indices: [UInt32]) throws {
        guard indices.count * 4 <= buffer(19).length else { throw VivoHybridExecutionError.resourceLimit("publication request overflow") }
        if !indices.isEmpty {
            indices.withUnsafeBytes { buffer(19).contents().copyMemory(from: $0.baseAddress!, byteCount: $0.count) }
        }
    }
    func status() -> VivoHybridGPUStatus {
        let p = buffer(17).contents().assumingMemoryBound(to: UInt32.self)
        return .init(flags: p[0], firstInvalidLane: p[1] == .max ? nil : p[1],
                     firstInvalidReaction: p[2] == .max ? nil : p[2],
                     unfinishedExactLanes: p[3], maximumExactEvents: p[4])
    }
    func publications(requests: [VivoHybridPublicationRequest]) -> [VivoHybridPublication] {
        let p = buffer(20).contents().assumingMemoryBound(to: UInt32.self)
        return requests.enumerated().map { i, request in
            let mode = VivoHybridGPUAuthority(rawValue: p[i * 4 + 2])!
            let count: UInt32? = mode == .deterministicRK2 ? nil : p[i * 4 + 1]
            return .init(speciesIndex: request.speciesIndex, laneIndex: request.laneIndex,
                         authority: mode, value: count.map(Double.init) ?? Double(Float(bitPattern: p[i * 4])),
                         exactCount: count)
        }
    }
}
