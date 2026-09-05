import CryptoKit
import Foundation
import Metal
import NumiVivoShaders

public struct VivoMetalTargetScreenDescription: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let authoritativeLikelihoodFingerprint: VivoFingerprint
    public let configuration: VivoMetalTargetScreenConfiguration
    public let executableFingerprint: VivoFingerprint
    public let shaderFingerprints: [String: VivoFingerprint]
    public let operatingSystem: String
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let cases: [VivoMetalTargetCaseDescription]
    public let retainedBufferBytes: Int
    public let method: String
    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

private struct VivoScreenInputCommand {
    var lanes: UInt32
    var drugSpecies: UInt32
    var competitorSpecies: UInt32
    var baseExposureM: Float
}
private struct VivoScreenObservationCommand {
    var targetSpecies: SIMD4<UInt32>
    var data: SIMD4<Float>
    var shape: SIMD4<UInt32>
}

/// A reusable OFFLINE likelihood workspace over existing ProgramPack F1 kernels.
/// No generic simulation state is committed here, and no second reaction law is
/// implemented. One arena per calibration case, not per particle. All parameter,
/// input, assay and score buffers persist across SMC proposal batches.
public actor VivoMetalTargetLikelihoodScreen {
    public nonisolated let description: VivoMetalTargetScreenDescription
    public nonisolated let screeningPolicy: VivoPosteriorScreeningPolicy
    private let prepared: VivoPreparedTargetPosterior
    private let policy: VivoMetalTargetScreenConfiguration
    private let queue: MTLCommandQueue
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let exposePipeline: MTLComputePipelineState
    private let observePipeline: MTLComputePipelineState
    private let storage: [CaseStorage]
    private var inFlight = false
    private var poisoned = false

    private final class CaseStorage: @unchecked Sendable {
        let plan: VivoMetalTargetCasePlan
        let arena: VivoMetalArena
        let initial: MTLBuffer
        let parameters: MTLBuffer
        let inputs: MTLBuffer
        let noise: MTLBuffer
        let scores: MTLBuffer
        let failures: MTLBuffer
        let retainedBytes: Int
        init(plan: VivoMetalTargetCasePlan, device: MTLDevice, queue: MTLCommandQueue) throws {
            self.plan = plan
            let n = plan.description.capacity
            arena = try VivoMetalArena(device: device, commandQueue: queue, pack: plan.pack,
                configuration: plan.configuration, couplingCapacity: n, publicationCapacity: 1)
            func shared(_ length: Int, _ label: String) throws -> MTLBuffer {
                guard let buffer = device.makeBuffer(length: max(length, 4), options: .storageModeShared) else {
                    throw VivoRuntimeError.allocationFailed(label)
                }
                buffer.label = "NumiVivo.TargetScreen." + label
                return buffer
            }
            initial = try shared(arena.capacities.stateElements * 4, "Initial")
            parameters = try shared(plan.parameterNames.count * n * 4, "Parameters")
            inputs = try shared(n * MemoryLayout<SIMD2<Float>>.stride, "Inputs")
            noise = try shared(n * MemoryLayout<SIMD4<Float>>.stride, "Noise")
            scores = try shared(n * MemoryLayout<SIMD2<Float>>.stride, "Scores")
            failures = try shared(n * 4, "Failures")
            let values = try plan.pack.initialState(laneCount: n)
            values.withUnsafeBytes { raw in initial.contents().copyMemory(from: raw.baseAddress!, byteCount: raw.count) }
            retainedBytes = arena.heap.size + [arena.commandBuffer, arena.statusBuffer, arena.eventBuffer,
                arena.couplingBuffer, arena.publicationRequestBuffer, arena.publicationOutputBuffer,
                arena.stateReadbackBuffer, initial, parameters, inputs, noise, scores, failures].reduce(0) { $0 + $1.length }
        }
    }

    public static func make(problem: VivoTargetPosteriorProblem,
                            configuration: VivoMetalTargetScreenConfiguration = .init(),
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoMetalTargetLikelihoodScreen {
        try configuration.validate()
        let prepared = try VivoPreparedTargetPosterior(problem)
        guard (1...32).contains(prepared.calibrationCases.count),
              MemoryLayout<VivoScreenInputCommand>.stride == 16,
              MemoryLayout<VivoScreenObservationCommand>.stride == 48 else {
            throw VivoPosteriorError.invalid("screen case capacity or Swift/Metal ABI mismatch")
        }
        let plans = try prepared.calibrationCases.map { try VivoMetalTargetCasePlan.make($0, prepared: prepared, policy: configuration) }
        guard plans.reduce(0, { $0 + $1.description.steps }) <= configuration.maximumTotalPlannedSteps else {
            throw VivoPosteriorError.budget("aggregate fixed-step screen budget")
        }
        // Conservative preallocation allowance for the fixed small target graph.
        // Driver/pipeline allocation is not included in the retained-buffer cap.
        let estimated = plans.reduce(0) { total, plan in
            total + 8 * plan.pack.data.count + 2_097_152
                + plan.description.capacity * (Int(plan.pack.runtimeContract.speciesCount) * 64
                    + plan.parameterNames.count * 16 + 512)
        }
        guard estimated <= configuration.maximumRetainedBufferBytes else {
            throw VivoPosteriorError.budget("estimated retained Metal workspace exceeds configured buffer allowance")
        }
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoRuntimeError.incompatibleDevice("screen requires unified-memory Metal and a command queue")
        }
        queue.label = "NumiVivo.TargetLikelihoodScreen"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in [NumiVivoKernel.clearStatus, .f1HeunPredict, .f1HeunCorrect, .validateShadow, .evaluateMonitors] {
            pipelines[kernel] = try await catalog.pipeline(kernel)
        }
        var hashes: [String: VivoFingerprint] = [:], scoringSource = ""
        for name in ["NumiVivoTargetLikelihood", "NumiVivoProgramPackRuntime"] {
            let bundle = NumiVivoShaderResources.bundle
            guard let url = bundle.url(forResource: name, withExtension: "metal")
                ?? bundle.url(forResource: name, withExtension: "metal", subdirectory: "Resources") else {
                throw NumiVivoShaderError.sourceResourceMissing
            }
            let bytes = try Data(contentsOf: url)
            guard bytes.count <= 4_194_304, let text = String(data: bytes, encoding: .utf8) else {
                throw VivoPosteriorError.invalid("shader resource representation")
            }
            hashes[name] = try VivoCanonicalJSON.fingerprint(bytes)
            if name == "NumiVivoTargetLikelihood" { scoringSource = text }
        }
        let options = MTLCompileOptions(); options.fastMathEnabled = false
        let library = try device.makeLibrary(source: scoringSource, options: options)
        guard let expose = library.makeFunction(name: "nvivo_target_likelihood_expose"),
              let observe = library.makeFunction(name: "nvivo_target_likelihood_observe") else {
            throw NumiVivoShaderError.functionMissing("target likelihood functions")
        }
        let exposePipeline = try device.makeComputePipelineState(function: expose)
        let observePipeline = try device.makeComputePipelineState(function: observe)
        var storage: [CaseStorage] = [], bytes = 0
        for plan in plans {
            let item = try CaseStorage(plan: plan, device: device, queue: queue)
            bytes += item.retainedBytes
            guard bytes <= configuration.maximumRetainedBufferBytes else {
                throw VivoPosteriorError.budget("retained Metal buffers exceed configured allowance")
            }
            storage.append(item)
        }
        let description = try VivoMetalTargetScreenDescription(schemaVersion: 1,
            authoritativeLikelihoodFingerprint: prepared.plan.likelihoodFingerprint, configuration: configuration,
            executableFingerprint: executableFingerprint(), shaderFingerprints: hashes,
            operatingSystem: ProcessInfo.processInfo.operatingSystemVersionString, deviceName: device.name,
            deviceRegistryID: device.registryID, cases: plans.map(\.description), retainedBufferBytes: bytes,
            method: "numivivo.f1-prior-fixed-independent-gaussian-screen.v1; fp32; fastMath=false")
        return try .init(description: description, prepared: prepared, policy: configuration, queue: queue,
                         pipelines: pipelines, expose: exposePipeline, observe: observePipeline, storage: storage)
    }

    private init(description: VivoMetalTargetScreenDescription, prepared: VivoPreparedTargetPosterior,
                 policy: VivoMetalTargetScreenConfiguration, queue: MTLCommandQueue,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline], expose: MTLComputePipelineState,
                 observe: MTLComputePipelineState, storage: [CaseStorage]) throws {
        self.description = description; self.prepared = prepared; self.policy = policy; self.queue = queue
        self.pipelines = pipelines; exposePipeline = expose; observePipeline = observe; self.storage = storage
        screeningPolicy = try .init(fingerprint: description.fingerprint(), maximumEvaluations: policy.maximumScreeningEvaluations)
    }

    public nonisolated func screen() -> VivoPosteriorScreen {
        .init(policy: screeningPolicy, evaluate: { try await self.evaluate($0) })
    }

    public func evaluate(_ candidates: [VivoPosteriorCandidate]) async throws -> [VivoPosteriorEvaluation] {
        guard !inFlight else { throw VivoPosteriorError.busy }
        guard !poisoned, candidates.count <= prepared.plan.configuration.particleCount,
              Set(candidates.map(\.ordinal)).count == candidates.count else {
            throw VivoPosteriorError.invalid("screen workspace failed, batch too large or duplicate candidate identity")
        }
        if candidates.isEmpty { return [] }
        for candidate in candidates {
            let expected = try VivoPosteriorCandidate(ordinal: candidate.ordinal,
                coordinates: candidate.coordinates, parameters: prepared.plan.parameters)
            guard expected == candidate else { throw VivoPosteriorError.invalid("screen candidate values/coordinates mismatch") }
        }
        inFlight = true
        defer { inFlight = false }
        var scores = [Double](repeating: 0, count: candidates.count), compensation = scores
        for item in storage {
            try Task.checkCancellation()
            try fill(item, candidates: candidates)
            let values = try await execute(item)
            for index in candidates.indices {
                let corrected = values[index] - compensation[index], next = scores[index] + corrected
                compensation[index] = (next - scores[index]) - corrected; scores[index] = next
                guard next.isFinite else { throw VivoPosteriorError.numerical("screen likelihood sum overflow") }
            }
        }
        return candidates.indices.map { .init(candidate: candidates[$0], logLikelihood: scores[$0]) }
    }

    private func fill(_ item: CaseStorage, candidates: [VivoPosteriorCandidate]) throws {
        let plan = item.plan, n = plan.description.capacity
        let parameters = item.parameters.contents().bindMemory(to: Float.self, capacity: plan.parameterNames.count * n)
        let inputs = item.inputs.contents().bindMemory(to: SIMD2<Float>.self, capacity: n)
        let noise = item.noise.contents().bindMemory(to: SIMD4<Float>.self, capacity: n)
        let bounds = try plan.pack.parameterMetadata()
        for lane in candidates.indices {
            let values = candidates[lane].values
            let experiment = try prepared.experiment(for: plan.item, values: values)
            let m = experiment.kinetics
            let table: [String: Double] = ["kon": m.association.value, "koff": m.dissociation.value,
                "kinact": m.inactivation.value, "kturnover": m.targetTurnover.value,
                "competitor-kon": m.competitor?.association.value ?? 0,
                "competitor-koff": m.competitor?.dissociation.value ?? 0]
            for (column, name) in plan.parameterNames.enumerated() {
                guard let value = table[name] else { throw VivoPosteriorError.invalid("screen parameter mapping") }
                let narrowed = try vivoScreenFloat(value, label: name)
                guard Double(narrowed) >= bounds[column].minimum, Double(narrowed) <= bounds[column].maximum else {
                    throw VivoPosteriorError.invalid("candidate outside compiled numerical bounds")
                }
                parameters[column * n + lane] = narrowed
            }
            let scaleBinding = prepared.problem.bindings.firstIndex { $0.field == .exposureScale && $0.caseIdentifiers.contains(plan.item.identifier) }
            let scale = try vivoScreenFloat(scaleBinding.map { values[$0] } ?? 1, label: "exposure scale")
            for knot in plan.item.experiment.exposure.knots {
                let result = Float(knot.unboundDrugM) * scale
                guard result.isFinite, result == 0 ? knot.unboundDrugM == 0 || scale == 0 : result >= Float.leastNormalMagnitude else {
                    throw VivoPosteriorError.invalid("FP32 exposure product overflow/underflow")
                }
            }
            inputs[lane] = SIMD2(scale, try vivoScreenFloat(m.competitor?.unboundConcentration.value ?? 0, label: "competitor concentration"))
            let assay = try prepared.problem.assayNoise(caseIdentifier: plan.item.identifier, values: values)
            noise[lane] = SIMD4(try vivoScreenFloat(assay.scale, label: "noise scale"),
                try vivoScreenFloat(assay.additionalStandardDeviation, label: "extra SD"),
                try vivoScreenFloat(assay.bias, label: "assay bias"), 0)
        }
        // Padding is computational only, never extra statistical particles.
        // Duplicate lane zero so inactive capacity cannot introduce new failures.
        for lane in candidates.count..<n {
            for column in plan.parameterNames.indices { parameters[column * n + lane] = parameters[column * n] }
            inputs[lane] = inputs[0]; noise[lane] = noise[0]
        }
    }

    private func execute(_ item: CaseStorage) async throws -> [Double] {
        let plan = item.plan, arena = item.arena, n = plan.description.capacity
        guard let reset = queue.makeCommandBuffer(), let blit = reset.makeBlitCommandEncoder() else {
            throw VivoRuntimeError.commandQueueUnavailable
        }
        blit.copy(from: item.initial, sourceOffset: 0, to: arena.baseState, destinationOffset: 0, size: item.initial.length)
        blit.copy(from: item.parameters, sourceOffset: 0, to: arena.parameterBuffer, destinationOffset: 0,
                  size: plan.parameterNames.count * n * 4)
        for buffer in [arena.temporalCandidate, item.scores, item.failures] {
            blit.fill(buffer: buffer, range: 0..<buffer.length, value: 0)
        }
        blit.endEncoding()
        try encode(.clearStatus, count: 1, command: reset, buffers: [arena.statusBuffer])
        try await complete(reset)
        var current = arena.baseState, next = arena.candidateState, stepIndex: UInt32 = 0
        var offset = 0
        while offset < plan.operations.count {
            try Task.checkCancellation()
            guard let command = queue.makeCommandBuffer() else { throw VivoRuntimeError.commandQueueUnavailable }
            command.label = "NumiVivo.TargetScreen.\(plan.item.identifier).\(offset)"
            let limit = min(offset + policy.operationsPerCommand, plan.operations.count)
            for op in plan.operations[offset..<limit] {
                switch op {
                case .exposure(let base):
                    var input = VivoScreenInputCommand(lanes: UInt32(n), drugSpecies: plan.drugIndex,
                        competitorSpecies: plan.competitorIndex ?? UInt32.max, baseExposureM: base)
                    try custom(exposePipeline, command: command, buffers: [current, item.inputs], uniform: &input, count: n)
                    let uniform = try VivoRuntimeCommandABI(pack: plan.pack, configuration: plan.configuration,
                        stepIndex: stepIndex, substepIndex: 0, dt: 0, absoluteTime: 0)
                    try encode(.validateShadow, count: arena.capacities.stateElements, command: command,
                               buffers: [current, arena.programBuffer, arena.statusBuffer], uniform: uniform)
                case .step(let dt, let after):
                    let uniform = try VivoRuntimeCommandABI(pack: plan.pack, configuration: plan.configuration,
                        stepIndex: stepIndex, substepIndex: 0, dt: dt, absoluteTime: after)
                    try encode(.f1HeunPredict, count: arena.capacities.stateElements, command: command,
                        buffers: [current, arena.stageState, arena.derivativeK1, arena.parameterBuffer, arena.programBuffer, arena.statusBuffer], uniform: uniform)
                    try encode(.f1HeunCorrect, count: arena.capacities.stateElements, command: command,
                        buffers: [current, arena.stageState, arena.derivativeK1, next, arena.parameterBuffer, arena.programBuffer, arena.statusBuffer], uniform: uniform)
                    try encode(.evaluateMonitors, count: n, command: command,
                        buffers: [next, arena.temporalCandidate, arena.parameterBuffer, arena.programBuffer, arena.eventBuffer, arena.statusBuffer], uniform: uniform)
                    try encode(.validateShadow, count: arena.capacities.stateElements, command: command,
                        buffers: [next, arena.programBuffer, arena.statusBuffer], uniform: uniform)
                    swap(&current, &next); stepIndex += 1
                case .observe(let index):
                    let observation = plan.item.observations[index]
                    let kind: UInt32
                    switch observation.observable {
                    case .drugOccupancy: kind = 0
                    case .covalentOccupancy: kind = 1
                    case .freeTargetFraction: kind = 2
                    case .competitorOccupancy: kind = 3
                    }
                    var uniform = VivoScreenObservationCommand(targetSpecies: plan.targetIndices,
                        data: SIMD4(Float(observation.value), Float(observation.standardDeviation ?? 0), 1e-4, 0),
                        shape: SIMD4(kind, UInt32(n), 0, 0))
                    try custom(observePipeline, command: command,
                        buffers: [current, item.noise, item.scores, item.failures], uniform: &uniform, count: n)
                }
            }
            try await complete(command)
            let status = arena.status()
            guard status.flags.rawValue == 0, status.requestedResponse == 0, status.eventDropped == 0 else {
                throw VivoPosteriorError.numerical("Metal screen state failed: flags=\(status.flags.rawValue), case=\(plan.item.identifier)")
            }
            offset = limit
        }
        let failures = item.failures.contents().assumingMemoryBound(to: UInt32.self)
        let scores = item.scores.contents().assumingMemoryBound(to: SIMD2<Float>.self)
        for lane in 0..<n {
            guard failures[lane] == 0, scores[lane].x.isFinite else {
                throw VivoPosteriorError.numerical("Metal observation failure in lane \(lane); no candidate was filtered out")
            }
        }
        return (0..<n).map { Double(scores[$0].x) }
    }

    private func encode(_ kernel: NumiVivoKernel, count: Int, command: MTLCommandBuffer,
                        buffers: [MTLBuffer], uniform: VivoRuntimeCommandABI? = nil) throws {
        guard let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoRuntimeError.pipelineEncodingFailed(kernel.rawValue)
        }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        if var uniform { encoder.setBytes(&uniform, length: MemoryLayout<VivoRuntimeCommandABI>.stride, index: buffers.count) }
        encoder.dispatchThreads(pipeline.gridSize(for: count), threadsPerThreadgroup: pipeline.threadgroupSize(for: count))
    }
    private func custom<T>(_ pipeline: MTLComputePipelineState, command: MTLCommandBuffer,
                           buffers: [MTLBuffer], uniform: inout T, count: Int) throws {
        guard let encoder = command.makeComputeCommandEncoder() else { throw VivoRuntimeError.commandQueueUnavailable }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        encoder.setBytes(&uniform, length: MemoryLayout<T>.stride, index: buffers.count)
        encoder.dispatchThreads(.init(width: count, height: 1, depth: 1),
            threadsPerThreadgroup: .init(width: max(1, min(pipeline.threadExecutionWidth * 4, pipeline.maxTotalThreadsPerThreadgroup)), height: 1, depth: 1))
    }
    private func complete(_ command: MTLCommandBuffer) async throws {
        do {
            try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
                command.addCompletedHandler { finished in
                    if let error = finished.error { continuation.resume(throwing: VivoRuntimeError.commandFailed(String(describing: error))) }
                    else { continuation.resume(returning: ()) }
                }
                command.commit()
            }
        } catch { poisoned = true; throw error }
        // Cancellation never releases shared CPU buffers before GPU completion.
        try Task.checkCancellation()
    }

    private static func executableFingerprint() throws -> VivoFingerprint {
        guard let url = Bundle.main.executableURL else { throw VivoPosteriorError.invalid("executing binary URL unavailable") }
        let handle = try FileHandle(forReadingFrom: url)
        defer { try? handle.close() }
        var hash = SHA256(), size = 0
        while let data = try handle.read(upToCount: 1_048_576), !data.isEmpty {
            size += data.count
            guard size <= 536_870_912 else { throw VivoPosteriorError.budget("executable identity snapshot size") }
            hash.update(data: data)
        }
        guard size > 0 else { throw VivoPosteriorError.invalid("empty executing-binary snapshot") }
        return try .init(bytes: Array(hash.finalize()))
    }
}
