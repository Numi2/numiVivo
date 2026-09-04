import Foundation
import Metal
import NumiVivoShaders

public struct VivoPartitionRuntimeConfiguration: Codable, Equatable, Sendable {
    public var timeStep: Float
    public var minimumTimeStep: Float
    public var maximumTimeStep: Float
    public var maximumAttempts: UInt32
    public var negativeTolerance: Float

    public init(
        timeStep: Float,
        minimumTimeStep: Float = 1e-8,
        maximumTimeStep: Float = 3_600,
        maximumAttempts: UInt32 = 16,
        negativeTolerance: Float = 1e-9
    ) {
        self.timeStep = timeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.maximumAttempts = maximumAttempts
        self.negativeTolerance = negativeTolerance
    }

    public func validate() throws {
        guard timeStep.isFinite,
              minimumTimeStep.isFinite,
              maximumTimeStep.isFinite,
              negativeTolerance.isFinite,
              minimumTimeStep > 0,
              maximumTimeStep >= minimumTimeStep,
              timeStep >= minimumTimeStep,
              timeStep <= maximumTimeStep,
              maximumAttempts > 0,
              negativeTolerance >= 0 else {
            throw VivoPartitionRuntimeError.invalidConfiguration
        }
    }
}

public struct VivoPartitionStatusFlags: OptionSet, Codable, Sendable {
    public let rawValue: UInt32
    public init(rawValue: UInt32) { self.rawValue = rawValue }

    public static let needsSmallerStep = Self(rawValue: 1 << 0)
    public static let nonFinite = Self(rawValue: 1 << 1)
    public static let negative = Self(rawValue: 1 << 2)
    public static let aboveBound = Self(rawValue: 1 << 3)

    public static let terminalInvalid: Self = [.nonFinite]
}

public struct VivoPartitionRuntimeStatus: Codable, Equatable, Sendable {
    public let flags: VivoPartitionStatusFlags
    public let invalidElementCount: UInt32
    public let firstInvalidElement: UInt32?
    public let maximumConcentration: Float
    public let maximumNegativeMagnitude: Float
}

public enum VivoPartitionStepDisposition: String, Codable, Sendable {
    case committed
    case committedWithReducedStep
    case rejected
}

public struct VivoPartitionStepCertificate: Codable, Equatable, Sendable {
    public let disposition: VivoPartitionStepDisposition
    public let modelFingerprint: String
    public let deviceName: String
    public let registryID: UInt64
    public let stepIndex: UInt32
    public let requestedStep: Float
    public let stabilityBound: Float
    public let acceptedStep: Float?
    public let attemptCount: UInt32
    public let timeBefore: Float
    public let timeAfter: Float
    public let status: VivoPartitionRuntimeStatus
    public let totalAmountBefore: [Double]?
    public let totalAmountAfter: [Double]?
}

public struct VivoPartitionStateSnapshot: Codable, Equatable, Sendable {
    public let modelFingerprint: String
    public let stepIndex: UInt32
    public let absoluteTime: Float
    public let compartmentCount: UInt32
    public let analyteCount: UInt32
    public let concentrations: [Float]

    public func concentration(analyte: UInt32, compartment: UInt32) -> Float? {
        guard analyte < analyteCount, compartment < compartmentCount else { return nil }
        return concentrations[Int(analyte * compartmentCount + compartment)]
    }
}

public enum VivoPartitionRuntimeError: Error, LocalizedError, Sendable {
    case invalidConfiguration
    case invalidInitialState(String)
    case deviceUnavailable
    case commandQueueUnavailable
    case shaderResourceUnavailable
    case shaderCompilation(String)
    case pipelineUnavailable(String)
    case allocationFailed(String)
    case commandFailed(String)
    case arithmeticOverflow

    public var errorDescription: String? {
        switch self {
        case .invalidConfiguration: return "invalid partition runtime configuration"
        case .invalidInitialState(let reason): return "invalid partition initial state: \(reason)"
        case .deviceUnavailable: return "no Metal device is available"
        case .commandQueueUnavailable: return "Metal command queue is unavailable"
        case .shaderResourceUnavailable: return "NumiVivo partition Metal source is unavailable"
        case .shaderCompilation(let reason): return "partition Metal compilation failed: \(reason)"
        case .pipelineUnavailable(let name): return "partition pipeline is unavailable: \(name)"
        case .allocationFailed(let name): return "partition Metal allocation failed: \(name)"
        case .commandFailed(let reason): return "partition Metal command failed: \(reason)"
        case .arithmeticOverflow: return "partition runtime arithmetic overflow"
        }
    }
}

private struct VivoPartitionRuntimeCommandABI {
    var compartmentCount: UInt32
    var analyteCount: UInt32
    var edgeCount: UInt32
    var stepIndex: UInt32
    var deltaTime: Float
    var absoluteTime: Float
    var negativeTolerance: Float
    var attemptIndex: UInt32
}

private struct VivoPartitionRuntimeStatusABI {
    var flags: UInt32
    var invalidElementCount: UInt32
    var firstInvalidElement: UInt32
    var maximumConcentrationBits: UInt32
    var maximumNegativeMagnitudeBits: UInt32
    var reserved0: UInt32
    var reserved1: UInt32
    var reserved2: UInt32
}

/// Transactional reversible-partition runtime. The gather formulation computes
/// equal and opposite amount transfer for each edge and therefore preserves each
/// analyte's total amount up to floating-point integration error.
public actor VivoPhysiologicalPartitionRuntime {
    public nonisolated let model: VivoPhysiologicalPartitionModel
    public nonisolated let configuration: VivoPartitionRuntimeConfiguration
    public nonisolated let deviceName: String
    public nonisolated let registryID: UInt64

    private let device: MTLDevice
    private let queue: MTLCommandQueue
    private let compiled: VivoCompiledPartitionModel
    private let clearPipeline: MTLComputePipelineState
    private let stagePipeline: MTLComputePipelineState
    private let finalizePipeline: MTLComputePipelineState
    private let validatePipeline: MTLComputePipelineState
    private let compartmentBuffer: MTLBuffer
    private let analyteBuffer: MTLBuffer
    private let edgeBuffer: MTLBuffer
    private var stateA: MTLBuffer
    private var stateB: MTLBuffer
    private let stageState: MTLBuffer
    private let firstDerivative: MTLBuffer
    private let statusBuffer: MTLBuffer
    private let readbackBuffer: MTLBuffer
    private let elementCount: Int
    private let byteCount: Int
    private var authoritativeIsA = true
    private var absoluteTime: Float = 0
    private var nextStepIndex: UInt32 = 0

    public static func make(
        model: VivoPhysiologicalPartitionModel,
        configuration: VivoPartitionRuntimeConfiguration,
        initialConcentrations: [Float],
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoPhysiologicalPartitionRuntime {
        try model.validate()
        try configuration.validate()
        let compiled = try VivoCompiledPartitionModel(model: model)
        let elements = model.compartments.count.multipliedReportingOverflow(by: model.analytes.count)
        guard !elements.overflow else { throw VivoPartitionRuntimeError.arithmeticOverflow }
        let bytes = elements.partialValue.multipliedReportingOverflow(by: MemoryLayout<Float>.stride)
        guard !bytes.overflow, bytes.partialValue > 0 else { throw VivoPartitionRuntimeError.arithmeticOverflow }

        let expanded: [Float]
        if initialConcentrations.count == model.analytes.count {
            var values = [Float](repeating: 0, count: elements.partialValue)
            for analyte in model.analytes.indices {
                let value = initialConcentrations[analyte]
                let definition = model.analytes[analyte]
                guard value.isFinite,
                      value >= definition.minimumConcentration,
                      value <= definition.maximumConcentration else {
                    throw VivoPartitionRuntimeError.invalidInitialState("analyte \(analyte) is outside its bounds")
                }
                for compartment in model.compartments.indices {
                    values[analyte * model.compartments.count + compartment] = value
                }
            }
            expanded = values
        } else if initialConcentrations.count == elements.partialValue {
            for analyte in model.analytes.indices {
                let definition = model.analytes[analyte]
                for compartment in model.compartments.indices {
                    let value = initialConcentrations[analyte * model.compartments.count + compartment]
                    guard value.isFinite,
                          value >= definition.minimumConcentration,
                          value <= definition.maximumConcentration else {
                        throw VivoPartitionRuntimeError.invalidInitialState("analyte \(analyte), compartment \(compartment) is outside bounds")
                    }
                }
            }
            expanded = initialConcentrations
        } else {
            throw VivoPartitionRuntimeError.invalidInitialState(
                "expected analyteCount or analyteCount × compartmentCount values"
            )
        }

        let device = try requestedDevice ?? MTLCreateSystemDefaultDevice().orThrowPartition(VivoPartitionRuntimeError.deviceUnavailable)
        guard let queue = device.makeCommandQueue() else {
            throw VivoPartitionRuntimeError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.Partition.Queue"
        guard let sourceURL = NumiVivoShaderResources.bundle.url(
            forResource: "NumiVivoPartitionKernels",
            withExtension: "metal"
        ) else {
            throw VivoPartitionRuntimeError.shaderResourceUnavailable
        }
        let source = String(decoding: try Data(contentsOf: sourceURL), as: UTF8.self)
        let library: MTLLibrary
        do { library = try device.makeLibrary(source: source, options: nil) }
        catch { throw VivoPartitionRuntimeError.shaderCompilation(error.localizedDescription) }

        func pipeline(_ name: String) throws -> MTLComputePipelineState {
            guard let function = library.makeFunction(name: name) else {
                throw VivoPartitionRuntimeError.pipelineUnavailable(name)
            }
            do { return try device.makeComputePipelineState(function: function) }
            catch { throw VivoPartitionRuntimeError.shaderCompilation("\(name): \(error.localizedDescription)") }
        }

        let compartmentBuffer = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.compartments,
            label: "NumiVivo.Partition.Compartments"
        )
        let analyteBuffer = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.analytes,
            label: "NumiVivo.Partition.Analytes"
        )
        let edgeFallback = VivoPartitionEdgeABI(
            analyteIndex: 0,
            sourceCompartment: 0,
            targetCompartment: 0,
            partitionCoefficient: 1,
            clearance: 0,
            sourceUnboundFraction: 1,
            targetUnboundFraction: 1
        )
        let edgeBuffer = try makePrivateBufferAllowingEmpty(
            device: device,
            queue: queue,
            values: compiled.edges,
            fallback: edgeFallback,
            label: "NumiVivo.Partition.Edges"
        )
        let stateA = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: expanded,
            label: "NumiVivo.Partition.StateA"
        )
        guard let stateB = device.makeBuffer(length: bytes.partialValue, options: .storageModePrivate),
              let stage = device.makeBuffer(length: bytes.partialValue, options: .storageModePrivate),
              let derivative = device.makeBuffer(length: bytes.partialValue, options: .storageModePrivate),
              let status = device.makeBuffer(length: MemoryLayout<VivoPartitionRuntimeStatusABI>.stride, options: .storageModeShared),
              let readback = device.makeBuffer(length: bytes.partialValue, options: .storageModeShared) else {
            throw VivoPartitionRuntimeError.allocationFailed("transaction state")
        }
        stateB.label = "NumiVivo.Partition.StateB"
        stage.label = "NumiVivo.Partition.Stage"
        derivative.label = "NumiVivo.Partition.FirstDerivative"
        status.label = "NumiVivo.Partition.Status"
        readback.label = "NumiVivo.Partition.Readback"

        return VivoPhysiologicalPartitionRuntime(
            model: model,
            configuration: configuration,
            device: device,
            queue: queue,
            compiled: compiled,
            clearPipeline: try pipeline("numivivo_partition::nvivo_partition_clear_status"),
            stagePipeline: try pipeline("numivivo_partition::nvivo_partition_stage"),
            finalizePipeline: try pipeline("numivivo_partition::nvivo_partition_finalize"),
            validatePipeline: try pipeline("numivivo_partition::nvivo_partition_validate"),
            compartmentBuffer: compartmentBuffer,
            analyteBuffer: analyteBuffer,
            edgeBuffer: edgeBuffer,
            stateA: stateA,
            stateB: stateB,
            stageState: stage,
            firstDerivative: derivative,
            statusBuffer: status,
            readbackBuffer: readback,
            elementCount: elements.partialValue,
            byteCount: bytes.partialValue
        )
    }

    private init(
        model: VivoPhysiologicalPartitionModel,
        configuration: VivoPartitionRuntimeConfiguration,
        device: MTLDevice,
        queue: MTLCommandQueue,
        compiled: VivoCompiledPartitionModel,
        clearPipeline: MTLComputePipelineState,
        stagePipeline: MTLComputePipelineState,
        finalizePipeline: MTLComputePipelineState,
        validatePipeline: MTLComputePipelineState,
        compartmentBuffer: MTLBuffer,
        analyteBuffer: MTLBuffer,
        edgeBuffer: MTLBuffer,
        stateA: MTLBuffer,
        stateB: MTLBuffer,
        stageState: MTLBuffer,
        firstDerivative: MTLBuffer,
        statusBuffer: MTLBuffer,
        readbackBuffer: MTLBuffer,
        elementCount: Int,
        byteCount: Int
    ) {
        self.model = model
        self.configuration = configuration
        self.device = device
        self.queue = queue
        self.compiled = compiled
        self.clearPipeline = clearPipeline
        self.stagePipeline = stagePipeline
        self.finalizePipeline = finalizePipeline
        self.validatePipeline = validatePipeline
        self.compartmentBuffer = compartmentBuffer
        self.analyteBuffer = analyteBuffer
        self.edgeBuffer = edgeBuffer
        self.stateA = stateA
        self.stateB = stateB
        self.stageState = stageState
        self.firstDerivative = firstDerivative
        self.statusBuffer = statusBuffer
        self.readbackBuffer = readbackBuffer
        self.elementCount = elementCount
        self.byteCount = byteCount
        self.deviceName = device.name
        self.registryID = device.registryID
    }

    public func time() -> Float { absoluteTime }
    public func stepIndex() -> UInt32 { nextStepIndex }
    public func conservativeMaximumStep() -> Float {
        max(configuration.minimumTimeStep, min(configuration.maximumTimeStep, compiled.conservativeMaximumStep))
    }

    public func step(
        deltaTime requestedStep: Float? = nil,
        permitAdaptiveReduction: Bool = true,
        certifyAmountConservation: Bool = false
    ) async throws -> VivoPartitionStepCertificate {
        let requested = requestedStep ?? configuration.timeStep
        guard requested.isFinite,
              requested >= configuration.minimumTimeStep,
              requested <= configuration.maximumTimeStep else {
            throw VivoPartitionRuntimeError.invalidConfiguration
        }
        let stability = conservativeMaximumStep()
        var candidate = requested
        if candidate > stability {
            guard permitAdaptiveReduction else {
                return rejectedCertificate(
                    requested: requested,
                    stability: stability,
                    attempts: 0,
                    status: .init(flags: [.needsSmallerStep], invalidElementCount: 0, firstInvalidElement: nil, maximumConcentration: 0, maximumNegativeMagnitude: 0)
                )
            }
            candidate = stability
        }

        let beforeAmounts = certifyAmountConservation ? try await totalAmounts() : nil
        var last = VivoPartitionRuntimeStatus(
            flags: [],
            invalidElementCount: 0,
            firstInvalidElement: nil,
            maximumConcentration: 0,
            maximumNegativeMagnitude: 0
        )
        for attempt in 0..<configuration.maximumAttempts {
            last = try await executeAttempt(deltaTime: candidate, attempt: attempt)
            if !last.flags.intersection(.terminalInvalid).isEmpty {
                return rejectedCertificate(requested: requested, stability: stability, attempts: attempt + 1, status: last)
            }
            if last.flags.contains(.needsSmallerStep) {
                let reduced = candidate * 0.5
                guard permitAdaptiveReduction,
                      attempt + 1 < configuration.maximumAttempts,
                      reduced >= configuration.minimumTimeStep,
                      reduced < candidate else {
                    return rejectedCertificate(requested: requested, stability: stability, attempts: attempt + 1, status: last)
                }
                candidate = reduced
                continue
            }

            authoritativeIsA.toggle()
            let before = absoluteTime
            absoluteTime += candidate
            let completedStep = nextStepIndex
            nextStepIndex &+= 1
            let afterAmounts = certifyAmountConservation ? try await totalAmounts() : nil
            return .init(
                disposition: candidate == requested ? .committed : .committedWithReducedStep,
                modelFingerprint: model.fingerprint,
                deviceName: deviceName,
                registryID: registryID,
                stepIndex: completedStep,
                requestedStep: requested,
                stabilityBound: stability,
                acceptedStep: candidate,
                attemptCount: attempt + 1,
                timeBefore: before,
                timeAfter: absoluteTime,
                status: last,
                totalAmountBefore: beforeAmounts,
                totalAmountAfter: afterAmounts
            )
        }
        return rejectedCertificate(requested: requested, stability: stability, attempts: configuration.maximumAttempts, status: last)
    }

    public func snapshot() async throws -> VivoPartitionStateSnapshot {
        let values = try await readAuthoritativeState()
        return .init(
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            absoluteTime: absoluteTime,
            compartmentCount: UInt32(model.compartments.count),
            analyteCount: UInt32(model.analytes.count),
            concentrations: values
        )
    }

    public func totalAmounts() async throws -> [Double] {
        let state = try await readAuthoritativeState()
        let compartments = model.compartments.count
        return model.analytes.indices.map { analyte in
            var amount: Double = 0
            for compartment in model.compartments.indices {
                amount += Double(state[analyte * compartments + compartment]) *
                    Double(model.compartments[compartment].volumeCubicMetres)
            }
            return amount
        }
    }

    private var authoritativeState: MTLBuffer { authoritativeIsA ? stateA : stateB }
    private var candidateState: MTLBuffer { authoritativeIsA ? stateB : stateA }

    private func executeAttempt(deltaTime: Float, attempt: UInt32) async throws -> VivoPartitionRuntimeStatus {
        guard let command = queue.makeCommandBuffer() else {
            throw VivoPartitionRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Partition.Step.\(nextStepIndex).Attempt.\(attempt)"
        var runtimeCommand = VivoPartitionRuntimeCommandABI(
            compartmentCount: UInt32(model.compartments.count),
            analyteCount: UInt32(model.analytes.count),
            edgeCount: UInt32(model.edges.count),
            stepIndex: nextStepIndex,
            deltaTime: deltaTime,
            absoluteTime: absoluteTime + deltaTime,
            negativeTolerance: configuration.negativeTolerance,
            attemptIndex: attempt
        )

        try encode(pipeline: clearPipeline, commandBuffer: command, count: 1, label: "Partition.Clear") { encoder in
            encoder.setBuffer(statusBuffer, offset: 0, index: 0)
        }
        try encode(pipeline: stagePipeline, commandBuffer: command, count: elementCount, label: "Partition.Stage") { encoder in
            encoder.setBuffer(authoritativeState, offset: 0, index: 0)
            encoder.setBuffer(stageState, offset: 0, index: 1)
            encoder.setBuffer(firstDerivative, offset: 0, index: 2)
            encoder.setBuffer(compartmentBuffer, offset: 0, index: 3)
            encoder.setBuffer(edgeBuffer, offset: 0, index: 4)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPartitionRuntimeCommandABI>.stride, index: 5)
        }
        try encode(pipeline: finalizePipeline, commandBuffer: command, count: elementCount, label: "Partition.Finalize") { encoder in
            encoder.setBuffer(authoritativeState, offset: 0, index: 0)
            encoder.setBuffer(stageState, offset: 0, index: 1)
            encoder.setBuffer(firstDerivative, offset: 0, index: 2)
            encoder.setBuffer(candidateState, offset: 0, index: 3)
            encoder.setBuffer(compartmentBuffer, offset: 0, index: 4)
            encoder.setBuffer(edgeBuffer, offset: 0, index: 5)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPartitionRuntimeCommandABI>.stride, index: 6)
        }
        try encode(pipeline: validatePipeline, commandBuffer: command, count: elementCount, label: "Partition.Validate") { encoder in
            encoder.setBuffer(candidateState, offset: 0, index: 0)
            encoder.setBuffer(analyteBuffer, offset: 0, index: 1)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPartitionRuntimeCommandABI>.stride, index: 2)
            encoder.setBuffer(statusBuffer, offset: 0, index: 3)
        }
        try await complete(command)

        let raw = statusBuffer.contents().load(as: VivoPartitionRuntimeStatusABI.self)
        return .init(
            flags: .init(rawValue: raw.flags),
            invalidElementCount: raw.invalidElementCount,
            firstInvalidElement: raw.firstInvalidElement == UInt32.max ? nil : raw.firstInvalidElement,
            maximumConcentration: Float(bitPattern: raw.maximumConcentrationBits),
            maximumNegativeMagnitude: Float(bitPattern: raw.maximumNegativeMagnitudeBits)
        )
    }

    private func readAuthoritativeState() async throws -> [Float] {
        guard let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoPartitionRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Partition.Readback"
        blit.copy(from: authoritativeState, sourceOffset: 0, to: readbackBuffer, destinationOffset: 0, size: byteCount)
        blit.endEncoding()
        try await complete(command)
        let pointer = readbackBuffer.contents().assumingMemoryBound(to: Float.self)
        return Array(UnsafeBufferPointer(start: pointer, count: elementCount))
    }

    private func rejectedCertificate(
        requested: Float,
        stability: Float,
        attempts: UInt32,
        status: VivoPartitionRuntimeStatus
    ) -> VivoPartitionStepCertificate {
        .init(
            disposition: .rejected,
            modelFingerprint: model.fingerprint,
            deviceName: deviceName,
            registryID: registryID,
            stepIndex: nextStepIndex,
            requestedStep: requested,
            stabilityBound: stability,
            acceptedStep: nil,
            attemptCount: attempts,
            timeBefore: absoluteTime,
            timeAfter: absoluteTime,
            status: status,
            totalAmountBefore: nil,
            totalAmountAfter: nil
        )
    }

    private func encode(
        pipeline: MTLComputePipelineState,
        commandBuffer: MTLCommandBuffer,
        count: Int,
        label: String,
        bindings: (MTLComputeCommandEncoder) -> Void
    ) throws {
        guard count > 0, let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoPartitionRuntimeError.commandQueueUnavailable
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline)
        bindings(encoder)
        let width = Self.threadgroupWidth(pipeline)
        encoder.dispatchThreads(
            MTLSize(width: count, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
        )
        encoder.endEncoding()
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(throwing: VivoPartitionRuntimeError.commandFailed(
                        completed.error?.localizedDescription ?? "status \(completed.status.rawValue)"
                    ))
                }
            }
            command.commit()
        }
    }

    private static func threadgroupWidth(_ pipeline: MTLComputePipelineState) -> Int {
        let width = max(1, pipeline.threadExecutionWidth)
        let groups = max(1, min(8, pipeline.maxTotalThreadsPerThreadgroup / width))
        return width * groups
    }

    private static func makePrivateBuffer<T>(
        device: MTLDevice,
        queue: MTLCommandQueue,
        values: [T],
        label: String
    ) throws -> MTLBuffer {
        guard !values.isEmpty else { throw VivoPartitionRuntimeError.allocationFailed("\(label) is empty") }
        return try makePrivateBufferAllowingEmpty(
            device: device,
            queue: queue,
            values: values,
            fallback: values[0],
            label: label
        )
    }

    private static func makePrivateBufferAllowingEmpty<T>(
        device: MTLDevice,
        queue: MTLCommandQueue,
        values: [T],
        fallback: T,
        label: String
    ) throws -> MTLBuffer {
        let uploaded = values.isEmpty ? [fallback] : values
        let size = uploaded.count.multipliedReportingOverflow(by: MemoryLayout<T>.stride)
        guard !size.overflow else { throw VivoPartitionRuntimeError.arithmeticOverflow }
        guard let staging = device.makeBuffer(length: size.partialValue, options: .storageModeShared),
              let destination = device.makeBuffer(length: size.partialValue, options: .storageModePrivate),
              let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoPartitionRuntimeError.allocationFailed(label)
        }
        uploaded.withUnsafeBytes { bytes in
            staging.contents().copyMemory(from: bytes.baseAddress!, byteCount: size.partialValue)
        }
        staging.label = "\(label).Staging"
        destination.label = label
        blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: size.partialValue)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoPartitionRuntimeError.commandFailed(command.error?.localizedDescription ?? "\(label) upload failed")
        }
        return destination
    }
}

private extension Optional {
    func orThrowPartition(_ error: @autoclosure () -> Error) throws -> Wrapped {
        guard let value = self else { throw error() }
        return value
    }
}
