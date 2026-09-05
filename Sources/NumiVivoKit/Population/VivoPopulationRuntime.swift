import Foundation
import Metal
import NumiVivoShaders

extension VivoPopulationPhenotypeABI: @unchecked Sendable {}
extension VivoPopulationTransitionABI: @unchecked Sendable {}

public struct VivoPopulationRuntimeConfiguration: Codable, Equatable, Sendable {
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
        negativeTolerance: Float = 1e-6
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
            throw VivoPopulationRuntimeError.invalidConfiguration
        }
    }
}

public struct VivoPopulationStatusFlags: OptionSet, Codable, Sendable {
    public let rawValue: UInt32
    public init(rawValue: UInt32) { self.rawValue = rawValue }

    public static let needsSmallerStep = Self(rawValue: 1 << 0)
    public static let nonFinite = Self(rawValue: 1 << 1)
    public static let negativeDensity = Self(rawValue: 1 << 2)
    public static let densityOverflow = Self(rawValue: 1 << 3)
    public static let invalidTopology = Self(rawValue: 1 << 4)

    public static let terminalInvalid: Self = [.nonFinite, .invalidTopology]
}

public struct VivoPopulationRuntimeStatus: Codable, Equatable, Sendable {
    public let flags: VivoPopulationStatusFlags
    public let invalidElementCount: UInt32
    public let firstInvalidElement: UInt32?
    public let maximumDensity: Float
    public let maximumNegativeMagnitude: Float
}

public enum VivoPopulationStepDisposition: String, Codable, Sendable {
    case committed
    case committedWithReducedStep
    case rejected
}

public struct VivoPopulationStepCertificate: Codable, Equatable, Sendable {
    public let disposition: VivoPopulationStepDisposition
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
    public let status: VivoPopulationRuntimeStatus
}

public struct VivoPopulationStateSnapshot: Codable, Equatable, Sendable {
    public let modelFingerprint: String
    public let stepIndex: UInt32
    public let absoluteTime: Float
    public let phenotypeCount: UInt32
    public let grid: VivoPopulationGrid
    public let densities: [Float]

    public func density(phenotype: UInt32, voxel: UInt32) -> Float? {
        guard let voxelCount = grid.voxelCount,
              phenotype < phenotypeCount,
              voxel < voxelCount else { return nil }
        return densities[Int(phenotype * voxelCount + voxel)]
    }
}

public enum VivoPopulationRuntimeError: Error, LocalizedError, Sendable {
    case invalidConfiguration
    case invalidInitialState(String)
    case invalidRegulatorState(String)
    case invalidVelocity(String)
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
        case .invalidConfiguration: return "invalid population runtime configuration"
        case .invalidInitialState(let reason): return "invalid initial population state: \(reason)"
        case .invalidRegulatorState(let reason): return "invalid regulator field state: \(reason)"
        case .invalidVelocity(let reason): return "invalid population velocity field: \(reason)"
        case .deviceUnavailable: return "no Metal device is available"
        case .commandQueueUnavailable: return "Metal command queue is unavailable"
        case .shaderResourceUnavailable: return "NumiVivo population Metal source is unavailable"
        case .shaderCompilation(let reason): return "population Metal compilation failed: \(reason)"
        case .pipelineUnavailable(let name): return "population pipeline is unavailable: \(name)"
        case .allocationFailed(let name): return "population Metal allocation failed: \(name)"
        case .commandFailed(let reason): return "population Metal command failed: \(reason)"
        case .arithmeticOverflow: return "population runtime arithmetic overflow"
        }
    }
}

private struct VivoPopulationRuntimeCommandABI {
    var phenotypeCount: UInt32
    var voxelCount: UInt32
    var regulatorFieldCount: UInt32
    var transitionCount: UInt32
    var width: UInt32
    var height: UInt32
    var depth: UInt32
    var boundaryMode: UInt32
    var spacingX: Float
    var spacingY: Float
    var spacingZ: Float
    var deltaTime: Float
    var negativeTolerance: Float
    var absoluteTime: Float
    var stepIndex: UInt32
    var attemptIndex: UInt32
}

private struct VivoPopulationRuntimeStatusABI {
    var flags: UInt32
    var invalidElementCount: UInt32
    var firstInvalidElement: UInt32
    var maximumDensityBits: UInt32
    var minimumDensityMagnitudeBits: UInt32
    var reserved0: UInt32
    var reserved1: UInt32
    var reserved2: UInt32
}

/// Continuous multicellular population runtime with logistic proliferation,
/// death, phenotype transitions, dense inter-population interactions, diffusion,
/// and upwind advection. RK2 candidates remain private until validation succeeds.
public actor VivoPopulationRuntime {
    public nonisolated let model: VivoPopulationModel
    public nonisolated let configuration: VivoPopulationRuntimeConfiguration
    public nonisolated let deviceName: String
    public nonisolated let registryID: UInt64

    private let device: MTLDevice
    private let queue: MTLCommandQueue
    private let compiled: VivoPopulationCompiledModel
    private let clearPipeline: MTLComputePipelineState
    private let stagePipeline: MTLComputePipelineState
    private let finalizePipeline: MTLComputePipelineState
    private let validatePipeline: MTLComputePipelineState
    private let phenotypeBuffer: MTLBuffer
    private let transitionBuffer: MTLBuffer
    private let interactionBuffer: MTLBuffer
    private let regulatorBuffer: MTLBuffer
    private let velocityBuffer: MTLBuffer
    private var stateA: MTLBuffer
    private var stateB: MTLBuffer
    private let stageState: MTLBuffer
    private let firstDerivative: MTLBuffer
    private let statusBuffer: MTLBuffer
    private let readbackBuffer: MTLBuffer
    private let stateElementCount: Int
    private let stateByteCount: Int
    private let regulatorElementCount: Int
    private let voxelCount: UInt32
    private var regulatorMirror: [Float]
    private var velocityMirror: [SIMD4<Float>]
    private var maximumVelocity = SIMD3<Float>(repeating: 0)
    private var authoritativeIsA = true
    private var absoluteTime: Float = 0
    private var nextStepIndex: UInt32 = 0

    public static func make(
        model: VivoPopulationModel,
        configuration: VivoPopulationRuntimeConfiguration,
        initialDensities: [Float],
        initialRegulatorFields: [Float]? = nil,
        initialVelocity: [SIMD3<Float>]? = nil,
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoPopulationRuntime {
        try model.validate()
        try configuration.validate()
        let compiled = try VivoPopulationCompiledModel(model: model)
        guard let voxelCount = model.grid.voxelCount else {
            throw VivoPopulationRuntimeError.arithmeticOverflow
        }
        let phenotypeCount = model.phenotypes.count
        let stateElements = phenotypeCount.multipliedReportingOverflow(by: Int(voxelCount))
        guard !stateElements.overflow else { throw VivoPopulationRuntimeError.arithmeticOverflow }
        let stateBytes = stateElements.partialValue.multipliedReportingOverflow(by: MemoryLayout<Float>.stride)
        guard !stateBytes.overflow, stateBytes.partialValue > 0 else {
            throw VivoPopulationRuntimeError.arithmeticOverflow
        }

        let expandedDensities: [Float]
        if initialDensities.count == phenotypeCount {
            var values = [Float](repeating: 0, count: stateElements.partialValue)
            for phenotype in 0..<phenotypeCount {
                guard initialDensities[phenotype].isFinite,
                      initialDensities[phenotype] >= model.phenotypes[phenotype].minimumDensity,
                      initialDensities[phenotype] <= model.phenotypes[phenotype].maximumDensity else {
                    throw VivoPopulationRuntimeError.invalidInitialState("phenotype \(phenotype) is outside its bounds")
                }
                for voxel in 0..<Int(voxelCount) {
                    values[phenotype * Int(voxelCount) + voxel] = initialDensities[phenotype]
                }
            }
            expandedDensities = values
        } else if initialDensities.count == stateElements.partialValue {
            for phenotype in 0..<phenotypeCount {
                let bounds = model.phenotypes[phenotype]
                let base = phenotype * Int(voxelCount)
                for voxel in 0..<Int(voxelCount) {
                    let value = initialDensities[base + voxel]
                    guard value.isFinite, value >= bounds.minimumDensity, value <= bounds.maximumDensity else {
                        throw VivoPopulationRuntimeError.invalidInitialState("density at phenotype \(phenotype), voxel \(voxel) is invalid")
                    }
                }
            }
            expandedDensities = initialDensities
        } else {
            throw VivoPopulationRuntimeError.invalidInitialState(
                "expected phenotypeCount or phenotypeCount × voxelCount values"
            )
        }

        let fieldElementsResult = model.regulatorFields.count.multipliedReportingOverflow(by: Int(voxelCount))
        guard !fieldElementsResult.overflow else { throw VivoPopulationRuntimeError.arithmeticOverflow }
        let fieldElements = fieldElementsResult.partialValue
        let regulatorValues = try validatedRegulators(
            model: model,
            values: initialRegulatorFields ?? [Float](repeating: 0, count: fieldElements),
            voxelCount: Int(voxelCount)
        )
        let velocityValues = try validatedVelocity(
            initialVelocity ?? [SIMD3<Float>](repeating: .zero, count: Int(voxelCount)),
            voxelCount: Int(voxelCount)
        )

        let device = try requestedDevice ?? MTLCreateSystemDefaultDevice().orThrowPopulation(VivoPopulationRuntimeError.deviceUnavailable)
        guard let queue = device.makeCommandQueue() else {
            throw VivoPopulationRuntimeError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.Population.Queue"
        guard let sourceURL = NumiVivoShaderResources.bundle.url(
            forResource: "NumiVivoPopulationKernels",
            withExtension: "metal"
        ) else {
            throw VivoPopulationRuntimeError.shaderResourceUnavailable
        }
        let source = String(decoding: try Data(contentsOf: sourceURL), as: UTF8.self)
        let library: MTLLibrary
        do { library = try await device.makeLibrary(source: source, options: nil) }
        catch { throw VivoPopulationRuntimeError.shaderCompilation(error.localizedDescription) }

        func pipeline(_ name: String) throws -> MTLComputePipelineState {
            guard let function = library.makeFunction(name: name) else {
                throw VivoPopulationRuntimeError.pipelineUnavailable(name)
            }
            do { return try device.makeComputePipelineState(function: function) }
            catch { throw VivoPopulationRuntimeError.shaderCompilation("\(name): \(error.localizedDescription)") }
        }

        let phenotypeBuffer = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.phenotypes,
            label: "NumiVivo.Population.Phenotypes"
        )
        let transitionBuffer = try makePrivateBufferAllowingEmpty(
            device: device,
            queue: queue,
            values: compiled.transitions,
            fallback: VivoPopulationTransitionABI(
                sourcePhenotype: 0,
                destinationPhenotype: 0,
                regulatorField: UInt32.max,
                mode: VivoPopulationTransitionMode.constitutive.rawValue,
                baseRate: 0,
                maximumRegulatedRate: 0,
                threshold: 1,
                hillCoefficient: 1
            ),
            label: "NumiVivo.Population.Transitions"
        )
        let interactionBuffer = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.interactionMatrix,
            label: "NumiVivo.Population.Interactions"
        )
        let regulatorBuffer = try makePrivateBufferAllowingEmpty(
            device: device,
            queue: queue,
            values: regulatorValues,
            fallback: Float.zero,
            label: "NumiVivo.Population.Regulators"
        )
        let velocityBuffer = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: velocityValues,
            label: "NumiVivo.Population.Velocity"
        )
        let stateA = try makePrivateBuffer(
            device: device,
            queue: queue,
            values: expandedDensities,
            label: "NumiVivo.Population.StateA"
        )
        guard let stateB = device.makeBuffer(length: stateBytes.partialValue, options: .storageModePrivate),
              let stageState = device.makeBuffer(length: stateBytes.partialValue, options: .storageModePrivate),
              let firstDerivative = device.makeBuffer(length: stateBytes.partialValue, options: .storageModePrivate),
              let status = device.makeBuffer(length: MemoryLayout<VivoPopulationRuntimeStatusABI>.stride, options: .storageModeShared),
              let readback = device.makeBuffer(length: stateBytes.partialValue, options: .storageModeShared) else {
            throw VivoPopulationRuntimeError.allocationFailed("transaction state")
        }
        stateB.label = "NumiVivo.Population.StateB"
        stageState.label = "NumiVivo.Population.Stage"
        firstDerivative.label = "NumiVivo.Population.FirstDerivative"
        status.label = "NumiVivo.Population.Status"
        readback.label = "NumiVivo.Population.Readback"

        return VivoPopulationRuntime(
            model: model,
            configuration: configuration,
            device: device,
            queue: queue,
            compiled: compiled,
            clearPipeline: try pipeline("numivivo_population::nvivo_population_clear_status"),
            stagePipeline: try pipeline("numivivo_population::nvivo_population_stage"),
            finalizePipeline: try pipeline("numivivo_population::nvivo_population_finalize"),
            validatePipeline: try pipeline("numivivo_population::nvivo_population_validate"),
            phenotypeBuffer: phenotypeBuffer,
            transitionBuffer: transitionBuffer,
            interactionBuffer: interactionBuffer,
            regulatorBuffer: regulatorBuffer,
            velocityBuffer: velocityBuffer,
            stateA: stateA,
            stateB: stateB,
            stageState: stageState,
            firstDerivative: firstDerivative,
            statusBuffer: status,
            readbackBuffer: readback,
            stateElementCount: stateElements.partialValue,
            stateByteCount: stateBytes.partialValue,
            regulatorElementCount: fieldElements,
            voxelCount: voxelCount,
            regulatorMirror: regulatorValues,
            velocityMirror: velocityValues
        )
    }

    private init(
        model: VivoPopulationModel,
        configuration: VivoPopulationRuntimeConfiguration,
        device: MTLDevice,
        queue: MTLCommandQueue,
        compiled: VivoPopulationCompiledModel,
        clearPipeline: MTLComputePipelineState,
        stagePipeline: MTLComputePipelineState,
        finalizePipeline: MTLComputePipelineState,
        validatePipeline: MTLComputePipelineState,
        phenotypeBuffer: MTLBuffer,
        transitionBuffer: MTLBuffer,
        interactionBuffer: MTLBuffer,
        regulatorBuffer: MTLBuffer,
        velocityBuffer: MTLBuffer,
        stateA: MTLBuffer,
        stateB: MTLBuffer,
        stageState: MTLBuffer,
        firstDerivative: MTLBuffer,
        statusBuffer: MTLBuffer,
        readbackBuffer: MTLBuffer,
        stateElementCount: Int,
        stateByteCount: Int,
        regulatorElementCount: Int,
        voxelCount: UInt32,
        regulatorMirror: [Float],
        velocityMirror: [SIMD4<Float>]
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
        self.phenotypeBuffer = phenotypeBuffer
        self.transitionBuffer = transitionBuffer
        self.interactionBuffer = interactionBuffer
        self.regulatorBuffer = regulatorBuffer
        self.velocityBuffer = velocityBuffer
        self.stateA = stateA
        self.stateB = stateB
        self.stageState = stageState
        self.firstDerivative = firstDerivative
        self.statusBuffer = statusBuffer
        self.readbackBuffer = readbackBuffer
        self.stateElementCount = stateElementCount
        self.stateByteCount = stateByteCount
        self.regulatorElementCount = regulatorElementCount
        self.voxelCount = voxelCount
        self.regulatorMirror = regulatorMirror
        self.velocityMirror = velocityMirror
        self.deviceName = device.name
        self.registryID = device.registryID
        self.maximumVelocity = Self.maximumVelocity(velocityMirror)
    }

    public func time() -> Float { absoluteTime }
    public func stepIndex() -> UInt32 { nextStepIndex }

    public func regulatorFields() -> [Float] { regulatorMirror }

    public func setRegulatorFields(_ values: [Float]) throws {
        regulatorMirror = try Self.validatedRegulators(
            model: model,
            values: values,
            voxelCount: Int(voxelCount)
        )
        if regulatorElementCount > 0 {
            try replacePrivateBuffer(regulatorMirror, destination: regulatorBuffer, label: "RegulatorUpdate")
        }
    }

    public func setRegulatorField(id: String, values: [Float]) throws {
        guard let field = model.regulatorFields.firstIndex(where: { $0.id == id }) else {
            throw VivoPopulationRuntimeError.invalidRegulatorState("unknown field \(id)")
        }
        guard values.count == Int(voxelCount) else {
            throw VivoPopulationRuntimeError.invalidRegulatorState("field \(id) has incorrect voxel count")
        }
        let definition = model.regulatorFields[field]
        guard values.allSatisfy({
            $0.isFinite && $0 >= definition.minimum && $0 <= definition.maximum
        }) else {
            throw VivoPopulationRuntimeError.invalidRegulatorState("field \(id) violates bounds")
        }
        let start = field * Int(voxelCount)
        regulatorMirror.replaceSubrange(start..<(start + Int(voxelCount)), with: values)
        try replacePrivateBuffer(regulatorMirror, destination: regulatorBuffer, label: "RegulatorFieldUpdate")
    }

    public func setVelocity(_ values: [SIMD3<Float>]) throws {
        velocityMirror = try Self.validatedVelocity(values, voxelCount: Int(voxelCount))
        maximumVelocity = Self.maximumVelocity(velocityMirror)
        try replacePrivateBuffer(velocityMirror, destination: velocityBuffer, label: "VelocityUpdate")
    }

    public func conservativeMaximumStep() -> Float {
        let grid = model.grid
        let advection = maximumVelocity.x / grid.spacingXMetres +
            maximumVelocity.y / grid.spacingYMetres +
            maximumVelocity.z / grid.spacingZMetres
        let total = compiled.conservativeKineticRate + compiled.conservativeDiffusionRate + advection
        guard total.isFinite, total > 0 else { return configuration.maximumTimeStep }
        return max(configuration.minimumTimeStep, min(configuration.maximumTimeStep, 0.35 / total))
    }

    public func step(
        deltaTime requestedStep: Float? = nil,
        permitAdaptiveReduction: Bool = true
    ) async throws -> VivoPopulationStepCertificate {
        let requested = requestedStep ?? configuration.timeStep
        guard requested.isFinite,
              requested >= configuration.minimumTimeStep,
              requested <= configuration.maximumTimeStep else {
            throw VivoPopulationRuntimeError.invalidConfiguration
        }
        let stability = conservativeMaximumStep()
        var candidate = requested
        if candidate > stability {
            guard permitAdaptiveReduction else {
                return rejectedCertificate(
                    requested: requested,
                    stability: stability,
                    attempts: 0,
                    status: .init(
                        flags: [.needsSmallerStep],
                        invalidElementCount: 0,
                        firstInvalidElement: nil,
                        maximumDensity: 0,
                        maximumNegativeMagnitude: 0
                    )
                )
            }
            candidate = stability
        }

        var last = VivoPopulationRuntimeStatus(
            flags: [],
            invalidElementCount: 0,
            firstInvalidElement: nil,
            maximumDensity: 0,
            maximumNegativeMagnitude: 0
        )
        for attempt in 0..<configuration.maximumAttempts {
            last = try await executeAttempt(deltaTime: candidate, attempt: attempt)
            if !last.flags.intersection(.terminalInvalid).isEmpty {
                return rejectedCertificate(
                    requested: requested,
                    stability: stability,
                    attempts: attempt + 1,
                    status: last
                )
            }
            if last.flags.contains(.needsSmallerStep) {
                let reduced = candidate * 0.5
                guard permitAdaptiveReduction,
                      attempt + 1 < configuration.maximumAttempts,
                      reduced >= configuration.minimumTimeStep,
                      reduced < candidate else {
                    return rejectedCertificate(
                        requested: requested,
                        stability: stability,
                        attempts: attempt + 1,
                        status: last
                    )
                }
                candidate = reduced
                continue
            }

            authoritativeIsA.toggle()
            let before = absoluteTime
            absoluteTime += candidate
            let completedStep = nextStepIndex
            nextStepIndex &+= 1
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
                status: last
            )
        }
        return rejectedCertificate(
            requested: requested,
            stability: stability,
            attempts: configuration.maximumAttempts,
            status: last
        )
    }

    public func snapshot() async throws -> VivoPopulationStateSnapshot {
        guard let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoPopulationRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Population.Snapshot"
        blit.copy(
            from: authoritativeState,
            sourceOffset: 0,
            to: readbackBuffer,
            destinationOffset: 0,
            size: stateByteCount
        )
        blit.endEncoding()
        try await complete(command)
        let pointer = readbackBuffer.contents().assumingMemoryBound(to: Float.self)
        return .init(
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            absoluteTime: absoluteTime,
            phenotypeCount: UInt32(model.phenotypes.count),
            grid: model.grid,
            densities: Array(UnsafeBufferPointer(start: pointer, count: stateElementCount))
        )
    }

    private var authoritativeState: MTLBuffer { authoritativeIsA ? stateA : stateB }
    private var candidateState: MTLBuffer { authoritativeIsA ? stateB : stateA }

    private func executeAttempt(deltaTime: Float, attempt: UInt32) async throws -> VivoPopulationRuntimeStatus {
        guard let command = queue.makeCommandBuffer() else {
            throw VivoPopulationRuntimeError.commandQueueUnavailable
        }
        command.label = "NumiVivo.Population.Step.\(nextStepIndex).Attempt.\(attempt)"
        var runtimeCommand = VivoPopulationRuntimeCommandABI(
            phenotypeCount: UInt32(model.phenotypes.count),
            voxelCount: voxelCount,
            regulatorFieldCount: UInt32(model.regulatorFields.count),
            transitionCount: UInt32(model.transitions.count),
            width: model.grid.width,
            height: model.grid.height,
            depth: model.grid.depth,
            boundaryMode: model.grid.boundary.rawValue,
            spacingX: model.grid.spacingXMetres,
            spacingY: model.grid.spacingYMetres,
            spacingZ: model.grid.spacingZMetres,
            deltaTime: deltaTime,
            negativeTolerance: configuration.negativeTolerance,
            absoluteTime: absoluteTime + deltaTime,
            stepIndex: nextStepIndex,
            attemptIndex: attempt
        )

        try encode(pipeline: clearPipeline, commandBuffer: command, count: 1, label: "Population.Clear") { encoder in
            encoder.setBuffer(statusBuffer, offset: 0, index: 0)
        }
        try encode(pipeline: stagePipeline, commandBuffer: command, count: stateElementCount, label: "Population.Stage") { encoder in
            encoder.setBuffer(authoritativeState, offset: 0, index: 0)
            encoder.setBuffer(stageState, offset: 0, index: 1)
            encoder.setBuffer(firstDerivative, offset: 0, index: 2)
            encoder.setBuffer(phenotypeBuffer, offset: 0, index: 3)
            encoder.setBuffer(transitionBuffer, offset: 0, index: 4)
            encoder.setBuffer(interactionBuffer, offset: 0, index: 5)
            encoder.setBuffer(regulatorBuffer, offset: 0, index: 6)
            encoder.setBuffer(velocityBuffer, offset: 0, index: 7)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPopulationRuntimeCommandABI>.stride, index: 8)
        }
        try encode(pipeline: finalizePipeline, commandBuffer: command, count: stateElementCount, label: "Population.Finalize") { encoder in
            encoder.setBuffer(authoritativeState, offset: 0, index: 0)
            encoder.setBuffer(stageState, offset: 0, index: 1)
            encoder.setBuffer(firstDerivative, offset: 0, index: 2)
            encoder.setBuffer(candidateState, offset: 0, index: 3)
            encoder.setBuffer(phenotypeBuffer, offset: 0, index: 4)
            encoder.setBuffer(transitionBuffer, offset: 0, index: 5)
            encoder.setBuffer(interactionBuffer, offset: 0, index: 6)
            encoder.setBuffer(regulatorBuffer, offset: 0, index: 7)
            encoder.setBuffer(velocityBuffer, offset: 0, index: 8)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPopulationRuntimeCommandABI>.stride, index: 9)
        }
        try encode(pipeline: validatePipeline, commandBuffer: command, count: stateElementCount, label: "Population.Validate") { encoder in
            encoder.setBuffer(candidateState, offset: 0, index: 0)
            encoder.setBuffer(phenotypeBuffer, offset: 0, index: 1)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoPopulationRuntimeCommandABI>.stride, index: 2)
            encoder.setBuffer(statusBuffer, offset: 0, index: 3)
        }

        try await complete(command)
        let raw = statusBuffer.contents().load(as: VivoPopulationRuntimeStatusABI.self)
        return .init(
            flags: .init(rawValue: raw.flags),
            invalidElementCount: raw.invalidElementCount,
            firstInvalidElement: raw.firstInvalidElement == UInt32.max ? nil : raw.firstInvalidElement,
            maximumDensity: Float(bitPattern: raw.maximumDensityBits),
            maximumNegativeMagnitude: Float(bitPattern: raw.minimumDensityMagnitudeBits)
        )
    }

    private func rejectedCertificate(
        requested: Float,
        stability: Float,
        attempts: UInt32,
        status: VivoPopulationRuntimeStatus
    ) -> VivoPopulationStepCertificate {
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
            status: status
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
            throw VivoPopulationRuntimeError.commandQueueUnavailable
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

    private func replacePrivateBuffer<T>(_ values: [T], destination: MTLBuffer, label: String) throws {
        let byteCount = values.count.multipliedReportingOverflow(by: MemoryLayout<T>.stride)
        guard !byteCount.overflow, byteCount.partialValue <= destination.length else {
            throw VivoPopulationRuntimeError.arithmeticOverflow
        }
        guard let staging = device.makeBuffer(length: max(1, byteCount.partialValue), options: .storageModeShared),
              let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoPopulationRuntimeError.allocationFailed(label)
        }
        if byteCount.partialValue > 0 {
            values.withUnsafeBytes { bytes in
                staging.contents().copyMemory(from: bytes.baseAddress!, byteCount: byteCount.partialValue)
            }
            blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: byteCount.partialValue)
        }
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoPopulationRuntimeError.commandFailed(command.error?.localizedDescription ?? "\(label) upload failed")
        }
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(throwing: VivoPopulationRuntimeError.commandFailed(
                        completed.error?.localizedDescription ?? "status \(completed.status.rawValue)"
                    ))
                }
            }
            command.commit()
        }
    }

    private static func validatedRegulators(
        model: VivoPopulationModel,
        values: [Float],
        voxelCount: Int
    ) throws -> [Float] {
        let expected = model.regulatorFields.count.multipliedReportingOverflow(by: voxelCount)
        guard !expected.overflow, values.count == expected.partialValue else {
            throw VivoPopulationRuntimeError.invalidRegulatorState("field array has incorrect length")
        }
        for (fieldIndex, field) in model.regulatorFields.enumerated() {
            let base = fieldIndex * voxelCount
            for voxel in 0..<voxelCount {
                let value = values[base + voxel]
                guard value.isFinite, value >= field.minimum, value <= field.maximum else {
                    throw VivoPopulationRuntimeError.invalidRegulatorState("field \(field.id), voxel \(voxel) violates bounds")
                }
            }
        }
        return values
    }

    private static func validatedVelocity(_ values: [SIMD3<Float>], voxelCount: Int) throws -> [SIMD4<Float>] {
        guard values.count == voxelCount else {
            throw VivoPopulationRuntimeError.invalidVelocity("velocity array has incorrect length")
        }
        return try values.map { value in
            guard value.x.isFinite, value.y.isFinite, value.z.isFinite else {
                throw VivoPopulationRuntimeError.invalidVelocity("velocity contains a non-finite component")
            }
            return SIMD4(value.x, value.y, value.z, 0)
        }
    }

    private static func maximumVelocity(_ values: [SIMD4<Float>]) -> SIMD3<Float> {
        values.reduce(into: SIMD3<Float>(repeating: 0)) { result, value in
            result.x = max(result.x, abs(value.x))
            result.y = max(result.y, abs(value.y))
            result.z = max(result.z, abs(value.z))
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
        guard !values.isEmpty else { throw VivoPopulationRuntimeError.allocationFailed("\(label) is empty") }
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
        let byteCount = uploaded.count.multipliedReportingOverflow(by: MemoryLayout<T>.stride)
        guard !byteCount.overflow else { throw VivoPopulationRuntimeError.arithmeticOverflow }
        guard let staging = device.makeBuffer(length: byteCount.partialValue, options: .storageModeShared),
              let destination = device.makeBuffer(length: byteCount.partialValue, options: .storageModePrivate),
              let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoPopulationRuntimeError.allocationFailed(label)
        }
        uploaded.withUnsafeBytes { bytes in
            staging.contents().copyMemory(from: bytes.baseAddress!, byteCount: byteCount.partialValue)
        }
        staging.label = "\(label).Staging"
        destination.label = label
        blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: byteCount.partialValue)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoPopulationRuntimeError.commandFailed(command.error?.localizedDescription ?? "\(label) upload failed")
        }
        return destination
    }
}

private extension Optional {
    func orThrowPopulation(_ error: @autoclosure () -> Error) throws -> Wrapped {
        guard let value = self else { throw error() }
        return value
    }
}
