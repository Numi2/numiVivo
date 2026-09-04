import CryptoKit
import Foundation
import Metal
import NumiVivoShaders

public enum VivoSSAPropensityLaw: UInt32, Codable, CaseIterable, Sendable {
    case zeroOrder = 0
    case firstOrder = 1
    case secondOrderDistinct = 2
    case secondOrderSame = 3
    case hillActivation = 4
    case hillRepression = 5
}

public struct VivoSSAStoichiometricChange: Codable, Hashable, Sendable {
    public let speciesIndex: UInt32
    public let delta: Int32

    public init(speciesIndex: UInt32, delta: Int32) {
        self.speciesIndex = speciesIndex
        self.delta = delta
    }
}

public struct VivoSSAReaction: Codable, Hashable, Sendable {
    public let id: String
    public let law: VivoSSAPropensityLaw
    public let reactantA: UInt32?
    public let reactantB: UInt32?
    public let rate: Float
    public let parameter1: Float
    public let parameter2: Float
    public let parameter3: Float
    public let changes: [VivoSSAStoichiometricChange]
    public let critical: Bool

    public init(
        id: String,
        law: VivoSSAPropensityLaw,
        reactantA: UInt32? = nil,
        reactantB: UInt32? = nil,
        rate: Float,
        parameter1: Float = 0,
        parameter2: Float = 0,
        parameter3: Float = 0,
        changes: [VivoSSAStoichiometricChange],
        critical: Bool = false
    ) {
        self.id = id
        self.law = law
        self.reactantA = reactantA
        self.reactantB = reactantB
        self.rate = rate
        self.parameter1 = parameter1
        self.parameter2 = parameter2
        self.parameter3 = parameter3
        self.changes = changes
        self.critical = critical
    }
}

public struct VivoExactSSAModel: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let species: [String]
    public let reactions: [VivoSSAReaction]
    public let fingerprint: String

    public init(species: [String], reactions: [VivoSSAReaction]) throws {
        let unsigned = VivoExactSSAModel(
            schemaVersion: 1,
            species: species,
            reactions: reactions,
            fingerprint: ""
        )
        try unsigned.validate()
        self = .init(
            schemaVersion: unsigned.schemaVersion,
            species: unsigned.species,
            reactions: unsigned.reactions,
            fingerprint: try Self.fingerprint(unsigned)
        )
    }

    private init(
        schemaVersion: UInt32,
        species: [String],
        reactions: [VivoSSAReaction],
        fingerprint: String
    ) {
        self.schemaVersion = schemaVersion
        self.species = species
        self.reactions = reactions
        self.fingerprint = fingerprint
    }

    public func validate() throws {
        guard schemaVersion == 1 else { throw VivoExactSSAError.unsupportedModelVersion(schemaVersion) }
        guard !species.isEmpty, species.count <= Int(UInt32.max) else {
            throw VivoExactSSAError.invalidModel("species table is empty or exceeds the ABI")
        }
        guard !reactions.isEmpty, reactions.count <= Int(UInt32.max) else {
            throw VivoExactSSAError.invalidModel("reaction table is empty or exceeds the ABI")
        }
        guard Set(species).count == species.count,
              species.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 256 }) else {
            throw VivoExactSSAError.invalidModel("species identifiers must be unique, non-empty and bounded")
        }
        guard Set(reactions.map(\.id)).count == reactions.count else {
            throw VivoExactSSAError.invalidModel("reaction identifiers must be unique")
        }

        for (index, reaction) in reactions.enumerated() {
            guard !reaction.id.isEmpty, reaction.id.utf8.count <= 256 else {
                throw VivoExactSSAError.invalidReaction(index, "invalid identifier")
            }
            guard reaction.rate.isFinite, reaction.rate >= 0,
                  reaction.parameter1.isFinite,
                  reaction.parameter2.isFinite,
                  reaction.parameter3.isFinite else {
                throw VivoExactSSAError.invalidReaction(index, "kinetic parameters must be finite and rate must be non-negative")
            }
            let requiredReactants: Int
            switch reaction.law {
            case .zeroOrder: requiredReactants = 0
            case .firstOrder, .secondOrderSame, .hillActivation, .hillRepression: requiredReactants = 1
            case .secondOrderDistinct: requiredReactants = 2
            }
            if requiredReactants >= 1, reaction.reactantA == nil {
                throw VivoExactSSAError.invalidReaction(index, "propensity law requires reactantA")
            }
            if requiredReactants == 2, reaction.reactantB == nil {
                throw VivoExactSSAError.invalidReaction(index, "propensity law requires reactantB")
            }
            for reactant in [reaction.reactantA, reaction.reactantB].compactMap({ $0 }) {
                guard reactant < UInt32(species.count) else {
                    throw VivoExactSSAError.invalidReaction(index, "reactant index is outside the species table")
                }
            }
            if reaction.law == .secondOrderDistinct,
               reaction.reactantA == reaction.reactantB {
                throw VivoExactSSAError.invalidReaction(index, "secondOrderDistinct requires different reactants")
            }
            if reaction.law == .secondOrderSame,
               reaction.reactantB != nil,
               reaction.reactantB != reaction.reactantA {
                throw VivoExactSSAError.invalidReaction(index, "secondOrderSame accepts only one reactant identity")
            }
            if reaction.law == .hillActivation || reaction.law == .hillRepression {
                guard reaction.parameter1 > 0, reaction.parameter2 >= 1 else {
                    throw VivoExactSSAError.invalidReaction(index, "Hill laws require threshold > 0 and coefficient >= 1")
                }
            }
            guard !reaction.changes.isEmpty else {
                throw VivoExactSSAError.invalidReaction(index, "reaction has no state change")
            }
            var seen = Set<UInt32>()
            for change in reaction.changes {
                guard change.speciesIndex < UInt32(species.count), change.delta != 0 else {
                    throw VivoExactSSAError.invalidReaction(index, "stoichiometric change has invalid species or zero delta")
                }
                guard seen.insert(change.speciesIndex).inserted else {
                    throw VivoExactSSAError.invalidReaction(index, "duplicate species changes must be combined before compilation")
                }
            }
        }
    }

    private static func fingerprint(_ model: VivoExactSSAModel) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return SHA256.hash(data: try encoder.encode(model)).map { String(format: "%02x", $0) }.joined()
    }
}

public struct VivoExactSSAConfiguration: Codable, Equatable, Sendable {
    public var laneCount: UInt32
    public var timeStep: Float
    public var minimumTimeStep: Float
    public var maximumTimeStep: Float
    public var maximumEventsPerLane: UInt32
    public var maximumAttempts: UInt32
    public var seed: UInt64

    public init(
        laneCount: UInt32,
        timeStep: Float,
        minimumTimeStep: Float = 1e-7,
        maximumTimeStep: Float = 60,
        maximumEventsPerLane: UInt32 = 4_096,
        maximumAttempts: UInt32 = 16,
        seed: UInt64 = 0x4e554d495649564f
    ) {
        self.laneCount = laneCount
        self.timeStep = timeStep
        self.minimumTimeStep = minimumTimeStep
        self.maximumTimeStep = maximumTimeStep
        self.maximumEventsPerLane = maximumEventsPerLane
        self.maximumAttempts = maximumAttempts
        self.seed = seed
    }

    public func validate() throws {
        guard laneCount > 0,
              timeStep.isFinite,
              minimumTimeStep.isFinite,
              maximumTimeStep.isFinite,
              minimumTimeStep > 0,
              maximumTimeStep >= minimumTimeStep,
              timeStep >= minimumTimeStep,
              timeStep <= maximumTimeStep,
              maximumEventsPerLane > 0,
              maximumAttempts > 0 else {
            throw VivoExactSSAError.invalidConfiguration
        }
    }
}

public struct VivoExactSSAStatusFlags: OptionSet, Codable, Sendable {
    public let rawValue: UInt32

    public init(rawValue: UInt32) { self.rawValue = rawValue }

    public static let needsSmallerStep = Self(rawValue: 1 << 0)
    public static let invalidPropensity = Self(rawValue: 1 << 1)
    public static let negativeCount = Self(rawValue: 1 << 2)
    public static let countOverflow = Self(rawValue: 1 << 3)
    public static let invalidModel = Self(rawValue: 1 << 4)

    public static let invalid: Self = [.invalidPropensity, .negativeCount, .countOverflow, .invalidModel]
}

public struct VivoExactSSAStatus: Codable, Equatable, Sendable {
    public let flags: VivoExactSSAStatusFlags
    public let lanesNeedingSmallerStep: UInt32
    public let invalidLaneCount: UInt32
    public let totalEventCount: UInt32
    public let maximumEventsObserved: UInt32
    public let firstInvalidLane: UInt32?
    public let firstInvalidReaction: UInt32?
}

public enum VivoExactSSAStepDisposition: String, Codable, Sendable {
    case committed
    case committedWithReducedStep
    case rejected
}

public struct VivoExactSSAStepCertificate: Codable, Equatable, Sendable {
    public let disposition: VivoExactSSAStepDisposition
    public let modelFingerprint: String
    public let deviceName: String
    public let registryID: UInt64
    public let stepIndex: UInt32
    public let requestedStep: Float
    public let acceptedStep: Float?
    public let attemptCount: UInt32
    public let timeBefore: Float
    public let timeAfter: Float
    public let status: VivoExactSSAStatus
}

public struct VivoExactSSAStateSnapshot: Codable, Equatable, Sendable {
    public let modelFingerprint: String
    public let stepIndex: UInt32
    public let absoluteTime: Float
    public let speciesCount: UInt32
    public let laneCount: UInt32
    public let counts: [UInt32]

    public func count(species: UInt32, lane: UInt32) -> UInt32? {
        guard species < speciesCount, lane < laneCount else { return nil }
        return counts[Int(species * laneCount + lane)]
    }
}

public enum VivoExactSSAError: Error, LocalizedError, Sendable {
    case unsupportedModelVersion(UInt32)
    case invalidModel(String)
    case invalidReaction(Int, String)
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
        case .unsupportedModelVersion(let version): return "unsupported exact-SSA model version \(version)"
        case .invalidModel(let reason): return "invalid exact-SSA model: \(reason)"
        case .invalidReaction(let index, let reason): return "invalid exact-SSA reaction \(index): \(reason)"
        case .invalidConfiguration: return "invalid exact-SSA runtime configuration"
        case .invalidInitialState(let reason): return "invalid exact-SSA initial state: \(reason)"
        case .deviceUnavailable: return "no Metal device is available"
        case .commandQueueUnavailable: return "Metal command queue is unavailable"
        case .shaderResourceUnavailable: return "NumiVivo exact-SSA Metal source is unavailable"
        case .shaderCompilation(let reason): return "exact-SSA Metal compilation failed: \(reason)"
        case .pipelineUnavailable(let name): return "exact-SSA pipeline is unavailable: \(name)"
        case .allocationFailed(let name): return "Metal allocation failed: \(name)"
        case .commandFailed(let reason): return "Metal command failed: \(reason)"
        case .arithmeticOverflow: return "exact-SSA allocation arithmetic overflow"
        }
    }
}

private struct VivoSSAReactionABI {
    var law: UInt32
    var reactantA: UInt32
    var reactantB: UInt32
    var stoichiometryOffset: UInt32
    var stoichiometryCount: UInt32
    var flags: UInt32
    var reserved0: UInt32 = 0
    var reserved1: UInt32 = 0
    var rate: Float
    var parameter1: Float
    var parameter2: Float
    var parameter3: Float
}

private struct VivoSSAStoichiometryABI {
    var speciesIndex: UInt32
    var delta: Int32
}

private struct VivoSSARuntimeCommandABI {
    var laneCount: UInt32
    var speciesCount: UInt32
    var reactionCount: UInt32
    var maximumEventsPerLane: UInt32
    var stepIndex: UInt32
    var attemptIndex: UInt32
    var reserved0: UInt32 = 0
    var reserved1: UInt32 = 0
    var deltaTime: Float
    var absoluteTime: Float
    var seed: UInt64
}

private struct VivoSSARuntimeStatusABI {
    var flags: UInt32
    var lanesNeedingSmallerStep: UInt32
    var invalidLaneCount: UInt32
    var totalEventCount: UInt32
    var maximumEventsObserved: UInt32
    var firstInvalidLane: UInt32
    var firstInvalidReaction: UInt32
    var reserved: UInt32
}

private struct VivoExactSSACompiledModel {
    let reactions: [VivoSSAReactionABI]
    let stoichiometry: [VivoSSAStoichiometryABI]

    init(model: VivoExactSSAModel) throws {
        var reactionRecords: [VivoSSAReactionABI] = []
        var stoichiometryRecords: [VivoSSAStoichiometryABI] = []
        reactionRecords.reserveCapacity(model.reactions.count)

        for reaction in model.reactions {
            guard stoichiometryRecords.count <= Int(UInt32.max) - reaction.changes.count else {
                throw VivoExactSSAError.arithmeticOverflow
            }
            let offset = UInt32(stoichiometryRecords.count)
            stoichiometryRecords.append(contentsOf: reaction.changes.map {
                VivoSSAStoichiometryABI(speciesIndex: $0.speciesIndex, delta: $0.delta)
            })
            reactionRecords.append(.init(
                law: reaction.law.rawValue,
                reactantA: reaction.reactantA ?? UInt32.max,
                reactantB: reaction.reactantB ?? UInt32.max,
                stoichiometryOffset: offset,
                stoichiometryCount: UInt32(reaction.changes.count),
                flags: reaction.critical ? 1 : 0,
                rate: reaction.rate,
                parameter1: reaction.parameter1,
                parameter2: reaction.parameter2,
                parameter3: reaction.parameter3
            ))
        }
        reactions = reactionRecords
        stoichiometry = stoichiometryRecords
    }
}

/// Exact Gillespie direct-method runtime. One GPU thread owns one independent
/// stochastic lane, so molecular updates require no atomics and the random stream
/// is independent of thread scheduling. Authoritative counts remain in private
/// Metal memory and are swapped only after a complete candidate passes validation.
public actor VivoExactSSARuntime {
    public nonisolated let model: VivoExactSSAModel
    public nonisolated let configuration: VivoExactSSAConfiguration
    public nonisolated let deviceName: String
    public nonisolated let registryID: UInt64

    private let device: MTLDevice
    private let queue: MTLCommandQueue
    private let advancePipeline: MTLComputePipelineState
    private let clearPipeline: MTLComputePipelineState
    private let validatePipeline: MTLComputePipelineState
    private let reactionBuffer: MTLBuffer
    private let stoichiometryBuffer: MTLBuffer
    private var stateA: MTLBuffer
    private var stateB: MTLBuffer
    private let statusBuffer: MTLBuffer
    private let eventCountBuffer: MTLBuffer
    private let readbackBuffer: MTLBuffer
    private let stateByteCount: Int
    private var authoritativeIsA = true
    private var nextStepIndex: UInt32 = 0
    private var absoluteTime: Float = 0

    public static func make(
        model: VivoExactSSAModel,
        configuration: VivoExactSSAConfiguration,
        initialCounts: [UInt32],
        device requestedDevice: MTLDevice? = nil
    ) async throws -> VivoExactSSARuntime {
        try model.validate()
        try configuration.validate()
        let device = try requestedDevice ?? MTLCreateSystemDefaultDevice().orThrow(VivoExactSSAError.deviceUnavailable)
        guard let queue = device.makeCommandQueue() else {
            throw VivoExactSSAError.commandQueueUnavailable
        }
        queue.label = "NumiVivo.ExactSSA.Queue"

        let sourceURL = NumiVivoShaderResources.bundle.url(
            forResource: "NumiVivoExactSSAKernels",
            withExtension: "metal"
        )
        guard let sourceURL else { throw VivoExactSSAError.shaderResourceUnavailable }
        let source = String(decoding: try Data(contentsOf: sourceURL), as: UTF8.self)
        let library: MTLLibrary
        do {
            library = try device.makeLibrary(source: source, options: nil)
        } catch {
            throw VivoExactSSAError.shaderCompilation(error.localizedDescription)
        }

        func pipeline(_ name: String) throws -> MTLComputePipelineState {
            guard let function = library.makeFunction(name: name) else {
                throw VivoExactSSAError.pipelineUnavailable(name)
            }
            do { return try device.makeComputePipelineState(function: function) }
            catch { throw VivoExactSSAError.shaderCompilation("\(name): \(error.localizedDescription)") }
        }

        let compiled = try VivoExactSSACompiledModel(model: model)
        let speciesCount = model.species.count
        let laneCount = Int(configuration.laneCount)
        let elements = speciesCount.multipliedReportingOverflow(by: laneCount)
        guard !elements.overflow else { throw VivoExactSSAError.arithmeticOverflow }
        let bytes = elements.partialValue.multipliedReportingOverflow(by: MemoryLayout<UInt32>.stride)
        guard !bytes.overflow, bytes.partialValue > 0 else { throw VivoExactSSAError.arithmeticOverflow }

        let expandedCounts: [UInt32]
        if initialCounts.count == speciesCount {
            var values = [UInt32](repeating: 0, count: elements.partialValue)
            for species in 0..<speciesCount {
                for lane in 0..<laneCount {
                    values[species * laneCount + lane] = initialCounts[species]
                }
            }
            expandedCounts = values
        } else if initialCounts.count == elements.partialValue {
            expandedCounts = initialCounts
        } else {
            throw VivoExactSSAError.invalidInitialState(
                "expected speciesCount or speciesCount × laneCount values"
            )
        }

        let reactionBuffer = try Self.makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.reactions,
            label: "NumiVivo.ExactSSA.Reactions"
        )
        let stoichiometryBuffer = try Self.makePrivateBuffer(
            device: device,
            queue: queue,
            values: compiled.stoichiometry,
            label: "NumiVivo.ExactSSA.Stoichiometry"
        )
        let stateA = try Self.makePrivateBuffer(
            device: device,
            queue: queue,
            values: expandedCounts,
            label: "NumiVivo.ExactSSA.StateA"
        )
        guard let stateB = device.makeBuffer(length: bytes.partialValue, options: .storageModePrivate),
              let status = device.makeBuffer(length: MemoryLayout<VivoSSARuntimeStatusABI>.stride, options: .storageModeShared),
              let events = device.makeBuffer(length: laneCount * MemoryLayout<UInt32>.stride, options: .storageModeShared),
              let readback = device.makeBuffer(length: bytes.partialValue, options: .storageModeShared) else {
            throw VivoExactSSAError.allocationFailed("transaction or diagnostics buffer")
        }
        stateB.label = "NumiVivo.ExactSSA.StateB"
        status.label = "NumiVivo.ExactSSA.Status"
        events.label = "NumiVivo.ExactSSA.EventCounts"
        readback.label = "NumiVivo.ExactSSA.Readback"

        return VivoExactSSARuntime(
            model: model,
            configuration: configuration,
            device: device,
            queue: queue,
            advancePipeline: try pipeline("nvivo_exact_ssa::nvivo_exact_ssa_advance"),
            clearPipeline: try pipeline("nvivo_exact_ssa::nvivo_exact_ssa_clear_status"),
            validatePipeline: try pipeline("nvivo_exact_ssa::nvivo_exact_ssa_validate_counts"),
            reactionBuffer: reactionBuffer,
            stoichiometryBuffer: stoichiometryBuffer,
            stateA: stateA,
            stateB: stateB,
            statusBuffer: status,
            eventCountBuffer: events,
            readbackBuffer: readback,
            stateByteCount: bytes.partialValue
        )
    }

    private init(
        model: VivoExactSSAModel,
        configuration: VivoExactSSAConfiguration,
        device: MTLDevice,
        queue: MTLCommandQueue,
        advancePipeline: MTLComputePipelineState,
        clearPipeline: MTLComputePipelineState,
        validatePipeline: MTLComputePipelineState,
        reactionBuffer: MTLBuffer,
        stoichiometryBuffer: MTLBuffer,
        stateA: MTLBuffer,
        stateB: MTLBuffer,
        statusBuffer: MTLBuffer,
        eventCountBuffer: MTLBuffer,
        readbackBuffer: MTLBuffer,
        stateByteCount: Int
    ) {
        self.model = model
        self.configuration = configuration
        self.device = device
        self.queue = queue
        self.advancePipeline = advancePipeline
        self.clearPipeline = clearPipeline
        self.validatePipeline = validatePipeline
        self.reactionBuffer = reactionBuffer
        self.stoichiometryBuffer = stoichiometryBuffer
        self.stateA = stateA
        self.stateB = stateB
        self.statusBuffer = statusBuffer
        self.eventCountBuffer = eventCountBuffer
        self.readbackBuffer = readbackBuffer
        self.stateByteCount = stateByteCount
        self.deviceName = device.name
        self.registryID = device.registryID
    }

    public func time() -> Float { absoluteTime }
    public func stepIndex() -> UInt32 { nextStepIndex }

    public func step(
        deltaTime requestedStep: Float? = nil,
        permitAdaptiveReduction: Bool = true
    ) async throws -> VivoExactSSAStepCertificate {
        let requested = requestedStep ?? configuration.timeStep
        guard requested.isFinite,
              requested >= configuration.minimumTimeStep,
              requested <= configuration.maximumTimeStep else {
            throw VivoExactSSAError.invalidConfiguration
        }

        var candidateStep = requested
        var lastStatus = VivoExactSSAStatus(
            flags: [],
            lanesNeedingSmallerStep: 0,
            invalidLaneCount: 0,
            totalEventCount: 0,
            maximumEventsObserved: 0,
            firstInvalidLane: nil,
            firstInvalidReaction: nil
        )

        for attempt in 0..<configuration.maximumAttempts {
            lastStatus = try await executeAttempt(deltaTime: candidateStep, attempt: attempt)
            let invalid = !lastStatus.flags.intersection(.invalid).isEmpty
            if invalid {
                return certificate(
                    disposition: .rejected,
                    requested: requested,
                    accepted: nil,
                    attempts: attempt + 1,
                    status: lastStatus
                )
            }
            if lastStatus.flags.contains(.needsSmallerStep) {
                let reduced = candidateStep * 0.5
                guard permitAdaptiveReduction,
                      attempt + 1 < configuration.maximumAttempts,
                      reduced >= configuration.minimumTimeStep,
                      reduced < candidateStep else {
                    return certificate(
                        disposition: .rejected,
                        requested: requested,
                        accepted: nil,
                        attempts: attempt + 1,
                        status: lastStatus
                    )
                }
                candidateStep = reduced
                continue
            }

            authoritativeIsA.toggle()
            let before = absoluteTime
            absoluteTime += candidateStep
            let completedStep = nextStepIndex
            nextStepIndex &+= 1
            return VivoExactSSAStepCertificate(
                disposition: candidateStep == requested ? .committed : .committedWithReducedStep,
                modelFingerprint: model.fingerprint,
                deviceName: deviceName,
                registryID: registryID,
                stepIndex: completedStep,
                requestedStep: requested,
                acceptedStep: candidateStep,
                attemptCount: attempt + 1,
                timeBefore: before,
                timeAfter: absoluteTime,
                status: lastStatus
            )
        }

        return certificate(
            disposition: .rejected,
            requested: requested,
            accepted: nil,
            attempts: configuration.maximumAttempts,
            status: lastStatus
        )
    }

    public func snapshot() async throws -> VivoExactSSAStateSnapshot {
        guard let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoExactSSAError.commandQueueUnavailable
        }
        command.label = "NumiVivo.ExactSSA.Snapshot"
        blit.copy(
            from: authoritativeState,
            sourceOffset: 0,
            to: readbackBuffer,
            destinationOffset: 0,
            size: stateByteCount
        )
        blit.endEncoding()
        try await complete(command)
        let pointer = readbackBuffer.contents().assumingMemoryBound(to: UInt32.self)
        let count = stateByteCount / MemoryLayout<UInt32>.stride
        let values = Array(UnsafeBufferPointer(start: pointer, count: count))
        return .init(
            modelFingerprint: model.fingerprint,
            stepIndex: nextStepIndex,
            absoluteTime: absoluteTime,
            speciesCount: UInt32(model.species.count),
            laneCount: configuration.laneCount,
            counts: values
        )
    }

    public func eventCounts() -> [UInt32] {
        let pointer = eventCountBuffer.contents().assumingMemoryBound(to: UInt32.self)
        return Array(UnsafeBufferPointer(start: pointer, count: Int(configuration.laneCount)))
    }

    private var authoritativeState: MTLBuffer { authoritativeIsA ? stateA : stateB }
    private var candidateState: MTLBuffer { authoritativeIsA ? stateB : stateA }

    private func executeAttempt(deltaTime: Float, attempt: UInt32) async throws -> VivoExactSSAStatus {
        guard let command = queue.makeCommandBuffer() else {
            throw VivoExactSSAError.commandQueueUnavailable
        }
        command.label = "NumiVivo.ExactSSA.Step.\(nextStepIndex).Attempt.\(attempt)"

        guard let blit = command.makeBlitCommandEncoder() else {
            throw VivoExactSSAError.commandQueueUnavailable
        }
        blit.copy(
            from: authoritativeState,
            sourceOffset: 0,
            to: candidateState,
            destinationOffset: 0,
            size: stateByteCount
        )
        blit.endEncoding()

        try encode(
            pipeline: clearPipeline,
            commandBuffer: command,
            threadCount: 1,
            label: "NumiVivo.ExactSSA.Clear"
        ) { encoder in
            encoder.setBuffer(statusBuffer, offset: 0, index: 0)
        }

        var runtimeCommand = VivoSSARuntimeCommandABI(
            laneCount: configuration.laneCount,
            speciesCount: UInt32(model.species.count),
            reactionCount: UInt32(model.reactions.count),
            maximumEventsPerLane: configuration.maximumEventsPerLane,
            stepIndex: nextStepIndex,
            attemptIndex: attempt,
            deltaTime: deltaTime,
            absoluteTime: absoluteTime + deltaTime,
            seed: configuration.seed
        )

        try encode(
            pipeline: advancePipeline,
            commandBuffer: command,
            threadCount: Int(configuration.laneCount),
            label: "NumiVivo.ExactSSA.Advance"
        ) { encoder in
            encoder.setBuffer(candidateState, offset: 0, index: 0)
            encoder.setBuffer(reactionBuffer, offset: 0, index: 1)
            encoder.setBuffer(stoichiometryBuffer, offset: 0, index: 2)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoSSARuntimeCommandABI>.stride, index: 3)
            encoder.setBuffer(statusBuffer, offset: 0, index: 4)
            encoder.setBuffer(eventCountBuffer, offset: 0, index: 5)
        }

        let stateElements = stateByteCount / MemoryLayout<UInt32>.stride
        try encode(
            pipeline: validatePipeline,
            commandBuffer: command,
            threadCount: stateElements,
            label: "NumiVivo.ExactSSA.Validate"
        ) { encoder in
            encoder.setBuffer(candidateState, offset: 0, index: 0)
            encoder.setBytes(&runtimeCommand, length: MemoryLayout<VivoSSARuntimeCommandABI>.stride, index: 1)
            encoder.setBuffer(statusBuffer, offset: 0, index: 2)
        }

        try await complete(command)
        let raw = statusBuffer.contents().load(as: VivoSSARuntimeStatusABI.self)
        return .init(
            flags: .init(rawValue: raw.flags),
            lanesNeedingSmallerStep: raw.lanesNeedingSmallerStep,
            invalidLaneCount: raw.invalidLaneCount,
            totalEventCount: raw.totalEventCount,
            maximumEventsObserved: raw.maximumEventsObserved,
            firstInvalidLane: raw.firstInvalidLane == UInt32.max ? nil : raw.firstInvalidLane,
            firstInvalidReaction: raw.firstInvalidReaction == UInt32.max ? nil : raw.firstInvalidReaction
        )
    }

    private func encode(
        pipeline: MTLComputePipelineState,
        commandBuffer: MTLCommandBuffer,
        threadCount: Int,
        label: String,
        bindings: (MTLComputeCommandEncoder) -> Void
    ) throws {
        guard threadCount > 0,
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw VivoExactSSAError.commandQueueUnavailable
        }
        encoder.label = label
        encoder.setComputePipelineState(pipeline)
        bindings(encoder)
        let width = Self.threadgroupWidth(pipeline)
        encoder.dispatchThreads(
            MTLSize(width: threadCount, height: 1, depth: 1),
            threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1)
        )
        encoder.endEncoding()
    }

    private func certificate(
        disposition: VivoExactSSAStepDisposition,
        requested: Float,
        accepted: Float?,
        attempts: UInt32,
        status: VivoExactSSAStatus
    ) -> VivoExactSSAStepCertificate {
        .init(
            disposition: disposition,
            modelFingerprint: model.fingerprint,
            deviceName: deviceName,
            registryID: registryID,
            stepIndex: nextStepIndex,
            requestedStep: requested,
            acceptedStep: accepted,
            attemptCount: attempts,
            timeBefore: absoluteTime,
            timeAfter: absoluteTime,
            status: status
        )
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { completed in
                if completed.status == .completed {
                    continuation.resume()
                } else {
                    continuation.resume(throwing: VivoExactSSAError.commandFailed(
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
        guard !values.isEmpty else { throw VivoExactSSAError.allocationFailed("\(label) is empty") }
        let byteCount = values.count.multipliedReportingOverflow(by: MemoryLayout<T>.stride)
        guard !byteCount.overflow else { throw VivoExactSSAError.arithmeticOverflow }
        guard let staging = device.makeBuffer(length: byteCount.partialValue, options: .storageModeShared),
              let destination = device.makeBuffer(length: byteCount.partialValue, options: .storageModePrivate),
              let command = queue.makeCommandBuffer(),
              let blit = command.makeBlitCommandEncoder() else {
            throw VivoExactSSAError.allocationFailed(label)
        }
        values.withUnsafeBytes { bytes in
            staging.contents().copyMemory(from: bytes.baseAddress!, byteCount: byteCount.partialValue)
        }
        staging.label = "\(label).Staging"
        destination.label = label
        blit.copy(from: staging, sourceOffset: 0, to: destination, destinationOffset: 0, size: byteCount.partialValue)
        blit.endEncoding()
        command.commit()
        command.waitUntilCompleted()
        guard command.status == .completed else {
            throw VivoExactSSAError.commandFailed(command.error?.localizedDescription ?? "initial upload failed")
        }
        return destination
    }
}

private extension Optional {
    func orThrow(_ error: @autoclosure () -> Error) throws -> Wrapped {
        guard let value = self else { throw error() }
        return value
    }
}
