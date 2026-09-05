import Foundation
import Metal
import NumiVivoShaders

public enum VivoMDRuntimeError: Error, Sendable, CustomStringConvertible {
    case unsupported([String])
    case metal(String)
    public var description: String {
        switch self {
        case .unsupported(let values): return "MD execution is unsupported: " + values.joined(separator: "; ")
        case .metal(let value): return "MD Metal runtime: \(value)"
        }
    }
}

public struct VivoMDStepCertificate: Codable, Sendable, Equatable {
    public let systemFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let stepIndex: UInt64
    public let timeBeforePS: Double
    public let timeAfterPS: Double
    public let committed: Bool
    public let statusFlags: UInt32
    public let firstViolationParticle: UInt32?
    public let violationCount: UInt32
}

public actor VivoMDMetalRuntime {
    public nonisolated let system: VivoClassicalSystem
    public nonisolated let configuration: VivoMDConfiguration
    public nonisolated let systemFingerprint: VivoFingerprint
    public nonisolated let configurationFingerprint: VivoFingerprint
    public nonisolated let deviceName: String
    public nonisolated let deviceRegistryID: UInt64

    private let queue: MTLCommandQueue
    private let packed: VivoMDPackedSystem
    private let arena: VivoMDGPUArena
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let periodicCell: VivoPeriodicCell?
    private var acceptedStep: UInt64 = 0
    private var acceptedTimePS: Double = 0
    private var inFlight = false
    private var acceptedForceIsCurrent = false

    public static func make(system: VivoClassicalSystem,
                            initialState: VivoClassicalInitialState,
                            configuration: VivoMDConfiguration,
                            initialVelocitiesNMPerPS: [VivoVector3D]? = nil,
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        let report = try VivoMDCapabilityAnalyzer.analyze(system: system, initialState: initialState, configuration: configuration)
        var blockers = report.blockers
        if configuration.electrostatics == .pme { blockers.append("PME reciprocal-space execution is not installed yet") }
        if configuration.ensemble != .nve || configuration.thermostat != .none || configuration.barostat != .none {
            blockers.append("current Wave B stage executes NVE only")
        }
        if let cell = initialState.periodicCell {
            try VivoMDMetalABI.validateMinimumImageCutoff(configuration.cutoffNM, cell: cell)
        }
        guard blockers.isEmpty else { throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted()) }
        let packed = try VivoMDSystemPacker.pack(system)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoMDRuntimeError.metal("Apple-silicon unified memory and command queue are required")
        }
        queue.label = "NumiVivo.MD.Queue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let required: [NumiVivoKernel] = [
            .mdClearForce, .mdClearStatus, .mdBonded, .mdNonbondedDirect,
            .mdHalfKick, .mdDrift, .mdConstraintPosition, .mdConstraintVelocity,
            .mdValidateConstraints, .mdKinetic, .mdValidate
        ]
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in required { pipelines[kernel] = try await catalog.pipeline(kernel) }
        let velocities = initialVelocitiesNMPerPS ?? [VivoVector3D](repeating: .zero, count: system.particles.count)
        let arena = try await VivoMDGPUArena.make(device: device, queue: queue, packed: packed,
                                                  initial: initialState, velocities: velocities)
        return try .init(system: system, configuration: configuration, packed: packed,
                         queue: queue, arena: arena, pipelines: pipelines,
                         periodicCell: initialState.periodicCell, device: device)
    }

    private init(system: VivoClassicalSystem, configuration: VivoMDConfiguration,
                 packed: VivoMDPackedSystem, queue: MTLCommandQueue,
                 arena: VivoMDGPUArena, pipelines: [NumiVivoKernel: NumiVivoPipeline],
                 periodicCell: VivoPeriodicCell?, device: MTLDevice) throws {
        self.system = system; self.configuration = configuration; self.packed = packed
        self.queue = queue; self.arena = arena; self.pipelines = pipelines; self.periodicCell = periodicCell
        systemFingerprint = packed.systemFingerprint
        configurationFingerprint = try configuration.fingerprint()
        deviceName = device.name; deviceRegistryID = device.registryID
    }

    public func step() async throws -> VivoMDStepCertificate {
        guard !inFlight else { throw VivoMDRuntimeError.metal("MD operation already in flight") }
        guard acceptedStep < UInt64.max else { throw VivoMDRuntimeError.metal("MD step index overflow") }
        inFlight = true; defer { inFlight = false }
        let before = acceptedTimePS, abi = try VivoMDMetalABI.command(packed: packed, configuration: configuration, cell: periodicCell)
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("step command unavailable")
        }
        command.label = "NumiVivo.MD.step.\(acceptedStep)"
        let bytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0, to: arena.candidatePosition, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0, to: arena.candidateVelocity, destinationOffset: 0, size: bytes)
        blit.endEncoding()

        try encodeStatusClear(command)
        try encodeForces(command, position: arena.candidatePosition, abi: abi)
        try encode(.mdHalfKick, command: command,
                   buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: abi)
        try encode(.mdDrift, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi: abi)
        if !packed.constraints.isEmpty { try encodePositionConstraints(command, abi: abi) }
        try encodeForces(command, position: arena.candidatePosition, abi: abi)
        try encode(.mdHalfKick, command: command,
                   buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: abi)
        if !packed.constraints.isEmpty {
            try encodeVelocityConstraints(command, abi: abi)
            try encode(.mdValidateConstraints, command: command,
                       buffers: [arena.candidatePosition, arena.candidateVelocity,
                                 arena.constraints, arena.constraintOffsets, arena.constraintIncidence, arena.status], abi: abi)
        }
        try encode(.mdValidate, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity, arena.forceEnergy, arena.status], abi: abi)
        try await complete(command)

        let status = arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee
        let first = status.firstParticle == UInt32.max ? nil : status.firstParticle
        guard status.flags == 0 else {
            acceptedForceIsCurrent = false
            return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                         deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                         stepIndex: acceptedStep, timeBeforePS: before, timeAfterPS: before,
                         committed: false, statusFlags: status.flags,
                         firstViolationParticle: first, violationCount: status.violationCount)
        }
        arena.commit(); acceptedStep += 1; acceptedTimePS += configuration.timeStepPS
        acceptedForceIsCurrent = true
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                     deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                     stepIndex: acceptedStep - 1, timeBeforePS: before, timeAfterPS: acceptedTimePS,
                     committed: true, statusFlags: 0, firstViolationParticle: nil, violationCount: 0)
    }

    public func snapshot() async throws -> VivoMDStateSnapshot {
        guard !inFlight else { throw VivoMDRuntimeError.metal("cannot snapshot during MD operation") }
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("snapshot command unavailable")
        }
        let bytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0, to: arena.positionReadback, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0, to: arena.velocityReadback, destinationOffset: 0, size: bytes)
        blit.endEncoding(); try await complete(command)
        let pp = arena.positionReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let vv = arena.velocityReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
        positions.reserveCapacity(arena.particleCount); velocities.reserveCapacity(arena.particleCount)
        for i in 0..<arena.particleCount {
            positions.append(.init(Double(pp[i].x), Double(pp[i].y), Double(pp[i].z)))
            velocities.append(.init(Double(vv[i].x), Double(vv[i].y), Double(vv[i].z)))
        }
        let value = VivoMDStateSnapshot(systemFingerprint: systemFingerprint,
                                        configurationFingerprint: configurationFingerprint,
                                        stepIndex: acceptedStep, timePS: acceptedTimePS,
                                        positionsNM: positions, velocitiesNMPerPS: velocities,
                                        periodicCell: periodicCell)
        try value.validate(particleCount: arena.particleCount); return value
    }

    private func encodePositionConstraints(_ command: MTLCommandBuffer, abi: VivoMDMetalCommand) throws {
        var source = arena.candidatePosition, destination = arena.positionScratch
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintPosition, command: command,
                       buffers: [source, destination, arena.dynamics,
                                 arena.constraints, arena.constraintOffsets, arena.constraintIncidence, arena.status], abi: abi)
            swap(&source, &destination)
        }
        if source !== arena.candidatePosition { try copy(command, from: source, to: arena.candidatePosition) }
    }

    private func encodeVelocityConstraints(_ command: MTLCommandBuffer, abi: VivoMDMetalCommand) throws {
        var source = arena.candidateVelocity, destination = arena.velocityScratch
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintVelocity, command: command,
                       buffers: [arena.candidatePosition, source, destination, arena.dynamics,
                                 arena.constraints, arena.constraintOffsets, arena.constraintIncidence, arena.status], abi: abi)
            swap(&source, &destination)
        }
        if source !== arena.candidateVelocity { try copy(command, from: source, to: arena.candidateVelocity) }
    }

    private func encodeForces(_ command: MTLCommandBuffer, position: MTLBuffer, abi: VivoMDMetalCommand) throws {
        try encode(.mdClearForce, command: command, buffers: [arena.forceEnergy], abi: abi)
        try encode(.mdBonded, command: command,
                   buffers: [position, arena.forceEnergy, arena.bonds, arena.bondOffsets, arena.bondIncidence,
                             arena.angles, arena.angleOffsets, arena.angleIncidence,
                             arena.torsions, arena.torsionOffsets, arena.torsionIncidence, arena.status], abi: abi)
        try encode(.mdNonbondedDirect, command: command,
                   buffers: [position, arena.forceEnergy, arena.dynamics, arena.typeIndices, arena.pairMatrix,
                             arena.exceptions, arena.exceptionOffsets, arena.exceptionPartners,
                             arena.exceptionIndices, arena.status], abi: abi)
    }

    private func encodeStatusClear(_ command: MTLCommandBuffer) throws {
        guard let pipeline = pipelines[.mdClearStatus], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("status-clear encoder unavailable")
        }
        encoder.setComputePipelineState(pipeline.state); encoder.setBuffer(arena.status, offset: 0, index: 0)
        encoder.dispatchThreads(.init(width: 1, height: 1, depth: 1), threadsPerThreadgroup: .init(width: 1, height: 1, depth: 1)); encoder.endEncoding()
    }

    private func encode(_ kernel: NumiVivoKernel, command: MTLCommandBuffer,
                        buffers: [MTLBuffer], abi: VivoMDMetalCommand) throws {
        guard let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("encoder unavailable for \(kernel.rawValue)")
        }
        encoder.label = kernel.rawValue; encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        var value = abi; encoder.setBytes(&value, length: MemoryLayout<VivoMDMetalCommand>.stride, index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount)); encoder.endEncoding()
    }

    private func copy(_ command: MTLCommandBuffer, from: MTLBuffer, to: MTLBuffer) throws {
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("constraint copy encoder unavailable") }
        blit.copy(from: from, sourceOffset: 0, to: to, destinationOffset: 0,
                  size: arena.particleCount * MemoryLayout<SIMD4<Float>>.stride); blit.endEncoding()
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
            command.addCompletedHandler { value in
                if let error = value.error { continuation.resume(throwing: VivoMDRuntimeError.metal(String(describing: error))) }
                else { continuation.resume(returning: ()) }
            }
            command.commit()
        }
    }
}
