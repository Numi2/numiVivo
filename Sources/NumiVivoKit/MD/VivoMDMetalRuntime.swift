import Foundation
import Metal
import NumiVivoShaders

public enum VivoMDRuntimeError: Error, Sendable, CustomStringConvertible {
    case unsupported([String])
    case metal(String)
    case rejected(flags: UInt32, firstParticle: UInt32, violations: UInt32)

    public var description: String {
        switch self {
        case .unsupported(let reasons): return "MD execution is unsupported: " + reasons.joined(separator: "; ")
        case .metal(let message): return "MD Metal runtime: \(message)"
        case .rejected(let flags, let particle, let count):
            return "MD candidate rejected (flags=\(flags), firstParticle=\(particle), violations=\(count))"
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
    public let violationCount: UInt32
}

private struct VivoMDCommandABI {
    var particleCount: UInt32
    var typeCount: UInt32
    var electrostatics: UInt32
    var periodic: UInt32
    var dtPS: Float
    var cutoffNM: Float
    var coulombPrefactor: Float
    var reactionFieldK: Float
    var reactionFieldC: Float
    var minimumDistanceNM: Float
    var reserved0: Float = 0
    var reserved1: Float = 0
    var cellA: SIMD4<Float>
    var cellB: SIMD4<Float>
    var cellC: SIMD4<Float>
    var reciprocalA: SIMD4<Float>
    var reciprocalB: SIMD4<Float>
    var reciprocalC: SIMD4<Float>
}

private struct VivoMDStatusABI {
    var flags: UInt32
    var firstParticle: UInt32
    var violationCount: UInt32
    var reserved: UInt32
}

private final class VivoMDMetalArena: @unchecked Sendable {
    let device: MTLDevice
    let particleCount: Int
    let positionA: MTLBuffer
    let positionB: MTLBuffer
    let velocityA: MTLBuffer
    let velocityB: MTLBuffer
    let forceEnergy: MTLBuffer
    let kinetic: MTLBuffer
    let status: MTLBuffer
    let dynamics: MTLBuffer
    let typeIndices: MTLBuffer
    let pairMatrix: MTLBuffer
    let bonds: MTLBuffer
    let bondOffsets: MTLBuffer
    let bondIncidence: MTLBuffer
    let angles: MTLBuffer
    let angleOffsets: MTLBuffer
    let angleIncidence: MTLBuffer
    let torsions: MTLBuffer
    let torsionOffsets: MTLBuffer
    let torsionIncidence: MTLBuffer
    let exceptions: MTLBuffer
    let exceptionOffsets: MTLBuffer
    let exceptionPartners: MTLBuffer
    let exceptionIndices: MTLBuffer
    let positionReadback: MTLBuffer
    let velocityReadback: MTLBuffer
    private(set) var acceptedIsA = true

    var acceptedPosition: MTLBuffer { acceptedIsA ? positionA : positionB }
    var candidatePosition: MTLBuffer { acceptedIsA ? positionB : positionA }
    var acceptedVelocity: MTLBuffer { acceptedIsA ? velocityA : velocityB }
    var candidateVelocity: MTLBuffer { acceptedIsA ? velocityB : velocityA }
    func commit() { acceptedIsA.toggle() }

    static func make(device: MTLDevice, queue: MTLCommandQueue,
                     packed: VivoMDPackedSystem,
                     initial: VivoClassicalInitialState,
                     velocities: [VivoVector3D]) async throws -> VivoMDMetalArena {
        let count = Int(packed.particleCount)
        guard initial.positionsNM.count == count, velocities.count == count else {
            throw VivoMDRuntimeError.metal("initial position/velocity shape mismatch")
        }
        let positions = try initial.positionsNM.map { try float4($0, label: "initial position") }
        let velocity4 = try velocities.map { try float4($0, label: "initial velocity") }
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("cannot create initialization command")
        }
        var staging: [MTLBuffer] = []
        func privateBuffer<T>(_ values: [T], label: String) throws -> MTLBuffer {
            let byteCount = max(values.count * MemoryLayout<T>.stride, 16)
            guard let destination = device.makeBuffer(length: byteCount, options: .storageModePrivate) else {
                throw VivoMDRuntimeError.metal("allocation failed for \(label)")
            }
            destination.label = label
            if !values.isEmpty {
                let source = values.withUnsafeBytes { raw -> MTLBuffer? in
                    guard let base = raw.baseAddress else { return nil }
                    return device.makeBuffer(bytes: base, length: raw.count, options: .storageModeShared)
                }
                guard let source else { throw VivoMDRuntimeError.metal("staging allocation failed for \(label)") }
                staging.append(source)
                blit.copy(from: source, sourceOffset: 0, to: destination, destinationOffset: 0,
                          size: values.count * MemoryLayout<T>.stride)
            }
            return destination
        }
        let positionA = try privateBuffer(positions, label: "NumiVivo.MD.positionA")
        let positionB = try privateBuffer(positions, label: "NumiVivo.MD.positionB")
        let velocityA = try privateBuffer(velocity4, label: "NumiVivo.MD.velocityA")
        let velocityB = try privateBuffer(velocity4, label: "NumiVivo.MD.velocityB")
        let dynamics = try privateBuffer(packed.particleDynamics, label: "NumiVivo.MD.dynamics")
        let typeIndices = try privateBuffer(packed.particleTypeIndices, label: "NumiVivo.MD.typeIndices")
        let pairMatrix = try privateBuffer(packed.typePairC12C6, label: "NumiVivo.MD.pairMatrix")
        let bonds = try privateBuffer(packed.bonds, label: "NumiVivo.MD.bonds")
        let bondOffsets = try privateBuffer(packed.bondOffsets, label: "NumiVivo.MD.bondOffsets")
        let bondIncidence = try privateBuffer(packed.bondIncidence, label: "NumiVivo.MD.bondIncidence")
        let angles = try privateBuffer(packed.angles, label: "NumiVivo.MD.angles")
        let angleOffsets = try privateBuffer(packed.angleOffsets, label: "NumiVivo.MD.angleOffsets")
        let angleIncidence = try privateBuffer(packed.angleIncidence, label: "NumiVivo.MD.angleIncidence")
        let torsions = try privateBuffer(packed.torsions, label: "NumiVivo.MD.torsions")
        let torsionOffsets = try privateBuffer(packed.torsionOffsets, label: "NumiVivo.MD.torsionOffsets")
        let torsionIncidence = try privateBuffer(packed.torsionIncidence, label: "NumiVivo.MD.torsionIncidence")
        let exceptions = try privateBuffer(packed.pairExceptions, label: "NumiVivo.MD.exceptions")
        let exceptionOffsets = try privateBuffer(packed.exceptionOffsets, label: "NumiVivo.MD.exceptionOffsets")
        let exceptionPartners = try privateBuffer(packed.exceptionPartners, label: "NumiVivo.MD.exceptionPartners")
        let exceptionIndices = try privateBuffer(packed.exceptionIndices, label: "NumiVivo.MD.exceptionIndices")
        blit.endEncoding()
        try await complete(command)
        _ = staging
        guard let forceEnergy = device.makeBuffer(length: max(count * MemoryLayout<SIMD4<Float>>.stride, 16), options: .storageModePrivate),
              let kinetic = device.makeBuffer(length: max(count * MemoryLayout<Float>.stride, 16), options: .storageModePrivate),
              let status = device.makeBuffer(length: MemoryLayout<VivoMDStatusABI>.stride, options: .storageModeShared),
              let positionReadback = device.makeBuffer(length: max(count * MemoryLayout<SIMD4<Float>>.stride, 16), options: .storageModeShared),
              let velocityReadback = device.makeBuffer(length: max(count * MemoryLayout<SIMD4<Float>>.stride, 16), options: .storageModeShared) else {
            throw VivoMDRuntimeError.metal("dynamic/readback buffer allocation failed")
        }
        return VivoMDMetalArena(device: device, particleCount: count,
                                positionA: positionA, positionB: positionB,
                                velocityA: velocityA, velocityB: velocityB,
                                forceEnergy: forceEnergy, kinetic: kinetic, status: status,
                                dynamics: dynamics, typeIndices: typeIndices, pairMatrix: pairMatrix,
                                bonds: bonds, bondOffsets: bondOffsets, bondIncidence: bondIncidence,
                                angles: angles, angleOffsets: angleOffsets, angleIncidence: angleIncidence,
                                torsions: torsions, torsionOffsets: torsionOffsets, torsionIncidence: torsionIncidence,
                                exceptions: exceptions, exceptionOffsets: exceptionOffsets,
                                exceptionPartners: exceptionPartners, exceptionIndices: exceptionIndices,
                                positionReadback: positionReadback, velocityReadback: velocityReadback)
    }

    private init(device: MTLDevice, particleCount: Int,
                 positionA: MTLBuffer, positionB: MTLBuffer, velocityA: MTLBuffer, velocityB: MTLBuffer,
                 forceEnergy: MTLBuffer, kinetic: MTLBuffer, status: MTLBuffer,
                 dynamics: MTLBuffer, typeIndices: MTLBuffer, pairMatrix: MTLBuffer,
                 bonds: MTLBuffer, bondOffsets: MTLBuffer, bondIncidence: MTLBuffer,
                 angles: MTLBuffer, angleOffsets: MTLBuffer, angleIncidence: MTLBuffer,
                 torsions: MTLBuffer, torsionOffsets: MTLBuffer, torsionIncidence: MTLBuffer,
                 exceptions: MTLBuffer, exceptionOffsets: MTLBuffer, exceptionPartners: MTLBuffer,
                 exceptionIndices: MTLBuffer, positionReadback: MTLBuffer, velocityReadback: MTLBuffer) {
        self.device = device; self.particleCount = particleCount
        self.positionA = positionA; self.positionB = positionB
        self.velocityA = velocityA; self.velocityB = velocityB
        self.forceEnergy = forceEnergy; self.kinetic = kinetic; self.status = status
        self.dynamics = dynamics; self.typeIndices = typeIndices; self.pairMatrix = pairMatrix
        self.bonds = bonds; self.bondOffsets = bondOffsets; self.bondIncidence = bondIncidence
        self.angles = angles; self.angleOffsets = angleOffsets; self.angleIncidence = angleIncidence
        self.torsions = torsions; self.torsionOffsets = torsionOffsets; self.torsionIncidence = torsionIncidence
        self.exceptions = exceptions; self.exceptionOffsets = exceptionOffsets
        self.exceptionPartners = exceptionPartners; self.exceptionIndices = exceptionIndices
        self.positionReadback = positionReadback; self.velocityReadback = velocityReadback
    }

    private static func float4(_ value: VivoVector3D, label: String) throws -> SIMD4<Float> {
        guard value.isFinite,
              abs(value.x) <= Double(Float.greatestFiniteMagnitude),
              abs(value.y) <= Double(Float.greatestFiniteMagnitude),
              abs(value.z) <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoMDRuntimeError.metal("\(label) is outside FP32 range")
        }
        return .init(Float(value.x), Float(value.y), Float(value.z), 0)
    }
    private static func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { value in
                if let error = value.error { continuation.resume(throwing: VivoMDRuntimeError.metal(String(describing: error))) }
                else { continuation.resume(returning: ()) }
            }
            command.commit()
        }
    }
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
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let periodicCell: VivoPeriodicCell?
    private let arena: VivoMDMetalArena
    private var acceptedStep: UInt64 = 0
    private var acceptedTimePS: Double = 0
    private var inFlight = false

    public static func make(system: VivoClassicalSystem,
                            initialState: VivoClassicalInitialState,
                            configuration: VivoMDConfiguration,
                            initialVelocitiesNMPerPS: [VivoVector3D]? = nil,
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        let report = try VivoMDCapabilityAnalyzer.analyze(system: system, initialState: initialState, configuration: configuration)
        var blockers = report.blockers
        if configuration.electrostatics == .pme { blockers.append("PME reciprocal-space execution is not yet present in the first Wave B kernel set") }
        if !system.constraints.isEmpty { blockers.append("constraint projection is not yet present in the first Wave B kernel set") }
        if configuration.ensemble != .nve || configuration.thermostat != .none || configuration.barostat != .none {
            blockers.append("the first Wave B kernel set currently executes NVE only")
        }
        guard blockers.isEmpty else { throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted()) }
        let packed = try VivoMDSystemPacker.pack(system)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoMDRuntimeError.metal("Apple-silicon unified memory and a command queue are required")
        }
        queue.label = "NumiVivo.MD.Queue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        let kernels: [NumiVivoKernel] = [.mdClear, .mdBonded, .mdNonbondedDirect, .mdHalfKick, .mdDrift, .mdKinetic, .mdValidate]
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in kernels { pipelines[kernel] = try await catalog.pipeline(kernel) }
        let velocities = initialVelocitiesNMPerPS ?? [VivoVector3D](repeating: .zero, count: system.particles.count)
        let arena = try await VivoMDMetalArena.make(device: device, queue: queue, packed: packed,
                                                    initial: initialState, velocities: velocities)
        return try VivoMDMetalRuntime(system: system, configuration: configuration, packed: packed,
                                      queue: queue, arena: arena, periodicCell: initialState.periodicCell,
                                      device: device)
    }

    private init(system: VivoClassicalSystem, configuration: VivoMDConfiguration,
                 packed: VivoMDPackedSystem, queue: MTLCommandQueue,
                 arena: VivoMDMetalArena, periodicCell: VivoPeriodicCell?, device: MTLDevice) throws {
        self.system = system; self.configuration = configuration; self.packed = packed
        self.queue = queue; self.arena = arena; self.periodicCell = periodicCell
        systemFingerprint = packed.systemFingerprint
        configurationFingerprint = try configuration.fingerprint()
        deviceName = device.name; deviceRegistryID = device.registryID
        let catalog = try? NumiVivoPipelineCatalog(device: device)
        _ = catalog
        // Pipelines are installed by factory immediately below through the local helper.
        self.pipelines = [:]
        throw VivoMDRuntimeError.metal("internal factory initializer must install pipelines")
    }

    private init(system: VivoClassicalSystem, configuration: VivoMDConfiguration,
                 packed: VivoMDPackedSystem, queue: MTLCommandQueue, arena: VivoMDMetalArena,
                 periodicCell: VivoPeriodicCell?, device: MTLDevice,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline]) throws {
        self.system = system; self.configuration = configuration; self.packed = packed
        self.queue = queue; self.arena = arena; self.periodicCell = periodicCell
        self.pipelines = pipelines
        systemFingerprint = packed.systemFingerprint
        configurationFingerprint = try configuration.fingerprint()
        deviceName = device.name; deviceRegistryID = device.registryID
    }

    public func step() async throws -> VivoMDStepCertificate {
        guard !inFlight else { throw VivoMDRuntimeError.metal("MD operation already in flight") }
        guard acceptedStep < UInt64.max else { throw VivoMDRuntimeError.metal("MD step index overflow") }
        inFlight = true
        defer { inFlight = false }
        let timeBefore = acceptedTimePS
        let commandABI = try makeCommand()
        guard let command = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("command buffer unavailable") }
        command.label = "NumiVivo.MD.step.\(acceptedStep)"
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("blit encoder unavailable") }
        let vectorBytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0, to: arena.candidatePosition, destinationOffset: 0, size: vectorBytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0, to: arena.candidateVelocity, destinationOffset: 0, size: vectorBytes)
        blit.endEncoding()
        try encodeForces(command, position: arena.candidatePosition, abi: commandABI)
        try encode(.mdHalfKick, command: command, buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: commandABI)
        try encode(.mdDrift, command: command, buffers: [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi: commandABI)
        try encodeForces(command, position: arena.candidatePosition, abi: commandABI)
        try encode(.mdHalfKick, command: command, buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: commandABI)
        try encode(.mdValidate, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity, arena.forceEnergy, arena.status], abi: commandABI)
        try await complete(command)
        let status = arena.status.contents().assumingMemoryBound(to: VivoMDStatusABI.self).pointee
        guard status.flags == 0 else {
            return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                         deviceName: deviceName, deviceRegistryID: deviceRegistryID, stepIndex: acceptedStep,
                         timeBeforePS: timeBefore, timeAfterPS: timeBefore, committed: false,
                         statusFlags: status.flags, violationCount: status.violationCount)
        }
        arena.commit()
        acceptedStep += 1
        acceptedTimePS += configuration.timeStepPS
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                     deviceName: deviceName, deviceRegistryID: deviceRegistryID, stepIndex: acceptedStep - 1,
                     timeBeforePS: timeBefore, timeAfterPS: acceptedTimePS, committed: true,
                     statusFlags: 0, violationCount: 0)
    }

    public func snapshot() async throws -> VivoMDStateSnapshot {
        guard !inFlight else { throw VivoMDRuntimeError.metal("cannot snapshot while MD operation is in flight") }
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("snapshot command unavailable")
        }
        let bytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0, to: arena.positionReadback, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0, to: arena.velocityReadback, destinationOffset: 0, size: bytes)
        blit.endEncoding()
        try await complete(command)
        let p = arena.positionReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let v = arena.velocityReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
        positions.reserveCapacity(arena.particleCount); velocities.reserveCapacity(arena.particleCount)
        for i in 0..<arena.particleCount {
            positions.append(.init(Double(p[i].x), Double(p[i].y), Double(p[i].z)))
            velocities.append(.init(Double(v[i].x), Double(v[i].y), Double(v[i].z)))
        }
        let result = VivoMDStateSnapshot(systemFingerprint: systemFingerprint,
                                         configurationFingerprint: configurationFingerprint,
                                         stepIndex: acceptedStep, timePS: acceptedTimePS,
                                         positionsNM: positions, velocitiesNMPerPS: velocities,
                                         periodicCell: periodicCell)
        try result.validate(particleCount: arena.particleCount)
        return result
    }

    private func encodeForces(_ command: MTLCommandBuffer, position: MTLBuffer, abi: VivoMDCommandABI) throws {
        try encode(.mdClear, command: command, buffers: [arena.forceEnergy, arena.status], abi: abi)
        try encode(.mdBonded, command: command,
                   buffers: [position, arena.forceEnergy,
                             arena.bonds, arena.bondOffsets, arena.bondIncidence,
                             arena.angles, arena.angleOffsets, arena.angleIncidence,
                             arena.torsions, arena.torsionOffsets, arena.torsionIncidence,
                             arena.status], abi: abi)
        try encode(.mdNonbondedDirect, command: command,
                   buffers: [position, arena.forceEnergy, arena.dynamics, arena.typeIndices,
                             arena.pairMatrix, arena.exceptions, arena.exceptionOffsets,
                             arena.exceptionPartners, arena.exceptionIndices, arena.status], abi: abi)
    }

    private func encode(_ kernel: NumiVivoKernel, command: MTLCommandBuffer,
                        buffers: [MTLBuffer], abi: VivoMDCommandABI) throws {
        guard let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("pipeline or encoder unavailable for \(kernel.rawValue)")
        }
        encoder.label = kernel.rawValue
        encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        var commandABI = abi
        encoder.setBytes(&commandABI, length: MemoryLayout<VivoMDCommandABI>.stride, index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
        encoder.endEncoding()
    }

    private func makeCommand() throws -> VivoMDCommandABI {
        let electrostatics: UInt32
        switch configuration.electrostatics {
        case .cutoff: electrostatics = 0
        case .reactionField: electrostatics = 1
        case .pme: throw VivoMDRuntimeError.unsupported(["PME kernel is not installed"])
        }
        let cutoff = Float(configuration.cutoffNM)
        let eps = configuration.relativeDielectric
        let prefactor = Float(138.935456 / eps)
        var krf: Float = 0, crf: Float = 0
        if configuration.electrostatics == .reactionField {
            let er = configuration.reactionFieldDielectric
            krf = Float((er - eps) / (2 * er + eps) / pow(configuration.cutoffNM, 3))
            crf = Float(3 * er / (2 * er + eps) / configuration.cutoffNM)
        }
        let cell = periodicCell
        let reciprocal = try cell.map(reciprocalCell)
        func f4(_ value: VivoVector3D?) -> SIMD4<Float> {
            guard let value else { return .zero }
            return .init(Float(value.x), Float(value.y), Float(value.z), 0)
        }
        return .init(particleCount: packed.particleCount, typeCount: packed.typeCount,
                     electrostatics: electrostatics, periodic: cell == nil ? 0 : 1,
                     dtPS: Float(configuration.timeStepPS), cutoffNM: cutoff,
                     coulombPrefactor: prefactor, reactionFieldK: krf, reactionFieldC: crf,
                     minimumDistanceNM: 1e-5,
                     cellA: f4(cell?.a), cellB: f4(cell?.b), cellC: f4(cell?.c),
                     reciprocalA: f4(reciprocal?.0), reciprocalB: f4(reciprocal?.1), reciprocalC: f4(reciprocal?.2))
    }

    private func reciprocalCell(_ cell: VivoPeriodicCell) throws -> (VivoVector3D,VivoVector3D,VivoVector3D) {
        let determinant = cell.a.dot(cell.b.cross(cell.c))
        guard determinant.isFinite, abs(determinant) > 1e-18 else {
            throw VivoMDRuntimeError.metal("periodic cell is singular")
        }
        return (cell.b.cross(cell.c) / determinant,
                cell.c.cross(cell.a) / determinant,
                cell.a.cross(cell.b) / determinant)
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation { continuation in
            command.addCompletedHandler { value in
                if let error = value.error { continuation.resume(throwing: VivoMDRuntimeError.metal(String(describing: error))) }
                else { continuation.resume(returning: ()) }
            }
            command.commit()
        }
    }
}
