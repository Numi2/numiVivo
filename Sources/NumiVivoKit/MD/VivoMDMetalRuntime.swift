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
    private let neighborGrid: VivoMDNeighborGrid?
    private let pmeEngine: VivoPMEEngine?
    private var acceptedStep: UInt64
    private var acceptedTimePS: Double
    private var inFlight = false

    public static func make(system: VivoClassicalSystem,
                            initialState: VivoClassicalInitialState,
                            configuration: VivoMDConfiguration,
                            initialVelocitiesNMPerPS: [VivoVector3D]? = nil,
                            device requestedDevice: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        try await makeInternal(system: system, initialState: initialState,
                               configuration: configuration,
                               velocities: initialVelocitiesNMPerPS,
                               startStep: 0,
                               startTimePS: initialState.sourceTimePS ?? 0,
                               device: requestedDevice)
    }

    public static func restore(system: VivoClassicalSystem,
                               configuration: VivoMDConfiguration,
                               checkpoint: VivoMDCheckpoint,
                               device requestedDevice: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        let systemFingerprint = try system.fingerprint()
        let configurationFingerprint = try configuration.fingerprint()
        try checkpoint.validate(particleCount: system.particles.count)
        guard checkpoint.systemFingerprint == systemFingerprint,
              checkpoint.configurationFingerprint == configurationFingerprint else {
            throw VivoArtifactValidationError.incompatible(
                "MD checkpoint system or configuration identity differs from requested runtime"
            )
        }
        let initial = VivoClassicalInitialState(systemFingerprint: systemFingerprint,
                                                positionsNM: checkpoint.positionsNM,
                                                periodicCell: checkpoint.periodicCell,
                                                sourceTimePS: checkpoint.timePS)
        return try await makeInternal(system: system, initialState: initial,
                                      configuration: configuration,
                                      velocities: checkpoint.velocitiesNMPerPS,
                                      startStep: checkpoint.acceptedStep,
                                      startTimePS: checkpoint.timePS,
                                      device: requestedDevice)
    }

    private static func makeInternal(system: VivoClassicalSystem,
                                     initialState: VivoClassicalInitialState,
                                     configuration: VivoMDConfiguration,
                                     velocities: [VivoVector3D]?,
                                     startStep: UInt64,
                                     startTimePS: Double,
                                     device requestedDevice: MTLDevice?) async throws -> VivoMDMetalRuntime {
        let report = try VivoMDCapabilityAnalyzer.analyze(system: system,
                                                          initialState: initialState,
                                                          configuration: configuration)
        var blockers = report.blockers
        if configuration.electrostatics == .pme, !configuration.resolvedNeighborListEnabled {
            blockers.append("PME real-space execution requires the bounded neighbor-list path")
        }
        if configuration.ensemble == .npt || configuration.barostat != .none {
            blockers.append("NPT barostat execution is not installed yet")
        }
        if configuration.thermostat == .velocityRescale {
            blockers.append("stochastic velocity-rescale thermostat is not installed yet")
        }
        if configuration.ensemble == .nvt, configuration.thermostat != .langevinMiddle {
            blockers.append("current NVT execution requires langevinMiddle")
        }
        if configuration.ensemble == .nve, configuration.thermostat != .none {
            blockers.append("NVE cannot use a thermostat")
        }
        if let cell = initialState.periodicCell {
            try VivoMDMetalABI.validateMinimumImageCutoff(configuration.cutoffNM, cell: cell)
            if configuration.resolvedNeighborListEnabled {
                try VivoMDMetalABI.validateNeighborRadius(configuration.cutoffNM + configuration.neighborSkinNM,
                                                          cell: cell)
            }
        }
        guard blockers.isEmpty else {
            throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted())
        }
        guard startTimePS.isFinite, startTimePS >= 0 else {
            throw VivoMDRuntimeError.metal("initial MD clock is invalid")
        }

        let packed = try VivoMDSystemPacker.pack(system)
        let device = try requestedDevice ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoMDRuntimeError.metal("Apple-silicon unified memory and command queue are required")
        }
        queue.label = "NumiVivo.MD.Queue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        var required: [NumiVivoKernel] = [
            .mdClearForce, .mdClearStatus,
            .mdUpdateVirtualPosition, .mdUpdateVirtualVelocity, .mdRedistributeVirtualForce,
            .mdBonded, .mdBuildNeighborList, .mdValidateNeighborDisplacement,
            .mdNonbondedNeighbor, .mdNonbondedDirect,
            .mdHalfKick, .mdDrift, .mdLangevin,
            .mdConstraintPosition, .mdConstraintVelocity, .mdValidateConstraints,
            .mdKinetic, .mdValidate
        ]
        if configuration.electrostatics == .pme {
            required.append(contentsOf: [.mdPMERealSpaceNeighbor, .mdPMEExceptionCorrection])
        }
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in required { pipelines[kernel] = try await catalog.pipeline(kernel) }
        let initialVelocities = velocities ?? [VivoVector3D](repeating: .zero,
                                                             count: system.particles.count)
        let arena = try await VivoMDGPUArena.make(device: device, queue: queue,
                                                  packed: packed,
                                                  initial: initialState,
                                                  velocities: initialVelocities,
                                                  configuration: configuration)
        let neighborGrid: VivoMDNeighborGrid?
        if configuration.resolvedNeighborListEnabled, let cell = initialState.periodicCell {
            neighborGrid = try await VivoMDNeighborGrid.make(
                device: device, catalog: catalog,
                particleCount: packed.particleCount, cell: cell,
                neighborRadiusNM: configuration.cutoffNM + configuration.neighborSkinNM,
                neighborCapacity: arena.neighborCapacity
            )
        } else { neighborGrid = nil }

        let pmeEngine: VivoPMEEngine?
        if configuration.electrostatics == .pme {
            guard let cell = initialState.periodicCell else {
                throw VivoMDRuntimeError.unsupported(["PME requires a periodic cell"])
            }
            let totalCharge = system.particles.reduce(0.0) { $0 + $1.chargeE }
            pmeEngine = try await VivoPMEEngine.make(device: device, catalog: catalog,
                                                     particleCount: packed.particleCount,
                                                     totalChargeE: totalCharge,
                                                     cell: cell,
                                                     configuration: configuration)
        } else { pmeEngine = nil }

        return try .init(system: system, configuration: configuration,
                         packed: packed, queue: queue, arena: arena,
                         pipelines: pipelines, periodicCell: initialState.periodicCell,
                         neighborGrid: neighborGrid, pmeEngine: pmeEngine,
                         device: device, acceptedStep: startStep,
                         acceptedTimePS: startTimePS)
    }

    private init(system: VivoClassicalSystem,
                 configuration: VivoMDConfiguration,
                 packed: VivoMDPackedSystem,
                 queue: MTLCommandQueue,
                 arena: VivoMDGPUArena,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline],
                 periodicCell: VivoPeriodicCell?,
                 neighborGrid: VivoMDNeighborGrid?,
                 pmeEngine: VivoPMEEngine?,
                 device: MTLDevice,
                 acceptedStep: UInt64,
                 acceptedTimePS: Double) throws {
        self.system = system; self.configuration = configuration; self.packed = packed
        self.queue = queue; self.arena = arena; self.pipelines = pipelines
        self.periodicCell = periodicCell; self.neighborGrid = neighborGrid; self.pmeEngine = pmeEngine
        self.acceptedStep = acceptedStep; self.acceptedTimePS = acceptedTimePS
        systemFingerprint = packed.systemFingerprint
        configurationFingerprint = try configuration.fingerprint()
        deviceName = device.name; deviceRegistryID = device.registryID
    }

    public func step() async throws -> VivoMDStepCertificate {
        try requireIdle()
        guard acceptedStep < UInt64.max else { throw VivoMDRuntimeError.metal("MD step index overflow") }
        inFlight = true; defer { inFlight = false }
        let before = acceptedTimePS
        let abi = try VivoMDMetalABI.command(packed: packed, configuration: configuration,
                                             cell: periodicCell, stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("step command unavailable")
        }
        command.label = "NumiVivo.MD.step.\(acceptedStep)"
        let bytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0,
                  to: arena.candidatePosition, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0,
                  to: arena.candidateVelocity, destinationOffset: 0, size: bytes)
        blit.endEncoding()

        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.candidatePosition,
                                  velocity: arena.candidateVelocity, abi: abi)
        if configuration.resolvedNeighborListEnabled {
            try buildNeighborList(command, position: arena.candidatePosition, abi: abi)
        }
        try encodeForces(command, position: arena.candidatePosition, abi: abi)
        try encode(.mdHalfKick, command: command,
                   buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: abi)
        if !packed.constraints.isEmpty { try encodeVelocityConstraints(command, abi: abi) }

        if configuration.ensemble == .nvt {
            var half = abi; half.dtPS *= 0.5
            try encode(.mdDrift, command: command,
                       buffers: [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi: half)
            if !packed.constraints.isEmpty { try encodePositionConstraints(command, abi: abi) }
            try encode(.mdLangevin, command: command,
                       buffers: [arena.candidateVelocity, arena.dynamics], abi: abi)
            if !packed.constraints.isEmpty { try encodeVelocityConstraints(command, abi: abi) }
            try encode(.mdDrift, command: command,
                       buffers: [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi: half)
            if !packed.constraints.isEmpty { try encodePositionConstraints(command, abi: abi) }
        } else {
            try encode(.mdDrift, command: command,
                       buffers: [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi: abi)
            if !packed.constraints.isEmpty { try encodePositionConstraints(command, abi: abi) }
        }

        try updateVirtualPosition(command, position: arena.candidatePosition, abi: abi)
        if configuration.resolvedNeighborListEnabled {
            try encode(.mdValidateNeighborDisplacement, command: command,
                       buffers: [arena.candidatePosition, arena.neighborReferencePosition, arena.status], abi: abi)
        }
        try encodeForces(command, position: arena.candidatePosition, abi: abi)
        try encode(.mdHalfKick, command: command,
                   buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: abi)
        if !packed.constraints.isEmpty {
            try encodeVelocityConstraints(command, abi: abi)
            try encode(.mdValidateConstraints, command: command,
                       buffers: [arena.candidatePosition, arena.candidateVelocity,
                                 arena.constraints, arena.constraintOffsets,
                                 arena.constraintIncidence, arena.status], abi: abi)
        }
        try updateVirtualVelocity(command, velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdValidate, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity, arena.forceEnergy, arena.status], abi: abi)
        try await complete(command)

        let status = arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee
        let first = status.firstParticle == UInt32.max ? nil : status.firstParticle
        guard status.flags == 0 else {
            return .init(systemFingerprint: systemFingerprint,
                         configurationFingerprint: configurationFingerprint,
                         deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                         stepIndex: acceptedStep, timeBeforePS: before, timeAfterPS: before,
                         committed: false, statusFlags: status.flags,
                         firstViolationParticle: first, violationCount: status.violationCount)
        }
        arena.commit(); acceptedStep += 1; acceptedTimePS += configuration.timeStepPS
        return .init(systemFingerprint: systemFingerprint,
                     configurationFingerprint: configurationFingerprint,
                     deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                     stepIndex: acceptedStep - 1, timeBeforePS: before,
                     timeAfterPS: acceptedTimePS, committed: true,
                     statusFlags: 0, firstViolationParticle: nil, violationCount: 0)
    }

    public func snapshot() async throws -> VivoMDStateSnapshot {
        try requireIdle()
        let values = try await readAcceptedState()
        let result = VivoMDStateSnapshot(systemFingerprint: systemFingerprint,
                                         configurationFingerprint: configurationFingerprint,
                                         stepIndex: acceptedStep, timePS: acceptedTimePS,
                                         positionsNM: values.positions,
                                         velocitiesNMPerPS: values.velocities,
                                         periodicCell: periodicCell)
        try result.validate(particleCount: arena.particleCount)
        return result
    }

    public func checkpoint() async throws -> VivoMDCheckpoint {
        let state = try await snapshot()
        let value = VivoMDCheckpoint(systemFingerprint: systemFingerprint,
                                     configurationFingerprint: configurationFingerprint,
                                     acceptedStep: acceptedStep, timePS: acceptedTimePS,
                                     positionsNM: state.positionsNM,
                                     velocitiesNMPerPS: state.velocitiesNMPerPS,
                                     periodicCell: periodicCell)
        try value.validate(particleCount: arena.particleCount)
        return value
    }

    public func observables() async throws -> VivoMDObservables {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        let abi = try VivoMDMetalABI.command(packed: packed, configuration: configuration,
                                             cell: periodicCell, stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer() else {
            throw VivoMDRuntimeError.metal("observables command unavailable")
        }
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.acceptedPosition,
                                  velocity: arena.acceptedVelocity, abi: abi)
        if configuration.resolvedNeighborListEnabled {
            try buildNeighborList(command, position: arena.acceptedPosition, abi: abi)
        }
        try encodeForces(command, position: arena.acceptedPosition, abi: abi)
        try encode(.mdKinetic, command: command,
                   buffers: [arena.acceptedVelocity, arena.dynamics, arena.kinetic], abi: abi)
        guard let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("observables readback encoder unavailable")
        }
        let vectorBytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        let scalarBytes = arena.particleCount * MemoryLayout<Float>.stride
        blit.copy(from: arena.forceEnergy, sourceOffset: 0,
                  to: arena.forceEnergyReadback, destinationOffset: 0, size: vectorBytes)
        blit.copy(from: arena.kinetic, sourceOffset: 0,
                  to: arena.kineticReadback, destinationOffset: 0, size: scalarBytes)
        blit.endEncoding(); try await complete(command)
        let status = arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee
        guard status.flags == 0 else {
            throw VivoMDRuntimeError.metal("accepted-state force evaluation failed with status \(status.flags)")
        }
        let energies = arena.forceEnergyReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let kinetics = arena.kineticReadback.contents().assumingMemoryBound(to: Float.self)
        var potential = 0.0, kinetic = 0.0
        for index in 0..<arena.particleCount {
            potential += Double(energies[index].w); kinetic += Double(kinetics[index])
        }
        let massive = UInt64(system.particles.lazy.filter { $0.massDa > 0 }.count)
        let multiplied = massive.multipliedReportingOverflow(by: 3)
        guard !multiplied.overflow else { throw VivoMDRuntimeError.metal("degrees-of-freedom count overflow") }
        let constraintCount = UInt64(system.constraints.count)
        let dof = multiplied.partialValue > constraintCount ? multiplied.partialValue - constraintCount : 0
        let temperature = dof > 0 ? 2 * kinetic / (Double(dof) * 0.00831446261815324) : nil
        guard potential.isFinite, kinetic.isFinite, temperature?.isFinite != false else {
            throw VivoMDRuntimeError.metal("accepted observables are non-finite")
        }
        return .init(systemFingerprint: systemFingerprint,
                     configurationFingerprint: configurationFingerprint,
                     stepIndex: acceptedStep, timePS: acceptedTimePS,
                     potentialEnergyKJPerMol: potential, kineticEnergyKJPerMol: kinetic,
                     temperatureK: temperature, degreesOfFreedom: dof)
    }

    private func normalizeVirtualSites(_ command: MTLCommandBuffer,
                                       position: MTLBuffer, velocity: MTLBuffer,
                                       abi: VivoMDMetalCommand) throws {
        guard !packed.virtualSites.isEmpty else { return }
        try updateVirtualPosition(command, position: position, abi: abi)
        try updateVirtualVelocity(command, velocity: velocity, abi: abi)
    }

    private func updateVirtualPosition(_ command: MTLCommandBuffer,
                                       position: MTLBuffer,
                                       abi: VivoMDMetalCommand) throws {
        guard !packed.virtualSites.isEmpty else { return }
        try encode(.mdUpdateVirtualPosition, command: command,
                   buffers: [position, arena.virtualSites,
                             arena.virtualSiteIndexByParticle, arena.status], abi: abi)
    }

    private func updateVirtualVelocity(_ command: MTLCommandBuffer,
                                       velocity: MTLBuffer,
                                       abi: VivoMDMetalCommand) throws {
        guard !packed.virtualSites.isEmpty else { return }
        try encode(.mdUpdateVirtualVelocity, command: command,
                   buffers: [velocity, arena.virtualSites,
                             arena.virtualSiteIndexByParticle], abi: abi)
    }

    private func buildNeighborList(_ command: MTLCommandBuffer,
                                   position: MTLBuffer,
                                   abi: VivoMDMetalCommand) throws {
        guard abi.neighborCapacity == arena.neighborCapacity else {
            throw VivoMDRuntimeError.metal("neighbor-list ABI capacity differs from arena allocation")
        }
        if let neighborGrid {
            guard neighborGrid.plan.neighborCapacity == arena.neighborCapacity else {
                throw VivoMDRuntimeError.metal("spatial-grid neighbor capacity differs from arena allocation")
            }
            try neighborGrid.encode(commandBuffer: command, positions: position,
                                    neighborCounts: arena.neighborCounts,
                                    neighborIndices: arena.neighborIndices,
                                    referencePositions: arena.neighborReferencePosition,
                                    status: arena.status)
        } else {
            try encode(.mdBuildNeighborList, command: command,
                       buffers: [position, arena.neighborCounts, arena.neighborIndices,
                                 arena.neighborReferencePosition, arena.status], abi: abi)
        }
    }

    private func encodePositionConstraints(_ command: MTLCommandBuffer,
                                           abi: VivoMDMetalCommand) throws {
        var source = arena.candidatePosition, destination = arena.positionScratch
        var sourceIsCandidate = true
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintPosition, command: command,
                       buffers: [source, destination, arena.dynamics,
                                 arena.constraints, arena.constraintOffsets,
                                 arena.constraintIncidence, arena.status], abi: abi)
            swap(&source, &destination); sourceIsCandidate.toggle()
        }
        if !sourceIsCandidate { try copy(command, from: source, to: arena.candidatePosition) }
    }

    private func encodeVelocityConstraints(_ command: MTLCommandBuffer,
                                           abi: VivoMDMetalCommand) throws {
        var source = arena.candidateVelocity, destination = arena.velocityScratch
        var sourceIsCandidate = true
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintVelocity, command: command,
                       buffers: [arena.candidatePosition, source, destination,
                                 arena.dynamics, arena.constraints,
                                 arena.constraintOffsets, arena.constraintIncidence,
                                 arena.status], abi: abi)
            swap(&source, &destination); sourceIsCandidate.toggle()
        }
        if !sourceIsCandidate { try copy(command, from: source, to: arena.candidateVelocity) }
    }

    private func encodeForces(_ command: MTLCommandBuffer,
                              position: MTLBuffer,
                              abi: VivoMDMetalCommand) throws {
        try encode(.mdClearForce, command: command, buffers: [arena.forceEnergy], abi: abi)
        try encode(.mdBonded, command: command,
                   buffers: [position, arena.forceEnergy,
                             arena.bonds, arena.bondOffsets, arena.bondIncidence,
                             arena.angles, arena.angleOffsets, arena.angleIncidence,
                             arena.torsions, arena.torsionOffsets,
                             arena.torsionIncidence, arena.status], abi: abi)

        if configuration.electrostatics == .pme {
            guard let pmeEngine else { throw VivoMDRuntimeError.metal("PME engine is absent") }
            try encode(.mdPMERealSpaceNeighbor, command: command,
                       buffers: [position, arena.forceEnergy, arena.dynamics,
                                 arena.typeIndices, arena.pairMatrix,
                                 arena.exceptions, arena.exceptionOffsets,
                                 arena.exceptionPartners, arena.exceptionIndices,
                                 arena.neighborCounts, arena.neighborIndices,
                                 arena.status], abi: abi)
            try pmeEngine.encodeReciprocal(commandBuffer: command,
                                           positions: position,
                                           dynamics: arena.dynamics,
                                           forceEnergy: arena.forceEnergy,
                                           status: arena.status)
            if !packed.pairExceptions.isEmpty {
                try encode(.mdPMEExceptionCorrection, command: command,
                           buffers: [position, arena.dynamics, arena.forceEnergy,
                                     arena.exceptions, arena.exceptionOffsets,
                                     arena.exceptionPartners, arena.exceptionIndices,
                                     arena.status], abi: abi)
            }
        } else if configuration.resolvedNeighborListEnabled {
            try encode(.mdNonbondedNeighbor, command: command,
                       buffers: [position, arena.forceEnergy, arena.dynamics,
                                 arena.typeIndices, arena.pairMatrix,
                                 arena.exceptions, arena.exceptionOffsets,
                                 arena.exceptionPartners, arena.exceptionIndices,
                                 arena.neighborCounts, arena.neighborIndices,
                                 arena.status], abi: abi)
        } else {
            try encode(.mdNonbondedDirect, command: command,
                       buffers: [position, arena.forceEnergy, arena.dynamics,
                                 arena.typeIndices, arena.pairMatrix,
                                 arena.exceptions, arena.exceptionOffsets,
                                 arena.exceptionPartners, arena.exceptionIndices,
                                 arena.status], abi: abi)
        }
        if !packed.virtualSites.isEmpty {
            try encode(.mdRedistributeVirtualForce, command: command,
                       buffers: [arena.forceEnergy, arena.virtualSites,
                                 arena.virtualParentOffsets,
                                 arena.virtualParentIncidence, arena.status], abi: abi)
        }
    }

    private func readAcceptedState() async throws
    -> (positions: [VivoVector3D], velocities: [VivoVector3D]) {
        inFlight = true; defer { inFlight = false }
        guard let command = queue.makeCommandBuffer() else {
            throw VivoMDRuntimeError.metal("state readback command unavailable")
        }
        let abi = try VivoMDMetalABI.command(packed: packed, configuration: configuration,
                                             cell: periodicCell, stepIndex: acceptedStep)
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.acceptedPosition,
                                  velocity: arena.acceptedVelocity, abi: abi)
        guard let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("state readback encoder unavailable")
        }
        let bytes = arena.particleCount * MemoryLayout<SIMD4<Float>>.stride
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0,
                  to: arena.positionReadback, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0,
                  to: arena.velocityReadback, destinationOffset: 0, size: bytes)
        blit.endEncoding(); try await complete(command)
        let status = arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee
        guard status.flags == 0 else {
            throw VivoMDRuntimeError.metal("derived-site normalization failed with status \(status.flags)")
        }
        let p = arena.positionReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let v = arena.velocityReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
        positions.reserveCapacity(arena.particleCount); velocities.reserveCapacity(arena.particleCount)
        for index in 0..<arena.particleCount {
            positions.append(.init(Double(p[index].x), Double(p[index].y), Double(p[index].z)))
            velocities.append(.init(Double(v[index].x), Double(v[index].y), Double(v[index].z)))
        }
        return (positions, velocities)
    }

    private func encodeStatusClear(_ command: MTLCommandBuffer) throws {
        guard let pipeline = pipelines[.mdClearStatus],
              let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("status-clear encoder unavailable")
        }
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(arena.status, offset: 0, index: 0)
        encoder.dispatchThreads(.init(width: 1, height: 1, depth: 1),
                                threadsPerThreadgroup: .init(width: 1, height: 1, depth: 1))
        encoder.endEncoding()
    }

    private func encode(_ kernel: NumiVivoKernel,
                        command: MTLCommandBuffer,
                        buffers: [MTLBuffer],
                        abi: VivoMDMetalCommand) throws {
        guard let pipeline = pipelines[kernel],
              let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("encoder unavailable for \(kernel.rawValue)")
        }
        encoder.label = kernel.rawValue; encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        var value = abi
        encoder.setBytes(&value, length: MemoryLayout<VivoMDMetalCommand>.stride,
                         index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
        encoder.endEncoding()
    }

    private func copy(_ command: MTLCommandBuffer, from: MTLBuffer, to: MTLBuffer) throws {
        guard let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("copy encoder unavailable")
        }
        blit.copy(from: from, sourceOffset: 0, to: to, destinationOffset: 0,
                  size: arena.particleCount * MemoryLayout<SIMD4<Float>>.stride)
        blit.endEncoding()
    }

    private func requireIdle() throws {
        guard !inFlight else { throw VivoMDRuntimeError.metal("MD operation already in flight") }
    }

    private func complete(_ command: MTLCommandBuffer) async throws {
        try await withCheckedThrowingContinuation {
            (continuation: CheckedContinuation<Void, Error>) in
            command.addCompletedHandler { value in
                if let error = value.error {
                    continuation.resume(throwing: VivoMDRuntimeError.metal(String(describing: error)))
                } else { continuation.resume(returning: ()) }
            }
            command.commit()
        }
    }
}
