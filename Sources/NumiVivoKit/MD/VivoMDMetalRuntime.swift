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
    public let barostat: VivoMDBarostatMoveCertificate?
}

private struct VivoMDMinPositionCommand {
    var particleCount: UInt32
    var stepScale: Float
    var maxDisplacementNM: Float
    var reserved: Float = 0
}
private struct VivoMDMinReduceCommand {
    var sourceCount: UInt32
    var reserved0: UInt32 = 0
    var reserved1: UInt32 = 0
    var reserved2: UInt32 = 0
}
private struct VivoMDReductionWorkspace {
    let termsA: MTLBuffer
    let termsB: MTLBuffer
    let scalar: MTLBuffer
}

public actor VivoMDMetalRuntime {
    public nonisolated let system: VivoClassicalSystem
    public nonisolated let configuration: VivoMDConfiguration
    public nonisolated let systemFingerprint: VivoFingerprint
    public nonisolated let configurationFingerprint: VivoFingerprint
    public nonisolated let deviceName: String
    public nonisolated let deviceRegistryID: UInt64

    private let device: MTLDevice
    private let catalog: NumiVivoPipelineCatalog
    private let queue: MTLCommandQueue
    private let packed: VivoMDPackedSystem
    private let arena: VivoMDGPUArena
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let totalChargeE: Double
    private let barostatEngine: VivoMDBarostatEngine?
    private var periodicCell: VivoPeriodicCell?
    private var neighborGrid: VivoMDNeighborGrid?
    private var pmeEngine: VivoPMEEngine?
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
        let sf = try system.fingerprint(), cf = try configuration.fingerprint()
        try checkpoint.validate(particleCount: system.particles.count)
        guard checkpoint.systemFingerprint == sf, checkpoint.configurationFingerprint == cf else {
            throw VivoArtifactValidationError.incompatible("MD checkpoint system or configuration identity differs from requested runtime")
        }
        let initial = VivoClassicalInitialState(systemFingerprint: sf,
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
        if configuration.thermostat == .velocityRescale {
            blockers.append("stochastic velocity-rescale thermostat is not installed")
        }
        if configuration.ensemble != .nve, configuration.thermostat != .langevinMiddle {
            blockers.append("current thermostatted NVT/NPT execution requires langevinMiddle")
        }
        if configuration.ensemble == .nve, configuration.thermostat != .none {
            blockers.append("NVE cannot use a thermostat")
        }
        if let cell = initialState.periodicCell {
            try VivoMDMetalABI.validateMinimumImageCutoff(configuration.cutoffNM, cell: cell)
            if configuration.resolvedNeighborListEnabled {
                try VivoMDMetalABI.validateNeighborRadius(configuration.cutoffNM + configuration.neighborSkinNM, cell: cell)
            }
        }
        guard blockers.isEmpty else { throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted()) }
        guard startTimePS.isFinite, startTimePS >= 0 else { throw VivoMDRuntimeError.metal("initial MD clock is invalid") }

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
            .mdMinimizeTerms, .mdMinimizeReduceStage, .mdMinimizePosition, .mdZeroVelocity,
            .mdHalfKick, .mdDrift, .mdLangevin,
            .mdConstraintPosition, .mdConstraintVelocity, .mdValidateConstraints,
            .mdKinetic, .mdValidate
        ]
        if configuration.electrostatics == .pme {
            required += [.mdPMERealSpaceNeighbor, .mdPMEExceptionCorrection]
        }
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for kernel in required { pipelines[kernel] = try await catalog.pipeline(kernel) }

        let initialVelocities = velocities ?? [VivoVector3D](repeating: .zero, count: system.particles.count)
        let arena = try await VivoMDGPUArena.make(device: device, queue: queue, packed: packed,
                                                  initial: initialState, velocities: initialVelocities,
                                                  configuration: configuration)
        let totalCharge = system.particles.reduce(0.0) { $0 + $1.chargeE }
        let grid = try await makeNeighborGrid(device: device, catalog: catalog, packed: packed,
                                              arena: arena, configuration: configuration,
                                              cell: initialState.periodicCell)
        let pme = try await makePME(device: device, catalog: catalog, packed: packed,
                                    totalChargeE: totalCharge, configuration: configuration,
                                    cell: initialState.periodicCell)
        let barostat: VivoMDBarostatEngine? = configuration.ensemble == .npt
            ? try await VivoMDBarostatEngine.make(device: device, catalog: catalog, system: system)
            : nil
        return try .init(system: system, configuration: configuration,
                         packed: packed, device: device, catalog: catalog,
                         queue: queue, arena: arena, pipelines: pipelines,
                         periodicCell: initialState.periodicCell,
                         neighborGrid: grid, pmeEngine: pme,
                         totalChargeE: totalCharge, barostatEngine: barostat,
                         acceptedStep: startStep, acceptedTimePS: startTimePS)
    }

    private static func makeNeighborGrid(device: MTLDevice,
                                         catalog: NumiVivoPipelineCatalog,
                                         packed: VivoMDPackedSystem,
                                         arena: VivoMDGPUArena,
                                         configuration: VivoMDConfiguration,
                                         cell: VivoPeriodicCell?) async throws -> VivoMDNeighborGrid? {
        guard configuration.resolvedNeighborListEnabled, let cell else { return nil }
        return try await VivoMDNeighborGrid.make(device: device, catalog: catalog,
                                                  particleCount: packed.particleCount,
                                                  cell: cell,
                                                  neighborRadiusNM: configuration.cutoffNM + configuration.neighborSkinNM,
                                                  neighborCapacity: arena.neighborCapacity)
    }

    private static func makePME(device: MTLDevice,
                                catalog: NumiVivoPipelineCatalog,
                                packed: VivoMDPackedSystem,
                                totalChargeE: Double,
                                configuration: VivoMDConfiguration,
                                cell: VivoPeriodicCell?) async throws -> VivoPMEEngine? {
        guard configuration.electrostatics == .pme else { return nil }
        guard let cell else { throw VivoMDRuntimeError.unsupported(["PME requires a periodic cell"]) }
        return try await VivoPMEEngine.make(device: device, catalog: catalog,
                                            particleCount: packed.particleCount,
                                            totalChargeE: totalChargeE,
                                            cell: cell, configuration: configuration)
    }

    private init(system: VivoClassicalSystem,
                 configuration: VivoMDConfiguration,
                 packed: VivoMDPackedSystem,
                 device: MTLDevice,
                 catalog: NumiVivoPipelineCatalog,
                 queue: MTLCommandQueue,
                 arena: VivoMDGPUArena,
                 pipelines: [NumiVivoKernel: NumiVivoPipeline],
                 periodicCell: VivoPeriodicCell?,
                 neighborGrid: VivoMDNeighborGrid?,
                 pmeEngine: VivoPMEEngine?,
                 totalChargeE: Double,
                 barostatEngine: VivoMDBarostatEngine?,
                 acceptedStep: UInt64,
                 acceptedTimePS: Double) throws {
        self.system = system; self.configuration = configuration; self.packed = packed
        self.device = device; self.catalog = catalog; self.queue = queue; self.arena = arena
        self.pipelines = pipelines; self.periodicCell = periodicCell
        self.neighborGrid = neighborGrid; self.pmeEngine = pmeEngine
        self.totalChargeE = totalChargeE; self.barostatEngine = barostatEngine
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
        let stepIndex = acceptedStep
        let abi = try currentCommand(stepIndex: stepIndex)
        guard let command = queue.makeCommandBuffer(), let blit = command.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("step command unavailable")
        }
        command.label = "NumiVivo.MD.step.\(stepIndex)"
        let bytes = vectorBytes
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0, to: arena.candidatePosition, destinationOffset: 0, size: bytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0, to: arena.candidateVelocity, destinationOffset: 0, size: bytes)
        blit.endEncoding()
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        if configuration.resolvedNeighborListEnabled { try buildNeighborList(command, position: arena.candidatePosition, abi: abi) }
        try encodeForces(command, position: arena.candidatePosition, abi: abi)
        try encode(.mdHalfKick, command: command,
                   buffers: [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi: abi)
        if !packed.constraints.isEmpty { try encodeVelocityConstraints(command, abi: abi) }

        if configuration.ensemble != .nve {
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
        let status = readStatus()
        let first = status.firstParticle == UInt32.max ? nil : status.firstParticle
        guard status.flags == 0 else {
            return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                         deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                         stepIndex: stepIndex, timeBeforePS: before, timeAfterPS: before,
                         committed: false, statusFlags: status.flags,
                         firstViolationParticle: first, violationCount: status.violationCount,
                         barostat: nil)
        }
        arena.commit(); acceptedStep += 1; acceptedTimePS += configuration.timeStepPS
        let barostatMove: VivoMDBarostatMoveCertificate?
        if configuration.ensemble == .npt, acceptedStep % UInt64(configuration.barostatInterval) == 0 {
            barostatMove = try await attemptBarostatMove()
        } else { barostatMove = nil }
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                     deviceName: deviceName, deviceRegistryID: deviceRegistryID,
                     stepIndex: stepIndex, timeBeforePS: before, timeAfterPS: acceptedTimePS,
                     committed: true, statusFlags: 0, firstViolationParticle: nil,
                     violationCount: 0, barostat: barostatMove)
    }

    public func minimize(_ minimization: VivoMDMinimizationConfiguration = .init()) async throws -> VivoMDMinimizationCertificate {
        try minimization.validate(); try requireIdle(); inFlight = true; defer { inFlight = false }
        let workspace = try makeReductionWorkspace(label: "min")
        var stepScale = minimization.initialStepScale
        var attempted: UInt32 = 0, accepted: UInt32 = 0, rejected: UInt32 = 0
        var baseline = try await evaluateState(position: arena.acceptedPosition,
                                               velocity: arena.acceptedVelocity,
                                               workspace: workspace)
        let initialEnergy = baseline.energy
        while attempted < minimization.maximumIterations,
              baseline.maximumForce > minimization.forceToleranceKJPerMolNM,
              stepScale >= minimization.minimumStepScale {
            attempted += 1
            let abi = try currentCommand(stepIndex: acceptedStep)
            guard let proposal = queue.makeCommandBuffer(), let blit = proposal.makeBlitCommandEncoder() else {
                throw VivoMDRuntimeError.metal("minimization proposal command unavailable")
            }
            blit.copy(from: arena.acceptedPosition, sourceOffset: 0,
                      to: arena.candidatePosition, destinationOffset: 0, size: vectorBytes)
            blit.endEncoding()
            try encodeStatusClear(proposal)
            try encode(.mdZeroVelocity, command: proposal, buffers: [arena.candidateVelocity], abi: abi)
            try encodeMinimizationPosition(proposal, abi: abi,
                                           stepScale: stepScale,
                                           maximumDisplacementNM: minimization.maximumDisplacementNM)
            if !packed.constraints.isEmpty { try encodePositionConstraints(proposal, abi: abi) }
            try updateVirtualPosition(proposal, position: arena.candidatePosition, abi: abi)
            if !packed.constraints.isEmpty {
                try encode(.mdValidateConstraints, command: proposal,
                           buffers: [arena.candidatePosition, arena.candidateVelocity,
                                     arena.constraints, arena.constraintOffsets,
                                     arena.constraintIncidence, arena.status], abi: abi)
            }
            try await complete(proposal)
            if readStatus().flags != 0 {
                rejected += 1; stepScale *= minimization.rejectedStepShrink
                baseline = try await evaluateState(position: arena.acceptedPosition,
                                                   velocity: arena.acceptedVelocity,
                                                   workspace: workspace)
                continue
            }
            let candidate = try await evaluateState(position: arena.candidatePosition,
                                                    velocity: arena.candidateVelocity,
                                                    workspace: workspace)
            if candidate.energy < baseline.energy {
                arena.commit(); accepted += 1
                stepScale = min(minimization.maximumStepScale, stepScale * minimization.acceptedStepGrowth)
                baseline = candidate
            } else {
                rejected += 1; stepScale *= minimization.rejectedStepShrink
                baseline = try await evaluateState(position: arena.acceptedPosition,
                                                   velocity: arena.acceptedVelocity,
                                                   workspace: workspace)
            }
        }
        let final = try await evaluateState(position: arena.acceptedPosition,
                                            velocity: arena.acceptedVelocity,
                                            workspace: workspace)
        let abi = try currentCommand(stepIndex: acceptedStep)
        guard let zero = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("minimization velocity-reset command unavailable") }
        try encode(.mdZeroVelocity, command: zero, buffers: [arena.acceptedVelocity], abi: abi)
        try updateVirtualVelocity(zero, velocity: arena.acceptedVelocity, abi: abi)
        try await complete(zero)
        return .init(systemFingerprint: systemFingerprint,
                     configurationFingerprint: configurationFingerprint,
                     converged: final.maximumForce <= minimization.forceToleranceKJPerMolNM,
                     attemptedIterations: attempted, acceptedIterations: accepted,
                     rejectedIterations: rejected,
                     initialPotentialEnergyKJPerMol: initialEnergy,
                     finalPotentialEnergyKJPerMol: final.energy,
                     finalMaximumForceKJPerMolNM: final.maximumForce,
                     finalStepScale: stepScale)
    }

    public func snapshot() async throws -> VivoMDStateSnapshot {
        try requireIdle(); let values = try await readAcceptedState()
        let result = VivoMDStateSnapshot(systemFingerprint: systemFingerprint,
                                         configurationFingerprint: configurationFingerprint,
                                         stepIndex: acceptedStep, timePS: acceptedTimePS,
                                         positionsNM: values.positions,
                                         velocitiesNMPerPS: values.velocities,
                                         periodicCell: periodicCell)
        try result.validate(particleCount: arena.particleCount); return result
    }

    public func checkpoint() async throws -> VivoMDCheckpoint {
        let state = try await snapshot()
        let value = VivoMDCheckpoint(systemFingerprint: systemFingerprint,
                                     configurationFingerprint: configurationFingerprint,
                                     acceptedStep: acceptedStep, timePS: acceptedTimePS,
                                     positionsNM: state.positionsNM,
                                     velocitiesNMPerPS: state.velocitiesNMPerPS,
                                     periodicCell: periodicCell)
        try value.validate(particleCount: arena.particleCount); return value
    }

    public func observables() async throws -> VivoMDObservables {
        try requireIdle(); inFlight = true; defer { inFlight = false }
        let workspace = try makeReductionWorkspace(label: "obs")
        let potential = try await evaluateState(position: arena.acceptedPosition,
                                                velocity: arena.acceptedVelocity,
                                                workspace: workspace).energy
        let abi = try currentCommand(stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("kinetic command unavailable") }
        try encode(.mdKinetic, command: command,
                   buffers: [arena.acceptedVelocity, arena.dynamics, arena.kinetic], abi: abi)
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("kinetic readback unavailable") }
        let scalarBytes = arena.particleCount * MemoryLayout<Float>.stride
        blit.copy(from: arena.kinetic, sourceOffset: 0,
                  to: arena.kineticReadback, destinationOffset: 0, size: scalarBytes)
        blit.endEncoding(); try await complete(command)
        let kinetics = arena.kineticReadback.contents().assumingMemoryBound(to: Float.self)
        var kinetic = 0.0
        for index in 0..<arena.particleCount { kinetic += Double(kinetics[index]) }
        let massive = UInt64(system.particles.lazy.filter { $0.massDa > 0 }.count)
        let product = massive.multipliedReportingOverflow(by: 3)
        guard !product.overflow else { throw VivoMDRuntimeError.metal("degrees-of-freedom count overflow") }
        let constraints = UInt64(system.constraints.count)
        let dof = product.partialValue > constraints ? product.partialValue - constraints : 0
        let temperature = dof > 0 ? 2 * kinetic / (Double(dof) * 0.00831446261815324) : nil
        guard potential.isFinite, kinetic.isFinite, temperature?.isFinite != false else {
            throw VivoMDRuntimeError.metal("accepted observables are non-finite")
        }
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
                     stepIndex: acceptedStep, timePS: acceptedTimePS,
                     potentialEnergyKJPerMol: potential, kineticEnergyKJPerMol: kinetic,
                     temperatureK: temperature, degreesOfFreedom: dof)
    }

    private func attemptBarostatMove() async throws -> VivoMDBarostatMoveCertificate {
        guard let engine = barostatEngine, let oldCell = periodicCell,
              let temperature = configuration.targetTemperatureK,
              let pressureBar = configuration.targetPressureBar else {
            throw VivoMDRuntimeError.metal("NPT barostat state is incomplete")
        }
        let oldVolume = cellVolume(oldCell)
        let uProposal = deterministicUniform(seed: configuration.randomSeed,
                                             step: acceptedStep, stream: 0x4e505456)
        let uAccept = deterministicUniform(seed: configuration.randomSeed,
                                           step: acceptedStep, stream: 0x4e505441)
        let dx = (2 * uProposal - 1) * configuration.resolvedBarostatMaximumLogVolumeStep
        let volumeRatio = exp(dx)
        let scale = exp(dx / 3)
        let proposedCell = VivoPeriodicCell(a: oldCell.a * scale,
                                            b: oldCell.b * scale,
                                            c: oldCell.c * scale)
        let proposedVolume = oldVolume * volumeRatio
        guard proposedCell.isValid, proposedVolume.isFinite, proposedVolume > 0 else {
            return rejectedBarostat(oldVolume: oldVolume, dx: dx)
        }
        do {
            try VivoMDMetalABI.validateMinimumImageCutoff(configuration.cutoffNM, cell: proposedCell)
            if configuration.resolvedNeighborListEnabled {
                try VivoMDMetalABI.validateNeighborRadius(configuration.cutoffNM + configuration.neighborSkinNM,
                                                          cell: proposedCell)
            }
        } catch {
            return rejectedBarostat(oldVolume: oldVolume, dx: dx)
        }

        let workspace = try makeReductionWorkspace(label: "npt")
        let oldState = try await evaluateState(position: arena.acceptedPosition,
                                               velocity: arena.acceptedVelocity,
                                               workspace: workspace)
        guard let copy = queue.makeCommandBuffer(), let blit = copy.makeBlitCommandEncoder() else {
            throw VivoMDRuntimeError.metal("barostat copy command unavailable")
        }
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0,
                  to: arena.candidatePosition, destinationOffset: 0, size: vectorBytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0,
                  to: arena.candidateVelocity, destinationOffset: 0, size: vectorBytes)
        blit.endEncoding()
        try engine.encodeProposal(commandBuffer: copy,
                                  positions: arena.candidatePosition,
                                  dynamics: arena.dynamics,
                                  cell: oldCell, scale: scale)
        try await complete(copy)

        let oldGrid = neighborGrid, oldPME = pmeEngine
        let proposedGrid: VivoMDNeighborGrid?
        let proposedPME: VivoPMEEngine?
        do {
            proposedGrid = try await Self.makeNeighborGrid(device: device, catalog: catalog,
                                                           packed: packed, arena: arena,
                                                           configuration: configuration,
                                                           cell: proposedCell)
            proposedPME = try await Self.makePME(device: device, catalog: catalog,
                                                 packed: packed, totalChargeE: totalChargeE,
                                                 configuration: configuration,
                                                 cell: proposedCell)
        } catch {
            return rejectedBarostat(oldVolume: oldVolume, dx: dx)
        }
        periodicCell = proposedCell; neighborGrid = proposedGrid; pmeEngine = proposedPME
        do {
            try await validateCandidateConstraints()
            let proposedState = try await evaluateState(position: arena.candidatePosition,
                                                        velocity: arena.candidateVelocity,
                                                        workspace: workspace)
            let pressureWork = pressureBar * (proposedVolume - oldVolume) * 0.0602214076
            let beta = 1 / (0.00831446261815324 * temperature)
            let moleculeCount = Double(engine.plan.componentCount)
            let logAcceptance = -beta * (proposedState.energy - oldState.energy + pressureWork)
                + (moleculeCount + 1) * dx
            let accepted = log(uAccept) < min(0, logAcceptance)
            if accepted {
                arena.commit()
                return .init(stepIndex: acceptedStep, attempted: true, accepted: true,
                             volumeBeforeNM3: oldVolume, volumeAfterNM3: proposedVolume,
                             logVolumeDelta: dx, logAcceptanceRatio: logAcceptance)
            }
            periodicCell = oldCell; neighborGrid = oldGrid; pmeEngine = oldPME
            return .init(stepIndex: acceptedStep, attempted: true, accepted: false,
                         volumeBeforeNM3: oldVolume, volumeAfterNM3: oldVolume,
                         logVolumeDelta: dx, logAcceptanceRatio: logAcceptance)
        } catch {
            periodicCell = oldCell; neighborGrid = oldGrid; pmeEngine = oldPME
            return rejectedBarostat(oldVolume: oldVolume, dx: dx)
        }
    }

    private func rejectedBarostat(oldVolume: Double, dx: Double) -> VivoMDBarostatMoveCertificate {
        .init(stepIndex: acceptedStep, attempted: true, accepted: false,
              volumeBeforeNM3: oldVolume, volumeAfterNM3: oldVolume,
              logVolumeDelta: dx, logAcceptanceRatio: -Double.greatestFiniteMagnitude)
    }

    private func validateCandidateConstraints() async throws {
        guard !packed.constraints.isEmpty else { return }
        let abi = try currentCommand(stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("barostat validation command unavailable") }
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.candidatePosition,
                                  velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdValidateConstraints, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity,
                             arena.constraints, arena.constraintOffsets,
                             arena.constraintIncidence, arena.status], abi: abi)
        try encode(.mdValidate, command: command,
                   buffers: [arena.candidatePosition, arena.candidateVelocity,
                             arena.forceEnergy, arena.status], abi: abi)
        try await complete(command)
        guard readStatus().flags == 0 else { throw VivoMDRuntimeError.metal("barostat candidate violates constraints") }
    }

    private func evaluateState(position: MTLBuffer,
                               velocity: MTLBuffer,
                               workspace: VivoMDReductionWorkspace) async throws
    -> (energy: Double, maximumForce: Double) {
        let abi = try currentCommand(stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("state evaluation command unavailable") }
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: position, velocity: velocity, abi: abi)
        if configuration.resolvedNeighborListEnabled { try buildNeighborList(command, position: position, abi: abi) }
        try encodeForces(command, position: position, abi: abi)
        try encode(.mdValidate, command: command,
                   buffers: [position, velocity, arena.forceEnergy, arena.status], abi: abi)
        try encode(.mdMinimizeTerms, command: command,
                   buffers: [arena.forceEnergy, workspace.termsA], abi: abi)
        var current = workspace.termsA, scratch = workspace.termsB
        var count = UInt32(arena.particleCount)
        while count > 1 {
            let output = (count + 1) / 2
            try encodeMinimizationReduce(command, source: current,
                                         destination: scratch,
                                         sourceCount: count, outputCount: output)
            swap(&current, &scratch); count = output
        }
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("state reduction readback unavailable") }
        blit.copy(from: current, sourceOffset: 0,
                  to: workspace.scalar, destinationOffset: 0,
                  size: MemoryLayout<SIMD2<Float>>.stride)
        blit.endEncoding(); try await complete(command)
        guard readStatus().flags == 0 else { throw VivoMDRuntimeError.metal("state force evaluation rejected") }
        let value = workspace.scalar.contents().assumingMemoryBound(to: SIMD2<Float>.self).pointee
        guard value.x.isFinite, value.y.isFinite, value.y >= 0 else { throw VivoMDRuntimeError.metal("state reduction is non-finite") }
        return (Double(value.x), Double(value.y))
    }

    private func makeReductionWorkspace(label: String) throws -> VivoMDReductionWorkspace {
        let bytes = max(arena.particleCount * MemoryLayout<SIMD2<Float>>.stride, 16)
        guard let a = device.makeBuffer(length: bytes, options: .storageModePrivate),
              let b = device.makeBuffer(length: bytes, options: .storageModePrivate),
              let scalar = device.makeBuffer(length: 16, options: .storageModeShared) else {
            throw VivoMDRuntimeError.metal("\(label) reduction allocation failed")
        }
        a.label = "NumiVivo.MD.\(label).reduceA"; b.label = "NumiVivo.MD.\(label).reduceB"
        scalar.label = "NumiVivo.MD.\(label).scalar"
        return .init(termsA: a, termsB: b, scalar: scalar)
    }

    private func encodeMinimizationPosition(_ command: MTLCommandBuffer,
                                            abi: VivoMDMetalCommand,
                                            stepScale: Double,
                                            maximumDisplacementNM: Double) throws {
        guard let pipeline = pipelines[.mdMinimizePosition], let encoder = command.makeComputeCommandEncoder(),
              stepScale <= Double(Float.greatestFiniteMagnitude),
              maximumDisplacementNM <= Double(Float.greatestFiniteMagnitude) else {
            throw VivoMDRuntimeError.metal("minimization position encoder or FP32 parameters unavailable")
        }
        var md = abi
        var min = VivoMDMinPositionCommand(particleCount: packed.particleCount,
                                           stepScale: Float(stepScale),
                                           maxDisplacementNM: Float(maximumDisplacementNM))
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(arena.candidatePosition, offset: 0, index: 0)
        encoder.setBuffer(arena.forceEnergy, offset: 0, index: 1)
        encoder.setBuffer(arena.dynamics, offset: 0, index: 2)
        encoder.setBytes(&md, length: MemoryLayout<VivoMDMetalCommand>.stride, index: 3)
        encoder.setBytes(&min, length: MemoryLayout<VivoMDMinPositionCommand>.stride, index: 4)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
        encoder.endEncoding()
    }

    private func encodeMinimizationReduce(_ command: MTLCommandBuffer,
                                          source: MTLBuffer,
                                          destination: MTLBuffer,
                                          sourceCount: UInt32,
                                          outputCount: UInt32) throws {
        guard let pipeline = pipelines[.mdMinimizeReduceStage], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("minimization reduction encoder unavailable")
        }
        var reduce = VivoMDMinReduceCommand(sourceCount: sourceCount)
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(source, offset: 0, index: 0); encoder.setBuffer(destination, offset: 0, index: 1)
        encoder.setBytes(&reduce, length: MemoryLayout<VivoMDMinReduceCommand>.stride, index: 2)
        encoder.dispatchThreads(pipeline.gridSize(for: Int(outputCount)),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: Int(outputCount)))
        encoder.endEncoding()
    }

    private func normalizeVirtualSites(_ command: MTLCommandBuffer,
                                       position: MTLBuffer,
                                       velocity: MTLBuffer,
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
                   buffers: [position, arena.virtualSites, arena.virtualSiteIndexByParticle, arena.status], abi: abi)
    }
    private func updateVirtualVelocity(_ command: MTLCommandBuffer,
                                       velocity: MTLBuffer,
                                       abi: VivoMDMetalCommand) throws {
        guard !packed.virtualSites.isEmpty else { return }
        try encode(.mdUpdateVirtualVelocity, command: command,
                   buffers: [velocity, arena.virtualSites, arena.virtualSiteIndexByParticle], abi: abi)
    }

    private func buildNeighborList(_ command: MTLCommandBuffer,
                                   position: MTLBuffer,
                                   abi: VivoMDMetalCommand) throws {
        guard abi.neighborCapacity == arena.neighborCapacity else { throw VivoMDRuntimeError.metal("neighbor-list ABI capacity differs from arena allocation") }
        if let neighborGrid {
            guard neighborGrid.plan.neighborCapacity == arena.neighborCapacity else { throw VivoMDRuntimeError.metal("spatial-grid neighbor capacity differs from arena allocation") }
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

    private func readAcceptedState() async throws -> (positions: [VivoVector3D], velocities: [VivoVector3D]) {
        inFlight = true; defer { inFlight = false }
        let abi = try currentCommand(stepIndex: acceptedStep)
        guard let command = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("state readback command unavailable") }
        try encodeStatusClear(command)
        try normalizeVirtualSites(command, position: arena.acceptedPosition,
                                  velocity: arena.acceptedVelocity, abi: abi)
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("state readback encoder unavailable") }
        blit.copy(from: arena.acceptedPosition, sourceOffset: 0,
                  to: arena.positionReadback, destinationOffset: 0, size: vectorBytes)
        blit.copy(from: arena.acceptedVelocity, sourceOffset: 0,
                  to: arena.velocityReadback, destinationOffset: 0, size: vectorBytes)
        blit.endEncoding(); try await complete(command)
        guard readStatus().flags == 0 else { throw VivoMDRuntimeError.metal("derived-site normalization failed") }
        let p = arena.positionReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let v = arena.velocityReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
        positions.reserveCapacity(arena.particleCount); velocities.reserveCapacity(arena.particleCount)
        for i in 0..<arena.particleCount {
            positions.append(.init(Double(p[i].x), Double(p[i].y), Double(p[i].z)))
            velocities.append(.init(Double(v[i].x), Double(v[i].y), Double(v[i].z)))
        }
        return (positions, velocities)
    }

    private func currentCommand(stepIndex: UInt64) throws -> VivoMDMetalCommand {
        try VivoMDMetalABI.command(packed: packed, configuration: configuration,
                                   cell: periodicCell, stepIndex: stepIndex)
    }
    private var vectorBytes: Int { arena.particleCount * MemoryLayout<SIMD4<Float>>.stride }
    private func readStatus() -> VivoMDMetalStatus {
        arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee
    }
    private func cellVolume(_ cell: VivoPeriodicCell) -> Double {
        abs(cell.a.dot(cell.b.cross(cell.c)))
    }
    private func deterministicUniform(seed: UInt64, step: UInt64, stream: UInt64) -> Double {
        var z = seed ^ (step &* 0x9E3779B97F4A7C15) ^ stream
        z &+= 0x9E3779B97F4A7C15
        z = (z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) &* 0x94D049BB133111EB
        z ^= z >> 31
        return (Double(z >> 11) + 0.5) * 0x1.0p-53
    }

    private func encodeStatusClear(_ command: MTLCommandBuffer) throws {
        guard let pipeline = pipelines[.mdClearStatus], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("status-clear encoder unavailable")
        }
        encoder.setComputePipelineState(pipeline.state); encoder.setBuffer(arena.status, offset: 0, index: 0)
        encoder.dispatchThreads(.init(width: 1, height: 1, depth: 1),
                                threadsPerThreadgroup: .init(width: 1, height: 1, depth: 1))
        encoder.endEncoding()
    }
    private func encode(_ kernel: NumiVivoKernel,
                        command: MTLCommandBuffer,
                        buffers: [MTLBuffer],
                        abi: VivoMDMetalCommand) throws {
        guard let pipeline = pipelines[kernel], let encoder = command.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("encoder unavailable for \(kernel.rawValue)")
        }
        encoder.label = kernel.rawValue; encoder.setComputePipelineState(pipeline.state)
        for (index, buffer) in buffers.enumerated() { encoder.setBuffer(buffer, offset: 0, index: index) }
        var value = abi
        encoder.setBytes(&value, length: MemoryLayout<VivoMDMetalCommand>.stride, index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount),
                                threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
        encoder.endEncoding()
    }
    private func copy(_ command: MTLCommandBuffer, from: MTLBuffer, to: MTLBuffer) throws {
        guard let blit = command.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("copy encoder unavailable") }
        blit.copy(from: from, sourceOffset: 0, to: to, destinationOffset: 0, size: vectorBytes); blit.endEncoding()
    }
    private func requireIdle() throws {
        guard !inFlight else { throw VivoMDRuntimeError.metal("MD operation already in flight") }
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
