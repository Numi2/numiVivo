import Foundation
import Metal
import NumiVivoShaders

public enum VivoMDRuntimeError: Error, Sendable, CustomStringConvertible {
    case unsupported([String])
    case metal(String)
    case candidateRejected(UInt32)
    public var description: String {
        switch self {
        case .unsupported(let reasons): return "MD execution is unsupported: " + reasons.joined(separator: "; ")
        case .metal(let reason): return "MD Metal runtime: " + reason
        case .candidateRejected(let flags): return "MD numerical candidate rejected with status \(flags)"
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

private struct VivoMDCellResources {
    let cell: VivoPeriodicCell?
    let neighbors: VivoMDNeighborGrid?
    let pme: VivoPMEEngine?
    let command: VivoMDMetalCommand
}
private struct VivoMDMinPositionCommand {
    var particleCount: UInt32
    var stepScale: Float
    var maxDisplacementNM: Float
    var reserved: Float = 0
}

/// Positions, velocities, periodic geometry and cell-dependent plans have one
/// owner. A dynamics step publishes none of them before its last fallible work.
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
    private let work: VivoMDWorkBuffers
    private let pipelines: [NumiVivoKernel: NumiVivoPipeline]
    private let totalChargeE: Double
    private let barostatEngine: VivoMDBarostatEngine?
    private var cellResources: VivoMDCellResources
    private var acceptedStep: UInt64
    private var acceptedTimePS: Double
    private var inFlight = false
    private var deviceFailure: String?

    public static func make(system: VivoClassicalSystem, initialState: VivoClassicalInitialState,
                            configuration: VivoMDConfiguration,
                            initialVelocitiesNMPerPS: [VivoVector3D]? = nil,
                            device: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        try await makeInternal(system: system, initial: initialState, configuration: configuration,
            velocities: initialVelocitiesNMPerPS, step: 0, time: initialState.sourceTimePS ?? 0, device: device)
    }
    public static func restore(system: VivoClassicalSystem, configuration: VivoMDConfiguration,
                               checkpoint: VivoMDCheckpoint, device: MTLDevice? = nil) async throws -> VivoMDMetalRuntime {
        try checkpoint.validate(particleCount: system.particles.count)
        let systemID = try system.fingerprint()
        let configurationID = try configuration.fingerprint()
        guard checkpoint.systemFingerprint == systemID, checkpoint.configurationFingerprint == configurationID else {
            throw VivoArtifactValidationError.incompatible("MD restore system/configuration identity mismatch")
        }
        let initial = VivoClassicalInitialState(systemFingerprint: systemID, positionsNM: checkpoint.positionsNM,
                                                periodicCell: checkpoint.periodicCell, sourceTimePS: checkpoint.timePS)
        return try await makeInternal(system: system, initial: initial, configuration: configuration,
            velocities: checkpoint.velocitiesNMPerPS, step: checkpoint.acceptedStep, time: checkpoint.timePS, device: device)
    }
    private static func makeInternal(system: VivoClassicalSystem, initial: VivoClassicalInitialState,
                                     configuration: VivoMDConfiguration, velocities: [VivoVector3D]?,
                                     step: UInt64, time: Double, device requested: MTLDevice?) async throws -> VivoMDMetalRuntime {
        let report = try VivoMDCapabilityAnalyzer.analyze(system: system, initialState: initial, configuration: configuration)
        guard report.executable else { throw VivoMDRuntimeError.unsupported(report.blockers) }
        guard time.isFinite, time >= 0 else { throw VivoMDRuntimeError.metal("invalid MD clock") }
        let packed = try VivoMDSystemPacker.pack(system)
        let device = try requested ?? VivoMetalDeviceSelector.productionDevice()
        guard device.hasUnifiedMemory, let queue = device.makeCommandQueue() else {
            throw VivoMDRuntimeError.metal("Apple-silicon unified memory and a command queue are required")
        }
        queue.label = "NumiVivo.MD.Queue"
        let catalog = try NumiVivoPipelineCatalog(device: device)
        var names: [NumiVivoKernel] = [
            .mdClearForce, .mdClearStatus, .mdUpdateVirtualPosition, .mdUpdateVirtualVelocity, .mdRedistributeVirtualForce,
            .mdBonded, .mdBuildNeighborList, .mdValidateNeighborDisplacement, .mdNonbondedNeighbor, .mdNonbondedDirect,
            .mdMinimizeTerms, .mdMinimizeReduceStage, .mdMinimizePosition, .mdZeroVelocity,
            .mdHalfKick, .mdDrift, .mdLangevin, .mdConstraintPosition, .mdConstraintVelocity, .mdValidateConstraints, .mdValidate,
            .mdDriftConstraintImpulse, .mdProjectedForceSeed, .mdValidateProjectedForce, .mdObservationTerms, .mdSumPairReduce
        ]
        if configuration.electrostatics == .pme { names += [.mdPMERealSpaceNeighbor, .mdPMEExceptionCorrection] }
        var pipelines: [NumiVivoKernel: NumiVivoPipeline] = [:]
        for name in names { pipelines[name] = try await catalog.pipeline(name) }
        let arena = try await VivoMDGPUArena.make(device: device, queue: queue, packed: packed, initial: initial,
            velocities: velocities ?? [VivoVector3D](repeating: .zero, count: system.particles.count), configuration: configuration)
        let work = try VivoMDWorkBuffers(device: device, particles: system.particles)
        let totalCharge = packed.particleDynamics.reduce(0.0) { $0 + Double($1.z) }
        let cell = try await makeCellResources(device: device, catalog: catalog, packed: packed, arena: arena,
            configuration: configuration, totalCharge: totalCharge, cell: initial.periodicCell)
        let barostat: VivoMDBarostatEngine?
        if configuration.ensemble == .npt {
            barostat = try await VivoMDBarostatEngine.make(device: device, catalog: catalog, system: system)
        } else { barostat = nil }
        let runtime = try VivoMDMetalRuntime(system: system, configuration: configuration, device: device, catalog: catalog,
            queue: queue, packed: packed, arena: arena, work: work, pipelines: pipelines, totalCharge: totalCharge,
            cell: cell, barostat: barostat, step: step, time: time)
        try await runtime.establishDerivedState()
        return runtime
    }
    private static func makeCellResources(device: MTLDevice, catalog: NumiVivoPipelineCatalog,
        packed: VivoMDPackedSystem, arena: VivoMDGPUArena, configuration: VivoMDConfiguration,
        totalCharge: Double, cell: VivoPeriodicCell?) async throws -> VivoMDCellResources {
        if let cell { try VivoMDExecutionPreflight.validateCell(cell, configuration: configuration) }
        let grid: VivoMDNeighborGrid?
        if configuration.resolvedNeighborListEnabled, let cell {
            grid = try await VivoMDNeighborGrid.make(device: device, catalog: catalog,
                particleCount: packed.particleCount, cell: cell,
                neighborRadiusNM: configuration.cutoffNM + configuration.neighborSkinNM,
                neighborCapacity: arena.neighborCapacity)
        } else { grid = nil }
        let pme: VivoPMEEngine?
        if configuration.electrostatics == .pme {
            guard let cell else { throw VivoMDRuntimeError.unsupported(["PME requires a periodic cell"]) }
            pme = try await VivoPMEEngine.make(device: device, catalog: catalog, particleCount: packed.particleCount,
                totalChargeE: totalCharge, cell: cell, configuration: configuration)
        } else { pme = nil }
        let command = try VivoMDMetalABI.command(packed: packed, configuration: configuration, cell: cell, stepIndex: 0)
        return .init(cell: cell, neighbors: grid, pme: pme, command: command)
    }
    private init(system: VivoClassicalSystem, configuration: VivoMDConfiguration, device: MTLDevice,
        catalog: NumiVivoPipelineCatalog, queue: MTLCommandQueue, packed: VivoMDPackedSystem, arena: VivoMDGPUArena,
        work: VivoMDWorkBuffers, pipelines: [NumiVivoKernel: NumiVivoPipeline], totalCharge: Double,
        cell: VivoMDCellResources, barostat: VivoMDBarostatEngine?, step: UInt64, time: Double) throws {
        self.system = system; self.configuration = configuration; self.device = device; self.catalog = catalog
        self.queue = queue; self.packed = packed; self.arena = arena; self.work = work; self.pipelines = pipelines
        totalChargeE = totalCharge; cellResources = cell; barostatEngine = barostat
        acceptedStep = step; acceptedTimePS = time
        systemFingerprint = packed.systemFingerprint; configurationFingerprint = try configuration.fingerprint()
        deviceName = device.name; deviceRegistryID = device.registryID
    }

    public func step() async throws -> VivoMDStepCertificate {
        try reserve(); defer { inFlight = false }
        try Task.checkCancellation()
        let nextTime = acceptedTimePS + configuration.timeStepPS
        guard acceptedStep < UInt64.max, nextTime.isFinite, nextTime > acceptedTimePS else {
            throw VivoMDRuntimeError.metal("MD clock/step cannot advance")
        }
        let phase = cellResources, abi = command(for: phase)
        let buffer = try makeCommand("step")
        try copyAccepted(into: buffer); try clear(buffer)
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try neighbors(buffer, position: arena.candidatePosition, phase: phase, abi: abi)
        try forces(buffer, position: arena.candidatePosition, phase: phase, abi: abi)
        try encode(.mdHalfKick, buffer, [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi)
        try projectVelocity(buffer, position: arena.candidatePosition, source: arena.candidateVelocity,
                            scratch: arena.velocityScratch, dynamics: arena.dynamics, abi: abi)
        if configuration.ensemble == .nve {
            try constrainedDrift(buffer, abi: abi)
        } else {
            var half = abi; half.dtPS *= 0.5
            try constrainedDrift(buffer, abi: half)
            try encode(.mdLangevin, buffer, [arena.candidateVelocity, arena.dynamics], abi)
            try projectVelocity(buffer, position: arena.candidatePosition, source: arena.candidateVelocity,
                                scratch: arena.velocityScratch, dynamics: arena.dynamics, abi: abi)
            try constrainedDrift(buffer, abi: half)
        }
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try neighbors(buffer, position: arena.candidatePosition, phase: phase, abi: abi)
        try forces(buffer, position: arena.candidatePosition, phase: phase, abi: abi)
        try encode(.mdHalfKick, buffer, [arena.candidateVelocity, arena.forceEnergy, arena.dynamics], abi)
        try projectVelocity(buffer, position: arena.candidatePosition, source: arena.candidateVelocity,
                            scratch: arena.velocityScratch, dynamics: arena.dynamics, abi: abi)
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try validateCandidate(buffer, abi: abi)
        try await complete(buffer); try Task.checkCancellation()
        let status = readStatus()
        guard status.flags == 0 else {
            return certificate(committed: false, status: status, barostat: nil, timeAfter: acceptedTimePS)
        }
        var finalPhase = phase
        var pressure: VivoMDBarostatMoveCertificate?
        if configuration.ensemble == .npt, (acceptedStep + 1) % UInt64(configuration.barostatInterval) == 0 {
            let result = try await pressureMove(phase: phase)
            finalPhase = result.phase; pressure = result.certificate
        }
        try Task.checkCancellation()
        let result = certificate(committed: true,
            status: .init(flags: 0, firstParticle: .max, violationCount: 0, reserved: 0),
            barostat: pressure, timeAfter: nextTime)
        // No suspension or fallible operation between state and clock publication.
        arena.commit(); cellResources = finalPhase; acceptedStep += 1; acceptedTimePS = nextTime
        return result
    }

    public func thermalize(temperatureK: Double, seed: UInt64) async throws -> VivoMDCheckpoint {
        try reserve(); defer { inFlight = false }; try Task.checkCancellation()
        guard temperatureK.isFinite, temperatureK > 0, Float(temperatureK).isFinite, Float(temperatureK) > 0 else {
            throw VivoMDRuntimeError.metal("invalid initialization temperature")
        }
        var abi = command(for: cellResources)
        abi.langevinA = 0; abi.targetTemperatureK = Float(temperatureK)
        let namespace = seed ^ 0x544845524d414c31
        abi.seedLow = UInt32(truncatingIfNeeded: namespace); abi.seedHigh = UInt32(truncatingIfNeeded: namespace >> 32)
        let buffer = try makeCommand("thermalize")
        try copyAccepted(into: buffer); try clear(buffer)
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdLangevin, buffer, [arena.candidateVelocity, arena.dynamics], abi)
        try projectVelocity(buffer, position: arena.candidatePosition, source: arena.candidateVelocity,
                            scratch: arena.velocityScratch, dynamics: arena.dynamics, abi: abi)
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdClearForce, buffer, [arena.forceEnergy], abi)
        try validateCandidate(buffer, abi: abi)
        try await complete(buffer); try ensureNumericalSuccess(); try Task.checkCancellation()
        // Complete and validate even the requested readback before changing owners.
        let state = try await readSnapshotReserved(position: arena.candidatePosition, velocity: arena.candidateVelocity)
        let result = try checkpoint(from: state)
        try Task.checkCancellation(); arena.commit()
        return result
    }

    public func minimize(_ settings: VivoMDMinimizationConfiguration = .init()) async throws -> VivoMDMinimizationCertificate {
        try settings.validate(); try reserve(); defer { inFlight = false }; try Task.checkCancellation()
        let phase = cellResources, abi = command(for: phase)
        let start = try makeCommand("minimize.initialProjection")
        try copyAccepted(into: start); try clear(start)
        try encode(.mdZeroVelocity, start, [arena.candidateVelocity], abi)
        try projectPosition(start, abi: abi)
        try normalize(start, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdClearForce, start, [arena.forceEnergy], abi); try validateCandidate(start, abi: abi)
        try await complete(start); try ensureNumericalSuccess()
        var baseline = try await evaluate(position: arena.candidatePosition, velocity: arena.candidateVelocity, phase: phase, gradient: true)
        try Task.checkCancellation(); arena.commit()
        let initialEnergy = baseline.energy
        var scale = settings.initialStepScale
        var attempts: UInt32 = 0, accepted: UInt32 = 0, rejected: UInt32 = 0
        while attempts < settings.maximumIterations,
              baseline.maximumForce > settings.forceToleranceKJPerMolNM, scale >= settings.minimumStepScale {
            try Task.checkCancellation(); attempts += 1
            let trial = try makeCommand("minimize.trial")
            try copyAccepted(into: trial); try clear(trial)
            try minimizePosition(trial, abi: abi, scale: scale, maximum: settings.maximumDisplacementNM)
            try projectPosition(trial, abi: abi)
            try normalize(trial, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
            try validateCandidate(trial, abi: abi)
            try await complete(trial)
            var candidate: (energy: Double, maximumForce: Double)?
            if readStatus().flags == 0 {
                do {
                    candidate = try await evaluate(position: arena.candidatePosition, velocity: arena.candidateVelocity,
                                                   phase: phase, gradient: true)
                } catch VivoMDRuntimeError.candidateRejected(_) { candidate = nil }
            }
            try Task.checkCancellation()
            if let candidate, candidate.energy < baseline.energy {
                arena.commit(); accepted += 1; baseline = candidate
                scale = min(settings.maximumStepScale, scale * settings.acceptedStepGrowth)
            } else {
                rejected += 1; scale *= settings.rejectedStepShrink
                baseline = try await evaluate(position: arena.acceptedPosition, velocity: arena.acceptedVelocity,
                                               phase: phase, gradient: true)
            }
        }
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
            converged: baseline.maximumForce <= settings.forceToleranceKJPerMolNM,
            attemptedIterations: attempts, acceptedIterations: accepted, rejectedIterations: rejected,
            initialPotentialEnergyKJPerMol: initialEnergy, finalPotentialEnergyKJPerMol: baseline.energy,
            finalMaximumForceKJPerMolNM: baseline.maximumForce, finalStepScale: scale)
    }

    public func snapshot() async throws -> VivoMDStateSnapshot {
        try reserve(); defer { inFlight = false }
        return try await readSnapshotReserved()
    }
    public func checkpoint() async throws -> VivoMDCheckpoint {
        let state = try await snapshot()
        return try checkpoint(from: state)
    }
    public func observables() async throws -> VivoMDObservables {
        try reserve(); defer { inFlight = false }
        return try await observationsReserved()
    }
    public func sample(includeObservables: Bool) async throws -> (state: VivoMDStateSnapshot, observables: VivoMDObservables?) {
        try reserve(); defer { inFlight = false }
        let state = try await readSnapshotReserved()
        let observations = includeObservables ? try await observationsReserved() : nil
        return (state, observations)
    }

    private func pressureMove(phase: VivoMDCellResources) async throws
    -> (phase: VivoMDCellResources, certificate: VivoMDBarostatMoveCertificate) {
        guard let engine = barostatEngine, let oldCell = phase.cell, let temperature = configuration.targetTemperatureK,
              let pressure = configuration.targetPressureBar else { throw VivoMDRuntimeError.metal("incomplete NPT state") }
        let step = acceptedStep + 1, oldVolume = oldCell.volumeNM3
        let deltaLogVolume = (2 * uniform(seed: configuration.randomSeed, step: step, stream: 0x4e505456) - 1)
            * configuration.resolvedBarostatMaximumLogVolumeStep
        let scale = exp(deltaLogVolume / 3)
        let proposedCell = VivoPeriodicCell(a: oldCell.a * scale, b: oldCell.b * scale, c: oldCell.c * scale)
        func rejection(_ reason: String) -> VivoMDBarostatMoveCertificate {
            .init(stepIndex: step, attempted: true, accepted: false,
                  volumeBeforeNM3: oldVolume, volumeAfterNM3: oldVolume,
                  logVolumeDelta: deltaLogVolume, logAcceptanceRatio: -Double.greatestFiniteMagnitude,
                  rejectionReason: reason)
        }
        do { try VivoMDExecutionPreflight.validateCell(proposedCell, configuration: configuration) }
        catch { return (phase, rejection("proposalOutsideCellContract")) }
        let proposedVolume = proposedCell.volumeNM3
        guard proposedVolume.isFinite, proposedVolume > 0 else { return (phase, rejection("invalidVolume")) }
        let proposalPhase = try await Self.makeCellResources(device: device, catalog: catalog, packed: packed,
            arena: arena, configuration: configuration, totalCharge: totalChargeE, cell: proposedCell)
        let baseline = try await evaluate(position: arena.candidatePosition, velocity: arena.candidateVelocity, phase: phase)
        let proposal = try makeCommand("npt.proposal")
        try copy(proposal, arena.candidatePosition, arena.positionScratch)
        try copy(proposal, arena.candidateVelocity, arena.velocityScratch)
        try clear(proposal)
        try engine.encodeProposal(commandBuffer: proposal, positions: arena.candidatePosition, dynamics: arena.dynamics,
                                  cell: oldCell, scale: scale, status: arena.status)
        let newABI = command(for: proposalPhase)
        try normalize(proposal, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: newABI)
        try validateCandidate(proposal, abi: newABI)
        try await complete(proposal)
        var result: VivoMDBarostatMoveCertificate
        if readStatus().flags != 0 {
            result = rejection("molecularGeometryOrConstraints")
        } else {
            do {
                let trial = try await evaluate(position: arena.candidatePosition, velocity: arena.candidateVelocity, phase: proposalPhase)
                let score = -(trial.energy - baseline.energy + pressure * (proposedVolume - oldVolume) * 0.0602214076)
                    / (0.00831446261815324 * temperature)
                    + (Double(engine.plan.componentCount) + 1) * deltaLogVolume
                guard score.isFinite else { throw VivoMDRuntimeError.candidateRejected(1) }
                let take = log(uniform(seed: configuration.randomSeed, step: step, stream: 0x4e505441)) < min(0, score)
                result = .init(stepIndex: step, attempted: true, accepted: take, volumeBeforeNM3: oldVolume,
                    volumeAfterNM3: take ? proposedVolume : oldVolume,
                    logVolumeDelta: deltaLogVolume, logAcceptanceRatio: score)
                if take { return (proposalPhase, result) }
            } catch VivoMDRuntimeError.candidateRejected(_) { result = rejection("invalidTrialEnergy") }
        }
        let restore = try makeCommand("npt.discard")
        try copy(restore, arena.positionScratch, arena.candidatePosition)
        try copy(restore, arena.velocityScratch, arena.candidateVelocity)
        try await complete(restore)
        return (phase, result)
    }

    private func evaluate(position: MTLBuffer, velocity: MTLBuffer, phase: VivoMDCellResources, gradient: Bool = false) async throws
    -> (energy: Double, maximumForce: Double) {
        let abi = command(for: phase), buffer = try makeCommand("evaluate")
        try clear(buffer); try neighbors(buffer, position: position, phase: phase, abi: abi)
        try forces(buffer, position: position, phase: phase, abi: abi)
        try encode(.mdValidate, buffer, [position, velocity, arena.forceEnergy, arena.status], abi)
        let direction: MTLBuffer
        if gradient {
            try encode(.mdProjectedForceSeed, buffer, [arena.forceEnergy, arena.dynamics, work.directionA], abi)
            try projectVelocity(buffer, position: position, source: work.directionA, scratch: work.directionB,
                                dynamics: work.unitDynamics, abi: abi)
            if !packed.constraints.isEmpty {
                try encode(.mdValidateProjectedForce, buffer, [position, work.directionA, arena.constraints,
                    arena.constraintOffsets, arena.constraintIncidence, arena.status], abi)
            }
            direction = work.directionA
        } else { direction = arena.forceEnergy }
        try encode(.mdMinimizeTerms, buffer, [arena.forceEnergy, direction, arena.dynamics, work.termsA], abi)
        try reduce(buffer, sumBoth: false)
        try await complete(buffer); try ensureNumericalSuccess()
        let value = work.scalar.contents().assumingMemoryBound(to: SIMD2<Float>.self).pointee
        guard value.x.isFinite, value.y.isFinite, value.y >= 0 else { throw VivoMDRuntimeError.candidateRejected(1) }
        return (Double(value.x), Double(value.y))
    }
    private func observationsReserved() async throws -> VivoMDObservables {
        let abi = command(for: cellResources), buffer = try makeCommand("observables")
        try clear(buffer); try neighbors(buffer, position: arena.acceptedPosition, phase: cellResources, abi: abi)
        try forces(buffer, position: arena.acceptedPosition, phase: cellResources, abi: abi)
        try encode(.mdValidate, buffer, [arena.acceptedPosition, arena.acceptedVelocity, arena.forceEnergy, arena.status], abi)
        try encode(.mdObservationTerms, buffer, [arena.forceEnergy, arena.acceptedVelocity, arena.dynamics, work.termsA], abi)
        try reduce(buffer, sumBoth: true); try await complete(buffer); try ensureNumericalSuccess()
        let value = work.scalar.contents().assumingMemoryBound(to: SIMD2<Float>.self).pointee
        let potential = Double(value.x), kinetic = Double(value.y)
        let massive = UInt64(system.particles.filter { $0.massDa > 0 }.count)
        let constraints = UInt64(system.constraints.count)
        guard massive * 3 > constraints else { throw VivoMDRuntimeError.metal("nonpositive constrained degrees of freedom") }
        let dof = massive * 3 - constraints
        let temperature = 2 * kinetic / (Double(dof) * 0.00831446261815324)
        guard potential.isFinite, kinetic.isFinite, kinetic >= 0, temperature.isFinite else { throw VivoMDRuntimeError.candidateRejected(1) }
        return .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
            stepIndex: acceptedStep, timePS: acceptedTimePS, potentialEnergyKJPerMol: potential,
            kineticEnergyKJPerMol: kinetic, temperatureK: temperature, degreesOfFreedom: dof)
    }
    private func reduce(_ buffer: MTLCommandBuffer, sumBoth: Bool) throws {
        var current = work.termsA, scratch = work.termsB, count = UInt32(arena.particleCount)
        while count > 1 {
            let output = count / 2 + count % 2
            let kernel: NumiVivoKernel = sumBoth ? .mdSumPairReduce : .mdMinimizeReduceStage
            guard let pipeline = pipelines[kernel], let encoder = buffer.makeComputeCommandEncoder() else {
                throw VivoMDRuntimeError.metal("reduction encoder unavailable")
            }
            encoder.setComputePipelineState(pipeline.state)
            encoder.setBuffer(current, offset: 0, index: 0); encoder.setBuffer(scratch, offset: 0, index: 1)
            var command = SIMD4<UInt32>(count, 0, 0, 0)
            encoder.setBytes(&command, length: 16, index: 2)
            encoder.dispatchThreads(pipeline.gridSize(for: Int(output)), threadsPerThreadgroup: pipeline.threadgroupSize(for: Int(output)))
            encoder.endEncoding(); swap(&current, &scratch); count = output
        }
        guard let blit = buffer.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("reduction readback unavailable") }
        blit.copy(from: current, sourceOffset: 0, to: work.scalar, destinationOffset: 0, size: 8); blit.endEncoding()
    }

    private func constrainedDrift(_ buffer: MTLCommandBuffer, abi: VivoMDMetalCommand) throws {
        try encode(.mdDrift, buffer, [arena.candidatePosition, arena.candidateVelocity, arena.dynamics], abi)
        if !packed.constraints.isEmpty {
            try copy(buffer, arena.candidatePosition, work.referencePosition)
            try projectPosition(buffer, abi: abi)
            try encode(.mdDriftConstraintImpulse, buffer, [arena.candidatePosition, work.referencePosition,
                arena.candidateVelocity, arena.dynamics, arena.status], abi)
        }
    }
    private func projectPosition(_ buffer: MTLCommandBuffer, abi: VivoMDMetalCommand) throws {
        guard !packed.constraints.isEmpty else { return }
        var source = arena.candidatePosition, destination = arena.positionScratch
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintPosition, buffer, [source, destination, arena.dynamics, arena.constraints,
                arena.constraintOffsets, arena.constraintIncidence, arena.status], abi)
            swap(&source, &destination)
        }
        if configuration.maximumConstraintIterations % 2 != 0 { try copy(buffer, source, arena.candidatePosition) }
    }
    private func projectVelocity(_ buffer: MTLCommandBuffer, position: MTLBuffer, source original: MTLBuffer,
                                  scratch: MTLBuffer, dynamics: MTLBuffer, abi: VivoMDMetalCommand) throws {
        guard !packed.constraints.isEmpty else { return }
        var source = original, destination = scratch
        for _ in 0..<configuration.maximumConstraintIterations {
            try encode(.mdConstraintVelocity, buffer, [position, source, destination, dynamics, arena.constraints,
                arena.constraintOffsets, arena.constraintIncidence, arena.status], abi)
            swap(&source, &destination)
        }
        if configuration.maximumConstraintIterations % 2 != 0 { try copy(buffer, source, original) }
    }
    private func validateCandidate(_ buffer: MTLCommandBuffer, abi: VivoMDMetalCommand) throws {
        if !packed.constraints.isEmpty {
            try encode(.mdValidateConstraints, buffer, [arena.candidatePosition, arena.candidateVelocity, arena.constraints,
                arena.constraintOffsets, arena.constraintIncidence, arena.status], abi)
        }
        try encode(.mdValidate, buffer, [arena.candidatePosition, arena.candidateVelocity, arena.forceEnergy, arena.status], abi)
    }
    private func normalize(_ buffer: MTLCommandBuffer, position: MTLBuffer, velocity: MTLBuffer, abi: VivoMDMetalCommand) throws {
        guard !packed.virtualSites.isEmpty else { return }
        try encode(.mdUpdateVirtualPosition, buffer, [position, arena.virtualSites, arena.virtualSiteIndexByParticle, arena.status], abi)
        try encode(.mdUpdateVirtualVelocity, buffer, [velocity, arena.virtualSites, arena.virtualSiteIndexByParticle], abi)
    }
    private func neighbors(_ buffer: MTLCommandBuffer, position: MTLBuffer, phase: VivoMDCellResources, abi: VivoMDMetalCommand) throws {
        guard configuration.resolvedNeighborListEnabled else { return }
        guard abi.neighborCapacity == arena.neighborCapacity else { throw VivoMDRuntimeError.metal("neighbor ABI/allocation mismatch") }
        if let grid = phase.neighbors {
            try grid.encode(commandBuffer: buffer, positions: position, neighborCounts: arena.neighborCounts,
                neighborIndices: arena.neighborIndices, referencePositions: arena.neighborReferencePosition, status: arena.status)
        } else {
            try encode(.mdBuildNeighborList, buffer, [position, arena.neighborCounts, arena.neighborIndices,
                arena.neighborReferencePosition, arena.status], abi)
        }
    }
    private func forces(_ buffer: MTLCommandBuffer, position: MTLBuffer, phase: VivoMDCellResources, abi: VivoMDMetalCommand) throws {
        try encode(.mdClearForce, buffer, [arena.forceEnergy], abi)
        try encode(.mdBonded, buffer, [position, arena.forceEnergy, arena.bonds, arena.bondOffsets, arena.bondIncidence,
            arena.angles, arena.angleOffsets, arena.angleIncidence, arena.torsions, arena.torsionOffsets, arena.torsionIncidence, arena.status], abi)
        var pairBuffers = [position, arena.forceEnergy, arena.dynamics, arena.typeIndices, arena.pairMatrix,
                           arena.exceptions, arena.exceptionOffsets, arena.exceptionPartners, arena.exceptionIndices]
        if configuration.resolvedNeighborListEnabled { pairBuffers += [arena.neighborCounts, arena.neighborIndices] }
        pairBuffers.append(arena.status)
        if configuration.electrostatics == .pme {
            guard let pme = phase.pme else { throw VivoMDRuntimeError.metal("PME engine absent") }
            try encode(.mdPMERealSpaceNeighbor, buffer, pairBuffers, abi)
            try pme.encodeReciprocal(commandBuffer: buffer, positions: position, dynamics: arena.dynamics,
                                    forceEnergy: arena.forceEnergy, status: arena.status)
            if !packed.pairExceptions.isEmpty {
                try encode(.mdPMEExceptionCorrection, buffer, [position, arena.dynamics, arena.forceEnergy, arena.exceptions,
                    arena.exceptionOffsets, arena.exceptionPartners, arena.exceptionIndices, arena.status], abi)
            }
        } else {
            try encode(configuration.resolvedNeighborListEnabled ? .mdNonbondedNeighbor : .mdNonbondedDirect, buffer, pairBuffers, abi)
        }
        if !packed.virtualSites.isEmpty {
            try encode(.mdRedistributeVirtualForce, buffer, [arena.forceEnergy, arena.virtualSites,
                arena.virtualParentOffsets, arena.virtualParentIncidence, arena.status], abi)
        }
    }
    private func minimizePosition(_ buffer: MTLCommandBuffer, abi: VivoMDMetalCommand, scale: Double, maximum: Double) throws {
        guard Float(scale).isFinite, Float(scale) > 0, Float(maximum).isFinite, Float(maximum) > 0,
              let pipeline = pipelines[.mdMinimizePosition], let encoder = buffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("minimization scale or encoder invalid")
        }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state)
        encoder.setBuffer(arena.candidatePosition, offset: 0, index: 0)
        encoder.setBuffer(work.directionA, offset: 0, index: 1); encoder.setBuffer(arena.dynamics, offset: 0, index: 2)
        var md = abi
        var command = VivoMDMinPositionCommand(particleCount: packed.particleCount, stepScale: Float(scale), maxDisplacementNM: Float(maximum))
        encoder.setBytes(&md, length: MemoryLayout<VivoMDMetalCommand>.stride, index: 3)
        encoder.setBytes(&command, length: 16, index: 4)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount), threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
    }

    private func establishDerivedState() async throws {
        try reserve(); defer { inFlight = false }
        let buffer = try makeCommand("initializeDerived"), abi = command(for: cellResources)
        try copyAccepted(into: buffer); try clear(buffer)
        try normalize(buffer, position: arena.candidatePosition, velocity: arena.candidateVelocity, abi: abi)
        try encode(.mdClearForce, buffer, [arena.forceEnergy], abi)
        try encode(.mdValidate, buffer, [arena.candidatePosition, arena.candidateVelocity, arena.forceEnergy, arena.status], abi)
        try await complete(buffer); try ensureNumericalSuccess(); arena.commit()
    }
    private func readSnapshotReserved(position: MTLBuffer? = nil, velocity: MTLBuffer? = nil) async throws -> VivoMDStateSnapshot {
        let buffer = try makeCommand("snapshot")
        try copy(buffer, position ?? arena.acceptedPosition, arena.positionReadback)
        try copy(buffer, velocity ?? arena.acceptedVelocity, arena.velocityReadback)
        try await complete(buffer)
        let positionsPointer = arena.positionReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        let velocitiesPointer = arena.velocityReadback.contents().assumingMemoryBound(to: SIMD4<Float>.self)
        var positions: [VivoVector3D] = [], velocities: [VivoVector3D] = []
        positions.reserveCapacity(arena.particleCount); velocities.reserveCapacity(arena.particleCount)
        for i in 0..<arena.particleCount {
            let p = positionsPointer[i], v = velocitiesPointer[i]
            positions.append(.init(Double(p.x), Double(p.y), Double(p.z)))
            velocities.append(.init(Double(v.x), Double(v.y), Double(v.z)))
        }
        let state = VivoMDStateSnapshot(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint,
            stepIndex: acceptedStep, timePS: acceptedTimePS, positionsNM: positions, velocitiesNMPerPS: velocities, periodicCell: cellResources.cell)
        try state.validate(particleCount: arena.particleCount)
        return state
    }
    private func checkpoint(from state: VivoMDStateSnapshot) throws -> VivoMDCheckpoint {
        let result = VivoMDCheckpoint(systemFingerprint: state.systemFingerprint, configurationFingerprint: state.configurationFingerprint,
            acceptedStep: state.stepIndex, timePS: state.timePS, positionsNM: state.positionsNM,
            velocitiesNMPerPS: state.velocitiesNMPerPS, periodicCell: state.periodicCell)
        try result.validate(particleCount: arena.particleCount)
        return result
    }
    private func certificate(committed: Bool, status: VivoMDMetalStatus,
                             barostat: VivoMDBarostatMoveCertificate?, timeAfter: Double) -> VivoMDStepCertificate {
        .init(systemFingerprint: systemFingerprint, configurationFingerprint: configurationFingerprint, deviceName: deviceName,
            deviceRegistryID: deviceRegistryID, stepIndex: acceptedStep, timeBeforePS: acceptedTimePS, timeAfterPS: timeAfter,
            committed: committed, statusFlags: status.flags,
            firstViolationParticle: status.firstParticle == .max ? nil : status.firstParticle,
            violationCount: status.violationCount, barostat: barostat)
    }
    private func command(for phase: VivoMDCellResources) -> VivoMDMetalCommand {
        var value = phase.command
        value.stepLow = UInt32(truncatingIfNeeded: acceptedStep); value.stepHigh = UInt32(truncatingIfNeeded: acceptedStep >> 32)
        return value
    }
    private func uniform(seed: UInt64, step: UInt64, stream: UInt64) -> Double {
        var value = seed ^ (step &* 0x9E3779B97F4A7C15) ^ stream
        value &+= 0x9E3779B97F4A7C15
        value = (value ^ (value >> 30)) &* 0xBF58476D1CE4E5B9
        value = (value ^ (value >> 27)) &* 0x94D049BB133111EB
        value ^= value >> 31
        return (Double(value >> 12) + 0.5) * 0x1.0p-52
    }
    private func reserve() throws {
        guard !inFlight else { throw VivoMDRuntimeError.metal("MD operation already in flight") }
        if let deviceFailure { throw VivoMDRuntimeError.metal("runtime stopped after device failure: " + deviceFailure) }
        inFlight = true
    }
    private func makeCommand(_ label: String) throws -> MTLCommandBuffer {
        guard let buffer = queue.makeCommandBuffer() else { throw VivoMDRuntimeError.metal("command buffer unavailable") }
        buffer.label = "NumiVivo.MD.\(acceptedStep).\(label)"
        return buffer
    }
    private func copyAccepted(into buffer: MTLCommandBuffer) throws {
        try copy(buffer, arena.acceptedPosition, arena.candidatePosition)
        try copy(buffer, arena.acceptedVelocity, arena.candidateVelocity)
    }
    private func copy(_ buffer: MTLCommandBuffer, _ source: MTLBuffer, _ destination: MTLBuffer) throws {
        guard let blit = buffer.makeBlitCommandEncoder() else { throw VivoMDRuntimeError.metal("blit encoder unavailable") }
        blit.copy(from: source, sourceOffset: 0, to: destination, destinationOffset: 0, size: arena.particleCount * 16)
        blit.endEncoding()
    }
    private func clear(_ buffer: MTLCommandBuffer) throws {
        guard let pipeline = pipelines[.mdClearStatus], let encoder = buffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("clear encoder unavailable")
        }
        encoder.setComputePipelineState(pipeline.state); encoder.setBuffer(arena.status, offset: 0, index: 0)
        encoder.dispatchThreads(.init(width: 1, height: 1, depth: 1), threadsPerThreadgroup: .init(width: 1, height: 1, depth: 1))
        encoder.endEncoding()
    }
    private func encode(_ kernel: NumiVivoKernel, _ buffer: MTLCommandBuffer,
                        _ buffers: [MTLBuffer], _ abi: VivoMDMetalCommand) throws {
        guard let pipeline = pipelines[kernel], let encoder = buffer.makeComputeCommandEncoder() else {
            throw VivoMDRuntimeError.metal("missing encoder for \(kernel.rawValue)")
        }
        defer { encoder.endEncoding() }
        encoder.setComputePipelineState(pipeline.state); encoder.label = kernel.rawValue
        for (index, resource) in buffers.enumerated() { encoder.setBuffer(resource, offset: 0, index: index) }
        var value = abi
        encoder.setBytes(&value, length: MemoryLayout<VivoMDMetalCommand>.stride, index: buffers.count)
        encoder.dispatchThreads(pipeline.gridSize(for: arena.particleCount), threadsPerThreadgroup: pipeline.threadgroupSize(for: arena.particleCount))
    }
    private func readStatus() -> VivoMDMetalStatus { arena.status.contents().assumingMemoryBound(to: VivoMDMetalStatus.self).pointee }
    private func ensureNumericalSuccess() throws {
        let status = readStatus()
        if status.flags != 0 { throw VivoMDRuntimeError.candidateRejected(status.flags) }
    }
    private func complete(_ buffer: MTLCommandBuffer) async throws {
        do {
            try await withCheckedThrowingContinuation { (continuation: CheckedContinuation<Void, Error>) in
                buffer.addCompletedHandler { completed in
                    if completed.status != .completed {
                        continuation.resume(throwing: VivoMDRuntimeError.metal(
                            completed.error.map { String(describing: $0) } ?? "command did not complete"))
                    } else { continuation.resume(returning: ()) }
                }
                buffer.commit()
            }
        } catch { deviceFailure = String(describing: error); throw error }
    }
}
