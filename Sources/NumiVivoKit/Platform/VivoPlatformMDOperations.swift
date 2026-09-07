import Foundation

public struct VivoWorkflowMDConfiguration: Codable, Sendable, Equatable {
    public var dynamics: VivoMDConfiguration
    public var steps: Int
    public var maximumParticles: Int
    public var maximumMeshPoints: Int
    public init(dynamics: VivoMDConfiguration, steps: Int, maximumParticles: Int = 100_000,
                maximumMeshPoints: Int = 4_194_304) {
        self.dynamics = dynamics; self.steps = steps; self.maximumParticles = maximumParticles
        self.maximumMeshPoints = maximumMeshPoints
    }
    public func validate() throws {
        try dynamics.validate()
        guard (1...100_000).contains(steps), (1...1_000_000).contains(maximumParticles),
              (1...67_108_864).contains(maximumMeshPoints) else { throw VivoChemistryError.invalid("bounded workflow MD stage settings") }
        // NPT changes the reciprocal mesh allocation during a stage. The existing
        // protocol runner owns that adaptive path; a fixed-memory recipe stage
        // must not claim a constant allocation reservation for it.
        guard dynamics.ensemble != .npt else { throw VivoChemistryError.unsupported("use md-protocol-run for volume-adaptive NPT; workflow MD segments currently support NVE/NVT") }
    }
}
public struct VivoWorkflowMDStageResult: Codable, Sendable, Equatable {
    public let configuration: VivoWorkflowMDConfiguration
    public let start: VivoMDCheckpoint
    public let end: VivoMDCheckpoint
    public let observables: VivoMDObservables
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let allocationReservationBytes: Int
    public let validationScope: String
}
public enum VivoPlatformMDOperations {
    private static let scope = "accepted native Metal segment; checkpoint source, clock, shape and finite-observable validation; cached trajectories are integrity-checked, not numerically replayed; not equilibrium, long-time drift or cross-device reproducibility qualification"
    private struct Inputs {
        let system: VivoClassicalSystem
        let initial: VivoClassicalInitialState
        let checkpoint: VivoMDCheckpoint?
        let configuration: VivoWorkflowMDConfiguration
        let reservation: Int
    }
    private static func prepare(_ cfg: VivoJSONValue, _ data: [String: Data], _ budget: VivoChemistryBudget, resume: Bool) throws -> Inputs {
        let c = try VivoPlatformOperations.decode(VivoWorkflowMDConfiguration.self, cfg)
        try c.validate(); try budget.validate()
        guard Set(data.keys) == Set(["system", resume ? "checkpoint" : "state"]) else { throw VivoChemistryError.invalid("MD workflow input slots") }
        let system = try VivoPlatformOperations.input(VivoClassicalSystem.self, "system", data)
        try VivoClassicalSystemValidator.validate(system)
        guard system.particles.count <= c.maximumParticles else { throw VivoChemistryError.resourceLimit("MD workflow particle capacity") }
        let initial: VivoClassicalInitialState, checkpoint: VivoMDCheckpoint?
        if resume {
            let value = try VivoPlatformOperations.input(VivoMDCheckpoint.self, "checkpoint", data)
            try value.validate(particleCount: system.particles.count)
            guard value.systemFingerprint == (try system.fingerprint()), value.configurationFingerprint == (try c.dynamics.fingerprint()),
                  value.acceptedStep <= UInt64.max - UInt64(c.steps) else { throw VivoChemistryError.invalid("MD workflow restart source/configuration/step") }
            checkpoint = value
            initial = .init(systemFingerprint: value.systemFingerprint, positionsNM: value.positionsNM,
                            periodicCell: value.periodicCell, sourceTimePS: value.timePS)
        } else {
            initial = try VivoPlatformOperations.input(VivoClassicalInitialState.self, "state", data); checkpoint = nil
        }
        let support = try VivoMDCapabilityAnalyzer.analyze(system: system, initialState: initial, configuration: c.dynamics)
        guard support.executable else { throw VivoMDRuntimeError.unsupported(support.blockers) }
        let n = system.particles.count, types = Set(system.particles.map(\.typeIdentifier)).count
        var memory = 0
        func reserve(_ count: Int, _ stride: Int) throws {
            let (bytes, overflow) = count.multipliedReportingOverflow(by: stride)
            guard !overflow, bytes >= 0, bytes <= budget.maximumBytes - memory else { throw VivoChemistryError.resourceLimit("MD workflow allocation reservation") }
            memory += bytes
        }
        // Conservative packed/host/device tables, staging and scratch panels.
        try reserve(n, 2048)
        try reserve(types * types, 32)
        let neighbors = min(max(n - 1, 1), Int(c.dynamics.resolvedMaximumNeighborsPerParticle))
        try reserve(n, neighbors * 16)
        let terms = system.bonds.count + system.angles.count + system.torsions.count + system.constraints.count
            + system.nonbondedExceptions.count + (system.linearVirtualSites?.count ?? 0)
        try reserve(terms, 512)
        if c.dynamics.electrostatics == .pme, let cell = initial.periodicCell {
            if c.dynamics.pmeGridDimensions == nil {
                for vector in [cell.a, cell.b, cell.c] {
                    let points = vector.norm / c.dynamics.resolvedPMEGridSpacingNM
                    guard points.isFinite, points <= 65_536 else { throw VivoChemistryError.resourceLimit("MD workflow PME axis capacity") }
                }
            }
            let mesh = try VivoPMEPlan.make(cell: cell, cutoffNM: c.dynamics.cutoffNM,
                tolerance: c.dynamics.resolvedPMETolerance, targetGridSpacingNM: c.dynamics.resolvedPMEGridSpacingNM,
                fixedGridDimensions: c.dynamics.pmeGridDimensions)
            guard mesh.gridPointCount <= UInt64(c.maximumMeshPoints) else { throw VivoChemistryError.resourceLimit("MD workflow reciprocal mesh capacity") }
            try reserve(Int(mesh.gridPointCount), 256)
        }
        return .init(system: system, initial: initial, checkpoint: checkpoint, configuration: c, reservation: memory)
    }
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        [false, true].map { resume in
            let kinds = ["system": "vivo.classical-system", resume ? "checkpoint" : "state": resume ? "vivo.md-checkpoint" : "vivo.classical-initial-state"]
            let outputs: [VivoChemistryTaskOutput] = [.init(name: "checkpoint", kind: "vivo.md-checkpoint"), .init(name: "result", kind: "vivo.workflow-md-stage")]
            let operation = VivoChemistryOperation(identifier: resume ? "vivo.platform.md-continue" : "vivo.platform.md-start",
                version: "1", implementationFingerprint: id, numericalBackend: "metal-fp32", outputs: outputs,
                execute: { cfg, data, budget in
                    let p = try prepare(cfg, data, budget, resume: resume)
                    let runtime: VivoMDMetalRuntime
                    if let checkpoint = p.checkpoint {
                        runtime = try await .restore(system: p.system, configuration: p.configuration.dynamics, checkpoint: checkpoint)
                    } else { runtime = try await .make(system: p.system, initialState: p.initial, configuration: p.configuration.dynamics) }
                    let start = try await runtime.checkpoint()
                    for _ in 0..<p.configuration.steps {
                        try Task.checkCancellation()
                        let candidate = try await runtime.step()
                        guard candidate.committed else { throw VivoChemistryError.convergence("MD segment rejected; no completed checkpoint published") }
                    }
                    let end = try await runtime.checkpoint(), observables = try await runtime.observables()
                    let result = VivoWorkflowMDStageResult(configuration: p.configuration, start: start, end: end,
                        observables: observables, deviceName: runtime.deviceName, deviceRegistryID: runtime.deviceRegistryID,
                        allocationReservationBytes: p.reservation, validationScope: scope)
                    return ["checkpoint": try VivoCanonicalJSON.encode(end), "result": try VivoCanonicalJSON.encode(result)]
                }, validateOutputs: { cfg, data, output, budget in
                    let p = try prepare(cfg, data, budget, resume: resume)
                    guard Set(output.keys) == Set(outputs.map(\.name)) else { throw VivoChemistryError.invalid("MD segment output slots") }
                    let result = try VivoPlatformOperations.input(VivoWorkflowMDStageResult.self, "result", output)
                    let end = try VivoPlatformOperations.input(VivoMDCheckpoint.self, "checkpoint", output)
                    try result.start.validate(particleCount: p.system.particles.count); try end.validate(particleCount: p.system.particles.count)
                    let system = try p.system.fingerprint(), config = try p.configuration.dynamics.fingerprint()
                    let firstStep = p.checkpoint?.acceptedStep ?? 0, firstTime = p.initial.sourceTimePS ?? 0
                    var time = firstTime
                    for _ in 0..<p.configuration.steps { time += p.configuration.dynamics.timeStepPS }
                    let obs = result.observables
                    guard result.configuration == p.configuration, result.end == end,
                          result.start.systemFingerprint == system, result.start.configurationFingerprint == config,
                          result.start.acceptedStep == firstStep, result.start.timePS == firstTime,
                          p.checkpoint == nil || result.start == p.checkpoint,
                          end.systemFingerprint == system, end.configurationFingerprint == config,
                          end.acceptedStep == firstStep + UInt64(p.configuration.steps), end.timePS == time,
                          end.periodicCell == result.start.periodicCell,
                          obs.systemFingerprint == system, obs.configurationFingerprint == config, obs.stepIndex == end.acceptedStep,
                          obs.timePS == time, obs.potentialEnergyKJPerMol.isFinite, obs.kineticEnergyKJPerMol.isFinite,
                          obs.kineticEnergyKJPerMol >= 0, obs.temperatureK?.isFinite != false,
                          obs.totalEnergyKJPerMol == obs.potentialEnergyKJPerMol + obs.kineticEnergyKJPerMol,
                          !result.deviceName.isEmpty, result.allocationReservationBytes == p.reservation, result.validationScope == scope else {
                        throw VivoChemistryError.invalid("MD segment checkpoint, physical source, clock or observable binding")
                    }
                })
            return VivoWorkflowDefinition(operation: operation, inputKinds: kinds,
                summary: resume ? "Continue an accepted NVE/NVT Metal checkpoint without resetting its state or random namespace." : "Start a bounded NVE/NVT Metal segment and publish an accepted checkpoint.",
                validationScope: scope, validateConfiguration: { try VivoPlatformOperations.decode(VivoWorkflowMDConfiguration.self, $0).validate() })
        }
    }
}
