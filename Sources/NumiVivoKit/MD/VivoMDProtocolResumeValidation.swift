import Foundation

/// Validates the persisted protocol before allocating a device runtime. Hashes
/// protect bytes; these checks establish the relationships between those bytes.
enum VivoMDProtocolResumeValidation {
    struct Prefix: Sendable {
        let cursor: VivoMDProtocolCheckpoint
        let checkpoint: VivoMDCheckpoint
        let lastSampledStep: UInt64?
        let lastObservedStep: UInt64?
        let trajectoryValidation: VivoMDTrajectoryValidation?
    }

    static func load(system: VivoClassicalSystem, plan: VivoMDProtocolPlan,
                     store: VivoArtifactStore, checkpoint: VivoFingerprint) async throws -> Prefix {
        try plan.validate()
        guard plan.systemFingerprint == (try system.fingerprint()) else { throw invalid("system identity") }
        let cursor = try await read(VivoMDProtocolCheckpoint.self, checkpoint, "md-protocol-checkpoint", store)
        try cursor.validate(plan: plan)
        let initial = try await read(VivoClassicalInitialState.self, cursor.initialStateArtifact,
                                     "classical-initial-state", store, maximumBytes: 512 * 1024 * 1024)
        try initial.validate(particleCount: system.particles.count)
        guard initial.systemFingerprint == plan.systemFingerprint else { throw invalid("initial-state identity") }

        var previous: VivoMDCheckpoint?
        for (index, hash) in cursor.priorStageReports.enumerated() {
            try Task.checkCancellation()
            let report = try await read(VivoMDProtocolStageReport.self, hash, "md-stage-report", store)
            let entry = try await state(report.entryCheckpoint, system, plan.stages[index], store)
            let exit = try await state(report.exitCheckpoint, system, plan.stages[index], store)
            try await validateEntry(entry, transition: report.transition, previous: previous,
                                    index: index, initial: initial, system: system, plan: plan, store: store)
            try validateReport(report, index: index, successful: true, entry: entry, exit: exit, plan: plan)
            _ = try await validateOutputs(stage: plan.stages[index], completed: report.committedSteps,
                entry: entry, current: exit, manifest: report.trajectoryManifest,
                observationTail: report.observationTail, observationCount: report.observationCount,
                allowPending: false, system: system, plan: plan, store: store)
            previous = exit
        }

        let stage = plan.stages[cursor.stageIndex]
        let entry = try await state(cursor.entryCheckpoint, system, stage, store)
        let current = try await state(cursor.currentCheckpoint, system, stage, store)
        try validateProgress(entry: entry, current: current, completed: cursor.completedStepsInStage, stage: stage)
        try await validateEntry(entry, transition: cursor.transition, previous: previous,
                                index: cursor.stageIndex, initial: initial, system: system, plan: plan, store: store)
        if let hash = cursor.activeStageReport {
            let report = try await read(VivoMDProtocolStageReport.self, hash, "md-stage-report", store)
            guard report.entryCheckpoint == cursor.entryCheckpoint, report.exitCheckpoint == cursor.currentCheckpoint,
                  report.transition == cursor.transition, report.committedSteps == cursor.completedStepsInStage,
                  report.trajectoryManifest == cursor.trajectoryManifest,
                  report.observationTail == cursor.observationTail, report.observationCount == cursor.observationCount else {
                throw invalid("active report and cursor disagree")
            }
            try validateReport(report, index: cursor.stageIndex, successful: cursor.phase == .stageFinished,
                               entry: entry, exit: current, plan: plan)
        }
        let outputs = try await validateOutputs(stage: stage, completed: cursor.completedStepsInStage,
            entry: entry, current: current, manifest: cursor.trajectoryManifest,
            observationTail: cursor.observationTail, observationCount: cursor.observationCount,
            allowPending: cursor.phase == .running, system: system, plan: plan, store: store)
        return .init(cursor: cursor, checkpoint: current,
                     lastSampledStep: outputs.sample, lastObservedStep: outputs.observation,
                     trajectoryValidation: outputs.trajectory)
    }

    private static func state(_ hash: VivoFingerprint, _ system: VivoClassicalSystem,
                              _ stage: VivoMDProtocolStage, _ store: VivoArtifactStore) async throws -> VivoMDCheckpoint {
        let value = try await read(VivoMDCheckpoint.self, hash, "md-checkpoint", store,
                                   maximumBytes: 512 * 1024 * 1024)
        try value.validate(particleCount: system.particles.count)
        guard value.systemFingerprint == (try system.fingerprint()),
              value.configurationFingerprint == (try stage.configuration.fingerprint()) else {
            throw invalid("checkpoint system/configuration identity")
        }
        return value
    }

    private static func validateProgress(entry: VivoMDCheckpoint, current: VivoMDCheckpoint,
                                         completed: UInt64, stage: VivoMDProtocolStage) throws {
        guard completed <= stage.steps, completed <= UInt64.max - entry.acceptedStep,
              current.acceptedStep == entry.acceptedStep + completed,
              current.timePS >= entry.timePS,
              stage.configuration.ensemble == .npt || current.periodicCell == entry.periodicCell else {
            throw invalid("checkpoint progress, clock or cell")
        }
        if stage.kind == .minimization || completed == 0 {
            guard current.timePS == entry.timePS else { throw invalid("nonadvancing stage advanced physical time") }
        } else {
            guard current.timePS > entry.timePS else { throw invalid("dynamics did not advance physical time") }
            // The runtime accumulates one FP64 addition per committed step.
            // For positive operands, the standard gamma_n bound covers those
            // n roundings; two more cover this comparison's product and sum.
            // The factor two covers rounded bound evaluation and converting the
            // rounded expected sum back to an upper bound on the exact sum.
            // n <= 10^10, so gamma_(n+2) is below 1.2e-6. No O(n) clock replay.
            let count = Double(completed)
            let expected = entry.timePS + count * stage.configuration.timeStepPS
            let roundoff = (count + 2) * (Double.ulpOfOne / 2)
            let tolerance = 2 * (roundoff / (1 - roundoff) * max(expected, Double.leastNormalMagnitude))
            guard expected.isFinite, tolerance.isFinite, abs(current.timePS - expected) <= tolerance else {
                throw invalid("dynamics clock disagrees with its accepted steps and timestep")
            }
        }
    }

    private static func validateEntry(_ entry: VivoMDCheckpoint, transition: VivoFingerprint?,
                                      previous: VivoMDCheckpoint?, index: Int,
                                      initial: VivoClassicalInitialState, system: VivoClassicalSystem,
                                      plan: VivoMDProtocolPlan, store: VivoArtifactStore) async throws {
        let stage = plan.stages[index]
        let source: VivoMDCheckpoint
        if index == 0 {
            guard transition == nil, previous == nil, entry.acceptedStep == 0,
                  entry.timePS == (initial.sourceTimePS ?? 0), entry.periodicCell == initial.periodicCell else {
                throw invalid("first-stage source boundary")
            }
            source = .init(systemFingerprint: plan.systemFingerprint,
                configurationFingerprint: try stage.configuration.fingerprint(), acceptedStep: 0,
                timePS: initial.sourceTimePS ?? 0,
                positionsNM: initial.positionsNM.map { .init(Double(Float($0.x)), Double(Float($0.y)), Double(Float($0.z))) },
                velocitiesNMPerPS: Array(repeating: .zero, count: system.particles.count), periodicCell: initial.periodicCell)
        } else {
            guard let previous, let transition else { throw invalid("stage transition is missing") }
            let record = try await read(VivoMDStageTransition.self, transition, "md-stage-transition", store)
            let expected = try VivoMDStageTransfer.prepare(checkpoint: previous, source: plan.stages[index - 1].configuration,
                destination: stage.configuration, particleCount: system.particles.count,
                velocityInitialization: stage.velocityInitialization, thermalizationSeed: stage.thermalizationSeed)
            guard record == expected else { throw invalid("stage transition disagrees with the preceding accepted exit") }
            source = expected.destinationCheckpoint
            guard entry.acceptedStep == source.acceptedStep, entry.timePS == source.timePS,
                  entry.periodicCell == source.periodicCell else { throw invalid("stage entry reset the clock or cell") }
        }
        // Device setup reconstructs dependent sites. Only physical coordinates
        // are inputs; their identity cannot be excused by virtual-site rounding.
        for particle in system.particles where particle.role != .virtualSite {
            let i = Int(particle.index)
            guard entry.positionsNM[i] == source.positionsNM[i] else { throw invalid("stage entry changed physical positions") }
            if stage.velocityInitialization != .maxwellBoltzmann {
                guard entry.velocitiesNMPerPS[i] == source.velocitiesNMPerPS[i] else {
                    throw invalid("stage entry violated its velocity policy")
                }
            }
        }
    }

    private static func validateReport(_ report: VivoMDProtocolStageReport, index: Int, successful: Bool,
                                       entry: VivoMDCheckpoint, exit: VivoMDCheckpoint, plan: VivoMDProtocolPlan) throws {
        let stage = plan.stages[index]
        guard report.schema == "numivivo.org/md-stage-report/v1", report.planFingerprint == (try plan.fingerprint()),
              report.stageIndex == index, report.stageIdentifier == stage.identifier, report.successful == successful,
              report.startTimePS == entry.timePS, report.endTimePS == exit.timePS else { throw invalid("stage report binding") }
        try validateProgress(entry: entry, current: exit, completed: report.committedSteps, stage: stage)
        if stage.kind == .minimization {
            // Rejected trials re-evaluate the accepted geometry; unordered PME
            // sums do not promise bitwise energy monotonicity for that readback.
            // A rejected Double scale multiplication may also underflow to zero,
            // terminating the optimizer without advancing its accepted state.
            guard let settings = stage.minimization, let certificate = report.minimization, report.rejected == nil,
                  certificate.systemFingerprint == entry.systemFingerprint,
                  certificate.configurationFingerprint == entry.configurationFingerprint,
                  UInt64(certificate.attemptedIterations) == UInt64(certificate.acceptedIterations) + UInt64(certificate.rejectedIterations),
                  certificate.attemptedIterations <= settings.maximumIterations,
                  certificate.initialPotentialEnergyKJPerMol.isFinite, certificate.finalPotentialEnergyKJPerMol.isFinite,
                  certificate.finalMaximumForceKJPerMolNM.isFinite, certificate.finalMaximumForceKJPerMolNM >= 0,
                  certificate.finalStepScale.isFinite, certificate.finalStepScale >= 0,
                  certificate.finalStepScale > 0 || certificate.rejectedIterations > 0,
                  certificate.finalStepScale <= settings.maximumStepScale,
                  certificate.converged == (certificate.finalMaximumForceKJPerMolNM <= settings.forceToleranceKJPerMolNM),
                  successful == (certificate.converged || !stage.requireMinimizationConvergence) else {
                throw invalid("minimization certificate or required convergence gate")
            }
        } else {
            guard report.minimization == nil else { throw invalid("dynamics report contains a minimizer") }
            if successful {
                guard report.committedSteps == stage.steps, report.rejected == nil else { throw invalid("unfinished dynamics report") }
            } else {
                guard let rejection = report.rejected, !rejection.committed, rejection.statusFlags != 0,
                      report.committedSteps < stage.steps,
                      rejection.systemFingerprint == exit.systemFingerprint,
                      rejection.configurationFingerprint == exit.configurationFingerprint,
                      rejection.stepIndex == exit.acceptedStep, rejection.timeBeforePS == exit.timePS,
                      rejection.timeAfterPS == exit.timePS else { throw invalid("rejection certificate binding") }
            }
        }
    }

    /// A publication failure may leave only the CURRENT accepted sample or
    /// observation pending. That state can be reconstructed without replaying
    /// dynamics; an older hole cannot and must never be silently skipped.
    private static func validateOutputs(stage: VivoMDProtocolStage, completed: UInt64,
        entry: VivoMDCheckpoint, current: VivoMDCheckpoint, manifest hash: VivoFingerprint?,
        observationTail: VivoFingerprint?, observationCount: UInt64, allowPending: Bool,
        system: VivoClassicalSystem, plan: VivoMDProtocolPlan, store: VivoArtifactStore) async throws
        -> (sample: UInt64?, observation: UInt64?, trajectory: VivoMDTrajectoryValidation?) {
        func expectedLast(count: UInt64, interval: UInt64?) throws -> UInt64? {
            guard let interval else {
                guard count == 0 else { throw invalid("output for a disabled schedule") }; return nil
            }
            let expected = completed / interval
            let pendingCurrent = allowPending && completed > 0 && completed % interval == 0 && count == expected - 1
            guard count == expected || pendingCurrent else { throw invalid("output schedule has an unrecoverable gap") }
            return count == 0 ? nil : entry.acceptedStep + count * interval
        }
        var sample: UInt64?
        var trajectoryValidation: VivoMDTrajectoryValidation?
        if let hash {
            guard let interval = stage.sampleEvery else { throw invalid("trajectory for an unsampled stage") }
            let manifest = try await read(VivoMDTrajectoryManifest.self, hash, "md-trajectory-manifest", store, maximumBytes: 64 * 1024)
            try manifest.validate()
            sample = try expectedLast(count: manifest.frameCount, interval: interval)
            guard manifest.systemFingerprint == entry.systemFingerprint,
                  manifest.configurationFingerprint == entry.configurationFingerprint,
                  manifest.particleCount == UInt32(system.particles.count),
                  manifest.includeVelocities == plan.trajectoryIncludesVelocities,
                  manifest.sealed == !allowPending, manifest.lastStep == sample,
                  manifest.firstStep == (sample == nil ? nil : entry.acceptedStep + interval),
                  manifest.firstTimePS.map({ $0 > entry.timePS }) ?? true,
                  manifest.lastTimePS.map({ $0 <= current.timePS }) ?? true else { throw invalid("trajectory prefix binding") }
            let reader = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: hash)
            trajectoryValidation = try await reader.validate(scope: .restart)
        } else if stage.sampleEvery != nil { throw invalid("sampled stage has no durable trajectory prefix") }
        let observation = try expectedLast(count: observationCount, interval: stage.observablesEvery)
        guard (observation == nil) == (observationTail == nil) else { throw invalid("observation count/tail") }
        if let tail = observationTail {
            let link = try await read(VivoMDObservationLink.self, tail, "md-observation-link", store)
            let value = link.observation
            let degrees = UInt64(system.particles.filter { $0.massDa > 0 }.count) * 3
            guard link.schema == VivoMDObservationLink.schemaID, link.stageIdentifier == stage.identifier,
                  link.ordinal == observationCount - 1, (link.previous == nil) == (observationCount == 1),
                  value.systemFingerprint == entry.systemFingerprint,
                  value.configurationFingerprint == entry.configurationFingerprint,
                  value.stepIndex == observation, value.timePS > entry.timePS, value.timePS <= current.timePS,
                  value.potentialEnergyKJPerMol.isFinite, value.kineticEnergyKJPerMol.isFinite, value.kineticEnergyKJPerMol >= 0,
                  value.totalEnergyKJPerMol.isFinite, value.totalEnergyKJPerMol == value.potentialEnergyKJPerMol + value.kineticEnergyKJPerMol,
                  degrees > UInt64(system.constraints.count), value.degreesOfFreedom == degrees - UInt64(system.constraints.count),
                  let temperature = value.temperatureK, temperature.isFinite, temperature >= 0 else { throw invalid("observation prefix binding") }
        }
        return (sample, observation, trajectoryValidation)
    }

    private static func read<T: Decodable & Sendable>(_ type: T.Type, _ hash: VivoFingerprint,
        _ kind: String, _ store: VivoArtifactStore, maximumBytes: UInt64 = 2 * 1024 * 1024) async throws -> T {
        let descriptor = try await store.descriptor(for: hash)
        guard descriptor.kind == kind, descriptor.byteCount <= maximumBytes else { throw invalid("artifact type or size") }
        let data = try await store.data(for: hash, maximumBytes: Int(maximumBytes))
        guard UInt64(data.count) == descriptor.byteCount else { throw invalid("artifact descriptor byte count") }
        return try VivoCanonicalJSON.decode(type, from: data)
    }

    private static func invalid(_ detail: String) -> VivoArtifactValidationError {
        .incompatible("MD protocol restart: " + detail)
    }
}
