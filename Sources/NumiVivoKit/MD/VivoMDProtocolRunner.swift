import Foundation
@preconcurrency import Metal

/// A run owns its runtime exclusively; callers receive immutable artifact
/// identities rather than an arena they can mutate between protocol stages.
public actor VivoMDProtocolRunner {
    public nonisolated let runID: UUID
    public nonisolated let checkpointReferenceName: String
    public nonisolated let plan: VivoMDProtocolPlan
    public nonisolated let planFingerprint: VivoFingerprint
    private let system: VivoClassicalSystem
    private let store: VivoArtifactStore
    private let device: MTLDevice
    private var cursor: VivoMDProtocolCheckpoint
    private var runtime: VivoMDMetalRuntime?
    private var trajectory: VivoMDTrajectoryArchiveWriter?
    private var latestDurable: VivoStoredArtifact?
    private var busy = false

    public static func start(system: VivoClassicalSystem, initialState: VivoClassicalInitialState,
                             plan: VivoMDProtocolPlan, store: VivoArtifactStore,
                             device: MTLDevice? = nil) async throws -> VivoMDProtocolRunner {
        let device = try device ?? VivoMetalDeviceSelector.productionDevice()
        let capabilities = try plan.validate(system: system, initialState: initialState)
        let blockers = capabilities.flatMap(\.blockers)
        guard blockers.isEmpty else { throw VivoMDRuntimeError.unsupported(Array(Set(blockers)).sorted()) }
        try Task.checkCancellation()
        let storedPlan = try await put(plan, kind: "md-protocol", store: store)
        let initial = try await put(initialState, kind: "classical-initial-state", store: store)
        let stage = plan.stages[0]
        let runtime = try await VivoMDMetalRuntime.make(system: system, initialState: initialState,
                                                        configuration: stage.configuration, device: device)
        if stage.velocityInitialization == .maxwellBoltzmann {
            _ = try await runtime.thermalize(temperatureK: stage.configuration.targetTemperatureK!,
                                             seed: stage.thermalizationSeed!)
        }
        let entry = try await put(runtime.checkpoint(), kind: "md-checkpoint", store: store)
        let writer = try makeWriter(stage: stage, plan: plan, store: store, particleCount: system.particles.count)
        let runID = UUID()
        let cursor = VivoMDProtocolCheckpoint(schema: VivoMDProtocolCheckpoint.schemaID,
            numericalContract: VivoMDExecutionIdentity.current, runID: runID, resumedFrom: nil,
            planFingerprint: storedPlan.fingerprint, systemFingerprint: plan.systemFingerprint,
            initialStateArtifact: initial.fingerprint, stageIndex: 0, phase: .running, completedStepsInStage: 0,
            entryCheckpoint: entry.fingerprint, currentCheckpoint: entry.fingerprint, transition: nil,
            trajectoryManifest: nil, observationTail: nil, observationCount: 0,
            priorStageReports: [], activeStageReport: nil)
        let runner = VivoMDProtocolRunner(system: system, plan: plan, store: store, device: device,
            cursor: cursor, runtime: runtime, trajectory: writer, latestDurable: nil)
        try await runner.persist()
        return runner
    }

    /// Resume forks a new named run from a verified immutable prefix. It never
    /// changes the original run reference or reinitializes a running stage.
    public static func resume(system: VivoClassicalSystem, plan: VivoMDProtocolPlan,
                              store: VivoArtifactStore, checkpoint: VivoFingerprint,
                              device: MTLDevice? = nil) async throws -> VivoMDProtocolRunner {
        let device = try device ?? VivoMetalDeviceSelector.productionDevice()
        try plan.validate()
        guard plan.systemFingerprint == (try system.fingerprint()) else {
            throw VivoArtifactValidationError.incompatible("resume system differs from protocol")
        }
        var cursor = try await read(VivoMDProtocolCheckpoint.self, fingerprint: checkpoint,
                                    kind: "md-protocol-checkpoint", store: store)
        try cursor.validate(plan: plan)
        let current = try await read(VivoMDCheckpoint.self, fingerprint: cursor.currentCheckpoint,
                                     kind: "md-checkpoint", store: store, maximumBytes: 512 * 1024 * 1024)
        let entry = try await read(VivoMDCheckpoint.self, fingerprint: cursor.entryCheckpoint,
                                   kind: "md-checkpoint", store: store, maximumBytes: 512 * 1024 * 1024)
        try current.validate(particleCount: system.particles.count)
        try entry.validate(particleCount: system.particles.count)
        let stage = plan.stages[cursor.stageIndex]
        let configurationID = try stage.configuration.fingerprint()
        guard current.systemFingerprint == plan.systemFingerprint, entry.systemFingerprint == plan.systemFingerprint,
              current.configurationFingerprint == configurationID, entry.configurationFingerprint == configurationID,
              cursor.completedStepsInStage <= UInt64.max - entry.acceptedStep,
              current.acceptedStep == entry.acceptedStep + cursor.completedStepsInStage,
              current.timePS >= entry.timePS else {
            throw VivoArtifactValidationError.incompatible("protocol cursor and MD checkpoint disagree")
        }
        if stage.kind == .minimization, current.timePS != entry.timePS {
            throw VivoArtifactValidationError.invalid("minimization cursor advanced physical time")
        }
        for (index, hash) in cursor.priorStageReports.enumerated() {
            let report = try await read(VivoMDProtocolStageReport.self, fingerprint: hash,
                                        kind: "md-stage-report", store: store)
            guard report.planFingerprint == cursor.planFingerprint, report.stageIndex == index,
                  report.stageIdentifier == plan.stages[index].identifier, report.successful else {
                throw VivoArtifactValidationError.incompatible("completed-stage report disagrees with protocol")
            }
        }
        if let hash = cursor.activeStageReport {
            let report = try await read(VivoMDProtocolStageReport.self, fingerprint: hash,
                                        kind: "md-stage-report", store: store)
            guard report.planFingerprint == cursor.planFingerprint, report.stageIndex == cursor.stageIndex,
                  report.stageIdentifier == stage.identifier, report.exitCheckpoint == cursor.currentCheckpoint,
                  report.committedSteps == cursor.completedStepsInStage,
                  report.successful == (cursor.phase == .stageFinished) else {
                throw VivoArtifactValidationError.incompatible("active-stage report disagrees with cursor")
            }
        }
        if let tail = cursor.observationTail {
            let observation = try await read(VivoMDObservationLink.self, fingerprint: tail,
                                             kind: "md-observation-link", store: store)
            guard observation.schema == VivoMDObservationLink.schemaID,
                  observation.ordinal == cursor.observationCount - 1,
                  observation.stageIdentifier == stage.identifier,
                  observation.observation.systemFingerprint == plan.systemFingerprint,
                  observation.observation.configurationFingerprint == configurationID,
                  observation.observation.stepIndex <= current.acceptedStep,
                  observation.observation.timePS <= current.timePS else {
                throw VivoArtifactValidationError.incompatible("observation prefix disagrees with restart state")
            }
        }
        var writer: VivoMDTrajectoryArchiveWriter?
        if let hash = cursor.trajectoryManifest {
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: hash)
            let manifest = archive.manifest
            guard manifest.systemFingerprint == plan.systemFingerprint,
                  manifest.configurationFingerprint == configurationID,
                  manifest.particleCount == UInt32(system.particles.count),
                  manifest.includeVelocities == plan.trajectoryIncludesVelocities,
                  manifest.lastStep.map({ $0 <= current.acceptedStep }) ?? true,
                  manifest.lastTimePS.map({ $0 <= current.timePS }) ?? true,
                  manifest.sealed == (cursor.phase != .running) else {
                throw VivoArtifactValidationError.incompatible("trajectory prefix disagrees with restart cursor")
            }
            if cursor.phase == .running {
                writer = try await .resume(store: store, manifest: hash, targetChunkBytes: plan.trajectoryChunkBytes)
            }
        } else if stage.sampleEvery != nil {
            throw VivoArtifactValidationError.invalid("sampled stage cursor has no durable trajectory prefix")
        }
        let runtime: VivoMDMetalRuntime?
        if cursor.phase == .running {
            runtime = try await .restore(system: system, configuration: stage.configuration,
                                          checkpoint: current, device: device)
        } else { runtime = nil }
        cursor.runID = UUID(); cursor.resumedFrom = checkpoint
        let stored = try await store.descriptor(for: checkpoint)
        let runner = VivoMDProtocolRunner(system: system, plan: plan, store: store, device: device,
            cursor: cursor, runtime: runtime, trajectory: writer, latestDurable: stored)
        try await runner.persist()
        return runner
    }

    private init(system: VivoClassicalSystem, plan: VivoMDProtocolPlan, store: VivoArtifactStore,
                 device: MTLDevice, cursor: VivoMDProtocolCheckpoint, runtime: VivoMDMetalRuntime?,
                 trajectory: VivoMDTrajectoryArchiveWriter?, latestDurable: VivoStoredArtifact?) {
        self.system = system; self.plan = plan; self.store = store; self.device = device
        self.cursor = cursor; self.runtime = runtime; self.trajectory = trajectory; self.latestDurable = latestDurable
        runID = cursor.runID; planFingerprint = cursor.planFingerprint
        checkpointReferenceName = "md-\(cursor.runID.uuidString.lowercased())-checkpoint"
    }

    public func run() async throws -> VivoMDProtocolRunReceipt {
        guard !busy else { throw VivoArtifactValidationError.invalid("MD protocol runner is already active") }
        busy = true; defer { busy = false }
        do {
            while true {
                try Task.checkCancellation()
                if cursor.phase == .blocked { return receipt(.rejected, diagnostic: "stage gate blocked; inspect its immutable report") }
                if cursor.phase == .stageFinished {
                    if cursor.stageIndex + 1 == plan.stages.count { return receipt(.completed) }
                    try await advance()
                    continue
                }
                guard let runtime else { throw VivoArtifactValidationError.invalid("running protocol has no MD authority") }
                let stage = plan.stages[cursor.stageIndex]
                if stage.kind == .minimization {
                    let result = try await runtime.minimize(stage.minimization!)
                    let success = result.converged || !stage.requireMinimizationConvergence
                    try await finishStage(success: success, minimization: result, rejected: nil)
                    continue
                }
                while cursor.completedStepsInStage < stage.steps {
                    try Task.checkCancellation()
                    let result = try await runtime.step()
                    guard result.committed else {
                        try await finishStage(success: false, minimization: nil, rejected: result)
                        return receipt(.rejected, diagnostic: "MD candidate rejected; no automatic retry or timestep change")
                    }
                    cursor.completedStepsInStage += 1
                    let ordinal = cursor.completedStepsInStage
                    let sampleDue = stage.sampleEvery.map { ordinal % $0 == 0 } ?? false
                    let observationDue = stage.observablesEvery.map { ordinal % $0 == 0 } ?? false
                    if sampleDue {
                        let pair = try await runtime.sample(includeObservables: observationDue)
                        guard let trajectory else { throw VivoArtifactValidationError.invalid("trajectory writer missing") }
                        try await trajectory.append(pair.state)
                        if let observation = pair.observables { try await record(observation) }
                    } else if observationDue {
                        try await record(runtime.observables())
                    }
                    if ordinal % stage.checkpointEvery == 0 { try await persist() }
                }
                try await finishStage(success: true, minimization: nil, rejected: nil)
            }
        } catch {
            let disposition: VivoMDProtocolDisposition = error is CancellationError ? .cancelled : .failed
            var persistence: String?
            // A healthy runtime can export its last accepted boundary even after
            // task cancellation. A poisoned GPU runtime keeps the prior durable
            // cursor; never fabricate a checkpoint of the failed command.
            do { try await persist() } catch { persistence = String(describing: error) }
            return receipt(disposition, diagnostic: String(describing: error), persistence: persistence)
        }
    }

    private func advance() async throws {
        guard cursor.phase == .stageFinished, let previousReport = cursor.activeStageReport else {
            throw VivoArtifactValidationError.invalid("cannot transition an unfinished MD stage")
        }
        let nextIndex = cursor.stageIndex + 1
        let previous = plan.stages[cursor.stageIndex], next = plan.stages[nextIndex]
        let state = try await Self.read(VivoMDCheckpoint.self, fingerprint: cursor.currentCheckpoint,
                                        kind: "md-checkpoint", store: store, maximumBytes: 512 * 1024 * 1024)
        let transfer = try VivoMDStageTransfer.prepare(checkpoint: state, source: previous.configuration,
            destination: next.configuration, particleCount: system.particles.count,
            velocityInitialization: next.velocityInitialization, thermalizationSeed: next.thermalizationSeed)
        let transition = try await Self.put(transfer, kind: "md-stage-transition", store: store)
        // Release the old arena before creating the next configuration. The
        // durable checkpoint is the authority during this fallible transition.
        runtime = nil
        let nextRuntime = try await VivoMDMetalRuntime.restore(system: system, configuration: next.configuration,
            checkpoint: transfer.destinationCheckpoint, device: device)
        if next.velocityInitialization == .maxwellBoltzmann {
            _ = try await nextRuntime.thermalize(temperatureK: next.configuration.targetTemperatureK!, seed: next.thermalizationSeed!)
        }
        let entry = try await Self.put(nextRuntime.checkpoint(), kind: "md-checkpoint", store: store)
        let writer = try Self.makeWriter(stage: next, plan: plan, store: store, particleCount: system.particles.count)
        cursor.priorStageReports.append(previousReport)
        cursor.stageIndex = nextIndex; cursor.phase = .running; cursor.completedStepsInStage = 0
        cursor.entryCheckpoint = entry.fingerprint; cursor.currentCheckpoint = entry.fingerprint
        cursor.transition = transition.fingerprint; cursor.trajectoryManifest = nil
        cursor.observationTail = nil; cursor.observationCount = 0; cursor.activeStageReport = nil
        runtime = nextRuntime; trajectory = writer
        // Resume from this entry skips initialization; it is never sampled twice.
        try await persist()
    }

    private func record(_ observation: VivoMDObservables) async throws {
        guard cursor.observationCount < .max, observation.systemFingerprint == cursor.systemFingerprint,
              observation.configurationFingerprint == (try plan.stages[cursor.stageIndex].configuration.fingerprint()),
              observation.potentialEnergyKJPerMol.isFinite, observation.kineticEnergyKJPerMol.isFinite,
              observation.timePS.isFinite else { throw VivoArtifactValidationError.invalid("invalid protocol observation") }
        let link = VivoMDObservationLink(schema: VivoMDObservationLink.schemaID,
            stageIdentifier: plan.stages[cursor.stageIndex].identifier, ordinal: cursor.observationCount,
            previous: cursor.observationTail, observation: observation)
        let artifact = try await Self.put(link, kind: "md-observation-link", store: store)
        cursor.observationTail = artifact.fingerprint; cursor.observationCount += 1
    }

    private func finishStage(success: Bool, minimization: VivoMDMinimizationCertificate?,
                             rejected: VivoMDStepCertificate?) async throws {
        guard let runtime else { throw VivoArtifactValidationError.invalid("cannot finish without an MD state") }
        let state = try await runtime.checkpoint()
        let exit = try await Self.put(state, kind: "md-checkpoint", store: store)
        let entry = try await Self.read(VivoMDCheckpoint.self, fingerprint: cursor.entryCheckpoint,
                                        kind: "md-checkpoint", store: store, maximumBytes: 512 * 1024 * 1024)
        var sealedHash: VivoFingerprint?
        if let trajectory {
            // Build an immutable sealed view without sealing/mutating the writer
            // first. A later failed report write can still flush a valid running
            // prefix during error handling; no half-finalized cursor is emitted.
            let prefix = try await trajectory.snapshot()
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: prefix.fingerprint)
            let m = archive.manifest
            let sealed = VivoMDTrajectoryManifest(schema: m.schema, systemFingerprint: m.systemFingerprint,
                configurationFingerprint: m.configurationFingerprint, particleCount: m.particleCount,
                includeVelocities: m.includeVelocities, frameCount: m.frameCount, chunkCount: m.chunkCount,
                tail: m.tail, firstStep: m.firstStep, lastStep: m.lastStep,
                firstTimePS: m.firstTimePS, lastTimePS: m.lastTimePS, sealed: true)
            let artifact = try await store.put(data: VivoCanonicalJSON.encode(sealed), kind: "md-trajectory-manifest",
                                               mediaType: "application/vnd.numivivo.md-trajectory-manifest+json")
            sealedHash = artifact.fingerprint
        }
        let report = VivoMDProtocolStageReport(schema: "numivivo.org/md-stage-report/v1",
            planFingerprint: planFingerprint, stageIdentifier: plan.stages[cursor.stageIndex].identifier,
            stageIndex: cursor.stageIndex, successful: success, entryCheckpoint: cursor.entryCheckpoint,
            exitCheckpoint: exit.fingerprint, transition: cursor.transition, committedSteps: cursor.completedStepsInStage,
            startTimePS: entry.timePS, endTimePS: state.timePS, trajectoryManifest: sealedHash,
            observationTail: cursor.observationTail, observationCount: cursor.observationCount,
            minimization: minimization, rejected: rejected)
        let stored = try await Self.put(report, kind: "md-stage-report", store: store)
        cursor.currentCheckpoint = exit.fingerprint; cursor.activeStageReport = stored.fingerprint
        cursor.trajectoryManifest = sealedHash; cursor.phase = success ? .stageFinished : .blocked
        trajectory = nil
        try await persist()
    }

    private func persist() async throws {
        if let runtime {
            let checkpoint = try await runtime.checkpoint()
            let artifact = try await Self.put(checkpoint, kind: "md-checkpoint", store: store)
            cursor.currentCheckpoint = artifact.fingerprint
        }
        if let trajectory {
            let prefix = try await trajectory.snapshot()
            cursor.trajectoryManifest = prefix.fingerprint
        }
        try cursor.validate(plan: plan)
        let artifact = try await Self.put(cursor, kind: "md-protocol-checkpoint", store: store)
        latestDurable = artifact
        _ = try await store.setReference(checkpointReferenceName, to: artifact)
    }
    private func receipt(_ disposition: VivoMDProtocolDisposition, diagnostic: String? = nil,
                         persistence: String? = nil) -> VivoMDProtocolRunReceipt {
        .init(schema: "numivivo.org/md-protocol-receipt/v1", runID: runID, planFingerprint: planFingerprint,
              disposition: disposition, completedStages: cursor.stageIndex + (cursor.phase == .stageFinished ? 1 : 0),
              latestDurableCheckpoint: latestDurable, checkpointReference: checkpointReferenceName,
              diagnostic: diagnostic, persistenceDiagnostic: persistence)
    }
    private static func makeWriter(stage: VivoMDProtocolStage, plan: VivoMDProtocolPlan,
                                   store: VivoArtifactStore, particleCount: Int) throws -> VivoMDTrajectoryArchiveWriter? {
        guard stage.sampleEvery != nil else { return nil }
        return try .init(store: store, systemFingerprint: plan.systemFingerprint,
                         configurationFingerprint: stage.configuration.fingerprint(), particleCount: UInt32(particleCount),
                         includeVelocities: plan.trajectoryIncludesVelocities, targetChunkBytes: plan.trajectoryChunkBytes)
    }
    private static func put<T: Encodable & Sendable>(_ value: T, kind: String,
                                                    store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind,
                            mediaType: "application/vnd.numivivo.\(kind)+json")
    }
    private static func read<T: Decodable & Sendable>(_ type: T.Type, fingerprint: VivoFingerprint, kind: String,
                                                     store: VivoArtifactStore, maximumBytes: UInt64 = 2 * 1024 * 1024) async throws -> T {
        let descriptor = try await store.descriptor(for: fingerprint)
        guard descriptor.kind == kind, descriptor.byteCount <= maximumBytes else {
            throw VivoArtifactValidationError.incompatible("unexpected type/size for protocol artifact \(fingerprint.hex)")
        }
        return try VivoCanonicalJSON.decode(type, from: await store.data(for: fingerprint))
    }
}
