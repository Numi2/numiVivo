import Foundation
import Testing
@testable import NumiVivoKit

@Suite(.serialized) struct MDProtocolResumeTests {
    private struct Fixture {
        let root: URL
        let store: VivoArtifactStore
        let system: VivoClassicalSystem
        let initial: VivoClassicalInitialState
        var plan: VivoMDProtocolPlan
        let entry: VivoMDCheckpoint
        var cursor: VivoMDProtocolCheckpoint
    }

    private func configuration() -> VivoMDConfiguration {
        .init(timeStepPS: 0.0001, electrostatics: .cutoff, ensemble: .nve,
              thermostat: .none, targetTemperatureK: nil, frictionPerPS: nil, neighborListEnabled: false)
    }
    private func put<T: Encodable & Sendable>(_ value: T, _ kind: String,
                                              _ store: VivoArtifactStore) async throws -> VivoFingerprint {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind,
                            mediaType: "application/vnd.numivivo.\(kind)+json").fingerprint
    }
    private func read<T: Decodable & Sendable>(_ type: T.Type, _ hash: VivoFingerprint,
                                               _ store: VivoArtifactStore) async throws -> T {
        try VivoCanonicalJSON.decode(type, from: await store.data(for: hash))
    }
    private func fixture(minimization: Bool = false, twoStages: Bool = false, dynamicsSteps: UInt64 = 4) async throws -> Fixture {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-md-resume-\(UUID().uuidString)")
        let store = try VivoArtifactStore(rootURL: root)
        let identity = try VivoCanonicalJSON.fingerprint(Data("MD protocol resume regression".utf8))
        let system = VivoClassicalSystem(identifier: "resume-harmonic-pair", structureFingerprint: identity,
            particles: (0..<2).map { .init(index: UInt32($0), atomIndex: UInt32($0), typeIdentifier: "H",
                                          massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0) },
            bonds: [.init(a: 0, b: 1, lengthNM: 0.1, forceConstant: 1000)])
        let initial = VivoClassicalInitialState(systemFingerprint: try system.fingerprint(),
            positionsNM: [.init(-0.0625, 0, 0), .init(0.0625, 0, 0)])
        let stage = VivoMDProtocolStage(identifier: "prepare", kind: minimization ? .minimization : .dynamics,
            configuration: configuration(), steps: minimization ? 0 : dynamicsSteps,
            minimization: minimization ? .init(maximumIterations: 8, forceToleranceKJPerMolNM: 1) : nil,
            velocityInitialization: .zero, sampleEvery: minimization ? nil : 1,
            observablesEvery: minimization ? nil : 1, checkpointEvery: 1)
        var stages = [stage]
        if twoStages {
            stages.append(.init(identifier: "production", kind: .dynamics, configuration: configuration(),
                steps: 4, velocityInitialization: .preserve, sampleEvery: 1, observablesEvery: 1, checkpointEvery: 1))
        }
        let plan = VivoMDProtocolPlan(identifier: "resume-regression", systemFingerprint: try system.fingerprint(),
                                      stages: stages, trajectoryIncludesVelocities: true, trajectoryChunkBytes: 512)
        let entry = VivoMDCheckpoint(systemFingerprint: try system.fingerprint(), configurationFingerprint: try configuration().fingerprint(),
            acceptedStep: 0, timePS: 0, positionsNM: initial.positionsNM,
            velocitiesNMPerPS: [.zero, .zero], periodicCell: nil)
        let entryHash = try await put(entry, "md-checkpoint", store)
        let planHash = try await put(plan, "md-protocol", store)
        let initialHash = try await put(initial, "classical-initial-state", store)
        let manifest: VivoFingerprint?
        if minimization { manifest = nil }
        else {
            let writer = try VivoMDTrajectoryArchiveWriter(store: store, systemFingerprint: entry.systemFingerprint,
                configurationFingerprint: entry.configurationFingerprint, particleCount: 2, includeVelocities: true,
                targetChunkBytes: plan.trajectoryChunkBytes)
            manifest = try await writer.snapshot().fingerprint
        }
        let cursor = VivoMDProtocolCheckpoint(schema: VivoMDProtocolCheckpoint.schemaID,
            numericalContract: VivoMDExecutionIdentity.current, runID: UUID(), resumedFrom: nil,
            planFingerprint: planHash, systemFingerprint: plan.systemFingerprint, initialStateArtifact: initialHash,
            stageIndex: 0, phase: .running, completedStepsInStage: 0,
            entryCheckpoint: entryHash, currentCheckpoint: entryHash, transition: nil, trajectoryManifest: manifest,
            observationTail: nil, observationCount: 0, priorStageReports: [], activeStageReport: nil)
        return .init(root: root, store: store, system: system, initial: initial, plan: plan, entry: entry, cursor: cursor)
    }
    private func load(_ fixture: Fixture) async throws -> VivoMDProtocolResumeValidation.Prefix {
        let hash = try await put(fixture.cursor, "md-protocol-checkpoint", fixture.store)
        return try await VivoMDProtocolResumeValidation.load(system: fixture.system, plan: fixture.plan,
                                                              store: fixture.store, checkpoint: hash)
    }
    private func rejects(_ operation: () async throws -> Void) async {
        do { try await operation(); Issue.record("Expected an inconsistent MD protocol prefix to be rejected") }
        catch { }
    }
    private func minimizationReport(_ f: Fixture, converged: Bool = true, successful: Bool = true) -> VivoMDProtocolStageReport {
        .init(schema: "numivivo.org/md-stage-report/v1", planFingerprint: f.cursor.planFingerprint,
            stageIdentifier: f.plan.stages[0].identifier, stageIndex: 0, successful: successful,
            entryCheckpoint: f.cursor.entryCheckpoint, exitCheckpoint: f.cursor.currentCheckpoint,
            transition: nil, committedSteps: 0, startTimePS: 0, endTimePS: 0, trajectoryManifest: nil,
            observationTail: nil, observationCount: 0,
            minimization: .init(systemFingerprint: f.entry.systemFingerprint, configurationFingerprint: f.entry.configurationFingerprint,
                converged: converged, attemptedIterations: 0, acceptedIterations: 0, rejectedIterations: 0,
                initialPotentialEnergyKJPerMol: 1, finalPotentialEnergyKJPerMol: 1,
                finalMaximumForceKJPerMolNM: converged ? 0 : 2, finalStepScale: 1e-5), rejected: nil)
    }
    private func altered<T: Codable>(_ value: T, _ change: (inout [String: Any]) -> Void) throws -> T {
        var object = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(value)) as? [String: Any])
        change(&object)
        return try VivoCanonicalJSON.decode(T.self, from: JSONSerialization.data(withJSONObject: object, options: .sortedKeys))
    }

    @Test func requiredMinimizationGateCannotBeBypassedByARehashedReport() async throws {
        var f = try await fixture(minimization: true); defer { try? FileManager.default.removeItem(at: f.root) }
        f.cursor.phase = .stageFinished
        let valid = minimizationReport(f)
        f.cursor.activeStageReport = try await put(valid, "md-stage-report", f.store)
        _ = try await load(f)

        let invalidReports = [
            minimizationReport(f, converged: false, successful: true),
            try altered(valid) { $0.removeValue(forKey: "minimization") },
            try altered(valid) { $0["schema"] = "unrelated-schema" },
            try altered(valid) { $0["observationCount"] = 1 },
            try altered(valid) { $0["startTimePS"] = 1 }
        ]
        for report in invalidReports {
            f.cursor.activeStageReport = try await put(report, "md-stage-report", f.store)
            await rejects { _ = try await load(f) }
        }
        // A legitimate blocked prefix remains inspectable and cannot advance.
        f.cursor.phase = .blocked
        f.cursor.activeStageReport = try await put(minimizationReport(f, converged: false, successful: false), "md-stage-report", f.store)
        #expect(try await load(f).cursor.phase == .blocked)
    }

    @Test func priorStageAndTransitionMustBindTheActualEntryState() async throws {
        var f = try await fixture(minimization: true, twoStages: true)
        defer { try? FileManager.default.removeItem(at: f.root) }
        let report = minimizationReport(f)
        f.cursor.priorStageReports = [try await put(report, "md-stage-report", f.store)]
        f.cursor.stageIndex = 1
        let transfer = try VivoMDStageTransfer.prepare(checkpoint: f.entry, source: f.plan.stages[0].configuration,
            destination: f.plan.stages[1].configuration, particleCount: 2)
        f.cursor.transition = try await put(transfer, "md-stage-transition", f.store)
        let writer = try VivoMDTrajectoryArchiveWriter(store: f.store, systemFingerprint: f.entry.systemFingerprint,
            configurationFingerprint: f.entry.configurationFingerprint, particleCount: 2, includeVelocities: true)
        f.cursor.trajectoryManifest = try await writer.snapshot().fingerprint
        _ = try await load(f)
        let validTransition = f.cursor.transition
        f.cursor.transition = nil
        await rejects { _ = try await load(f) }
        f.cursor.transition = try await put(try altered(transfer) { $0["velocityInitialization"] = "zero" }, "md-stage-transition", f.store)
        await rejects { _ = try await load(f) }
        f.cursor.transition = validTransition
        var moved = f.entry; moved.positionsNM[0].x = -0.125
        let movedHash = try await put(moved, "md-checkpoint", f.store)
        f.cursor.entryCheckpoint = movedHash; f.cursor.currentCheckpoint = movedHash
        await rejects { _ = try await load(f) }
    }

    @Test func blockedMinimizerPreservesRecomputedEnergyAndUnderflowTermination() async throws {
        var f = try await fixture(minimization: true)
        defer { try? FileManager.default.removeItem(at: f.root) }
        f.cursor.phase = .blocked
        var report = minimizationReport(f, converged: false, successful: false)
        report = try altered(report) {
            var certificate = $0["minimization"] as! [String: Any]
            certificate["attemptedIterations"] = 1; certificate["rejectedIterations"] = 1
            certificate["finalStepScale"] = 5e-6
            certificate["finalPotentialEnergyKJPerMol"] = Double(Float(1).nextUp)
            $0["minimization"] = certificate
        }
        // After rejection, minimize() re-evaluates the unchanged accepted
        // geometry. PME charge spreading uses unordered FP32 atomic sums, so
        // resume must not invent exact energy monotonicity for that readback.
        // This is a certificate contract test, not a measured PME-noise result.
        f.cursor.activeStageReport = try await put(report, "md-stage-report", f.store)
        #expect(try await load(f).cursor.phase == .blocked)

        f.plan.stages[0].minimization?.initialStepScale = .leastNonzeroMagnitude
        f.plan.stages[0].minimization?.minimumStepScale = .leastNonzeroMagnitude
        f.cursor.planFingerprint = try await put(f.plan, "md-protocol", f.store)
        let settings = try #require(f.plan.stages[0].minimization)
        try settings.validate()
        // These are precisely the Double operands in the rejected-step branch.
        let exhaustedScale = settings.initialStepScale * settings.rejectedStepShrink
        #expect(exhaustedScale == 0)
        report = try altered(minimizationReport(f, converged: false, successful: false)) {
            var certificate = $0["minimization"] as! [String: Any]
            certificate["attemptedIterations"] = 1; certificate["rejectedIterations"] = 1
            certificate["finalStepScale"] = exhaustedScale
            $0["minimization"] = certificate
        }
        f.cursor.activeStageReport = try await put(report, "md-stage-report", f.store)
        #expect(try await load(f).cursor.phase == .blocked)
        let impossible = try altered(report) {
            var certificate = $0["minimization"] as! [String: Any]
            certificate["attemptedIterations"] = 0; certificate["rejectedIterations"] = 0
            $0["minimization"] = certificate
        }
        f.cursor.activeStageReport = try await put(impossible, "md-stage-report", f.store)
        await rejects { _ = try await load(f) }
    }

    @Test func onlyTheCurrentAcceptedPublicationMayBePending() async throws {
        var f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        var current = f.entry; current.acceptedStep = 1; current.timePS = configuration().timeStepPS
        f.cursor.completedStepsInStage = 1
        f.cursor.currentCheckpoint = try await put(current, "md-checkpoint", f.store)
        let pending = try await load(f)
        #expect(pending.lastSampledStep == nil && pending.lastObservedStep == nil)
        current.timePS = 1
        f.cursor.currentCheckpoint = try await put(current, "md-checkpoint", f.store)
        await rejects { _ = try await load(f) }
        current.timePS = configuration().timeStepPS
        current.acceptedStep = 2; current.timePS += configuration().timeStepPS
        f.cursor.completedStepsInStage = 2
        f.cursor.currentCheckpoint = try await put(current, "md-checkpoint", f.store)
        await rejects { _ = try await load(f) }
    }

    @Test func clockValidationAcceptsRepeatedFP64AccumulationAndRejectsARehashedTimeChange() async throws {
        var f = try await fixture(dynamicsSteps: 1_000_000)
        defer { try? FileManager.default.removeItem(at: f.root) }
        f.plan.stages[0].sampleEvery = nil; f.plan.stages[0].observablesEvery = nil
        f.cursor.planFingerprint = try await put(f.plan, "md-protocol", f.store)
        f.cursor.trajectoryManifest = nil
        var current = f.entry
        let timestep = f.plan.stages[0].configuration.timeStepPS
        for _ in 0..<f.plan.stages[0].steps { current.timePS += timestep }
        current.acceptedStep = f.plan.stages[0].steps
        #expect(current.timePS != Double(current.acceptedStep) * timestep)
        f.cursor.completedStepsInStage = current.acceptedStep
        f.cursor.currentCheckpoint = try await put(current, "md-checkpoint", f.store)
        _ = try await load(f)
        current.timePS += 0.01
        f.cursor.currentCheckpoint = try await put(current, "md-checkpoint", f.store)
        await rejects { _ = try await load(f) }
    }

    @Test func nativeMinimizationAndThermalStageRetainTheirVerifiedTransitions() async throws {
        let f = try await fixture(minimization: true, twoStages: true)
        defer { try? FileManager.default.removeItem(at: f.root) }
        var plan = f.plan
        // A declared loose threshold is enough to exercise the actual native
        // minimizer gate without claiming a physically equilibrated structure.
        plan.stages[0].minimization?.forceToleranceKJPerMolNM = 100
        plan.stages[1].configuration = .init(timeStepPS: 0.0001, electrostatics: .cutoff,
                                            neighborListEnabled: false, randomSeed: 19)
        plan.stages[1].velocityInitialization = .maxwellBoltzmann
        plan.stages[1].thermalizationSeed = 73
        let runner = try await VivoMDProtocolRunner.start(system: f.system, initialState: f.initial, plan: plan, store: f.store)
        let receipt = try await runner.run()
        #expect(receipt.disposition == .completed && receipt.completedStages == 2)
        let durable = try #require(receipt.latestDurableCheckpoint)
        let verified = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: plan,
                                                                    store: f.store, checkpoint: durable.fingerprint)
        #expect(verified.cursor.priorStageReports.count == 1 && verified.cursor.stageIndex == 1)
        let prior = try await read(VivoMDProtocolStageReport.self, #require(verified.cursor.priorStageReports.first), f.store)
        #expect(prior.minimization?.converged == true && prior.successful)
        let transition = try await read(VivoMDStageTransition.self, #require(verified.cursor.transition), f.store)
        #expect(transition.velocityInitialization == .maxwellBoltzmann && transition.thermalizationSeed == 73)
        #expect(verified.checkpoint.velocitiesNMPerPS.contains { $0 != .zero })
        let resumed = try await VivoMDProtocolRunner.resume(system: f.system, plan: plan, store: f.store, checkpoint: durable.fingerprint)
        let resumedReceipt = try await resumed.run()
        #expect(resumedReceipt.disposition == .completed && resumedReceipt.completedStages == 2)
        let resumedHash = try #require(resumedReceipt.latestDurableCheckpoint).fingerprint
        let resumedPrefix = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: plan,
                                                                         store: f.store, checkpoint: resumedHash)
        #expect(resumedPrefix.checkpoint == verified.checkpoint)
        #expect(resumedPrefix.cursor.transition == verified.cursor.transition)
        #expect(try await f.store.reference(runner.checkpointReferenceName).artifact.fingerprint == durable.fingerprint)

        // Exercise the unsuccessful native gate too: one tiny minimizer step
        // cannot remove this harmonic pair's initial finite bond error.
        plan.stages[0].minimization?.maximumIterations = 1
        plan.stages[0].minimization?.forceToleranceKJPerMolNM = 1e-12
        let blocked = try await VivoMDProtocolRunner.start(system: f.system, initialState: f.initial, plan: plan, store: f.store)
        let blockedReceipt = try await blocked.run()
        #expect(blockedReceipt.disposition == .rejected && blockedReceipt.completedStages == 0)
        let blockedHash = try #require(blockedReceipt.latestDurableCheckpoint).fingerprint
        let blockedResume = try await VivoMDProtocolRunner.resume(system: f.system, plan: plan, store: f.store, checkpoint: blockedHash)
        let stillBlocked = try await blockedResume.run()
        #expect(stillBlocked.disposition == .rejected && stillBlocked.completedStages == 0)
    }

    @Test func nativeResumeRepairsPartialPublicationsWithoutSkippingOrDuplicatingSteps() async throws {
        // A one-step stage also exercises recovery immediately before stage
        // finalization, when there is no next dynamics iteration to repair it.
        for requestedSteps in [UInt64(1), UInt64(4)] {
            let f = try await fixture(dynamicsSteps: requestedSteps); defer { try? FileManager.default.removeItem(at: f.root) }
            let runtime = try await VivoMDMetalRuntime.make(system: f.system, initialState: f.initial, configuration: configuration())
            #expect(try await runtime.step().committed)
            let accepted = try await runtime.checkpoint(), pair = try await runtime.sample(includeObservables: true)
            let acceptedHash = try await put(accepted, "md-checkpoint", f.store)
            let firstObservation = try #require(pair.observables)
            let observationHash = try await put(VivoMDObservationLink(schema: VivoMDObservationLink.schemaID,
                stageIdentifier: "prepare", ordinal: 0, previous: nil, observation: firstObservation), "md-observation-link", f.store)
            let writer = try VivoMDTrajectoryArchiveWriter(store: f.store, systemFingerprint: accepted.systemFingerprint,
                configurationFingerprint: accepted.configurationFingerprint, particleCount: 2, includeVelocities: true,
                targetChunkBytes: f.plan.trajectoryChunkBytes)
            try await writer.append(pair.state)
            let sampledManifest = try await writer.snapshot().fingerprint
            var expectedFrames = [VivoMDArchiveFrame(snapshot: pair.state, includeVelocities: true)]
            for _ in 1..<requestedSteps {
                #expect(try await runtime.step().committed)
                expectedFrames.append(VivoMDArchiveFrame(snapshot: try await runtime.snapshot(), includeVelocities: true))
            }
            let expected = try await runtime.checkpoint()

            // These are the actual durable states left by a failure before append,
            // between append and observation, or after both publications. Each uses
            // the same real Metal checkpoint and immutable prefixes as the runner.
            for (samplePresent, observationPresent) in [(false, false), (true, false), (false, true), (true, true)] {
                var cursor = f.cursor
                cursor.completedStepsInStage = 1; cursor.currentCheckpoint = acceptedHash
                if samplePresent { cursor.trajectoryManifest = sampledManifest }
                if observationPresent { cursor.observationCount = 1; cursor.observationTail = observationHash }
                let sourceHash = try await put(cursor, "md-protocol-checkpoint", f.store)
                let original = try await f.store.data(for: sourceHash)
                let runner = try await VivoMDProtocolRunner.resume(system: f.system, plan: f.plan, store: f.store, checkpoint: sourceHash)
                let receipt = try await runner.run()
                #expect(receipt.disposition == .completed)
                let durable = try #require(receipt.latestDurableCheckpoint)
                let prefix = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: f.plan,
                                                                           store: f.store, checkpoint: durable.fingerprint)
                #expect(prefix.checkpoint == expected)
                #expect(prefix.cursor.completedStepsInStage == requestedSteps && prefix.cursor.observationCount == requestedSteps)
                #expect(try await f.store.data(for: sourceHash) == original)
                let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: #require(prefix.cursor.trajectoryManifest))
                var actualFrames: [VivoMDArchiveFrame] = []
                for hash in try await archive.index() { actualFrames += try await archive.readChunk(hash) }
                #expect(actualFrames == expectedFrames)
                var tail = prefix.cursor.observationTail
                var steps: [UInt64] = []
                while let hash = tail {
                    let observation = try await read(VivoMDObservationLink.self, hash, f.store)
                    steps.append(observation.observation.stepIndex); tail = observation.previous
                }
                #expect(steps == Array((1...requestedSteps).reversed()))
            }
        }
    }

    @Test func actualObservationWriteFailurePreservesARecoverableAcceptedBoundary() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let reference = try await VivoMDMetalRuntime.make(system: f.system, initialState: f.initial, configuration: configuration())
        var expectedFrames: [VivoMDArchiveFrame] = []
        var expectedObservations: [VivoMDObservables] = []
        var firstCheckpoint: VivoMDCheckpoint?
        for step in 1...4 {
            #expect(try await reference.step().committed)
            let pair = try await reference.sample(includeObservables: true)
            expectedFrames.append(.init(snapshot: pair.state, includeVelocities: true))
            expectedObservations.append(try #require(pair.observables))
            if step == 1 { firstCheckpoint = try await reference.checkpoint() }
        }
        let finalCheckpoint = try await reference.checkpoint()
        let firstLink = VivoMDObservationLink(schema: VivoMDObservationLink.schemaID, stageIdentifier: "prepare",
            ordinal: 0, previous: nil, observation: try #require(expectedObservations.first))
        let hash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(firstLink)).hex
        let sentinel = f.root.appendingPathComponent("objects/sha256/\(hash.prefix(2))/\(hash.dropFirst(2).prefix(2))/\(hash)")
        // Rooted storage refuses a directory at an immutable object destination:
        // linkat sees EEXIST and the existing-object regular-file check fails.
        // This exercises the real production write/error/persist path without
        // adding an injectable failure seam to the runtime or artifact store.
        try #require(!FileManager.default.fileExists(atPath: sentinel.path))
        try FileManager.default.createDirectory(at: sentinel, withIntermediateDirectories: true)
        let runner = try await VivoMDProtocolRunner.start(system: f.system, initialState: f.initial, plan: f.plan, store: f.store)
        let failed = try await runner.run()
        #expect(failed.disposition == .failed && failed.diagnostic != nil && failed.persistenceDiagnostic == nil)
        let durable = try #require(failed.latestDurableCheckpoint)
        let savedBytes = try await f.store.data(for: durable.fingerprint)
        let prefix = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: f.plan,
                                                                  store: f.store, checkpoint: durable.fingerprint)
        #expect(prefix.checkpoint == firstCheckpoint)
        #expect(prefix.cursor.phase == .running && prefix.cursor.completedStepsInStage == 1)
        #expect(prefix.cursor.observationCount == 0 && prefix.lastSampledStep == 1 && prefix.lastObservedStep == nil)
        #expect(try await f.store.reference(runner.checkpointReferenceName).artifact.fingerprint == durable.fingerprint)
        try #require(FileManager.default.contentsOfDirectory(atPath: sentinel.path).isEmpty)
        try FileManager.default.removeItem(at: sentinel)

        // Both retrying this same actor and forking its failed durable prefix
        // must publish the missing observation once and then resume dynamics.
        let sameRunnerReceipt = try await runner.run()
        let restarted = try await VivoMDProtocolRunner.resume(system: f.system, plan: f.plan,
                                                               store: f.store, checkpoint: durable.fingerprint)
        let restartedReceipt = try await restarted.run()
        for receipt in [sameRunnerReceipt, restartedReceipt] {
            #expect(receipt.disposition == .completed)
            let finalHash = try #require(receipt.latestDurableCheckpoint).fingerprint
            let repaired = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: f.plan,
                                                                         store: f.store, checkpoint: finalHash)
            #expect(repaired.checkpoint == finalCheckpoint && repaired.cursor.observationCount == 4)
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: #require(repaired.cursor.trajectoryManifest))
            var frames: [VivoMDArchiveFrame] = []
            for link in try await archive.index() { frames += try await archive.readChunk(link) }
            #expect(frames == expectedFrames)
            var observations: [VivoMDObservables] = [], tail = repaired.cursor.observationTail
            while let hash = tail {
                let link = try await read(VivoMDObservationLink.self, hash, f.store)
                observations.append(link.observation); tail = link.previous
            }
            #expect(Array(observations.reversed()) == expectedObservations)
        }
        #expect(try await f.store.data(for: durable.fingerprint) == savedBytes)
    }

    @Test func actualStageReportWriteFailureCannotResetAnUnsuccessfulMinimizer() async throws {
        let f = try await fixture(minimization: true, twoStages: true)
        defer { try? FileManager.default.removeItem(at: f.root) }
        var plan = f.plan
        plan.stages[0].minimization?.maximumIterations = 1
        plan.stages[0].minimization?.forceToleranceKJPerMolNM = 1e-12
        let settings = try #require(plan.stages[0].minimization)
        let reference = try await VivoMDMetalRuntime.make(system: f.system, initialState: f.initial,
                                                          configuration: plan.stages[0].configuration)
        let entry = try await reference.checkpoint()
        let outcome = try await reference.minimize(settings)
        let exit = try await reference.checkpoint()
        #expect(!outcome.converged && outcome.attemptedIterations == 1 && outcome.acceptedIterations == 1)
        #expect(exit.positionsNM != entry.positionsNM && exit.acceptedStep == entry.acceptedStep)
        let expectedReport = VivoMDProtocolStageReport(schema: "numivivo.org/md-stage-report/v1",
            planFingerprint: try plan.fingerprint(), stageIdentifier: "prepare", stageIndex: 0, successful: false,
            entryCheckpoint: try entry.fingerprint(), exitCheckpoint: try exit.fingerprint(), transition: nil,
            committedSteps: 0, startTimePS: entry.timePS, endTimePS: exit.timePS, trajectoryManifest: nil,
            observationTail: nil, observationCount: 0, minimization: outcome, rejected: nil)
        let reportHash = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(expectedReport))
        let hex = reportHash.hex
        let sentinel = f.root.appendingPathComponent("objects/sha256/\(hex.prefix(2))/\(hex.dropFirst(2).prefix(2))/\(hex)")
        try #require(!FileManager.default.fileExists(atPath: sentinel.path))
        try FileManager.default.createDirectory(at: sentinel, withIntermediateDirectories: true)
        let runner = try await VivoMDProtocolRunner.start(system: f.system, initialState: f.initial, plan: plan, store: f.store)
        let initialReference = try await f.store.reference(runner.checkpointReferenceName)
        let initialBytes = try await f.store.data(for: initialReference.artifact.fingerprint)

        // The first attempt advances accepted geometry but fails its required
        // gate. Repeated storage failures must keep THAT terminal outcome; they
        // must not expose that geometry with a new running minimizer budget.
        for _ in 0..<2 {
            let failed = try await runner.run()
            #expect(failed.disposition == .failed && failed.persistenceDiagnostic != nil)
            #expect(failed.latestDurableCheckpoint == initialReference.artifact)
            #expect(try await f.store.reference(runner.checkpointReferenceName) == initialReference)
            let durable = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: plan, store: f.store,
                                                                       checkpoint: initialReference.artifact.fingerprint)
            #expect(durable.cursor.phase == .running && durable.cursor.activeStageReport == nil)
            #expect(durable.checkpoint == entry && durable.checkpoint != exit)
        }

        // A process-level resume can only restart the old durable entry, not
        // continue the unpublished failed geometry with another iteration budget.
        let resumed = try await VivoMDProtocolRunner.resume(system: f.system, plan: plan, store: f.store,
                                                            checkpoint: initialReference.artifact.fingerprint)
        let resumedInitial = try await f.store.reference(resumed.checkpointReferenceName)
        let resumedFailure = try await resumed.run()
        #expect(resumedFailure.disposition == .failed && resumedFailure.persistenceDiagnostic != nil)
        #expect(resumedFailure.latestDurableCheckpoint == resumedInitial.artifact)
        try #require(FileManager.default.contentsOfDirectory(atPath: sentinel.path).isEmpty)
        try FileManager.default.removeItem(at: sentinel)

        for pendingRunner in [runner, resumed] {
            let repaired = try await pendingRunner.run()
            #expect(repaired.disposition == .rejected && repaired.completedStages == 0)
            let hash = try #require(repaired.latestDurableCheckpoint).fingerprint
            let verified = try await VivoMDProtocolResumeValidation.load(system: f.system, plan: plan,
                                                                        store: f.store, checkpoint: hash)
            #expect(verified.cursor.phase == .blocked && verified.cursor.activeStageReport == reportHash)
            #expect(verified.checkpoint == exit)
            let storedReport = try await read(VivoMDProtocolStageReport.self, reportHash, f.store)
            #expect(storedReport == expectedReport && storedReport.minimization == outcome)
            let blockedResume = try await VivoMDProtocolRunner.resume(system: f.system, plan: plan, store: f.store, checkpoint: hash)
            let blockedReceipt = try await blockedResume.run()
            #expect(blockedReceipt.disposition == .rejected && blockedReceipt.completedStages == 0)
        }
        #expect(try await f.store.data(for: initialReference.artifact.fingerprint) == initialBytes)
    }
}
