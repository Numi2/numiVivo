import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

/// Synthetic integrity fixtures in real rooted storage. No test creates a native
/// runtime or presents these declared scalar measurements as sampled MD evidence.
@Suite(.serialized) struct MolecularSamplingArchiveReaderTests: Sendable {
    private struct Fixture: Sendable {
        let root: URL
        let store: VivoArtifactStore
        let request: VivoMolecularSamplingRunRequest
        let requestID: VivoFingerprint
        let cursor: VivoMolecularSamplingCursor
        let artifact: VivoStoredArtifact
        let checkpoints: [VivoMDCheckpoint]
        let diagnostic: VivoMolecularSamplingResult?
    }

    private func put<T: Encodable>(_ value: T, kind: String, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind, mediaType: "application/json")
    }

    private func directory() throws -> URL {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-sampling-reader-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        return root
    }

    private func randomValue(replica: Int, ordinal: Int) -> Double {
        var x = UInt64(ordinal + 1) &+ UInt64(replica + 1) &* 0x9e3779b97f4a7c15
        x = (x ^ (x >> 30)) &* 0xbf58476d1ce4e5b9
        x = (x ^ (x >> 27)) &* 0x94d049bb133111eb
        x ^= x >> 31
        return Double(Float(0.125 + Double(x >> 11) / 9_007_199_254_740_992 * 0.125))
    }

    private func fixture(blocks: Int = 2, framesPerBlock: Int = 2, passingLatest: Bool = false,
                         timeStepPS: Double = 1.0 / 256, acceptedCell: VivoPeriodicCell? = nil) async throws -> Fixture {
        let root = try directory()
        do {
            let store = try VivoArtifactStore(rootURL: root, createIfNeeded: false)
            let hydrogen = try #require(VivoElement.from(symbol: "H"))
            let initialCell: VivoPeriodicCell? = acceptedCell == nil ? nil : .init(a: .init(4, 0, 0), b: .init(0, 4, 0), c: .init(0, 0, 4))
            var structure = VivoMolecularStructure(identifier: "synthetic-sampling-reader", atoms: [
                .init(index: 0, name: "H1", element: hydrogen), .init(index: 1, name: "H2", element: hydrogen)
            ], bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: [.zero, .init(0.125, 0, 0)])])
            structure.periodicCell = initialCell
            let system = VivoClassicalSystem(identifier: "synthetic-sampling-reader",
                structureFingerprint: try VivoStructureCodec.fingerprint(structure), particles: [
                    .init(index: 0, atomIndex: 0, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
                    .init(index: 1, atomIndex: 1, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
                ])
            let systemID = try system.fingerprint(), seeds: [UInt64] = [17, 29]
            let states = seeds.indices.map { i in
                VivoClassicalInitialState(systemFingerprint: systemID, positionsNM: [.zero, .init(0.125, 0, 0)],
                    periodicCell: initialCell, sourceTimePS: i == 0 ? 0.25 : 0.75)
            }
            let md = VivoMDConfiguration(timeStepPS: timeStepPS, cutoffNM: 0.5, neighborSkinNM: 0,
                electrostatics: .cutoff, ensemble: acceptedCell == nil ? .nvt : .npt,
                barostat: acceptedCell == nil ? .none : .monteCarloIsotropic,
                targetPressureBar: acceptedCell == nil ? nil : 1, neighborListEnabled: false)
            let convergence = VivoMolecularSamplingConfiguration(minimumRetainedFramesPerReplica: passingLatest ? 1024 : 64,
                minimumReplicas: 2, maximumRHat: 1.1, minimumEffectiveSamplesPerReplica: 20,
                maximumAutocorrelationLag: 64)
            let request = VivoMolecularSamplingRunRequest(structure: structure, system: system,
                initialStates: states, replicaSeeds: seeds, md: md, contextIdentifier: "synthetic integrity fixture",
                observables: [.init(identifier: "distance", kind: .distance(atomA: 0, atomB: 1), maximumMeanStandardError: 100)],
                convergence: convergence, equilibrationSteps: 2, stepsPerBlock: UInt64(framesPerBlock * 2),
                sampleEvery: 2, maximumBlocks: max(2, blocks), requiredConsecutivePasses: 2, trajectoryChunkBytes: 1024 * 1024)
            try request.validate()
            let requestID = try await put(request, kind: "molecular-sampling-run", store: store).fingerprint
            var replicas: [VivoMolecularReplicaCursor] = [], writers: [VivoMDTrajectoryArchiveWriter] = []
            for i in seeds.indices {
                var cfg = md; cfg.randomSeed = seeds[i]
                let initialID = try await put(states[i], kind: "classical-initial-state", store: store).fingerprint
                replicas.append(.init(series: .init(identifier: "replica-\(i)", sourceFingerprint: initialID,
                    configuration: cfg, steps: [], timesPS: [], valuesByObservable: [[]]),
                    mdCheckpoint: nil, trajectory: nil, minimization: nil))
                writers.append(try .init(store: store, systemFingerprint: systemID,
                    configurationFingerprint: cfg.fingerprint(), particleCount: 2,
                    includeVelocities: false, targetChunkBytes: request.trajectoryChunkBytes))
            }
            var cursor = VivoMolecularSamplingCursor(schema: VivoMolecularSamplingCursor.schema,
                requestFingerprint: requestID, completedBlocks: 0, consecutivePasses: 0, replicas: replicas, diagnostics: [])
            var checkpoints: [VivoMDCheckpoint] = [], latest: VivoMolecularSamplingResult?
            var clocks = states.map { $0.sourceTimePS! }
            for i in clocks.indices {
                for _ in 0..<request.equilibrationSteps { clocks[i] += md.timeStepPS }
            }
            for block in 0..<blocks {
                checkpoints = []
                for i in seeds.indices {
                    let cfg = cursor.replicas[i].series.configuration
                    var final: VivoMDStateSnapshot?
                    for offset in 0..<framesPerBlock {
                        let ordinal = block * framesPerBlock + offset
                        let step = request.equilibrationSteps + UInt64(ordinal + 1) * request.sampleEvery
                        for _ in 0..<request.sampleEvery { clocks[i] += md.timeStepPS }
                        let distance = passingLatest ? randomValue(replica: i, ordinal: ordinal) : 0.125 + Double(ordinal % 3) / 128
                        let state = VivoMDStateSnapshot(systemFingerprint: systemID,
                            configurationFingerprint: try cfg.fingerprint(), stepIndex: step,
                            timePS: clocks[i],
                            positionsNM: [.zero, .init(distance, 0, 0)],
                            velocitiesNMPerPS: [.init(0.25, -0.125, 0.0625), .init(-0.25, 0.125, -0.0625)], periodicCell: acceptedCell)
                        try await writers[i].append(state); final = state
                        cursor.replicas[i].series.steps.append(step)
                        cursor.replicas[i].series.timesPS.append(state.timePS)
                        cursor.replicas[i].series.valuesByObservable[0].append(distance)
                    }
                    let state = try #require(final)
                    let checkpoint = VivoMDCheckpoint(systemFingerprint: systemID, configurationFingerprint: state.configurationFingerprint,
                        acceptedStep: state.stepIndex, timePS: state.timePS, positionsNM: state.positionsNM,
                        velocitiesNMPerPS: state.velocitiesNMPerPS, periodicCell: state.periodicCell)
                    try checkpoint.validate(particleCount: 2)
                    cursor.replicas[i].mdCheckpoint = try await put(checkpoint, kind: "md-checkpoint", store: store).fingerprint
                    let manifest = try await writers[i].snapshot()
                    cursor.replicas[i].trajectory = manifest.fingerprint
                    cursor.replicas[i].series.sourceFingerprint = manifest.fingerprint
                    checkpoints.append(checkpoint)
                }
                let analysis = VivoMolecularSamplingRequest(structureFingerprint: system.structureFingerprint,
                    systemFingerprint: systemID, contextIdentifier: request.contextIdentifier,
                    observables: request.observables, replicas: cursor.replicas.map(\.series), configuration: convergence)
                let diagnostic = try VivoMolecularSampling.analyze(analysis)
                cursor.completedBlocks += 1
                cursor.consecutivePasses = diagnostic.converged ? cursor.consecutivePasses + 1 : 0
                cursor.diagnostics.append(try await put(diagnostic, kind: "molecular-sampling-result", store: store).fingerprint)
                latest = diagnostic
            }
            let artifact = try await put(cursor, kind: "molecular-sampling-checkpoint", store: store)
            return .init(root: root, store: store, request: request, requestID: requestID,
                cursor: cursor, artifact: artifact, checkpoints: checkpoints, diagnostic: latest)
        } catch {
            try? FileManager.default.removeItem(at: root)
            throw error
        }
    }

    private func rejects(_ body: () async throws -> Void) async {
        do { try await body(); Issue.record("invalid sampling archive selection was accepted") } catch {}
    }

    private func selection(_ f: Fixture, cursor: VivoMolecularSamplingCursor? = nil,
                           limits: VivoMolecularSamplingReadLimits = .init()) async throws -> VivoValidatedMolecularSamplingSelection {
        let id: VivoFingerprint
        if let cursor { id = try await put(cursor, kind: "molecular-sampling-checkpoint", store: f.store).fingerprint }
        else { id = f.artifact.fingerprint }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: id, limits: limits)
        return try await reader.selectReplica(index: 0)
    }

    private func receipt(_ f: Fixture, status: VivoMolecularSamplingRunStatus,
                         checkpoint: VivoFingerprint? = nil) async throws -> VivoFingerprint {
        let value = VivoMolecularSamplingRunReceipt(schema: VivoMolecularSamplingRunReceipt.schema,
            requestFingerprint: f.requestID, status: status, completedBlocks: f.cursor.completedBlocks,
            consecutivePasses: f.cursor.consecutivePasses, checkpoint: checkpoint ?? f.artifact.fingerprint,
            trajectoryManifests: f.cursor.replicas.compactMap(\.trajectory), sampling: f.diagnostic,
            diagnostic: status == .cancelled || status == .rejected ? "synthetic termination record" : nil,
            interpretation: VivoMolecularSampling.interpretation)
        return try await put(value, kind: "molecular-sampling-receipt", store: f.store).fingerprint
    }

    @Test func selectedReplicaRetainsExactPayloadSeedAndSourceClock() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        let inspection = reader.inspect()
        #expect(inspection.completedBlocks == 2 && inspection.consecutivePasses == 0)
        #expect(!inspection.declaredObservableCriteriaSatisfied && inspection.replicas.count == 2)
        for index in 0..<2 {
            let selected = try await reader.selectReplica(index: index, validation: .allPayloads)
            let originalID = try #require(f.cursor.replicas[index].mdCheckpoint)
            #expect(selected.store === f.store)
            #expect(selected.requestFingerprint == f.requestID && selected.samplingCheckpointFingerprint == f.artifact.fingerprint)
            #expect(selected.replicaIndex == index && selected.replicaIdentifier == "replica-\(index)")
            #expect(selected.replicaSeed == f.request.replicaSeeds[index])
            #expect(selected.configuration.randomSeed == f.request.replicaSeeds[index])
            #expect(selected.mdCheckpoint == f.checkpoints[index] && selected.mdCheckpointFingerprint == originalID)
            #expect(selected.mdCheckpointData == (try await f.store.data(for: originalID)))
            #expect(selected.mdCheckpoint.acceptedStep == 10)
            #expect(selected.mdCheckpoint.timePS == f.request.initialStates[index].sourceTimePS! + 10.0 / 256)
            #expect(selected.trajectoryValidation.scope == .allPayloads && selected.trajectoryValidation.verifiedPayloads == 2)
            #expect(selected.trajectoryManifest.lastStep == selected.mdCheckpoint.acceptedStep)
            #expect(selected.sourceReceipt == nil && selected.sourceReceiptFingerprint == nil && selected.sourceTermination == .notRecorded)
        }
        for index in [-1, 2] { await rejects { _ = try await reader.selectReplica(index: index) } }
    }

    @Test func selectionRetainsOriginalNoncanonicalCheckpointBytes() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let checkpoint = f.checkpoints[0]
        var original = Data(" \n".utf8)
        original.append(try VivoCanonicalJSON.encode(checkpoint))
        original.append(Data("\n\t ".utf8))
        let stored = try await f.store.put(data: original, kind: "md-checkpoint", mediaType: "application/json")
        let canonicalID = try checkpoint.fingerprint()
        try #require(stored.fingerprint != canonicalID)
        var cursor = f.cursor
        cursor.replicas[0].mdCheckpoint = stored.fingerprint
        let selected = try await selection(f, cursor: cursor)
        #expect(selected.mdCheckpointData == original && selected.mdCheckpointFingerprint == stored.fingerprint)
        #expect(selected.mdCheckpoint == checkpoint)
        #expect(try selected.mdCheckpoint.fingerprint() == canonicalID)
        #expect(selected.trajectoryManifestFingerprint == f.cursor.replicas[0].trajectory)
    }

    @Test func nonbinaryTimestepPreservesRepeatedlyAccumulatedSourceClock() async throws {
        let f = try await fixture(framesPerBlock: 32, timeStepPS: 0.003)
        defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        var differsFromClosedForm = false
        for index in 0..<2 {
            let selected = try await reader.selectReplica(index: index)
            let sourceTime = try #require(f.request.initialStates[index].sourceTimePS)
            var accumulated = sourceTime
            for _ in 0..<selected.mdCheckpoint.acceptedStep { accumulated += 0.003 }
            #expect(selected.mdCheckpoint.acceptedStep == 130)
            #expect(selected.mdCheckpoint.timePS == accumulated)
            #expect(selected.mdCheckpointData == (try await f.store.data(for: selected.mdCheckpointFingerprint)))
            differsFromClosedForm = differsFromClosedForm || accumulated != sourceTime + Double(selected.mdCheckpoint.acceptedStep) * 0.003
        }
        #expect(differsFromClosedForm)
    }

    @Test func zeroBlocksNeverPromoteInitialStateOrOrphanArtifacts() async throws {
        let f = try await fixture(blocks: 0); defer { try? FileManager.default.removeItem(at: f.root) }
        var cfg = f.request.md; cfg.randomSeed = f.request.replicaSeeds[0]
        let cp = VivoMDCheckpoint(systemFingerprint: try f.request.system.fingerprint(), configurationFingerprint: try cfg.fingerprint(),
            acceptedStep: 4, timePS: 0.25 + 4.0 / 256, positionsNM: [.zero, .init(0.125, 0, 0)],
            velocitiesNMPerPS: [.init(0.25, 0, 0), .init(-0.25, 0, 0)], periodicCell: nil)
        let orphanCheckpoint = try await put(cp, kind: "md-checkpoint", store: f.store)
        let writer = try VivoMDTrajectoryArchiveWriter(store: f.store, systemFingerprint: cp.systemFingerprint,
            configurationFingerprint: cp.configurationFingerprint, particleCount: 2, includeVelocities: false)
        try await writer.append(.init(systemFingerprint: cp.systemFingerprint, configurationFingerprint: cp.configurationFingerprint,
            stepIndex: cp.acceptedStep, timePS: cp.timePS, positionsNM: cp.positionsNM, velocitiesNMPerPS: cp.velocitiesNMPerPS, periodicCell: nil))
        let orphanTrajectory = try await writer.snapshot()
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        #expect(reader.inspect().completedBlocks == 0)
        #expect(reader.inspect().replicas.allSatisfy { $0.mdCheckpointFingerprint == nil })
        await rejects { _ = try await reader.selectReplica(index: 0) }
        var publishedOrphan = f.cursor
        publishedOrphan.replicas[0].mdCheckpoint = orphanCheckpoint.fingerprint
        publishedOrphan.replicas[0].trajectory = orphanTrajectory.fingerprint
        publishedOrphan.replicas[0].series.sourceFingerprint = orphanTrajectory.fingerprint
        let invalidOrphan = try await put(publishedOrphan, kind: "molecular-sampling-checkpoint", store: f.store)
        await rejects {
            _ = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: invalidOrphan.fingerprint)
        }
        var wrongSource = f.cursor
        wrongSource.replicas[0].series.sourceFingerprint = f.requestID
        let invalidSource = try await put(wrongSource, kind: "molecular-sampling-checkpoint", store: f.store)
        await rejects {
            _ = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: invalidSource.fingerprint)
        }
    }

    @Test func acceptedNPTCellIsAuthoritativeAndReceivesCurrentCapabilityChecks() async throws {
        let accepted = VivoPeriodicCell(a: .init(4.25, 0, 0), b: .init(0, 4.25, 0), c: .init(0, 0, 4.25))
        let f = try await fixture(acceptedCell: accepted)
        defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await selection(f)
        #expect(selected.configuration.ensemble == .npt)
        #expect(selected.mdCheckpoint.periodicCell == accepted)
        #expect(selected.mdCheckpoint.periodicCell != f.request.initialStates[0].periodicCell)

        let skew = VivoPeriodicCell(a: .init(4.25, 0, 0), b: .init(0.5, 4.25, 0), c: .init(0, 0, 4.25))
        let unsupported = try await fixture(acceptedCell: skew)
        defer { try? FileManager.default.removeItem(at: unsupported.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: unsupported.store, checkpoint: unsupported.artifact.fingerprint)
        let manifest = try #require(unsupported.cursor.replicas[0].trajectory)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: unsupported.store, manifest: manifest)
        let tail = try #require(archive.manifest.tail)
        let frames = try await archive.readChunk(tail)
        #expect(frames.last?.periodicCell == skew && unsupported.checkpoints[0].periodicCell == skew)
        await rejects { _ = try await reader.selectReplica(index: 0) }
    }

    @Test(arguments: ["missing-contract", "old-contract", "configuration", "system", "step", "negative-time", "non-fp32", "tail-position", "tail-cell"])
    func changedCheckpointCannotAcquireSelectionAuthority(change: String) async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        var cp = f.checkpoints[0]
        switch change {
        case "missing-contract": cp.numericalContract = nil
        case "old-contract": cp.numericalContract = "numivivo.org/md-metal-numerics/v2"
        case "configuration": cp.configurationFingerprint = f.checkpoints[1].configurationFingerprint
        case "system": cp.systemFingerprint = f.requestID
        case "step": cp.acceptedStep -= 1
        case "negative-time": cp.timePS = -0.25
        case "non-fp32": cp.positionsNM[0].x = 0.1
        case "tail-position": cp.positionsNM[0].x = 0.125
        case "tail-cell": cp.periodicCell = .init(a: .init(4, 0, 0), b: .init(0, 4, 0), c: .init(0, 0, 4))
        default: Issue.record("unknown test mutation"); return
        }
        var cursor = f.cursor
        cursor.replicas[0].mdCheckpoint = try await put(cp, kind: "md-checkpoint", store: f.store).fingerprint
        let changed = cursor
        await rejects { _ = try await selection(f, cursor: changed) }
    }

    @Test(arguments: ["seed", "interior-step", "interior-time", "scalar-shape", "identifier"])
    func scalarPrefixCannotDetachFromReplicaSchedule(change: String) async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        var cursor = f.cursor
        switch change {
        case "seed": cursor.replicas[0].series.configuration.randomSeed = 999
        case "interior-step": cursor.replicas[0].series.steps[0] += 1
        case "interior-time": cursor.replicas[0].series.timesPS[0] += 1.0 / 512
        case "scalar-shape": cursor.replicas[0].series.valuesByObservable[0].removeLast()
        case "identifier": cursor.replicas[0].series.identifier = "replica-1"
        default: Issue.record("unknown test mutation"); return
        }
        let changed = cursor
        await rejects { _ = try await selection(f, cursor: changed) }
    }

    @Test func unrelatedValidPassingDiagnosticCannotForgeConsecutivePasses() async throws {
        let f = try await fixture(framesPerBlock: 512, passingLatest: true)
        defer { try? FileManager.default.removeItem(at: f.root) }
        try #require(f.diagnostic?.converged == true && f.cursor.consecutivePasses == 1)
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        #expect(!reader.inspect().declaredObservableCriteriaSatisfied && reader.inspect().diagnosticEvaluations == 2)
        await rejects { _ = try await selection(f, limits: .init(maximumDiagnosticEvaluations: 1)) }
        await rejects { _ = try await selection(f, limits: .init(maximumAutocorrelationWork: 1)) }
        var forged = f.cursor
        forged.diagnostics[0] = forged.diagnostics[1]
        forged.consecutivePasses = 2
        let changed = forged
        await rejects { _ = try await selection(f, cursor: changed) }
        let forgedID = try await put(changed, kind: "molecular-sampling-checkpoint", store: f.store).fingerprint
        // At maximumBlocks even the previous buggy fast path cannot allocate a
        // runtime: it would simply return a falsely converged receipt. Exercise
        // the public runner's shared admission before that terminal fast path.
        await rejects {
            _ = try await VivoMolecularSamplingRunner.run(f.request, store: f.store, resumeFrom: forgedID)
        }
    }

    @Test func selectedPayloadScopeNeverOverstatesOlderCoordinateVerification() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        let manifest = try #require(f.cursor.replicas[0].trajectory)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: manifest)
        let links = try await archive.index()
        try #require(links.count == 2)
        let first = try await archive.readLink(links[0])
        let firstDescriptor = try await f.store.descriptor(for: first.payload)
        let firstPath = f.root.appendingPathComponent(firstDescriptor.objectPath)
        let firstOriginal = try Data(contentsOf: firstPath)
        var corrupt = firstOriginal; corrupt[0] ^= 1
        try corrupt.write(to: firstPath)
        let selected = try await reader.selectReplica(index: 0)
        #expect(selected.trajectoryValidation.scope == .restart && selected.trajectoryValidation.verifiedPayloads == 1)
        await rejects { _ = try await reader.selectReplica(index: 0, validation: .allPayloads) }
        try firstOriginal.write(to: firstPath)
        let last = try await archive.readLink(links[1])
        let lastDescriptor = try await f.store.descriptor(for: last.payload)
        let lastPath = f.root.appendingPathComponent(lastDescriptor.objectPath)
        var finalCorruption = try Data(contentsOf: lastPath); finalCorruption[0] ^= 1
        try finalCorruption.write(to: lastPath)
        await rejects { _ = try await reader.selectReplica(index: 0) }
    }

    @Test func terminalReceiptsBindTheExactAcceptedPrefixWithoutPromotingCriteria() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        for status in [VivoMolecularSamplingRunStatus.cancelled, .rejected, .budgetExhausted] {
            let hash = try await receipt(f, status: status)
            let selected = try await reader.selectReplica(index: 0, receipt: hash)
            #expect(selected.sourceReceiptFingerprint == hash && selected.sourceReceipt?.status == status)
            #expect(!selected.inspection.declaredObservableCriteriaSatisfied)
            switch status {
            case .cancelled: #expect(selected.sourceTermination == .cancelled)
            case .rejected: #expect(selected.sourceTermination == .rejected)
            case .budgetExhausted: #expect(selected.sourceTermination == .budgetExhausted)
            case .converged: Issue.record("unexpected test status")
            }
        }
        let stale = try await receipt(f, status: .cancelled, checkpoint: f.requestID)
        await rejects { _ = try await reader.selectReplica(index: 0, receipt: stale) }
        let falseSuccess = try await receipt(f, status: .converged)
        await rejects { _ = try await reader.selectReplica(index: 0, receipt: falseSuccess) }
    }

    @Test func resumeValidatesEveryReplicaWhileSelectionVerifiesTheChosenPayload() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        let secondManifest = try #require(f.cursor.replicas[1].trajectory)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: secondManifest)
        let tail = try #require(archive.manifest.tail)
        let link = try await archive.readLink(tail)
        let descriptor = try await f.store.descriptor(for: link.payload)
        let path = f.root.appendingPathComponent(descriptor.objectPath)
        let original = try Data(contentsOf: path)
        var corrupt = original; corrupt[0] ^= 1
        try corrupt.write(to: path)
        _ = try await reader.selectReplica(index: 0)
        await rejects { _ = try await reader.validateForResume(request: f.request) }
        try original.write(to: path)
        #expect(try await reader.validateForResume(request: f.request).count == 2)
    }

    @Test func allReplicaResumeSharesOneCumulativeReadBudget() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        var metadata = UInt64(f.cursor.replicas.count) * 64 * 1024
        let ids = [f.artifact.fingerprint, f.requestID]
            + f.cursor.replicas.compactMap(\.mdCheckpoint)
            + f.cursor.replicas.compactMap(\.trajectory)
            + [try #require(f.cursor.diagnostics.last)]
        for id in ids { metadata += try await f.store.descriptor(for: id).byteCount }
        let manifestID = try #require(f.cursor.replicas[0].trajectory)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: manifestID)
        let manifest = archive.manifest
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: manifest.particleCount, includeVelocities: false)
        let wire = manifest.frameCount * UInt64(stride) + manifest.chunkCount * UInt64(VivoMDTrajectoryChunkCodec.headerBytes)
        let tailCeiling = min(wire, UInt64(VivoMDTrajectoryChunkCodec.maximumChunkBytes))
        // The reader documents these conservative admission charges. One
        // selection fits exactly; independently resetting that budget for each
        // resumed replica would incorrectly admit both trajectories.
        let oneTraversal = manifest.chunkCount * 64 * 1024 + 2 * tailCeiling + 64 * 1024
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint,
            limits: .init(maximumReadBytes: metadata + oneTraversal))
        _ = try await reader.selectReplica(index: 0)
        await rejects { _ = try await reader.validateForResume(request: f.request) }
        let sufficient = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint,
            limits: .init(maximumReadBytes: metadata + 2 * oneTraversal))
        #expect(try await sufficient.validateForResume(request: f.request).count == 2)
    }

    @Test func readingDoesNotPublishOrFollowReplacementRoot() async throws {
        let f = try await fixture()
        let moved = f.root.deletingLastPathComponent().appendingPathComponent(f.root.lastPathComponent + "-moved")
        defer { try? FileManager.default.removeItem(at: f.root); try? FileManager.default.removeItem(at: moved) }
        let name = VivoMolecularSamplingRunner.checkpointReferenceName(requestFingerprint: f.requestID)
        #expect(name == "molecular-sampling-" + f.requestID.hex + "-checkpoint")
        #expect(name != VivoMolecularSamplingRunner.checkpointReferenceName(requestFingerprint: f.artifact.fingerprint))
        _ = try await f.store.setReference(name, to: f.artifact)
        await rejects { _ = try await f.store.reference(name, maximumObjectBytes: 1) }
        await rejects { _ = try await f.store.reference(name, maximumObjectBytes: -1) }
        #expect(try await f.store.reference(name, maximumObjectBytes: Int(f.artifact.byteCount)).artifact == f.artifact)
        let before = try await f.store.list()
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        try FileManager.default.moveItem(at: f.root, to: moved)
        guard Darwin.mkdir(f.root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        let replacement = try VivoArtifactStore(rootURL: f.root, createIfNeeded: false)
        let selected = try await reader.selectReplica(index: 0)
        let selectedRootURL = await selected.store.rootURL
        let replacementRootURL = await replacement.rootURL
        #expect(selected.store === f.store && selectedRootURL == replacementRootURL)
        #expect(try await f.store.data(for: selected.mdCheckpointFingerprint) == selected.mdCheckpointData)
        #expect(try await f.store.list() == before)
        #expect(try await replacement.list().isEmpty)
        #expect(try await f.store.reference(name).artifact.fingerprint == f.artifact.fingerprint)
    }

    @Test func cancelledReaderAndSelectionDoNotPublish() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let before = try await f.store.list()
        let openTask = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        }
        do { _ = try await openTask.value; Issue.record("pre-cancelled reader opened") } catch is CancellationError {}
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.artifact.fingerprint)
        let selectTask = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await reader.selectReplica(index: 0)
        }
        do { _ = try await selectTask.value; Issue.record("pre-cancelled selection succeeded") } catch is CancellationError {}
        #expect(try await f.store.list() == before)
    }

    @Test func explicitAllocationScheduleAndTraversalLimitsAreEnforced() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let limits: [VivoMolecularSamplingReadLimits] = [
            .init(maximumRequestBytes: 1), .init(maximumCursorBytes: 1), .init(maximumCheckpointBytes: 1),
            .init(maximumDiagnosticBytes: 1), .init(maximumReadBytes: 1), .init(maximumReplicas: 1),
            .init(maximumParticles: 1), .init(maximumScalarElements: 7), .init(maximumAcceptedSteps: 9),
            .init(maximumTrajectoryChunks: 1)
        ]
        for limit in limits { await rejects { _ = try await selection(f, limits: limit) } }
        _ = try await selection(f, limits: .init(maximumReplicas: 2, maximumParticles: 2,
            maximumScalarElements: 8, maximumAcceptedSteps: 10, maximumTrajectoryChunks: 2))
    }
}
