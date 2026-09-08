import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

/// Small synthetic accepted prefixes in real rooted storage. These host tests
/// never create a native runtime or claim sampled/converged molecular evidence.
@Suite(.serialized) struct SamplingReactionSeedTests: Sendable {
    private struct Fixture: Sendable {
        let root: URL
        let store: VivoArtifactStore
        let request: VivoMolecularSamplingRunRequest
        let requestID: VivoFingerprint
        let cursor: VivoMolecularSamplingCursor
        let checkpointID: VivoFingerprint
        let diagnostic: VivoMolecularSamplingResult
        let implementation: VivoFingerprint
    }

    private func put<T: Encodable>(_ value: T, kind: String, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind, mediaType: "application/json")
    }

    private func fixture(noncanonicalCheckpoint: Bool = false, periodic: Bool = false) async throws -> Fixture {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-sampling-seed-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        do {
            let store = try VivoArtifactStore(rootURL: root, createIfNeeded: false)
            let cell: VivoPeriodicCell? = periodic ? .init(a: .init(2, 0, 0), b: .init(0, 2, 0), c: .init(0, 0, 2)) : nil
            let hydrogen = try #require(VivoElement.from(symbol: "H"))
            let structure = VivoMolecularStructure(identifier: "synthetic-sampling-export", atoms: [
                .init(index: 0, name: "H1", element: hydrogen), .init(index: 1, name: "H2", element: hydrogen)
            ], bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: [.init(0.125, 0, 0), .zero])])
            let system = VivoClassicalSystem(identifier: "synthetic-sampling-export",
                structureFingerprint: try VivoStructureCodec.fingerprint(structure), particles: [
                    .init(index: 0, atomIndex: 1, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
                    .init(index: 1, atomIndex: 0, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
                ])
            let systemID = try system.fingerprint(), seeds: [UInt64] = [17, 29]
            let states = seeds.indices.map { index in
                VivoClassicalInitialState(systemFingerprint: systemID, positionsNM: [.zero, .init(0.125, 0, 0)],
                    periodicCell: cell, sourceTimePS: index == 0 ? 0.25 : 0.75)
            }
            let md = VivoMDConfiguration(timeStepPS: 1.0 / 256, cutoffNM: 0.5, neighborSkinNM: 0,
                electrostatics: .cutoff, ensemble: .nvt, thermostat: .langevinMiddle,
                targetTemperatureK: 300, frictionPerPS: 1, neighborListEnabled: false)
            let request = VivoMolecularSamplingRunRequest(structure: structure, system: system,
                initialStates: states, replicaSeeds: seeds, md: md, contextIdentifier: "synthetic export integrity fixture",
                observables: [.init(identifier: "distance", kind: .distance(atomA: 0, atomB: 1), maximumMeanStandardError: 100)],
                convergence: .init(minimumRetainedFramesPerReplica: 64, minimumReplicas: 2,
                    maximumRHat: 1.1, minimumEffectiveSamplesPerReplica: 20, maximumAutocorrelationLag: 3),
                equilibrationSteps: 2, stepsPerBlock: 4, sampleEvery: 2, maximumBlocks: 2,
                requiredConsecutivePasses: 2, trajectoryChunkBytes: 1024 * 1024)
            try request.validate()
            let requestID = try await put(request, kind: "molecular-sampling-run", store: store).fingerprint
            var replicas: [VivoMolecularReplicaCursor] = []
            for index in seeds.indices {
                var configuration = md; configuration.randomSeed = seeds[index]
                let configurationID = try configuration.fingerprint()
                let writer = try VivoMDTrajectoryArchiveWriter(store: store, systemFingerprint: systemID,
                    configurationFingerprint: configurationID, particleCount: 2, includeVelocities: false,
                    targetChunkBytes: request.trajectoryChunkBytes)
                var steps: [UInt64] = [], times: [Double] = [], distances: [Double] = []
                var final: VivoMDStateSnapshot?, manifest: VivoStoredArtifact?
                for ordinal in 0..<2 {
                    let step = UInt64(4 + ordinal * 2)
                    let distance = 0.125 + Double(ordinal + index) / 128
                    let state = VivoMDStateSnapshot(systemFingerprint: systemID, configurationFingerprint: configurationID,
                        stepIndex: step, timePS: states[index].sourceTimePS! + Double(step) / 256,
                        positionsNM: [.zero, .init(distance, 0, 0)],
                        velocitiesNMPerPS: [.init(0.25, -0.125, 0.0625), .init(-0.25, 0.125, -0.0625)], periodicCell: cell)
                    try await writer.append(state)
                    // Force two valid chunks so allPayloads is stronger than restart.
                    manifest = try await writer.snapshot(); final = state
                    steps.append(step); times.append(state.timePS); distances.append(distance)
                }
                let state = try #require(final), trajectory = try #require(manifest)
                let checkpoint = VivoMDCheckpoint(systemFingerprint: systemID, configurationFingerprint: configurationID,
                    acceptedStep: state.stepIndex, timePS: state.timePS, positionsNM: state.positionsNM,
                    velocitiesNMPerPS: state.velocitiesNMPerPS, periodicCell: state.periodicCell)
                try checkpoint.validate(particleCount: 2)
                var bytes = try VivoCanonicalJSON.encode(checkpoint)
                if noncanonicalCheckpoint && index == 0 { bytes = Data(" \n".utf8) + bytes + Data("\n\t ".utf8) }
                let checkpointID = try await store.put(data: bytes, kind: "md-checkpoint", mediaType: "application/json").fingerprint
                replicas.append(.init(series: .init(identifier: "replica-\(index)", sourceFingerprint: trajectory.fingerprint,
                    configuration: configuration, steps: steps, timesPS: times, valuesByObservable: [distances]),
                    mdCheckpoint: checkpointID, trajectory: trajectory.fingerprint, minimization: nil))
            }
            let diagnostic = try VivoMolecularSampling.analyze(.init(structureFingerprint: system.structureFingerprint,
                systemFingerprint: systemID, contextIdentifier: request.contextIdentifier,
                observables: request.observables, replicas: replicas.map(\.series), configuration: request.convergence))
            let diagnosticID = try await put(diagnostic, kind: "molecular-sampling-result", store: store).fingerprint
            let cursor = VivoMolecularSamplingCursor(schema: VivoMolecularSamplingCursor.schema, requestFingerprint: requestID,
                completedBlocks: 1, consecutivePasses: 0, replicas: replicas, diagnostics: [diagnosticID])
            let checkpointID = try await put(cursor, kind: "molecular-sampling-checkpoint", store: store).fingerprint
            return .init(root: root, store: store, request: request, requestID: requestID, cursor: cursor,
                checkpointID: checkpointID, diagnostic: diagnostic,
                implementation: try VivoCanonicalJSON.fingerprint(Data("synthetic-export-implementation".utf8)))
        } catch {
            try? FileManager.default.removeItem(at: root)
            throw error
        }
    }

    private func select(_ f: Fixture, receipt: VivoFingerprint? = nil,
                        validation: VivoMolecularSamplingSelectionValidation = .restart) async throws -> VivoValidatedMolecularSamplingSelection {
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.checkpointID)
        return try await reader.selectReplica(index: 0, validation: validation, receipt: receipt)
    }

    private func rejects(_ body: () async throws -> Void) async {
        do { try await body(); Issue.record("invalid sampling reaction seed was accepted") } catch {}
    }

    private func inventory(_ root: URL) throws -> [String: VivoFingerprint] {
        let enumerator = try #require(FileManager.default.enumerator(at: root,
            includingPropertiesForKeys: [.isRegularFileKey]))
        var files: [String: VivoFingerprint] = [:]
        for case let file as URL in enumerator {
            if try file.resourceValues(forKeys: [.isRegularFileKey]).isRegularFile == true {
                let relative = String(file.path.dropFirst(root.path.count + 1))
                files[relative] = try VivoCanonicalJSON.fingerprint(Data(contentsOf: file))
            }
        }
        return files
    }

    private func corrupt(_ fingerprint: VivoFingerprint, fixture f: Fixture) async throws {
        let descriptor = try await f.store.descriptor(for: fingerprint)
        let path = f.root.appendingPathComponent(descriptor.objectPath)
        var data = try Data(contentsOf: path)
        try #require(!data.isEmpty)
        data[0] ^= 1
        try data.write(to: path)
    }

    private func source(_ f: Fixture) async throws -> VivoVerifiedMolecularSamplingExport {
        let publication = try await VivoMolecularSamplingExporter.publish(select(f), implementationFingerprint: f.implementation)
        return try await VivoMolecularSamplingExporter.verifiedExport(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
    }

    private func request(_ source: VivoVerifiedMolecularSamplingExport, connected: Bool = false) throws -> VivoSamplingReactionSeedRequest {
        .init(sourceExport: source.receipt,
            destination: try VivoReactionQualificationWorkflow.template(connected ? "h3-connected-rate" : "h2-minimum"),
            assignments: connected ? [
                .init(target: .connectedEndpoint(identifier: "H0-H1_plus_H2", componentIndex: 0), sourceAtomByNucleus: [0, 1]),
                .init(target: .connectedEndpoint(identifier: "H0_plus_H1-H2", componentIndex: 1), sourceAtomByNucleus: [0, 1])
            ] : [.init(target: .qualification, sourceAtomByNucleus: [0, 1])],
            modelTransferStatement: "Synthetic accepted classical prefix (mass 1 Da, target 300 K) supplies coordinates only. Fresh full-CI STO-3G qualification retains 1.008 Da and 298.15 K. These are different Hamiltonians and thermodynamic contexts; the prefix is not a QM ensemble.")
    }

    @Test func remapsChemicalAtomOrderAndPreservesEveryDestinationSetting() async throws {
        let f = try await fixture(noncanonicalCheckpoint: true); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f)
        var request = try request(source)
        request.assignments[0].sourceAtomByNucleus = [1, 0]
        guard case .qualify(var expected) = request.destination.calculation else { throw VivoChemistryError.invalid("test seed") }
        // Particle order is deliberately the reverse of chemical atom order.
        let nm = source.selection.mdCheckpoint.positionsNM
        expected.model.system.nuclei[0].positionBohr = SIMD3(nm[0].x, nm[0].y, nm[0].z) / VivoAtomicUnits.bohrInNM
        expected.model.system.nuclei[1].positionBohr = SIMD3(nm[1].x, nm[1].y, nm[1].z) / VivoAtomicUnits.bohrInNM
        let assembled = try VivoSamplingReactionSeed.assemble(request, source: source)
        #expect(assembled == .init(.qualify(request: expected)))
        #expect(!source.receipt.provenance.declaredObservableCriteriaSatisfied)
        let publication = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        let lineage = publication.receipt.provenance
        #expect(lineage.source == source.receipt.provenance)
        #expect(lineage.source.mdCheckpointFingerprint != lineage.source.canonicalCheckpointFingerprint)
        #expect(lineage.sourceConfiguration.targetTemperatureK == 300)
        #expect(lineage.transfers[0].sourceParticleByNucleus == [0, 1])
        #expect(lineage.transfers[0].sourceMassesDa == [1, 1])
        #expect(lineage.transfers[0].destinationMassesDa == [1.008, 1.008])
        #expect(lineage.transfers[0].destinationThermochemistry.temperatureK == 298.15)
        #expect(lineage.transfers[0].destinationSolver == .fullCI)
        #expect(lineage.modelTransferStatement == request.modelTransferStatement)
        #expect(lineage.interpretation == VivoSamplingReactionSeed.interpretation)
        let output = try #require(publication.receipt.outputs.first { $0.name == "request" })
        #expect(output.kind == "vivo.reaction-calculation-request")
        let data = try await VivoChemistryWorkflow(store: f.store).payload(artifact: output.artifact, expectedKind: output.kind)
        #expect(try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self, from: data) == assembled)
        #expect(try VivoCanonicalJSON.fingerprint(data) == lineage.assembledReactionRequestFingerprint)
        #expect(try await f.store.descriptor(for: output.artifact).kind == "chemistry-output")
    }

    @Test func bothConnectedEndpointSeedsChangeWhileSaddleAndIsolatedAtomsRemainExact() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), request = try request(source, connected: true)
        guard case .connectedReaction(var expected) = request.destination.calculation else { throw VivoChemistryError.invalid("test connected seed") }
        let nm = source.selection.mdCheckpoint.positionsNM
        let positions = [nm[1], nm[0]].map { SIMD3($0.x, $0.y, $0.z) / VivoAtomicUnits.bohrInNM }
        for (endpoint, component) in [(0, 0), (1, 1)] {
            for atom in 0..<2 { expected.endpoints[endpoint].components[component].qualification.model.system.nuclei[atom].positionBohr = positions[atom] }
        }
        #expect(try VivoSamplingReactionSeed.assemble(request, source: source) == .init(.connectedReaction(request: expected)))
        let publication = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        #expect(publication.receipt.provenance.transfers.count == 2)
        #expect(publication.receipt.provenance.transfers[0].assembledQualificationFingerprint == publication.receipt.provenance.transfers[1].assembledQualificationFingerprint)
        _ = try await VivoSamplingReactionSeed.verify(publication.receipt, store: f.store, implementationFingerprint: f.implementation)
    }

    @Test func rejectsIncompleteDuplicateOutOfRangeAndMismatchedElementMappingsBeforePublication() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), original = try request(source), before = try inventory(f.root)
        for indices in [[0], [0, 0], [-1, 1], [0, 2], [0, 1, 2]] {
            var altered = original; altered.assignments[0].sourceAtomByNucleus = indices
            await rejects { _ = try await VivoSamplingReactionSeed.publish(altered, source: source, implementationFingerprint: f.implementation) }
        }
        var element = original
        guard case .qualify(var q) = element.destination.calculation else { throw VivoChemistryError.invalid("test seed") }
        q.model.system.nuclei[0].atomicNumber = 2
        element.destination = .init(.qualify(request: q))
        await rejects { _ = try await VivoSamplingReactionSeed.publish(element, source: source, implementationFingerprint: f.implementation) }
        q.model.system.nuclei[0].atomicNumber = 1
        q.model.system.pointCharges = [.init(chargeE: 0.5, positionBohr: .init(5, 0, 0))]
        element.destination = .init(.qualify(request: q))
        await rejects { _ = try await VivoSamplingReactionSeed.publish(element, source: source, implementationFingerprint: f.implementation) }
        #expect(try inventory(f.root) == before)
    }

    @Test func exactVerifiedExportCannotBeReplacedByAnotherValidReplica() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f)
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: f.store, checkpoint: f.checkpointID)
        let alternate = try await VivoMolecularSamplingExporter.publish(reader.selectReplica(index: 1), implementationFingerprint: f.implementation)
        var mismatched = try request(source); mismatched.sourceExport = alternate.receipt
        let before = try inventory(f.root)
        await rejects { _ = try await VivoSamplingReactionSeed.publish(mismatched, source: source, implementationFingerprint: f.implementation) }
        #expect(try inventory(f.root) == before)
    }

    @Test func rejectsAbsentDuplicateWrongAndNonSeedTargets() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), original = try request(source, connected: true), before = try inventory(f.root)
        var bad: [VivoSamplingReactionSeedRequest] = []
        for target: VivoSamplingReactionSeedTarget in [.qualification, .connectedEndpoint(identifier: "missing", componentIndex: 0),
            .connectedEndpoint(identifier: "H0-H1_plus_H2", componentIndex: -1),
            .connectedEndpoint(identifier: "H0-H1_plus_H2", componentIndex: 2), .connectedSaddle] {
            var altered = original; altered.assignments = [.init(target: target, sourceAtomByNucleus: [0, 1])]; bad.append(altered)
        }
        var duplicate = original; duplicate.assignments.append(duplicate.assignments[0]); bad.append(duplicate)
        var empty = original; empty.assignments = []; bad.append(empty)
        var nonSeed = original; nonSeed.destination = try VivoReactionQualificationWorkflow.template("h2-global-embedding"); bad.append(nonSeed)
        for altered in bad { await rejects { _ = try await VivoSamplingReactionSeed.publish(altered, source: source, implementationFingerprint: f.implementation) } }
        #expect(try inventory(f.root) == before)
    }

    @Test func rejectsImplicitPeriodicEnvironmentRemoval() async throws {
        let f = try await fixture(periodic: true); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), request = try request(source), before = try inventory(f.root)
        try #require(source.selection.mdCheckpoint.periodicCell != nil)
        await rejects { _ = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation) }
        #expect(try inventory(f.root) == before)
    }

    @Test func immutableVerificationIgnoresMissingAndRepointedCachesWithoutWriting() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), request = try request(source)
        let first = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        let second = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        #expect(!first.reused && second.reused && first.receipt == second.receipt && first.artifact == second.artifact)
        let reference = "chemistry-task-" + first.receipt.taskFingerprint.hex
        let exportReference = "chemistry-task-" + source.receipt.taskFingerprint.hex
        try await f.store.removeReference(reference); try await f.store.removeReference(exportReference)
        let before = try inventory(f.root)
        let report = try await VivoSamplingReactionSeed.verify(first.receipt, store: f.store,
            implementationFingerprint: f.implementation, minimumValidation: .allPayloads)
        #expect(report.receiptFingerprint == first.artifact.fingerprint && report.verifiedOutputCount == 2)
        #expect(report.sourceVerification.validationScope == "allPayloads")
        #expect(report.sourceVerification.recordedValidationScope == "restart")
        #expect(try inventory(f.root) == before)
        let wrong = try await f.store.descriptor(for: f.requestID)
        _ = try await f.store.setReference(reference, to: wrong)
        let repointed = try inventory(f.root)
        _ = try await VivoSamplingReactionSeed.verify(first.receipt, store: f.store, implementationFingerprint: f.implementation)
        #expect(try inventory(f.root) == repointed)
        #expect(try await f.store.reference(reference).artifact == wrong)
    }

    @Test(arguments: ["source-tail", "request-output", "workflow-receipt", "assembly-input"])
    func corruptedSourceOrImmutableHandoffCannotVerifyAndIsNeverRegenerated(_ target: String) async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), request = try request(source)
        let publication = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        let victim: VivoFingerprint
        switch target {
        case "source-tail":
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: source.selection.trajectoryManifestFingerprint)
            victim = try await archive.readLink(try #require(archive.manifest.tail)).payload
        case "request-output": victim = try #require(publication.receipt.outputs.first { $0.name == "request" }).artifact
        case "workflow-receipt": victim = publication.receipt.workflowReceiptFingerprint
        default: victim = try #require(publication.receipt.task.inputs.first { $0.name == "assembly" }).artifact
        }
        try await corrupt(victim, fixture: f)
        let before = try inventory(f.root)
        await rejects { _ = try await VivoSamplingReactionSeed.verify(publication.receipt, store: f.store, implementationFingerprint: f.implementation) }
        #expect(try inventory(f.root) == before)
    }

    @Test func forgedReceiptAndWrongImplementationCannotRebindAValidWorkflow() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), request = try request(source)
        let publication = try await VivoSamplingReactionSeed.publish(request, source: source, implementationFingerprint: f.implementation)
        var json = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(publication.receipt)) as? [String: Any])
        var body = try #require(json["request"] as? [String: Any]); body["modelTransferStatement"] = "Altered model-transfer assertion"; json["request"] = body
        let altered = try VivoCanonicalJSON.decode(VivoSamplingReactionSeedReceipt.self, from: JSONSerialization.data(withJSONObject: json, options: [.sortedKeys]))
        let wrong = try VivoCanonicalJSON.fingerprint(Data("different-implementation".utf8)), before = try inventory(f.root)
        await rejects { _ = try await VivoSamplingReactionSeed.verify(altered, store: f.store, implementationFingerprint: f.implementation) }
        await rejects { _ = try await VivoSamplingReactionSeed.verify(publication.receipt, store: f.store, implementationFingerprint: wrong) }
        await rejects { _ = try await VivoSamplingReactionSeed.verify(publication.receipt, store: f.store, implementationFingerprint: f.implementation, samplingImplementationFingerprint: wrong) }
        #expect(try inventory(f.root) == before)
    }

    @Test func boundedStatementBudgetAndPreCancellationFailBeforeAnyWrite() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let source = try await source(f), original = try request(source), before = try inventory(f.root)
        for statement in [" \n\t ", String(repeating: "x", count: 8193), "hidden\u{0}context"] {
            var altered = original; altered.modelTransferStatement = statement
            await rejects { _ = try await VivoSamplingReactionSeed.publish(altered, source: source, implementationFingerprint: f.implementation) }
        }
        var tooSmall = original
        guard case .qualify(var q) = original.destination.calculation else { throw VivoChemistryError.invalid("test seed") }
        q.model.budget.maximumBytes = 1; tooSmall.destination = .init(.qualify(request: q))
        await rejects { _ = try await VivoSamplingReactionSeed.publish(tooSmall, source: source, implementationFingerprint: f.implementation) }
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await VivoSamplingReactionSeed.publish(original, source: source, implementationFingerprint: f.implementation)
        }
        do { _ = try await task.value; Issue.record("cancelled assembly published") } catch is CancellationError {}
        #expect(try inventory(f.root) == before)
    }
}
