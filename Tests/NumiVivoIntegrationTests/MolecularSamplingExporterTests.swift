import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

/// Small synthetic accepted prefixes in real rooted storage. These host tests
/// never create a native runtime or claim sampled/converged molecular evidence.
@Suite(.serialized) struct MolecularSamplingExporterTests: Sendable {
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

    private struct WorkflowReceipt: Codable {
        let schema: String
        let taskFingerprint: VivoFingerprint
        let outputs: [VivoChemistryOutputReceipt]
    }

    private func put<T: Encodable>(_ value: T, kind: String, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind, mediaType: "application/json")
    }

    private func fixture(noncanonicalCheckpoint: Bool = false) async throws -> Fixture {
        let root = FileManager.default.temporaryDirectory.appendingPathComponent("numivivo-sampling-export-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        do {
            let store = try VivoArtifactStore(rootURL: root, createIfNeeded: false)
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
                    sourceTimePS: index == 0 ? 0.25 : 0.75)
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
                        velocitiesNMPerPS: [.init(0.25, -0.125, 0.0625), .init(-0.25, 0.125, -0.0625)], periodicCell: nil)
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

    private func payload(_ name: String, publication: VivoMolecularSamplingExportPublication,
                         store: VivoArtifactStore) async throws -> Data {
        let output = try #require(publication.receipt.outputs.first { $0.name == name })
        return try await VivoChemistryWorkflow(store: store).payload(artifact: output.artifact, expectedKind: output.kind)
    }

    private func rejects(_ body: () async throws -> Void) async {
        do { try await body(); Issue.record("invalid sampling export was accepted") } catch {}
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

    @Test func publishesOriginalCheckpointAndMappedOutputsWithoutRestampingState() async throws {
        let f = try await fixture(noncanonicalCheckpoint: true); defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await select(f)
        let original = try await f.store.data(for: selected.mdCheckpointFingerprint)
        let originalDescriptor = try await f.store.descriptor(for: selected.mdCheckpointFingerprint)
        let canonicalID = try selected.mdCheckpoint.fingerprint()
        try #require(selected.mdCheckpointFingerprint != canonicalID)
        let publication = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        #expect(!publication.reused)
        let receipt = publication.receipt
        #expect(receipt.implementationFingerprint == f.implementation)
        #expect(try receipt.task.fingerprint() == receipt.taskFingerprint)
        #expect(receipt.provenance.mdCheckpointFingerprint == selected.mdCheckpointFingerprint)
        #expect(receipt.provenance.canonicalCheckpointFingerprint == canonicalID)
        #expect(receipt.provenance.samplingCheckpointFingerprint == f.checkpointID)
        #expect(receipt.provenance.requestFingerprint == f.requestID)
        #expect(receipt.provenance.replicaIndex == 0 && receipt.provenance.replicaIdentifier == "replica-0")
        #expect(receipt.provenance.replicaSeed == 17)
        #expect(receipt.provenance.sourceTermination == "notRecorded")
        #expect(!receipt.provenance.declaredObservableCriteriaSatisfied)
        #expect(receipt.provenance.sourceReceiptFingerprint == nil)
        let expectedKinds = ["checkpoint": "vivo.md-checkpoint", "source-structure": "vivo.molecular-structure-document",
            "system": "vivo.classical-system", "configuration": "vivo.md-configuration", "snapshot": "vivo.md-state-snapshot",
            "structure": "vivo.molecular-structure-document", "frame": "vivo.trajectory-frame",
            "mapping": "vivo.md-snapshot-mapping", "provenance": "vivo.molecular-sampling-selected-state"]
        #expect(Set(receipt.outputs.map(\.name)) == Set(expectedKinds.keys))
        #expect(Set(receipt.payloadFingerprints.keys) == Set(expectedKinds.keys))
        for output in receipt.outputs {
            #expect(output.kind == expectedKinds[output.name])
            let descriptor = try await f.store.descriptor(for: output.artifact)
            #expect(descriptor.kind == "chemistry-output")
            let bytes = try await payload(output.name, publication: publication, store: f.store)
            #expect(receipt.payloadFingerprints[output.name] == (try VivoCanonicalJSON.fingerprint(bytes)))
        }
        #expect(try await payload("checkpoint", publication: publication, store: f.store) == original)
        #expect(try await f.store.data(for: selected.mdCheckpointFingerprint) == original)
        #expect(try await f.store.descriptor(for: selected.mdCheckpointFingerprint) == originalDescriptor)
        #expect(originalDescriptor.kind == "md-checkpoint")
        let source = try VivoCanonicalJSON.decode(VivoMolecularStructureDocument.self,
            from: await payload("source-structure", publication: publication, store: f.store))
        #expect(source.structure == f.request.structure)
        #expect(source.structureFingerprint == f.request.system.structureFingerprint)
        #expect(try VivoCanonicalJSON.decode(VivoClassicalSystem.self,
            from: await payload("system", publication: publication, store: f.store)) == f.request.system)
        #expect(try VivoCanonicalJSON.decode(VivoMDConfiguration.self,
            from: await payload("configuration", publication: publication, store: f.store)) == selected.configuration)
        let snapshot = try VivoCanonicalJSON.decode(VivoMDStateSnapshot.self,
            from: await payload("snapshot", publication: publication, store: f.store))
        #expect(snapshot.systemFingerprint == selected.mdCheckpoint.systemFingerprint)
        #expect(snapshot.configurationFingerprint == selected.mdCheckpoint.configurationFingerprint)
        #expect(snapshot.stepIndex == 6 && snapshot.timePS == 0.25 + 6.0 / 256)
        #expect(snapshot.positionsNM == selected.mdCheckpoint.positionsNM)
        #expect(snapshot.velocitiesNMPerPS == selected.mdCheckpoint.velocitiesNMPerPS)
        #expect(snapshot.periodicCell == selected.mdCheckpoint.periodicCell)
        #expect(snapshot.potentialEnergyKJPerMol == nil && snapshot.kineticEnergyKJPerMol == nil && snapshot.temperatureK == nil)
        let mapping = try VivoCanonicalJSON.decode(VivoWorkflowSnapshotMapping.self,
            from: await payload("mapping", publication: publication, store: f.store))
        let frame = try VivoCanonicalJSON.decode(VivoTrajectoryFrame.self,
            from: await payload("frame", publication: publication, store: f.store))
        let document = try VivoCanonicalJSON.decode(VivoMolecularStructureDocument.self,
            from: await payload("structure", publication: publication, store: f.store))
        #expect(mapping.atomToParticle == [1, 0] && mapping.checkpointPayload == canonicalID)
        #expect(frame.positionsNM == [snapshot.positionsNM[1], snapshot.positionsNM[0]])
        #expect(frame.velocitiesNMPerPS == [snapshot.velocitiesNMPerPS[1], snapshot.velocitiesNMPerPS[0]])
        #expect(frame.step == snapshot.stepIndex && frame.timePS == snapshot.timePS)
        #expect(document.structure.conformers[0].positionsNM == frame.positionsNM)
        #expect(document.sourceFingerprint == canonicalID && document.structureFingerprint == mapping.snapshotStructure)
        #expect(try await payload("provenance", publication: publication, store: f.store) == VivoCanonicalJSON.encode(receipt.provenance))
        let verification = try await VivoMolecularSamplingExporter.verify(receipt, store: f.store,
            implementationFingerprint: f.implementation)
        #expect(verification.exportReceiptFingerprint == publication.artifact.fingerprint)
        #expect(verification.verifiedOutputCount == 9)
        #expect(verification.mdCheckpointFingerprint == selected.mdCheckpointFingerprint)
        #expect(verification.canonicalCheckpointFingerprint == canonicalID)
        #expect(verification.sourceTermination == "notRecorded" && !verification.declaredObservableCriteriaSatisfied)
    }

    @Test func repeatedPublicationReusesTheExactTaskReceiptAndOutputs() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await select(f)
        let first = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        let before = try inventory(f.root)
        let second = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        #expect(second.reused)
        #expect(first.artifact == second.artifact)
        #expect(try VivoCanonicalJSON.encode(first.receipt) == VivoCanonicalJSON.encode(second.receipt))
        #expect(first.receipt.taskFingerprint == second.receipt.taskFingerprint)
        #expect(first.receipt.workflowReceiptFingerprint == second.receipt.workflowReceiptFingerprint)
        #expect(first.receipt.outputs == second.receipt.outputs)
        #expect(try inventory(f.root) == before)
    }

    @Test func cachedExportStillRevalidatesTheSelectedSourceTail() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await select(f)
        let publication = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: selected.trajectoryManifestFingerprint)
        let tail = try await archive.readLink(try #require(archive.manifest.tail))
        try await corrupt(tail.payload, fixture: f)
        let before = try inventory(f.root)
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
                implementationFingerprint: f.implementation)
        }
        #expect(try inventory(f.root) == before)
    }

    @Test func corruptedWrappedOutputCannotVerify() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let publication = try await VivoMolecularSamplingExporter.publish(select(f), implementationFingerprint: f.implementation)
        let frame = try #require(publication.receipt.outputs.first { $0.name == "frame" })
        try await corrupt(frame.artifact, fixture: f)
        let before = try inventory(f.root)
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
                implementationFingerprint: f.implementation)
        }
        #expect(try inventory(f.root) == before)
    }

    @Test func forgedProvenanceAndWrongImplementationCannotVerify() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let publication = try await VivoMolecularSamplingExporter.publish(select(f), implementationFingerprint: f.implementation)
        var json = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(publication.receipt)) as? [String: Any])
        var provenance = try #require(json["provenance"] as? [String: Any])
        provenance["sourceTermination"] = "converged"
        provenance["declaredObservableCriteriaSatisfied"] = true
        json["provenance"] = provenance
        let forged = try VivoCanonicalJSON.decode(VivoMolecularSamplingExportReceipt.self,
            from: JSONSerialization.data(withJSONObject: json, options: [.sortedKeys]))
        let wrongImplementation = try VivoCanonicalJSON.fingerprint(Data("wrong-export-implementation".utf8))
        let alternate = try await VivoMolecularSamplingExporter.publish(select(f), implementationFingerprint: wrongImplementation)
        _ = try await VivoMolecularSamplingExporter.verify(alternate.receipt, store: f.store,
            implementationFingerprint: wrongImplementation)
        var substituted = try #require(JSONSerialization.jsonObject(with: VivoCanonicalJSON.encode(publication.receipt)) as? [String: Any])
        substituted["workflowReceiptFingerprint"] = try JSONSerialization.jsonObject(
            with: VivoCanonicalJSON.encode(alternate.receipt.workflowReceiptFingerprint), options: [.fragmentsAllowed])
        let rebound = try VivoCanonicalJSON.decode(VivoMolecularSamplingExportReceipt.self,
            from: JSONSerialization.data(withJSONObject: substituted, options: [.sortedKeys]))
        let before = try inventory(f.root)
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(forged, store: f.store, implementationFingerprint: f.implementation)
        }
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
                implementationFingerprint: wrongImplementation)
        }
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(rebound, store: f.store,
                implementationFingerprint: f.implementation)
        }
        #expect(try inventory(f.root) == before)
    }

    @Test func publicationUsesTheSelectedStoreAfterItsPathIsReplaced() async throws {
        let f = try await fixture()
        let moved = f.root.deletingLastPathComponent().appendingPathComponent(f.root.lastPathComponent + "-moved")
        defer { try? FileManager.default.removeItem(at: f.root); try? FileManager.default.removeItem(at: moved) }
        let selected = try await select(f)
        try FileManager.default.moveItem(at: f.root, to: moved)
        guard Darwin.mkdir(f.root.path, 0o700) == 0 else { throw POSIXError(POSIXErrorCode(rawValue: errno) ?? .EIO) }
        let replacement = try VivoArtifactStore(rootURL: f.root, createIfNeeded: false)
        let publication = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        #expect(selected.store === f.store)
        #expect(try await f.store.verify(publication.artifact.fingerprint))
        #expect(FileManager.default.fileExists(atPath: moved.appendingPathComponent(publication.artifact.objectPath).path))
        #expect(try await replacement.list().isEmpty)
        #expect(try inventory(f.root).isEmpty)
        _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
    }

    @Test func preCancelledPublicationCannotReturnSuccessOrWriteArtifacts() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await select(f), before = try inventory(f.root)
        let descriptors = try await f.store.list()
        let task = Task {
            withUnsafeCurrentTask { $0?.cancel() }
            return try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        }
        do { _ = try await task.value; Issue.record("pre-cancelled export returned success") } catch is CancellationError {}
        #expect(try await f.store.list() == descriptors)
        #expect(try inventory(f.root) == before)
    }

    @Test func cancelledPartialPrefixRetainsItsReceiptWithoutClaimingConvergence() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let sourceReceipt = VivoMolecularSamplingRunReceipt(schema: VivoMolecularSamplingRunReceipt.schema,
            requestFingerprint: f.requestID, status: .cancelled, completedBlocks: 1, consecutivePasses: 0,
            checkpoint: f.checkpointID, trajectoryManifests: f.cursor.replicas.compactMap(\.trajectory),
            sampling: f.diagnostic, diagnostic: "synthetic cancellation", interpretation: VivoMolecularSampling.interpretation)
        let receiptID = try await put(sourceReceipt, kind: "molecular-sampling-receipt", store: f.store).fingerprint
        let selected = try await select(f, receipt: receiptID)
        let publication = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        #expect(publication.receipt.provenance.sourceReceiptFingerprint == receiptID)
        #expect(publication.receipt.provenance.sourceTermination == "cancelled")
        #expect(!publication.receipt.provenance.declaredObservableCriteriaSatisfied)
        let verification = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
        #expect(verification.sourceTermination == "cancelled" && !verification.declaredObservableCriteriaSatisfied)
        #expect(verification.samplingCheckpointFingerprint == f.checkpointID)
        #expect(try await f.store.data(for: receiptID) == VivoCanonicalJSON.encode(sourceReceipt))
    }

    @Test func strongerVerificationChecksAllPayloadsWithoutRewritingReceiptScope() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let selected = try await select(f)
        try #require(selected.trajectoryManifest.chunkCount == 2)
        #expect(selected.trajectoryValidation.verifiedPayloads == 1)
        let publication = try await VivoMolecularSamplingExporter.publish(selected, implementationFingerprint: f.implementation)
        #expect(publication.receipt.provenance.validationScope == "restart")
        let before = try inventory(f.root)
        let verification = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation, minimumValidation: .allPayloads)
        #expect(verification.validationScope == "allPayloads")
        #expect(verification.recordedValidationScope == "restart")
        #expect(verification.verifiedPayloads == 2)
        #expect(verification.verifiedPayloadBytes > publication.receipt.provenance.verifiedPayloadBytes)
        #expect(verification.exportReceiptFingerprint == publication.artifact.fingerprint)
        #expect(publication.receipt.provenance.validationScope == "restart")
        #expect(try inventory(f.root) == before)
        let archive = try await VivoMDTrajectoryArchiveReader.open(store: f.store, manifest: selected.trajectoryManifestFingerprint)
        let links = try await archive.index()
        let first = try await archive.readLink(try #require(links.first))
        try await corrupt(first.payload, fixture: f)
        _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
        await rejects {
            _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
                implementationFingerprint: f.implementation, minimumValidation: .allPayloads)
        }
    }

    @Test func immutableReceiptVerificationSurvivesDeletedOrRepointedCacheWithoutWrites() async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let publication = try await VivoMolecularSamplingExporter.publish(select(f), implementationFingerprint: f.implementation)
        let reference = "chemistry-task-" + publication.receipt.taskFingerprint.hex
        try await f.store.removeReference(reference)
        let descriptors = try await f.store.list(), before = try inventory(f.root)
        let verification = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
        #expect(verification.exportReceiptFingerprint == publication.artifact.fingerprint)
        #expect(try await f.store.list() == descriptors)
        #expect(try inventory(f.root) == before)
        do { _ = try await f.store.reference(reference); Issue.record("verification republished a missing cache") }
        catch VivoArtifactStoreError.referenceMissing(_) {}
        let replacement = try await f.store.descriptor(for: f.requestID)
        _ = try await f.store.setReference(reference, to: replacement)
        let repointed = try inventory(f.root)
        _ = try await VivoMolecularSamplingExporter.verify(publication.receipt, store: f.store,
            implementationFingerprint: f.implementation)
        #expect(try inventory(f.root) == repointed)
        #expect(try await f.store.reference(reference).artifact == replacement)
    }

    @Test(arguments: ["receipt", "envelope"])
    func workflowVerificationRejectsOversizedWireBeforeReadingOrHashing(_ oversized: String) async throws {
        let f = try await fixture(); defer { try? FileManager.default.removeItem(at: f.root) }
        let expected = ["value": Data("1".utf8)]
        let operation = VivoChemistryOperation(identifier: "test.sampling-export-wire-limit", version: "1",
            implementationFingerprint: f.implementation, outputs: [.init(name: "value", kind: "synthetic-value")],
            execute: { _, _, _ in expected },
            validateOutputs: { _, _, payloads, _ in
                guard payloads == expected else { throw VivoChemistryError.invalid("synthetic workflow output differs") }
            })
        let task = VivoChemistryTask(operation: operation.identifier, version: operation.version,
            implementationFingerprint: f.implementation, inputs: [], configuration: .object([:]), outputs: operation.outputs,
            resources: .init(maximumInputBytes: 16, maximumOutputBytes: 16))
        let workflow = VivoChemistryWorkflow(store: f.store)
        let valid = try await workflow.run(task, using: operation)
        _ = try await workflow.verifyReceipt(valid.receiptFingerprint, task: task, using: operation)
        let wireID: VivoFingerprint, receiptID: VivoFingerprint
        if oversized == "receipt" {
            var data = try await f.store.data(for: valid.receiptFingerprint)
            data.append(Data(repeating: 0x20, count: 128 * 1024))
            wireID = try await f.store.put(data: data, kind: "chemistry-task-receipt",
                mediaType: "application/vnd.numivivo.chemistry-task-receipt+json").fingerprint
            receiptID = wireID
        } else {
            let output = try #require(valid.outputs.first)
            var data = try await f.store.data(for: output.artifact)
            data.append(Data(repeating: 0x20, count: 128 * 1024))
            wireID = try await f.store.put(data: data, kind: "chemistry-output",
                mediaType: "application/vnd.numivivo.chemistry-output+json").fingerprint
            let receipt = WorkflowReceipt(schema: "numivivo.org/chemistry-task-receipt/v1",
                taskFingerprint: valid.taskFingerprint,
                outputs: [.init(name: output.name, kind: output.kind, artifact: wireID)])
            receiptID = try await put(receipt, kind: "chemistry-task-receipt", store: f.store).fingerprint
        }
        // If verification reads/hashes before bounding wire size, this corruption
        // would raise integrityFailure instead of the rooted fstat limit error.
        try await corrupt(wireID, fixture: f)
        let before = try inventory(f.root)
        do {
            _ = try await workflow.verifyReceipt(receiptID, task: task, using: operation)
            Issue.record("oversized workflow wire object was accepted")
        } catch VivoRootedFileStore.Failure.exceededLimit(_) {}
        #expect(try inventory(f.root) == before)
    }
}
