import Foundation
import Darwin
import Testing
@testable import NumiVivoKit

/// Real Metal sampling and immutable-prefix continuation. The tiny, explicitly
/// synthetic harmonic model is a routing/identity fixture, not an equilibrium,
/// molecular accuracy, kinetics or performance qualification.
@Suite(.serialized) struct MolecularSamplingNativeBridgeTests: Sendable {
    private struct ReplicaEvidence: Codable {
        let index: Int
        let seed: UInt64
        let checkpoint: String
        let trajectory: String
        let acceptedStep: UInt64
        let timePS: Double
    }
    private struct Evidence: Codable {
        let schema = "numivivo.org/test-evidence/sampling-prefix/v1"
        let status = "success"
        let nativeStatus = "accepted-metal-sampling-and-exact-prefix-resume"
        let sourceCommit: String
        let numericalContract: String
        let model = "synthetic two-particle harmonic routing fixture; not physical H2 parameters"
        let store: String
        let requestFile: String
        let request: String
        let checkpoint: String
        let sourceReceipt: String
        let reference: String
        let resumedFrom: String
        let completedBlocks: Int
        let sourceTermination: String
        let declaredObservableCriteriaSatisfied: Bool
        let exactResumedCursor: Bool
        let replicas: [ReplicaEvidence]
    }

    @Test func acceptedReplicaPrefixResumesExactlyAndRetainsExportInputs() async throws {
        let environment = ProcessInfo.processInfo.environment
        let retained = environment["NUMIVIVO_TEST_ARTIFACTS"]
        let parent: URL
        if let retained {
            guard !retained.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
                throw VivoArtifactValidationError.invalid("NUMIVIVO_TEST_ARTIFACTS must name a directory")
            }
            parent = URL(fileURLWithPath: retained, isDirectory: true).standardizedFileURL
            try FileManager.default.createDirectory(at: parent, withIntermediateDirectories: true)
        } else {
            parent = FileManager.default.temporaryDirectory
        }
        let root = parent.appendingPathComponent("molecular-sampling-prefix-\(UUID().uuidString)")
        guard Darwin.mkdir(root.path, 0o700) == 0 else {
            throw NSError(domain: NSPOSIXErrorDomain, code: Int(errno), userInfo: [NSFilePathErrorKey: root.path])
        }
        defer { if retained == nil { try? FileManager.default.removeItem(at: root) } }
        print("NUMIVIVO_SAMPLING_PREFIX_ROOT=\(root.path)")
        let storeRoot = root.appendingPathComponent("store")
        let store = try VivoArtifactStore(rootURL: storeRoot)
        let hydrogen = try #require(VivoElement.from(symbol: "H"))
        let positions: [VivoVector3D] = [.zero, .init(0.125, 0, 0)]
        let structure = VivoMolecularStructure(identifier: "native-sampling-bridge-harmonic-fixture", atoms: [
            .init(index: 0, name: "H1", element: hydrogen), .init(index: 1, name: "H2", element: hydrogen)
        ], bonds: [.init(atomA: 0, atomB: 1)], conformers: [.init(positionsNM: positions)])
        let system = VivoClassicalSystem(identifier: "synthetic-harmonic-routing-parameters",
            structureFingerprint: try VivoStructureCodec.fingerprint(structure), particles: [
                .init(index: 0, atomIndex: 0, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0),
                .init(index: 1, atomIndex: 1, typeIdentifier: "H", massDa: 1, chargeE: 0, sigmaNM: 0, epsilonKJPerMol: 0)
            ], bonds: [.init(a: 0, b: 1, lengthNM: 0.125, forceConstant: 1000)])
        let systemID = try system.fingerprint()
        let initialStates = [0.125, 0.375].map {
            VivoClassicalInitialState(systemFingerprint: systemID, positionsNM: positions, sourceTimePS: $0)
        }
        let md = VivoMDConfiguration(timeStepPS: 1.0 / 1024, cutoffNM: 0.5, neighborSkinNM: 0,
            electrostatics: .cutoff, ensemble: .nvt, thermostat: .langevinMiddle,
            targetTemperatureK: 300, frictionPerPS: 1, neighborListEnabled: false)
        let request = VivoMolecularSamplingRunRequest(structure: structure, system: system,
            initialStates: initialStates, replicaSeeds: [17, 29], md: md,
            contextIdentifier: "native accepted-prefix routing fixture; declared criteria intentionally unmet",
            observables: [.init(identifier: "distance", kind: .distance(atomA: 0, atomB: 1), maximumMeanStandardError: 100)],
            convergence: .init(minimumRetainedFramesPerReplica: 64, minimumReplicas: 2,
                maximumRHat: 1.1, minimumEffectiveSamplesPerReplica: 20, maximumAutocorrelationLag: 64),
            equilibrationSteps: 4, stepsPerBlock: 8, sampleEvery: 4, maximumBlocks: 4,
            requiredConsecutivePasses: 2, trajectoryChunkBytes: 1024 * 1024)
        let requestFile = root.appendingPathComponent("request.json")
        try VivoCanonicalJSON.encode(request).write(to: requestFile, options: .withoutOverwriting)
        let complete = try await VivoMolecularSamplingRunner.run(request, store: store)
        try VivoCanonicalJSON.encode(complete).write(to: root.appendingPathComponent("native-receipt.json"), options: .withoutOverwriting)
        try #require(complete.status == .budgetExhausted, "native producer failed: \(complete.diagnostic ?? "no diagnostic")")
        try #require(complete.completedBlocks == 4 && complete.consecutivePasses == 0)
        try #require(complete.sampling?.converged == false)
        let receiptArtifact = try await store.put(data: VivoCanonicalJSON.encode(complete),
            kind: "molecular-sampling-receipt", mediaType: "application/json")
        let reference = VivoMolecularSamplingRunner.checkpointReferenceName(requestFingerprint: complete.requestFingerprint)
        try #require(try await store.reference(reference).artifact.fingerprint == complete.checkpoint)

        // Select an immutable B1 cursor actually published by this native run.
        // No synthetic cursor, timing race or mutable-reference rewind is used.
        let cursors = try await store.list(kind: "molecular-sampling-checkpoint")
        var firstBlock: VivoFingerprint?
        var zeroBlock: VivoFingerprint?
        for artifact in cursors {
            let cursor = try VivoCanonicalJSON.decode(VivoMolecularSamplingCursor.self,
                from: await store.data(for: artifact.fingerprint, maximumBytes: 1024 * 1024, verify: true))
            try #require(cursor.requestFingerprint == complete.requestFingerprint)
            if cursor.completedBlocks == 1 { try #require(firstBlock == nil); firstBlock = artifact.fingerprint }
            if cursor.completedBlocks == 0 { try #require(zeroBlock == nil); zeroBlock = artifact.fingerprint }
        }
        try #require(cursors.count == 5)
        let b0 = try await VivoMolecularSamplingArchiveReader.open(store: store, checkpoint: #require(zeroBlock))
        try #require(b0.inspect().completedBlocks == 0)
        var emptyRejected = false
        do { _ = try await b0.selectReplica(index: 0) } catch { emptyRejected = true }
        try #require(emptyRejected)
        let b1 = try #require(firstBlock)
        let resumed = try await VivoMolecularSamplingRunner.run(request, store: store, resumeFrom: b1)
        try VivoCanonicalJSON.encode(resumed).write(to: root.appendingPathComponent("resumed-receipt.json"), options: .withoutOverwriting)
        try #require(resumed == complete, "B1 continuation must reproduce the complete accepted cursor, diagnostics and manifests")
        let terminal = try await VivoMolecularSamplingRunner.run(request, store: store, resumeFrom: complete.checkpoint)
        try #require(terminal == complete)
        try #require(try await store.list(kind: "molecular-sampling-checkpoint").count == 5)
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: store, checkpoint: complete.checkpoint)
        try #require(!reader.inspect().declaredObservableCriteriaSatisfied)
        var replicas: [ReplicaEvidence] = []
        for index in initialStates.indices {
            let selected = try await reader.selectReplica(index: index, validation: .allPayloads, receipt: receiptArtifact.fingerprint)
            let checkpoint = selected.mdCheckpoint
            try #require(selected.sourceTermination == .budgetExhausted)
            try #require(selected.replicaSeed == request.replicaSeeds[index])
            try #require(selected.configuration.randomSeed == request.replicaSeeds[index])
            try #require(checkpoint.acceptedStep == 36)
            try #require(checkpoint.timePS == initialStates[index].sourceTimePS! + 36.0 / 1024)
            try #require(checkpoint.numericalContract == VivoMDExecutionIdentity.current)
            try #require(checkpoint.velocitiesNMPerPS.contains { $0 != .zero })
            try #require(selected.trajectoryManifest.frameCount == 8)
            try #require(selected.trajectoryValidation.scope == .allPayloads && selected.trajectoryValidation.verifiedPayloads == 4)
            try #require(try await store.data(for: selected.mdCheckpointFingerprint) == selected.mdCheckpointData)
            try selected.mdCheckpointData.write(to: root.appendingPathComponent("replica-\(index)-checkpoint.json"), options: .withoutOverwriting)
            replicas.append(.init(index: index, seed: selected.replicaSeed, checkpoint: selected.mdCheckpointFingerprint.hex,
                trajectory: selected.trajectoryManifestFingerprint.hex, acceptedStep: checkpoint.acceptedStep, timePS: checkpoint.timePS))
        }
        try #require(replicas[0].checkpoint != replicas[1].checkpoint)
        let evidence = Evidence(sourceCommit: environment["NUMIVIVO_TEST_SOURCE_COMMIT"] ?? "unrecorded",
            numericalContract: VivoMDExecutionIdentity.current, store: storeRoot.path, requestFile: requestFile.path,
            request: complete.requestFingerprint.hex, checkpoint: complete.checkpoint.hex,
            sourceReceipt: receiptArtifact.fingerprint.hex, reference: reference, resumedFrom: b1.hex,
            completedBlocks: complete.completedBlocks, sourceTermination: complete.status.rawValue,
            declaredObservableCriteriaSatisfied: false, exactResumedCursor: true, replicas: replicas)
        try VivoCanonicalJSON.encode(evidence).write(to: root.appendingPathComponent("sampling-prefix-receipt.json"), options: .withoutOverwriting)
    }
}
