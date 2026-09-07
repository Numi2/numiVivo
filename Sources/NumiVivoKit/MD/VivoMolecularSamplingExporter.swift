import Foundation

/// Lineage of one accepted replica at a completed cross-replica block boundary.
/// This is a selection record, not a new convergence or kinetic qualification.
public struct VivoMolecularSamplingExportProvenance: Codable, Sendable, Equatable {
    public let schema: String
    public let requestFingerprint: VivoFingerprint
    public let samplingCheckpointFingerprint: VivoFingerprint
    public let replicaIndex: Int
    public let replicaIdentifier: String
    public let replicaSeed: UInt64
    /// Hash of the original stored checkpoint bytes, including their formatting.
    public let mdCheckpointFingerprint: VivoFingerprint
    /// Hash used by the existing snapshot mapper for the canonical decoded state.
    public let canonicalCheckpointFingerprint: VivoFingerprint
    public let systemFingerprint: VivoFingerprint
    public let sourceStructureFingerprint: VivoFingerprint
    public let configurationFingerprint: VivoFingerprint
    public let numericalContract: String
    public let acceptedStep: UInt64
    public let timePS: Double
    public let periodicCell: VivoPeriodicCell?
    public let trajectoryManifestFingerprint: VivoFingerprint
    public let trajectoryFrameCount: UInt64
    public let trajectoryChunkCount: UInt64
    public let validationScope: String
    public let indexedChunks: UInt64
    public let verifiedPayloads: UInt64
    public let verifiedPayloadBytes: UInt64
    public let completedBlocks: Int
    public let consecutivePasses: Int
    public let requiredConsecutivePasses: Int
    public let declaredObservableCriteriaSatisfied: Bool
    public let diagnosticEvaluations: Int
    public let latestDiagnosticFingerprint: VivoFingerprint?
    public let latestDiagnostic: VivoMolecularSamplingResult?
    public let criteriaHistoryScope: String
    public let sourceReceiptFingerprint: VivoFingerprint?
    public let sourceTermination: String
    /// Producer reference retained without revalidating minimization history.
    public let minimizationFingerprint: VivoFingerprint?
    public let interpretation: String
}

/// Immutable output identities. Task configuration binds the pre-task provenance
/// digest; this final record binds task and output artifacts without a hash cycle.
public struct VivoMolecularSamplingExportReceipt: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/molecular-sampling-export/v1"
    public let schema: String
    public let implementationFingerprint: VivoFingerprint
    public let task: VivoChemistryTask
    public let taskFingerprint: VivoFingerprint
    public let workflowReceiptFingerprint: VivoFingerprint
    public let provenance: VivoMolecularSamplingExportProvenance
    public let outputs: [VivoChemistryOutputReceipt]
    /// Digests of decoded payload bytes, distinct from chemistry-output envelopes.
    public let payloadFingerprints: [String: VivoFingerprint]

    public func fingerprint() throws -> VivoFingerprint {
        try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoMolecularSamplingExportPublication: Sendable {
    public let receipt: VivoMolecularSamplingExportReceipt
    public let artifact: VivoStoredArtifact
    /// Cache reuse is an observation of this call, not part of immutable lineage.
    public let reused: Bool
}

public struct VivoMolecularSamplingExportVerification: Codable, Sendable, Equatable {
    public let exportReceiptFingerprint: VivoFingerprint
    public let taskFingerprint: VivoFingerprint
    public let samplingCheckpointFingerprint: VivoFingerprint
    public let mdCheckpointFingerprint: VivoFingerprint
    public let canonicalCheckpointFingerprint: VivoFingerprint
    public let sourceTermination: String
    public let declaredObservableCriteriaSatisfied: Bool
    public let recordedValidationScope: String
    /// Effective scope freshly checked, possibly stronger than the export scope.
    public let validationScope: String
    public let indexedChunks: UInt64
    public let verifiedPayloads: UInt64
    public let verifiedPayloadBytes: UInt64
    public let verifiedOutputCount: Int
}

public enum VivoMolecularSamplingExporter {
    public static let operationIdentifier = "vivo.platform.molecular-sampling-export"
    public static let operationVersion = "1"
    public static let outputs: [VivoChemistryTaskOutput] = [
        .init(name: "checkpoint", kind: "vivo.md-checkpoint"),
        .init(name: "source-structure", kind: "vivo.molecular-structure-document"),
        .init(name: "system", kind: "vivo.classical-system"),
        .init(name: "configuration", kind: "vivo.md-configuration"),
        .init(name: "snapshot", kind: "vivo.md-state-snapshot"),
        .init(name: "structure", kind: "vivo.molecular-structure-document"),
        .init(name: "frame", kind: "vivo.trajectory-frame"),
        .init(name: "mapping", kind: "vivo.md-snapshot-mapping"),
        .init(name: "provenance", kind: "vivo.molecular-sampling-selected-state")
    ].sorted { $0.name < $1.name }

    /// Publish through the chemistry workflow's existing output envelopes. The
    /// selection retains the exact rooted store actor; no alternate store or path
    /// can be supplied. Its evidence describes the earlier successful selection.
    public static func publish(_ selection: VivoValidatedMolecularSamplingSelection,
                               implementationFingerprint: VivoFingerprint,
                               budget: VivoChemistryBudget = .init()) async throws -> VivoMolecularSamplingExportPublication {
        try Task.checkCancellation()
        let prepared = try prepare(selection, implementationFingerprint: implementationFingerprint, budget: budget)
        let workflow = VivoChemistryWorkflow(store: selection.store)
        let result = try await workflow.run(prepared.task, using: operation(prepared))
        try Task.checkCancellation()
        let receipt = try makeReceipt(prepared, result: result)
        let artifact = try await selection.store.put(data: VivoCanonicalJSON.encode(receipt),
            kind: "molecular-sampling-export", mediaType: "application/json")
        try Task.checkCancellation()
        return .init(receipt: receipt, artifact: artifact, reused: result.reused)
    }

    /// Freshly validate the source sampling prefix and reconstruct every payload,
    /// then verify the exact immutable workflow receipt. Mutable cache references
    /// are neither read nor written, and missing artifacts are never regenerated.
    /// Ordinary chemistry receipt verification alone does not perform this graph
    /// traversal or reconstruct sampling diagnostics.
    ///
    /// Sampling limits bound one effective selected traversal. Strengthening a
    /// restart receipt additionally reads one manifest and tail link, each capped
    /// at 64 KiB, to reconstruct its recorded payload-byte claim. Workflow reads
    /// have separate wire bounds derived from the persisted task resources; pure
    /// reconstruction is admitted by that task's chemistry budget.
    public static func verify(_ receipt: VivoMolecularSamplingExportReceipt,
                              store: VivoArtifactStore,
                              implementationFingerprint: VivoFingerprint,
                              limits: VivoMolecularSamplingReadLimits = .init(),
                              minimumValidation: VivoMolecularSamplingSelectionValidation? = nil) async throws -> VivoMolecularSamplingExportVerification {
        try Task.checkCancellation()
        guard receipt.schema == VivoMolecularSamplingExportReceipt.schemaID,
              receipt.implementationFingerprint == implementationFingerprint,
              let recorded = VivoMolecularSamplingSelectionValidation(rawValue: receipt.provenance.validationScope) else {
            throw invalid("receipt schema, implementation or validation scope")
        }
        let effective: VivoMolecularSamplingSelectionValidation =
            recorded == .allPayloads || minimumValidation == .allPayloads ? .allPayloads : .restart
        let reader = try await VivoMolecularSamplingArchiveReader.open(store: store,
            checkpoint: receipt.provenance.samplingCheckpointFingerprint, limits: limits)
        let selection = try await reader.selectReplica(index: receipt.provenance.replicaIndex,
            validation: effective, receipt: receipt.provenance.sourceReceiptFingerprint)
        var recordedTailBytes: UInt64?
        if recorded != effective {
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: store,
                manifest: selection.trajectoryManifestFingerprint)
            guard let tail = archive.manifest.tail else { throw invalid("accepted trajectory has no tail") }
            recordedTailBytes = try await archive.readLink(tail).payloadBytes
        }
        let prepared = try prepare(selection, implementationFingerprint: implementationFingerprint,
            budget: receipt.task.resources.budget, recordedScope: recorded, recordedTailBytes: recordedTailBytes)
        guard receipt.task == prepared.task, receipt.taskFingerprint == (try prepared.task.fingerprint()),
              receipt.provenance == prepared.provenance,
              receipt.payloadFingerprints == (try payloadFingerprints(prepared.payloads)) else {
            throw invalid("reconstructed selection, task or payload identity differs")
        }
        let workflow = VivoChemistryWorkflow(store: store)
        let result = try await workflow.verifyReceipt(receipt.workflowReceiptFingerprint,
            task: prepared.task, using: operation(prepared))
        guard receipt == (try makeReceipt(prepared, result: result)) else {
            throw invalid("immutable workflow receipt or output identity differs")
        }
        try Task.checkCancellation()
        return .init(exportReceiptFingerprint: try receipt.fingerprint(), taskFingerprint: result.taskFingerprint,
            samplingCheckpointFingerprint: selection.samplingCheckpointFingerprint,
            mdCheckpointFingerprint: selection.mdCheckpointFingerprint,
            canonicalCheckpointFingerprint: prepared.provenance.canonicalCheckpointFingerprint,
            sourceTermination: selection.sourceTermination.rawValue,
            declaredObservableCriteriaSatisfied: selection.inspection.declaredObservableCriteriaSatisfied,
            recordedValidationScope: recorded.rawValue, validationScope: effective.rawValue,
            indexedChunks: selection.trajectoryValidation.indexedChunks,
            verifiedPayloads: selection.trajectoryValidation.verifiedPayloads,
            verifiedPayloadBytes: selection.trajectoryValidation.verifiedPayloadBytes,
            verifiedOutputCount: result.outputs.count)
    }

    private struct Prepared: Sendable {
        let task: VivoChemistryTask
        let provenance: VivoMolecularSamplingExportProvenance
        let payloads: [String: Data]
    }

    private static func prepare(_ selection: VivoValidatedMolecularSamplingSelection,
                                implementationFingerprint: VivoFingerprint, budget: VivoChemistryBudget,
                                recordedScope: VivoMolecularSamplingSelectionValidation? = nil,
                                recordedTailBytes: UInt64? = nil) throws -> Prepared {
        try Task.checkCancellation(); try budget.validate()
        let checkpoint = selection.mdCheckpoint, request = selection.request, inspection = selection.inspection
        _ = try budget.elements([request.system.particles.count, 3], simultaneousArrays: 16)
        guard try VivoCanonicalJSON.fingerprint(selection.mdCheckpointData) == selection.mdCheckpointFingerprint,
              checkpoint.numericalContract == VivoMDExecutionIdentity.current else {
            throw invalid("captured checkpoint identity")
        }
        let source = try VivoMolecularStructureDocument(structure: request.structure)
        let mapped = try VivoPlatformSnapshotOperations.map(source: source, system: request.system,
            checkpoint: checkpoint, budget: budget)
        let snapshot = VivoMDStateSnapshot(systemFingerprint: checkpoint.systemFingerprint,
            configurationFingerprint: checkpoint.configurationFingerprint, stepIndex: checkpoint.acceptedStep,
            timePS: checkpoint.timePS, positionsNM: checkpoint.positionsNM,
            velocitiesNMPerPS: checkpoint.velocitiesNMPerPS, periodicCell: checkpoint.periodicCell)
        try snapshot.validate(particleCount: request.system.particles.count)
        let validation = selection.trajectoryValidation
        guard let actualScope = VivoMolecularSamplingSelectionValidation(rawValue: validation.scope.rawValue) else {
            throw invalid("selection must include accepted coordinate payload")
        }
        let scope = recordedScope ?? actualScope
        let downgrading = scope == .restart && actualScope == .allPayloads
        guard scope == actualScope || (downgrading && recordedTailBytes != nil) else {
            throw invalid("recorded scope cannot exceed fresh evidence")
        }
        let provenance = VivoMolecularSamplingExportProvenance(
            schema: "numivivo.org/molecular-sampling-selected-state/v1",
            requestFingerprint: selection.requestFingerprint,
            samplingCheckpointFingerprint: selection.samplingCheckpointFingerprint,
            replicaIndex: selection.replicaIndex, replicaIdentifier: selection.replicaIdentifier,
            replicaSeed: selection.replicaSeed, mdCheckpointFingerprint: selection.mdCheckpointFingerprint,
            canonicalCheckpointFingerprint: mapped.mapping.checkpointPayload,
            systemFingerprint: checkpoint.systemFingerprint, sourceStructureFingerprint: source.structureFingerprint,
            configurationFingerprint: checkpoint.configurationFingerprint, numericalContract: VivoMDExecutionIdentity.current,
            acceptedStep: checkpoint.acceptedStep, timePS: checkpoint.timePS, periodicCell: checkpoint.periodicCell,
            trajectoryManifestFingerprint: selection.trajectoryManifestFingerprint,
            trajectoryFrameCount: selection.trajectoryManifest.frameCount,
            trajectoryChunkCount: selection.trajectoryManifest.chunkCount, validationScope: scope.rawValue,
            indexedChunks: validation.indexedChunks, verifiedPayloads: downgrading ? 1 : validation.verifiedPayloads,
            verifiedPayloadBytes: downgrading ? recordedTailBytes! : validation.verifiedPayloadBytes,
            completedBlocks: inspection.completedBlocks, consecutivePasses: inspection.consecutivePasses,
            requiredConsecutivePasses: inspection.requiredConsecutivePasses,
            declaredObservableCriteriaSatisfied: inspection.declaredObservableCriteriaSatisfied,
            diagnosticEvaluations: inspection.diagnosticEvaluations,
            latestDiagnosticFingerprint: inspection.latestDiagnosticFingerprint, latestDiagnostic: inspection.latestDiagnostic,
            criteriaHistoryScope: inspection.criteriaHistoryScope,
            sourceReceiptFingerprint: selection.sourceReceiptFingerprint, sourceTermination: selection.sourceTermination.rawValue,
            minimizationFingerprint: selection.minimizationFingerprint,
            interpretation: "selected accepted replica at the last complete cross-replica block; declared scalar criteria and recorded termination retained separately; no new equilibration, convergence, minimization-history, kinetic or electronic claim")
        let provenanceData = try VivoCanonicalJSON.encode(provenance)
        let payloads = [
            "checkpoint": selection.mdCheckpointData,
            "source-structure": try source.canonicalData(),
            "system": try VivoCanonicalJSON.encode(request.system),
            "configuration": try VivoCanonicalJSON.encode(selection.configuration),
            "snapshot": try VivoCanonicalJSON.encode(snapshot),
            "structure": try mapped.document.canonicalData(),
            "frame": try VivoCanonicalJSON.encode(mapped.frame),
            "mapping": try VivoCanonicalJSON.encode(mapped.mapping),
            "provenance": provenanceData
        ]
        var remaining = budget.maximumBytes
        for payload in payloads.values {
            guard payload.count <= remaining else { throw VivoChemistryError.resourceLimit("sampling export output byte budget") }
            remaining -= payload.count
        }
        var inputs: [VivoChemistryTaskInput] = [
            .init(name: "request", artifact: selection.requestFingerprint, kind: "molecular-sampling-run"),
            .init(name: "sampling-checkpoint", artifact: selection.samplingCheckpointFingerprint, kind: "molecular-sampling-checkpoint"),
            .init(name: "checkpoint", artifact: selection.mdCheckpointFingerprint, kind: "md-checkpoint")
        ]
        if let receipt = selection.sourceReceiptFingerprint {
            inputs.append(.init(name: "source-receipt", artifact: receipt, kind: "molecular-sampling-receipt"))
        }
        // A digest avoids UInt64 seed/step loss through the general JSON value
        // model's signed-integer/Double representation. The full typed lineage
        // remains an output payload and part of this export receipt.
        let task = VivoChemistryTask(operation: operationIdentifier, version: operationVersion,
            implementationFingerprint: implementationFingerprint, inputs: inputs,
            configuration: .object(["provenanceFingerprint": .string(try VivoCanonicalJSON.fingerprint(provenanceData).hex)]),
            outputs: outputs, resources: .init(budget: budget,
                maximumInputBytes: budget.maximumBytes, maximumOutputBytes: budget.maximumBytes))
        _ = try task.fingerprint(); try Task.checkCancellation()
        return Prepared(task: task, provenance: provenance, payloads: payloads)
    }

    private static func operation(_ prepared: Prepared) -> VivoChemistryOperation {
        let validate: @Sendable (VivoJSONValue, [String: Data], [String: Data], VivoChemistryBudget) throws -> Void = {
            configuration, inputs, payloads, budget in
            try Task.checkCancellation()
            guard configuration == prepared.task.configuration, budget == prepared.task.resources.budget,
                  Set(inputs.keys) == Set(prepared.task.inputs.map(\.name)), payloads == prepared.payloads else {
                throw invalid("workflow configuration, input names or reconstructed payloads differ")
            }
            for input in prepared.task.inputs {
                guard let data = inputs[input.name], try VivoCanonicalJSON.fingerprint(data) == input.artifact else {
                    throw invalid("source input byte identity differs")
                }
            }
        }
        return VivoChemistryOperation(identifier: operationIdentifier, version: operationVersion,
            implementationFingerprint: prepared.task.implementationFingerprint, outputs: outputs,
            execute: { configuration, inputs, budget in
                try validate(configuration, inputs, prepared.payloads, budget)
                return prepared.payloads
            }, validateOutputs: validate)
    }

    private static func payloadFingerprints(_ payloads: [String: Data]) throws -> [String: VivoFingerprint] {
        try payloads.mapValues { try VivoCanonicalJSON.fingerprint($0) }
    }

    private static func makeReceipt(_ prepared: Prepared, result: VivoChemistryTaskResult) throws -> VivoMolecularSamplingExportReceipt {
        .init(schema: VivoMolecularSamplingExportReceipt.schemaID,
            implementationFingerprint: prepared.task.implementationFingerprint,
            task: prepared.task, taskFingerprint: result.taskFingerprint,
            workflowReceiptFingerprint: result.receiptFingerprint, provenance: prepared.provenance,
            outputs: result.outputs, payloadFingerprints: try payloadFingerprints(prepared.payloads))
    }

    private static func invalid(_ message: String) -> VivoChemistryError {
        .invalid("molecular sampling export: " + message)
    }
}
