import Foundation

/// A seed request location, never an already qualified point or Hessian.
public enum VivoSamplingReactionSeedTarget: Codable, Sendable, Equatable, Hashable {
    case qualification
    case connectedSaddle
    case connectedEndpoint(identifier: String, componentIndex: Int)
}

public struct VivoSamplingReactionSeedAssignment: Codable, Sendable, Equatable {
    public var target: VivoSamplingReactionSeedTarget
    /// Source chemical atom index for each destination nucleus, in nucleus order.
    /// Each assignment covers its entire destination without duplicate atoms.
    public var sourceAtomByNucleus: [Int]
    public init(target: VivoSamplingReactionSeedTarget, sourceAtomByNucleus: [Int]) {
        self.target = target; self.sourceAtomByNucleus = sourceAtomByNucleus
    }
}

/// Explicit transfer between separately declared classical and electronic models.
/// Sampling convergence is neither required nor promoted to a QM ensemble claim.
public struct VivoSamplingReactionSeedRequest: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/sampling-reaction-seed-request/v1"
    public var schema: String
    public var sourceExport: VivoMolecularSamplingExportReceipt
    public var destination: VivoReactionCalculationRequest
    public var assignments: [VivoSamplingReactionSeedAssignment]
    public var modelTransferStatement: String
    public init(sourceExport: VivoMolecularSamplingExportReceipt, destination: VivoReactionCalculationRequest,
                assignments: [VivoSamplingReactionSeedAssignment], modelTransferStatement: String) {
        schema = Self.schemaID; self.sourceExport = sourceExport; self.destination = destination
        self.assignments = assignments; self.modelTransferStatement = modelTransferStatement
    }
}

public struct VivoSamplingReactionSeedTransfer: Codable, Sendable, Equatable {
    public let assignment: VivoSamplingReactionSeedAssignment
    public let sourceParticleByNucleus: [UInt32]
    public let sourceMassesDa: [Double]
    public let destinationMassesDa: [Double]
    public let destinationThermochemistry: VivoThermochemistryConfiguration
    public let destinationSolver: VivoNuclearSolver
    public let destinationBasisFingerprint: VivoFingerprint
    public let initialQualificationFingerprint: VivoFingerprint
    public let assembledQualificationFingerprint: VivoFingerprint
}

public struct VivoSamplingReactionSeedProvenance: Codable, Sendable, Equatable {
    public let schema: String
    public let assemblyRequestFingerprint: VivoFingerprint
    public let sourceExportReceiptFingerprint: VivoFingerprint
    public let source: VivoMolecularSamplingExportProvenance
    public let sourceConfiguration: VivoMDConfiguration
    public let initialReactionRequestFingerprint: VivoFingerprint
    public let assembledReactionRequestFingerprint: VivoFingerprint
    public let transfers: [VivoSamplingReactionSeedTransfer]
    public let modelTransferStatement: String
    public let interpretation: String
}

public struct VivoSamplingReactionSeedReceipt: Codable, Sendable, Equatable {
    public static let schemaID = "numivivo.org/sampling-reaction-seed/v1"
    public let schema: String
    public let implementationFingerprint: VivoFingerprint
    public let request: VivoSamplingReactionSeedRequest
    public let task: VivoChemistryTask
    public let taskFingerprint: VivoFingerprint
    public let workflowReceiptFingerprint: VivoFingerprint
    public let provenance: VivoSamplingReactionSeedProvenance
    public let outputs: [VivoChemistryOutputReceipt]
    public let payloadFingerprints: [String: VivoFingerprint]
    public func fingerprint() throws -> VivoFingerprint { try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self)) }
}

public struct VivoSamplingReactionSeedPublication: Sendable {
    public let receipt: VivoSamplingReactionSeedReceipt
    public let artifact: VivoStoredArtifact
    public let reused: Bool
}

/// A report only. Reusing it as a source of authority is not supported.
public struct VivoSamplingReactionSeedVerification: Codable, Sendable, Equatable {
    public let receiptFingerprint: VivoFingerprint
    public let taskFingerprint: VivoFingerprint
    public let assembledReactionRequestFingerprint: VivoFingerprint
    public let sourceVerification: VivoMolecularSamplingExportVerification
    public let verifiedOutputCount: Int
}

public enum VivoSamplingReactionSeed {
    public static let operationIdentifier = "vivo.platform.sampling-reaction-seed"
    public static let operationVersion = "1"
    public static let interpretation = "accepted classical coordinates used only as initial nuclear qualification seeds; explicit atom mapping and source/destination mass, temperature and model contexts retained; every destination model, electron count, basis, optimizer and thermochemistry setting is preserved; fresh stationary-point and reaction qualification remains required; no sampled QM ensemble, equilibrium population, barrier or kinetic claim"
    public static let outputs: [VivoChemistryTaskOutput] = [
        .init(name: "provenance", kind: "vivo.sampling-reaction-seed-provenance"),
        .init(name: "request", kind: "vivo.reaction-calculation-request")
    ]

    /// Coordinate-only host assembly; the returned request still requires the
    /// existing nuclear or connected-reaction qualification operation.
    public static func assemble(_ request: VivoSamplingReactionSeedRequest,
                                source: VivoVerifiedMolecularSamplingExport) throws -> VivoReactionCalculationRequest {
        try assembleAndDescribe(request, source: source).request
    }

    /// The validated source owns the rooted store. No alternate store/path may
    /// redirect publication, even if the original filesystem path is replaced.
    public static func publish(_ request: VivoSamplingReactionSeedRequest,
                               source: VivoVerifiedMolecularSamplingExport,
                               implementationFingerprint: VivoFingerprint) async throws -> VivoSamplingReactionSeedPublication {
        try Task.checkCancellation()
        let prepared = try prepare(request, source: source, implementationFingerprint: implementationFingerprint)
        let store = source.selection.store
        for input in prepared.task.inputs {
            guard let bytes = prepared.inputs[input.name] else { throw invalid("missing prepared input") }
            let artifact = try await store.put(data: bytes, kind: input.kind, mediaType: "application/json")
            guard artifact.fingerprint == input.artifact else { throw invalid("stored input identity differs") }
        }
        let result = try await VivoChemistryWorkflow(store: store).run(prepared.task, using: operation(prepared))
        try Task.checkCancellation()
        let receipt = try makeReceipt(prepared, result: result)
        let artifact = try await store.put(data: VivoCanonicalJSON.encode(receipt),
            kind: "sampling-reaction-seed", mediaType: "application/json")
        try Task.checkCancellation()
        return .init(receipt: receipt, artifact: artifact, reused: result.reused)
    }

    /// Reopens the source archive and exact export, reconstructs coordinates and
    /// every context, and verifies the immutable workflow receipt. It never reads
    /// mutable cache references or regenerates a missing input/output artifact.
    /// An omitted sampling implementation means the same implementation as this
    /// assembly; an older retained export requires its explicit expected identity.
    public static func verify(_ receipt: VivoSamplingReactionSeedReceipt, store: VivoArtifactStore,
                              implementationFingerprint: VivoFingerprint,
                              samplingImplementationFingerprint: VivoFingerprint? = nil,
                              limits: VivoMolecularSamplingReadLimits = .init(),
                              minimumValidation: VivoMolecularSamplingSelectionValidation? = nil) async throws -> VivoSamplingReactionSeedVerification {
        try Task.checkCancellation()
        guard receipt.schema == VivoSamplingReactionSeedReceipt.schemaID,
              receipt.implementationFingerprint == implementationFingerprint else { throw invalid("receipt schema or implementation") }
        let source = try await VivoMolecularSamplingExporter.verifiedExport(receipt.request.sourceExport, store: store,
            implementationFingerprint: samplingImplementationFingerprint ?? implementationFingerprint,
            limits: limits, minimumValidation: minimumValidation)
        let prepared = try prepare(receipt.request, source: source, implementationFingerprint: implementationFingerprint)
        guard receipt.task == prepared.task, receipt.taskFingerprint == (try prepared.task.fingerprint()),
              receipt.provenance == prepared.provenance,
              receipt.payloadFingerprints == (try payloadFingerprints(prepared.payloads)) else { throw invalid("reconstructed seed identity differs") }
        let result = try await VivoChemistryWorkflow(store: store).verifyReceipt(receipt.workflowReceiptFingerprint,
            task: prepared.task, using: operation(prepared))
        guard receipt == (try makeReceipt(prepared, result: result)) else { throw invalid("immutable receipt or output identity differs") }
        try Task.checkCancellation()
        return .init(receiptFingerprint: try receipt.fingerprint(), taskFingerprint: result.taskFingerprint,
            assembledReactionRequestFingerprint: prepared.provenance.assembledReactionRequestFingerprint,
            sourceVerification: source.verification, verifiedOutputCount: result.outputs.count)
    }

    private static func assembleAndDescribe(_ request: VivoSamplingReactionSeedRequest,
                                           source: VivoVerifiedMolecularSamplingExport) throws
        -> (request: VivoReactionCalculationRequest, transfers: [VivoSamplingReactionSeedTransfer]) {
        try Task.checkCancellation()
        guard request.schema == VivoSamplingReactionSeedRequest.schemaID,
              request.sourceExport == source.receipt,
              (1...33).contains(request.assignments.count),
              Set(request.assignments.map(\.target)).count == request.assignments.count,
              !request.modelTransferStatement.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty,
              request.modelTransferStatement.utf8.count <= 8192,
              !request.modelTransferStatement.unicodeScalars.contains(where: { CharacterSet.controlCharacters.contains($0) && $0 != "\n" && $0 != "\t" }) else {
            throw invalid("schema, exact verified source, unique assignments or bounded model-transfer statement")
        }
        let selection = source.selection, budget = request.destination.budget
        try budget.validate()
        guard selection.mdCheckpoint.periodicCell == nil else {
            throw VivoChemistryError.unsupported("sampling reaction seed: periodic coordinates require an explicit qualified unwrapping/partition interface; no unwrapping or environment removal is inferred")
        }
        switch request.destination.calculation {
        case .qualify, .connectedReaction: break
        default: throw VivoChemistryError.unsupported("sampling reaction seed: only fresh nuclear or connected-reaction seed requests can receive coordinates")
        }
        // Admit destination coordinate validation before its pairwise checks.
        let seeds: [VivoNuclearQualificationRequest]
        switch request.destination.calculation {
        case .qualify(let q): seeds = [q]
        case .connectedReaction(let r): seeds = [r.saddle] + r.endpoints.flatMap { $0.components.map(\.qualification) }
        default: throw invalid("unsupported destination")
        }
        var validationWork = 0
        for seed in seeds {
            let n = seed.model.system.nuclei.count
            let count = try budget.elements([n, n], simultaneousArrays: 1)
            let sum = validationWork.addingReportingOverflow(count)
            guard !sum.overflow, sum.partialValue <= budget.maximumOperatorApplications else {
                throw VivoChemistryError.resourceLimit("sampling reaction seed destination validation work")
            }
            validationWork = sum.partialValue
        }
        try request.destination.validate()
        let mapped = try VivoPlatformSnapshotOperations.map(source: .init(structure: selection.request.structure),
            system: selection.request.system, checkpoint: selection.mdCheckpoint, budget: budget)
        var transfers: [VivoSamplingReactionSeedTransfer] = []
        func replace(_ initial: VivoNuclearQualificationRequest, assignment: VivoSamplingReactionSeedAssignment) throws -> VivoNuclearQualificationRequest {
            let indices = assignment.sourceAtomByNucleus, atoms = mapped.document.structure.atoms
            guard indices.count == initial.model.system.nuclei.count,
                  Set(indices).count == indices.count,
                  indices.allSatisfy({ $0 >= 0 && $0 < atoms.count }) else { throw invalid("mapping must cover every destination nucleus with distinct source atoms") }
            var result = initial
            for (nucleus, atom) in indices.enumerated() {
                guard Int(atoms[atom].element.atomicNumber) == initial.model.system.nuclei[nucleus].atomicNumber else { throw invalid("mapped element differs") }
                let nm = mapped.frame.positionsNM[atom]
                result.model.system.nuclei[nucleus].positionBohr = SIMD3(nm.x, nm.y, nm.z) / VivoAtomicUnits.bohrInNM
            }
            try result.validate()
            let particles = indices.map { mapped.mapping.atomToParticle[$0] }
            transfers.append(.init(assignment: assignment, sourceParticleByNucleus: particles,
                sourceMassesDa: particles.map { selection.request.system.particles[Int($0)].massDa },
                destinationMassesDa: result.massesDa, destinationThermochemistry: result.thermochemistry,
                destinationSolver: result.model.solver,
                destinationBasisFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result.model.basis)),
                initialQualificationFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(initial)),
                assembledQualificationFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(result))))
            return result
        }
        let assembled: VivoReactionCalculationRequest
        switch request.destination.calculation {
        case .qualify(var q):
            guard request.assignments.count == 1, request.assignments[0].target == .qualification else { throw invalid("standalone qualification target differs") }
            q = try replace(q, assignment: request.assignments[0]); assembled = .init(.qualify(request: q))
        case .connectedReaction(var r):
            for assignment in request.assignments {
                switch assignment.target {
                case .qualification: throw invalid("standalone target is not in a connected reaction")
                case .connectedSaddle: r.saddle = try replace(r.saddle, assignment: assignment)
                case .connectedEndpoint(let identifier, let component):
                    guard let endpoint = r.endpoints.firstIndex(where: { $0.identifier == identifier }),
                          component >= 0, component < r.endpoints[endpoint].components.count else { throw invalid("endpoint identifier or component index is absent") }
                    r.endpoints[endpoint].components[component].qualification = try replace(r.endpoints[endpoint].components[component].qualification, assignment: assignment)
                }
            }
            assembled = .init(.connectedReaction(request: r))
        default: throw invalid("unsupported destination")
        }
        try assembled.validate(); try Task.checkCancellation()
        return (assembled, transfers)
    }

    private struct Prepared: Sendable {
        let request: VivoSamplingReactionSeedRequest
        let task: VivoChemistryTask
        let provenance: VivoSamplingReactionSeedProvenance
        let inputs: [String: Data]
        let payloads: [String: Data]
    }

    private static func prepare(_ request: VivoSamplingReactionSeedRequest, source: VivoVerifiedMolecularSamplingExport,
                                implementationFingerprint: VivoFingerprint) throws -> Prepared {
        let assembled = try assembleAndDescribe(request, source: source), budget = request.destination.budget
        let requestData = try VivoCanonicalJSON.encode(request)
        let sourceData = try VivoCanonicalJSON.encode(source.receipt)
        let requestID = try VivoCanonicalJSON.fingerprint(requestData)
        let sourceID = try VivoCanonicalJSON.fingerprint(sourceData)
        let assembledData = try VivoCanonicalJSON.encode(assembled.request)
        let provenance = VivoSamplingReactionSeedProvenance(schema: "numivivo.org/sampling-reaction-seed-provenance/v1",
            assemblyRequestFingerprint: requestID, sourceExportReceiptFingerprint: sourceID,
            source: source.receipt.provenance, sourceConfiguration: source.selection.configuration,
            initialReactionRequestFingerprint: try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(request.destination)),
            assembledReactionRequestFingerprint: try VivoCanonicalJSON.fingerprint(assembledData),
            transfers: assembled.transfers, modelTransferStatement: request.modelTransferStatement, interpretation: interpretation)
        let inputs = ["assembly": requestData, "sampling-export": sourceData]
        let payloads = ["request": assembledData, "provenance": try VivoCanonicalJSON.encode(provenance)]
        for group in [inputs, payloads] {
            var remaining = budget.maximumBytes
            for bytes in group.values {
                guard bytes.count <= remaining else { throw VivoChemistryError.resourceLimit("sampling reaction seed input/output byte budget") }
                remaining -= bytes.count
            }
        }
        let task = VivoChemistryTask(operation: operationIdentifier, version: operationVersion,
            implementationFingerprint: implementationFingerprint, inputs: [
                .init(name: "assembly", artifact: requestID, kind: "sampling-reaction-seed-request"),
                .init(name: "sampling-export", artifact: sourceID, kind: "molecular-sampling-export")
            ], configuration: .object([:]), outputs: outputs,
            resources: .init(budget: budget, maximumInputBytes: budget.maximumBytes, maximumOutputBytes: budget.maximumBytes))
        _ = try task.fingerprint()
        return .init(request: request, task: task, provenance: provenance, inputs: inputs, payloads: payloads)
    }

    private static func operation(_ prepared: Prepared) -> VivoChemistryOperation {
        let validate: @Sendable (VivoJSONValue, [String: Data], [String: Data], VivoChemistryBudget) throws -> Void = {
            configuration, inputs, payloads, budget in
            try Task.checkCancellation()
            guard configuration == prepared.task.configuration, budget == prepared.task.resources.budget,
                  inputs == prepared.inputs, payloads == prepared.payloads else { throw invalid("workflow source, coordinates or model context differs") }
        }
        return .init(identifier: operationIdentifier, version: operationVersion,
            implementationFingerprint: prepared.task.implementationFingerprint, outputs: outputs,
            execute: { configuration, inputs, budget in
                try validate(configuration, inputs, prepared.payloads, budget); return prepared.payloads
            }, validateOutputs: validate)
    }

    private static func payloadFingerprints(_ payloads: [String: Data]) throws -> [String: VivoFingerprint] {
        try payloads.mapValues { try VivoCanonicalJSON.fingerprint($0) }
    }
    private static func makeReceipt(_ prepared: Prepared, result: VivoChemistryTaskResult) throws -> VivoSamplingReactionSeedReceipt {
        .init(schema: VivoSamplingReactionSeedReceipt.schemaID, implementationFingerprint: prepared.task.implementationFingerprint,
            request: prepared.request, task: prepared.task, taskFingerprint: result.taskFingerprint,
            workflowReceiptFingerprint: result.receiptFingerprint, provenance: prepared.provenance,
            outputs: result.outputs, payloadFingerprints: try payloadFingerprints(prepared.payloads))
    }
    private static func invalid(_ message: String) -> VivoChemistryError { .invalid("sampling reaction seed: " + message) }
}
