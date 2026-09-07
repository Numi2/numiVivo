import Foundation

/// Read and analysis admission limits, independent of the producer's budgets.
/// Object limits are applied before file allocation. The aggregate read budget
/// includes conservative trajectory object ceilings: 64 KiB per index link and
/// the admitted wire payload size. Descriptors have the store's separate metadata
/// limit. This budget is not a measurement of bytes read or peak RSS.
public struct VivoMolecularSamplingReadLimits: Sendable {
    public let maximumRequestBytes: Int
    public let maximumCursorBytes: Int
    public let maximumCheckpointBytes: Int
    public let maximumDiagnosticBytes: Int
    public let maximumReadBytes: UInt64
    public let maximumReplicas: Int
    public let maximumParticles: Int
    public let maximumScalarElements: Int
    public let maximumAcceptedSteps: UInt64
    public let maximumDiagnosticEvaluations: Int
    public let maximumAutocorrelationWork: UInt64
    public let maximumTrajectoryChunks: UInt64?

    public init(maximumRequestBytes: Int = 128 * 1024 * 1024,
                maximumCursorBytes: Int = 128 * 1024 * 1024,
                maximumCheckpointBytes: Int = 128 * 1024 * 1024,
                maximumDiagnosticBytes: Int = 128 * 1024 * 1024,
                maximumReadBytes: UInt64 = 2 * 1024 * 1024 * 1024,
                maximumReplicas: Int = 64, maximumParticles: Int = 1_000_000,
                maximumScalarElements: Int = 2_000_000,
                maximumAcceptedSteps: UInt64 = 10_000_000_000,
                maximumDiagnosticEvaluations: Int = 101,
                maximumAutocorrelationWork: UInt64 = 1_000_000_000,
                maximumTrajectoryChunks: UInt64? = nil) {
        self.maximumRequestBytes = maximumRequestBytes; self.maximumCursorBytes = maximumCursorBytes
        self.maximumCheckpointBytes = maximumCheckpointBytes; self.maximumDiagnosticBytes = maximumDiagnosticBytes
        self.maximumReadBytes = maximumReadBytes; self.maximumReplicas = maximumReplicas
        self.maximumParticles = maximumParticles; self.maximumScalarElements = maximumScalarElements
        self.maximumAcceptedSteps = maximumAcceptedSteps; self.maximumDiagnosticEvaluations = maximumDiagnosticEvaluations
        self.maximumAutocorrelationWork = maximumAutocorrelationWork; self.maximumTrajectoryChunks = maximumTrajectoryChunks
    }

    fileprivate func validate() throws {
        guard maximumRequestBytes > 0, maximumCursorBytes > 0, maximumCheckpointBytes > 0,
              maximumDiagnosticBytes > 0, maximumReadBytes > 0, maximumReplicas > 0,
              maximumParticles > 0, maximumScalarElements > 0, maximumAcceptedSteps > 0,
              maximumDiagnosticEvaluations > 0, maximumAutocorrelationWork > 0 else {
            throw samplingArchiveInvalid("read and analysis limits must be positive")
        }
    }
}

/// Accepted-state evidence always includes the selected final coordinate chunk.
public enum VivoMolecularSamplingSelectionValidation: String, Sendable, Equatable {
    case restart, allPayloads
}

/// The cursor has no terminal status. Only an explicitly supplied, bound receipt
/// can distinguish a recorded cancellation/rejection from an active run/crash.
public enum VivoMolecularSamplingSourceTermination: String, Sendable, Equatable {
    case notRecorded, converged, budgetExhausted, rejected, cancelled
}

public struct VivoMolecularSamplingReplicaInspection: Sendable {
    public let index: Int
    public let identifier: String
    public let seed: UInt64
    public let mdCheckpointFingerprint: VivoFingerprint?
    public let trajectoryManifestFingerprint: VivoFingerprint?
}

public struct VivoMolecularSamplingInspection: Sendable {
    public let requestFingerprint: VivoFingerprint
    public let samplingCheckpointFingerprint: VivoFingerprint
    public let completedBlocks: Int
    public let consecutivePasses: Int
    public let requiredConsecutivePasses: Int
    public let declaredObservableCriteriaSatisfied: Bool
    public let latestDiagnosticFingerprint: VivoFingerprint?
    public let latestDiagnostic: VivoMolecularSamplingResult?
    public let replicas: [VivoMolecularSamplingReplicaInspection]
    public let diagnosticEvaluations: Int
    /// The exact current diagnostic and final pass suffix are reconstructed from
    /// bound scalar prefixes. This does not establish absence of a much earlier
    /// stopping opportunity, nor remeasure observables from coordinate payloads.
    public let criteriaHistoryScope: String
}

/// Immutable, non-decodable accepted-prefix authority. The original checkpoint
/// bytes and hash are retained; no reconstruction can restamp its MD contract.
/// The captured actor pins the rooted store even when its path is moved/reused.
/// Evidence applies at validation time, not after arbitrary filesystem changes.
public struct VivoValidatedMolecularSamplingSelection: Sendable {
    let store: VivoArtifactStore
    public let requestFingerprint: VivoFingerprint
    public let samplingCheckpointFingerprint: VivoFingerprint
    public let request: VivoMolecularSamplingRunRequest
    public let replicaIndex: Int
    public let replicaIdentifier: String
    public let replicaSeed: UInt64
    public let configuration: VivoMDConfiguration
    public let mdCheckpointFingerprint: VivoFingerprint
    public let mdCheckpointData: Data
    public let mdCheckpoint: VivoMDCheckpoint
    public let trajectoryManifestFingerprint: VivoFingerprint
    public let trajectoryManifest: VivoMDTrajectoryManifest
    public let trajectoryValidation: VivoMDTrajectoryValidation
    /// Retained producer reference only; minimization history is not revalidated.
    public let minimizationFingerprint: VivoFingerprint?
    public let inspection: VivoMolecularSamplingInspection
    public let sourceReceiptFingerprint: VivoFingerprint?
    public let sourceReceipt: VivoMolecularSamplingRunReceipt?
    public let sourceTermination: VivoMolecularSamplingSourceTermination

    fileprivate init(reader: VivoMolecularSamplingArchiveReader, index: Int,
                     checkpoint: VivoMDCheckpoint, checkpointData: Data,
                     archive: VivoMDTrajectoryArchiveReader, validation: VivoMDTrajectoryValidation,
                     receiptFingerprint: VivoFingerprint?, receipt: VivoMolecularSamplingRunReceipt?) {
        store = reader.store; requestFingerprint = reader.cursor.requestFingerprint
        samplingCheckpointFingerprint = reader.checkpointFingerprint; request = reader.request
        replicaIndex = index; replicaIdentifier = reader.cursor.replicas[index].series.identifier
        replicaSeed = reader.request.replicaSeeds[index]; configuration = reader.cursor.replicas[index].series.configuration
        mdCheckpointFingerprint = reader.cursor.replicas[index].mdCheckpoint!
        mdCheckpointData = checkpointData; mdCheckpoint = checkpoint
        trajectoryManifestFingerprint = archive.manifestFingerprint; trajectoryManifest = archive.manifest
        trajectoryValidation = validation; minimizationFingerprint = reader.cursor.replicas[index].minimization
        inspection = reader.inspect()
        sourceReceiptFingerprint = receiptFingerprint; sourceReceipt = receipt
        sourceTermination = receipt.map { VivoMolecularSamplingSourceTermination(rawValue: $0.status.rawValue)! } ?? .notRecorded
    }
}

/// Internal continuation evidence binds the captured exact checkpoint to its
/// verified archive. Runner handoff does not reread/restamp the checkpoint.
struct VivoMolecularSamplingReplicaRestart: Sendable {
    let checkpointFingerprint: VivoFingerprint
    let checkpoint: VivoMDCheckpoint
    let trajectoryValidation: VivoMDTrajectoryValidation
    fileprivate init(checkpointFingerprint: VivoFingerprint, checkpoint: VivoMDCheckpoint,
                     trajectoryValidation: VivoMDTrajectoryValidation) {
        self.checkpointFingerprint = checkpointFingerprint; self.checkpoint = checkpoint
        self.trajectoryValidation = trajectoryValidation
    }
}

/// Read-only admission shared by replica selection and production continuation.
/// open validates every replica's checkpoint/manifest/scalar metadata and the
/// declared-observable suffix. Selection additionally verifies the chosen index
/// chain and payload scope; runner continuation requires this for every replica.
public struct VivoMolecularSamplingArchiveReader: Sendable {
    let store: VivoArtifactStore
    public let checkpointFingerprint: VivoFingerprint
    fileprivate let request: VivoMolecularSamplingRunRequest
    let cursor: VivoMolecularSamplingCursor
    let latest: VivoMolecularSamplingResult?
    private let limits: VivoMolecularSamplingReadLimits
    private let metadataBudget: VivoMolecularSamplingReadBudget
    private let replicas: [Replica]
    private let diagnosticEvaluations: Int

    private struct Replica: Sendable {
        let checkpoint: VivoMDCheckpoint?
        let checkpointData: Data?
        let archive: VivoMDTrajectoryArchiveReader?
    }

    private init(store: VivoArtifactStore, fingerprint: VivoFingerprint,
                 request: VivoMolecularSamplingRunRequest, cursor: VivoMolecularSamplingCursor,
                 latest: VivoMolecularSamplingResult?, limits: VivoMolecularSamplingReadLimits,
                 budget: VivoMolecularSamplingReadBudget, replicas: [Replica], diagnosticEvaluations: Int) {
        self.store = store; checkpointFingerprint = fingerprint; self.request = request; self.cursor = cursor
        self.latest = latest; self.limits = limits; metadataBudget = budget
        self.replicas = replicas; self.diagnosticEvaluations = diagnosticEvaluations
    }

    public static func open(store: VivoArtifactStore, checkpoint fingerprint: VivoFingerprint,
                            limits: VivoMolecularSamplingReadLimits = .init()) async throws -> Self {
        try Task.checkCancellation(); try limits.validate()
        var budget = VivoMolecularSamplingReadBudget(maximum: limits.maximumReadBytes)
        let cursorData = try await budget.read(fingerprint, kind: "molecular-sampling-checkpoint",
                                              maximumBytes: limits.maximumCursorBytes, store: store)
        let cursor = try VivoCanonicalJSON.decode(VivoMolecularSamplingCursor.self, from: cursorData)
        guard cursor.schema == VivoMolecularSamplingCursor.schema else { throw samplingArchiveInvalid("cursor schema") }
        let requestData = try await budget.read(cursor.requestFingerprint, kind: "molecular-sampling-run",
                                               maximumBytes: limits.maximumRequestBytes, store: store)
        let request = try VivoCanonicalJSON.decode(VivoMolecularSamplingRunRequest.self, from: requestData)
        guard request.replicaSeeds.count <= limits.maximumReplicas, request.system.particles.count <= limits.maximumParticles,
              request.structure.atoms.count <= limits.maximumParticles,
              request.system.particles.count <= UInt32.max, request.structure.atoms.count <= UInt32.max else {
            throw samplingArchiveInvalid("replica/particle read budget")
        }
        try request.validate(); try Task.checkCancellation()
        guard cursor.completedBlocks >= 0, cursor.completedBlocks <= request.maximumBlocks,
              cursor.consecutivePasses >= 0, cursor.consecutivePasses <= cursor.completedBlocks,
              cursor.consecutivePasses <= request.requiredConsecutivePasses,
              cursor.replicas.count == request.replicaSeeds.count, cursor.diagnostics.count == cursor.completedBlocks else {
            throw samplingArchiveInvalid("block cursor or replica count")
        }
        let productionSteps = try samplingArchiveProduct(UInt64(cursor.completedBlocks), request.stepsPerBlock)
        let acceptedStep = try samplingArchiveSum(request.equilibrationSteps, productionSteps)
        guard cursor.completedBlocks == 0 || acceptedStep <= limits.maximumAcceptedSteps else {
            throw samplingArchiveInvalid("accepted-step validation budget")
        }
        let frames = productionSteps / request.sampleEvery
        let scalarCount = try samplingArchiveProduct(try samplingArchiveProduct(frames, UInt64(cursor.replicas.count)), UInt64(request.observables.count))
        guard scalarCount <= UInt64(limits.maximumScalarElements) else { throw samplingArchiveInvalid("scalar-element read budget") }
        let systemID = try request.system.fingerprint(), structureID = try VivoStructureCodec.fingerprint(request.structure)
        var replicas: [Replica] = []
        for (index, replica) in cursor.replicas.enumerated() {
            try Task.checkCancellation()
            var configuration = request.md; configuration.randomSeed = request.replicaSeeds[index]
            let series = replica.series, initial = request.initialStates[index]
            guard (initial.sourceTimePS ?? 0) >= 0,
                  series.identifier == "replica-\(index)", series.configuration == configuration,
                  UInt64(series.timesPS.count) == frames, series.steps.count == series.timesPS.count,
                  series.valuesByObservable.count == request.observables.count,
                  series.valuesByObservable.allSatisfy({ $0.count == series.timesPS.count && $0.allSatisfy(\.isFinite) }) else {
                throw samplingArchiveInvalid("replica identity, configuration or scalar prefix shape")
            }
            if cursor.completedBlocks == 0 {
                guard replica.mdCheckpoint == nil, replica.trajectory == nil, replica.minimization == nil,
                      series.sourceFingerprint == (try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(initial))) else {
                    throw samplingArchiveInvalid("block zero cannot own unpublished replica state")
                }
                replicas.append(.init(checkpoint: nil, checkpointData: nil, archive: nil)); continue
            }
            for sample in series.steps.indices {
                if sample % 1024 == 0 { try Task.checkCancellation() }
                let step = try samplingArchiveSum(request.equilibrationSteps,
                    samplingArchiveProduct(UInt64(sample) + 1, request.sampleEvery))
                guard series.steps[sample] == step else { throw samplingArchiveInvalid("scalar accepted-step schedule") }
                try validateClock(series.timesPS[sample], initialTime: initial.sourceTimePS ?? 0,
                                  steps: step, timeStep: configuration.timeStepPS)
            }
            guard let checkpointID = replica.mdCheckpoint, let trajectoryID = replica.trajectory,
                  series.sourceFingerprint == trajectoryID,
                  (replica.minimization != nil) == (request.minimization != nil) else {
                throw samplingArchiveInvalid("missing durable replica state or minimization reference")
            }
            let data = try await budget.read(checkpointID, kind: "md-checkpoint",
                                            maximumBytes: limits.maximumCheckpointBytes, store: store)
            let checkpoint = try VivoCanonicalJSON.decode(VivoMDCheckpoint.self, from: data)
            try checkpoint.validate(particleCount: request.system.particles.count)
            _ = try await budget.read(trajectoryID, kind: "md-trajectory-manifest", maximumBytes: 64 * 1024, store: store)
            try budget.reserve(64 * 1024) // The trajectory owner freshly hashes its own bounded manifest read.
            let archive = try await VivoMDTrajectoryArchiveReader.open(store: store, manifest: trajectoryID)
            let manifest = archive.manifest
            guard checkpoint.systemFingerprint == systemID,
                  checkpoint.configurationFingerprint == (try configuration.fingerprint()), checkpoint.acceptedStep == acceptedStep,
                  manifest.systemFingerprint == systemID, manifest.configurationFingerprint == checkpoint.configurationFingerprint,
                  manifest.particleCount == UInt32(request.system.particles.count), !manifest.includeVelocities,
                  manifest.frameCount == frames, manifest.firstStep == series.steps.first, manifest.firstTimePS == series.timesPS.first,
                  manifest.lastStep == checkpoint.acceptedStep, manifest.lastTimePS == checkpoint.timePS, !manifest.sealed,
                  series.steps.last == checkpoint.acceptedStep, series.timesPS.last == checkpoint.timePS,
                  configuration.ensemble == .npt || checkpoint.periodicCell == initial.periodicCell else {
                throw samplingArchiveInvalid("manifest, scalar prefix and exact checkpoint binding")
            }
            try validateClock(checkpoint.timePS, initialTime: initial.sourceTimePS ?? 0,
                              steps: acceptedStep, timeStep: configuration.timeStepPS)
            replicas.append(.init(checkpoint: checkpoint, checkpointData: data, archive: archive))
        }
        var latest: VivoMolecularSamplingResult?, evaluations = 0, admittedWork: UInt64 = 0
        if let latestID = cursor.diagnostics.last {
            let data = try await budget.read(latestID, kind: "molecular-sampling-result",
                                            maximumBytes: limits.maximumDiagnosticBytes, store: store)
            let recorded = try VivoCanonicalJSON.decode(VivoMolecularSamplingResult.self, from: data)
            // Historical manifest hashes legitimately change between blocks. The
            // numerical prefix analysis uses the currently bound source IDs; old
            // diagnostic pass flags are never authority for the current counter.
            func analyze(blocks: Int) throws -> VivoMolecularSamplingResult {
                try Task.checkCancellation()
                guard evaluations < limits.maximumDiagnosticEvaluations else { throw samplingArchiveInvalid("diagnostic evaluation budget") }
                let work = try samplingArchiveSum(admittedWork, UInt64(request.convergence.maximumAutocorrelationWork))
                guard work <= limits.maximumAutocorrelationWork else { throw samplingArchiveInvalid("diagnostic autocorrelation work budget") }
                admittedWork = work; evaluations += 1
                let count = Int(UInt64(blocks) * request.stepsPerBlock / request.sampleEvery)
                let series = cursor.replicas.map { replica -> VivoMolecularReplicaSeries in
                    var value = replica.series
                    value.steps = Array(value.steps.prefix(count)); value.timesPS = Array(value.timesPS.prefix(count))
                    value.valuesByObservable = value.valuesByObservable.map { Array($0.prefix(count)) }
                    return value
                }
                let analysis = VivoMolecularSamplingRequest(structureFingerprint: structureID, systemFingerprint: systemID,
                    contextIdentifier: request.contextIdentifier, observables: request.observables,
                    replicas: series, configuration: request.convergence)
                let result = try VivoMolecularSampling.analyze(analysis)
                try Task.checkCancellation(); return result
            }
            let current = try analyze(blocks: cursor.completedBlocks)
            guard recorded == current else { throw samplingArchiveInvalid("current diagnostic does not reconstruct from its bound scalar prefix") }
            latest = current
            var passes = current.converged ? 1 : 0
            if current.converged && cursor.completedBlocks > 1 {
                for blocks in stride(from: cursor.completedBlocks - 1, through: 1, by: -1) {
                    if !(try analyze(blocks: blocks)).converged { break }
                    passes += 1
                    guard passes <= request.requiredConsecutivePasses else { throw samplingArchiveInvalid("pass suffix exceeds the stopping threshold") }
                }
            }
            guard passes == cursor.consecutivePasses else { throw samplingArchiveInvalid("consecutive-pass suffix differs from the scalar prefixes") }
        }
        try Task.checkCancellation()
        return .init(store: store, fingerprint: fingerprint, request: request, cursor: cursor,
                     latest: latest, limits: limits, budget: budget, replicas: replicas, diagnosticEvaluations: evaluations)
    }

    public func inspect() -> VivoMolecularSamplingInspection {
        .init(requestFingerprint: cursor.requestFingerprint, samplingCheckpointFingerprint: checkpointFingerprint,
              completedBlocks: cursor.completedBlocks, consecutivePasses: cursor.consecutivePasses,
              requiredConsecutivePasses: request.requiredConsecutivePasses,
              declaredObservableCriteriaSatisfied: cursor.consecutivePasses >= request.requiredConsecutivePasses,
              latestDiagnosticFingerprint: cursor.diagnostics.last, latestDiagnostic: latest,
              replicas: cursor.replicas.enumerated().map { index, value in
                  .init(index: index, identifier: value.series.identifier, seed: request.replicaSeeds[index],
                        mdCheckpointFingerprint: value.mdCheckpoint, trajectoryManifestFingerprint: value.trajectory)
              }, diagnosticEvaluations: diagnosticEvaluations,
              criteriaHistoryScope: "current diagnostic and final consecutive-pass suffix; earlier stopping opportunities and scalar remeasurement are not verified")
    }

    public func selectReplica(index: Int, validation: VivoMolecularSamplingSelectionValidation = .restart,
                              receipt: VivoFingerprint? = nil) async throws -> VivoValidatedMolecularSamplingSelection {
        try Task.checkCancellation()
        var budget = metadataBudget
        let selected = try await validateReplica(index: index, scope: validation, budget: &budget)
        let sourceReceipt = try await readReceipt(receipt, budget: &budget)
        try Task.checkCancellation()
        return .init(reader: self, index: index, checkpoint: selected.checkpoint, checkpointData: selected.data,
                     archive: selected.archive, validation: selected.validation,
                     receiptFingerprint: receipt, receipt: sourceReceipt)
    }

    /// Immediate internal handoff; never accepts a substitute store or manifest.
    func validateForResume(request expected: VivoMolecularSamplingRunRequest) async throws -> [VivoMolecularSamplingReplicaRestart] {
        guard request == expected else { throw samplingArchiveInvalid("resume request differs from the captured request") }
        try Task.checkCancellation()
        if cursor.completedBlocks == 0 { return [] }
        var budget = metadataBudget, tokens: [VivoMolecularSamplingReplicaRestart] = []
        for index in replicas.indices {
            let replica = try await validateReplica(index: index, scope: .restart, budget: &budget)
            tokens.append(.init(checkpointFingerprint: cursor.replicas[index].mdCheckpoint!,
                                checkpoint: replica.checkpoint, trajectoryValidation: replica.validation))
        }
        try Task.checkCancellation(); return tokens
    }

    private func validateReplica(index: Int, scope: VivoMolecularSamplingSelectionValidation,
                                 budget: inout VivoMolecularSamplingReadBudget) async throws
        -> (checkpoint: VivoMDCheckpoint, data: Data, archive: VivoMDTrajectoryArchiveReader, validation: VivoMDTrajectoryValidation) {
        try Task.checkCancellation()
        guard replicas.indices.contains(index) else { throw samplingArchiveInvalid("selected replica index is out of range") }
        guard cursor.completedBlocks > 0, let checkpoint = replicas[index].checkpoint,
              let data = replicas[index].checkpointData, let archive = replicas[index].archive else {
            throw samplingArchiveInvalid("no accepted complete production block is available")
        }
        let initial = VivoClassicalInitialState(systemFingerprint: checkpoint.systemFingerprint,
            positionsNM: checkpoint.positionsNM,
            periodicCell: checkpoint.periodicCell, sourceTimePS: checkpoint.timePS)
        let capabilities = try VivoMDCapabilityAnalyzer.analyze(system: request.system, initialState: initial,
                                                               configuration: cursor.replicas[index].series.configuration)
        guard capabilities.executable else { throw VivoMDRuntimeError.unsupported(capabilities.blockers) }
        let manifest = archive.manifest
        if let limit = limits.maximumTrajectoryChunks, manifest.chunkCount > limit { throw samplingArchiveInvalid("trajectory chunk work budget") }
        let stride = try VivoMDTrajectoryChunkCodec.frameBytes(particleCount: manifest.particleCount, includeVelocities: false)
        let wireBytes = try samplingArchiveSum(samplingArchiveProduct(manifest.frameCount, UInt64(stride)),
            samplingArchiveProduct(manifest.chunkCount, UInt64(VivoMDTrajectoryChunkCodec.headerBytes)))
        let tailCeiling = min(wireBytes, UInt64(VivoMDTrajectoryChunkCodec.maximumChunkBytes))
        try budget.reserve(samplingArchiveProduct(manifest.chunkCount, 64 * 1024))
        try budget.reserve(scope == .allPayloads ? wireBytes : tailCeiling)
        // One fresh bounded tail read compares exact checkpoint coordinates. It
        // does not repeat the index walk, and its bytes are separately admitted.
        try budget.reserve(samplingArchiveSum(64 * 1024, tailCeiling))
        let verified = try await archive.validate(scope: scope == .allPayloads ? .allPayloads : .restart,
                                                  maximumChunks: limits.maximumTrajectoryChunks)
        guard let tail = manifest.tail, let frame = try await archive.readChunk(tail).last,
              frame.stepIndex == checkpoint.acceptedStep, frame.timePS == checkpoint.timePS,
              frame.positionsNM == checkpoint.positionsNM, frame.periodicCell == checkpoint.periodicCell else {
            throw samplingArchiveInvalid("selected final payload differs from the exact checkpoint")
        }
        try Task.checkCancellation(); return (checkpoint, data, archive, verified)
    }

    private func readReceipt(_ fingerprint: VivoFingerprint?, budget: inout VivoMolecularSamplingReadBudget) async throws -> VivoMolecularSamplingRunReceipt? {
        guard let fingerprint else { return nil }
        let data = try await budget.read(fingerprint, kind: "molecular-sampling-receipt",
                                        maximumBytes: limits.maximumDiagnosticBytes, store: store)
        let receipt = try VivoCanonicalJSON.decode(VivoMolecularSamplingRunReceipt.self, from: data)
        guard receipt.schema == VivoMolecularSamplingRunReceipt.schema, receipt.requestFingerprint == cursor.requestFingerprint,
              receipt.checkpoint == checkpointFingerprint, receipt.completedBlocks == cursor.completedBlocks,
              receipt.consecutivePasses == cursor.consecutivePasses, receipt.trajectoryManifests == cursor.replicas.compactMap(\.trajectory),
              receipt.sampling == latest, receipt.interpretation == VivoMolecularSampling.interpretation else {
            throw samplingArchiveInvalid("source receipt differs from the captured complete prefix")
        }
        switch receipt.status {
        case .converged:
            guard cursor.consecutivePasses >= request.requiredConsecutivePasses else { throw samplingArchiveInvalid("receipt falsely reports converged") }
        case .budgetExhausted:
            guard cursor.completedBlocks == request.maximumBlocks, cursor.consecutivePasses < request.requiredConsecutivePasses else {
                throw samplingArchiveInvalid("receipt falsely reports budget exhaustion")
            }
        case .cancelled, .rejected: break
        }
        return receipt
    }

    private static func validateClock(_ time: Double, initialTime: Double, steps: UInt64, timeStep: Double) throws {
        guard time.isFinite, time >= initialTime else { throw samplingArchiveInvalid("nonfinite or regressing scalar clock") }
        if steps == 0 {
            guard time == initialTime else { throw samplingArchiveInvalid("zero-step clock differs") }; return
        }
        // Same conservative repeated-FP64-addition bound as MD protocol resume.
        // nextUp covers conversion of any caller-admitted UInt64 step count.
        let count = Double(steps), upperCount = count.nextUp
        let expected = initialTime + count * timeStep
        let roundoff = (upperCount + 2) * (Double.ulpOfOne / 2)
        guard roundoff < 1 else { throw samplingArchiveInvalid("clock bound cannot admit this step count") }
        let tolerance = 2 * (roundoff / (1 - roundoff) * max(expected, Double.leastNormalMagnitude))
        guard time > initialTime, expected.isFinite, tolerance.isFinite, abs(time - expected) <= tolerance else {
            throw samplingArchiveInvalid("clock differs from its original source time and accepted-step schedule")
        }
    }
}

private struct VivoMolecularSamplingReadBudget: Sendable {
    let maximum: UInt64
    private var admitted: UInt64 = 0
    init(maximum: UInt64) { self.maximum = maximum }
    mutating func reserve(_ amount: UInt64) throws {
        let total = try samplingArchiveSum(admitted, amount)
        guard total <= maximum else { throw samplingArchiveInvalid("aggregate read budget") }
        admitted = total
    }
    mutating func read(_ fingerprint: VivoFingerprint, kind: String, maximumBytes: Int,
                       store: VivoArtifactStore) async throws -> Data {
        try Task.checkCancellation()
        let descriptor = try await store.descriptor(for: fingerprint)
        try Task.checkCancellation()
        guard descriptor.kind == kind, descriptor.byteCount <= UInt64(maximumBytes) else {
            throw samplingArchiveInvalid("artifact kind or object read budget")
        }
        try reserve(descriptor.byteCount)
        let data = try await store.data(for: fingerprint, maximumBytes: min(maximumBytes, Int(descriptor.byteCount)), verify: true)
        try Task.checkCancellation()
        guard UInt64(data.count) == descriptor.byteCount else { throw samplingArchiveInvalid("artifact descriptor byte count") }
        return data
    }
}

private func samplingArchiveInvalid(_ message: String) -> VivoArtifactValidationError {
    .invalid("molecular sampling archive: " + message)
}
private func samplingArchiveProduct(_ a: UInt64, _ b: UInt64) throws -> UInt64 {
    let result = a.multipliedReportingOverflow(by: b)
    guard !result.overflow else { throw samplingArchiveInvalid("count multiplication overflow") }; return result.partialValue
}
private func samplingArchiveSum(_ a: UInt64, _ b: UInt64) throws -> UInt64 {
    let result = a.addingReportingOverflow(b)
    guard !result.overflow else { throw samplingArchiveInvalid("count addition overflow") }; return result.partialValue
}
