import Foundation

public enum VivoProteinStressAction: String, Codable, Sendable { case source, stageStart, sample, rejected }
public enum VivoProteinStressDisposition: String, Codable, Sendable { case completed, rejected, cancelled, failed }

public struct VivoProteinStressJournalEntry: Codable, Sendable, Equatable {
    public var schema: String = "numivivo.org/protein-stress-journal/v1"
    public let request: VivoFingerprint
    public let previous: VivoFingerprint?
    public let action: VivoProteinStressAction
    public let stageIndex: Int
    public let stepsInStage: UInt64
    public let checkpoint: VivoFingerprint
    public let cumulativeProtocolWorkKJPerMol: Double
    public let frame: VivoProteinStressFrame
    public let rejection: VivoMDStepCertificate?
}

public struct VivoProteinStressReceipt: Codable, Sendable, Equatable {
    public var schema: String = "numivivo.org/protein-stress-receipt/v1"
    public let request: VivoFingerprint
    public let disposition: VivoProteinStressDisposition
    public let journalTail: VivoFingerprint?
    public let checkpointReferenceName: String
    public let deviceName: String?
    public let diagnostic: String?
    public let interpretation: String
}

/// Domain orchestration only: every numerical step, force kernel, constraint,
/// thermostat, RNG counter and checkpoint remains owned by VivoMDMetalRuntime.
/// A bounded scalar/checkpoint journal uses the existing immutable artifact store.
public enum VivoProteinStressRunner {
    public static let interpretation = "piecewise-constant harmonic/NVT protocol; sampled forces and structural retention are conditional simulation observables, not unfolding-force qualification, melting temperature, equilibrium free energy, or bulk material properties"

    public static func run(_ request: VivoProteinStressRequest, store: VivoArtifactStore,
                           resumeFrom: VivoFingerprint? = nil) async throws -> VivoProteinStressReceipt {
        let compiled = try VivoProteinStressCompilation(request), requestID = try request.fingerprint()
        _ = try await put(request, kind: "protein-stress-request", store: store)
        let reference = "protein-stress-" + UUID().uuidString.lowercased()
        var current: VivoProteinStressJournalEntry?, tail: VivoFingerprint?, state = request.sourceCheckpoint
        var deviceName: String?, runtime: VivoMDMetalRuntime?
        func receipt(_ disposition: VivoProteinStressDisposition, _ diagnostic: String? = nil) -> VivoProteinStressReceipt {
            .init(request: requestID, disposition: disposition, journalTail: tail, checkpointReferenceName: reference,
                  deviceName: deviceName, diagnostic: diagnostic, interpretation: interpretation)
        }
        func entry(action: VivoProteinStressAction, stage: Int, steps: UInt64, state: VivoMDCheckpoint,
                   work: Double, rejection: VivoMDStepCertificate? = nil) async throws -> VivoProteinStressJournalEntry {
            let checkpoint = try await put(state, kind: "md-checkpoint", store: store)
            return .init(request: requestID, previous: tail, action: action, stageIndex: stage, stepsInStage: steps,
                         checkpoint: checkpoint.fingerprint, cumulativeProtocolWorkKJPerMol: work,
                         frame: try compiled.frame(checkpoint: state, stage: stage >= 0 ? stage : nil), rejection: rejection)
        }
        func publish(_ next: VivoProteinStressJournalEntry) async throws {
            let artifact = try await put(next, kind: "protein-stress-journal", store: store)
            // The immutable artifact already exists even if reference publication
            // fails; return that exact durable tail, never an unpersisted state.
            current = next; tail = artifact.fingerprint
            _ = try await store.setReference(reference, to: artifact)
        }
        do {
            if let resumeFrom {
                let verified = try await verify(request, store: store, journalTail: resumeFrom)
                current = verified.last!.entry; tail = resumeFrom
                state = try await read(VivoMDCheckpoint.self, id: current!.checkpoint, kind: "md-checkpoint", store: store)
                _ = try await store.setReference(reference, to: store.descriptor(for: resumeFrom))
                if current!.action == .rejected { return receipt(.rejected, "rejected protocol cannot resume as a successful candidate") }
            } else {
                try await publish(entry(action: .source, stage: -1, steps: 0, state: state, work: 0))
            }
            while let cursor = current {
                try Task.checkCancellation()
                let mustAdvance = cursor.stageIndex == -1 || cursor.stepsInStage == request.stages[cursor.stageIndex].steps
                if mustAdvance {
                    let destination = cursor.stageIndex + 1
                    if destination == request.stages.count { return receipt(.completed) }
                    let transfer = try compiled.transfer(state, from: cursor.stageIndex >= 0 ? cursor.stageIndex : nil, to: destination)
                    let cumulative = cursor.cumulativeProtocolWorkKJPerMol + transfer.workKJPerMol
                    guard cumulative.isFinite else { throw VivoProteinStressError.invalid("cumulative protocol work overflow") }
                    runtime = nil
                    let next = try await VivoMDMetalRuntime.restore(system: request.system, configuration: request.configuration(for: destination),
                        checkpoint: transfer.checkpoint, forceProvider: compiled.forceProvider(stage: destination))
                    deviceName = next.deviceName; runtime = next; state = transfer.checkpoint
                    try await publish(entry(action: .stageStart, stage: destination, steps: 0, state: state, work: cumulative))
                    continue
                }
                if runtime == nil {
                    runtime = try await VivoMDMetalRuntime.restore(system: request.system, configuration: request.configuration(for: cursor.stageIndex),
                        checkpoint: state, forceProvider: compiled.forceProvider(stage: cursor.stageIndex))
                    deviceName = runtime!.deviceName
                }
                guard let runtime else { throw VivoProteinStressError.invalid("running protocol has no MD owner") }
                var steps = cursor.stepsInStage
                let stage = request.stages[cursor.stageIndex]
                // Complete the next globally scheduled sample interval. This also
                // keeps a resumed run on the original sample grid.
                let untilSample = request.sampleEvery - steps % request.sampleEvery
                let stop = steps + min(stage.steps - steps, untilSample)
                while steps < stop {
                    try Task.checkCancellation()
                    let certificate = try await runtime.step()
                    if !certificate.committed {
                        state = try await runtime.checkpoint()
                        try await publish(entry(action: .rejected, stage: cursor.stageIndex, steps: steps, state: state,
                            work: cursor.cumulativeProtocolWorkKJPerMol, rejection: certificate))
                        return receipt(.rejected, "MD candidate rejected; no retry, timestep adjustment, or replacement replica")
                    }
                    steps += 1
                }
                state = try await runtime.checkpoint()
                try await publish(entry(action: .sample, stage: cursor.stageIndex, steps: steps, state: state,
                    work: cursor.cumulativeProtocolWorkKJPerMol))
            }
            throw VivoProteinStressError.invalid("missing durable source boundary")
        } catch {
            // An unsaved accepted suffix is not promoted. The durable journal is
            // the only restart authority; cancellation is distinct from rejection.
            return receipt(error is CancellationError ? .cancelled : .failed, String(describing: error))
        }
    }

    public struct VerifiedEntry: Sendable {
        public let fingerprint: VivoFingerprint
        public let entry: VivoProteinStressJournalEntry
        public let acceptedStep: UInt64
        public let timePS: Double
    }

    /// Reconstructs every stored metric and declared transition, checks original
    /// request bytes and exact checkpoint identities/clocks, and rejects truncated
    /// or reordered chains. Integrity is not a proof of MD physical accuracy.
    public static func verify(_ request: VivoProteinStressRequest, store: VivoArtifactStore,
                              journalTail: VivoFingerprint) async throws -> [VerifiedEntry] {
        let compiled = try VivoProteinStressCompilation(request), requestID = try request.fingerprint()
        let storedRequest: VivoProteinStressRequest = try await read(VivoProteinStressRequest.self, id: requestID,
                                                                      kind: "protein-stress-request", store: store)
        guard storedRequest == request else { throw VivoProteinStressError.invalid("request bytes differ") }
        var reversed: [(VivoFingerprint, VivoProteinStressJournalEntry)] = [], next: VivoFingerprint? = journalTail
        var seen = Set<VivoFingerprint>()
        while let id = next {
            guard reversed.count < compiled.maximumJournalEntries, seen.insert(id).inserted else {
                throw VivoProteinStressError.invalid("journal traversal limit or cycle")
            }
            let value: VivoProteinStressJournalEntry = try await read(VivoProteinStressJournalEntry.self, id: id, kind: "protein-stress-journal", store: store)
            guard value.schema == "numivivo.org/protein-stress-journal/v1", value.request == requestID,
                  value.cumulativeProtocolWorkKJPerMol.isFinite else { throw VivoProteinStressError.invalid("journal identity or work") }
            reversed.append((id, value)); next = value.previous
        }
        var prior: VivoProteinStressJournalEntry?, priorState: VivoMDCheckpoint?, result: [VerifiedEntry] = []
        for (id, value) in reversed.reversed() {
            try Task.checkCancellation()
            let state: VivoMDCheckpoint = try await read(VivoMDCheckpoint.self, id: value.checkpoint, kind: "md-checkpoint", store: store)
            try state.validate(particleCount: compiled.massesDa.count)
            guard state.systemFingerprint == request.sourceCheckpoint.systemFingerprint,
                  state.periodicCell == request.sourceCheckpoint.periodicCell,
                  (state.positionPrecision ?? .fp32) == .fp32 else { throw VivoProteinStressError.invalid("checkpoint physical identity") }
            if let prior, let priorState {
                guard prior.action != .rejected, request.stages.indices.contains(value.stageIndex),
                      state.configurationFingerprint == (try compiled.configurationFingerprint(stage: value.stageIndex)) else {
                    throw VivoProteinStressError.invalid("journal continuation or stage configuration")
                }
                if value.action == .stageStart {
                    guard value.stageIndex == prior.stageIndex + 1, value.stepsInStage == 0, value.rejection == nil,
                          prior.stageIndex == -1 || prior.stepsInStage == request.stages[prior.stageIndex].steps else {
                        throw VivoProteinStressError.invalid("premature or reordered stage transition")
                    }
                    let transfer = try compiled.transfer(priorState, from: prior.stageIndex >= 0 ? prior.stageIndex : nil, to: value.stageIndex)
                    guard state == transfer.checkpoint,
                          value.cumulativeProtocolWorkKJPerMol == prior.cumulativeProtocolWorkKJPerMol + transfer.workKJPerMol else {
                        throw VivoProteinStressError.invalid("stage transfer changed physical state, RNG counter, or work")
                    }
                } else {
                    guard value.action == .sample || value.action == .rejected,
                          value.stageIndex == prior.stageIndex, value.stepsInStage >= prior.stepsInStage,
                          prior.stepsInStage < request.stages[value.stageIndex].steps,
                          value.stepsInStage <= request.stages[value.stageIndex].steps,
                          value.cumulativeProtocolWorkKJPerMol == prior.cumulativeProtocolWorkKJPerMol else {
                        throw VivoProteinStressError.invalid("invalid accepted/rejected progress")
                    }
                    let delta = value.stepsInStage - prior.stepsInStage
                    let expectedStop = prior.stepsInStage + min(request.stages[value.stageIndex].steps - prior.stepsInStage,
                                                                 request.sampleEvery - prior.stepsInStage % request.sampleEvery)
                    if value.action == .sample {
                        guard value.stepsInStage == expectedStop, delta > 0, value.rejection == nil else {
                            throw VivoProteinStressError.invalid("sample differs from declared grid")
                        }
                    } else {
                        guard value.stepsInStage < expectedStop, let rejected = value.rejection,
                              !rejected.committed, rejected.statusFlags != 0,
                              rejected.systemFingerprint == state.systemFingerprint,
                              rejected.configurationFingerprint == state.configurationFingerprint,
                              rejected.timeBeforePS == state.timePS, rejected.timeAfterPS == state.timePS else {
                            throw VivoProteinStressError.invalid("rejection certificate or accepted prefix")
                        }
                        if delta == 0, state != priorState { throw VivoProteinStressError.invalid("rejected step changed accepted state") }
                    }
                    var expectedTime = priorState.timePS
                    for _ in 0..<delta { expectedTime += request.sourceConfiguration.timeStepPS }
                    guard state.acceptedStep == priorState.acceptedStep + delta, state.timePS == expectedTime else {
                        throw VivoProteinStressError.invalid("accepted clock or random-step discontinuity")
                    }
                }
            } else {
                guard value.action == .source, value.stageIndex == -1, value.stepsInStage == 0, value.previous == nil,
                      value.cumulativeProtocolWorkKJPerMol == 0, value.rejection == nil, state == request.sourceCheckpoint else {
                    throw VivoProteinStressError.invalid("journal does not begin at the exact prepared source")
                }
            }
            let expected = try compiled.frame(checkpoint: state, stage: value.stageIndex >= 0 ? value.stageIndex : nil)
            guard value.frame == expected else { throw VivoProteinStressError.invalid("stored metrics do not reconstruct from checkpoint") }
            result.append(.init(fingerprint: id, entry: value, acceptedStep: state.acceptedStep, timePS: state.timePS))
            prior = value; priorState = state
        }
        guard !result.isEmpty else { throw VivoProteinStressError.invalid("empty journal") }
        return result
    }

    static func put<T: Encodable & Sendable>(_ value: T, kind: String, store: VivoArtifactStore) async throws -> VivoStoredArtifact {
        try await store.put(data: VivoCanonicalJSON.encode(value), kind: kind, mediaType: "application/json")
    }
    static func read<T: Decodable & Sendable>(_ type: T.Type, id: VivoFingerprint, kind: String, store: VivoArtifactStore) async throws -> T {
        let descriptor = try await store.descriptor(for: id)
        guard descriptor.kind == kind, descriptor.mediaType == "application/json" else {
            throw VivoProteinStressError.invalid("artifact kind differs from requested journal object")
        }
        let bytes = try await store.data(for: id, maximumBytes: 128 * 1024 * 1024)
        return try VivoCanonicalJSON.decode(type, from: bytes)
    }
}
