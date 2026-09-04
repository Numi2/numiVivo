import CryptoKit
import Foundation

public enum VivoCampaignProgressEventKind: String, Codable, CaseIterable, Sendable {
    case claimed
    case checkpointed
    case committed
    case rejected
    case failed
    case abandoned
}

public enum VivoCampaignProgressJobState: String, Codable, CaseIterable, Sendable {
    case pending
    case active
    case committed
    case rejected
    case failed
}

public struct VivoCampaignJobClaim: Codable, Equatable, Sendable {
    public let jobIndex: UInt64
    public let jobDigest: String
    public let workerID: String
    public let claimID: String
    public let attempt: UInt32
    public let eventIndex: UInt64
}

public struct VivoCampaignProgressEvent: Codable, Equatable, Sendable {
    public let index: UInt64
    public let kind: VivoCampaignProgressEventKind
    public let jobIndex: UInt64
    public let jobDigest: String
    public let workerID: String
    public let claimID: String
    public let attempt: UInt32
    public let resultFingerprint: String?
    public let certificateFingerprint: String?
    public let checkpointFingerprint: String?
    public let message: String?
    public let wallClockUnixSeconds: Double
    public let previousDigest: String
    public let digest: String
}

public struct VivoCampaignProgressSnapshot: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let campaignFingerprint: String
    public let campaignLedgerHead: String
    public let jobCount: UInt64
    public let events: [VivoCampaignProgressEvent]
    public let journalHead: String
    public let fingerprint: String
}

public struct VivoCampaignProgressSummary: Codable, Equatable, Sendable {
    public let pending: UInt64
    public let active: UInt64
    public let committed: UInt64
    public let rejected: UInt64
    public let failed: UInt64
    public let eventCount: UInt64
    public let journalHead: String
}

public enum VivoCampaignResumeError: Error, LocalizedError, Sendable {
    case invalidCampaign
    case invalidWorker
    case invalidFingerprint(String)
    case unknownJob(UInt64)
    case jobDigestMismatch(UInt64)
    case claimNotActive(UInt64)
    case claimMismatch(UInt64)
    case terminalJob(UInt64)
    case invalidSnapshot(String)
    case journalMismatch(UInt64, String)
    case sourceFile(String)
    case destinationFile(String)
    case arithmeticOverflow
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidCampaign: return "campaign manifest is invalid for resumption"
        case .invalidWorker: return "campaign worker identifier is invalid"
        case .invalidFingerprint(let subject): return "campaign resume fingerprint is invalid: \(subject)"
        case .unknownJob(let job): return "campaign resume references unknown job \(job)"
        case .jobDigestMismatch(let job): return "campaign resume job digest mismatch for job \(job)"
        case .claimNotActive(let job): return "campaign job \(job) does not have an active claim"
        case .claimMismatch(let job): return "campaign claim token does not own job \(job)"
        case .terminalJob(let job): return "campaign job \(job) is already terminal"
        case .invalidSnapshot(let reason): return "campaign progress snapshot is invalid: \(reason)"
        case .journalMismatch(let event, let reason): return "campaign journal mismatch at event \(event): \(reason)"
        case .sourceFile(let reason): return "campaign progress source file failed: \(reason)"
        case .destinationFile(let reason): return "campaign progress destination file failed: \(reason)"
        case .arithmeticOverflow: return "campaign progress arithmetic overflow"
        case .encoding(let reason): return "campaign progress encoding failed: \(reason)"
        }
    }
}

private struct VivoCampaignReplayState {
    struct Job {
        var state: VivoCampaignProgressJobState = .pending
        var activeClaim: VivoCampaignJobClaim?
        var attempts: UInt32 = 0
        var lastCheckpointFingerprint: String?
    }

    var jobs: [UInt64: Job]
    var journalHead: String
}

public struct VivoCampaignProgressJournal: Sendable {
    public let campaign: VivoCampaignManifest
    public private(set) var snapshot: VivoCampaignProgressSnapshot
    private var replay: VivoCampaignReplayState

    public init(campaign: VivoCampaignManifest) throws {
        try VivoValidatedArtifactLoader.validateCampaignManifest(campaign)
        guard campaign.jobs.count <= Int(UInt64.max) else {
            throw VivoCampaignResumeError.arithmeticOverflow
        }
        self.campaign = campaign
        self.replay = .init(
            jobs: Dictionary(uniqueKeysWithValues: campaign.jobs.map {
                ($0.index, VivoCampaignReplayState.Job())
            }),
            journalHead: Self.zeroDigest
        )
        self.snapshot = try Self.makeSnapshot(
            campaign: campaign,
            events: [],
            journalHead: Self.zeroDigest
        )
    }

    public init(
        campaign: VivoCampaignManifest,
        snapshot: VivoCampaignProgressSnapshot
    ) throws {
        try VivoValidatedArtifactLoader.validateCampaignManifest(campaign)
        let replay = try Self.replay(campaign: campaign, snapshot: snapshot)
        self.campaign = campaign
        self.snapshot = snapshot
        self.replay = replay
    }

    public func summary() -> VivoCampaignProgressSummary {
        var counts: [VivoCampaignProgressJobState: UInt64] = [:]
        for job in replay.jobs.values { counts[job.state, default: 0] += 1 }
        return .init(
            pending: counts[.pending, default: 0],
            active: counts[.active, default: 0],
            committed: counts[.committed, default: 0],
            rejected: counts[.rejected, default: 0],
            failed: counts[.failed, default: 0],
            eventCount: UInt64(snapshot.events.count),
            journalHead: snapshot.journalHead
        )
    }

    public func state(jobIndex: UInt64) throws -> VivoCampaignProgressJobState {
        guard let job = replay.jobs[jobIndex] else {
            throw VivoCampaignResumeError.unknownJob(jobIndex)
        }
        return job.state
    }

    public func activeClaims() -> [VivoCampaignJobClaim] {
        replay.jobs.values.compactMap(\.activeClaim).sorted { $0.jobIndex < $1.jobIndex }
    }

    public func nextPendingJobIndex() -> UInt64? {
        campaign.jobs.first(where: { replay.jobs[$0.index]?.state == .pending })?.index
    }

    public mutating func claimNext(
        workerID: String,
        wallClock: Date = Date()
    ) throws -> VivoCampaignJobClaim? {
        guard !workerID.isEmpty, workerID.utf8.count <= 256 else {
            throw VivoCampaignResumeError.invalidWorker
        }
        guard let campaignJob = campaign.jobs.first(where: {
            replay.jobs[$0.index]?.state == .pending
        }) else {
            return nil
        }
        guard var state = replay.jobs[campaignJob.index] else {
            throw VivoCampaignResumeError.unknownJob(campaignJob.index)
        }
        let nextAttempt = state.attempts.addingReportingOverflow(1)
        guard !nextAttempt.overflow else { throw VivoCampaignResumeError.arithmeticOverflow }
        let eventIndex = UInt64(snapshot.events.count)
        let claimID = Self.claimDigest(
            campaignFingerprint: campaign.fingerprint,
            jobDigest: campaignJob.digest,
            workerID: workerID,
            attempt: nextAttempt.partialValue,
            eventIndex: eventIndex
        )
        let event = try makeEvent(
            kind: .claimed,
            job: campaignJob,
            workerID: workerID,
            claimID: claimID,
            attempt: nextAttempt.partialValue,
            wallClock: wallClock
        )
        let claim = VivoCampaignJobClaim(
            jobIndex: campaignJob.index,
            jobDigest: campaignJob.digest,
            workerID: workerID,
            claimID: claimID,
            attempt: nextAttempt.partialValue,
            eventIndex: event.index
        )
        state.state = .active
        state.activeClaim = claim
        state.attempts = nextAttempt.partialValue
        replay.jobs[campaignJob.index] = state
        try append(event)
        return claim
    }

    public mutating func checkpoint(
        claim: VivoCampaignJobClaim,
        checkpointFingerprint: String,
        wallClock: Date = Date()
    ) throws {
        try validateFingerprint(checkpointFingerprint, subject: "checkpoint")
        let job = try activeCampaignJob(for: claim)
        let event = try makeEvent(
            kind: .checkpointed,
            job: job,
            workerID: claim.workerID,
            claimID: claim.claimID,
            attempt: claim.attempt,
            checkpointFingerprint: checkpointFingerprint,
            wallClock: wallClock
        )
        replay.jobs[job.index]?.lastCheckpointFingerprint = checkpointFingerprint
        try append(event)
    }

    public mutating func commit(
        claim: VivoCampaignJobClaim,
        resultFingerprint: String,
        certificateFingerprint: String,
        checkpointFingerprint: String? = nil,
        message: String? = nil,
        wallClock: Date = Date()
    ) throws {
        try terminate(
            claim: claim,
            kind: .committed,
            resultFingerprint: resultFingerprint,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClock: wallClock
        )
    }

    public mutating func reject(
        claim: VivoCampaignJobClaim,
        certificateFingerprint: String,
        checkpointFingerprint: String? = nil,
        message: String? = nil,
        wallClock: Date = Date()
    ) throws {
        try terminate(
            claim: claim,
            kind: .rejected,
            resultFingerprint: nil,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClock: wallClock
        )
    }

    public mutating func fail(
        claim: VivoCampaignJobClaim,
        message: String,
        checkpointFingerprint: String? = nil,
        wallClock: Date = Date()
    ) throws {
        try terminate(
            claim: claim,
            kind: .failed,
            resultFingerprint: nil,
            certificateFingerprint: nil,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClock: wallClock
        )
    }

    public mutating func abandon(
        claim: VivoCampaignJobClaim,
        message: String,
        wallClock: Date = Date()
    ) throws {
        let job = try activeCampaignJob(for: claim)
        let event = try makeEvent(
            kind: .abandoned,
            job: job,
            workerID: claim.workerID,
            claimID: claim.claimID,
            attempt: claim.attempt,
            message: message,
            wallClock: wallClock
        )
        replay.jobs[job.index]?.state = .pending
        replay.jobs[job.index]?.activeClaim = nil
        try append(event)
    }

    public mutating func abandonAllActive(
        message: String,
        wallClock: Date = Date()
    ) throws {
        for claim in activeClaims() {
            try abandon(claim: claim, message: message, wallClock: wallClock)
        }
    }

    private mutating func terminate(
        claim: VivoCampaignJobClaim,
        kind: VivoCampaignProgressEventKind,
        resultFingerprint: String?,
        certificateFingerprint: String?,
        checkpointFingerprint: String?,
        message: String?,
        wallClock: Date
    ) throws {
        if let resultFingerprint {
            try validateFingerprint(resultFingerprint, subject: "result")
        }
        if let certificateFingerprint {
            try validateFingerprint(certificateFingerprint, subject: "certificate")
        }
        if let checkpointFingerprint {
            try validateFingerprint(checkpointFingerprint, subject: "checkpoint")
        }
        if let message, message.utf8.count > 4096 {
            throw VivoCampaignResumeError.invalidSnapshot("event message exceeds 4096 bytes")
        }
        let job = try activeCampaignJob(for: claim)
        let event = try makeEvent(
            kind: kind,
            job: job,
            workerID: claim.workerID,
            claimID: claim.claimID,
            attempt: claim.attempt,
            resultFingerprint: resultFingerprint,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClock: wallClock
        )
        switch kind {
        case .committed: replay.jobs[job.index]?.state = .committed
        case .rejected: replay.jobs[job.index]?.state = .rejected
        case .failed: replay.jobs[job.index]?.state = .failed
        default: throw VivoCampaignResumeError.invalidSnapshot("invalid terminal event")
        }
        replay.jobs[job.index]?.activeClaim = nil
        if let checkpointFingerprint {
            replay.jobs[job.index]?.lastCheckpointFingerprint = checkpointFingerprint
        }
        try append(event)
    }

    private func activeCampaignJob(
        for claim: VivoCampaignJobClaim
    ) throws -> VivoCampaignJob {
        guard let campaignJob = campaign.jobs.first(where: { $0.index == claim.jobIndex }) else {
            throw VivoCampaignResumeError.unknownJob(claim.jobIndex)
        }
        guard campaignJob.digest == claim.jobDigest else {
            throw VivoCampaignResumeError.jobDigestMismatch(claim.jobIndex)
        }
        guard let state = replay.jobs[claim.jobIndex], state.state == .active,
              let active = state.activeClaim else {
            throw VivoCampaignResumeError.claimNotActive(claim.jobIndex)
        }
        guard active == claim else {
            throw VivoCampaignResumeError.claimMismatch(claim.jobIndex)
        }
        return campaignJob
    }

    private func makeEvent(
        kind: VivoCampaignProgressEventKind,
        job: VivoCampaignJob,
        workerID: String,
        claimID: String,
        attempt: UInt32,
        resultFingerprint: String? = nil,
        certificateFingerprint: String? = nil,
        checkpointFingerprint: String? = nil,
        message: String? = nil,
        wallClock: Date
    ) throws -> VivoCampaignProgressEvent {
        let index = UInt64(snapshot.events.count)
        let wallClockUnixSeconds = wallClock.timeIntervalSince1970
        guard wallClockUnixSeconds.isFinite else {
            throw VivoCampaignResumeError.invalidSnapshot("wall clock is non-finite")
        }
        let digest = Self.eventDigest(
            previous: snapshot.journalHead,
            index: index,
            kind: kind,
            jobIndex: job.index,
            jobDigest: job.digest,
            workerID: workerID,
            claimID: claimID,
            attempt: attempt,
            resultFingerprint: resultFingerprint,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClockUnixSeconds: wallClockUnixSeconds
        )
        return .init(
            index: index,
            kind: kind,
            jobIndex: job.index,
            jobDigest: job.digest,
            workerID: workerID,
            claimID: claimID,
            attempt: attempt,
            resultFingerprint: resultFingerprint,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message,
            wallClockUnixSeconds: wallClockUnixSeconds,
            previousDigest: snapshot.journalHead,
            digest: digest
        )
    }

    private mutating func append(_ event: VivoCampaignProgressEvent) throws {
        var events = snapshot.events
        events.append(event)
        snapshot = try Self.makeSnapshot(
            campaign: campaign,
            events: events,
            journalHead: event.digest
        )
        replay.journalHead = event.digest
    }

    private static func makeSnapshot(
        campaign: VivoCampaignManifest,
        events: [VivoCampaignProgressEvent],
        journalHead: String
    ) throws -> VivoCampaignProgressSnapshot {
        let unsigned = VivoCampaignProgressSnapshot(
            schemaVersion: 1,
            campaignFingerprint: campaign.fingerprint,
            campaignLedgerHead: campaign.ledgerHead,
            jobCount: UInt64(campaign.jobs.count),
            events: events,
            journalHead: journalHead,
            fingerprint: ""
        )
        return .init(
            schemaVersion: unsigned.schemaVersion,
            campaignFingerprint: unsigned.campaignFingerprint,
            campaignLedgerHead: unsigned.campaignLedgerHead,
            jobCount: unsigned.jobCount,
            events: unsigned.events,
            journalHead: unsigned.journalHead,
            fingerprint: try canonicalFingerprint(unsigned)
        )
    }

    private static func replay(
        campaign: VivoCampaignManifest,
        snapshot: VivoCampaignProgressSnapshot
    ) throws -> VivoCampaignReplayState {
        guard snapshot.schemaVersion == 1,
              snapshot.campaignFingerprint == campaign.fingerprint,
              snapshot.campaignLedgerHead == campaign.ledgerHead,
              snapshot.jobCount == UInt64(campaign.jobs.count),
              isFingerprint(snapshot.journalHead),
              isFingerprint(snapshot.fingerprint) else {
            throw VivoCampaignResumeError.invalidSnapshot("identity or version")
        }
        let unsigned = VivoCampaignProgressSnapshot(
            schemaVersion: snapshot.schemaVersion,
            campaignFingerprint: snapshot.campaignFingerprint,
            campaignLedgerHead: snapshot.campaignLedgerHead,
            jobCount: snapshot.jobCount,
            events: snapshot.events,
            journalHead: snapshot.journalHead,
            fingerprint: ""
        )
        let actualSnapshotFingerprint = try canonicalFingerprint(unsigned)
        guard actualSnapshotFingerprint == snapshot.fingerprint else {
            throw VivoCampaignResumeError.invalidFingerprint("progress snapshot")
        }

        let campaignJobs = Dictionary(uniqueKeysWithValues: campaign.jobs.map { ($0.index, $0) })
        var state = VivoCampaignReplayState(
            jobs: Dictionary(uniqueKeysWithValues: campaign.jobs.map {
                ($0.index, VivoCampaignReplayState.Job())
            }),
            journalHead: zeroDigest
        )
        var previous = zeroDigest

        for (expectedIndex, event) in snapshot.events.enumerated() {
            guard event.index == UInt64(expectedIndex),
                  event.previousDigest == previous,
                  let campaignJob = campaignJobs[event.jobIndex],
                  event.jobDigest == campaignJob.digest,
                  isFingerprint(event.digest),
                  isFingerprint(event.claimID),
                  event.attempt > 0,
                  !event.workerID.isEmpty,
                  event.workerID.utf8.count <= 256,
                  event.wallClockUnixSeconds.isFinite else {
                throw VivoCampaignResumeError.journalMismatch(
                    event.index,
                    "index, identity, claim or wall clock"
                )
            }
            let expectedDigest = eventDigest(
                previous: previous,
                index: event.index,
                kind: event.kind,
                jobIndex: event.jobIndex,
                jobDigest: event.jobDigest,
                workerID: event.workerID,
                claimID: event.claimID,
                attempt: event.attempt,
                resultFingerprint: event.resultFingerprint,
                certificateFingerprint: event.certificateFingerprint,
                checkpointFingerprint: event.checkpointFingerprint,
                message: event.message,
                wallClockUnixSeconds: event.wallClockUnixSeconds
            )
            guard expectedDigest == event.digest else {
                throw VivoCampaignResumeError.journalMismatch(event.index, "digest")
            }
            guard var job = state.jobs[event.jobIndex] else {
                throw VivoCampaignResumeError.unknownJob(event.jobIndex)
            }
            switch event.kind {
            case .claimed:
                guard job.state == .pending,
                      event.attempt == job.attempts + 1,
                      event.claimID == claimDigest(
                        campaignFingerprint: campaign.fingerprint,
                        jobDigest: event.jobDigest,
                        workerID: event.workerID,
                        attempt: event.attempt,
                        eventIndex: event.index
                      ) else {
                    throw VivoCampaignResumeError.journalMismatch(event.index, "claim transition")
                }
                job.state = .active
                job.attempts = event.attempt
                job.activeClaim = .init(
                    jobIndex: event.jobIndex,
                    jobDigest: event.jobDigest,
                    workerID: event.workerID,
                    claimID: event.claimID,
                    attempt: event.attempt,
                    eventIndex: event.index
                )
            case .checkpointed:
                guard job.state == .active,
                      job.activeClaim?.claimID == event.claimID,
                      let checkpoint = event.checkpointFingerprint,
                      isFingerprint(checkpoint) else {
                    throw VivoCampaignResumeError.journalMismatch(event.index, "checkpoint transition")
                }
                job.lastCheckpointFingerprint = checkpoint
            case .committed, .rejected, .failed:
                guard job.state == .active,
                      job.activeClaim?.claimID == event.claimID else {
                    throw VivoCampaignResumeError.journalMismatch(event.index, "terminal transition")
                }
                if event.kind == .committed {
                    guard let result = event.resultFingerprint,
                          let certificate = event.certificateFingerprint,
                          isFingerprint(result),
                          isFingerprint(certificate) else {
                        throw VivoCampaignResumeError.journalMismatch(event.index, "committed fingerprints")
                    }
                    job.state = .committed
                } else if event.kind == .rejected {
                    guard let certificate = event.certificateFingerprint,
                          isFingerprint(certificate) else {
                        throw VivoCampaignResumeError.journalMismatch(event.index, "rejection certificate")
                    }
                    job.state = .rejected
                } else {
                    guard let message = event.message, !message.isEmpty else {
                        throw VivoCampaignResumeError.journalMismatch(event.index, "failure message")
                    }
                    job.state = .failed
                }
                if let checkpoint = event.checkpointFingerprint {
                    guard isFingerprint(checkpoint) else {
                        throw VivoCampaignResumeError.journalMismatch(event.index, "terminal checkpoint")
                    }
                    job.lastCheckpointFingerprint = checkpoint
                }
                job.activeClaim = nil
            case .abandoned:
                guard job.state == .active,
                      job.activeClaim?.claimID == event.claimID else {
                    throw VivoCampaignResumeError.journalMismatch(event.index, "abandon transition")
                }
                job.state = .pending
                job.activeClaim = nil
            }
            state.jobs[event.jobIndex] = job
            previous = event.digest
        }
        guard previous == snapshot.journalHead else {
            throw VivoCampaignResumeError.invalidSnapshot("journal head")
        }
        state.journalHead = previous
        return state
    }

    private func validateFingerprint(_ value: String, subject: String) throws {
        guard Self.isFingerprint(value) else {
            throw VivoCampaignResumeError.invalidFingerprint(subject)
        }
    }

    private static func claimDigest(
        campaignFingerprint: String,
        jobDigest: String,
        workerID: String,
        attempt: UInt32,
        eventIndex: UInt64
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-campaign-claim-v1".utf8))
        hasher.update(data: Data(campaignFingerprint.utf8))
        hasher.update(data: Data(jobDigest.utf8))
        hasher.update(data: Data(workerID.utf8))
        update(&hasher, attempt)
        update(&hasher, eventIndex)
        return hex(hasher.finalize())
    }

    private static func eventDigest(
        previous: String,
        index: UInt64,
        kind: VivoCampaignProgressEventKind,
        jobIndex: UInt64,
        jobDigest: String,
        workerID: String,
        claimID: String,
        attempt: UInt32,
        resultFingerprint: String?,
        certificateFingerprint: String?,
        checkpointFingerprint: String?,
        message: String?,
        wallClockUnixSeconds: Double
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-campaign-progress-event-v1".utf8))
        hasher.update(data: Data(previous.utf8))
        update(&hasher, index)
        hasher.update(data: Data(kind.rawValue.utf8))
        update(&hasher, jobIndex)
        hasher.update(data: Data(jobDigest.utf8))
        hasher.update(data: Data(workerID.utf8))
        hasher.update(data: Data(claimID.utf8))
        update(&hasher, attempt)
        updateOptional(&hasher, resultFingerprint)
        updateOptional(&hasher, certificateFingerprint)
        updateOptional(&hasher, checkpointFingerprint)
        updateOptional(&hasher, message)
        update(&hasher, wallClockUnixSeconds.bitPattern)
        return hex(hasher.finalize())
    }

    private static func canonicalFingerprint<T: Encodable>(_ value: T) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        do {
            return hex(SHA256.hash(data: try encoder.encode(value)))
        } catch {
            throw VivoCampaignResumeError.encoding(error.localizedDescription)
        }
    }

    private static func updateOptional(_ hasher: inout SHA256, _ value: String?) {
        if let value {
            hasher.update(data: Data([1]))
            update(&hasher, UInt64(value.utf8.count))
            hasher.update(data: Data(value.utf8))
        } else {
            hasher.update(data: Data([0]))
        }
    }

    private static func update<T: FixedWidthInteger>(_ hasher: inout SHA256, _ value: T) {
        var littleEndian = value.littleEndian
        Swift.withUnsafeBytes(of: &littleEndian) { bytes in
            hasher.update(data: Data(bytes))
        }
    }

    private static func isFingerprint(_ value: String) -> Bool {
        value.utf8.count == 64 && value.allSatisfy {
            $0.isHexDigit && !$0.isUppercase
        }
    }

    private static func hex<S: Sequence>(_ bytes: S) -> String where S.Element == UInt8 {
        bytes.map { String(format: "%02x", $0) }.joined()
    }

    private static let zeroDigest = String(repeating: "0", count: 64)
}

public actor VivoCampaignResumeStore {
    public nonisolated let url: URL
    public nonisolated let campaignFingerprint: String
    private var journal: VivoCampaignProgressJournal

    public static func openOrCreate(
        campaign: VivoCampaignManifest,
        at url: URL,
        abandonActiveClaimsOnOpen: Bool = false
    ) throws -> VivoCampaignResumeStore {
        let manager = FileManager.default
        let journal: VivoCampaignProgressJournal
        if manager.fileExists(atPath: url.path) {
            let data: Data
            do {
                data = try VivoValidatedArtifactLoader.data(at: url)
            } catch {
                throw VivoCampaignResumeError.sourceFile(error.localizedDescription)
            }
            let snapshot: VivoCampaignProgressSnapshot
            do {
                snapshot = try JSONDecoder().decode(
                    VivoCampaignProgressSnapshot.self,
                    from: data
                )
            } catch {
                throw VivoCampaignResumeError.sourceFile(error.localizedDescription)
            }
            var resumed = try VivoCampaignProgressJournal(
                campaign: campaign,
                snapshot: snapshot
            )
            if abandonActiveClaimsOnOpen {
                try resumed.abandonAllActive(
                    message: "active claim recovered after process restart"
                )
                try persist(resumed.snapshot, to: url)
            }
            journal = resumed
        } else {
            journal = try VivoCampaignProgressJournal(campaign: campaign)
            try persist(journal.snapshot, to: url)
        }
        return VivoCampaignResumeStore(
            url: url,
            campaignFingerprint: campaign.fingerprint,
            journal: journal
        )
    }

    private init(
        url: URL,
        campaignFingerprint: String,
        journal: VivoCampaignProgressJournal
    ) {
        self.url = url
        self.campaignFingerprint = campaignFingerprint
        self.journal = journal
    }

    public func snapshot() -> VivoCampaignProgressSnapshot { journal.snapshot }
    public func summary() -> VivoCampaignProgressSummary { journal.summary() }
    public func activeClaims() -> [VivoCampaignJobClaim] { journal.activeClaims() }

    public func claimNext(workerID: String) throws -> VivoCampaignJobClaim? {
        var candidate = journal
        let claim = try candidate.claimNext(workerID: workerID)
        if claim != nil {
            try Self.persist(candidate.snapshot, to: url)
            journal = candidate
        }
        return claim
    }

    public func checkpoint(
        claim: VivoCampaignJobClaim,
        checkpointFingerprint: String
    ) throws {
        var candidate = journal
        try candidate.checkpoint(
            claim: claim,
            checkpointFingerprint: checkpointFingerprint
        )
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    public func commit(
        claim: VivoCampaignJobClaim,
        resultFingerprint: String,
        certificateFingerprint: String,
        checkpointFingerprint: String? = nil,
        message: String? = nil
    ) throws {
        var candidate = journal
        try candidate.commit(
            claim: claim,
            resultFingerprint: resultFingerprint,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message
        )
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    public func reject(
        claim: VivoCampaignJobClaim,
        certificateFingerprint: String,
        checkpointFingerprint: String? = nil,
        message: String? = nil
    ) throws {
        var candidate = journal
        try candidate.reject(
            claim: claim,
            certificateFingerprint: certificateFingerprint,
            checkpointFingerprint: checkpointFingerprint,
            message: message
        )
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    public func fail(
        claim: VivoCampaignJobClaim,
        message: String,
        checkpointFingerprint: String? = nil
    ) throws {
        var candidate = journal
        try candidate.fail(
            claim: claim,
            message: message,
            checkpointFingerprint: checkpointFingerprint
        )
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    public func abandon(
        claim: VivoCampaignJobClaim,
        message: String
    ) throws {
        var candidate = journal
        try candidate.abandon(claim: claim, message: message)
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    public func abandonAllActive(message: String) throws {
        var candidate = journal
        try candidate.abandonAllActive(message: message)
        try Self.persist(candidate.snapshot, to: url)
        journal = candidate
    }

    private static func persist(
        _ snapshot: VivoCampaignProgressSnapshot,
        to url: URL
    ) throws {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.prettyPrinted, .sortedKeys, .withoutEscapingSlashes]
        let data: Data
        do { data = try encoder.encode(snapshot) }
        catch { throw VivoCampaignResumeError.encoding(error.localizedDescription) }
        do {
            try FileManager.default.createDirectory(
                at: url.deletingLastPathComponent(),
                withIntermediateDirectories: true
            )
            try data.write(to: url, options: [.atomic])
        } catch {
            throw VivoCampaignResumeError.destinationFile(error.localizedDescription)
        }
    }
}
