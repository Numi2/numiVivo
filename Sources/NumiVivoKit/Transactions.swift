import Foundation

public struct VivoTransactionEpoch: Codable, Sendable, Hashable {
    public var id: UUID
    public var sequence: UInt64
    public var simulationTime: Double
    public var deltaTime: Double
    public var programContentFingerprint: String
    public var metadata: [String: String]

    public init(
        id: UUID = UUID(),
        sequence: UInt64,
        simulationTime: Double,
        deltaTime: Double,
        programContentFingerprint: String,
        metadata: [String: String] = [:]
    ) {
        self.id = id
        self.sequence = sequence
        self.simulationTime = simulationTime
        self.deltaTime = deltaTime
        self.programContentFingerprint = programContentFingerprint
        self.metadata = metadata
    }
}

public enum VivoPrepareVote: String, Codable, Sendable {
    case prepared
    case reject
    case substepRequired
    case reversibleShutdown
    case permanentShutdown
}

public struct VivoPrepareReceipt: Codable, Sendable, Hashable {
    public var participantID: String
    public var epochID: UUID
    public var vote: VivoPrepareVote
    public var token: Data
    public var preparedStateSHA256: String
    public var diagnosticFlags: UInt32
    public var requestedMaximumDeltaTime: Double?
    public var metadata: [String: String]

    public init(
        participantID: String,
        epochID: UUID,
        vote: VivoPrepareVote,
        token: Data = Data(),
        preparedStateSHA256: String,
        diagnosticFlags: UInt32 = 0,
        requestedMaximumDeltaTime: Double? = nil,
        metadata: [String: String] = [:]
    ) {
        self.participantID = participantID
        self.epochID = epochID
        self.vote = vote
        self.token = token
        self.preparedStateSHA256 = preparedStateSHA256
        self.diagnosticFlags = diagnosticFlags
        self.requestedMaximumDeltaTime = requestedMaximumDeltaTime
        self.metadata = metadata
    }
}

public struct VivoTransactionalParticipant: Sendable {
    public typealias Prepare = @Sendable (VivoTransactionEpoch) async throws -> VivoPrepareReceipt
    public typealias Commit = @Sendable (VivoTransactionEpoch, VivoPrepareReceipt) async throws -> Void
    public typealias Rollback = @Sendable (VivoTransactionEpoch, VivoPrepareReceipt) async -> Void

    public let id: String
    private let prepareOperation: Prepare
    private let commitOperation: Commit
    private let rollbackOperation: Rollback

    public init(
        id: String,
        prepare: @escaping Prepare,
        commit: @escaping Commit,
        rollback: @escaping Rollback
    ) {
        self.id = id
        self.prepareOperation = prepare
        self.commitOperation = commit
        self.rollbackOperation = rollback
    }

    func prepare(_ epoch: VivoTransactionEpoch) async throws -> VivoPrepareReceipt {
        var receipt = try await prepareOperation(epoch)
        guard receipt.participantID == id else {
            throw VivoTransactionError.invalidReceipt(
                participant: id,
                reason: "receipt participantID is '\(receipt.participantID)'"
            )
        }
        guard receipt.epochID == epoch.id else {
            throw VivoTransactionError.invalidReceipt(
                participant: id,
                reason: "receipt epoch does not match the active transaction"
            )
        }
        guard !receipt.preparedStateSHA256.isEmpty else {
            throw VivoTransactionError.invalidReceipt(
                participant: id,
                reason: "prepared-state fingerprint is empty"
            )
        }
        receipt.participantID = id
        return receipt
    }

    func commit(_ epoch: VivoTransactionEpoch, receipt: VivoPrepareReceipt) async throws {
        try await commitOperation(epoch, receipt)
    }

    func rollback(_ epoch: VivoTransactionEpoch, receipt: VivoPrepareReceipt) async {
        await rollbackOperation(epoch, receipt)
    }
}

public struct VivoTransactionOutcome: Codable, Sendable, Hashable {
    public var epoch: VivoTransactionEpoch
    public var committed: Bool
    public var receipts: [VivoPrepareReceipt]
    public var limitingDeltaTime: Double?
    public var terminalVote: VivoPrepareVote?
}

public struct VivoInDoubtTransaction: Codable, Sendable, Hashable {
    public var epoch: VivoTransactionEpoch
    public var committedParticipantIDs: [String]
    public var unresolvedReceipts: [VivoPrepareReceipt]
    public var failure: String
}

public enum VivoTransactionError: Error, Sendable, CustomStringConvertible {
    case duplicateParticipant(String)
    case invalidEpoch(String)
    case invalidReceipt(participant: String, reason: String)
    case preparationFailed(participant: String, reason: String)
    case rejected(VivoTransactionOutcome)
    case commitInDoubt(VivoInDoubtTransaction)

    public var description: String {
        switch self {
        case .duplicateParticipant(let participant):
            return "Transaction contains duplicate participant '\(participant)'."
        case .invalidEpoch(let reason):
            return "Transaction epoch is invalid: \(reason)."
        case .invalidReceipt(let participant, let reason):
            return "Participant '\(participant)' returned an invalid prepare receipt: \(reason)."
        case .preparationFailed(let participant, let reason):
            return "Participant '\(participant)' failed during prepare: \(reason)."
        case .rejected(let outcome):
            return "Transaction \(outcome.epoch.id) was not prepared by every participant."
        case .commitInDoubt(let record):
            return "Transaction \(record.epoch.id) entered an in-doubt commit state: \(record.failure)."
        }
    }
}

public actor VivoTransactionCoordinator {
    private var inDoubt: [UUID: VivoInDoubtTransaction] = [:]

    public init() {}

    public func execute(
        epoch: VivoTransactionEpoch,
        participants suppliedParticipants: [VivoTransactionalParticipant]
    ) async throws -> VivoTransactionOutcome {
        guard epoch.simulationTime.isFinite,
              epoch.deltaTime.isFinite,
              epoch.deltaTime > 0,
              !epoch.programContentFingerprint.isEmpty else {
            throw VivoTransactionError.invalidEpoch(
                "time and deltaTime must be finite, deltaTime must be positive, and the program fingerprint is required"
            )
        }
        let participants = suppliedParticipants.sorted { $0.id < $1.id }
        var known = Set<String>()
        for participant in participants {
            guard known.insert(participant.id).inserted else {
                throw VivoTransactionError.duplicateParticipant(participant.id)
            }
        }
        guard !participants.isEmpty else {
            return VivoTransactionOutcome(
                epoch: epoch,
                committed: true,
                receipts: [],
                limitingDeltaTime: nil,
                terminalVote: nil
            )
        }

        var receiptsByID: [String: VivoPrepareReceipt] = [:]
        do {
            try await withThrowingTaskGroup(
                of: (String, VivoPrepareReceipt).self
            ) { group in
                for participant in participants {
                    group.addTask {
                        (participant.id, try await participant.prepare(epoch))
                    }
                }
                while let (id, receipt) = try await group.next() {
                    receiptsByID[id] = receipt
                }
            }
        } catch {
            let prepared = participants.compactMap { participant -> (VivoTransactionalParticipant, VivoPrepareReceipt)? in
                guard let receipt = receiptsByID[participant.id] else { return nil }
                return (participant, receipt)
            }
            await rollback(epoch: epoch, prepared: prepared)
            if let transactionError = error as? VivoTransactionError {
                throw transactionError
            }
            throw VivoTransactionError.preparationFailed(
                participant: "unknown",
                reason: String(describing: error)
            )
        }

        let orderedReceipts = participants.compactMap { receiptsByID[$0.id] }
        let limitingDeltaTime = orderedReceipts
            .compactMap(\.requestedMaximumDeltaTime)
            .filter { $0.isFinite && $0 > 0 }
            .min()
        let terminalVote = Self.strongestVote(in: orderedReceipts)
        let allPrepared = orderedReceipts.count == participants.count &&
                          orderedReceipts.allSatisfy { $0.vote == .prepared }
        guard allPrepared else {
            let prepared = zip(participants, orderedReceipts).filter { $0.1.vote == .prepared }
            await rollback(epoch: epoch, prepared: Array(prepared))
            let outcome = VivoTransactionOutcome(
                epoch: epoch,
                committed: false,
                receipts: orderedReceipts,
                limitingDeltaTime: limitingDeltaTime,
                terminalVote: terminalVote
            )
            throw VivoTransactionError.rejected(outcome)
        }

        var committedIDs: [String] = []
        for participant in participants {
            guard let receipt = receiptsByID[participant.id] else {
                let record = VivoInDoubtTransaction(
                    epoch: epoch,
                    committedParticipantIDs: committedIDs,
                    unresolvedReceipts: orderedReceipts,
                    failure: "missing receipt during commit"
                )
                inDoubt[epoch.id] = record
                throw VivoTransactionError.commitInDoubt(record)
            }
            do {
                try await participant.commit(epoch, receipt: receipt)
                committedIDs.append(participant.id)
            } catch {
                let unresolved = participants
                    .filter { !committedIDs.contains($0.id) }
                    .compactMap { receiptsByID[$0.id] }
                let record = VivoInDoubtTransaction(
                    epoch: epoch,
                    committedParticipantIDs: committedIDs,
                    unresolvedReceipts: unresolved,
                    failure: "participant '\(participant.id)': \(error)"
                )
                inDoubt[epoch.id] = record
                throw VivoTransactionError.commitInDoubt(record)
            }
        }

        inDoubt.removeValue(forKey: epoch.id)
        return VivoTransactionOutcome(
            epoch: epoch,
            committed: true,
            receipts: orderedReceipts,
            limitingDeltaTime: limitingDeltaTime,
            terminalVote: nil
        )
    }

    public func inDoubtTransactions() -> [VivoInDoubtTransaction] {
        inDoubt.values.sorted { $0.epoch.sequence < $1.epoch.sequence }
    }

    public func resolveInDoubt(epochID: UUID) {
        inDoubt.removeValue(forKey: epochID)
    }

    private func rollback(
        epoch: VivoTransactionEpoch,
        prepared: [(VivoTransactionalParticipant, VivoPrepareReceipt)]
    ) async {
        for (participant, receipt) in prepared.sorted(by: { $0.0.id > $1.0.id }) {
            await participant.rollback(epoch, receipt: receipt)
        }
    }

    private static func strongestVote(in receipts: [VivoPrepareReceipt]) -> VivoPrepareVote? {
        let rank: [VivoPrepareVote: Int] = [
            .prepared: 0,
            .substepRequired: 1,
            .reject: 2,
            .reversibleShutdown: 3,
            .permanentShutdown: 4
        ]
        return receipts
            .map(\.vote)
            .filter { $0 != .prepared }
            .max { (rank[$0] ?? 0) < (rank[$1] ?? 0) }
    }
}
