import CryptoKit
import Foundation

/// Stable participant identities used by the integrated NumiLab runtime.
public enum VivoNumiLabParticipant: String, CaseIterable, Codable, Sendable {
    case numanX = "numanx"
    case numiTissue = "numitissue"
    case numiBrain = "numibrain"
    case numiVivo = "numivivo"
    case environment = "environment"
    case experiment = "experiment"
}

/// How a consumer applies one coupling channel to its staged state.
public enum VivoCouplingSemantics: String, Codable, Sendable {
    case replace
    case add
    case minimum
    case maximum
    case rate
    case impulse
    case event
}

/// Description of a versioned channel crossing simulator boundaries.
public struct VivoCouplingChannelDescriptor: Codable, Hashable, Sendable {
    public let id: String
    public let producer: VivoNumiLabParticipant
    public let consumer: VivoNumiLabParticipant
    public let unit: String
    public let semantics: VivoCouplingSemantics
    public let shape: [UInt32]
    public let lowerBound: Float?
    public let upperBound: Float?
    public let required: Bool
    public let schemaVersion: UInt32

    public init(
        id: String,
        producer: VivoNumiLabParticipant,
        consumer: VivoNumiLabParticipant,
        unit: String,
        semantics: VivoCouplingSemantics,
        shape: [UInt32] = [1],
        lowerBound: Float? = nil,
        upperBound: Float? = nil,
        required: Bool = true,
        schemaVersion: UInt32 = 1
    ) {
        self.id = id
        self.producer = producer
        self.consumer = consumer
        self.unit = unit
        self.semantics = semantics
        self.shape = shape
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.required = required
        self.schemaVersion = schemaVersion
    }

    public var elementCount: Int? {
        var count: UInt64 = 1
        for dimension in shape {
            guard dimension > 0 else { return nil }
            let next = count.multipliedReportingOverflow(by: UInt64(dimension))
            guard !next.overflow, next.partialValue <= UInt64(Int.max) else { return nil }
            count = next.partialValue
        }
        return Int(count)
    }

    public func validate() throws {
        guard !id.isEmpty, id.utf8.count <= 256 else {
            throw VivoNumiLabCouplingError.invalidChannel("channel id is empty or too long")
        }
        guard !unit.isEmpty, unit.utf8.count <= 128 else {
            throw VivoNumiLabCouplingError.invalidChannel("channel \(id) has an invalid unit")
        }
        guard !shape.isEmpty, shape.count <= 8, elementCount != nil else {
            throw VivoNumiLabCouplingError.invalidChannel("channel \(id) has an invalid shape")
        }
        if let lowerBound, let upperBound, lowerBound > upperBound {
            throw VivoNumiLabCouplingError.invalidChannel("channel \(id) lower bound exceeds upper bound")
        }
        guard producer != consumer else {
            throw VivoNumiLabCouplingError.invalidChannel("channel \(id) must cross participant boundaries")
        }
    }
}

/// One immutable channel payload. Values are always flattened in row-major order.
public struct VivoCouplingSample: Codable, Equatable, Sendable {
    public let channelID: String
    public let producer: VivoNumiLabParticipant
    public let schemaVersion: UInt32
    public let unit: String
    public let values: [Float]
    public let quality: Float
    public let sourceStep: UInt64
    public let sourceTime: Double

    public init(
        channelID: String,
        producer: VivoNumiLabParticipant,
        schemaVersion: UInt32 = 1,
        unit: String,
        values: [Float],
        quality: Float = 1,
        sourceStep: UInt64,
        sourceTime: Double
    ) {
        self.channelID = channelID
        self.producer = producer
        self.schemaVersion = schemaVersion
        self.unit = unit
        self.values = values
        self.quality = quality
        self.sourceStep = sourceStep
        self.sourceTime = sourceTime
    }
}

public struct VivoCouplingClock: Codable, Equatable, Sendable {
    public let epoch: UInt64
    public let step: UInt64
    public let time: Double
    public let deltaTime: Double

    public init(epoch: UInt64, step: UInt64, time: Double, deltaTime: Double) {
        self.epoch = epoch
        self.step = step
        self.time = time
        self.deltaTime = deltaTime
    }

    public func validate() throws {
        guard time.isFinite, deltaTime.isFinite, deltaTime > 0 else {
            throw VivoNumiLabCouplingError.invalidClock
        }
    }
}

/// Complete immutable input to one distributed shadow step.
public struct VivoCouplingTransaction: Codable, Equatable, Sendable {
    public let id: String
    public let clock: VivoCouplingClock
    public let inputs: [VivoCouplingSample]
    public let iteration: UInt32

    public init(
        id: String,
        clock: VivoCouplingClock,
        inputs: [VivoCouplingSample],
        iteration: UInt32 = 0
    ) {
        self.id = id
        self.clock = clock
        self.inputs = inputs
        self.iteration = iteration
    }
}

public enum VivoPreparationDisposition: String, Codable, Sendable {
    case prepared
    case needsSmallerStep
    case rejected
    case reversibleShutdown
    case permanentShutdown
}

/// Participant-owned candidate result. Outputs are not visible to committed readers.
public struct VivoPreparedCouplingState: Codable, Equatable, Sendable {
    public let participant: VivoNumiLabParticipant
    public let transactionID: String
    public let disposition: VivoPreparationDisposition
    public let maximumAcceptedStep: Double?
    public let outputs: [VivoCouplingSample]
    public let stateDigest: String
    public let diagnostics: [String]

    public init(
        participant: VivoNumiLabParticipant,
        transactionID: String,
        disposition: VivoPreparationDisposition,
        maximumAcceptedStep: Double? = nil,
        outputs: [VivoCouplingSample] = [],
        stateDigest: String,
        diagnostics: [String] = []
    ) {
        self.participant = participant
        self.transactionID = transactionID
        self.disposition = disposition
        self.maximumAcceptedStep = maximumAcceptedStep
        self.outputs = outputs
        self.stateDigest = stateDigest
        self.diagnostics = diagnostics
    }
}

/// A participant must make `stageCommit` fallible and `releaseCommit` infallible.
/// Until `releaseCommit`, committed readers must continue seeing the previous state.
public protocol VivoCoupledParticipant: Sendable {
    var participantID: VivoNumiLabParticipant { get }

    func prepare(_ transaction: VivoCouplingTransaction) async throws -> VivoPreparedCouplingState
    func stageCommit(transactionID: String) async throws
    func releaseCommit(transactionID: String) async
    func rollback(transactionID: String) async
}

public enum VivoJointStepDisposition: String, Codable, Sendable {
    case committed
    case needsSmallerStep
    case rejected
    case reversibleShutdown
    case permanentShutdown
}

public struct VivoJointStepCertificate: Codable, Equatable, Sendable {
    public let transactionID: String
    public let disposition: VivoJointStepDisposition
    public let clock: VivoCouplingClock
    public let requestedStep: Double
    public let maximumAcceptedStep: Double?
    public let participants: [VivoPreparedCouplingState]
    public let routedOutputs: [VivoCouplingSample]
    public let ledgerDigest: String

    public init(
        transactionID: String,
        disposition: VivoJointStepDisposition,
        clock: VivoCouplingClock,
        requestedStep: Double,
        maximumAcceptedStep: Double?,
        participants: [VivoPreparedCouplingState],
        routedOutputs: [VivoCouplingSample],
        ledgerDigest: String
    ) {
        self.transactionID = transactionID
        self.disposition = disposition
        self.clock = clock
        self.requestedStep = requestedStep
        self.maximumAcceptedStep = maximumAcceptedStep
        self.participants = participants
        self.routedOutputs = routedOutputs
        self.ledgerDigest = ledgerDigest
    }
}

public enum VivoNumiLabCouplingError: Error, LocalizedError, Sendable {
    case duplicateParticipant(VivoNumiLabParticipant)
    case missingParticipant(VivoNumiLabParticipant)
    case duplicateChannel(String)
    case invalidChannel(String)
    case invalidClock
    case invalidTransaction(String)
    case participantMismatch
    case undeclaredOutput(String)
    case invalidOutput(String)
    case stagingFailure(String)

    public var errorDescription: String? {
        switch self {
        case .duplicateParticipant(let participant): return "duplicate participant: \(participant.rawValue)"
        case .missingParticipant(let participant): return "missing participant: \(participant.rawValue)"
        case .duplicateChannel(let channel): return "duplicate coupling channel: \(channel)"
        case .invalidChannel(let reason): return "invalid coupling channel: \(reason)"
        case .invalidClock: return "invalid coupling clock"
        case .invalidTransaction(let reason): return "invalid transaction: \(reason)"
        case .participantMismatch: return "participant returned a result for another transaction or identity"
        case .undeclaredOutput(let channel): return "participant produced undeclared channel: \(channel)"
        case .invalidOutput(let reason): return "invalid coupling output: \(reason)"
        case .stagingFailure(let reason): return "joint commit staging failed: \(reason)"
        }
    }
}

/// Deterministic two-phase coordinator for NumanX, NumiTissue, NumiBrain and NumiVivo.
/// Participant order is fixed, output validation is central, and no participant
/// publishes a candidate before every participant has staged the same transaction.
public actor VivoNumiLabCouplingCoordinator {
    private let participants: [VivoNumiLabParticipant: any VivoCoupledParticipant]
    private let orderedParticipants: [VivoNumiLabParticipant]
    private let channels: [String: VivoCouplingChannelDescriptor]
    private var previousLedgerDigest: String

    public init(
        participants: [any VivoCoupledParticipant],
        channels: [VivoCouplingChannelDescriptor],
        genesisDigest: String = String(repeating: "0", count: 64)
    ) throws {
        var participantMap: [VivoNumiLabParticipant: any VivoCoupledParticipant] = [:]
        for participant in participants {
            guard participantMap[participant.participantID] == nil else {
                throw VivoNumiLabCouplingError.duplicateParticipant(participant.participantID)
            }
            participantMap[participant.participantID] = participant
        }

        var channelMap: [String: VivoCouplingChannelDescriptor] = [:]
        for channel in channels {
            try channel.validate()
            guard channelMap[channel.id] == nil else {
                throw VivoNumiLabCouplingError.duplicateChannel(channel.id)
            }
            guard participantMap[channel.producer] != nil || channel.producer == .environment || channel.producer == .experiment else {
                throw VivoNumiLabCouplingError.missingParticipant(channel.producer)
            }
            guard participantMap[channel.consumer] != nil else {
                throw VivoNumiLabCouplingError.missingParticipant(channel.consumer)
            }
            channelMap[channel.id] = channel
        }

        self.participants = participantMap
        self.orderedParticipants = participantMap.keys.sorted { $0.rawValue < $1.rawValue }
        self.channels = channelMap
        self.previousLedgerDigest = genesisDigest
    }

    public func makeTransaction(
        clock: VivoCouplingClock,
        externalInputs: [VivoCouplingSample],
        iteration: UInt32 = 0
    ) throws -> VivoCouplingTransaction {
        try clock.validate()
        try validate(samples: externalInputs, permitInternalProducers: false)
        let id = Self.transactionDigest(clock: clock, inputs: externalInputs, iteration: iteration)
        return VivoCouplingTransaction(id: id, clock: clock, inputs: externalInputs, iteration: iteration)
    }

    public func step(_ transaction: VivoCouplingTransaction) async throws -> VivoJointStepCertificate {
        try transaction.clock.validate()
        guard transaction.id == Self.transactionDigest(
            clock: transaction.clock,
            inputs: transaction.inputs,
            iteration: transaction.iteration
        ) else {
            throw VivoNumiLabCouplingError.invalidTransaction("transaction fingerprint does not match its contents")
        }
        try validate(samples: transaction.inputs, permitInternalProducers: true)

        var prepared: [VivoPreparedCouplingState] = []
        prepared.reserveCapacity(orderedParticipants.count)

        do {
            for participantID in orderedParticipants {
                guard let participant = participants[participantID] else {
                    throw VivoNumiLabCouplingError.missingParticipant(participantID)
                }
                let participantInputs = transaction.inputs.filter {
                    channels[$0.channelID]?.consumer == participantID
                }
                let participantTransaction = VivoCouplingTransaction(
                    id: transaction.id,
                    clock: transaction.clock,
                    inputs: participantInputs,
                    iteration: transaction.iteration
                )
                let candidate = try await participant.prepare(participantTransaction)
                guard candidate.participant == participantID,
                      candidate.transactionID == transaction.id else {
                    throw VivoNumiLabCouplingError.participantMismatch
                }
                try validate(samples: candidate.outputs, permitInternalProducers: true)
                prepared.append(candidate)
            }
        } catch {
            await rollback(prepared, transactionID: transaction.id)
            throw error
        }

        let terminal = prepared.map(\.disposition).max(by: {
            Self.dispositionRank($0) < Self.dispositionRank($1)
        }) ?? .prepared
        let maximumAcceptedStep = prepared.compactMap(\.maximumAcceptedStep).min()
        let outputs = prepared.flatMap(\.outputs).sorted(by: Self.sampleOrder)

        guard terminal == .prepared else {
            await rollback(prepared, transactionID: transaction.id)
            let disposition: VivoJointStepDisposition
            switch terminal {
            case .prepared: disposition = .committed
            case .needsSmallerStep: disposition = .needsSmallerStep
            case .rejected: disposition = .rejected
            case .reversibleShutdown: disposition = .reversibleShutdown
            case .permanentShutdown: disposition = .permanentShutdown
            }
            return makeCertificate(
                transaction: transaction,
                disposition: disposition,
                maximumAcceptedStep: maximumAcceptedStep,
                prepared: prepared,
                outputs: outputs,
                advanceLedger: false
            )
        }

        var staged: [VivoPreparedCouplingState] = []
        do {
            for candidate in prepared {
                guard let participant = participants[candidate.participant] else {
                    throw VivoNumiLabCouplingError.missingParticipant(candidate.participant)
                }
                try await participant.stageCommit(transactionID: transaction.id)
                staged.append(candidate)
            }
        } catch {
            await rollback(staged, transactionID: transaction.id)
            let unstaged = prepared.filter { candidate in
                !staged.contains(where: { $0.participant == candidate.participant })
            }
            await rollback(unstaged, transactionID: transaction.id)
            throw VivoNumiLabCouplingError.stagingFailure(error.localizedDescription)
        }

        // `releaseCommit` is deliberately nonthrowing by protocol. Once every
        // participant has staged, publication cannot be partially declined.
        for candidate in prepared {
            await participants[candidate.participant]?.releaseCommit(transactionID: transaction.id)
        }

        return makeCertificate(
            transaction: transaction,
            disposition: .committed,
            maximumAcceptedStep: transaction.clock.deltaTime,
            prepared: prepared,
            outputs: outputs,
            advanceLedger: true
        )
    }

    private func validate(samples: [VivoCouplingSample], permitInternalProducers: Bool) throws {
        var seen = Set<String>()
        for sample in samples {
            guard let descriptor = channels[sample.channelID] else {
                throw VivoNumiLabCouplingError.undeclaredOutput(sample.channelID)
            }
            guard descriptor.producer == sample.producer,
                  descriptor.schemaVersion == sample.schemaVersion,
                  descriptor.unit == sample.unit else {
                throw VivoNumiLabCouplingError.invalidOutput("channel \(sample.channelID) metadata mismatch")
            }
            guard permitInternalProducers || sample.producer == .environment || sample.producer == .experiment else {
                throw VivoNumiLabCouplingError.invalidOutput("external input claims internal producer \(sample.producer.rawValue)")
            }
            guard let expectedCount = descriptor.elementCount, sample.values.count == expectedCount else {
                throw VivoNumiLabCouplingError.invalidOutput("channel \(sample.channelID) element count mismatch")
            }
            guard sample.sourceTime.isFinite,
                  sample.quality.isFinite,
                  sample.quality >= 0,
                  sample.quality <= 1,
                  sample.values.allSatisfy(\.isFinite) else {
                throw VivoNumiLabCouplingError.invalidOutput("channel \(sample.channelID) contains non-finite data")
            }
            for value in sample.values {
                if let lower = descriptor.lowerBound, value < lower {
                    throw VivoNumiLabCouplingError.invalidOutput("channel \(sample.channelID) violates its lower bound")
                }
                if let upper = descriptor.upperBound, value > upper {
                    throw VivoNumiLabCouplingError.invalidOutput("channel \(sample.channelID) violates its upper bound")
                }
            }
            let identity = "\(sample.producer.rawValue):\(sample.channelID)"
            guard seen.insert(identity).inserted else {
                throw VivoNumiLabCouplingError.invalidOutput("duplicate channel sample \(identity)")
            }
        }
    }

    private func rollback(_ candidates: [VivoPreparedCouplingState], transactionID: String) async {
        for candidate in candidates.reversed() {
            await participants[candidate.participant]?.rollback(transactionID: transactionID)
        }
    }

    private func makeCertificate(
        transaction: VivoCouplingTransaction,
        disposition: VivoJointStepDisposition,
        maximumAcceptedStep: Double?,
        prepared: [VivoPreparedCouplingState],
        outputs: [VivoCouplingSample],
        advanceLedger: Bool
    ) -> VivoJointStepCertificate {
        let digest = Self.ledgerDigest(
            previous: previousLedgerDigest,
            transaction: transaction,
            disposition: disposition,
            prepared: prepared,
            outputs: outputs
        )
        if advanceLedger {
            previousLedgerDigest = digest
        }
        return VivoJointStepCertificate(
            transactionID: transaction.id,
            disposition: disposition,
            clock: transaction.clock,
            requestedStep: transaction.clock.deltaTime,
            maximumAcceptedStep: maximumAcceptedStep,
            participants: prepared,
            routedOutputs: outputs,
            ledgerDigest: digest
        )
    }

    private static func transactionDigest(
        clock: VivoCouplingClock,
        inputs: [VivoCouplingSample],
        iteration: UInt32
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-coupling-transaction-v1".utf8))
        update(&hasher, clock.epoch)
        update(&hasher, clock.step)
        update(&hasher, clock.time.bitPattern)
        update(&hasher, clock.deltaTime.bitPattern)
        update(&hasher, iteration)
        for sample in inputs.sorted(by: sampleOrder) {
            hasher.update(data: Data(sample.channelID.utf8))
            hasher.update(data: Data([0]))
            hasher.update(data: Data(sample.producer.rawValue.utf8))
            hasher.update(data: Data([0]))
            update(&hasher, sample.schemaVersion)
            hasher.update(data: Data(sample.unit.utf8))
            hasher.update(data: Data([0]))
            update(&hasher, sample.sourceStep)
            update(&hasher, sample.sourceTime.bitPattern)
            update(&hasher, sample.quality.bitPattern)
            for value in sample.values { update(&hasher, value.bitPattern) }
        }
        return hasher.finalize().map { String(format: "%02x", $0) }.joined()
    }

    private static func ledgerDigest(
        previous: String,
        transaction: VivoCouplingTransaction,
        disposition: VivoJointStepDisposition,
        prepared: [VivoPreparedCouplingState],
        outputs: [VivoCouplingSample]
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-coupling-ledger-v1".utf8))
        hasher.update(data: Data(previous.utf8))
        hasher.update(data: Data(transaction.id.utf8))
        hasher.update(data: Data(disposition.rawValue.utf8))
        for candidate in prepared.sorted(by: { $0.participant.rawValue < $1.participant.rawValue }) {
            hasher.update(data: Data(candidate.participant.rawValue.utf8))
            hasher.update(data: Data(candidate.disposition.rawValue.utf8))
            hasher.update(data: Data(candidate.stateDigest.utf8))
        }
        for sample in outputs.sorted(by: sampleOrder) {
            hasher.update(data: Data(sample.channelID.utf8))
            for value in sample.values { update(&hasher, value.bitPattern) }
        }
        return hasher.finalize().map { String(format: "%02x", $0) }.joined()
    }

    private static func dispositionRank(_ disposition: VivoPreparationDisposition) -> Int {
        switch disposition {
        case .prepared: return 0
        case .needsSmallerStep: return 1
        case .rejected: return 2
        case .reversibleShutdown: return 3
        case .permanentShutdown: return 4
        }
    }

    private static func sampleOrder(_ lhs: VivoCouplingSample, _ rhs: VivoCouplingSample) -> Bool {
        if lhs.producer != rhs.producer { return lhs.producer.rawValue < rhs.producer.rawValue }
        if lhs.channelID != rhs.channelID { return lhs.channelID < rhs.channelID }
        if lhs.sourceStep != rhs.sourceStep { return lhs.sourceStep < rhs.sourceStep }
        return lhs.sourceTime < rhs.sourceTime
    }

    private static func update<T: FixedWidthInteger>(_ hasher: inout SHA256, _ value: T) {
        var littleEndian = value.littleEndian
        withUnsafeBytes(of: &littleEndian) { bytes in
            hasher.update(data: Data(bytes))
        }
    }
}
