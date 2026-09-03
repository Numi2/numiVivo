import Foundation

struct NVivoHashUniformsSwift: Sendable {
    var byteCountLow: UInt32
    var byteCountHigh: UInt32
    var chunkBytes: UInt32
    var chunkCount: UInt32
    var outputWordOffset: UInt32
    var domain: UInt32
    var logicalStepLow: UInt32
    var logicalStepHigh: UInt32
}

public struct VivoPreparedStep: Sendable, Hashable {
    public var receipt: VivoPrepareReceipt
    public var requestedDeltaTime: Float
    public var executedSubsteps: Int
    public var projectedLogicalStep: UInt64
    public var projectedTime: Double
    public var currentStateVersion: UInt32
    public var diagnostics: VivoStepDiagnostics
    public var events: [VivoEvent]

    public init(
        receipt: VivoPrepareReceipt,
        requestedDeltaTime: Float,
        executedSubsteps: Int,
        projectedLogicalStep: UInt64,
        projectedTime: Double,
        currentStateVersion: UInt32,
        diagnostics: VivoStepDiagnostics,
        events: [VivoEvent]
    ) {
        self.receipt = receipt
        self.requestedDeltaTime = requestedDeltaTime
        self.executedSubsteps = executedSubsteps
        self.projectedLogicalStep = projectedLogicalStep
        self.projectedTime = projectedTime
        self.currentStateVersion = currentStateVersion
        self.diagnostics = diagnostics
        self.events = events
    }
}

public struct VivoTransactionalRuntimeSnapshot: Sendable, Hashable, Codable {
    public var committedLogicalStep: UInt64
    public var committedTime: Double
    public var stateVersion: UInt32
    public var pendingEpochID: UUID?
    public var terminalStatus: VivoStepStatus?

    public init(
        committedLogicalStep: UInt64,
        committedTime: Double,
        stateVersion: UInt32,
        pendingEpochID: UUID?,
        terminalStatus: VivoStepStatus?
    ) {
        self.committedLogicalStep = committedLogicalStep
        self.committedTime = committedTime
        self.stateVersion = stateVersion
        self.pendingEpochID = pendingEpochID
        self.terminalStatus = terminalStatus
    }
}

public enum VivoTransactionalRuntimeError: Error, Sendable, CustomStringConvertible {
    case invalidEpoch(String)
    case programFingerprintMismatch(expected: String, actual: String)
    case sequenceMismatch(expected: UInt64, actual: UInt64)
    case timeMismatch(expected: Double, actual: Double)
    case alreadyPrepared(UUID)
    case noPreparedTransaction
    case epochMismatch(expected: UUID, actual: UUID)
    case receiptParticipantMismatch(expected: String, actual: String)
    case receiptTokenMismatch
    case receiptFingerprintMismatch
    case preparedStateChanged
    case invalidPreparedPublication(UInt32)

    public var description: String {
        switch self {
        case .invalidEpoch(let message):
            return "Invalid transaction epoch: \(message)"
        case .programFingerprintMismatch(let expected, let actual):
            return "Transaction program fingerprint '\(actual)' does not match loaded ProgramPack '\(expected)'."
        case .sequenceMismatch(let expected, let actual):
            return "Transaction sequence \(actual) does not match required sequence \(expected)."
        case .timeMismatch(let expected, let actual):
            return "Transaction simulationTime \(actual) does not match committed time \(expected)."
        case .alreadyPrepared(let id):
            return "Transaction \(id) is already prepared and must be committed or rolled back first."
        case .noPreparedTransaction:
            return "No NumiVivo transaction is currently prepared."
        case .epochMismatch(let expected, let actual):
            return "Prepared epoch \(expected) does not match supplied epoch \(actual)."
        case .receiptParticipantMismatch(let expected, let actual):
            return "Receipt participant '\(actual)' does not match runtime participant '\(expected)'."
        case .receiptTokenMismatch:
            return "Prepared transaction receipt token does not match."
        case .receiptFingerprintMismatch:
            return "Prepared transaction receipt fingerprint does not match."
        case .preparedStateChanged:
            return "Prepared GPU state no longer matches its attestation."
        case .invalidPreparedPublication(let status):
            return "GPU publication status \(status) is not a valid prepared state."
        }
    }
}
