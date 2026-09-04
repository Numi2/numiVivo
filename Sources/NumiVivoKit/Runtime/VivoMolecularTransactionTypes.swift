import Foundation

public enum VivoPreparedMolecularDisposition: String, Codable, Sendable {
    case prepared
    case requiresSmallerStep
    case rejected
    case reversibleShutdown
    case permanentShutdown
}

public struct VivoPreparedMolecularStep: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoPreparedMolecularDisposition
    public let sourceFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let fidelity: VivoFidelity
    public let stepIndex: UInt32
    public let timeBefore: Float
    public let requestedTimeStep: Float
    public let candidateTimeStep: Float
    public let attemptCount: UInt32
    public let status: VivoRuntimeStatus
    public let events: [VivoEvent]
    public let publications: [Float]

    public var canCommit: Bool {
        disposition == .prepared && !status.blocksCommit
    }
}

public struct VivoMolecularTransactionCertificate: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoStepDisposition
    public let sourceFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let fidelity: VivoFidelity
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let stepIndex: UInt32
    public let timeBefore: Float
    public let timeAfter: Float
    public let requestedTimeStep: Float
    public let acceptedTimeStep: Float?
    public let attemptCount: UInt32
    public let status: VivoRuntimeStatus

    public var committed: Bool {
        disposition == .committed || disposition == .committedWithReducedStep
    }
}

public struct VivoMolecularTransactionResult: Codable, Sendable, Equatable {
    public let certificate: VivoMolecularTransactionCertificate
    public let events: [VivoEvent]
    public let publications: [Float]
}
