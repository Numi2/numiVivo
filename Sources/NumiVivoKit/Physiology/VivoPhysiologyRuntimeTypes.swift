import Foundation

public struct VivoPhysiologyRuntimeConfiguration: Codable, Sendable, Equatable {
    public var maximumSubsteps: UInt32
    public var boundTolerance: Float
    public var maximumAbsoluteDerivative: Float
    public var maximumTransformsPerStep: Int
    public var maximumPublicationsPerStep: Int
    public var privateHeapHeadroom: Double

    public init(
        maximumSubsteps: UInt32 = 12,
        boundTolerance: Float = 1e-6,
        maximumAbsoluteDerivative: Float = .greatestFiniteMagnitude,
        maximumTransformsPerStep: Int = 262_144,
        maximumPublicationsPerStep: Int = 65_536,
        privateHeapHeadroom: Double = 1.15
    ) {
        self.maximumSubsteps = maximumSubsteps
        self.boundTolerance = boundTolerance
        self.maximumAbsoluteDerivative = maximumAbsoluteDerivative
        self.maximumTransformsPerStep = maximumTransformsPerStep
        self.maximumPublicationsPerStep = maximumPublicationsPerStep
        self.privateHeapHeadroom = privateHeapHeadroom
    }

    public func validate(for model: PreparedVivoPhysiologyModel) throws {
        guard maximumSubsteps > 0,
              boundTolerance.isFinite, boundTolerance >= 0,
              maximumAbsoluteDerivative.isFinite, maximumAbsoluteDerivative > 0,
              maximumTransformsPerStep > 0,
              maximumPublicationsPerStep > 0,
              privateHeapHeadroom.isFinite, privateHeapHeadroom >= 1,
              model.schema == PreparedVivoPhysiologyModel.schema,
              model.environmentCount > 0,
              model.pairCount > 0,
              model.incidenceOffsets.count == Int(model.pairCount) + 1,
              model.clearances.count == Int(model.pairCount),
              model.initialState.count == Int(model.pairCount) * Int(model.environmentCount),
              model.minimumTimeStepSeconds > 0,
              model.preferredTimeStepSeconds >= model.minimumTimeStepSeconds,
              model.maximumTimeStepSeconds >= model.preferredTimeStepSeconds else {
            throw VivoRuntimeError.invalidConfiguration("physiology runtime configuration or prepared model is invalid")
        }
    }
}

public struct VivoPhysiologyPublicationRequest: Codable, Sendable, Equatable, Hashable {
    public let pairIndex: UInt32
    public let environmentIndex: UInt32
    public let flags: UInt32

    public init(pairIndex: UInt32, environmentIndex: UInt32, flags: UInt32 = 0) {
        self.pairIndex = pairIndex
        self.environmentIndex = environmentIndex
        self.flags = flags
    }
}

public struct VivoPhysiologyPrepareRequest: Codable, Sendable, Equatable {
    public var timeStepSeconds: Double?
    public var preUpdates: [VivoPhysiologyStateUpdate]
    public var postUpdates: [VivoPhysiologyStateUpdate]
    public var publications: [VivoPhysiologyPublicationRequest]
    public var permitAdaptiveReduction: Bool

    public init(
        timeStepSeconds: Double? = nil,
        preUpdates: [VivoPhysiologyStateUpdate] = [],
        postUpdates: [VivoPhysiologyStateUpdate] = [],
        publications: [VivoPhysiologyPublicationRequest] = [],
        permitAdaptiveReduction: Bool = true
    ) {
        self.timeStepSeconds = timeStepSeconds
        self.preUpdates = preUpdates
        self.postUpdates = postUpdates
        self.publications = publications
        self.permitAdaptiveReduction = permitAdaptiveReduction
    }
}

public enum VivoPreparedPhysiologyDisposition: String, Codable, Sendable {
    case prepared
    case requiresSmallerStep
    case rejected
}

public struct VivoPreparedPhysiologyStep: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoPreparedPhysiologyDisposition
    public let modelFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let timeBefore: Double
    public let requestedTimeStep: Double
    public let candidateTimeStep: Double
    public let attemptCount: UInt32
    public let status: VivoPhysiologyRuntimeStatus
    public let publications: [Float]
    public let appliedDoseIdentifiers: [String]
    public let nextBoundarySeconds: Double?

    public var canCommit: Bool {
        disposition == .prepared && !status.blocksCommit
    }
}

public enum VivoPhysiologyStepDisposition: String, Codable, Sendable {
    case committed
    case committedWithReducedStep
    case rejected
}

public struct VivoPhysiologyStepCertificate: Codable, Sendable, Equatable {
    public let transactionID: UUID
    public let disposition: VivoPhysiologyStepDisposition
    public let modelFingerprint: VivoFingerprint
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let stepIndex: UInt32
    public let timeBefore: Double
    public let timeAfter: Double
    public let requestedTimeStep: Double
    public let acceptedTimeStep: Double?
    public let attemptCount: UInt32
    public let status: VivoPhysiologyRuntimeStatus
    public let appliedDoseIdentifiers: [String]

    public var committed: Bool {
        disposition == .committed || disposition == .committedWithReducedStep
    }
}

public struct VivoPhysiologyStepResult: Codable, Sendable, Equatable {
    public let certificate: VivoPhysiologyStepCertificate
    public let publications: [Float]
}

public struct VivoPhysiologySnapshot: Codable, Sendable, Equatable {
    public let modelFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let absoluteTimeSeconds: Double
    public let pairCount: UInt32
    public let environmentCount: UInt32
    public let values: [Float]

    public func value(pairIndex: UInt32, environmentIndex: UInt32) -> Float? {
        guard pairIndex < pairCount, environmentIndex < environmentCount else { return nil }
        return values[Int(pairIndex * environmentCount + environmentIndex)]
    }
}

public struct VivoPhysiologyCheckpoint: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/physiology-checkpoint/v1"

    public let schema: String
    public let modelFingerprint: VivoFingerprint
    public let stepIndex: UInt32
    public let absoluteTimeSeconds: Double
    public let doseCursor: Int
    public let pairCount: UInt32
    public let environmentCount: UInt32
    public let stateFP32LE: Data

    public init(
        modelFingerprint: VivoFingerprint,
        stepIndex: UInt32,
        absoluteTimeSeconds: Double,
        doseCursor: Int,
        pairCount: UInt32,
        environmentCount: UInt32,
        stateFP32LE: Data
    ) {
        self.schema = Self.schema
        self.modelFingerprint = modelFingerprint
        self.stepIndex = stepIndex
        self.absoluteTimeSeconds = absoluteTimeSeconds
        self.doseCursor = doseCursor
        self.pairCount = pairCount
        self.environmentCount = environmentCount
        self.stateFP32LE = stateFP32LE
    }

    public func validate() throws {
        let count = UInt64(pairCount) * UInt64(environmentCount)
        let bytes = count.multipliedReportingOverflow(by: 4)
        guard schema == Self.schema,
              pairCount > 0,
              environmentCount > 0,
              absoluteTimeSeconds.isFinite,
              absoluteTimeSeconds >= 0,
              doseCursor >= 0,
              !bytes.overflow,
              bytes.partialValue <= UInt64(Int.max),
              stateFP32LE.count == Int(bytes.partialValue) else {
            throw VivoArtifactValidationError.invalid("physiology checkpoint is invalid")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}
