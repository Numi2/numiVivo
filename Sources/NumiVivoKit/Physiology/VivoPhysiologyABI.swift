import Foundation

public enum VivoPhysiologyUpdateMode: UInt32, Codable, Sendable, CaseIterable {
    case replace = 0
    case add = 1
    case rate = 2
    case minimum = 3
    case maximum = 4
}

public struct VivoPhysiologyStateUpdate: Codable, Sendable, Equatable, Hashable {
    public let pairIndex: UInt32
    public let environmentIndex: UInt32
    public let mode: VivoPhysiologyUpdateMode
    public let value: Float

    public init(
        pairIndex: UInt32,
        environmentIndex: UInt32,
        mode: VivoPhysiologyUpdateMode,
        value: Float
    ) {
        self.pairIndex = pairIndex
        self.environmentIndex = environmentIndex
        self.mode = mode
        self.value = value
    }
}

public struct VivoPhysiologyStateTransformABI: Sendable, Equatable {
    public static let replaceFlag: UInt32 = 1 << 0
    public static let minimumFlag: UInt32 = 1 << 1
    public static let maximumFlag: UInt32 = 1 << 2

    public var pairIndex: UInt32
    public var environmentIndex: UInt32
    public var flags: UInt32
    public var reserved: UInt32
    public var replacement: Float
    public var additiveDelta: Float
    public var minimum: Float
    public var maximum: Float

    public init(
        pairIndex: UInt32,
        environmentIndex: UInt32,
        flags: UInt32,
        replacement: Float,
        additiveDelta: Float,
        minimum: Float,
        maximum: Float
    ) {
        self.pairIndex = pairIndex
        self.environmentIndex = environmentIndex
        self.flags = flags
        self.reserved = 0
        self.replacement = replacement
        self.additiveDelta = additiveDelta
        self.minimum = minimum
        self.maximum = maximum
    }
}

public struct VivoPhysiologyPublicationRequestABI: Sendable, Equatable {
    public var pairIndex: UInt32
    public var environmentIndex: UInt32
    public var outputIndex: UInt32
    public var flags: UInt32

    public init(pairIndex: UInt32, environmentIndex: UInt32, outputIndex: UInt32, flags: UInt32 = 0) {
        self.pairIndex = pairIndex
        self.environmentIndex = environmentIndex
        self.outputIndex = outputIndex
        self.flags = flags
    }
}

public struct VivoPhysiologyRuntimeCommandABI: Sendable, Equatable {
    public var pairCount: UInt32
    public var environmentCount: UInt32
    public var stateElementCount: UInt32
    public var incidenceCount: UInt32

    public var preTransformCount: UInt32
    public var postTransformCount: UInt32
    public var publicationCount: UInt32
    public var stepIndex: UInt32

    public var appliedDoseCount: UInt32
    public var runtimeFlags: UInt32
    public var transactionWord0: UInt32
    public var transactionWord1: UInt32

    public var transactionWord2: UInt32
    public var transactionWord3: UInt32
    public var reservedWord0: UInt32
    public var reservedWord1: UInt32

    public var dt: Float
    public var absoluteTime: Float
    public var minimumTimeStep: Float
    public var maximumDerivative: Float

    public init(
        model: PreparedVivoPhysiologyModel,
        stepIndex: UInt32,
        dt: Float,
        absoluteTime: Float,
        transactionID: UUID,
        preTransformCount: Int,
        postTransformCount: Int,
        publicationCount: Int,
        appliedDoseCount: Int,
        maximumDerivative: Float
    ) throws {
        let stateCount = UInt64(model.pairCount) * UInt64(model.environmentCount)
        guard stateCount <= UInt64(UInt32.max),
              model.incidence.count <= Int(UInt32.max),
              preTransformCount >= 0, preTransformCount <= Int(UInt32.max),
              postTransformCount >= 0, postTransformCount <= Int(UInt32.max),
              publicationCount >= 0, publicationCount <= Int(UInt32.max),
              appliedDoseCount >= 0, appliedDoseCount <= Int(UInt32.max),
              dt.isFinite, dt > 0,
              absoluteTime.isFinite, absoluteTime >= 0,
              maximumDerivative.isFinite, maximumDerivative >= 0 else {
            throw VivoArtifactValidationError.invalid("physiology runtime command exceeds ABI or finite-value limits")
        }
        let words = Self.uuidWords(transactionID)
        self.pairCount = model.pairCount
        self.environmentCount = model.environmentCount
        self.stateElementCount = UInt32(stateCount)
        self.incidenceCount = UInt32(model.incidence.count)
        self.preTransformCount = UInt32(preTransformCount)
        self.postTransformCount = UInt32(postTransformCount)
        self.publicationCount = UInt32(publicationCount)
        self.stepIndex = stepIndex
        self.appliedDoseCount = UInt32(appliedDoseCount)
        self.runtimeFlags = 0
        self.transactionWord0 = words.0
        self.transactionWord1 = words.1
        self.transactionWord2 = words.2
        self.transactionWord3 = words.3
        self.reservedWord0 = 0
        self.reservedWord1 = 0
        self.dt = dt
        self.absoluteTime = absoluteTime
        self.minimumTimeStep = Float(model.minimumTimeStepSeconds)
        self.maximumDerivative = maximumDerivative
    }

    public static func validateMemoryLayout() throws {
        guard MemoryLayout<Self>.stride == 80,
              MemoryLayout<Self>.alignment == MemoryLayout<UInt32>.alignment,
              MemoryLayout<VivoPhysiologyStateTransformABI>.stride == 32,
              MemoryLayout<VivoPhysiologyPublicationRequestABI>.stride == 16,
              MemoryLayout<VivoPhysiologyIncidenceABI>.stride == 16,
              MemoryLayout<VivoPhysiologyClearanceABI>.stride == 16 else {
            throw VivoArtifactValidationError.incompatible("Swift physiology records do not match the Metal ABI")
        }
    }

    private static func uuidWords(_ uuid: UUID) -> (UInt32, UInt32, UInt32, UInt32) {
        var bytes = uuid.uuid
        return withUnsafeBytes(of: &bytes) { raw in
            let values = raw.bindMemory(to: UInt32.self)
            return (
                UInt32(littleEndian: values[0]),
                UInt32(littleEndian: values[1]),
                UInt32(littleEndian: values[2]),
                UInt32(littleEndian: values[3])
            )
        }
    }
}

public struct VivoPhysiologyRuntimeStatusABI: Sendable, Equatable {
    public var flags: UInt32
    public var violationCount: UInt32
    public var firstViolationPair: UInt32
    public var firstViolationEnvironment: UInt32

    public var maximumViolation: Float
    public var maximumAbsoluteDerivative: Float
    public var suggestedTimeStep: Float
    public var requestedResponse: UInt32

    public var appliedDoseCount: UInt32
    public var appliedPreTransformCount: UInt32
    public var appliedPostTransformCount: UInt32
    public var publicationCount: UInt32

    public var stepIndex: UInt32
    public var transactionWord0: UInt32
    public var transactionWord1: UInt32
    public var reserved: UInt32

    public init() {
        self.flags = 0
        self.violationCount = 0
        self.firstViolationPair = UInt32.max
        self.firstViolationEnvironment = UInt32.max
        self.maximumViolation = 0
        self.maximumAbsoluteDerivative = 0
        self.suggestedTimeStep = 0
        self.requestedResponse = 0
        self.appliedDoseCount = 0
        self.appliedPreTransformCount = 0
        self.appliedPostTransformCount = 0
        self.publicationCount = 0
        self.stepIndex = 0
        self.transactionWord0 = 0
        self.transactionWord1 = 0
        self.reserved = 0
    }

    public static func validateMemoryLayout() throws {
        guard MemoryLayout<Self>.stride == 64 else {
            throw VivoArtifactValidationError.incompatible("physiology status record does not match the Metal ABI")
        }
    }
}

public struct VivoPhysiologyRuntimeFlags: OptionSet, Codable, Sendable, Equatable, Hashable {
    public let rawValue: UInt32

    public init(rawValue: UInt32) {
        self.rawValue = rawValue
    }

    public static let nonFiniteState = Self(rawValue: 1 << 0)
    public static let belowMinimum = Self(rawValue: 1 << 1)
    public static let aboveMaximum = Self(rawValue: 1 << 2)
    public static let excessiveDerivative = Self(rawValue: 1 << 3)
    public static let invalidTransform = Self(rawValue: 1 << 4)
    public static let invalidIncidence = Self(rawValue: 1 << 5)
    public static let requiresSubstep = Self(rawValue: 1 << 6)
    public static let rejected = Self(rawValue: 1 << 7)
}

public struct VivoPhysiologyRuntimeStatus: Codable, Sendable, Equatable {
    public let flags: VivoPhysiologyRuntimeFlags
    public let violationCount: UInt32
    public let firstViolationPair: UInt32?
    public let firstViolationEnvironment: UInt32?
    public let maximumViolation: Float
    public let maximumAbsoluteDerivative: Float
    public let suggestedTimeStep: Float
    public let requestedResponse: UInt32
    public let appliedDoseCount: UInt32
    public let appliedPreTransformCount: UInt32
    public let appliedPostTransformCount: UInt32
    public let publicationCount: UInt32
    public let stepIndex: UInt32

    public init(raw: VivoPhysiologyRuntimeStatusABI) {
        self.flags = VivoPhysiologyRuntimeFlags(rawValue: raw.flags)
        self.violationCount = raw.violationCount
        self.firstViolationPair = raw.firstViolationPair == UInt32.max ? nil : raw.firstViolationPair
        self.firstViolationEnvironment = raw.firstViolationEnvironment == UInt32.max ? nil : raw.firstViolationEnvironment
        self.maximumViolation = raw.maximumViolation
        self.maximumAbsoluteDerivative = raw.maximumAbsoluteDerivative
        self.suggestedTimeStep = raw.suggestedTimeStep
        self.requestedResponse = raw.requestedResponse
        self.appliedDoseCount = raw.appliedDoseCount
        self.appliedPreTransformCount = raw.appliedPreTransformCount
        self.appliedPostTransformCount = raw.appliedPostTransformCount
        self.publicationCount = raw.publicationCount
        self.stepIndex = raw.stepIndex
    }

    public var blocksCommit: Bool {
        !flags.intersection([
            .nonFiniteState,
            .belowMinimum,
            .aboveMaximum,
            .invalidTransform,
            .invalidIncidence,
            .rejected
        ]).isEmpty || violationCount > 0
    }
}
