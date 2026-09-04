import CryptoKit
import Foundation
import Metal

public struct VivoAdaptiveCohortRuntimeLayout: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let representation: VivoStateRepresentation
    public let authoritativeElementsPerChunk: UInt64
    public let stateCopyCount: UInt32
    public let randomCounterBytes: UInt64
    public let delayedQueueBytes: UInt64
    public let topologyBytes: UInt64
    public let auxiliaryBytes: UInt64
    public let integerMigrationTolerance: Float

    public init(
        cohortID: UInt32,
        representation: VivoStateRepresentation,
        authoritativeElementsPerChunk: UInt64,
        stateCopyCount: UInt32 = 4,
        randomCounterBytes: UInt64 = 0,
        delayedQueueBytes: UInt64 = 0,
        topologyBytes: UInt64 = 0,
        auxiliaryBytes: UInt64 = 0,
        integerMigrationTolerance: Float = 1e-4
    ) {
        self.cohortID = cohortID
        self.representation = representation
        self.authoritativeElementsPerChunk = authoritativeElementsPerChunk
        self.stateCopyCount = stateCopyCount
        self.randomCounterBytes = randomCounterBytes
        self.delayedQueueBytes = delayedQueueBytes
        self.topologyBytes = topologyBytes
        self.auxiliaryBytes = auxiliaryBytes
        self.integerMigrationTolerance = integerMigrationTolerance
    }

    public var bytesPerElement: UInt64 {
        switch representation {
        case .continuousFP32, .discreteUInt32: return 4
        case .mixedContinuousDiscrete: return 8
        }
    }

    public func stateCopyBytes() throws -> UInt64 {
        let result = authoritativeElementsPerChunk.multipliedReportingOverflow(by: bytesPerElement)
        guard !result.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
        return result.partialValue
    }

    public func requiredStateBytes() throws -> UInt64 {
        let copyBytes = try stateCopyBytes()
        let result = copyBytes.multipliedReportingOverflow(by: UInt64(stateCopyCount))
        guard !result.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
        return result.partialValue
    }

    public func requiredScratchBytes(alignment: UInt64 = 256) throws -> UInt64 {
        var cursor: UInt64 = 0
        for length in [randomCounterBytes, delayedQueueBytes, topologyBytes, auxiliaryBytes] {
            cursor = try Self.align(cursor, to: alignment, cohortID: cohortID)
            let next = cursor.addingReportingOverflow(length)
            guard !next.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
            cursor = next.partialValue
        }
        return cursor
    }

    public func validate(against plan: VivoFidelityCohortPlan) throws {
        guard cohortID == plan.cohortID,
              authoritativeElementsPerChunk > 0,
              stateCopyCount >= 2,
              integerMigrationTolerance.isFinite,
              integerMigrationTolerance >= 0 else {
            throw VivoAdaptiveArenaError.invalidLayout(cohortID, "identity, element count, copy count or tolerance")
        }
        let state = try requiredStateBytes()
        let scratch = try requiredScratchBytes()
        guard state <= plan.stateBytesPerChunk else {
            throw VivoAdaptiveArenaError.invalidLayout(
                cohortID,
                "layout requires \(state) state bytes but plan permits \(plan.stateBytesPerChunk)"
            )
        }
        guard scratch <= plan.scratchBytesPerChunk else {
            throw VivoAdaptiveArenaError.invalidLayout(
                cohortID,
                "layout requires \(scratch) scratch bytes but plan permits \(plan.scratchBytesPerChunk)"
            )
        }
    }

    private static func align(_ value: UInt64, to alignment: UInt64, cohortID: UInt32) throws -> UInt64 {
        guard alignment > 0, alignment.nonzeroBitCount == 1 else {
            throw VivoAdaptiveArenaError.invalidAlignment
        }
        let mask = alignment - 1
        let sum = value.addingReportingOverflow(mask)
        guard !sum.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID) }
        return sum.partialValue & ~mask
    }
}

public struct VivoAdaptiveScratchOffsets: Codable, Equatable, Sendable {
    public let randomCounterOffset: UInt64
    public let delayedQueueOffset: UInt64
    public let topologyOffset: UInt64
    public let auxiliaryOffset: UInt64

    public init(layout: VivoAdaptiveCohortRuntimeLayout, alignment: UInt64 = 256) throws {
        func align(_ value: UInt64) throws -> UInt64 {
            guard alignment > 0, alignment.nonzeroBitCount == 1 else {
                throw VivoAdaptiveArenaError.invalidAlignment
            }
            let sum = value.addingReportingOverflow(alignment - 1)
            guard !sum.overflow else { throw VivoAdaptiveArenaError.arithmeticOverflow(layout.cohortID) }
            return sum.partialValue & ~(alignment - 1)
        }
        randomCounterOffset = 0
        delayedQueueOffset = try align(randomCounterOffset + layout.randomCounterBytes)
        topologyOffset = try align(delayedQueueOffset + layout.delayedQueueBytes)
        auxiliaryOffset = try align(topologyOffset + layout.topologyBytes)
    }
}

public struct VivoAdaptiveChunkDescriptor: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let chunkIndex: UInt32
    public let stateBytes: UInt64
    public let scratchBytes: UInt64
    public let stateCopyBytes: UInt64
    public let stateCopyCount: UInt32
    public let representation: VivoStateRepresentation
    public let scratchOffsets: VivoAdaptiveScratchOffsets

    public init(
        cohortID: UInt32,
        chunkIndex: UInt32,
        stateBytes: UInt64,
        scratchBytes: UInt64,
        stateCopyBytes: UInt64,
        stateCopyCount: UInt32,
        representation: VivoStateRepresentation,
        scratchOffsets: VivoAdaptiveScratchOffsets
    ) {
        self.cohortID = cohortID
        self.chunkIndex = chunkIndex
        self.stateBytes = stateBytes
        self.scratchBytes = scratchBytes
        self.stateCopyBytes = stateCopyBytes
        self.stateCopyCount = stateCopyCount
        self.representation = representation
        self.scratchOffsets = scratchOffsets
    }

    public func stateCopyOffset(_ index: UInt32) throws -> Int {
        guard index < stateCopyCount else {
            throw VivoAdaptiveArenaError.invalidStateCopy(cohortID, index)
        }
        let product = UInt64(index).multipliedReportingOverflow(by: stateCopyBytes)
        guard !product.overflow, product.partialValue <= UInt64(Int.max) else {
            throw VivoAdaptiveArenaError.arithmeticOverflow(cohortID)
        }
        return Int(product.partialValue)
    }
}

public struct VivoAdaptiveArenaChunkResources: @unchecked Sendable {
    public let descriptor: VivoAdaptiveChunkDescriptor
    public let state: MTLBuffer
    public let scratch: MTLBuffer

    public init(descriptor: VivoAdaptiveChunkDescriptor, state: MTLBuffer, scratch: MTLBuffer) {
        self.descriptor = descriptor
        self.state = state
        self.scratch = scratch
    }
}

public struct VivoAdaptiveArenaManifest: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let planFingerprint: String
    public let deviceName: String
    public let registryID: UInt64
    public let heapBytes: UInt64
    public let allocatedStateBytes: UInt64
    public let allocatedScratchBytes: UInt64
    public let chunkCount: UInt32
    public let chunks: [VivoAdaptiveChunkDescriptor]
    public let fingerprint: String

    public init(
        planFingerprint: String,
        deviceName: String,
        registryID: UInt64,
        heapBytes: UInt64,
        allocatedStateBytes: UInt64,
        allocatedScratchBytes: UInt64,
        chunkCount: UInt32,
        chunks: [VivoAdaptiveChunkDescriptor]
    ) throws {
        let unsigned = VivoAdaptiveArenaManifest(
            schemaVersion: 1,
            planFingerprint: planFingerprint,
            deviceName: deviceName,
            registryID: registryID,
            heapBytes: heapBytes,
            allocatedStateBytes: allocatedStateBytes,
            allocatedScratchBytes: allocatedScratchBytes,
            chunkCount: chunkCount,
            chunks: chunks,
            fingerprint: ""
        )
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        let digest = SHA256.hash(data: try encoder.encode(unsigned))
        self = .init(
            schemaVersion: unsigned.schemaVersion,
            planFingerprint: unsigned.planFingerprint,
            deviceName: unsigned.deviceName,
            registryID: unsigned.registryID,
            heapBytes: unsigned.heapBytes,
            allocatedStateBytes: unsigned.allocatedStateBytes,
            allocatedScratchBytes: unsigned.allocatedScratchBytes,
            chunkCount: unsigned.chunkCount,
            chunks: unsigned.chunks,
            fingerprint: digest.map { String(format: "%02x", $0) }.joined()
        )
    }

    private init(
        schemaVersion: UInt32,
        planFingerprint: String,
        deviceName: String,
        registryID: UInt64,
        heapBytes: UInt64,
        allocatedStateBytes: UInt64,
        allocatedScratchBytes: UInt64,
        chunkCount: UInt32,
        chunks: [VivoAdaptiveChunkDescriptor],
        fingerprint: String
    ) {
        self.schemaVersion = schemaVersion
        self.planFingerprint = planFingerprint
        self.deviceName = deviceName
        self.registryID = registryID
        self.heapBytes = heapBytes
        self.allocatedStateBytes = allocatedStateBytes
        self.allocatedScratchBytes = allocatedScratchBytes
        self.chunkCount = chunkCount
        self.chunks = chunks
        self.fingerprint = fingerprint
    }
}

public struct VivoAdaptiveMigrationStatusFlags: OptionSet, Codable, Sendable {
    public let rawValue: UInt32
    public init(rawValue: UInt32) { self.rawValue = rawValue }

    public static let invalidNumber = Self(rawValue: 1 << 0)
    public static let negativeCount = Self(rawValue: 1 << 1)
    public static let countOverflow = Self(rawValue: 1 << 2)
    public static let fractionalCount = Self(rawValue: 1 << 3)
    public static let floatPrecisionLoss = Self(rawValue: 1 << 4)
}

public struct VivoAdaptiveMigrationStatus: Codable, Equatable, Sendable {
    public let flags: VivoAdaptiveMigrationStatusFlags
    public let invalidElementCount: UInt32
    public let firstInvalidElement: UInt32?
    public let maximumError: Float
}

public struct VivoAdaptiveArenaMigrationRecord: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let chunkIndex: UInt32
    public let kind: VivoBufferMigrationKind
    public let sourceRepresentation: VivoStateRepresentation
    public let destinationRepresentation: VivoStateRepresentation
    public let copiedStateBytes: UInt64
    public let copiedRandomCounterBytes: UInt64
    public let copiedDelayedQueueBytes: UInt64
    public let rebuiltTopology: Bool
}

public struct VivoAdaptiveArenaMigrationCertificate: Codable, Equatable, Sendable {
    public let sourceManifestFingerprint: String
    public let destinationManifestFingerprint: String
    public let committed: Bool
    public let records: [VivoAdaptiveArenaMigrationRecord]
    public let status: VivoAdaptiveMigrationStatus
}

public enum VivoAdaptiveArenaError: Error, LocalizedError, Sendable {
    case invalidAlignment
    case invalidPlan(String)
    case missingLayout(UInt32)
    case duplicateLayout(UInt32)
    case invalidLayout(UInt32, String)
    case arithmeticOverflow(UInt32)
    case workingSetExceeded(required: UInt64, available: UInt64)
    case heapCreationFailed
    case allocationFailed(UInt32, UInt32, String)
    case missingChunk(UInt32, UInt32)
    case invalidStateCopy(UInt32, UInt32)
    case invalidUpload(UInt32, UInt32, String)
    case migrationUnsupported(UInt32, VivoStateRepresentation, VivoStateRepresentation)
    case migrationRejected(VivoAdaptiveMigrationStatus)
    case shaderResourceUnavailable
    case shaderCompilation(String)
    case pipelineUnavailable(String)
    case commandQueueUnavailable
    case commandFailed(String)

    public var errorDescription: String? {
        switch self {
        case .invalidAlignment: return "adaptive arena alignment is invalid"
        case .invalidPlan(let reason): return "invalid adaptive arena plan: \(reason)"
        case .missingLayout(let cohort): return "missing runtime layout for cohort \(cohort)"
        case .duplicateLayout(let cohort): return "duplicate runtime layout for cohort \(cohort)"
        case .invalidLayout(let cohort, let reason): return "invalid runtime layout for cohort \(cohort): \(reason)"
        case .arithmeticOverflow(let cohort): return "adaptive arena arithmetic overflow for cohort \(cohort)"
        case .workingSetExceeded(let required, let available): return "adaptive arena requires \(required) bytes but budget is \(available)"
        case .heapCreationFailed: return "Metal heap creation failed"
        case .allocationFailed(let cohort, let chunk, let resource): return "failed to allocate \(resource) for cohort \(cohort), chunk \(chunk)"
        case .missingChunk(let cohort, let chunk): return "missing adaptive arena chunk \(cohort):\(chunk)"
        case .invalidStateCopy(let cohort, let copy): return "invalid state copy \(copy) for cohort \(cohort)"
        case .invalidUpload(let cohort, let chunk, let reason): return "invalid upload for cohort \(cohort), chunk \(chunk): \(reason)"
        case .migrationUnsupported(let cohort, let source, let destination): return "unsupported migration for cohort \(cohort): \(source.rawValue) to \(destination.rawValue)"
        case .migrationRejected(let status): return "adaptive arena migration rejected with flags \(status.flags.rawValue)"
        case .shaderResourceUnavailable: return "NumiVivo migration Metal source is unavailable"
        case .shaderCompilation(let reason): return "adaptive arena shader compilation failed: \(reason)"
        case .pipelineUnavailable(let name): return "adaptive arena pipeline unavailable: \(name)"
        case .commandQueueUnavailable: return "adaptive arena Metal command queue is unavailable"
        case .commandFailed(let reason): return "adaptive arena Metal command failed: \(reason)"
        }
    }
}

public extension VivoAdaptiveFidelityPlan {
    /// Actual chunked state and scratch residency. This intentionally does not
    /// rely on `totalResidentBytes`, which may also include compiler-owned fixed
    /// tables and is not a substitute for multiplying per-chunk resources.
    var actualChunkedArenaBytes: UInt64? {
        var total: UInt64 = 0
        for cohort in cohorts {
            let perChunk = cohort.stateBytesPerChunk.addingReportingOverflow(cohort.scratchBytesPerChunk)
            guard !perChunk.overflow else { return nil }
            let allChunks = perChunk.partialValue.multipliedReportingOverflow(by: UInt64(cohort.chunkCount))
            guard !allChunks.overflow else { return nil }
            let next = total.addingReportingOverflow(allChunks.partialValue)
            guard !next.overflow else { return nil }
            total = next.partialValue
        }
        return total
    }
}
