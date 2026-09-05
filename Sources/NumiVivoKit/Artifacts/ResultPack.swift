import Foundation

public struct VivoMeasurementSample: Codable, Sendable, Equatable {
    public let measurementIdentifier: String
    public let replicateIndex: UInt32
    public let stepIndex: UInt32
    public let timeSeconds: Double
    public let values: [Double]
    public let unit: String
    public let laneIndices: [UInt32]?

    public init(
        measurementIdentifier: String,
        replicateIndex: UInt32,
        stepIndex: UInt32,
        timeSeconds: Double,
        values: [Double],
        unit: String,
        laneIndices: [UInt32]? = nil
    ) {
        self.measurementIdentifier = measurementIdentifier
        self.replicateIndex = replicateIndex
        self.stepIndex = stepIndex
        self.timeSeconds = timeSeconds
        self.values = values
        self.unit = unit
        self.laneIndices = laneIndices
    }

    public func validate() throws {
        guard !measurementIdentifier.isEmpty,
              timeSeconds.isFinite,
              timeSeconds >= 0,
              !unit.isEmpty,
              !values.isEmpty,
              values.allSatisfy(\.isFinite) else {
            throw VivoArtifactValidationError.invalid("measurement sample contains missing or non-finite data")
        }
        if let laneIndices {
            guard laneIndices.count == values.count,
                  Set(laneIndices).count == laneIndices.count else {
                throw VivoArtifactValidationError.invalid("measurement lane indices must be unique and match values")
            }
        }
    }
}

public struct VivoStepLedgerEntry: Codable, Sendable, Equatable {
    public let replicateIndex: UInt32
    public let stepIndex: UInt32
    public let disposition: VivoStepDisposition
    public let attemptedStep: Float
    public let acceptedStep: Float?
    public let attemptCount: UInt32
    public let timeBefore: Float
    public let timeAfter: Float
    public let runtimeFlags: UInt32
    public let violationCount: UInt32
    public let maximumViolation: Float
    public let maximumRate: Float
    public let eventCount: UInt32
    public let eventDropped: UInt32

    public init(replicateIndex: UInt32, certificate: VivoStepCertificate) {
        self.replicateIndex = replicateIndex
        self.stepIndex = certificate.stepIndex
        self.disposition = certificate.disposition
        self.attemptedStep = certificate.attemptedStep
        self.acceptedStep = certificate.acceptedStep
        self.attemptCount = certificate.attemptCount
        self.timeBefore = certificate.timeBefore
        self.timeAfter = certificate.timeAfter
        self.runtimeFlags = certificate.status.flags.rawValue
        self.violationCount = certificate.status.violationCount
        self.maximumViolation = certificate.status.maximumViolation
        self.maximumRate = certificate.status.maximumRate
        self.eventCount = certificate.status.eventCount
        self.eventDropped = certificate.status.eventDropped
    }
}

public struct VivoRunLedger: Codable, Sendable, Equatable {
    public let genesis: VivoFingerprint
    public private(set) var head: VivoFingerprint
    public private(set) var entryCount: UInt64
    public private(set) var entries: [VivoStepLedgerEntry]

    public init(genesis: VivoFingerprint) {
        self.genesis = genesis
        self.head = genesis
        self.entryCount = 0
        self.entries = []
    }

    @discardableResult
    public mutating func append(_ entry: VivoStepLedgerEntry) throws -> VivoFingerprint {
        let material = VivoLedgerHashMaterial(
            previous: head,
            index: entryCount,
            entry: entry
        )
        head = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(material))
        entries.append(entry)
        entryCount &+= 1
        return head
    }

    public func verify() throws -> Bool {
        var cursor = genesis
        for (index, entry) in entries.enumerated() {
            cursor = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(
                VivoLedgerHashMaterial(previous: cursor, index: UInt64(index), entry: entry)
            ))
        }
        return cursor == head && entryCount == UInt64(entries.count)
    }
}

private struct VivoLedgerHashMaterial: Codable {
    let previous: VivoFingerprint
    let index: UInt64
    let entry: VivoStepLedgerEntry
}

public struct VivoReplicateResult: Codable, Sendable, Equatable {
    public enum Termination: String, Codable, Sendable {
        case completedDuration
        case stopCondition
        case rejectedStep
        case reversibleShutdown
        case permanentShutdown
        case runtimeFailure
    }

    public let replicateIndex: UInt32
    public let seed: VivoRuntimeSeed
    public let termination: Termination
    public let terminationReason: String
    public let finalStepIndex: UInt32
    public let finalTimeSeconds: Double
    public let ledgerHead: VivoFingerprint
    public let committedSteps: UInt64
    public let rejectedSteps: UInt64
    public let measurementCount: UInt64
    public let checkpointFingerprints: [VivoFingerprint]

    public init(
        replicateIndex: UInt32,
        seed: VivoRuntimeSeed,
        termination: Termination,
        terminationReason: String,
        finalStepIndex: UInt32,
        finalTimeSeconds: Double,
        ledgerHead: VivoFingerprint,
        committedSteps: UInt64,
        rejectedSteps: UInt64,
        measurementCount: UInt64,
        checkpointFingerprints: [VivoFingerprint]
    ) {
        self.replicateIndex = replicateIndex
        self.seed = seed
        self.termination = termination
        self.terminationReason = terminationReason
        self.finalStepIndex = finalStepIndex
        self.finalTimeSeconds = finalTimeSeconds
        self.ledgerHead = ledgerHead
        self.committedSteps = committedSteps
        self.rejectedSteps = rejectedSteps
        self.measurementCount = measurementCount
        self.checkpointFingerprints = checkpointFingerprints
    }
}

public struct VivoResultPack: Codable, Sendable, Equatable {
    public let experimentFingerprint: VivoFingerprint
    public let programFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint
    public let couplingFingerprint: VivoFingerprint?
    public let sourceProgramFingerprint: VivoFingerprint
    public let fidelity: VivoFidelity
    public let runtimeVersion: String
    public let deviceName: String
    public let deviceRegistryID: UInt64
    public let startedAt: Date
    public let finishedAt: Date
    public let replicates: [VivoReplicateResult]
    public let measurements: [VivoMeasurementSample]
    public let events: [VivoEvent]
    /// Optional for decoding legacy result packs; new experiment runs always supply these.
    public let recordedEvents: [VivoRecordedEvent]?
    public let measurementSummaries: [VivoMeasurementSummary]?
    public let ledgers: [VivoRunLedger]
    public let limitations: [String]
    public let annotations: [String: String]

    public init(
        experimentFingerprint: VivoFingerprint,
        programFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint,
        couplingFingerprint: VivoFingerprint?,
        sourceProgramFingerprint: VivoFingerprint,
        fidelity: VivoFidelity,
        runtimeVersion: String,
        deviceName: String,
        deviceRegistryID: UInt64,
        startedAt: Date,
        finishedAt: Date,
        replicates: [VivoReplicateResult],
        measurements: [VivoMeasurementSample],
        measurementSummaries: [VivoMeasurementSummary]? = nil,
        events: [VivoEvent],
        recordedEvents: [VivoRecordedEvent]? = nil,
        ledgers: [VivoRunLedger],
        limitations: [String],
        annotations: [String: String] = [:]
    ) {
        self.experimentFingerprint = experimentFingerprint
        self.programFingerprint = programFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.sourceProgramFingerprint = sourceProgramFingerprint
        self.fidelity = fidelity
        self.runtimeVersion = runtimeVersion
        self.deviceName = deviceName
        self.deviceRegistryID = deviceRegistryID
        self.startedAt = startedAt
        self.finishedAt = finishedAt
        self.replicates = replicates
        self.measurements = measurements
        self.events = events
        self.recordedEvents = recordedEvents
        self.measurementSummaries = measurementSummaries
        self.ledgers = ledgers
        self.limitations = limitations
        self.annotations = annotations
    }

    public func validate() throws {
        guard finishedAt >= startedAt,
              !runtimeVersion.isEmpty,
              !deviceName.isEmpty,
              replicates.count == ledgers.count else {
            throw VivoArtifactValidationError.invalid("result metadata or replicate ledger count is invalid")
        }
        guard Set(replicates.map(\.replicateIndex)).count == replicates.count else {
            throw VivoArtifactValidationError.invalid("result contains duplicate replicate indices")
        }
        let indices = Set(replicates.map(\.replicateIndex))
        if let records = recordedEvents {
            guard records.map(\.event) == events, records.allSatisfy({ indices.contains($0.replicateIndex) }) else {
                throw VivoArtifactValidationError.invalid("replicate-tagged events differ from raw events or reference unknown replicates")
            }
        }
        for summary in measurementSummaries ?? [] {
            try summary.validate()
            guard indices.contains(summary.replicateIndex) else { throw VivoArtifactValidationError.invalid("summary references an unknown replicate") }
        }
        for sample in measurements {
            try sample.validate()
            guard indices.contains(sample.replicateIndex) else { throw VivoArtifactValidationError.invalid("sample references an unknown replicate") }
        }
        for ledger in ledgers where try !ledger.verify() {
            throw VivoArtifactValidationError.invalid("result contains an invalid step ledger chain")
        }
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }
}

public struct VivoCheckpointPack: Codable, Sendable, Equatable {
    public static let formatVersion: UInt32 = 1

    public let version: UInt32
    public let programFingerprint: VivoFingerprint
    public let sourceProgramFingerprint: VivoFingerprint
    public let hostContextFingerprint: VivoFingerprint
    public let experimentFingerprint: VivoFingerprint?
    public let couplingFingerprint: VivoFingerprint?
    public let fidelity: VivoFidelity
    public let seed: VivoRuntimeSeed
    public let stepIndex: UInt32
    public let absoluteTimeSeconds: Double
    public let speciesCount: UInt32
    public let laneCount: UInt32
    public let parameterCount: UInt32
    public let parameterStride: UInt32
    public let temporalStateCount: UInt32
    public let stateFP32LE: Data
    public let parametersFP32LE: Data
    public let temporalStateFP32LE: Data
    public let transportFP32LE: Data
    public let ledgerHead: VivoFingerprint

    public init(
        programFingerprint: VivoFingerprint,
        sourceProgramFingerprint: VivoFingerprint,
        hostContextFingerprint: VivoFingerprint,
        experimentFingerprint: VivoFingerprint?,
        couplingFingerprint: VivoFingerprint?,
        fidelity: VivoFidelity,
        seed: VivoRuntimeSeed,
        stepIndex: UInt32,
        absoluteTimeSeconds: Double,
        speciesCount: UInt32,
        laneCount: UInt32,
        parameterCount: UInt32,
        parameterStride: UInt32,
        temporalStateCount: UInt32,
        stateFP32LE: Data,
        parametersFP32LE: Data,
        temporalStateFP32LE: Data,
        transportFP32LE: Data,
        ledgerHead: VivoFingerprint
    ) {
        self.version = Self.formatVersion
        self.programFingerprint = programFingerprint
        self.sourceProgramFingerprint = sourceProgramFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentFingerprint = experimentFingerprint
        self.couplingFingerprint = couplingFingerprint
        self.fidelity = fidelity
        self.seed = seed
        self.stepIndex = stepIndex
        self.absoluteTimeSeconds = absoluteTimeSeconds
        self.speciesCount = speciesCount
        self.laneCount = laneCount
        self.parameterCount = parameterCount
        self.parameterStride = parameterStride
        self.temporalStateCount = temporalStateCount
        self.stateFP32LE = stateFP32LE
        self.parametersFP32LE = parametersFP32LE
        self.temporalStateFP32LE = temporalStateFP32LE
        self.transportFP32LE = transportFP32LE
        self.ledgerHead = ledgerHead
    }

    public func validate() throws {
        guard version == Self.formatVersion,
              absoluteTimeSeconds.isFinite,
              absoluteTimeSeconds >= 0,
              speciesCount > 0,
              laneCount > 0,
              parameterStride > 0 else {
            throw VivoArtifactValidationError.invalid("checkpoint header is invalid")
        }
        let stateElements = try checkedProduct(UInt64(speciesCount), UInt64(laneCount), label: "checkpoint state")
        let parameterElements = try checkedProduct(UInt64(parameterCount), UInt64(parameterStride), label: "checkpoint parameters")
        let temporalElements = try checkedProduct(
            try checkedProduct(UInt64(temporalStateCount), UInt64(laneCount), label: "checkpoint temporal state"),
            2,
            label: "checkpoint temporal float2"
        )
        try requireByteCount(stateFP32LE, elements: stateElements, label: "state")
        try requireByteCount(parametersFP32LE, elements: parameterElements, label: "parameters")
        try requireByteCount(temporalStateFP32LE, elements: temporalElements, label: "temporal state")
        try requireByteCount(transportFP32LE, elements: UInt64(speciesCount) * 4, label: "transport")
    }

    public func fingerprint() throws -> VivoFingerprint {
        try validate()
        return try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(self))
    }

    private func requireByteCount(_ data: Data, elements: UInt64, label: String) throws {
        let expected = try checkedProduct(elements, 4, label: "\(label) bytes")
        guard expected <= UInt64(Int.max), data.count == Int(expected) else {
            throw VivoArtifactValidationError.invalid("checkpoint \(label) byte count does not match its shape")
        }
    }

    private func checkedProduct(_ lhs: UInt64, _ rhs: UInt64, label: String) throws -> UInt64 {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else {
            throw VivoArtifactValidationError.invalid("\(label) size overflow")
        }
        return result.partialValue
    }
}

public enum VivoLittleEndianFP32 {
    public static func encode(_ values: [Float]) -> Data {
        var data = Data(capacity: values.count * MemoryLayout<UInt32>.stride)
        for value in values {
            var bits = value.bitPattern.littleEndian
            withUnsafeBytes(of: &bits) { data.append(contentsOf: $0) }
        }
        return data
    }

    public static func decode(_ data: Data) throws -> [Float] {
        guard data.count.isMultiple(of: MemoryLayout<UInt32>.stride) else {
            throw VivoArtifactValidationError.invalid("FP32 little-endian payload length is not divisible by four")
        }
        var values: [Float] = []
        values.reserveCapacity(data.count / 4)
        var offset = 0
        while offset < data.count {
            let bits = data.withUnsafeBytes { raw -> UInt32 in
                var value: UInt32 = 0
                memcpy(&value, raw.baseAddress!.advanced(by: offset), 4)
                return UInt32(littleEndian: value)
            }
            values.append(Float(bitPattern: bits))
            offset += 4
        }
        return values
    }
}
