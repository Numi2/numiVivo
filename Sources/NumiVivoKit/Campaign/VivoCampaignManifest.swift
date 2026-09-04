import CryptoKit
import Foundation

public enum VivoCampaignSamplingMode: String, Codable, CaseIterable, Sendable {
    case cartesian
    case latinHypercube
    case halton
}

public enum VivoCampaignAxisTransform: String, Codable, CaseIterable, Sendable {
    case linear
    case logarithmic10
}

public struct VivoCampaignArtifactReferences: Codable, Equatable, Sendable {
    public let programPackFingerprint: String
    public let mechanismPackFingerprint: String
    public let hostContextFingerprint: String
    public let experimentDefinitionFingerprint: String?
    public let couplingPackFingerprint: String?

    public init(
        programPackFingerprint: String,
        mechanismPackFingerprint: String,
        hostContextFingerprint: String,
        experimentDefinitionFingerprint: String? = nil,
        couplingPackFingerprint: String? = nil
    ) {
        self.programPackFingerprint = programPackFingerprint
        self.mechanismPackFingerprint = mechanismPackFingerprint
        self.hostContextFingerprint = hostContextFingerprint
        self.experimentDefinitionFingerprint = experimentDefinitionFingerprint
        self.couplingPackFingerprint = couplingPackFingerprint
    }

    public func validate() throws {
        for (name, value) in [
            ("programPackFingerprint", programPackFingerprint),
            ("mechanismPackFingerprint", mechanismPackFingerprint),
            ("hostContextFingerprint", hostContextFingerprint)
        ] {
            guard Self.isFingerprint(value) else {
                throw VivoCampaignError.invalidDefinition("\(name) is not a SHA-256 fingerprint")
            }
        }
        for (name, value) in [
            ("experimentDefinitionFingerprint", experimentDefinitionFingerprint),
            ("couplingPackFingerprint", couplingPackFingerprint)
        ] {
            if let value, !Self.isFingerprint(value) {
                throw VivoCampaignError.invalidDefinition("\(name) is not a SHA-256 fingerprint")
            }
        }
    }

    private static func isFingerprint(_ value: String) -> Bool {
        value.utf8.count == 64 && value.allSatisfy { $0.isHexDigit && !$0.isUppercase }
    }
}

public struct VivoCampaignAxis: Codable, Equatable, Sendable {
    public let parameterID: String
    public let unit: String
    public let transform: VivoCampaignAxisTransform
    public let explicitValues: [Double]?
    public let lowerBound: Double?
    public let upperBound: Double?
    public let sampleCount: UInt32?

    public init(
        parameterID: String,
        unit: String,
        transform: VivoCampaignAxisTransform = .linear,
        explicitValues: [Double]
    ) {
        self.parameterID = parameterID
        self.unit = unit
        self.transform = transform
        self.explicitValues = explicitValues
        self.lowerBound = nil
        self.upperBound = nil
        self.sampleCount = UInt32(exactly: explicitValues.count)
    }

    public init(
        parameterID: String,
        unit: String,
        transform: VivoCampaignAxisTransform = .linear,
        lowerBound: Double,
        upperBound: Double,
        sampleCount: UInt32
    ) {
        self.parameterID = parameterID
        self.unit = unit
        self.transform = transform
        self.explicitValues = nil
        self.lowerBound = lowerBound
        self.upperBound = upperBound
        self.sampleCount = sampleCount
    }

    public func validate() throws {
        guard !parameterID.isEmpty,
              parameterID.utf8.count <= 256,
              !unit.isEmpty,
              unit.utf8.count <= 128 else {
            throw VivoCampaignError.invalidAxis(parameterID, "identifier or unit")
        }
        if let values = explicitValues {
            guard lowerBound == nil,
                  upperBound == nil,
                  !values.isEmpty,
                  values.count <= Int(UInt32.max),
                  values.allSatisfy(\.isFinite),
                  Set(values.map(\.bitPattern)).count == values.count else {
                throw VivoCampaignError.invalidAxis(parameterID, "explicit value set")
            }
            if transform == .logarithmic10, values.contains(where: { $0 <= 0 }) {
                throw VivoCampaignError.invalidAxis(parameterID, "logarithmic values must be positive")
            }
            return
        }
        guard let lowerBound,
              let upperBound,
              let sampleCount,
              lowerBound.isFinite,
              upperBound.isFinite,
              lowerBound < upperBound,
              sampleCount > 0 else {
            throw VivoCampaignError.invalidAxis(parameterID, "bounds or sample count")
        }
        if transform == .logarithmic10, lowerBound <= 0 {
            throw VivoCampaignError.invalidAxis(parameterID, "logarithmic lower bound must be positive")
        }
    }

    func cartesianValues() throws -> [Double] {
        if let explicitValues { return explicitValues }
        guard let lowerBound, let upperBound, let sampleCount else {
            throw VivoCampaignError.invalidAxis(parameterID, "missing range")
        }
        if sampleCount == 1 { return [midpoint(lowerBound, upperBound)] }
        return (0..<sampleCount).map { index in
            interpolate(
                fraction: Double(index) / Double(sampleCount - 1),
                lower: lowerBound,
                upper: upperBound
            )
        }
    }

    func sample(fraction: Double) throws -> Double {
        if let values = explicitValues {
            let clamped = min(max(fraction, 0), 1.nextDown)
            let index = min(values.count - 1, Int(clamped * Double(values.count)))
            return values[index]
        }
        guard let lowerBound, let upperBound else {
            throw VivoCampaignError.invalidAxis(parameterID, "missing range")
        }
        return interpolate(fraction: fraction, lower: lowerBound, upper: upperBound)
    }

    private func interpolate(fraction: Double, lower: Double, upper: Double) -> Double {
        switch transform {
        case .linear:
            return lower + min(max(fraction, 0), 1) * (upper - lower)
        case .logarithmic10:
            let low = log10(lower)
            let high = log10(upper)
            return pow(10, low + min(max(fraction, 0), 1) * (high - low))
        }
    }

    private func midpoint(_ lower: Double, _ upper: Double) -> Double {
        transform == .linear ? 0.5 * (lower + upper) : sqrt(lower * upper)
    }
}

public struct VivoCampaignDefinition: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let name: String
    public let artifacts: VivoCampaignArtifactReferences
    public let sampling: VivoCampaignSamplingMode
    public let axes: [VivoCampaignAxis]
    public let candidateCount: UInt32?
    public let replicateCount: UInt32
    public let baseSeed: UInt64
    public let maximumJobs: UInt64
    public let tags: [String: String]

    public init(
        name: String,
        artifacts: VivoCampaignArtifactReferences,
        sampling: VivoCampaignSamplingMode,
        axes: [VivoCampaignAxis],
        candidateCount: UInt32? = nil,
        replicateCount: UInt32 = 1,
        baseSeed: UInt64,
        maximumJobs: UInt64 = 1_000_000,
        tags: [String: String] = [:]
    ) {
        self.schemaVersion = 1
        self.name = name
        self.artifacts = artifacts
        self.sampling = sampling
        self.axes = axes
        self.candidateCount = candidateCount
        self.replicateCount = replicateCount
        self.baseSeed = baseSeed
        self.maximumJobs = maximumJobs
        self.tags = tags
    }

    public func validate() throws {
        guard schemaVersion == 1,
              !name.isEmpty,
              name.utf8.count <= 256,
              !axes.isEmpty,
              axes.count <= 256,
              Set(axes.map(\.parameterID)).count == axes.count,
              replicateCount > 0,
              maximumJobs > 0,
              tags.count <= 4_096,
              tags.allSatisfy({ $0.key.utf8.count <= 256 && $0.value.utf8.count <= 4_096 }) else {
            throw VivoCampaignError.invalidDefinition("identity, axes, replicates, limits or tags")
        }
        try artifacts.validate()
        for axis in axes { try axis.validate() }
        switch sampling {
        case .cartesian:
            if candidateCount != nil {
                throw VivoCampaignError.invalidDefinition("cartesian campaigns derive candidate count from axes")
            }
        case .latinHypercube, .halton:
            guard let candidateCount, candidateCount > 0 else {
                throw VivoCampaignError.invalidDefinition("sampled campaigns require a positive candidateCount")
            }
        }
    }
}

public struct VivoCampaignAssignment: Codable, Hashable, Sendable {
    public let parameterID: String
    public let unit: String
    public let value: Double
}

public struct VivoCampaignCandidate: Codable, Equatable, Sendable {
    public let index: UInt32
    public let assignments: [VivoCampaignAssignment]
    public let fingerprint: String
}

public struct VivoCampaignJob: Codable, Equatable, Sendable {
    public let index: UInt64
    public let candidateIndex: UInt32
    public let replicateIndex: UInt32
    public let seed: UInt64
    public let candidateFingerprint: String
    public let previousJobDigest: String
    public let digest: String
}

public struct VivoCampaignManifest: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let definitionFingerprint: String
    public let name: String
    public let artifacts: VivoCampaignArtifactReferences
    public let sampling: VivoCampaignSamplingMode
    public let candidates: [VivoCampaignCandidate]
    public let jobs: [VivoCampaignJob]
    public let tags: [String: String]
    public let ledgerHead: String
    public let fingerprint: String
}

public enum VivoCampaignError: Error, LocalizedError, Sendable {
    case invalidDefinition(String)
    case invalidAxis(String, String)
    case arithmeticOverflow
    case jobLimitExceeded(UInt64)
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidDefinition(let reason): return "invalid campaign definition: \(reason)"
        case .invalidAxis(let parameter, let reason): return "invalid campaign axis \(parameter): \(reason)"
        case .arithmeticOverflow: return "campaign expansion arithmetic overflow"
        case .jobLimitExceeded(let count): return "campaign expands to \(count) jobs, exceeding its limit"
        case .encoding(let reason): return "campaign encoding failed: \(reason)"
        }
    }
}

public struct VivoCampaignCompiler: Sendable {
    public init() {}

    public func compile(_ definition: VivoCampaignDefinition) throws -> VivoCampaignManifest {
        try definition.validate()
        let definitionFingerprint = try fingerprint(definition)
        let assignments = try candidateAssignments(definition)
        guard assignments.count <= Int(UInt32.max) else {
            throw VivoCampaignError.arithmeticOverflow
        }

        var candidates: [VivoCampaignCandidate] = []
        candidates.reserveCapacity(assignments.count)
        for (index, assignment) in assignments.enumerated() {
            let ordered = assignment.sorted { $0.parameterID < $1.parameterID }
            candidates.append(.init(
                index: UInt32(index),
                assignments: ordered,
                fingerprint: try fingerprint(ordered)
            ))
        }

        let jobCount = UInt64(candidates.count).multipliedReportingOverflow(by: UInt64(definition.replicateCount))
        guard !jobCount.overflow else { throw VivoCampaignError.arithmeticOverflow }
        guard jobCount.partialValue <= definition.maximumJobs else {
            throw VivoCampaignError.jobLimitExceeded(jobCount.partialValue)
        }
        guard jobCount.partialValue <= UInt64(Int.max) else {
            throw VivoCampaignError.arithmeticOverflow
        }

        var jobs: [VivoCampaignJob] = []
        jobs.reserveCapacity(Int(jobCount.partialValue))
        var previous = String(repeating: "0", count: 64)
        var jobIndex: UInt64 = 0
        for candidate in candidates {
            for replicate in 0..<definition.replicateCount {
                let seed = derivedSeed(
                    baseSeed: definition.baseSeed,
                    definitionFingerprint: definitionFingerprint,
                    candidateFingerprint: candidate.fingerprint,
                    replicate: replicate
                )
                let digest = jobDigest(
                    previous: previous,
                    index: jobIndex,
                    candidate: candidate.fingerprint,
                    replicate: replicate,
                    seed: seed
                )
                jobs.append(.init(
                    index: jobIndex,
                    candidateIndex: candidate.index,
                    replicateIndex: replicate,
                    seed: seed,
                    candidateFingerprint: candidate.fingerprint,
                    previousJobDigest: previous,
                    digest: digest
                ))
                previous = digest
                jobIndex &+= 1
            }
        }

        let unsigned = VivoCampaignManifest(
            schemaVersion: 1,
            definitionFingerprint: definitionFingerprint,
            name: definition.name,
            artifacts: definition.artifacts,
            sampling: definition.sampling,
            candidates: candidates,
            jobs: jobs,
            tags: definition.tags,
            ledgerHead: previous,
            fingerprint: ""
        )
        return .init(
            schemaVersion: unsigned.schemaVersion,
            definitionFingerprint: unsigned.definitionFingerprint,
            name: unsigned.name,
            artifacts: unsigned.artifacts,
            sampling: unsigned.sampling,
            candidates: unsigned.candidates,
            jobs: unsigned.jobs,
            tags: unsigned.tags,
            ledgerHead: unsigned.ledgerHead,
            fingerprint: try fingerprint(unsigned)
        )
    }

    private func candidateAssignments(_ definition: VivoCampaignDefinition) throws -> [[VivoCampaignAssignment]] {
        switch definition.sampling {
        case .cartesian:
            var result: [[VivoCampaignAssignment]] = [[]]
            for axis in definition.axes.sorted(by: { $0.parameterID < $1.parameterID }) {
                let values = try axis.cartesianValues()
                let multiplied = result.count.multipliedReportingOverflow(by: values.count)
                guard !multiplied.overflow, UInt64(multiplied.partialValue) <= definition.maximumJobs else {
                    throw VivoCampaignError.jobLimitExceeded(UInt64.max)
                }
                var expanded: [[VivoCampaignAssignment]] = []
                expanded.reserveCapacity(multiplied.partialValue)
                for prefix in result {
                    for value in values {
                        var candidate = prefix
                        candidate.append(.init(parameterID: axis.parameterID, unit: axis.unit, value: value))
                        expanded.append(candidate)
                    }
                }
                result = expanded
            }
            return result

        case .latinHypercube:
            guard let count = definition.candidateCount else {
                throw VivoCampaignError.invalidDefinition("candidateCount")
            }
            let axes = definition.axes.sorted(by: { $0.parameterID < $1.parameterID })
            var result = [[VivoCampaignAssignment]](repeating: [], count: Int(count))
            for (axisIndex, axis) in axes.enumerated() {
                let permutation = deterministicPermutation(
                    count: Int(count),
                    seed: definition.baseSeed,
                    discriminator: UInt64(axisIndex)
                )
                for candidate in 0..<Int(count) {
                    // Centered strata avoid unrecorded pseudorandom jitter.
                    let fraction = (Double(permutation[candidate]) + 0.5) / Double(count)
                    result[candidate].append(.init(
                        parameterID: axis.parameterID,
                        unit: axis.unit,
                        value: try axis.sample(fraction: fraction)
                    ))
                }
            }
            return result

        case .halton:
            guard let count = definition.candidateCount else {
                throw VivoCampaignError.invalidDefinition("candidateCount")
            }
            let axes = definition.axes.sorted(by: { $0.parameterID < $1.parameterID })
            guard axes.count <= Self.primes.count else {
                throw VivoCampaignError.invalidDefinition("Halton campaigns support at most \(Self.primes.count) axes")
            }
            return try (0..<Int(count)).map { candidate in
                try axes.enumerated().map { axisIndex, axis in
                    .init(
                        parameterID: axis.parameterID,
                        unit: axis.unit,
                        value: try axis.sample(
                            fraction: radicalInverse(index: UInt64(candidate + 1), base: Self.primes[axisIndex])
                        )
                    )
                }
            }
        }
    }

    private func deterministicPermutation(count: Int, seed: UInt64, discriminator: UInt64) -> [Int] {
        guard count > 1 else { return Array(0..<count) }
        return (0..<count).sorted { left, right in
            let leftHash = orderingHash(seed: seed, discriminator: discriminator, value: UInt64(left))
            let rightHash = orderingHash(seed: seed, discriminator: discriminator, value: UInt64(right))
            return leftHash == rightHash ? left < right : leftHash < rightHash
        }
    }

    private func orderingHash(seed: UInt64, discriminator: UInt64, value: UInt64) -> String {
        var hasher = SHA256()
        update(&hasher, seed)
        update(&hasher, discriminator)
        update(&hasher, value)
        return hasher.finalize().map { String(format: "%02x", $0) }.joined()
    }

    private func derivedSeed(
        baseSeed: UInt64,
        definitionFingerprint: String,
        candidateFingerprint: String,
        replicate: UInt32
    ) -> UInt64 {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-campaign-seed-v1".utf8))
        update(&hasher, baseSeed)
        hasher.update(data: Data(definitionFingerprint.utf8))
        hasher.update(data: Data(candidateFingerprint.utf8))
        update(&hasher, replicate)
        let digest = Array(hasher.finalize())
        return digest.prefix(8).enumerated().reduce(UInt64(0)) { result, pair in
            result | (UInt64(pair.element) << UInt64(pair.offset * 8))
        }
    }

    private func jobDigest(
        previous: String,
        index: UInt64,
        candidate: String,
        replicate: UInt32,
        seed: UInt64
    ) -> String {
        var hasher = SHA256()
        hasher.update(data: Data("numivivo-campaign-job-v1".utf8))
        hasher.update(data: Data(previous.utf8))
        update(&hasher, index)
        hasher.update(data: Data(candidate.utf8))
        update(&hasher, replicate)
        update(&hasher, seed)
        return hasher.finalize().map { String(format: "%02x", $0) }.joined()
    }

    private func radicalInverse(index: UInt64, base: UInt32) -> Double {
        var value = index
        var denominator = Double(base)
        var result = 0.0
        while value > 0 {
            result += Double(value % UInt64(base)) / denominator
            value /= UInt64(base)
            denominator *= Double(base)
        }
        return result
    }

    private func fingerprint<T: Encodable>(_ value: T) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        do {
            return SHA256.hash(data: try encoder.encode(value)).map { String(format: "%02x", $0) }.joined()
        } catch {
            throw VivoCampaignError.encoding(error.localizedDescription)
        }
    }

    private func update<T: FixedWidthInteger>(_ hasher: inout SHA256, _ value: T) {
        var littleEndian = value.littleEndian
        withUnsafeBytes(of: &littleEndian) { hasher.update(data: Data($0)) }
    }

    private static let primes: [UInt32] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
        131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
        197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
        271, 277, 281, 283, 293, 307, 311, 313
    ]
}
