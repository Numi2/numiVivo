import CryptoKit
import Foundation

public struct VivoCampaignJobResourceEstimate: Codable, Equatable, Sendable {
    public let jobIndex: UInt64
    public let jobDigest: String
    public let topologyFingerprint: String
    public let mode: VivoAdaptiveFidelityMode
    public let stateBytes: UInt64
    public let scratchBytes: UInt64
    public let immutableTableBytes: UInt64
    public let expectedLogicalSteps: UInt64
    public let preferredLaneMultiple: UInt32

    public init(
        jobIndex: UInt64,
        jobDigest: String,
        topologyFingerprint: String,
        mode: VivoAdaptiveFidelityMode,
        stateBytes: UInt64,
        scratchBytes: UInt64,
        immutableTableBytes: UInt64,
        expectedLogicalSteps: UInt64,
        preferredLaneMultiple: UInt32 = 1
    ) {
        self.jobIndex = jobIndex
        self.jobDigest = jobDigest
        self.topologyFingerprint = topologyFingerprint
        self.mode = mode
        self.stateBytes = stateBytes
        self.scratchBytes = scratchBytes
        self.immutableTableBytes = immutableTableBytes
        self.expectedLogicalSteps = expectedLogicalSteps
        self.preferredLaneMultiple = preferredLaneMultiple
    }

    public var mutableBytes: UInt64? {
        let result = stateBytes.addingReportingOverflow(scratchBytes)
        return result.overflow ? nil : result.partialValue
    }

    public func validate() throws {
        guard jobDigest.utf8.count == 64,
              topologyFingerprint.utf8.count == 64,
              stateBytes > 0,
              expectedLogicalSteps > 0,
              preferredLaneMultiple > 0,
              mutableBytes != nil else {
            throw VivoCampaignBatchError.invalidEstimate(jobIndex)
        }
    }
}

public struct VivoCampaignBatchPolicy: Codable, Equatable, Sendable {
    public var workingSetFraction: Double
    public var immutableCacheFraction: Double
    public var minimumJobsPerBatch: UInt32
    public var maximumJobsPerBatch: UInt32
    public var maximumConcurrentBatches: UInt32
    public var reserveBytes: UInt64
    public var permitMixedFidelity: Bool

    public init(
        workingSetFraction: Double = 0.78,
        immutableCacheFraction: Double = 0.15,
        minimumJobsPerBatch: UInt32 = 1,
        maximumJobsPerBatch: UInt32 = 65_536,
        maximumConcurrentBatches: UInt32 = 2,
        reserveBytes: UInt64 = 512 * 1_024 * 1_024,
        permitMixedFidelity: Bool = false
    ) {
        self.workingSetFraction = workingSetFraction
        self.immutableCacheFraction = immutableCacheFraction
        self.minimumJobsPerBatch = minimumJobsPerBatch
        self.maximumJobsPerBatch = maximumJobsPerBatch
        self.maximumConcurrentBatches = maximumConcurrentBatches
        self.reserveBytes = reserveBytes
        self.permitMixedFidelity = permitMixedFidelity
    }

    public func validate() throws {
        guard workingSetFraction > 0,
              workingSetFraction <= 0.95,
              immutableCacheFraction >= 0,
              immutableCacheFraction < workingSetFraction,
              minimumJobsPerBatch > 0,
              maximumJobsPerBatch >= minimumJobsPerBatch,
              maximumConcurrentBatches > 0 else {
            throw VivoCampaignBatchError.invalidPolicy
        }
    }
}

public struct VivoCampaignBatch: Codable, Equatable, Sendable {
    public let batchIndex: UInt32
    public let topologyFingerprint: String
    public let modes: [VivoAdaptiveFidelityMode]
    public let jobIndices: [UInt64]
    public let jobDigests: [String]
    public let immutableTableBytes: UInt64
    public let mutableBytes: UInt64
    public let totalResidentBytes: UInt64
    public let maximumExpectedLogicalSteps: UInt64
    public let preferredLaneMultiple: UInt32
    public let fingerprint: String
}

public struct VivoCampaignExecutionWave: Codable, Equatable, Sendable {
    public let waveIndex: UInt32
    public let batchIndices: [UInt32]
    public let totalResidentBytes: UInt64
}

public struct VivoCampaignBatchPlan: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let campaignFingerprint: String
    public let device: VivoAppleGPUExecutionProfile
    public let policy: VivoCampaignBatchPolicy
    public let usableWorkingSetBytes: UInt64
    public let batches: [VivoCampaignBatch]
    public let waves: [VivoCampaignExecutionWave]
    public let fingerprint: String
}

public enum VivoCampaignBatchError: Error, LocalizedError, Sendable {
    case invalidPolicy
    case invalidCampaign
    case invalidEstimate(UInt64)
    case duplicateEstimate(UInt64)
    case missingEstimate(UInt64)
    case unknownJob(UInt64)
    case workingSetUnavailable
    case jobExceedsWorkingSet(UInt64, UInt64, UInt64)
    case arithmeticOverflow
    case encoding(String)

    public var errorDescription: String? {
        switch self {
        case .invalidPolicy: return "campaign batch policy is invalid"
        case .invalidCampaign: return "campaign manifest is invalid"
        case .invalidEstimate(let job): return "resource estimate for job \(job) is invalid"
        case .duplicateEstimate(let job): return "resource estimate for job \(job) is duplicated"
        case .missingEstimate(let job): return "resource estimate for job \(job) is missing"
        case .unknownJob(let job): return "resource estimate references unknown job \(job)"
        case .workingSetUnavailable: return "Apple GPU working-set size is unavailable"
        case .jobExceedsWorkingSet(let job, let required, let available): return "job \(job) requires \(required) bytes; \(available) are available"
        case .arithmeticOverflow: return "campaign batch planning arithmetic overflow"
        case .encoding(let reason): return "campaign batch plan encoding failed: \(reason)"
        }
    }
}

/// Deterministic first-fit-decreasing scheduler. Jobs are grouped by topology so
/// immutable ProgramPack and mechanism tables are loaded once per batch. Unless
/// explicitly permitted, a batch also contains one numerical authority only.
public struct VivoCampaignBatchPlanner: Sendable {
    public let policy: VivoCampaignBatchPolicy

    public init(policy: VivoCampaignBatchPolicy = .init()) throws {
        try policy.validate()
        self.policy = policy
    }

    public func plan(
        campaign: VivoCampaignManifest,
        estimates: [VivoCampaignJobResourceEstimate],
        device: VivoAppleGPUExecutionProfile
    ) throws -> VivoCampaignBatchPlan {
        guard campaign.schemaVersion == 1,
              !campaign.fingerprint.isEmpty,
              campaign.jobs.map(\.index).sorted() == Array(Set(campaign.jobs.map(\.index))).sorted() else {
            throw VivoCampaignBatchError.invalidCampaign
        }
        guard device.recommendedWorkingSetBytes > policy.reserveBytes else {
            throw VivoCampaignBatchError.workingSetUnavailable
        }

        var estimateMap: [UInt64: VivoCampaignJobResourceEstimate] = [:]
        let knownJobs = Set(campaign.jobs.map(\.index))
        for estimate in estimates {
            try estimate.validate()
            guard knownJobs.contains(estimate.jobIndex) else {
                throw VivoCampaignBatchError.unknownJob(estimate.jobIndex)
            }
            guard estimateMap[estimate.jobIndex] == nil else {
                throw VivoCampaignBatchError.duplicateEstimate(estimate.jobIndex)
            }
            estimateMap[estimate.jobIndex] = estimate
        }
        for job in campaign.jobs where estimateMap[job.index] == nil {
            throw VivoCampaignBatchError.missingEstimate(job.index)
        }

        let fractionalBudget = UInt64(
            (Double(device.recommendedWorkingSetBytes) * policy.workingSetFraction).rounded(.down)
        )
        guard fractionalBudget > policy.reserveBytes else {
            throw VivoCampaignBatchError.workingSetUnavailable
        }
        let usable = fractionalBudget - policy.reserveBytes

        struct GroupKey: Hashable, Comparable {
            let topology: String
            let mode: String

            static func < (lhs: Self, rhs: Self) -> Bool {
                lhs.topology == rhs.topology ? lhs.mode < rhs.mode : lhs.topology < rhs.topology
            }
        }
        let grouped = Dictionary(grouping: estimates) { estimate in
            GroupKey(
                topology: estimate.topologyFingerprint,
                mode: policy.permitMixedFidelity ? "mixed" : estimate.mode.rawValue
            )
        }

        struct MutableBatch {
            let topology: String
            var modes: Set<VivoAdaptiveFidelityMode>
            var jobs: [VivoCampaignJobResourceEstimate]
            var immutableBytes: UInt64
            var mutableBytes: UInt64
            var maximumSteps: UInt64
            var laneMultiple: UInt32

            var residentBytes: UInt64 { immutableBytes + mutableBytes }
        }

        var mutableBatches: [MutableBatch] = []
        for key in grouped.keys.sorted() {
            let jobs = grouped[key]!.sorted { left, right in
                let leftBytes = left.mutableBytes ?? 0
                let rightBytes = right.mutableBytes ?? 0
                if leftBytes != rightBytes { return leftBytes > rightBytes }
                return left.jobIndex < right.jobIndex
            }
            for job in jobs {
                guard let jobMutable = job.mutableBytes else {
                    throw VivoCampaignBatchError.invalidEstimate(job.jobIndex)
                }
                let jobTotal = job.immutableTableBytes.addingReportingOverflow(jobMutable)
                guard !jobTotal.overflow else { throw VivoCampaignBatchError.arithmeticOverflow }
                guard jobTotal.partialValue <= usable else {
                    throw VivoCampaignBatchError.jobExceedsWorkingSet(job.jobIndex, jobTotal.partialValue, usable)
                }

                var selected: Int?
                var bestRemainder = UInt64.max
                for index in mutableBatches.indices {
                    let batch = mutableBatches[index]
                    guard batch.topology == key.topology,
                          batch.jobs.count < Int(policy.maximumJobsPerBatch),
                          policy.permitMixedFidelity || batch.modes == [job.mode] else {
                        continue
                    }
                    let immutable = max(batch.immutableBytes, job.immutableTableBytes)
                    let mutable = batch.mutableBytes.addingReportingOverflow(jobMutable)
                    guard !mutable.overflow else { throw VivoCampaignBatchError.arithmeticOverflow }
                    let resident = immutable.addingReportingOverflow(mutable.partialValue)
                    guard !resident.overflow, resident.partialValue <= usable else { continue }
                    let remainder = usable - resident.partialValue
                    if remainder < bestRemainder {
                        bestRemainder = remainder
                        selected = index
                    }
                }

                if let index = selected {
                    mutableBatches[index].jobs.append(job)
                    mutableBatches[index].modes.insert(job.mode)
                    mutableBatches[index].immutableBytes = max(mutableBatches[index].immutableBytes, job.immutableTableBytes)
                    mutableBatches[index].mutableBytes += jobMutable
                    mutableBatches[index].maximumSteps = max(mutableBatches[index].maximumSteps, job.expectedLogicalSteps)
                    mutableBatches[index].laneMultiple = leastCommonMultiple(
                        mutableBatches[index].laneMultiple,
                        job.preferredLaneMultiple
                    )
                } else {
                    mutableBatches.append(.init(
                        topology: key.topology,
                        modes: [job.mode],
                        jobs: [job],
                        immutableBytes: job.immutableTableBytes,
                        mutableBytes: jobMutable,
                        maximumSteps: job.expectedLogicalSteps,
                        laneMultiple: job.preferredLaneMultiple
                    ))
                }
            }
        }

        var batches: [VivoCampaignBatch] = []
        batches.reserveCapacity(mutableBatches.count)
        for (index, batch) in mutableBatches.enumerated() {
            guard index <= Int(UInt32.max) else { throw VivoCampaignBatchError.arithmeticOverflow }
            let orderedJobs = batch.jobs.sorted { $0.jobIndex < $1.jobIndex }
            let resident = batch.immutableBytes.addingReportingOverflow(batch.mutableBytes)
            guard !resident.overflow else { throw VivoCampaignBatchError.arithmeticOverflow }
            let unsigned = VivoCampaignBatch(
                batchIndex: UInt32(index),
                topologyFingerprint: batch.topology,
                modes: batch.modes.sorted { $0.rawValue < $1.rawValue },
                jobIndices: orderedJobs.map(\.jobIndex),
                jobDigests: orderedJobs.map(\.jobDigest),
                immutableTableBytes: batch.immutableBytes,
                mutableBytes: batch.mutableBytes,
                totalResidentBytes: resident.partialValue,
                maximumExpectedLogicalSteps: batch.maximumSteps,
                preferredLaneMultiple: batch.laneMultiple,
                fingerprint: ""
            )
            batches.append(.init(
                batchIndex: unsigned.batchIndex,
                topologyFingerprint: unsigned.topologyFingerprint,
                modes: unsigned.modes,
                jobIndices: unsigned.jobIndices,
                jobDigests: unsigned.jobDigests,
                immutableTableBytes: unsigned.immutableTableBytes,
                mutableBytes: unsigned.mutableBytes,
                totalResidentBytes: unsigned.totalResidentBytes,
                maximumExpectedLogicalSteps: unsigned.maximumExpectedLogicalSteps,
                preferredLaneMultiple: unsigned.preferredLaneMultiple,
                fingerprint: try fingerprint(unsigned)
            ))
        }

        let waves = try scheduleWaves(batches: batches, usableBytes: usable)
        let unsigned = VivoCampaignBatchPlan(
            schemaVersion: 1,
            campaignFingerprint: campaign.fingerprint,
            device: device,
            policy: policy,
            usableWorkingSetBytes: usable,
            batches: batches,
            waves: waves,
            fingerprint: ""
        )
        return .init(
            schemaVersion: unsigned.schemaVersion,
            campaignFingerprint: unsigned.campaignFingerprint,
            device: unsigned.device,
            policy: unsigned.policy,
            usableWorkingSetBytes: unsigned.usableWorkingSetBytes,
            batches: unsigned.batches,
            waves: unsigned.waves,
            fingerprint: try fingerprint(unsigned)
        )
    }

    private func scheduleWaves(
        batches: [VivoCampaignBatch],
        usableBytes: UInt64
    ) throws -> [VivoCampaignExecutionWave] {
        var unscheduled = batches.sorted {
            if $0.totalResidentBytes != $1.totalResidentBytes {
                return $0.totalResidentBytes > $1.totalResidentBytes
            }
            return $0.batchIndex < $1.batchIndex
        }
        var waves: [VivoCampaignExecutionWave] = []
        while !unscheduled.isEmpty {
            var selected: [VivoCampaignBatch] = []
            var selectedIndices: [Int] = []
            var bytes: UInt64 = 0
            for (index, batch) in unscheduled.enumerated() {
                guard selected.count < Int(policy.maximumConcurrentBatches) else { break }
                let next = bytes.addingReportingOverflow(batch.totalResidentBytes)
                guard !next.overflow else { throw VivoCampaignBatchError.arithmeticOverflow }
                if next.partialValue <= usableBytes {
                    bytes = next.partialValue
                    selected.append(batch)
                    selectedIndices.append(index)
                }
            }
            guard !selected.isEmpty else {
                throw VivoCampaignBatchError.workingSetUnavailable
            }
            guard waves.count <= Int(UInt32.max) else { throw VivoCampaignBatchError.arithmeticOverflow }
            waves.append(.init(
                waveIndex: UInt32(waves.count),
                batchIndices: selected.map(\.batchIndex).sorted(),
                totalResidentBytes: bytes
            ))
            for index in selectedIndices.reversed() { unscheduled.remove(at: index) }
        }
        return waves
    }

    private func leastCommonMultiple(_ left: UInt32, _ right: UInt32) -> UInt32 {
        let divisor = greatestCommonDivisor(left, right)
        let reduced = left / max(divisor, 1)
        let result = UInt64(reduced) * UInt64(right)
        return result > UInt64(UInt32.max) ? UInt32.max : UInt32(result)
    }

    private func greatestCommonDivisor(_ left: UInt32, _ right: UInt32) -> UInt32 {
        var a = left
        var b = right
        while b != 0 {
            let remainder = a % b
            a = b
            b = remainder
        }
        return a
    }

    private func fingerprint<T: Encodable>(_ value: T) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        do {
            return SHA256.hash(data: try encoder.encode(value)).map { String(format: "%02x", $0) }.joined()
        } catch {
            throw VivoCampaignBatchError.encoding(error.localizedDescription)
        }
    }
}
