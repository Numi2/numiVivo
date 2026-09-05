import CryptoKit
import Foundation
#if canImport(Metal)
@preconcurrency import Metal
#endif

public enum VivoAdaptiveFidelityMode: String, Codable, CaseIterable, Sendable {
    case logic
    case deterministicRK2
    case tauLeap
    case exactSSA
    case spatialSplit
    case tissueCoupled
}

public enum VivoStateRepresentation: String, Codable, Sendable {
    case continuousFP32
    case discreteUInt32
    case mixedContinuousDiscrete
}

public struct VivoAppleGPUExecutionProfile: Codable, Equatable, Sendable {
    public let deviceName: String
    public let registryID: UInt64
    public let threadExecutionWidth: UInt32
    public let maximumThreadsPerThreadgroup: UInt32
    public let recommendedWorkingSetBytes: UInt64
    public let hasUnifiedMemory: Bool
    public let supportsDynamicLibraries: Bool
    public let supportsRaytracing: Bool

    public init(
        deviceName: String,
        registryID: UInt64,
        threadExecutionWidth: UInt32,
        maximumThreadsPerThreadgroup: UInt32,
        recommendedWorkingSetBytes: UInt64,
        hasUnifiedMemory: Bool,
        supportsDynamicLibraries: Bool,
        supportsRaytracing: Bool
    ) {
        self.deviceName = deviceName
        self.registryID = registryID
        self.threadExecutionWidth = max(1, threadExecutionWidth)
        self.maximumThreadsPerThreadgroup = max(1, maximumThreadsPerThreadgroup)
        self.recommendedWorkingSetBytes = recommendedWorkingSetBytes
        self.hasUnifiedMemory = hasUnifiedMemory
        self.supportsDynamicLibraries = supportsDynamicLibraries
        self.supportsRaytracing = supportsRaytracing
    }
}

#if canImport(Metal)
public extension VivoAppleGPUExecutionProfile {
    static func inspect(device: MTLDevice, representativePipeline: MTLComputePipelineState? = nil) -> Self {
        let executionWidth = representativePipeline?.threadExecutionWidth ?? 32
        let maximumThreads = representativePipeline?.maxTotalThreadsPerThreadgroup ?? 256
        return .init(
            deviceName: device.name,
            registryID: device.registryID,
            threadExecutionWidth: UInt32(clamping: executionWidth),
            maximumThreadsPerThreadgroup: UInt32(clamping: maximumThreads),
            recommendedWorkingSetBytes: UInt64(device.recommendedMaxWorkingSetSize),
            hasUnifiedMemory: device.hasUnifiedMemory,
            supportsDynamicLibraries: device.supportsDynamicLibraries,
            supportsRaytracing: device.supportsRaytracing
        )
    }
}
#endif

/// Compiler and runtime statistics for one structurally homogeneous reaction cohort.
public struct VivoMechanismWorkload: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let laneCount: UInt64
    public let speciesCount: UInt32
    public let reactionCount: UInt32
    public let parameterCount: UInt32
    public let minimumRepresentedCount: Double
    public let medianRepresentedCount: Double
    public let maximumPropensityPerSecond: Double
    public let expectedFiringsPerLanePerStep: Double
    public let stiffnessRatio: Double
    public let criticalReactionFraction: Double
    public let relativeParameterUncertainty: Double
    public let requestedRelativeError: Double
    public let stateRepresentation: VivoStateRepresentation
    public let requiresSpatialTransport: Bool
    public let requiresTissueCoupling: Bool
    public let topologyMayChange: Bool
    public let containsDelayedReactions: Bool
    public let containsIrreversibleState: Bool

    public init(
        cohortID: UInt32,
        laneCount: UInt64,
        speciesCount: UInt32,
        reactionCount: UInt32,
        parameterCount: UInt32,
        minimumRepresentedCount: Double,
        medianRepresentedCount: Double,
        maximumPropensityPerSecond: Double,
        expectedFiringsPerLanePerStep: Double,
        stiffnessRatio: Double,
        criticalReactionFraction: Double,
        relativeParameterUncertainty: Double,
        requestedRelativeError: Double,
        stateRepresentation: VivoStateRepresentation,
        requiresSpatialTransport: Bool = false,
        requiresTissueCoupling: Bool = false,
        topologyMayChange: Bool = false,
        containsDelayedReactions: Bool = false,
        containsIrreversibleState: Bool = false
    ) {
        self.cohortID = cohortID
        self.laneCount = laneCount
        self.speciesCount = speciesCount
        self.reactionCount = reactionCount
        self.parameterCount = parameterCount
        self.minimumRepresentedCount = minimumRepresentedCount
        self.medianRepresentedCount = medianRepresentedCount
        self.maximumPropensityPerSecond = maximumPropensityPerSecond
        self.expectedFiringsPerLanePerStep = expectedFiringsPerLanePerStep
        self.stiffnessRatio = stiffnessRatio
        self.criticalReactionFraction = criticalReactionFraction
        self.relativeParameterUncertainty = relativeParameterUncertainty
        self.requestedRelativeError = requestedRelativeError
        self.stateRepresentation = stateRepresentation
        self.requiresSpatialTransport = requiresSpatialTransport
        self.requiresTissueCoupling = requiresTissueCoupling
        self.topologyMayChange = topologyMayChange
        self.containsDelayedReactions = containsDelayedReactions
        self.containsIrreversibleState = containsIrreversibleState
    }

    public func validate() throws {
        guard laneCount > 0, speciesCount > 0 else {
            throw VivoAdaptiveFidelityError.invalidWorkload(cohortID, "laneCount and speciesCount must be positive")
        }
        let finite = [
            minimumRepresentedCount,
            medianRepresentedCount,
            maximumPropensityPerSecond,
            expectedFiringsPerLanePerStep,
            stiffnessRatio,
            criticalReactionFraction,
            relativeParameterUncertainty,
            requestedRelativeError
        ].allSatisfy(\.isFinite)
        guard finite else {
            throw VivoAdaptiveFidelityError.invalidWorkload(cohortID, "workload contains a non-finite statistic")
        }
        guard minimumRepresentedCount >= 0,
              medianRepresentedCount >= minimumRepresentedCount,
              maximumPropensityPerSecond >= 0,
              expectedFiringsPerLanePerStep >= 0,
              stiffnessRatio >= 1,
              (0...1).contains(criticalReactionFraction),
              relativeParameterUncertainty >= 0,
              requestedRelativeError > 0,
              requestedRelativeError <= 1 else {
            throw VivoAdaptiveFidelityError.invalidWorkload(cohortID, "workload statistics are outside valid bounds")
        }
    }
}

public struct VivoAdaptiveFidelityPolicy: Codable, Equatable, Sendable {
    public var exactCountThreshold: Double
    public var exactCriticalFractionThreshold: Double
    public var exactMaximumReactions: UInt32
    public var exactMaximumExpectedFirings: Double
    public var deterministicCountThreshold: Double
    public var deterministicUncertaintyFloor: Double
    public var deterministicMaximumRequestedError: Double
    public var tauMaximumExpectedFirings: Double
    public var stiffRatioThreshold: Double
    public var workingSetFraction: Double
    public var minimumLanesPerChunk: UInt32
    public var maximumLanesPerChunk: UInt32
    public var preferredThreadgroupsPerChunk: UInt32
    public var stateBufferMultiplicity: UInt32
    public var extraScratchFraction: Double

    public init(
        exactCountThreshold: Double = 64,
        exactCriticalFractionThreshold: Double = 0.2,
        exactMaximumReactions: UInt32 = 2_048,
        exactMaximumExpectedFirings: Double = 24,
        deterministicCountThreshold: Double = 20_000,
        deterministicUncertaintyFloor: Double = 0.02,
        deterministicMaximumRequestedError: Double = 0.05,
        tauMaximumExpectedFirings: Double = 4_096,
        stiffRatioThreshold: Double = 1_000,
        workingSetFraction: Double = 0.72,
        minimumLanesPerChunk: UInt32 = 256,
        maximumLanesPerChunk: UInt32 = 1_048_576,
        preferredThreadgroupsPerChunk: UInt32 = 256,
        stateBufferMultiplicity: UInt32 = 4,
        extraScratchFraction: Double = 0.25
    ) {
        self.exactCountThreshold = exactCountThreshold
        self.exactCriticalFractionThreshold = exactCriticalFractionThreshold
        self.exactMaximumReactions = exactMaximumReactions
        self.exactMaximumExpectedFirings = exactMaximumExpectedFirings
        self.deterministicCountThreshold = deterministicCountThreshold
        self.deterministicUncertaintyFloor = deterministicUncertaintyFloor
        self.deterministicMaximumRequestedError = deterministicMaximumRequestedError
        self.tauMaximumExpectedFirings = tauMaximumExpectedFirings
        self.stiffRatioThreshold = stiffRatioThreshold
        self.workingSetFraction = workingSetFraction
        self.minimumLanesPerChunk = minimumLanesPerChunk
        self.maximumLanesPerChunk = maximumLanesPerChunk
        self.preferredThreadgroupsPerChunk = preferredThreadgroupsPerChunk
        self.stateBufferMultiplicity = stateBufferMultiplicity
        self.extraScratchFraction = extraScratchFraction
    }

    public func validate() throws {
        guard exactCountThreshold >= 0,
              (0...1).contains(exactCriticalFractionThreshold),
              exactMaximumReactions > 0,
              exactMaximumExpectedFirings > 0,
              deterministicCountThreshold >= exactCountThreshold,
              deterministicUncertaintyFloor >= 0,
              deterministicMaximumRequestedError > 0,
              tauMaximumExpectedFirings >= exactMaximumExpectedFirings,
              stiffRatioThreshold >= 1,
              workingSetFraction > 0,
              workingSetFraction <= 0.95,
              minimumLanesPerChunk > 0,
              maximumLanesPerChunk >= minimumLanesPerChunk,
              preferredThreadgroupsPerChunk > 0,
              stateBufferMultiplicity >= 2,
              extraScratchFraction >= 0 else {
            throw VivoAdaptiveFidelityError.invalidPolicy
        }
    }
}

public enum VivoBufferMigrationKind: String, Codable, Sendable {
    case none
    case reinterpretInPlace
    case continuousToDiscrete
    case discreteToContinuous
    case resizeAndCopy
    case rebuildSpatialTopology
    case rebuildTissueCoupling
}

public struct VivoBufferMigrationPlan: Codable, Equatable, Sendable {
    public let kind: VivoBufferMigrationKind
    public let sourceMode: VivoAdaptiveFidelityMode?
    public let destinationMode: VivoAdaptiveFidelityMode
    public let preserveRandomCounter: Bool
    public let preserveDelayedQueue: Bool
    public let requiresTransactionBarrier: Bool
    public let reason: String
}

public struct VivoFidelityCohortPlan: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let mode: VivoAdaptiveFidelityMode
    public let laneCount: UInt64
    public let lanesPerChunk: UInt32
    public let chunkCount: UInt32
    public let threadsPerThreadgroup: UInt32
    public let threadgroupsPerChunk: UInt32
    public let stateBytesPerChunk: UInt64
    public let scratchBytesPerChunk: UInt64
    public let totalResidentBytes: UInt64
    public let recommendedMaximumStep: Double
    public let requiresPrivateHeap: Bool
    public let migration: VivoBufferMigrationPlan
    public let rationale: [String]
}

public struct VivoAdaptiveFidelityPlan: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let device: VivoAppleGPUExecutionProfile
    public let workingSetBudgetBytes: UInt64
    public let cohorts: [VivoFidelityCohortPlan]
    public let totalResidentBytes: UInt64
    public let fingerprint: String
}

public enum VivoAdaptiveFidelityError: Error, LocalizedError, Sendable {
    case invalidPolicy
    case invalidWorkload(UInt32, String)
    case arithmeticOverflow(UInt32)
    case workingSetUnavailable
    case workingSetExceeded(required: UInt64, available: UInt64)

    public var errorDescription: String? {
        switch self {
        case .invalidPolicy: return "adaptive fidelity policy is invalid"
        case .invalidWorkload(let cohort, let reason): return "invalid workload for cohort \(cohort): \(reason)"
        case .arithmeticOverflow(let cohort): return "allocation arithmetic overflow for cohort \(cohort)"
        case .workingSetUnavailable: return "Apple GPU working-set size is unavailable"
        case .workingSetExceeded(let required, let available): return "planned resident bytes \(required) exceed budget \(available)"
        }
    }
}

/// Converts biological and numerical workload statistics into explicit execution,
/// allocation and state-migration decisions. The plan is deterministic for the
/// same workload, policy and GPU profile.
public struct VivoAdaptiveFidelityPlanner: Sendable {
    public let policy: VivoAdaptiveFidelityPolicy

    public init(policy: VivoAdaptiveFidelityPolicy = .init()) throws {
        try policy.validate()
        self.policy = policy
    }

    public func plan(
        workloads: [VivoMechanismWorkload],
        device: VivoAppleGPUExecutionProfile,
        previous: VivoAdaptiveFidelityPlan? = nil
    ) throws -> VivoAdaptiveFidelityPlan {
        guard device.recommendedWorkingSetBytes > 0 else {
            throw VivoAdaptiveFidelityError.workingSetUnavailable
        }
        for workload in workloads { try workload.validate() }

        let budgetDouble = Double(device.recommendedWorkingSetBytes) * policy.workingSetFraction
        guard budgetDouble.isFinite, budgetDouble > 0, budgetDouble <= Double(UInt64.max) else {
            throw VivoAdaptiveFidelityError.workingSetUnavailable
        }
        let budget = UInt64(budgetDouble.rounded(.down))
        let previousByCohort = Dictionary(uniqueKeysWithValues: (previous?.cohorts ?? []).map { ($0.cohortID, $0) })

        var plans: [VivoFidelityCohortPlan] = []
        var total: UInt64 = 0
        for workload in workloads.sorted(by: { $0.cohortID < $1.cohortID }) {
            let modeSelection = selectMode(workload)
            let mode = modeSelection.mode
            let elementBytes = bytesPerStateElement(workload.stateRepresentation)
            let fixedBytes = try fixedTableBytes(workload)
            let bytesPerLane = try checkedMultiply(
                UInt64(workload.speciesCount),
                UInt64(elementBytes),
                cohortID: workload.cohortID
            )
            let stateMultiplicity = UInt64(policy.stateBufferMultiplicity)
            let privateStatePerLane = try checkedMultiply(bytesPerLane, stateMultiplicity, cohortID: workload.cohortID)
            let modeScratchPerLane = scratchBytesPerLane(mode: mode, workload: workload)
            let perLane = try checkedAdd(privateStatePerLane, modeScratchPerLane, cohortID: workload.cohortID)

            let fairShare = max(
                UInt64(policy.minimumLanesPerChunk) * max(perLane, 1),
                budget / UInt64(max(workloads.count, 1))
            )
            let availableForLanes = fairShare > fixedBytes ? fairShare - fixedBytes : perLane
            let rawLanes = perLane == 0 ? workload.laneCount : availableForLanes / perLane
            let boundedLanes = min(
                workload.laneCount,
                UInt64(policy.maximumLanesPerChunk),
                max(UInt64(policy.minimumLanesPerChunk), rawLanes)
            )
            let executionWidth = UInt64(max(1, device.threadExecutionWidth))
            let alignedLanes = max(
                executionWidth,
                min(workload.laneCount, (boundedLanes / executionWidth) * executionWidth)
            )
            guard alignedLanes <= UInt64(UInt32.max) else {
                throw VivoAdaptiveFidelityError.arithmeticOverflow(workload.cohortID)
            }
            let lanesPerChunk = UInt32(alignedLanes)
            let chunks64 = (workload.laneCount + alignedLanes - 1) / alignedLanes
            guard chunks64 <= UInt64(UInt32.max) else {
                throw VivoAdaptiveFidelityError.arithmeticOverflow(workload.cohortID)
            }
            let chunkCount = UInt32(chunks64)

            let stateBytes = try checkedMultiply(privateStatePerLane, alignedLanes, cohortID: workload.cohortID)
            let scratchBase = try checkedMultiply(modeScratchPerLane, alignedLanes, cohortID: workload.cohortID)
            let extraScratch = UInt64((Double(stateBytes) * policy.extraScratchFraction).rounded(.up))
            let scratchBytes = try checkedAdd(scratchBase, extraScratch, cohortID: workload.cohortID)
            let resident = try checkedAdd(
                fixedBytes,
                try checkedAdd(stateBytes, scratchBytes, cohortID: workload.cohortID),
                cohortID: workload.cohortID
            )
            total = try checkedAdd(total, resident, cohortID: workload.cohortID)

            let threads = preferredThreads(device: device)
            let groups = min(
                policy.preferredThreadgroupsPerChunk,
                UInt32((alignedLanes + UInt64(threads) - 1) / UInt64(threads))
            )
            let previousMode = previousByCohort[workload.cohortID]?.mode
            let migration = migrationPlan(
                from: previousMode,
                to: mode,
                workload: workload
            )
            let maximumStep = recommendedMaximumStep(workload: workload, mode: mode)
            plans.append(.init(
                cohortID: workload.cohortID,
                mode: mode,
                laneCount: workload.laneCount,
                lanesPerChunk: lanesPerChunk,
                chunkCount: chunkCount,
                threadsPerThreadgroup: threads,
                threadgroupsPerChunk: max(1, groups),
                stateBytesPerChunk: stateBytes,
                scratchBytesPerChunk: scratchBytes,
                totalResidentBytes: resident,
                recommendedMaximumStep: maximumStep,
                requiresPrivateHeap: true,
                migration: migration,
                rationale: modeSelection.rationale
            ))
        }

        guard total <= budget else {
            throw VivoAdaptiveFidelityError.workingSetExceeded(required: total, available: budget)
        }

        let unsigned = VivoAdaptiveFidelityPlan(
            schemaVersion: 1,
            device: device,
            workingSetBudgetBytes: budget,
            cohorts: plans,
            totalResidentBytes: total,
            fingerprint: ""
        )
        let fingerprint = try Self.fingerprint(unsigned)
        return .init(
            schemaVersion: unsigned.schemaVersion,
            device: unsigned.device,
            workingSetBudgetBytes: unsigned.workingSetBudgetBytes,
            cohorts: unsigned.cohorts,
            totalResidentBytes: unsigned.totalResidentBytes,
            fingerprint: fingerprint
        )
    }

    private func selectMode(_ workload: VivoMechanismWorkload) -> (mode: VivoAdaptiveFidelityMode, rationale: [String]) {
        if workload.requiresTissueCoupling {
            return (.tissueCoupled, ["cohort exchanges state with the tissue-scale runtime"])
        }
        if workload.requiresSpatialTransport {
            return (.spatialSplit, ["cohort requires spatial transport or gradients"])
        }
        if workload.reactionCount == 0 {
            return (.logic, ["cohort has no kinetic reactions"])
        }

        let isDiscrete = workload.stateRepresentation != .continuousFP32
        let exactByCount = workload.minimumRepresentedCount <= policy.exactCountThreshold
        let exactByCriticality = workload.criticalReactionFraction >= policy.exactCriticalFractionThreshold
        let exactAffordable = workload.reactionCount <= policy.exactMaximumReactions &&
            workload.expectedFiringsPerLanePerStep <= policy.exactMaximumExpectedFirings
        if isDiscrete && (exactByCount || exactByCriticality) && exactAffordable {
            var reasons = ["discrete molecular state requires event-resolved execution"]
            if exactByCount { reasons.append("minimum molecular count is below exact threshold") }
            if exactByCriticality { reasons.append("critical reaction fraction exceeds exact threshold") }
            return (.exactSSA, reasons)
        }

        let deterministicByCount = workload.medianRepresentedCount >= policy.deterministicCountThreshold
        let deterministicErrorAcceptable = workload.requestedRelativeError >= policy.deterministicMaximumRequestedError
        let modelUncertaintyDominates = workload.relativeParameterUncertainty >= policy.deterministicUncertaintyFloor
        if deterministicByCount && deterministicErrorAcceptable && modelUncertaintyDominates {
            return (
                .deterministicRK2,
                [
                    "molecular populations are large enough for a continuous state",
                    "requested numerical error is below declared model uncertainty"
                ]
            )
        }

        if isDiscrete && workload.expectedFiringsPerLanePerStep <= policy.tauMaximumExpectedFirings {
            return (.tauLeap, ["discrete state is outside exact-SSA cost envelope but remains leap-compatible"])
        }

        return (
            .deterministicRK2,
            workload.stiffnessRatio >= policy.stiffRatioThreshold
                ? ["large event volume requires continuous integration", "stiffness requires runtime step restriction"]
                : ["large event volume requires continuous integration"]
        )
    }

    private func migrationPlan(
        from source: VivoAdaptiveFidelityMode?,
        to destination: VivoAdaptiveFidelityMode,
        workload: VivoMechanismWorkload
    ) -> VivoBufferMigrationPlan {
        guard let source else {
            return .init(
                kind: destination == .spatialSplit ? .rebuildSpatialTopology : destination == .tissueCoupled ? .rebuildTissueCoupling : .resizeAndCopy,
                sourceMode: nil,
                destinationMode: destination,
                preserveRandomCounter: true,
                preserveDelayedQueue: workload.containsDelayedReactions,
                requiresTransactionBarrier: true,
                reason: "initial allocation"
            )
        }
        guard source != destination else {
            return .init(
                kind: workload.topologyMayChange ? .resizeAndCopy : .none,
                sourceMode: source,
                destinationMode: destination,
                preserveRandomCounter: true,
                preserveDelayedQueue: workload.containsDelayedReactions,
                requiresTransactionBarrier: workload.topologyMayChange,
                reason: workload.topologyMayChange ? "topology changed within the same fidelity mode" : "fidelity and topology are unchanged"
            )
        }

        let discrete: Set<VivoAdaptiveFidelityMode> = [.tauLeap, .exactSSA]
        let kind: VivoBufferMigrationKind
        if destination == .spatialSplit { kind = .rebuildSpatialTopology }
        else if destination == .tissueCoupled { kind = .rebuildTissueCoupling }
        else if discrete.contains(source) && !discrete.contains(destination) { kind = .discreteToContinuous }
        else if !discrete.contains(source) && discrete.contains(destination) { kind = .continuousToDiscrete }
        else { kind = .reinterpretInPlace }

        return .init(
            kind: kind,
            sourceMode: source,
            destinationMode: destination,
            preserveRandomCounter: true,
            preserveDelayedQueue: workload.containsDelayedReactions,
            requiresTransactionBarrier: true,
            reason: "fidelity changed from \(source.rawValue) to \(destination.rawValue)"
        )
    }

    private func preferredThreads(device: VivoAppleGPUExecutionProfile) -> UInt32 {
        let width = max(1, device.threadExecutionWidth)
        let preferred = width * 8
        return min(max(width, preferred), device.maximumThreadsPerThreadgroup)
    }

    private func bytesPerStateElement(_ representation: VivoStateRepresentation) -> UInt32 {
        switch representation {
        case .continuousFP32, .discreteUInt32: return 4
        case .mixedContinuousDiscrete: return 8
        }
    }

    private func scratchBytesPerLane(mode: VivoAdaptiveFidelityMode, workload: VivoMechanismWorkload) -> UInt64 {
        switch mode {
        case .logic: return UInt64(workload.speciesCount) / 8 + 32
        case .deterministicRK2: return UInt64(workload.speciesCount) * 8 + UInt64(workload.reactionCount) * 4
        case .tauLeap: return UInt64(workload.reactionCount) * 8 + 64
        case .exactSSA: return UInt64(workload.reactionCount) * 4 + 64
        case .spatialSplit: return UInt64(workload.speciesCount) * 12 + UInt64(workload.reactionCount) * 4
        case .tissueCoupled: return UInt64(workload.speciesCount) * 16 + UInt64(workload.reactionCount) * 4
        }
    }

    private func fixedTableBytes(_ workload: VivoMechanismWorkload) throws -> UInt64 {
        let species = try checkedMultiply(UInt64(workload.speciesCount), 32, cohortID: workload.cohortID)
        let reactions = try checkedMultiply(UInt64(workload.reactionCount), 80, cohortID: workload.cohortID)
        let parameters = try checkedMultiply(UInt64(workload.parameterCount), 48, cohortID: workload.cohortID)
        return try checkedAdd(species, try checkedAdd(reactions, parameters, cohortID: workload.cohortID), cohortID: workload.cohortID)
    }

    private func recommendedMaximumStep(workload: VivoMechanismWorkload, mode: VivoAdaptiveFidelityMode) -> Double {
        let rateBound = workload.maximumPropensityPerSecond > 0 ? 0.1 / workload.maximumPropensityPerSecond : 60
        let eventBound = workload.expectedFiringsPerLanePerStep > 0 ? 1 / workload.expectedFiringsPerLanePerStep : 60
        switch mode {
        case .logic: return 60
        case .exactSSA: return min(60, max(1e-9, eventBound * policy.exactMaximumExpectedFirings))
        case .tauLeap: return min(60, max(1e-9, rateBound * 10))
        case .deterministicRK2: return min(60, max(1e-9, rateBound))
        case .spatialSplit, .tissueCoupled: return min(60, max(1e-9, rateBound * 0.5))
        }
    }

    private func checkedMultiply(_ lhs: UInt64, _ rhs: UInt64, cohortID: UInt32) throws -> UInt64 {
        let result = lhs.multipliedReportingOverflow(by: rhs)
        guard !result.overflow else { throw VivoAdaptiveFidelityError.arithmeticOverflow(cohortID) }
        return result.partialValue
    }

    private func checkedAdd(_ lhs: UInt64, _ rhs: UInt64, cohortID: UInt32) throws -> UInt64 {
        let result = lhs.addingReportingOverflow(rhs)
        guard !result.overflow else { throw VivoAdaptiveFidelityError.arithmeticOverflow(cohortID) }
        return result.partialValue
    }

    private static func fingerprint(_ plan: VivoAdaptiveFidelityPlan) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        let data = try encoder.encode(plan)
        return SHA256.hash(data: data).map { String(format: "%02x", $0) }.joined()
    }
}
