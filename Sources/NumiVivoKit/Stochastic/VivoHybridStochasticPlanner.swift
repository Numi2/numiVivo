import CryptoKit
import Foundation

public enum VivoHybridExecutionMode: String, Codable, CaseIterable, Sendable {
    case exactSSA
    case tauLeap
    case deterministicRK2
    case spatialSplit
}

public struct VivoHybridReactionDescriptor: Codable, Equatable, Sendable {
    public let reactionIndex: UInt32
    public let speciesIndices: [UInt32]
    public let critical: Bool
    public let spatial: Bool
    public let minimumReactantCount: Double
    public let medianReactantCount: Double
    public let expectedFiringsPerStep: Double
    public let maximumPropensityPerSecond: Double
    public let stiffnessRatio: Double

    public init(
        reactionIndex: UInt32,
        speciesIndices: [UInt32],
        critical: Bool,
        spatial: Bool = false,
        minimumReactantCount: Double,
        medianReactantCount: Double,
        expectedFiringsPerStep: Double,
        maximumPropensityPerSecond: Double,
        stiffnessRatio: Double = 1
    ) {
        self.reactionIndex = reactionIndex
        self.speciesIndices = speciesIndices
        self.critical = critical
        self.spatial = spatial
        self.minimumReactantCount = minimumReactantCount
        self.medianReactantCount = medianReactantCount
        self.expectedFiringsPerStep = expectedFiringsPerStep
        self.maximumPropensityPerSecond = maximumPropensityPerSecond
        self.stiffnessRatio = stiffnessRatio
    }

    public func validate(speciesCount: UInt32) throws {
        guard !speciesIndices.isEmpty,
              Set(speciesIndices).count == speciesIndices.count,
              speciesIndices.allSatisfy({ $0 < speciesCount }),
              minimumReactantCount.isFinite,
              medianReactantCount.isFinite,
              expectedFiringsPerStep.isFinite,
              maximumPropensityPerSecond.isFinite,
              stiffnessRatio.isFinite,
              minimumReactantCount >= 0,
              medianReactantCount >= minimumReactantCount,
              expectedFiringsPerStep >= 0,
              maximumPropensityPerSecond >= 0,
              stiffnessRatio >= 1 else {
            throw VivoHybridPlanningError.invalidReaction(reactionIndex)
        }
    }
}

public struct VivoHybridStochasticPolicy: Codable, Equatable, Sendable {
    public var exactEnterCount: Double
    public var exactExitCount: Double
    public var exactMaximumExpectedFirings: Double
    public var exactMaximumReactionsPerCohort: UInt32
    public var deterministicEnterCount: Double
    public var deterministicExitCount: Double
    public var deterministicMinimumExpectedFirings: Double
    public var tauRelativeErrorTarget: Double
    public var tauMaximumExpectedFirings: Double
    public var stiffnessStepFactor: Double
    public var maximumExactEventsPerStep: UInt32

    public init(
        exactEnterCount: Double = 48,
        exactExitCount: Double = 96,
        exactMaximumExpectedFirings: Double = 32,
        exactMaximumReactionsPerCohort: UInt32 = 2_048,
        deterministicEnterCount: Double = 25_000,
        deterministicExitCount: Double = 12_500,
        deterministicMinimumExpectedFirings: Double = 2_048,
        tauRelativeErrorTarget: Double = 0.03,
        tauMaximumExpectedFirings: Double = 16_384,
        stiffnessStepFactor: Double = 0.1,
        maximumExactEventsPerStep: UInt32 = 8_192
    ) {
        self.exactEnterCount = exactEnterCount
        self.exactExitCount = exactExitCount
        self.exactMaximumExpectedFirings = exactMaximumExpectedFirings
        self.exactMaximumReactionsPerCohort = exactMaximumReactionsPerCohort
        self.deterministicEnterCount = deterministicEnterCount
        self.deterministicExitCount = deterministicExitCount
        self.deterministicMinimumExpectedFirings = deterministicMinimumExpectedFirings
        self.tauRelativeErrorTarget = tauRelativeErrorTarget
        self.tauMaximumExpectedFirings = tauMaximumExpectedFirings
        self.stiffnessStepFactor = stiffnessStepFactor
        self.maximumExactEventsPerStep = maximumExactEventsPerStep
    }

    public func validate() throws {
        guard exactEnterCount >= 0,
              exactExitCount > exactEnterCount,
              exactMaximumExpectedFirings > 0,
              exactMaximumReactionsPerCohort > 0,
              deterministicEnterCount > deterministicExitCount,
              deterministicExitCount > exactExitCount,
              deterministicMinimumExpectedFirings > 0,
              tauRelativeErrorTarget > 0,
              tauRelativeErrorTarget < 1,
              tauMaximumExpectedFirings >= exactMaximumExpectedFirings,
              stiffnessStepFactor > 0,
              stiffnessStepFactor <= 1,
              maximumExactEventsPerStep > 0 else {
            throw VivoHybridPlanningError.invalidPolicy
        }
    }
}

public struct VivoHybridCohortPlan: Codable, Equatable, Sendable {
    public let cohortID: UInt32
    public let mode: VivoHybridExecutionMode
    public let reactionIndices: [UInt32]
    public let speciesIndices: [UInt32]
    public let recommendedMaximumStep: Double
    public let tauRelativeErrorTarget: Double?
    public let maximumExactEventsPerStep: UInt32?
    public let requiresStateMigration: Bool
    public let previousMode: VivoHybridExecutionMode?
    public let rationale: [String]
}

public struct VivoHybridStochasticPlan: Codable, Equatable, Sendable {
    public let schemaVersion: UInt32
    public let speciesCount: UInt32
    public let cohorts: [VivoHybridCohortPlan]
    public let uncoupledSpecies: [UInt32]
    public let fingerprint: String

    public func cohort(owning species: UInt32) -> VivoHybridCohortPlan? {
        cohorts.first(where: { $0.speciesIndices.contains(species) })
    }
}

public enum VivoHybridPlanningError: Error, LocalizedError, Sendable {
    case invalidPolicy
    case invalidSpeciesCount
    case invalidReaction(UInt32)
    case duplicateReaction(UInt32)
    case arithmeticOverflow
    case inconsistentPreviousPlan

    public var errorDescription: String? {
        switch self {
        case .invalidPolicy: return "hybrid stochastic policy is invalid"
        case .invalidSpeciesCount: return "hybrid stochastic species count must be positive"
        case .invalidReaction(let index): return "invalid hybrid reaction descriptor \(index)"
        case .duplicateReaction(let index): return "duplicate hybrid reaction index \(index)"
        case .arithmeticOverflow: return "hybrid stochastic planner arithmetic overflow"
        case .inconsistentPreviousPlan: return "previous hybrid plan is inconsistent with the current species table"
        }
    }
}

/// Partitions the reaction hypergraph into independent components before choosing
/// a numerical method. Reactions sharing any dynamic species are never assigned
/// to different authorities, preventing exact, leap and continuous solvers from
/// independently consuming the same molecular pool.
public struct VivoHybridStochasticPlanner: Sendable {
    public let policy: VivoHybridStochasticPolicy

    public init(policy: VivoHybridStochasticPolicy = .init()) throws {
        try policy.validate()
        self.policy = policy
    }

    public func plan(
        speciesCount: UInt32,
        reactions: [VivoHybridReactionDescriptor],
        previous: VivoHybridStochasticPlan? = nil
    ) throws -> VivoHybridStochasticPlan {
        guard speciesCount > 0 else { throw VivoHybridPlanningError.invalidSpeciesCount }
        if let previous, previous.speciesCount != speciesCount {
            throw VivoHybridPlanningError.inconsistentPreviousPlan
        }
        var seenReactions = Set<UInt32>()
        for reaction in reactions {
            try reaction.validate(speciesCount: speciesCount)
            guard seenReactions.insert(reaction.reactionIndex).inserted else {
                throw VivoHybridPlanningError.duplicateReaction(reaction.reactionIndex)
            }
        }

        var unionFind = UnionFind(count: reactions.count)
        var ownerBySpecies: [UInt32: Int] = [:]
        for (reactionPosition, reaction) in reactions.enumerated() {
            for species in reaction.speciesIndices {
                if let owner = ownerBySpecies[species] {
                    unionFind.union(reactionPosition, owner)
                } else {
                    ownerBySpecies[species] = reactionPosition
                }
            }
        }

        var componentMembers: [Int: [Int]] = [:]
        for position in reactions.indices {
            componentMembers[unionFind.find(position), default: []].append(position)
        }
        let previousBySpecies: [UInt32: VivoHybridExecutionMode] = previousModeMap(previous)

        let orderedComponents = componentMembers.values.sorted { left, right in
            let leftReaction = left.map { reactions[$0].reactionIndex }.min() ?? 0
            let rightReaction = right.map { reactions[$0].reactionIndex }.min() ?? 0
            return leftReaction < rightReaction
        }

        var cohortPlans: [VivoHybridCohortPlan] = []
        cohortPlans.reserveCapacity(orderedComponents.count)
        for (cohortPosition, memberPositions) in orderedComponents.enumerated() {
            guard cohortPosition <= Int(UInt32.max) else { throw VivoHybridPlanningError.arithmeticOverflow }
            let members = memberPositions.map { reactions[$0] }
            let componentSpecies = Array(Set(members.flatMap(\.speciesIndices))).sorted()
            let previousModes = Set(componentSpecies.compactMap { previousBySpecies[$0] })
            let previousMode = previousModes.count == 1 ? previousModes.first : nil
            let selection = selectMode(members: members, previousMode: previousMode)
            let maximumStep = recommendedMaximumStep(members: members, mode: selection.mode)
            cohortPlans.append(.init(
                cohortID: UInt32(cohortPosition),
                mode: selection.mode,
                reactionIndices: members.map(\.reactionIndex).sorted(),
                speciesIndices: componentSpecies,
                recommendedMaximumStep: maximumStep,
                tauRelativeErrorTarget: selection.mode == .tauLeap ? policy.tauRelativeErrorTarget : nil,
                maximumExactEventsPerStep: selection.mode == .exactSSA ? policy.maximumExactEventsPerStep : nil,
                requiresStateMigration: previousMode != nil && previousMode != selection.mode,
                previousMode: previousMode,
                rationale: selection.rationale
            ))
        }

        let ownedSpecies = Set(cohortPlans.flatMap(\.speciesIndices))
        let uncoupled = (0..<speciesCount).filter { !ownedSpecies.contains($0) }
        let unsigned = VivoHybridStochasticPlan(
            schemaVersion: 1,
            speciesCount: speciesCount,
            cohorts: cohortPlans,
            uncoupledSpecies: uncoupled,
            fingerprint: ""
        )
        return .init(
            schemaVersion: unsigned.schemaVersion,
            speciesCount: unsigned.speciesCount,
            cohorts: unsigned.cohorts,
            uncoupledSpecies: unsigned.uncoupledSpecies,
            fingerprint: try fingerprint(unsigned)
        )
    }

    private func selectMode(
        members: [VivoHybridReactionDescriptor],
        previousMode: VivoHybridExecutionMode?
    ) -> (mode: VivoHybridExecutionMode, rationale: [String]) {
        if members.contains(where: \.spatial) {
            return (.spatialSplit, ["reaction component contains a spatial transport dependency"])
        }

        let minimumCount = members.map(\.minimumReactantCount).min() ?? 0
        let medianCount = median(members.map(\.medianReactantCount))
        let expectedFirings = members.reduce(0) { $0 + $1.expectedFiringsPerStep }
        let hasCritical = members.contains(where: \.critical)
        let exactThreshold = previousMode == .exactSSA ? policy.exactExitCount : policy.exactEnterCount
        let exactAffordable = members.count <= Int(policy.exactMaximumReactionsPerCohort) &&
            expectedFirings <= policy.exactMaximumExpectedFirings
        if (hasCritical || minimumCount <= exactThreshold), exactAffordable {
            var reasons = ["component is inside exact-event cost bounds"]
            if hasCritical { reasons.append("component contains a critical reaction") }
            if minimumCount <= exactThreshold { reasons.append("a reactant count is below the exact threshold") }
            return (.exactSSA, reasons)
        }

        let deterministicThreshold = previousMode == .deterministicRK2
            ? policy.deterministicExitCount
            : policy.deterministicEnterCount
        if medianCount >= deterministicThreshold,
           expectedFirings >= policy.deterministicMinimumExpectedFirings,
           !hasCritical {
            return (
                .deterministicRK2,
                [
                    "component molecular populations are above the continuous-state threshold",
                    "event volume exceeds the deterministic-entry threshold"
                ]
            )
        }

        if expectedFirings <= policy.tauMaximumExpectedFirings {
            return (.tauLeap, ["component is leap-compatible but outside the exact-event cost envelope"])
        }

        return (
            .deterministicRK2,
            ["event volume exceeds the bounded tau-leap envelope"]
        )
    }

    private func recommendedMaximumStep(
        members: [VivoHybridReactionDescriptor],
        mode: VivoHybridExecutionMode
    ) -> Double {
        let maximumRate = members.map(\.maximumPropensityPerSecond).max() ?? 0
        let maximumStiffness = members.map(\.stiffnessRatio).max() ?? 1
        let rateBound = maximumRate > 0 ? policy.stiffnessStepFactor / maximumRate : 60
        let stiffnessBound = maximumStiffness > 1 ? 1 / maximumStiffness : 60
        switch mode {
        case .exactSSA: return min(60, max(1e-12, rateBound * Double(policy.maximumExactEventsPerStep)))
        case .tauLeap: return min(60, max(1e-12, rateBound / policy.tauRelativeErrorTarget))
        case .deterministicRK2: return min(60, max(1e-12, min(rateBound, stiffnessBound)))
        case .spatialSplit: return min(60, max(1e-12, min(rateBound, stiffnessBound) * 0.5))
        }
    }

    private func previousModeMap(_ plan: VivoHybridStochasticPlan?) -> [UInt32: VivoHybridExecutionMode] {
        guard let plan else { return [:] }
        var result: [UInt32: VivoHybridExecutionMode] = [:]
        for cohort in plan.cohorts {
            for species in cohort.speciesIndices { result[species] = cohort.mode }
        }
        return result
    }

    private func median(_ values: [Double]) -> Double {
        guard !values.isEmpty else { return 0 }
        let sorted = values.sorted()
        let middle = sorted.count / 2
        if sorted.count.isMultiple(of: 2) {
            return 0.5 * (sorted[middle - 1] + sorted[middle])
        }
        return sorted[middle]
    }

    private func fingerprint(_ plan: VivoHybridStochasticPlan) throws -> String {
        let encoder = JSONEncoder()
        encoder.outputFormatting = [.sortedKeys, .withoutEscapingSlashes]
        return SHA256.hash(data: try encoder.encode(plan)).map { String(format: "%02x", $0) }.joined()
    }
}

private struct UnionFind {
    private var parent: [Int]
    private var rank: [UInt8]

    init(count: Int) {
        parent = Array(0..<count)
        rank = [UInt8](repeating: 0, count: count)
    }

    mutating func find(_ element: Int) -> Int {
        if parent[element] != element {
            parent[element] = find(parent[element])
        }
        return parent[element]
    }

    mutating func union(_ left: Int, _ right: Int) {
        let leftRoot = find(left)
        let rightRoot = find(right)
        guard leftRoot != rightRoot else { return }
        if rank[leftRoot] < rank[rightRoot] {
            parent[leftRoot] = rightRoot
        } else if rank[leftRoot] > rank[rightRoot] {
            parent[rightRoot] = leftRoot
        } else {
            parent[rightRoot] = leftRoot
            rank[leftRoot] &+= 1
        }
    }
}
