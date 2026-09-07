import Foundation

public struct VivoRefinementCriterion: Codable, Sendable, Equatable {
    public var identifier: String
    public var metricUnit: String
    public var observableFingerprint: String
    public var maximumMetric: Double
    public var requiredConfirmations: Int
    public var requiresDisjointConfirmationSources: Bool
    public init(identifier: String, metricUnit: String, observableFingerprint: String, maximumMetric: Double, requiredConfirmations: Int = 1,
                requiresDisjointConfirmationSources: Bool = true) {
        self.identifier = identifier; self.metricUnit = metricUnit; self.observableFingerprint = observableFingerprint; self.maximumMetric = maximumMetric
        self.requiredConfirmations = requiredConfirmations; self.requiresDisjointConfirmationSources = requiresDisjointConfirmationSources
    }
}
public enum VivoRefinementActionRole: String, Codable, Sendable { case discovery, confirmation }
public struct VivoRefinementActionProposal: Codable, Sendable, Equatable {
    public var identifier: String
    public var criterionIdentifier: String
    public var role: VivoRefinementActionRole
    public var prerequisites: [String]
    public var costClass: String
    public var declaredWorkUnits: Int
    public var initialEstimatedSeconds: Double
    /// A declared ranking prior, not observed uncertainty reduction or a bound.
    public var expectedMetricReduction: Double
    public init(identifier: String, criterionIdentifier: String, role: VivoRefinementActionRole,
                prerequisites: [String] = [], costClass: String, declaredWorkUnits: Int,
                initialEstimatedSeconds: Double, expectedMetricReduction: Double) {
        self.identifier = identifier; self.criterionIdentifier = criterionIdentifier; self.role = role
        self.prerequisites = prerequisites; self.costClass = costClass; self.declaredWorkUnits = declaredWorkUnits
        self.initialEstimatedSeconds = initialEstimatedSeconds; self.expectedMetricReduction = expectedMetricReduction
    }
}
/// Created by native evidence adapters, never accepted merely because a user
/// supplies a small scalar. Workflow readers bind this to the producing task.
public struct VivoRefinementMetricEvidence: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/refinement-metric-evidence/v1"
    public let schema: String
    public let metricUnit: String
    public let observableFingerprint: String
    public let value: Double?
    public let nativeChecksPassed: Bool
    public let evidenceIdentifier: String
    public let executionSourceIdentifiers: [String]
    public let interpretation: String
    public init(metricUnit: String, observableFingerprint: String, value: Double?, nativeChecksPassed: Bool, evidenceIdentifier: String,
                executionSourceIdentifiers: [String], interpretation: String) {
        schema = Self.schema; self.metricUnit = metricUnit; self.observableFingerprint = observableFingerprint; self.value = value
        self.nativeChecksPassed = nativeChecksPassed; self.evidenceIdentifier = evidenceIdentifier
        self.executionSourceIdentifiers = executionSourceIdentifiers; self.interpretation = interpretation
    }
    public func validate() throws {
        guard schema == Self.schema, !metricUnit.isEmpty, !interpretation.isEmpty,
              observableFingerprint.count == 64, observableFingerprint.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
              value == nil || (value!.isFinite && value! >= 0), !nativeChecksPassed || value != nil,
              evidenceIdentifier.count == 64, evidenceIdentifier.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
              Set(executionSourceIdentifiers).count == executionSourceIdentifiers.count,
              executionSourceIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 512 }) else {
            throw VivoChemistryError.invalid("refinement metric evidence shape")
        }
    }
}
public struct VivoRefinementActionRecord: Codable, Sendable, Equatable {
    public let actionIdentifier: String
    public let metric: VivoRefinementMetricEvidence?
    public let elapsedSeconds: Double
    /// Cached/replay-heavy execution must not calibrate a fresh calculation.
    public let entirelyUncachedExecution: Bool
    public let failure: String?
}
public struct VivoRefinementActionScore: Codable, Sendable, Equatable {
    public let actionIdentifier: String
    public let estimatedSeconds: Double
    public let usedMeasuredCostClass: Bool
    public let expectedNormalizedBenefit: Double
    public let priority: Double
}
public struct VivoRefinementCriterionState: Codable, Sendable, Equatable {
    public let criterionIdentifier: String
    public let qualifyingConfirmations: [String]
    public let confirmationFailed: Bool
    public let satisfied: Bool
}
public struct VivoRefinementDecision: Codable, Sendable, Equatable {
    public let criteria: [VivoRefinementCriterionState]
    public let rankedActions: [VivoRefinementActionScore]
    public var selectedActionIdentifier: String? { rankedActions.first?.actionIdentifier }
    public var allCriteriaSatisfied: Bool { !criteria.isEmpty && criteria.allSatisfy(\.satisfied) }
}

/// Allocation policy only. The workflow executor supplies reconstructed native
/// observations. Timings affect order, not acceptance tolerances or evidence.
public enum VivoAdaptiveRefinementPolicy {
    public static func validate(criteria: [VivoRefinementCriterion], proposals: [VivoRefinementActionProposal]) throws {
        guard (1...32).contains(criteria.count), (1...128).contains(proposals.count),
              Set(criteria.map(\.identifier)).count == criteria.count,
              Set(proposals.map(\.identifier)).count == proposals.count else { throw VivoChemistryError.invalid("adaptive criterion/action identities") }
        let criteriaIDs = Set(criteria.map(\.identifier)), actionIDs = Set(proposals.map(\.identifier))
        for c in criteria {
            guard !c.identifier.isEmpty, c.identifier.utf8.count <= 256, !c.metricUnit.isEmpty,
                  c.observableFingerprint.count == 64, c.observableFingerprint.utf8.allSatisfy({ (48...57).contains($0) || (97...102).contains($0) }),
                  c.maximumMetric.isFinite, c.maximumMetric > 0, (1...16).contains(c.requiredConfirmations),
                  proposals.filter({ $0.criterionIdentifier == c.identifier && $0.role == .confirmation }).count >= c.requiredConfirmations else {
                throw VivoChemistryError.invalid("adaptive criterion needs an explicit tolerance and enough confirmation actions")
            }
        }
        for a in proposals {
            guard !a.identifier.isEmpty, a.identifier.utf8.count <= 256, criteriaIDs.contains(a.criterionIdentifier),
                  !a.costClass.isEmpty, a.costClass.utf8.count <= 256, a.declaredWorkUnits > 0,
                  a.declaredWorkUnits <= 1_000_000_000_000,
                  a.initialEstimatedSeconds.isFinite, a.initialEstimatedSeconds > 0,
                  a.expectedMetricReduction.isFinite, a.expectedMetricReduction >= 0,
                  Set(a.prerequisites).count == a.prerequisites.count,
                  Set(a.prerequisites).isSubset(of: actionIDs), !a.prerequisites.contains(a.identifier) else {
                throw VivoChemistryError.invalid("adaptive action costs, benefit prior or prerequisites")
            }
        }
        var visited = Set<String>()
        while visited.count < proposals.count {
            let ready = proposals.filter { !visited.contains($0.identifier) && Set($0.prerequisites).isSubset(of: visited) }
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("adaptive action dependency cycle") }
            visited.formUnion(ready.map(\.identifier))
        }
    }
    public static func decide(criteria: [VivoRefinementCriterion], proposals: [VivoRefinementActionProposal],
                              records: [VivoRefinementActionRecord], remainingDeclaredWorkUnits: Int) throws -> VivoRefinementDecision {
        try validate(criteria: criteria,proposals: proposals)
        let actions = Dictionary(uniqueKeysWithValues: proposals.map { ($0.identifier,$0) })
        guard Set(records.map(\.actionIdentifier)).count == records.count, remainingDeclaredWorkUnits >= 0 else {
            throw VivoChemistryError.invalid("adaptive execution ledger")
        }
        var priorSuccess = Set<String>(), confirmationStarted = Set<String>()
        var priorDiscovery: [String: VivoRefinementMetricEvidence] = [:]
        for record in records {
            guard let action = actions[record.actionIdentifier], record.elapsedSeconds.isFinite, record.elapsedSeconds >= 0,
                  Set(action.prerequisites).isSubset(of: priorSuccess),
                  (record.failure == nil) == (record.metric != nil) else { throw VivoChemistryError.invalid("adaptive observation ledger or execution order") }
            try record.metric?.validate()
            let criterion = criteria.first { $0.identifier == action.criterionIdentifier }!
            if action.role == .discovery {
                guard !confirmationStarted.contains(criterion.identifier) else {
                    throw VivoChemistryError.invalid("confirmation evidence cannot be reused for discovery tuning")
                }
                priorDiscovery[criterion.identifier] = record.metric
            } else {
                let hasDiscovery = proposals.contains { $0.criterionIdentifier == criterion.identifier && $0.role == .discovery }
                let previous = priorDiscovery[criterion.identifier]
                guard !hasDiscovery || (previous?.nativeChecksPassed == true && (previous?.value ?? .infinity) <= criterion.maximumMetric) else {
                    throw VivoChemistryError.invalid("confirmation ran before discovery met its criterion")
                }
                confirmationStarted.insert(criterion.identifier)
            }
            if record.failure == nil { priorSuccess.insert(record.actionIdentifier) }
        }
        let done = Set(records.map(\.actionIdentifier)), succeeded = Set(records.filter { $0.failure == nil }.map(\.actionIdentifier))
        var states: [VivoRefinementCriterionState] = [], scores: [VivoRefinementActionScore] = []
        for criterion in criteria {
            let related = records.filter { actions[$0.actionIdentifier]!.criterionIdentifier == criterion.identifier }
            for record in related {
                if let metric = record.metric, metric.metricUnit != criterion.metricUnit || metric.observableFingerprint != criterion.observableFingerprint {
                    throw VivoChemistryError.invalid("adaptive observations mix incompatible metric units")
                }
            }
            let discovery = related.filter { actions[$0.actionIdentifier]!.role == .discovery }
            let confirmations = related.filter { actions[$0.actionIdentifier]!.role == .confirmation }
            var sources = Set(discovery.flatMap { $0.metric?.executionSourceIdentifiers ?? [] }), evidence = Set(discovery.compactMap { $0.metric?.evidenceIdentifier })
            var qualified: [String] = [], failed = false
            for record in confirmations {
                guard let metric = record.metric else { failed = true; continue }
                let independent = !criterion.requiresDisjointConfirmationSources ||
                    (!metric.executionSourceIdentifiers.isEmpty && sources.isDisjoint(with: metric.executionSourceIdentifiers))
                let unique = evidence.insert(metric.evidenceIdentifier).inserted
                if metric.nativeChecksPassed, let value = metric.value, value <= criterion.maximumMetric, independent, unique {
                    qualified.append(record.actionIdentifier)
                } else { failed = true }
                sources.formUnion(metric.executionSourceIdentifiers)
            }
            let satisfied = !failed && qualified.count >= criterion.requiredConfirmations
            states.append(.init(criterionIdentifier: criterion.identifier,qualifyingConfirmations: qualified,
                                confirmationFailed: failed,satisfied: satisfied))
            if satisfied || failed { continue }
            let hasDiscovery = proposals.contains { $0.criterionIdentifier == criterion.identifier && $0.role == .discovery }
            let last = discovery.last?.metric
            let readyToConfirm = !hasDiscovery || (last?.nativeChecksPassed == true && (last?.value ?? .infinity) <= criterion.maximumMetric)
            for action in proposals where action.criterionIdentifier == criterion.identifier && !done.contains(action.identifier) {
                guard action.declaredWorkUnits <= remainingDeclaredWorkUnits,
                      Set(action.prerequisites).isSubset(of: succeeded),
                      readyToConfirm ? action.role == .confirmation : action.role == .discovery else { continue }
                let timings = records.filter { $0.failure == nil && $0.entirelyUncachedExecution && $0.elapsedSeconds > 0 &&
                    actions[$0.actionIdentifier]!.costClass == action.costClass }
                let measuredUnits = timings.reduce(0.0) { $0+Double(actions[$1.actionIdentifier]!.declaredWorkUnits) }
                let seconds = measuredUnits > 0
                    ? timings.reduce(0.0) { $0+$1.elapsedSeconds }/measuredUnits*Double(action.declaredWorkUnits)
                    : action.initialEstimatedSeconds
                let benefit = action.role == .confirmation ? 1.0 : max(1e-12,action.expectedMetricReduction/criterion.maximumMetric)
                let priority = benefit/max(1e-9,seconds)
                guard seconds.isFinite, seconds > 0, benefit.isFinite, priority.isFinite else {
                    throw VivoChemistryError.invalid("adaptive action score overflow")
                }
                scores.append(.init(actionIdentifier: action.identifier,estimatedSeconds: seconds,
                    usedMeasuredCostClass: measuredUnits > 0,expectedNormalizedBenefit: benefit,priority: priority))
            }
        }
        scores.sort { $0.priority == $1.priority ? $0.actionIdentifier < $1.actionIdentifier : $0.priority > $1.priority }
        return .init(criteria: states,rankedActions: scores)
    }
}
