import Foundation

public enum VivoRefinementActionRole: String, Codable, Sendable { case discovery, confirmation }
public struct VivoRefinementCriterion: Codable, Sendable, Equatable {
    public var identifier: String
    public var metricUnit: String
    public var observableFingerprint: String
    public var maximumMetric: Double
    public var requiredConfirmations: Int
    public var requiresDisjointConfirmationSources: Bool
    public init(identifier: String,metricUnit: String,observableFingerprint: String,maximumMetric: Double,
                requiredConfirmations: Int = 1,requiresDisjointConfirmationSources: Bool = true) {
        self.identifier = identifier; self.metricUnit = metricUnit; self.observableFingerprint = observableFingerprint; self.maximumMetric = maximumMetric
        self.requiredConfirmations = requiredConfirmations; self.requiresDisjointConfirmationSources = requiresDisjointConfirmationSources
    }
}
public struct VivoRefinementActionProposal: Codable, Sendable, Equatable {
    public var identifier: String
    public var criterionIdentifier: String
    public var role: VivoRefinementActionRole
    public var prerequisites: [String]
    public var costClass: String
    public var declaredWorkUnits: Int
    public var initialEstimatedSeconds: Double
    /// A declared prioritization prior, never acceptance evidence.
    public var expectedMetricReduction: Double
    public init(identifier: String,criterionIdentifier: String,role: VivoRefinementActionRole,
                prerequisites: [String] = [],costClass: String,declaredWorkUnits: Int,
                initialEstimatedSeconds: Double,expectedMetricReduction: Double) {
        self.identifier = identifier; self.criterionIdentifier = criterionIdentifier; self.role = role
        self.prerequisites = prerequisites; self.costClass = costClass; self.declaredWorkUnits = declaredWorkUnits
        self.initialEstimatedSeconds = initialEstimatedSeconds; self.expectedMetricReduction = expectedMetricReduction
    }
}
/// Produced by a registered native evidence adapter, not by accepting an
/// arbitrary scalar in the CLI. The workflow binds/reconstructs its source.
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
    public init(metricUnit: String,observableFingerprint: String,value: Double?,nativeChecksPassed: Bool,evidenceIdentifier: String,
                executionSourceIdentifiers: [String],interpretation: String) {
        schema = Self.schema; self.metricUnit = metricUnit; self.observableFingerprint = observableFingerprint; self.value = value; self.nativeChecksPassed = nativeChecksPassed
        self.evidenceIdentifier = evidenceIdentifier; self.executionSourceIdentifiers = executionSourceIdentifiers; self.interpretation = interpretation
    }
}
public struct VivoRefinementActionRecord: Codable, Sendable, Equatable {
    public let actionIdentifier: String
    public let metric: VivoRefinementMetricEvidence?
    public let elapsedSeconds: Double
    /// Cached task validation is not a measurement of fresh execution cost.
    public let entirelyUncachedExecution: Bool
    public let failure: String?
}
public struct VivoRefinementPriority: Codable, Sendable, Equatable {
    public let actionIdentifier: String
    public let expectedSeconds: Double
    public let priority: Double
    public let costObservationCount: Int
}
public struct VivoRefinementCriterionState: Codable, Sendable, Equatable {
    public let criterionIdentifier: String
    public let latestDiscoveryMetric: Double?
    public let discoverySatisfied: Bool
    public let acceptedConfirmations: Int
    public let confirmationFailed: Bool
    public let satisfied: Bool
}
public struct VivoAdaptiveRefinementDecision: Codable, Sendable, Equatable {
    public let rankedActions: [VivoRefinementPriority]
    public let criteria: [VivoRefinementCriterionState]
    public let allCriteriaSatisfied: Bool
    public let failedConfirmation: Bool
}
public enum VivoAdaptiveRefinementPolicy {
    private static func name(_ text: String) -> Bool { !text.isEmpty && text.utf8.count <= 512 && !text.contains("\u{0}") }
    private static func hash(_ text: String) -> Bool { text.count == 64 && text.allSatisfy { "0123456789abcdef".contains($0) } }
    public static func validate(criteria: [VivoRefinementCriterion],actions: [VivoRefinementActionProposal]) throws {
        guard (1...32).contains(criteria.count), (1...128).contains(actions.count),
              Set(criteria.map(\.identifier)).count == criteria.count, Set(actions.map(\.identifier)).count == actions.count else {
            throw VivoChemistryError.invalid("adaptive campaign criterion/action capacity or duplicates")
        }
        let ids = Set(actions.map(\.identifier)), criterionIDs = Set(criteria.map(\.identifier))
        for criterion in criteria {
            guard name(criterion.identifier), name(criterion.metricUnit), hash(criterion.observableFingerprint), criterion.maximumMetric.isFinite, criterion.maximumMetric > 0,
                  (1...16).contains(criterion.requiredConfirmations), actions.contains(where: { $0.criterionIdentifier == criterion.identifier && $0.role == .discovery }),
                  actions.filter({ $0.criterionIdentifier == criterion.identifier && $0.role == .confirmation }).count >= criterion.requiredConfirmations else {
                throw VivoChemistryError.invalid("adaptive criterion has no valid discovery/confirmation contract")
            }
        }
        for action in actions {
            guard name(action.identifier), criterionIDs.contains(action.criterionIdentifier), name(action.costClass),
                  action.declaredWorkUnits > 0, action.declaredWorkUnits <= 1_000_000_000_000,
                  action.initialEstimatedSeconds.isFinite, action.initialEstimatedSeconds > 0,
                  action.expectedMetricReduction.isFinite, action.expectedMetricReduction >= 0,
                  Set(action.prerequisites).count == action.prerequisites.count,
                  Set(action.prerequisites).isSubset(of: ids), !action.prerequisites.contains(action.identifier) else {
                throw VivoChemistryError.invalid("adaptive action prerequisite, prior or work declaration")
            }
        }
        var remaining = actions, completed = Set<String>()
        while !remaining.isEmpty {
            let ready = remaining.filter { Set($0.prerequisites).isSubset(of: completed) }
            guard !ready.isEmpty else { throw VivoChemistryError.invalid("adaptive action dependency cycle") }
            completed.formUnion(ready.map(\.identifier)); remaining.removeAll { completed.contains($0.identifier) }
        }
    }
    public static func decide(criteria: [VivoRefinementCriterion],actions: [VivoRefinementActionProposal],
                              records: [VivoRefinementActionRecord],remainingDeclaredWork: Int) throws -> VivoAdaptiveRefinementDecision {
        try validate(criteria: criteria,actions: actions)
        let actionsByID = Dictionary(uniqueKeysWithValues: actions.map { ($0.identifier,$0) })
        guard remainingDeclaredWork >= 0, Set(records.map(\.actionIdentifier)).count == records.count else {
            throw VivoChemistryError.invalid("adaptive repeated execution record or remaining work")
        }
        var completed = Set<String>()
        for record in records {
            guard let action = actionsByID[record.actionIdentifier], Set(action.prerequisites).isSubset(of: completed),
                  record.elapsedSeconds.isFinite, record.elapsedSeconds >= 0,
                  (record.failure != nil) != (record.metric != nil) else {
                throw VivoChemistryError.invalid("adaptive execution order, metric/failure or timing")
            }
            if let metric = record.metric {
                let criterion = criteria.first { $0.identifier == action.criterionIdentifier }!
                guard metric.schema == VivoRefinementMetricEvidence.schema, metric.metricUnit == criterion.metricUnit,
                      metric.observableFingerprint == criterion.observableFingerprint,
                      hash(metric.evidenceIdentifier), name(metric.interpretation),
                      metric.value == nil || (metric.value!.isFinite && metric.value! >= 0),
                      Set(metric.executionSourceIdentifiers).count == metric.executionSourceIdentifiers.count,
                      metric.executionSourceIdentifiers.allSatisfy(name) else {
                    throw VivoChemistryError.invalid("adaptive metric identity, units or values")
                }
            }
            completed.insert(record.actionIdentifier)
        }
        var states: [VivoRefinementCriterionState] = []
        for criterion in criteria {
            var latest: Double?, ready = false, confirmations = 0, failed = false
            var discoverySources = Set<String>(), confirmedSources = Set<String>(), evidence = Set<String>()
            for record in records where actionsByID[record.actionIdentifier]!.criterionIdentifier == criterion.identifier {
                let action = actionsByID[record.actionIdentifier]!, metric = record.metric
                let acceptable = metric?.nativeChecksPassed == true && (metric?.value).map { $0 <= criterion.maximumMetric } == true
                if action.role == .discovery {
                    guard confirmations == 0, !failed else { throw VivoChemistryError.invalid("confirmation evidence cannot be reused for discovery tuning") }
                    latest = metric?.value; ready = acceptable
                    if let metric { discoverySources.formUnion(metric.executionSourceIdentifiers) }
                } else {
                    guard ready else { throw VivoChemistryError.invalid("confirmation ran before discovery met its criterion") }
                    let sources = Set(metric?.executionSourceIdentifiers ?? [])
                    let independent = !criterion.requiresDisjointConfirmationSources ||
                        (!sources.isEmpty && sources.isDisjoint(with: discoverySources.union(confirmedSources)))
                    let unique = metric.map { evidence.insert($0.evidenceIdentifier).inserted } ?? false
                    if !acceptable || !independent || !unique { failed = true }
                    else { confirmations += 1; confirmedSources.formUnion(sources) }
                }
            }
            states.append(.init(criterionIdentifier: criterion.identifier,latestDiscoveryMetric: latest,discoverySatisfied: ready,
                acceptedConfirmations: confirmations,confirmationFailed: failed,satisfied: ready && !failed && confirmations >= criterion.requiredConfirmations))
        }
        var priorities: [VivoRefinementPriority] = []
        for action in actions where !completed.contains(action.identifier) && Set(action.prerequisites).isSubset(of: completed) && action.declaredWorkUnits <= remainingDeclaredWork {
            let state = states.first { $0.criterionIdentifier == action.criterionIdentifier }!, criterion = criteria.first { $0.identifier == action.criterionIdentifier }!
            if state.satisfied || state.confirmationFailed || (action.role == .confirmation && !state.discoverySatisfied) ||
                (action.role == .discovery && state.discoverySatisfied) { continue }
            let observations = records.filter { $0.failure == nil && $0.entirelyUncachedExecution && $0.elapsedSeconds > 0 && actionsByID[$0.actionIdentifier]!.costClass == action.costClass }
            let ratios = observations.map { $0.elapsedSeconds/Double(actionsByID[$0.actionIdentifier]!.declaredWorkUnits) }.sorted()
            let seconds = ratios.isEmpty ? action.initialEstimatedSeconds : ratios[ratios.count/2]*Double(action.declaredWorkUnits)
            let benefit = action.role == .confirmation ? 1 : max(action.expectedMetricReduction/criterion.maximumMetric,1e-12)
            let priority = benefit/max(seconds,1e-9)
            guard seconds.isFinite, priority.isFinite else { throw VivoChemistryError.invalid("adaptive priority arithmetic") }
            priorities.append(.init(actionIdentifier: action.identifier,expectedSeconds: seconds,priority: priority,costObservationCount: ratios.count))
        }
        priorities.sort { $0.priority == $1.priority ? $0.actionIdentifier < $1.actionIdentifier : $0.priority > $1.priority }
        return .init(rankedActions: priorities,criteria: states,allCriteriaSatisfied: states.allSatisfy(\.satisfied),failedConfirmation: states.contains(where: \.confirmationFailed))
    }
}
