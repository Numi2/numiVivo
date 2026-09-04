import Foundation

/// Creates an explicit executable authority plan from the actual reaction graph.
/// Unspecified components use exact SSA; computational cost never silently
/// promotes a low-count component to a deterministic approximation.
public enum VivoHybridExplicitPlanBuilder {
    public static func make(model: VivoExactSSAModel,
                            reactionAuthorities: [String: VivoHybridExecutionMode] = [:],
                            maximumStep: Double = 0.01) throws -> VivoHybridStochasticPlan {
        try model.validate()
        guard maximumStep.isFinite, maximumStep > 0 else {
            throw VivoHybridExecutionError.invalidPlan("maximum step must be positive and finite")
        }
        let identifiers = Set(model.reactions.map(\.id))
        guard reactionAuthorities.keys.allSatisfy({ identifiers.contains($0) }) else {
            throw VivoHybridExecutionError.invalidPlan("authority override names an absent reaction")
        }
        var parents = Array(model.species.indices)
        var ranks = [UInt8](repeating: 0, count: parents.count)
        func root(_ value: Int) -> Int {
            var result = value
            while parents[result] != result { result = parents[result] }
            return result
        }
        func unite(_ left: Int, _ right: Int) {
            let a = root(left), b = root(right)
            if a == b { return }
            if ranks[a] < ranks[b] { parents[a] = b }
            else {
                parents[b] = a
                if ranks[a] == ranks[b] { ranks[a] += 1 }
            }
        }
        var dependencies: [[UInt32]] = []
        for reaction in model.reactions {
            let members = Array(Set(reaction.changes.map(\.speciesIndex) +
                                    [reaction.reactantA, reaction.reactantB].compactMap { $0 })).sorted()
            guard let first = members.first else { throw VivoHybridExecutionError.invalidPlan("empty reaction") }
            for member in members.dropFirst() { unite(Int(first), Int(member)) }
            dependencies.append(members)
        }
        var groups: [Int: [UInt32]] = [:]
        for index in model.reactions.indices {
            groups[root(Int(dependencies[index][0])), default: []].append(UInt32(index))
        }
        let ordered = groups.values.sorted { $0[0] < $1[0] }
        var cohorts: [VivoHybridCohortPlan] = []
        var owned = Set<UInt32>()
        for (position, reactions) in ordered.enumerated() {
            let species = Array(Set(reactions.flatMap { dependencies[Int($0)] })).sorted()
            owned.formUnion(species)
            let requested = Set(reactions.compactMap { reactionAuthorities[model.reactions[Int($0)].id] })
            guard requested.count <= 1 else {
                throw VivoHybridExecutionError.invalidPlan("connected reactions have conflicting authorities; propensity dependencies cannot be split")
            }
            let mode = requested.first ?? .exactSSA
            guard mode != .spatialSplit else {
                throw VivoHybridExecutionError.invalidPlan("spatial authority needs an explicit transport participant")
            }
            cohorts.append(.init(cohortID: UInt32(position), mode: mode,
                                  reactionIndices: reactions, speciesIndices: species,
                                  recommendedMaximumStep: maximumStep,
                                  tauRelativeErrorTarget: nil, maximumExactEventsPerStep: nil,
                                  requiresStateMigration: false, previousMode: nil,
                                  rationale: [requested.isEmpty ? "explicit default: exact SSA" : "caller-selected component authority",
                                              "ownership includes all state changes and propensity inputs"]))
        }
        let unsigned = VivoHybridStochasticPlan(schemaVersion: 1, speciesCount: UInt32(model.species.count),
                                                cohorts: cohorts, uncoupledSpecies: model.species.indices.map(UInt32.init).filter { !owned.contains($0) },
                                                fingerprint: "")
        let fingerprint = try VivoCanonicalJSON.fingerprint(VivoCanonicalJSON.encode(unsigned)).hex
        let plan = VivoHybridStochasticPlan(schemaVersion: unsigned.schemaVersion, speciesCount: unsigned.speciesCount,
                                            cohorts: cohorts, uncoupledSpecies: unsigned.uncoupledSpecies, fingerprint: fingerprint)
        _ = try VivoHybridExecutionCompiler.compile(model: model, plan: plan)
        return plan
    }
}
