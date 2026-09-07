import Foundation

public enum VivoAdaptiveCampaignExamples {
    /// Independent draws from an exactly enumerable four-state distribution.
    /// This is a statistical software fixture, not a molecular trajectory.
    public static func equilibriumRequest(seed: UInt64,samplesPerOrigin: Int = 256) -> VivoMBARTargetRefinementRequest {
        let energy = [0.0,log(2.0),log(3.0),log(4.0)]
        var rng = VivoSplitMix64(state: seed), samples: [VivoMBARTargetSample] = []
        for origin in 0..<2 {
            let weights = energy.map { origin == 0 ? 1.0 : exp(-$0) }, total = weights.reduce(0,+)
            for sample in 0..<samplesPerOrigin {
                var selector = rng.unit()*total, microstate = weights.count-1
                for i in weights.indices { selector -= weights[i]; if selector <= 0 { microstate = i; break } }
                let identity = "finite-state-seed-\(seed)-origin-\(origin)-draw-\(sample)"
                samples.append(.init(identifier: identity,originStateIndex: origin,independentBlockIdentifier: identity,
                    coordinate: microstate < 2 ? 0 : 1,sampledReducedPotentials: [0,energy[microstate]],targetReducedPotentials: [0,energy[microstate]]))
            }
        }
        var cfg = VivoMBARTargetConfiguration(bootstrapReplicates: 32,minimumBinEffectiveSamples: 10,
                                              maximumBinWeight: 0.1,maximumPrimitiveElements: 50_000_000)
        cfg.bootstrapSeed = seed ^ 0xB0057
        return .init(sampledStateIdentifiers: ["finite-state-uniform","finite-state-inverse-multiplicity"],
            targetIdentifiers: ["reference","refined"],matchingSampledStateIndices: [0,1],binEdges: [-0.5,0.5,1.5],
            coordinateUnit: "dimensionless finite-state partition",samples: samples,configuration: cfg)
    }
    public static func make(_ name: String) throws -> VivoAdaptiveCampaignRequest {
        switch name {
        case "equilibrium-correction":
            let discovery = equilibriumRequest(seed: 0xC01D), confirmation = equilibriumRequest(seed: 0xC02D)
            let context = try VivoAdaptiveMetricContext.targetProfile(discovery).hex
            let criterion = VivoRefinementCriterion(identifier: "finite-state-profile",metricUnit: "dimensionless-bin-free-energy-standard-deviation",
                observableFingerprint: context,maximumMetric: 0.15)
            func action(_ identifier: String,_ role: VivoRefinementActionRole,_ input: VivoMBARTargetRefinementRequest,
                        prerequisites: [String]) throws -> VivoAdaptiveCampaignAction {
                let budget = VivoChemistryBudget(maximumOperatorApplications: input.configuration.maximumPrimitiveElements)
                let recipe = VivoWorkflowRecipe(identifier: identifier,artifacts: [
                    .init(identifier: "request",source: .json(kind: "vivo.mbar-target-refinement-request",payload: try VivoPlatformOperations.json(input)))],
                    nodes: [.init(identifier: "correct",operation: "vivo.platform.mbar-target-refinement",inputs: ["request":.artifact(identifier: "request")],
                        resources: .init(budget: budget))],outputs: [.init(name: "metric",node: "correct",port: "metric"),.init(name: "result",node: "correct",port: "result")])
                return .init(proposal: .init(identifier: identifier,criterionIdentifier: criterion.identifier,role: role,prerequisites: prerequisites,
                    costClass: "native-target-mbar-fixture",declaredWorkUnits: 50_000_000,initialEstimatedSeconds: 1,expectedMetricReduction: 0.1),recipe: recipe)
            }
            return try .init(identifier: "finite-state-equilibrium-refinement-regression",criteria: [criterion],actions: [
                action("discovery",.discovery,discovery,prerequisites: []),action("confirmation",.confirmation,confirmation,prerequisites: ["discovery"])],
                maximumActions: 2,maximumDeclaredWorkUnits: 100_000_000)
        case "electronic-crosscheck":
            let direct = try VivoPropertyRefinementExamples.algebraicPath(), selected = try VivoPropertyRefinementExamples.algebraicPath(selectedCI: true)
            let criterion = try VivoRefinementCriterion(identifier: "algebraic-electronic-profile",metricUnit: "normalized-electronic-profile-sensitivity",
                observableFingerprint: VivoAdaptiveMetricContext.electronic(request: direct).hex,maximumMetric: 1,
                requiresDisjointConfirmationSources: false)
            func action(_ identifier: String,_ role: VivoRefinementActionRole,_ input: VivoPropertyDirectedSpaceRequest,
                        prerequisites: [String]) throws -> VivoAdaptiveCampaignAction {
                let recipe = VivoWorkflowRecipe(identifier: identifier,artifacts: [
                    .init(identifier: "request",source: .json(kind: "vivo.property-directed-space-request",payload: try VivoPlatformOperations.json(input)))],
                    nodes: [.init(identifier: "refine",operation: "vivo.native.property-directed-space",version: "2",inputs: ["request":.artifact(identifier: "request")],resources: .init(budget: input.budget)),
                        .init(identifier: "metric",operation: "vivo.platform.metric-electronic-refinement",inputs: ["result":.output(node: "refine",port: "refinement")],resources: .init(budget: input.budget))],
                    outputs: [.init(name: "metric",node: "metric",port: "metric"),.init(name: "result",node: "refine",port: "refinement")])
                return .init(proposal: .init(identifier: identifier,criterionIdentifier: criterion.identifier,role: role,prerequisites: prerequisites,
                    costClass: identifier,declaredWorkUnits: input.budget.maximumOperatorApplications,initialEstimatedSeconds: 1,expectedMetricReduction: 1),recipe: recipe)
            }
            return try .init(identifier: "algebraic-two-solver-crosscheck-not-independent-geometries",criteria: [criterion],actions: [
                action("direct",.discovery,direct,prerequisites: []),action("selected",.confirmation,selected,prerequisites: ["direct"])],
                maximumActions: 2,maximumDeclaredWorkUnits: 2*direct.budget.maximumOperatorApplications)
        default: throw VivoChemistryError.invalid("campaign template: equilibrium-correction or electronic-crosscheck")
        }
    }
}
