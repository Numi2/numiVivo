import Foundation

public extension VivoWorkflowTemplates {
    /// End-to-end algebraic software example. It neither samples a molecule nor
    /// qualifies a rate. Each native stage retains its physical/provenance checks.
    static func propertyRefinement() throws -> VivoWorkflowRecipe {
        let request = try VivoPropertyRefinementExamples.algebraicPath()
        let resources = VivoChemistryResourceContract(budget: request.budget,
            maximumInputBytes: request.budget.maximumBytes,maximumOutputBytes: request.budget.maximumBytes)
        let solver = VivoSelectedCISolverRequest(configuration: .init(maximumDeterminants: 3,
            selectionBatchSize: 1,minimumSelectionContributionHartree: 0,pt2ToleranceHartree: 1e-10,
            eigenResidualTolerance: 1e-12,maximumDavidsonSubspace: 8,
            maximumExternalDeterminants: 32,fullResidualTolerance: 1e-9))
        return try .init(identifier: "property-refinement",artifacts: [
            .init(identifier: "request",source: .json(kind: "vivo.property-directed-space-request",payload: VivoPlatformOperations.json(request)))
        ],nodes: [
            .init(identifier: "refine",operation: "vivo.native.property-directed-space",
                inputs: ["request":.artifact(identifier: "request")],resources: resources),
            .init(identifier: "anchor",operation: "vivo.native.refined-anchor-hamiltonian",
                inputs: ["refinement":.output(node: "refine",port: "refinement")],
                configuration: VivoPlatformOperations.json(VivoRefinedHamiltonianSelection(pointIdentifier: "barrier-point")),resources: resources),
            .init(identifier: "selected",operation: "vivo.native.selected-ci",version: "2",
                inputs: ["hamiltonian":.output(node: "anchor",port: "hamiltonian")],
                configuration: VivoPlatformOperations.json(solver),resources: resources),
            .init(identifier: "state",operation: "vivo.platform.selected-ci-state",
                inputs: ["electronic":.output(node: "selected",port: "electronic"),
                         "hamiltonian":.output(node: "anchor",port: "hamiltonian")],resources: resources),
            .init(identifier: "information",operation: "vivo.native.selective-orbital-information",
                inputs: ["state":.output(node: "state",port: "state")],
                configuration: VivoPlatformOperations.json(VivoOrbitalInformationSelection(orbitals: [0,1],pairs: [.init(0,1)])),resources: resources)
        ],outputs: [
            .init(name: "refinement",node: "refine",port: "refinement"),
            .init(name: "anchor",node: "anchor",port: "hamiltonian"),
            .init(name: "electronic",node: "selected",port: "electronic"),
            .init(name: "information",node: "information",port: "information")
        ])
    }
}
