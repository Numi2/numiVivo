import Foundation

/// Thin typed adapters to the authoritative correlation/refinement operations.
/// A completed exploratory node is not necessarily a qualified calculation.
/// Downstream anchor/state operations enforce their stronger input conditions.
public enum VivoPlatformRefinementOperations {
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        let native = VivoCorrelationRefinementOperations.self
        func definition(_ operation: VivoChemistryOperation, _ inputs: [String: String], _ summary: String,
                        _ configure: @escaping @Sendable (VivoJSONValue) throws -> Void) -> VivoWorkflowDefinition {
            .init(operation: operation,inputKinds: inputs,summary: summary,
                validationScope: "source-bound native numerical reconstruction; bounded probes and candidate-pool sensitivity are not chemical accuracy or kinetic qualification",
                validateConfiguration: configure)
        }
        var definitions: [VivoWorkflowDefinition] = [
            definition(native.selectedCI(implementationFingerprint: id),["hamiltonian":"vivo.embedded-hamiltonian"],
                       "Bounded selected CI with complete connected residuals and explicit acceptance/probe policy.") { cfg in
                // Reuse the same strict top-level schema dispatch as chemistry-solve.
                let dispatched = try VivoNativeSolverDispatch.decode(VivoCanonicalJSON.encode(cfg),implementationFingerprint: id)
                guard dispatched.operation.identifier == "vivo.native.selected-ci" else {
                    throw VivoChemistryError.invalid("selected-CI recipe requires the selected-CI solver schema")
                }
            },
            definition(native.orbitalInformation(implementationFingerprint: id),["state":"vivo.ci-state"],
                       "Requested orbital marginals, entropy and pair mutual information with explicit unmeasured coverage.") { cfg in
                _ = try VivoPlatformOperations.decode(VivoOrbitalInformationSelection.self,cfg)
            },
            definition(native.propertyDirectedSpace(implementationFingerprint: id),["request":"vivo.property-directed-space-request"],
                       "Reaction-balanced nested electronic-space exploration, candidate-union checks and held-out confirmation.",VivoPlatformOperations.empty),
            definition(native.refinedHamiltonian(implementationFingerprint: id),["refinement":"vivo.property-directed-space-result"],
                       "Reconstruct one declared anchor and bind its Hamiltonian to established candidate-pool sensitivity.") { cfg in
                let selection = try VivoPlatformOperations.decode(VivoRefinedHamiltonianSelection.self,cfg)
                guard !selection.pointIdentifier.isEmpty, selection.pointIdentifier.utf8.count <= 256 else {
                    throw VivoChemistryError.invalid("refined-anchor point identifier")
                }
            }
        ]
        definitions.append(VivoPlatformOperations.pure(identifier: "vivo.platform.selected-ci-state",id: id,
            inputs: ["electronic":"vivo.selected-ci-result","hamiltonian":"vivo.embedded-hamiltonian"],
            outputs: [.init(name: "state",kind: "vivo.ci-state")],
            summary: "Extract a residual-converged selected state only after validating its physical Hamiltonian binding.",
            configure: VivoPlatformOperations.empty,calculate: { _,inputs,budget in
                let result = try VivoPlatformOperations.input(VivoSelectedCIResult.self,"electronic",inputs)
                let h = try VivoPlatformOperations.input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs)
                try VivoSelectedCI.validate(result,hamiltonian: h,budget: budget)
                guard result.converged else { throw VivoChemistryError.convergence("unconverged probe cannot supply an accepted CI-state node") }
                return ["state": try VivoCanonicalJSON.encode(result.state)]
            }))
        return definitions
    }
}
