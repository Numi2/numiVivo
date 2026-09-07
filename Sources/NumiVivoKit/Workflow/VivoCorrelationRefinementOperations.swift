import Foundation

/// Explicit distinction between an exploratory probe and an accepted eigenpair.
/// Neither policy asserts global-ground-root identity or chemical accuracy.
public struct VivoSelectedCISolverRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/selected-ci-solver-request/v1"
    public let schema: String
    public let configuration: VivoSelectedCIConfiguration
    public let initialDeterminants: [UInt64]?
    public let requireConverged: Bool
    public init(configuration: VivoSelectedCIConfiguration = .init(), initialDeterminants: [UInt64]? = nil,
                requireConverged: Bool = true) {
        schema = Self.schema; self.configuration = configuration
        self.initialDeterminants = initialDeterminants; self.requireConverged = requireConverged
    }
    public func validate(budget: VivoChemistryBudget) throws {
        guard schema == Self.schema else { throw VivoChemistryError.unsupported("selected-CI request schema") }
        try configuration.validate(budget: budget)
        if let seeds = initialDeterminants {
            guard !seeds.isEmpty, seeds.count <= configuration.maximumDeterminants,
                  Set(seeds).count == seeds.count else { throw VivoChemistryError.invalid("selected-CI request seed set") }
        }
    }
    func expectedSeeds(for h: VivoEmbeddedHamiltonian) -> [UInt64] {
        if let seeds = initialDeterminants { return seeds.sorted() }
        var reference: UInt64 = 0
        for p in 0..<h.alphaElectrons { reference |= UInt64(1) << (2*p) }
        for p in 0..<h.betaElectrons { reference |= UInt64(1) << (2*p+1) }
        return [reference]
    }
}

/// Operations share the existing content-addressed workflow and validation-on-
/// cache-read contract. A receipt authenticates an output, not a success claim.
public enum VivoCorrelationRefinementOperations {
    public static func selectedCI(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: "vivo.native.selected-ci",version: "2",implementationFingerprint: id,
            outputs: [.init(name: "electronic",kind: "vivo.selected-ci-result")],execute: { cfg,inputs,budget in
                guard Set(inputs.keys) == ["hamiltonian"], let data = inputs["hamiltonian"] else {
                    throw VivoChemistryError.invalid("selected-CI input slots")
                }
                let h = try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from: data)
                let request = try VivoCanonicalJSON.decode(VivoSelectedCISolverRequest.self,from: VivoCanonicalJSON.encode(cfg))
                try request.validate(budget: budget)
                let result = try VivoSelectedCI.solve(h,configuration: request.configuration,
                                                      initialDeterminants: request.initialDeterminants,budget: budget)
                if request.requireConverged && !result.converged {
                    throw VivoChemistryError.convergence("selected CI: \(result.termination.rawValue); no accepted eigenpair receipt")
                }
                return ["electronic": try VivoCanonicalJSON.encode(result)]
            },validateOutputs: { cfg,inputs,outputs,budget in
                guard Set(inputs.keys) == ["hamiltonian"], Set(outputs.keys) == ["electronic"],
                      let input = inputs["hamiltonian"], let output = outputs["electronic"] else {
                    throw VivoChemistryError.invalid("selected-CI output slots")
                }
                let request = try VivoCanonicalJSON.decode(VivoSelectedCISolverRequest.self,from: VivoCanonicalJSON.encode(cfg))
                let h = try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from: input)
                let result = try VivoCanonicalJSON.decode(VivoSelectedCIResult.self,from: output)
                try request.validate(budget: budget); try h.validate(budget: budget)
                guard h.orbitalCount <= 31, result.configuration == request.configuration,
                      result.seedDeterminants == request.expectedSeeds(for: h),
                      !request.requireConverged || result.converged else {
                    throw VivoChemistryError.invalid("selected-CI result/request or convergence-policy mismatch")
                }
                try VivoSelectedCI.validate(result,hamiltonian: h,budget: budget)
            })
    }
    public static func orbitalInformation(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: "vivo.native.selective-orbital-information",version: "1",implementationFingerprint: id,
            outputs: [.init(name: "information",kind: "vivo.selective-orbital-information")],execute: { cfg,inputs,budget in
                guard Set(inputs.keys) == ["state"], let data = inputs["state"] else {
                    throw VivoChemistryError.invalid("orbital-information input slots")
                }
                let state = try VivoCanonicalJSON.decode(VivoCIState.self,from: data)
                let selection = try VivoCanonicalJSON.decode(VivoOrbitalInformationSelection.self,from: VivoCanonicalJSON.encode(cfg))
                let result = try VivoSelectiveOrbitalInformation.analyze(state,selection: selection,budget: budget)
                return ["information": try VivoCanonicalJSON.encode(result)]
            },validateOutputs: { cfg,inputs,outputs,budget in
                guard Set(inputs.keys) == ["state"], Set(outputs.keys) == ["information"],
                      let input = inputs["state"], let output = outputs["information"] else {
                    throw VivoChemistryError.invalid("orbital-information output slots")
                }
                let state = try VivoCanonicalJSON.decode(VivoCIState.self,from: input)
                let selection = try VivoCanonicalJSON.decode(VivoOrbitalInformationSelection.self,from: VivoCanonicalJSON.encode(cfg))
                let result = try VivoCanonicalJSON.decode(VivoSelectiveOrbitalInformationResult.self,from: output)
                try VivoSelectiveOrbitalInformation.validate(result,state: state,selection: selection,budget: budget)
            })
    }
    /// Persists incomplete exploration as such. It never converts a budget limit
    /// or failed holdout into a sensitivity-established result.
    public static func propertyDirectedSpace(implementationFingerprint id: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier: "vivo.native.property-directed-space",version: "2",implementationFingerprint: id,
            outputs: [.init(name: "refinement",kind: "vivo.property-directed-space-result")],execute: { cfg,inputs,budget in
                guard cfg == .object([:]), Set(inputs.keys) == ["request"], let input = inputs["request"] else {
                    throw VivoChemistryError.invalid("property-refinement input/configuration slots")
                }
                let request = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceRequest.self,from: input)
                guard request.budget == budget else { throw VivoChemistryError.invalid("refinement task/request budgets differ") }
                return ["refinement": try VivoCanonicalJSON.encode(VivoPropertyDirectedSpace.run(request))]
            },validateOutputs: { cfg,inputs,outputs,budget in
                guard cfg == .object([:]), Set(inputs.keys) == ["request"], Set(outputs.keys) == ["refinement"],
                      let input = inputs["request"], let output = outputs["refinement"] else {
                    throw VivoChemistryError.invalid("property-refinement output/configuration slots")
                }
                let request = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceRequest.self,from: input)
                guard request.budget == budget else { throw VivoChemistryError.invalid("refinement task/request budgets differ") }
                let result = try VivoCanonicalJSON.decode(VivoPropertyDirectedSpaceResult.self,from: output)
                try VivoPropertyDirectedSpace.validate(result,request: request)
            })
    }
}
