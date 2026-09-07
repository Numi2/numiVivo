import Foundation

public struct VivoReactiveTrainingRequest: Codable, Sendable, Equatable {
    public let authority: VivoNuclearPotentialDefinition
    public let baseline: VivoNuclearPotentialDefinition
    public let labels: [VivoReactiveTrainingLabel]
    public let heldOutGroups: [String]
    public let configuration: VivoReactiveSurrogateConfiguration
    public init(authority: VivoNuclearPotentialDefinition, baseline: VivoNuclearPotentialDefinition,
                labels: [VivoReactiveTrainingLabel], heldOutGroups: [String], configuration: VivoReactiveSurrogateConfiguration) {
        self.authority = authority; self.baseline = baseline; self.labels = labels
        self.heldOutGroups = heldOutGroups; self.configuration = configuration
    }
}
public struct VivoRingPolymerWorkflowRequest: Codable, Sendable, Equatable {
    public let potential: VivoNuclearMetalSpecification
    public let checkpoint: VivoRingPolymerCheckpoint
    public let sweeps: Int
    public init(potential: VivoNuclearMetalSpecification, checkpoint: VivoRingPolymerCheckpoint, sweeps: Int) {
        self.potential = potential; self.checkpoint = checkpoint; self.sweeps = sweeps
    }
    func admit(_ budget: VivoChemistryBudget) throws {
        let reserved = try potential.reservationBytes(budget: budget)
        try checkpoint.validate()
        let c = checkpoint.configuration
        guard checkpoint.definition == (try potential.definition()), (1...100_000).contains(sweeps),
              c.maximumPrimitiveWork <= budget.maximumOperatorApplications,
              Double(reserved)+Double(sweeps+16)*Double(c.beadCount)*Double(checkpoint.definition.atomIndices.count)*32 <= Double(budget.maximumBytes) else {
            throw VivoChemistryError.resourceLimit("ring workflow identity, stored coordinates or numerical budget")
        }
    }
}
public struct VivoReactiveSamplingWorkflowRequest: Codable, Sendable, Equatable {
    public let authority: VivoNuclearMetalSpecification
    public let baseline: VivoNuclearMetalSpecification
    public let model: VivoReactiveSurrogateModel
    public let checkpoint: VivoReactiveSamplingCheckpoint
    public let sweeps: Int
    public let metalInference: Bool
    public init(authority: VivoNuclearMetalSpecification, baseline: VivoNuclearMetalSpecification,
                model: VivoReactiveSurrogateModel, checkpoint: VivoReactiveSamplingCheckpoint, sweeps: Int, metalInference: Bool = true) {
        self.authority = authority; self.baseline = baseline; self.model = model; self.checkpoint = checkpoint
        self.sweeps = sweeps; self.metalInference = metalInference
    }
    func admit(_ budget: VivoChemistryBudget) throws {
        try checkpoint.validate(model: model)
        let a = try authority.reservationBytes(budget: budget), b = try baseline.reservationBytes(budget: budget)
        let backend = metalInference ? "metal-fp32-rbf-analytic-derivative/v1" : "native-fp64-rbf/v1"
        guard model.payload.authorityDefinition == (try authority.definition()),
              model.payload.baselineDefinition == (try baseline.definition()), (1...100_000).contains(sweeps),
              checkpoint.backendProfile == backend, model.payload.qualification.passed,
              Double(a)+Double(b)+Double(checkpoint.configuration.maximumStoredCoordinateElements)*8 <= Double(budget.maximumBytes),
              Double(sweeps)*Double(checkpoint.configuration.integrationSteps+1) <= Double(budget.maximumOperatorApplications) else {
            throw VivoChemistryError.resourceLimit("reactive workflow model, profile, storage or propagation budget")
        }
    }
}
public struct VivoReactiveLabelGeometry: Codable, Sendable, Equatable {
    public let identifier: String
    public let sourceGroup: String
    public let positionsNM: [VivoVector3D]
    public init(identifier: String, sourceGroup: String, positionsNM: [VivoVector3D]) {
        self.identifier = identifier; self.sourceGroup = sourceGroup; self.positionsNM = positionsNM
    }
}
public struct VivoReactiveLabelWorkflowRequest: Codable, Sendable, Equatable {
    public let authority: VivoNuclearMetalSpecification
    public let baseline: VivoNuclearMetalSpecification
    public let geometries: [VivoReactiveLabelGeometry]
    public init(authority: VivoNuclearMetalSpecification, baseline: VivoNuclearMetalSpecification, geometries: [VivoReactiveLabelGeometry]) {
        self.authority = authority; self.baseline = baseline; self.geometries = geometries
    }
    func admit(_ budget: VivoChemistryBudget) throws {
        let a = try authority.reservationBytes(budget: budget), b = try baseline.reservationBytes(budget: budget)
        let x = try authority.definition(), y = try baseline.definition()
        guard (1...4096).contains(geometries.count), Set(geometries.map(\.identifier)).count == geometries.count,
              x.atomIndices == y.atomIndices, x.particleIndices == y.particleIndices, x.massesDa == y.massesDa, x.periodicCell == y.periodicCell,
              geometries.allSatisfy({ !$0.identifier.isEmpty && !$0.sourceGroup.isEmpty && $0.positionsNM.count == x.atomIndices.count && $0.positionsNM.allSatisfy(\.isFinite) }),
              Double(a)+Double(b)+Double(geometries.count)*Double(x.atomIndices.count)*256 <= Double(budget.maximumBytes),
              2*geometries.count <= budget.maximumOperatorApplications else {
            throw VivoChemistryError.resourceLimit("reactive authority label mapping, source groups or admission")
        }
    }
}

public enum VivoPlatformPrecisionSamplingOperations {
    private static let scope = "Native complete-Hamiltonian equilibrium sampling with fixed source/precision/model identities; cached samples validate checkpoint, RNG/energy bookkeeping and finite values, not a new GPU/electronic replay. Not a physical trajectory, independent-sample claim or chemical qualification."
    public static func definitions(implementationFingerprint id: VivoFingerprint) -> [VivoWorkflowDefinition] {
        let pure = [
            VivoPlatformOperations.pure(identifier: "vivo.platform.barrier-tunnelling",id: id,
                inputs: ["request":"vivo.barrier-tunnelling-request"],outputs: [.init(name: "correction",kind: "vivo.barrier-tunnelling-result")],
                summary: "Wigner small-correction or converged one-dimensional Eckart tunnelling, separate from classical kappa.",configure: VivoPlatformOperations.empty,
                calculate: { _,input,budget in
                    let r = try VivoPlatformOperations.input(VivoBarrierTunnellingRequest.self,"request",input)
                    guard r.maximumEvaluations <= budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("tunnelling workflow budget") }
                    return ["correction":try VivoCanonicalJSON.encode(VivoBarrierTunnelling.calculate(r))]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.global-kinetic-uncertainty",id: id,
                inputs: ["request":"vivo.global-kinetic-uncertainty-request"],outputs: [.init(name: "prediction",kind: "vivo.global-joint-uncertainty-result")],
                summary: "Full nonlinear kinetic propagation of supplied joint draws and explicit model probabilities, preserving unresolved sources.",configure: VivoPlatformOperations.empty,
                calculate: { _,input,budget in
                    let r = try VivoPlatformOperations.input(VivoGlobalKineticUncertaintyRequest.self,"request",input)
                    guard r.ensemble.maximumPrimitiveWork <= budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("joint kinetic workflow budget") }
                    _ = try budget.elements([r.ensemble.draws.count,max(1,r.parameters.count)],simultaneousArrays: 4)
                    return ["prediction":try VivoCanonicalJSON.encode(VivoGlobalKineticUncertainty.calculate(r))]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.reactive-surrogate-train",id: id,
                inputs: ["request":"vivo.reactive-training-request"],outputs: [.init(name: "model",kind: "vivo.reactive-surrogate-model")],
                summary: "Deterministic energy-and-force delta fit with strictly group-separated held-out validation.",configure: VivoPlatformOperations.empty,
                calculate: { _,input,budget in
                    let r = try VivoPlatformOperations.input(VivoReactiveTrainingRequest.self,"request",input)
                    guard r.configuration.maximumPrimitiveWork <= budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("reactive training workflow budget") }
                    _ = try budget.elements([r.labels.count,max(1,r.authority.atomIndices.count),16],simultaneousArrays: 4)
                    return ["model":try VivoCanonicalJSON.encode(VivoReactiveDeltaSurrogate.train(authority: r.authority,baseline: r.baseline,
                        labels: r.labels,heldOutGroups: r.heldOutGroups,configuration: r.configuration))]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.ring-polymer-convergence",id: id,
                inputs: ["request":"vivo.ring-polymer-convergence-request"],outputs: [.init(name:"assessment",kind:"vivo.ring-polymer-convergence-result")],
                summary: "Independent-chain diagnostics and consecutive finite-bead convergence for one declared equilibrium observable.",configure: VivoPlatformOperations.empty,
                calculate: { _,input,budget in
                    let r=try VivoPlatformOperations.input(VivoRingPolymerConvergenceRequest.self,"request",input)
                    guard r.configuration.maximumPrimitiveWork <= budget.maximumOperatorApplications else {
                        throw VivoChemistryError.resourceLimit("ring convergence workflow budget")
                    }
                    _=try budget.elements([r.chains.count,max(1,r.chains.map { $0.run.observations.count }.max() ?? 1)],simultaneousArrays:6)
                    return ["assessment":try VivoCanonicalJSON.encode(VivoRingPolymerConvergence.assess(r))]
                }),
            VivoPlatformOperations.pure(identifier: "vivo.platform.reactive-surrogate-coverage",id: id,
                inputs: ["request":"vivo.reactive-surrogate-coverage-request"],outputs: [.init(name:"assessment",kind:"vivo.reactive-surrogate-coverage-result")],
                summary: "Grouped external-domain coverage of a frozen surrogate against authoritative complete-Hamiltonian labels.",configure: VivoPlatformOperations.empty,
                calculate: { _,input,budget in
                    let r=try VivoPlatformOperations.input(VivoReactiveSurrogateCoverageRequest.self,"request",input)
                    guard r.configuration.maximumPrimitiveWork <= budget.maximumOperatorApplications else {
                        throw VivoChemistryError.resourceLimit("reactive coverage workflow budget")
                    }
                    _=try budget.elements([r.labels.count,max(1,r.model.payload.authorityDefinition.atomIndices.count)],simultaneousArrays:8)
                    return ["assessment":try VivoCanonicalJSON.encode(VivoReactiveSurrogateCoverage.assess(r))]
                })]
        return pure+[ring(id),reactive(id),labels(id)]
    }
    private static func ring(_ id: VivoFingerprint) -> VivoWorkflowDefinition {
        let kinds = ["request":"vivo.ring-polymer-workflow-request"]
        let outputs = [VivoChemistryTaskOutput(name: "run",kind: "vivo.ring-polymer-equilibrium-run"),.init(name: "checkpoint",kind: "vivo.ring-polymer-equilibrium-checkpoint")]
        let op = VivoChemistryOperation(identifier: "vivo.platform.ring-polymer-sample",version: "1",implementationFingerprint: id,
            numericalBackend: "metal-fp32-with-fp64-nuclear-hmc",outputs: outputs,
            execute: { cfg,input,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys) else { throw VivoChemistryError.invalid("ring workflow ports") }
                let r = try VivoPlatformOperations.input(VivoRingPolymerWorkflowRequest.self,"request",input); try r.admit(budget)
                let potential = try await r.potential.make(budget: budget)
                let run = try await VivoRingPolymerSampling.run(potential: potential,checkpoint: r.checkpoint,sweeps: r.sweeps)
                try run.validate()
                return ["run":try VivoCanonicalJSON.encode(run),"checkpoint":try VivoCanonicalJSON.encode(run.end)]
            },validateOutputs: { cfg,input,output,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys), Set(output.keys) == Set(outputs.map(\.name)) else { throw VivoChemistryError.invalid("ring validated ports") }
                let r = try VivoPlatformOperations.input(VivoRingPolymerWorkflowRequest.self,"request",input); try r.admit(budget)
                let run = try VivoPlatformOperations.input(VivoRingPolymerRun.self,"run",output)
                let cp = try VivoPlatformOperations.input(VivoRingPolymerCheckpoint.self,"checkpoint",output)
                try run.validate()
                guard run.start == r.checkpoint, run.end == cp, run.end.sweep-run.start.sweep == UInt64(r.sweeps) else { throw VivoChemistryError.invalid("ring workflow continuation binding") }
            })
        return .init(operation: op,inputKinds: kinds,summary: "Resumable fixed-cell quantum-nuclear equilibrium sampling through the existing complete MD/QM-MM force authority.",validationScope: scope,validateConfiguration: VivoPlatformOperations.empty)
    }
    private static func reactive(_ id: VivoFingerprint) -> VivoWorkflowDefinition {
        let kinds = ["request":"vivo.reactive-sampling-workflow-request"]
        let outputs = [VivoChemistryTaskOutput(name: "run",kind: "vivo.reactive-equilibrium-run"),.init(name: "checkpoint",kind: "vivo.reactive-equilibrium-checkpoint")]
        let op = VivoChemistryOperation(identifier: "vivo.platform.reactive-surrogate-sample",version: "1",implementationFingerprint: id,
            numericalBackend: "metal-fp32-with-fp64-authority-corrected-hmc",outputs: outputs,
            execute: { cfg,input,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys) else { throw VivoChemistryError.invalid("reactive workflow ports") }
                let r = try VivoPlatformOperations.input(VivoReactiveSamplingWorkflowRequest.self,"request",input); try r.admit(budget)
                let authority = try await r.authority.make(budget: budget), baseline = try await r.baseline.make(budget: budget)
                let backend: VivoReactiveEnergyForceBackend? = r.metalInference ? try .metal(model: r.model,maximumBatchSize: 1,maximumBytes: min(budget.maximumBytes,1_048_576)) : nil
                let run = try await VivoReactiveSurrogateSampling.run(model: r.model,authority: authority,baseline: baseline,
                    checkpoint: r.checkpoint,sweeps: r.sweeps,backend: backend)
                try run.validate(model: r.model)
                return ["run":try VivoCanonicalJSON.encode(run),"checkpoint":try VivoCanonicalJSON.encode(run.end)]
            },validateOutputs: { cfg,input,output,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys), Set(output.keys) == Set(outputs.map(\.name)) else { throw VivoChemistryError.invalid("reactive validated ports") }
                let r = try VivoPlatformOperations.input(VivoReactiveSamplingWorkflowRequest.self,"request",input); try r.admit(budget)
                let run = try VivoPlatformOperations.input(VivoReactiveSamplingRun.self,"run",output)
                let cp = try VivoPlatformOperations.input(VivoReactiveSamplingCheckpoint.self,"checkpoint",output)
                try run.validate(model: r.model)
                guard run.start == r.checkpoint, run.end == cp, run.end.sweep-run.start.sweep == UInt64(r.sweeps) else { throw VivoChemistryError.invalid("reactive workflow continuation binding") }
            })
        return .init(operation: op,inputKinds: kinds,summary: "Frozen reactive-surrogate proposals with complete authoritative endpoint Metropolis correction and new acquisition labels.",validationScope: scope,validateConfiguration: VivoPlatformOperations.empty)
    }
    private static func labels(_ id: VivoFingerprint) -> VivoWorkflowDefinition {
        let kinds = ["request":"vivo.reactive-label-workflow-request"]
        let outputs = [VivoChemistryTaskOutput(name: "labels",kind: "vivo.reactive-training-labels")]
        let op = VivoChemistryOperation(identifier: "vivo.platform.reactive-surrogate-label",version: "1",implementationFingerprint: id,
            numericalBackend: "metal-fp32-with-fp64-BO",outputs: outputs,
            execute: { cfg,input,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys) else { throw VivoChemistryError.invalid("reactive label ports") }
                let r = try VivoPlatformOperations.input(VivoReactiveLabelWorkflowRequest.self,"request",input); try r.admit(budget)
                let a = try await r.authority.make(budget: budget), b = try await r.baseline.make(budget: budget)
                var labels: [VivoReactiveTrainingLabel] = []
                for g in r.geometries {
                    let exact = try await a.checked(g.positionsNM), base = try await b.checked(g.positionsNM)
                    labels.append(.init(identifier: g.identifier,sourceGroup: g.sourceGroup,authority: exact,baseline: base))
                }
                return ["labels":try VivoCanonicalJSON.encode(labels)]
            },validateOutputs: { cfg,input,output,budget in
                try VivoPlatformOperations.empty(cfg)
                guard Set(input.keys) == Set(kinds.keys), Set(output.keys) == ["labels"] else { throw VivoChemistryError.invalid("validated label ports") }
                let r = try VivoPlatformOperations.input(VivoReactiveLabelWorkflowRequest.self,"request",input); try r.admit(budget)
                let labels = try VivoPlatformOperations.input([VivoReactiveTrainingLabel].self,"labels",output)
                guard labels.count == r.geometries.count else { throw VivoChemistryError.invalid("label count") }
                for (g,label) in zip(r.geometries,labels) {
                    guard label.identifier == g.identifier, label.sourceGroup == g.sourceGroup else { throw VivoChemistryError.invalid("label group identity") }
                    try label.authority.validate(definition: r.authority.definition(),positionsNM: g.positionsNM)
                    try label.baseline.validate(definition: r.baseline.definition(),positionsNM: g.positionsNM)
                }
            })
        return .init(operation: op,inputKinds: kinds,summary: "Full-Hamiltonian energy/force labels at identical mapped geometries; group provenance is preserved for held-out fitting.",validationScope: scope,validateConfiguration: VivoPlatformOperations.empty)
    }
}
