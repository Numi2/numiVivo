import Foundation

public enum VivoElectronicWorkflowCalculation: Codable, Sendable, Equatable {
    case hartreeFock(configuration: VivoSCFConfiguration)
    case lda(configuration: VivoLDAConfiguration)
    case cpcmRHF(configuration: VivoCPCMRHFConfiguration)
    case correlated(reference: VivoSCFConfiguration, solver: VivoManyBodySolverRequest)
}
public struct VivoElectronicWorkflowRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/electronic-workflow/v1"
    public let schema: String
    public let system: VivoElectronicSystem
    public let basis: VivoGaussianBasis
    public let calculation: VivoElectronicWorkflowCalculation
    public let budget: VivoChemistryBudget
    public init(system:VivoElectronicSystem,basis:VivoGaussianBasis,calculation:VivoElectronicWorkflowCalculation,
                budget:VivoChemistryBudget = .init()) {
        schema=Self.schema; self.system=system; self.basis=basis; self.calculation=calculation; self.budget=budget
    }
    public func validate() throws {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("electronic workflow schema") }
        try budget.validate(); try system.validate(); try basis.validate(nucleusCount:system.nuclei.count)
        switch calculation {
        case .hartreeFock(let c): try c.validate()
        case .lda(let c): try c.validate(atomCount:system.nuclei.count); try c.cpcm?.validate(system:system)
        case .cpcmRHF(let c): try c.validate(system:system)
        case .correlated(let reference,let solver):
            try reference.validate()
            guard reference.reference == .restricted, system.alphaElectrons==system.betaElectrons else {
                throw VivoChemistryError.unsupported("molecular correlated preparation currently requires a closed-shell RHF reference; open-spin embedded Hamiltonians can use VivoManyBodySolver directly")
            }
            switch solver {
            case .ccsd(let c): try c.validate()
            case .casscf(_,let c): try c.validate()
            case .mp2(let gap): guard gap.isFinite, gap>0 else { throw VivoChemistryError.invalid("MP2 denominator threshold") }
            default: break
            }
        }
    }
}
public struct VivoElectronicWorkflowPlan: Sendable {
    public let nodes: [VivoChemistryDAGNode]
    public let resultNode: String
    public let resultOutput: String
    public let resultKind: String
}
public enum VivoElectronicWorkflowPlanner {
    /// The actual executable/backend fingerprint is supplied by the application;
    /// arbitrary version labels are not treated as an implementation identity.
    public static func plan(_ request:VivoElectronicWorkflowRequest, systemArtifact:VivoFingerprint,
                            basisArtifact:VivoFingerprint, implementationFingerprint:VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        try request.validate()
        let resources=VivoChemistryResourceContract(budget:request.budget,
            maximumInputBytes:request.budget.maximumBytes,maximumOutputBytes:request.budget.maximumBytes)
        func json<T:Encodable>(_ value:T) throws -> VivoJSONValue {
            try VivoCanonicalJSON.decode(VivoJSONValue.self,from:VivoCanonicalJSON.encode(value))
        }
        let ops=VivoElectronicWorkflowOperations.self, id=implementationFingerprint
        var nodes=[VivoChemistryDAGNode(identifier:"integrals",operation:ops.integrals(implementationFingerprint:id),
            inputs:[.artifact(name:"system",fingerprint:systemArtifact,kind:"vivo.electronic-system"),
                    .artifact(name:"basis",fingerprint:basisArtifact,kind:"vivo.gaussian-basis")],
            configuration:.object([:]),resources:resources)]
        let ao=VivoChemistryDAGInput.output(name:"integrals",node:"integrals",output:"integrals",kind:"vivo.ao-integrals")
        switch request.calculation {
        case .hartreeFock(let cfg):
            nodes.append(.init(identifier:"result",operation:ops.hartreeFock(implementationFingerprint:id),
                inputs:[ao],configuration:try json(cfg),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"scf",resultKind:"vivo.hartree-fock")
        case .lda(let cfg):
            nodes.append(.init(identifier:"result",operation:ops.lda(implementationFingerprint:id),
                inputs:[ao],configuration:try json(cfg),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"lda",resultKind:"vivo.lda")
        case .cpcmRHF(let cfg):
            nodes.append(.init(identifier:"result",operation:ops.cpcmRHF(implementationFingerprint:id),
                inputs:[ao],configuration:try json(cfg),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"scf",resultKind:"vivo.cpcm-rhf")
        case .correlated(let reference,let solver):
            nodes.append(.init(identifier:"reference",operation:ops.hartreeFock(implementationFingerprint:id),
                inputs:[ao],configuration:try json(reference),resources:resources))
            nodes.append(.init(identifier:"hamiltonian",operation:ops.hamiltonian(implementationFingerprint:id),
                inputs:[ao,.output(name:"scf",node:"reference",output:"scf",kind:"vivo.hartree-fock")],
                configuration:try json(reference),resources:resources))
            nodes.append(.init(identifier:"result",operation:ops.manyBody(implementationFingerprint:id),
                inputs:[.output(name:"hamiltonian",node:"hamiltonian",output:"hamiltonian",kind:"vivo.embedded-hamiltonian")],
                configuration:try json(solver),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"manyBody",resultKind:"vivo.many-body-result")
        }
    }
}
