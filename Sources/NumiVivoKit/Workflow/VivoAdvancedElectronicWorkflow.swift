import Foundation

public enum VivoAdvancedElectronicCalculation:Codable,Sendable,Equatable {
    case densityFittedMP2(configuration:VivoRIWorkflowConfiguration)
    case densityFittedCorrelated(configuration:VivoRIWorkflowConfiguration,space:VivoRIActiveSpace,solver:VivoAdvancedManyBodyRequest)
    case densityFittedEmbedding(configuration:VivoRIWorkflowConfiguration,space:VivoRIActiveSpace,embedding:VivoECCWorkflowConfiguration)
    case densityFittedDMET(configuration:VivoRIWorkflowConfiguration,space:VivoRIActiveSpace,embedding:VivoECCDMETConfiguration)
    case smoothCPCM(configuration:VivoSmoothSCFWorkflowConfiguration)
    case correlated(reference:VivoSCFConfiguration,solver:VivoAdvancedManyBodyRequest)
    case correlatedEmbedding(configuration:VivoECCWorkflowConfiguration)
    case eccDMET(configuration:VivoECCDMETConfiguration)
}
/// Both electronic request schemas use the same artifact store, DAG scheduler
/// and CLI commands. Added RI cases preserve decoding of existing requests.
public struct VivoAdvancedElectronicWorkflowRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/advanced-electronic-workflow/v1"
    public let schema:String
    public let system:VivoElectronicSystem
    public let basis:VivoGaussianBasis
    public let calculation:VivoAdvancedElectronicCalculation
    public let budget:VivoChemistryBudget
    public init(system:VivoElectronicSystem,basis:VivoGaussianBasis,calculation:VivoAdvancedElectronicCalculation,budget:VivoChemistryBudget = .init()) {
        schema=Self.schema;self.system=system;self.basis=basis;self.calculation=calculation;self.budget=budget
    }
    public func validate() throws {
        guard schema==Self.schema else { throw VivoChemistryError.invalid("advanced electronic workflow schema") }
        try system.validate();try basis.validate(nucleusCount:system.nuclei.count);try budget.validate()
        func validateRI(_ c:VivoRIWorkflowConfiguration) throws {
            try c.auxiliaryBasis.validate(nucleusCount:system.nuclei.count);try c.scf.validate()
            guard c.scf.reference == .restricted,system.alphaElectrons==system.betaElectrons,
                  c.fitting.metricEigenvalueThreshold.isFinite,c.fitting.metricEigenvalueThreshold>0,
                  (1...8192).contains(c.fitting.maximumAuxiliaryFunctions),c.fitting.pairBlockSize>0 else {
                throw VivoChemistryError.unsupported("RI workflows require a closed-shell RHF reference and valid fitting configuration")
            }
        }
        switch calculation {
        case .densityFittedMP2(let c):try validateRI(c)
        case .densityFittedCorrelated(let c,let space,_),
             .densityFittedEmbedding(let c,let space,_),
             .densityFittedDMET(let c,let space,_):
            try validateRI(c)
            let count=try VivoGaussianIntegralEngine.expanded(system:system,basis:basis,budget:budget).count
            try space.validate(orbitalCount:count,occupied:system.alphaElectrons)
            let m=space.activeOrbitals.count
            _ = try budget.elements([m,m,m,m],simultaneousArrays:4)
        case .smoothCPCM(let c):try c.scf.validate()
        case .correlated(let reference,_):
            try reference.validate()
            guard reference.reference == .restricted,system.alphaElectrons==system.betaElectrons else {
                throw VivoChemistryError.unsupported("correlated molecular preparation uses RHF; direct Hamiltonian solvers preserve their own spin contract")
            }
        case .correlatedEmbedding, .eccDMET:
            guard system.alphaElectrons==system.betaElectrons else {
                throw VivoChemistryError.unsupported("correlated moment matching currently requires a non-spin-polarized reference")
            }
        }
    }
}
public enum VivoAdvancedElectronicWorkflowPlanner {
    public static func plan(_ request:VivoAdvancedElectronicWorkflowRequest,systemArtifact:VivoFingerprint,basisArtifact:VivoFingerprint,
                            implementationFingerprint id:VivoFingerprint) throws -> VivoElectronicWorkflowPlan {
        try request.validate()
        let resources=VivoChemistryResourceContract(budget:request.budget,maximumInputBytes:request.budget.maximumBytes,maximumOutputBytes:request.budget.maximumBytes)
        func json<T:Encodable>(_ value:T) throws -> VivoJSONValue {
            try VivoCanonicalJSON.decode(VivoJSONValue.self,from:VivoCanonicalJSON.encode(value))
        }
        let inputs:[VivoChemistryDAGInput]=[.artifact(name:"system",fingerprint:systemArtifact,kind:"vivo.electronic-system"),
            .artifact(name:"basis",fingerprint:basisArtifact,kind:"vivo.gaussian-basis")]
        let ops=VivoAdvancedChemistryOperations.self,base=VivoElectronicWorkflowOperations.self
        var nodes:[VivoChemistryDAGNode]=[]
        let riConfiguration:VivoRIWorkflowConfiguration?,space:VivoRIActiveSpace?
        switch request.calculation {
        case .densityFittedMP2(let c):riConfiguration=c;space=nil
        case .densityFittedCorrelated(let c,let s,_),.densityFittedEmbedding(let c,let s,_),.densityFittedDMET(let c,let s,_):
            riConfiguration=c;space=s
        default:riConfiguration=nil;space=nil
        }
        if let c=riConfiguration {
            nodes.append(.init(identifier:"integrals",operation:ops.ri(implementationFingerprint:id),inputs:inputs,configuration:try json(c),resources:resources))
            nodes.append(.init(identifier:"reference",operation:ops.riHF(implementationFingerprint:id),
                inputs:[.output(name:"integrals",node:"integrals",output:"integrals",kind:"vivo.ri-integrals")],configuration:try json(c.scf),resources:resources))
            let reference=VivoChemistryDAGInput.output(name:"reference",node:"reference",output:"reference",kind:"vivo.ri-hf")
            if let space {
                nodes.append(.init(identifier:"hamiltonian",operation:ops.riHamiltonian(implementationFingerprint:id),
                    inputs:[reference],configuration:try json(space),resources:resources))
            } else {
                nodes.append(.init(identifier:"result",operation:ops.riMP2(implementationFingerprint:id),
                    inputs:[reference],configuration:.object([:]),resources:resources))
                return .init(nodes:nodes,resultNode:"result",resultOutput:"mp2",resultKind:"vivo.ri-mp2")
            }
        } else {
            nodes.append(.init(identifier:"integrals",operation:base.integrals(implementationFingerprint:id),inputs:inputs,configuration:.object([:]),resources:resources))
            let ao=VivoChemistryDAGInput.output(name:"integrals",node:"integrals",output:"integrals",kind:"vivo.ao-integrals")
            if case .smoothCPCM(let c)=request.calculation {
                nodes.append(.init(identifier:"result",operation:ops.smoothSCF(implementationFingerprint:id),inputs:[ao],configuration:try json(c),resources:resources))
                return .init(nodes:nodes,resultNode:"result",resultOutput:"scf",resultKind:"vivo.smooth-cpcm-scf")
            }
            let reference:VivoSCFConfiguration
            if case .correlated(let r,_)=request.calculation { reference=r } else { reference = .init() }
            nodes.append(.init(identifier:"reference",operation:base.hartreeFock(implementationFingerprint:id),inputs:[ao],configuration:try json(reference),resources:resources))
            nodes.append(.init(identifier:"hamiltonian",operation:base.hamiltonian(implementationFingerprint:id),
                inputs:[ao,.output(name:"scf",node:"reference",output:"scf",kind:"vivo.hartree-fock")],configuration:try json(reference),resources:resources))
        }
        // The downstream operations are deliberately identical for dense and RI
        // preparation. No duplicate correlated, ECC or DMET solver is introduced.
        let h=VivoChemistryDAGInput.output(name:"hamiltonian",node:"hamiltonian",output:"hamiltonian",kind:"vivo.embedded-hamiltonian")
        switch request.calculation {
        case .correlated(_,let solver),.densityFittedCorrelated(_,_,let solver):
            nodes.append(.init(identifier:"result",operation:ops.manyBody(implementationFingerprint:id),inputs:[h],configuration:try json(solver),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"manyBody",resultKind:"vivo.advanced-many-body-result")
        case .correlatedEmbedding(let c),.densityFittedEmbedding(_,_,let c):
            nodes.append(.init(identifier:"result",operation:ops.ecc(implementationFingerprint:id),inputs:[h],configuration:try json(c),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"embedding",resultKind:"vivo.ecc-selected-moment-result")
        case .eccDMET(let c),.densityFittedDMET(_,_,let c):
            nodes.append(.init(identifier:"result",operation:ops.eccDMET(implementationFingerprint:id),inputs:[h],configuration:try json(c),resources:resources))
            return .init(nodes:nodes,resultNode:"result",resultOutput:"embedding",resultKind:"vivo.ecc-dmet-result")
        default:throw VivoChemistryError.invalid("unresolved advanced calculation")
        }
    }
}
