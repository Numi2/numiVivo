import Foundation

public enum VivoReactionCalculation:Codable,Sendable,Equatable {
    case qualify(request:VivoNuclearQualificationRequest)
    case connectedReaction(request:VivoConnectedReactionRequest)
    case residualBarrier(request:VivoResidualBarrierRequest)
    case globalEmbedding(request:VivoVariationalEmbeddingRequest)
    case eccSolventClosure(request:VivoECCSolventClosureRequest)
    case connectivity(request:VivoReactionConnectivityRequest)
    case transitionStateTheory(request:VivoTransitionStateTheoryRequest)
    case reproductionPreflight(package:VivoReproductionPackage)
    case barrierConvergence(request:VivoBarrierConvergenceRequest)
    case correlatedSolvent(request:VivoCorrelatedSolventRequest)
    case solvatedPath(request:VivoSolvatedECCPathRequest)
    case harmonicBarrier(saddle:VivoNuclearQualifiedPoint,reactants:[VivoNuclearQualifiedPoint])
    case descent(saddle:VivoNuclearQualifiedPoint,configuration:VivoNuclearDescentConfiguration)
}
public struct VivoReactionCalculationRequest:Codable,Sendable,Equatable {
    public static let schema="numivivo.org/reaction-calculation/v1"
    public let schema:String
    public let calculation:VivoReactionCalculation
    public init(_ calculation:VivoReactionCalculation) {schema=Self.schema;self.calculation=calculation}
    public var budget:VivoChemistryBudget {
        switch calculation {
        case .qualify(let r):return r.model.budget
        case .connectedReaction(let r):return r.saddle.model.budget
        case .residualBarrier(let r):return r.baseline.budget
        case .globalEmbedding(let r):return r.molecule.budget
        case .eccSolventClosure(let r):return r.molecule.budget
        case .connectivity(let r):return r.saddle.request.model.budget
        case .transitionStateTheory(let r):return r.connectivity.request.saddle.request.model.budget
        case .reproductionPreflight(let r):return r.budget
        case .barrierConvergence(let r):return r.budget
        case .correlatedSolvent(let r):return r.budget
        case .solvatedPath(let r):return r.path.budget
        case .harmonicBarrier(let saddle,_),.descent(let saddle,_):return saddle.request.model.budget
        }
    }
    public func validate() throws {
        guard schema==Self.schema else {throw VivoChemistryError.invalid("unknown reaction-calculation schema")}
        try budget.validate()
        switch calculation {
        case .qualify(let r):try r.validate()
        case .connectedReaction(let r):try r.validate()
        case .residualBarrier(let r):try r.validate()
        case .globalEmbedding(let r):try r.validate()
        case .eccSolventClosure(let r):try r.validate()
        case .connectivity(let r):
            try VivoReactionConnectivity.validateRequest(r)
            guard r.endpoints.allSatisfy({ $0.components.allSatisfy { $0.point.request.model.budget==budget } }) else {
                throw VivoChemistryError.invalid("connectivity endpoint resource contracts differ from the task")
            }
        case .transitionStateTheory(let r):
            try r.validate()
            guard r.connectivity.request.endpoints.allSatisfy({ $0.components.allSatisfy { $0.point.request.model.budget==budget } }) else {
                throw VivoChemistryError.invalid("TST endpoint resource contracts differ from the saddle task")
            }
        case .reproductionPreflight(let r):_ = try VivoReproductionPreflight.inspect(r)
        case .barrierConvergence(let r):try r.validate()
        case .correlatedSolvent(let r):try r.validate()
        case .solvatedPath(let r):try r.validate()
        case .harmonicBarrier(let saddle,let reactants):
            try saddle.request.validate()
            guard (1...16).contains(reactants.count),reactants.allSatisfy({$0.request.model.budget==budget}) else {
                throw VivoChemistryError.invalid("barrier input population or resource contracts")
            }
        case .descent(let saddle,_):try saddle.request.validate()
        }
    }
}
public enum VivoReactionCalculationResult:Codable,Sendable,Equatable {
    case qualified(point:VivoNuclearQualifiedPoint)
    case connectedReaction(result:VivoTransitionStateTheoryResult)
    case residualBarrier(result:VivoResidualBarrierResult)
    case globalEmbedding(result:VivoVariationalEmbeddingResult)
    case eccSolventClosure(result:VivoECCSolventClosureResult)
    case connectivity(result:VivoReactionConnectivityResult)
    case transitionStateTheory(result:VivoTransitionStateTheoryResult)
    case reproductionPreflight(result:VivoReproductionReadiness)
    case barrierConvergence(result:VivoBarrierConvergenceResult)
    case correlatedSolvent(result:VivoCorrelatedSolventResult)
    case solvatedPath(result:VivoSolvatedECCPathResult)
    case harmonicBarrier(result:VivoHarmonicBarrierEstimate)
    case descent(result:VivoNuclearDescentResult)
}
public enum VivoReactionQualificationWorkflow {
    public static func execute(_ request:VivoReactionCalculationRequest) throws -> VivoReactionCalculationResult {
        try request.validate()
        switch request.calculation {
        case .qualify(let r):return .qualified(point:try VivoNuclearQualification.run(r))
        case .connectedReaction(let r):return .connectedReaction(result:try VivoConnectedReaction.run(r))
        case .residualBarrier(let r):return .residualBarrier(result:try VivoResidualBarrierCampaign.run(r))
        case .globalEmbedding(let r):return .globalEmbedding(result:try VivoVariationalEmbedding.run(r))
        case .eccSolventClosure(let r):return .eccSolventClosure(result:try VivoECCSolventClosure.run(r))
        case .connectivity(let r):
            let result=try VivoReactionConnectivity.run(r)
            guard result.converged else { throw VivoChemistryError.convergence("mapped connection did not pass independent refinement") }
            return .connectivity(result:result)
        case .transitionStateTheory(let r):return .transitionStateTheory(result:try VivoTransitionStateTheory.estimate(r))
        case .reproductionPreflight(let r):return .reproductionPreflight(result:try VivoReproductionPreflight.inspect(r))
        case .barrierConvergence(let r):return .barrierConvergence(result:try VivoBarrierConvergence.run(r))
        case .correlatedSolvent(let r):return .correlatedSolvent(result:try VivoCorrelatedSolvation.solve(r))
        case .solvatedPath(let r):return .solvatedPath(result:try VivoSolvatedECCPath.solve(r))
        case .harmonicBarrier(let saddle,let reactants):return .harmonicBarrier(result:try VivoHarmonicBarrier.estimate(saddle:saddle,reactants:reactants))
        case .descent(let saddle,let cfg):return .descent(result:try VivoNuclearDescent.trace(saddle,configuration:cfg))
        }
    }
    public static func validate(_ result:VivoReactionCalculationResult,request:VivoReactionCalculationRequest) throws {
        try request.validate()
        switch (request.calculation,result) {
        case (.qualify(let r),.qualified(let point)):try VivoNuclearQualification.validate(point,request:r)
        case (.connectedReaction(let r),.connectedReaction(let result)):try VivoConnectedReaction.validate(result,request:r)
        case (.residualBarrier(let r),.residualBarrier(let result)):try VivoResidualBarrierCampaign.validate(result,request:r)
        case (.globalEmbedding(let r),.globalEmbedding(let result)):try VivoVariationalEmbedding.validate(result,request:r)
        case (.eccSolventClosure(let r),.eccSolventClosure(let result)):try VivoECCSolventClosure.validate(result,request:r)
        case (.connectivity(let r),.connectivity(let result)):try VivoReactionConnectivity.validate(result,request:r)
        case (.transitionStateTheory(let r),.transitionStateTheory(let result)):try VivoTransitionStateTheory.validate(result,request:r)
        case (.reproductionPreflight(let r),.reproductionPreflight(let result)):
            guard result == (try VivoReproductionPreflight.inspect(r)) else { throw VivoChemistryError.invalid("reproduction preflight result binding") }
        case (.barrierConvergence(let r),.barrierConvergence(let report)):try VivoBarrierConvergence.validate(report,request:r)
        case (.correlatedSolvent(let r),.correlatedSolvent(let result)):try VivoCorrelatedSolvation.validate(result,request:r)
        case (.solvatedPath(let r),.solvatedPath(let path)):try VivoSolvatedECCPath.validate(path,request:r)
        case (.harmonicBarrier(let saddle,let reactants),.harmonicBarrier(let barrier)):
            let rebuilt=try VivoHarmonicBarrier.estimate(saddle:saddle,reactants:reactants)
            guard barrier==rebuilt else {throw VivoChemistryError.invalid("harmonic barrier reconstruction differs")}
        case (.descent(let saddle,let cfg),.descent(let result)):
            let rebuilt=try VivoNuclearDescent.trace(saddle,configuration:cfg)
            guard rebuilt==result else {throw VivoChemistryError.invalid("descent reconstruction differs")}
        default:throw VivoChemistryError.invalid("reaction output has a different calculation than its request")
        }
    }
    public static func operation(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.reaction-qualification",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"result",kind:"vivo.reaction-calculation-result")],execute:{ cfg,inputs,budget in
                guard cfg == .object([:]),Set(inputs.keys)==Set(["request"]),let data=inputs["request"] else {throw VivoChemistryError.invalid("reaction operation slots")}
                let request=try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self,from:data)
                guard request.budget==budget else {throw VivoChemistryError.invalid("reaction budget binding")}
                return ["result":try VivoCanonicalJSON.encode(execute(request))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                guard cfg == .object([:]),Set(inputs.keys)==Set(["request"]),Set(outputs.keys)==Set(["result"]),
                      let data=inputs["request"],let output=outputs["result"] else {throw VivoChemistryError.invalid("reaction validation slots")}
                let request=try VivoCanonicalJSON.decode(VivoReactionCalculationRequest.self,from:data)
                guard request.budget==budget else {throw VivoChemistryError.invalid("reaction validation budget")}
                try validate(VivoCanonicalJSON.decode(VivoReactionCalculationResult.self,from:output),request:request)
            })
    }
    public static func template(_ name:String) throws -> VivoReactionCalculationRequest {
        if name=="h3-connected-rate" {
            func seed(_ name:String) throws -> VivoNuclearQualificationRequest {
                guard case .qualify(let request)=try template(name).calculation else { throw VivoChemistryError.invalid("nuclear seed template") }
                return request
            }
            let saddle=try seed("h3-saddle"),h2=try seed("h2-minimum"),atom=try seed("h-atom")
            return .init(.connectedReaction(request:.init(atomIdentifiers:["H0","H1","H2"],saddle:saddle,
                endpoints:[.init(identifier:"H0-H1_plus_H2",components:[.init(atomIndices:[0,1],qualification:h2),.init(atomIndices:[2],qualification:atom)]),
                           .init(identifier:"H0_plus_H1-H2",components:[.init(atomIndices:[0],qualification:atom),.init(atomIndices:[1,2],qualification:h2)])],
                reactantEndpointIdentifier:"H0-H1_plus_H2",connectivity:.init(initialDisplacementMassWeighted:0.015,
                    stepMassWeighted:0.04,endpointMaximumGradient:1e-7,comparisonSamples:128))))
        }

        if name=="h3-residual-barrier" { return .init(.residualBarrier(request:VivoResidualBarrierCampaign.template())) }
        if name=="h2-global-embedding" { return .init(.globalEmbedding(request:VivoVariationalEmbedding.template())) }
        if name=="h2-ecc-solvent-closure" { return .init(.eccSolventClosure(request:VivoECCSolventClosure.hydrogenControl())) }
        if name=="paper-michael-inputs" { return .init(.reproductionPreflight(package:.init(target:.acrylamideMethanethiolate,resultIdentifier:"supply-exact-paper-figure-table-method"))) }
        if name=="paper-btk-inputs" { return .init(.reproductionPreflight(package:.init(target:.btkSnapshot,resultIdentifier:"supply-exact-snapshot-method-result"))) }
        if name=="h3-barrier-convergence" {return .init(.barrierConvergence(request:VivoBarrierBenchmarks.hydrogenExchange631G()))}
        if name=="h3-barrier-convergence-ensemble" {return .init(.barrierConvergence(request:VivoBarrierBenchmarks.hydrogenExchange631G(ensemble:true)))}
        if name=="h2-solvated-path" {return .init(.solvatedPath(request:.init(path:try VivoMolecularECCPath.hydrogenStretchTemplate(),solvent:.init(dielectricConstant:4,angularPoints:50))))}
        if name=="h2-equilibrium-cpcm" {
            let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
                .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
            return .init(.correlatedSolvent(request:.init(system:system,basis:.hydrogenSTO3G(nucleusIndices:[0,1]),
                solvent:.init(dielectricConstant:4,angularPoints:50))))
        }
        if name=="h2-equilibrium-minimum" {
            guard case .qualify(var request)=try template("h2-minimum").calculation else {
                throw VivoChemistryError.invalid("internal nuclear template binding")
            }
            request.model.solver = .equilibriumFullCI
            request.model.solvent = .init(dielectricConstant:4,angularPoints:50)
            request.model.correlatedSolventConfiguration = .init(densityTolerance:1e-10,
                potentialToleranceHartree:1e-10,energyToleranceHartree:1e-12,ciResidualTolerance:1e-13)
            return .init(.qualify(request:request))
        }
        let n:Int,positions:[SIMD3<Double>],na:Int,nb:Int,operation:VivoNuclearOperation,kind:VivoStationaryKind,symmetry:Int,degeneracy:Int
        switch name {
        case "h2-minimum":n=2;positions=[.init(0,0,-0.75),.init(0,0,0.75)];na=1;nb=1;operation = .minimize;kind = .minimum;symmetry=2;degeneracy=1
        case "h3-saddle":n=3;positions=[.init(0,0,-1.85),.zero,.init(0,0,1.85)];na=2;nb=1;operation = .refineSaddle;kind = .firstOrderSaddle;symmetry=2;degeneracy=2
        case "h-atom":n=1;positions=[.zero];na=1;nb=0;operation = .characterize;kind = .minimum;symmetry=1;degeneracy=2
        default:throw VivoChemistryError.invalid("unknown reaction template")
        }
        let system=VivoElectronicSystem(nuclei:positions.enumerated().map{.init(atomicNumber:1,positionBohr:$0.element,structureAtomIndex:UInt32($0.offset))},alphaElectrons:na,betaElectrons:nb)
        let request=VivoNuclearQualificationRequest(model:.init(system:system,basis:.hydrogenSTO3G(nucleusIndices:Array(0..<n)),solver:.fullCI),
            massesDa:[Double](repeating:1.008,count:n),operation:operation,kind:kind,
            thermochemistry:.init(rotationalSymmetryNumber:symmetry,electronicDegeneracy:degeneracy))
        return .init(.qualify(request:request))
    }
}
