import Foundation

public enum VivoReactionCalculation:Codable,Sendable,Equatable {
    case qualify(request:VivoNuclearQualificationRequest)
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
        case .solvatedPath(let r):return r.path.budget
        case .harmonicBarrier(let saddle,_),.descent(let saddle,_):return saddle.request.model.budget
        }
    }
    public func validate() throws {
        guard schema==Self.schema else {throw VivoChemistryError.invalid("unknown reaction-calculation schema")}
        try budget.validate()
        switch calculation {
        case .qualify(let r):try r.validate()
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
    case solvatedPath(result:VivoSolvatedECCPathResult)
    case harmonicBarrier(result:VivoHarmonicBarrierEstimate)
    case descent(result:VivoNuclearDescentResult)
}
public enum VivoReactionQualificationWorkflow {
    public static func execute(_ request:VivoReactionCalculationRequest) throws -> VivoReactionCalculationResult {
        try request.validate()
        switch request.calculation {
        case .qualify(let r):return .qualified(point:try VivoNuclearQualification.run(r))
        case .solvatedPath(let r):return .solvatedPath(result:try VivoSolvatedECCPath.solve(r))
        case .harmonicBarrier(let saddle,let reactants):return .harmonicBarrier(result:try VivoHarmonicBarrier.estimate(saddle:saddle,reactants:reactants))
        case .descent(let saddle,let cfg):return .descent(result:try VivoNuclearDescent.trace(saddle,configuration:cfg))
        }
    }
    public static func validate(_ result:VivoReactionCalculationResult,request:VivoReactionCalculationRequest) throws {
        try request.validate()
        switch (request.calculation,result) {
        case (.qualify(let r),.qualified(let point)):try VivoNuclearQualification.validate(point,request:r)
        case (.solvatedPath(let r),.solvatedPath(let path)):try VivoSolvatedECCPath.validate(path,request:r)
        case (.harmonicBarrier(let saddle,let reactants),.harmonicBarrier(let barrier)):
            let rebuilt=try VivoHarmonicBarrier.estimate(saddle:saddle,reactants:reactants)
            guard barrier==rebuilt else {throw VivoChemistryError.invalid("harmonic barrier reconstruction differs")}
        case (.descent(let saddle,let cfg),.descent(let result)):
            // Deterministic local calculation: rerun the branch integration, not
            // merely a check of the reported endpoint flag.
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
        if name=="h2-solvated-path" {return .init(.solvatedPath(request:.init(path:try VivoMolecularECCPath.hydrogenStretchTemplate(),solvent:.init(dielectricConstant:4,angularPoints:50))))}
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
