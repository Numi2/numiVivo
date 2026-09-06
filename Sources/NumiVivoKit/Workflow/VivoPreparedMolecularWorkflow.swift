import Foundation

public enum VivoPreparedMolecularCalculation: Codable, Sendable, Equatable {
    case prepare(request:VivoMolecularPreparationRequest)
    case sampling(request:VivoMolecularSamplingRequest)
    case rate(request:VivoPreparedReactionRateRequest,kinetics:VivoCovalentKineticPack?)
    public var budget:VivoChemistryBudget {
        switch self {
        case .rate(let r,_):return r.transitionState.connectivity.request.saddle.request.model.budget
        default:return .init()
        }
    }
}
public enum VivoPreparedMolecularCalculationResult: Codable, Sendable, Equatable {
    case prepared(result:VivoMolecularPreparationResult)
    case sampling(result:VivoMolecularSamplingResult)
    case rate(result:VivoPreparedReactionRateResult,kinetics:VivoCovalentKineticPack?)
}

/// Pure preparation, sampling analysis and rate derivation use the same verified
/// task/DAG executor as electronic chemistry. Adaptive MD itself retains its
/// existing exact checkpoints and append-only coordinate archive.
public enum VivoPreparedMolecularWorkflow {
    public static func execute(_ calculation:VivoPreparedMolecularCalculation) throws -> VivoPreparedMolecularCalculationResult {
        switch calculation {
        case .prepare(let r):return .prepared(result:try VivoMolecularPreparation.prepare(r))
        case .sampling(let r):return .sampling(result:try VivoMolecularSampling.analyze(r))
        case .rate(let r,let model):
            let result=try VivoPreparedReactionRate.calculate(r)
            let updated:VivoCovalentKineticPack?
            if let model {
                updated=try VivoTransitionStateDerivation.calculate(result.kineticRequest).applyingInactivation(to:model)
            } else { updated=nil }
            return .rate(result:result,kinetics:updated)
        }
    }
    public static func operation(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.prepared-molecular-calculation",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"result",kind:"vivo.prepared-molecular-result")],execute:{ cfg,inputs,budget in
                guard cfg == .object([:]),Set(inputs.keys)==Set(["request"]),let data=inputs["request"] else {
                    throw VivoChemistryError.invalid("prepared molecular task input contract")
                }
                let request=try VivoCanonicalJSON.decode(VivoPreparedMolecularCalculation.self,from:data)
                guard request.budget==budget else { throw VivoChemistryError.invalid("prepared molecular budget binding") }
                return ["result":try VivoCanonicalJSON.encode(execute(request))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                guard cfg == .object([:]),Set(inputs.keys)==Set(["request"]),Set(outputs.keys)==Set(["result"]),
                      let data=inputs["request"],let result=outputs["result"] else {
                    throw VivoChemistryError.invalid("prepared molecular output contract")
                }
                let request=try VivoCanonicalJSON.decode(VivoPreparedMolecularCalculation.self,from:data)
                guard request.budget==budget else { throw VivoChemistryError.invalid("prepared molecular validation budget") }
                let actual=try VivoCanonicalJSON.decode(VivoPreparedMolecularCalculationResult.self,from:result)
                guard actual == (try execute(request)) else {
                    throw VivoChemistryError.invalid("prepared molecular output does not reconstruct from its source artifacts")
                }
            })
    }
}
