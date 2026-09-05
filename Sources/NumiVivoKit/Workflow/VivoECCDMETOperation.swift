import Foundation

public extension VivoAdvancedChemistryOperations {
    static func eccDMET(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.ecc-dmet",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"embedding",kind:"vivo.ecc-dmet-result")],execute:{cfg,inputs,budget in
                guard Set(inputs.keys)==["hamiltonian"],let data=inputs["hamiltonian"] else{throw VivoChemistryError.invalid("ECC-DMET input slots")}
                let h=try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from:data)
                let settings=try VivoCanonicalJSON.decode(VivoECCDMETConfiguration.self,from:VivoCanonicalJSON.encode(cfg))
                let result=try VivoECCDMET.solve(h,configuration:settings,budget:budget)
                guard result.converged else{throw VivoChemistryError.convergence("ECC-DMET: \(result.termination)")}
                return ["embedding":try VivoCanonicalJSON.encode(result)]
            },validateOutputs:{cfg,inputs,outputs,budget in
                guard Set(inputs.keys)==["hamiltonian"],Set(outputs.keys)==["embedding"],
                      let input=inputs["hamiltonian"],let output=outputs["embedding"] else{throw VivoChemistryError.invalid("ECC-DMET output slots")}
                let h=try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from:input)
                let result=try VivoCanonicalJSON.decode(VivoECCDMETResult.self,from:output)
                let settings=try VivoCanonicalJSON.decode(VivoECCDMETConfiguration.self,from:VivoCanonicalJSON.encode(cfg))
                try VivoECCDMET.validate(result,hamiltonian:h,configuration:settings,budget:budget)
            })
    }
}
