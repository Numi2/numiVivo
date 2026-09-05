import Foundation

public enum VivoNativeChemistryOperations {
    /// Fingerprint must identify the actual build/kernel/compiler/backend set.
    /// It is not a human version string and is never filled with a fake constant.
    public static func hartreeFock(implementationFingerprint: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.gaussian-hartree-fock",version:"1",implementationFingerprint:implementationFingerprint,
              outputs:[.init(name:"integrals",kind:"vivo.ao-integrals"),.init(name:"scf",kind:"vivo.hartree-fock")],
              execute:{ configuration,inputs,budget in
                guard Set(inputs.keys)==Set(["system","basis"]),let systemData=inputs["system"],let basisData=inputs["basis"] else {
                    throw VivoChemistryError.invalid("HF task requires exactly system and basis inputs")
                }
                let system=try VivoCanonicalJSON.decode(VivoElectronicSystem.self,from:systemData)
                let basis=try VivoCanonicalJSON.decode(VivoGaussianBasis.self,from:basisData)
                let cfg=try VivoCanonicalJSON.decode(VivoSCFConfiguration.self,from:VivoCanonicalJSON.encode(configuration))
                let ao=try VivoGaussianIntegralEngine.compute(system:system,basis:basis,budget:budget)
                let scf=try VivoHartreeFock.solve(system:system,integrals:ao,configuration:cfg,budget:budget)
                return ["integrals":try VivoCanonicalJSON.encode(ao),"scf":try VivoCanonicalJSON.encode(scf)]
              },validateOutputs:{ configuration,inputs,outputs,budget in
                guard let aoData=outputs["integrals"],let scfData=outputs["scf"] else { throw VivoChemistryError.invalid("HF output slots") }
                let ao=try VivoCanonicalJSON.decode(VivoAOIntegrals.self,from:aoData)
                let scf=try VivoCanonicalJSON.decode(VivoHartreeFockResult.self,from:scfData)
                guard let systemData=inputs["system"],let basisData=inputs["basis"] else { throw VivoChemistryError.invalid("HF input binding") }
                let system=try VivoCanonicalJSON.decode(VivoElectronicSystem.self,from:systemData)
                let basis=try VivoCanonicalJSON.decode(VivoGaussianBasis.self,from:basisData)
                let cfg=try VivoCanonicalJSON.decode(VivoSCFConfiguration.self,from:VivoCanonicalJSON.encode(configuration))
                guard ao.sourceBasis==basis else { throw VivoChemistryError.invalid("HF basis binding") }
                try VivoHartreeFock.validate(result:scf,system:system,integrals:ao,configuration:cfg,budget:budget)
              })
    }
    public static func configurationInteraction(implementationFingerprint: VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.embedded-ci",version:"1",implementationFingerprint:implementationFingerprint,
              outputs:[.init(name:"ci",kind:"vivo.ci-state")],execute:{ configuration,inputs,budget in
                guard Set(inputs.keys)==Set(["hamiltonian"]),let data=inputs["hamiltonian"] else { throw VivoChemistryError.invalid("CI task requires one Hamiltonian input") }
                let h=try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from:data)
                let method=try VivoCanonicalJSON.decode(VivoCIMethod.self,from:VivoCanonicalJSON.encode(configuration))
                return ["ci":try VivoCanonicalJSON.encode(VivoConfigurationInteraction.solve(h,method:method,budget:budget))]
              },validateOutputs:{ configuration,inputs,outputs,budget in
                guard let data=outputs["ci"] else { throw VivoChemistryError.invalid("CI output slot") }
                let ci=try VivoCanonicalJSON.decode(VivoCIResult.self,from:data);try ci.state.validate(budget:budget)
                guard let source=inputs["hamiltonian"] else { throw VivoChemistryError.invalid("CI source binding") }
                let h=try VivoCanonicalJSON.decode(VivoEmbeddedHamiltonian.self,from:source)
                let method=try VivoCanonicalJSON.decode(VivoCIMethod.self,from:VivoCanonicalJSON.encode(configuration))
                guard ci.method==method,ci.state.orbitalCount==h.orbitalCount,ci.state.alphaElectrons==h.alphaElectrons,
                      ci.state.betaElectrons==h.betaElectrons else { throw VivoChemistryError.invalid("CI sector/method binding") }
                let energy=try VivoCIDensityMatrices.compute(ci.state,budget:budget).energy(of:h)
                guard abs(energy-ci.energyHartree)<1e-8 else { throw VivoChemistryError.invalid("CI RDM energy binding") }
                guard ci.energyHartree.isFinite,ci.eigenResidualNorm.isFinite,ci.eigenResidualNorm>=0,ci.eigenResidualNorm<1e-8 else {
                    throw VivoChemistryError.invalid("CI output eigenpair")
                }
              })
    }
}
