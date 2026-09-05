import Foundation

public struct VivoRIWorkflowConfiguration:Codable,Sendable,Equatable {
    public let auxiliaryBasis:VivoGaussianBasis
    public let fitting:VivoRIConfiguration
    public let scf:VivoSCFConfiguration
    public init(auxiliaryBasis:VivoGaussianBasis,fitting:VivoRIConfiguration = .init(),scf:VivoSCFConfiguration = .init()) {
        self.auxiliaryBasis=auxiliaryBasis;self.fitting=fitting;self.scf=scf
    }
}
public struct VivoSmoothSCFWorkflowConfiguration:Codable,Sendable,Equatable {
    public let solvent:VivoSmoothCPCMConfiguration
    public let scf:VivoSCFConfiguration
    public init(solvent:VivoSmoothCPCMConfiguration = .init(),scf:VivoSCFConfiguration = .init()) { self.solvent=solvent;self.scf=scf }
}
public struct VivoECCWorkflowConfiguration:Codable,Sendable,Equatable {
    public let fragments:[VivoECCFragment]
    public let operators:[VivoECCCorrelationOperator]
    public let settings:VivoECCSelfConsistencyConfiguration
    public init(fragments:[VivoECCFragment],operators:[VivoECCCorrelationOperator],settings:VivoECCSelfConsistencyConfiguration = .init()) {
        self.fragments=fragments;self.operators=operators;self.settings=settings
    }
}
/// Versioned operations in the existing verified artifact DAG. No second store,
/// no numerical provider dependency, and no fallback to a different method.
public enum VivoAdvancedChemistryOperations {
    private static func decode<T:Decodable>(_ type:T.Type,_ cfg:VivoJSONValue) throws -> T {
        try VivoCanonicalJSON.decode(type,from:VivoCanonicalJSON.encode(cfg))
    }
    private static func input<T:Decodable>(_ type:T.Type,_ key:String,_ values:[String:Data]) throws -> T {
        guard let data=values[key] else { throw VivoChemistryError.invalid("missing advanced chemistry slot \(key)") }
        return try VivoCanonicalJSON.decode(type,from:data)
    }
    private static func slots(_ data:[String:Data],_ expected:[String]) throws {
        guard Set(data.keys)==Set(expected) else { throw VivoChemistryError.invalid("advanced chemistry slot contract") }
    }
    public static func manyBody(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.advanced-many-body",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"manyBody",kind:"vivo.advanced-many-body-result")],execute:{ cfg,inputs,budget in
                try slots(inputs,["hamiltonian"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs),request=try decode(VivoAdvancedManyBodyRequest.self,cfg)
                return ["manyBody":try VivoCanonicalJSON.encode(VivoManyBodySolver.solve(h,request:request,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["hamiltonian"]);try slots(outputs,["manyBody"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs),result=try input(VivoAdvancedManyBodyResult.self,"manyBody",outputs)
                try VivoManyBodySolver.validate(result,hamiltonian:h,request:decode(VivoAdvancedManyBodyRequest.self,cfg),budget:budget)
            })
    }
    public static func ri(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.ri-integrals",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"integrals",kind:"vivo.ri-integrals")],execute:{ cfg,inputs,budget in
                try slots(inputs,["system","basis"])
                let s=try input(VivoElectronicSystem.self,"system",inputs),b=try input(VivoGaussianBasis.self,"basis",inputs)
                let c=try decode(VivoRIWorkflowConfiguration.self,cfg)
                return ["integrals":try VivoCanonicalJSON.encode(VivoDensityFitting.compute(system:s,basis:b,auxiliaryBasis:c.auxiliaryBasis,configuration:c.fitting,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["system","basis"]);try slots(outputs,["integrals"])
                let r=try input(VivoRIIntegrals.self,"integrals",outputs),c=try decode(VivoRIWorkflowConfiguration.self,cfg)
                try r.oneElectron.validate(budget:budget);try r.factors.validate(budget:budget)
                guard r.oneElectron.system == (try input(VivoElectronicSystem.self,"system",inputs)),
                      r.oneElectron.basis == (try input(VivoGaussianBasis.self,"basis",inputs)),r.auxiliaryBasis==c.auxiliaryBasis,
                      r.configuration==c.fitting,r.factors.orbitalCount==r.oneElectron.count,r.whiteningDefect.isFinite,r.whiteningDefect<1e-5 else {
                    throw VivoChemistryError.invalid("RI output source/metric identity")
                }
            })
    }
    public static func riHF(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.ri-hf",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"reference",kind:"vivo.ri-hf")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals"])
                let ri=try input(VivoRIIntegrals.self,"integrals",inputs),settings=try decode(VivoSCFConfiguration.self,cfg)
                return ["reference":try VivoCanonicalJSON.encode(VivoFactorizedHartreeFock.solve(ri,configuration:settings,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals"]);try slots(outputs,["reference"])
                let r=try input(VivoRIHartreeFockResult.self,"reference",outputs),settings=try decode(VivoSCFConfiguration.self,cfg)
                guard r.source == (try input(VivoRIIntegrals.self,"integrals",inputs)),r.configuration==settings,
                      r.scf.finalCommutatorNorm.isFinite,r.scf.finalCommutatorNorm<=settings.commutatorTolerance else {
                    throw VivoChemistryError.invalid("RI-HF result binding")
                }
                // This operation is the closed-shell RI-MP2 preparation path.
                // Rebuilding its RI Fock in MP2 also verifies reference binding.
                _ = try VivoFactorizedMP2.solve(r,budget:budget)
            })
    }
    public static func riMP2(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.ri-mp2",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"mp2",kind:"vivo.ri-mp2")],execute:{ cfg,inputs,budget in
                try slots(inputs,["reference"])
                guard cfg == .object([:]) else { throw VivoChemistryError.invalid("RI-MP2 configuration") }
                return ["mp2":try VivoCanonicalJSON.encode(VivoFactorizedMP2.solve(input(VivoRIHartreeFockResult.self,"reference",inputs),budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["reference"]);try slots(outputs,["mp2"])
                guard cfg == .object([:]) else { throw VivoChemistryError.invalid("RI-MP2 configuration") }
                let r=try input(VivoMP2Result.self,"mp2",outputs),expected=try VivoFactorizedMP2.solve(input(VivoRIHartreeFockResult.self,"reference",inputs),budget:budget)
                guard r.totalEnergyHartree.isFinite,abs(r.totalEnergyHartree-expected.totalEnergyHartree)<1e-8,
                      abs(r.referenceEnergyHartree-expected.referenceEnergyHartree)<1e-8 else { throw VivoChemistryError.invalid("RI-MP2 result binding") }
            })
    }
    public static func smoothSCF(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.smooth-cpcm-scf",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"scf",kind:"vivo.smooth-cpcm-scf")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs),settings=try decode(VivoSmoothSCFWorkflowConfiguration.self,cfg)
                return ["scf":try VivoCanonicalJSON.encode(VivoSmoothSolventSCF.solve(system:ao.sourceSystem,integrals:ao,solvent:settings.solvent,scf:settings.scf,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals"]);try slots(outputs,["scf"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs),r=try input(VivoSmoothSolventSCFResult.self,"scf",outputs)
                let settings=try decode(VivoSmoothSCFWorkflowConfiguration.self,cfg),hf=r.scf
                try ao.validate(budget:budget);try settings.scf.validate()
                guard hf.reference==settings.scf.reference,hf.alphaElectrons==ao.sourceSystem.alphaElectrons,
                      hf.betaElectrons==ao.sourceSystem.betaElectrons,hf.alphaElectrons<=ao.count,hf.betaElectrons<=ao.count,
                      hf.finalCommutatorNorm.isFinite,hf.finalCommutatorNorm<=settings.scf.commutatorTolerance,
                      r.solvent.configuration==settings.solvent else { throw VivoChemistryError.invalid("smooth solvent SCF sector/settings") }
                for c in [hf.alphaCoefficients,hf.betaCoefficients] {
                    guard c.rows==ao.count,c.columns==ao.count,c.values.allSatisfy(\.isFinite),
                          try ao.overlap.congruence(c).adding(.identity(ao.count),scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("smooth SCF orbital metric") }
                }
                let da=VivoHartreeFock.density(hf.alphaCoefficients,occupied:hf.alphaElectrons),db=VivoHartreeFock.density(hf.betaCoefficients,occupied:hf.betaElectrons)
                guard try da.adding(hf.alphaDensity,scale:-1).frobeniusNorm<=settings.scf.densityTolerance,
                      try db.adding(hf.betaDensity,scale:-1).frobeniusNorm<=settings.scf.densityTolerance else { throw VivoChemistryError.invalid("smooth SCF density/orbital binding") }
                let op=try VivoSmoothCPCMOperator(system:ao.sourceSystem,basis:ao.sourceBasis,configuration:settings.solvent,budget:budget)
                let pcm=try op.evaluate(totalDensity:da.adding(db)),gas=VivoHartreeFock.focks(ao,da,db)
                let energy=VivoHartreeFock.energy(ao,da,db,gas.0,gas.1)+pcm.polarizationEnergyHartree
                let fa=try gas.0.adding(pcm.reactionPotentialMatrix),fb=try gas.1.adding(pcm.reactionPotentialMatrix)
                let errors=try VivoHartreeFock.error(fa,da,ao.overlap,hf.alphaCoefficients)+VivoHartreeFock.error(fb,db,ao.overlap,hf.betaCoefficients)
                guard hf.energyHartree.isFinite,abs(energy-hf.energyHartree)<max(1e-8,10*settings.scf.energyToleranceHartree),
                      abs(pcm.polarizationEnergyHartree-r.solvent.polarizationEnergyHartree)<1e-8,
                      errors.reduce(0.0,{hypot($0,$1)})<=1.01*settings.scf.commutatorTolerance else { throw VivoChemistryError.invalid("smooth SCF physical energy/residual binding") }
            })
    }
    public static func ecc(implementationFingerprint id:VivoFingerprint)->VivoChemistryOperation {
        .init(identifier:"vivo.native.ecc-selected-moment-cycle",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"embedding",kind:"vivo.ecc-selected-moment-result")],execute:{ cfg,inputs,budget in
                try slots(inputs,["hamiltonian"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs),settings=try decode(VivoECCWorkflowConfiguration.self,cfg)
                return ["embedding":try VivoCanonicalJSON.encode(VivoECCSelfConsistency.solve(h,fragments:settings.fragments,
                    operators:settings.operators,configuration:settings.settings,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["hamiltonian"]);try slots(outputs,["embedding"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs),settings=try decode(VivoECCWorkflowConfiguration.self,cfg)
                let r=try input(VivoECCSelfConsistencyResult.self,"embedding",outputs)
                try h.validate(budget:budget)
                guard r.converged,r.numberMatching.converged,r.correlationPotentialHartree.count==settings.operators.count,
                      r.correlationPotentialHartree.allSatisfy({$0.isFinite && abs($0)<=settings.settings.maximumPotentialHartree}),
                      r.momentResiduals.count==settings.operators.count,r.momentResiduals.allSatisfy(\.isFinite),
                      r.referenceMoments.count==settings.operators.count,r.fragmentMoments.count==settings.operators.count,
                      r.momentResiduals.reduce(0.0,{hypot($0,$1)})<=settings.settings.momentTolerance,
                      abs(r.numberMatching.populationResidual)<=settings.settings.numberMatching.populationTolerance,
                      r.numberMatching.targetFragmentPopulationSum==Double(h.alphaElectrons+h.betaElectrons),
                      r.numberMatching.states.map(\.identifier)==settings.fragments.map(\.identifier),
                      r.clusterCoefficients.count==settings.fragments.count else { throw VivoChemistryError.convergence("ECC moment/number feedback did not satisfy its declared contract") }
                for i in r.momentResiduals.indices {
                    guard r.referenceMoments[i].isFinite,r.fragmentMoments[i].isFinite,
                          abs(r.momentResiduals[i]-r.referenceMoments[i]+r.fragmentMoments[i])<1e-10 else { throw VivoChemistryError.invalid("ECC reported moment residual identity") }
                }
                for state in r.numberMatching.states { try state.biasedCI.state.validate(budget:budget) }
            })
    }
}
