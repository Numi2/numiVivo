import Foundation

/// Native operations over the existing verified artifact DAG. Geometry/basis
/// integrals are a separate reusable node: changing CCSD to CASCI need not
/// recompute the AO integrals or converged RHF reference.
public enum VivoElectronicWorkflowOperations {
    private static func decode<T:Decodable>(_ type:T.Type,_ value:VivoJSONValue) throws -> T {
        try VivoCanonicalJSON.decode(type,from:VivoCanonicalJSON.encode(value))
    }
    private static func input<T:Decodable>(_ type:T.Type,_ name:String,_ values:[String:Data]) throws -> T {
        guard let data=values[name] else { throw VivoChemistryError.invalid("missing native chemistry slot \(name)") }
        return try VivoCanonicalJSON.decode(type,from:data)
    }
    private static func slots(_ values:[String:Data],_ names:[String]) throws {
        guard Set(values.keys)==Set(names) else { throw VivoChemistryError.invalid("native chemistry input/output slot mismatch") }
    }
    public static func integrals(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.ao-integrals",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"integrals",kind:"vivo.ao-integrals")],execute:{ cfg,inputs,budget in
                try slots(inputs,["system","basis"])
                guard cfg == .object([:]) else { throw VivoChemistryError.invalid("integral operation accepts no unrecognized configuration") }
                let s=try input(VivoElectronicSystem.self,"system",inputs), b=try input(VivoGaussianBasis.self,"basis",inputs)
                return ["integrals":try VivoCanonicalJSON.encode(VivoGaussianIntegralEngine.compute(system:s,basis:b,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["system","basis"]); try slots(outputs,["integrals"])
                guard cfg == .object([:]) else { throw VivoChemistryError.invalid("integral configuration") }
                let ao=try input(VivoAOIntegrals.self,"integrals",outputs)
                try ao.validate(budget:budget)
                guard ao.sourceSystem == (try input(VivoElectronicSystem.self,"system",inputs)),
                      ao.sourceBasis == (try input(VivoGaussianBasis.self,"basis",inputs)) else { throw VivoChemistryError.invalid("AO integral input binding") }
            })
    }
    public static func hartreeFock(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.hf-from-integrals",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"scf",kind:"vivo.hartree-fock")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), settings=try decode(VivoSCFConfiguration.self,cfg)
                return ["scf":try VivoCanonicalJSON.encode(VivoHartreeFock.solve(system:ao.sourceSystem,integrals:ao,configuration:settings,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals"]); try slots(outputs,["scf"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), hf=try input(VivoHartreeFockResult.self,"scf",outputs)
                try VivoHartreeFock.validate(result:hf,system:ao.sourceSystem,integrals:ao,configuration:decode(VivoSCFConfiguration.self,cfg),budget:budget)
            })
    }
    public static func hamiltonian(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.rhf-hamiltonian",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"hamiltonian",kind:"vivo.embedded-hamiltonian")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals","scf"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), hf=try input(VivoHartreeFockResult.self,"scf",inputs)
                let settings=try decode(VivoSCFConfiguration.self,cfg)
                guard settings.reference == .restricted else { throw VivoChemistryError.unsupported("RHF spatial Hamiltonian preparation requires restricted orbitals") }
                try VivoHartreeFock.validate(result:hf,system:ao.sourceSystem,integrals:ao,configuration:settings,budget:budget)
                let h=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,alphaElectrons:hf.alphaElectrons,
                    betaElectrons:hf.betaElectrons,orbitalIdentifiers:(0..<ao.count).map { "rhf-orbital-\($0)" },
                    energyReference:"electronic; scalar inherited verbatim from AO integrals; no frozen core or thermal/standard-state correction",budget:budget)
                var provenance=h.provenance
                provenance["aoArtifactSHA256"]=try VivoCanonicalJSON.fingerprint(inputs["integrals"]!).hex
                provenance["referenceSHA256"]=try VivoCanonicalJSON.fingerprint(inputs["scf"]!).hex
                let bound=VivoEmbeddedHamiltonian(orbitalIdentifiers:h.orbitalIdentifiers,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
                    oneElectron:h.oneElectron,twoElectron:h.twoElectron,constantEnergyHartree:h.constantEnergyHartree,
                    energyReference:h.energyReference,provenance:provenance)
                return ["hamiltonian":try VivoCanonicalJSON.encode(bound)]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals","scf"]); try slots(outputs,["hamiltonian"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",outputs), ao=try input(VivoAOIntegrals.self,"integrals",inputs)
                let hf=try input(VivoHartreeFockResult.self,"scf",inputs), settings=try decode(VivoSCFConfiguration.self,cfg)
                try h.validate(budget:budget)
                try VivoHartreeFock.validate(result:hf,system:ao.sourceSystem,integrals:ao,configuration:settings,budget:budget)
                guard settings.reference == .restricted, h.alphaElectrons==hf.alphaElectrons, h.betaElectrons==hf.betaElectrons,
                      h.orbitalIdentifiers==(0..<ao.count).map({"rhf-orbital-\($0)"}), h.constantEnergyHartree==ao.constantEnergyHartree,
                      h.provenance["aoArtifactSHA256"] == (try VivoCanonicalJSON.fingerprint(inputs["integrals"]!).hex),
                      h.provenance["referenceSHA256"] == (try VivoCanonicalJSON.fingerprint(inputs["scf"]!).hex),
                      try h.oneElectron.adding(ao.coreHamiltonian.congruence(hf.alphaCoefficients),scale:-1).frobeniusNorm<1e-9 else {
                    throw VivoChemistryError.invalid("embedded Hamiltonian provenance/reference binding")
                }
            })
    }
    public static func manyBody(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.many-body",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"manyBody",kind:"vivo.many-body-result")],execute:{ cfg,inputs,budget in
                try slots(inputs,["hamiltonian"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs), request=try decode(VivoManyBodySolverRequest.self,cfg)
                return ["manyBody":try VivoCanonicalJSON.encode(VivoManyBodySolver.solve(h,request:request,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["hamiltonian"]); try slots(outputs,["manyBody"])
                let h=try input(VivoEmbeddedHamiltonian.self,"hamiltonian",inputs), result=try input(VivoManyBodySolverResult.self,"manyBody",outputs)
                try VivoManyBodySolver.validate(result,hamiltonian:h,request:decode(VivoManyBodySolverRequest.self,cfg),budget:budget)
            })
    }
    private static func checkDensity(_ p:VivoQMMatrix,_ c:VivoQMMatrix,ao:VivoAOIntegrals) throws {
        let n=ao.count
        guard p.rows==n, p.columns==n, c.rows==n, c.columns==n,
              p.values.allSatisfy(\.isFinite), c.values.allSatisfy(\.isFinite),
              try ao.overlap.congruence(c).adding(VivoQMMatrix.identity(n),scale:-1).frobeniusNorm<1e-7,
              try p.adding(VivoHartreeFock.density(c,occupied:ao.sourceSystem.alphaElectrons).scaled(2),scale:-1).frobeniusNorm<1e-7 else {
            throw VivoChemistryError.invalid("mean-field output AO-metric density/coefficients")
        }
    }
    private static func checkPCM(_ actual:VivoCPCMResult,_ expected:VivoCPCMResult) throws {
        guard actual.dielectricConstant==expected.dielectricConstant, actual.tesserae==expected.tesserae,
              actual.apparentSurfaceChargesE.count==expected.apparentSurfaceChargesE.count,
              actual.apparentSurfaceChargesE.allSatisfy(\.isFinite),
              abs(actual.polarizationEnergyHartree-expected.polarizationEnergyHartree)<1e-8,
              zip(actual.apparentSurfaceChargesE,expected.apparentSurfaceChargesE).reduce(0.0,{hypot($0,$1.0-$1.1)})<1e-7,
              try actual.reactionPotentialMatrix.adding(expected.reactionPotentialMatrix,scale:-1).frobeniusNorm<1e-7 else {
            throw VivoChemistryError.invalid("C-PCM output reaction field does not match final density/cavity")
        }
    }
    public static func lda(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.lda-from-integrals",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"lda",kind:"vivo.lda")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), settings=try decode(VivoLDAConfiguration.self,cfg)
                return ["lda":try VivoCanonicalJSON.encode(VivoRestrictedLDA.solve(system:ao.sourceSystem,integrals:ao,configuration:settings,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals"]); try slots(outputs,["lda"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), r=try input(VivoLDAResult.self,"lda",outputs)
                let settings=try decode(VivoLDAConfiguration.self,cfg)
                try ao.validate(budget:budget); try settings.validate(atomCount:ao.sourceSystem.nuclei.count)
                guard ao.sourceSystem.alphaElectrons==ao.sourceSystem.betaElectrons, ao.sourceSystem.alphaElectrons<=ao.count else { throw VivoChemistryError.invalid("restricted LDA output spin sector") }
                try checkDensity(r.totalDensity,r.coefficients,ao:ao)
                let grid=try VivoDFTQuadrature.build(system:ao.sourceSystem,configuration:settings.grid)
                let xc=try VivoLDAWorkspace(ao:ao,grid:grid,budget:budget).evaluate(r.totalDensity)
                let j=VivoRestrictedLDA.coulomb(ao,density:r.totalDensity)
                var fock=try ao.coreHamiltonian.adding(j).adding(xc.matrix)
                let gas=ao.constantEnergyHartree+zip(r.totalDensity.values,ao.coreHamiltonian.values).reduce(0.0) { $0+$1.0*$1.1 }
                    + 0.5*zip(r.totalDensity.values,j.values).reduce(0.0) { $0+$1.0*$1.1 }+xc.energy
                var energy=gas
                if let configuration=settings.cpcm {
                    guard let actual=r.cpcm else { throw VivoChemistryError.invalid("LDA result omits requested solvent") }
                    let expected=try VivoCPCM.evaluate(system:ao.sourceSystem,integrals:ao,totalDensity:r.totalDensity,configuration:configuration,budget:budget)
                    try checkPCM(actual,expected); energy += expected.polarizationEnergyHartree
                    fock=try fock.adding(expected.reactionPotentialMatrix)
                } else if r.cpcm != nil { throw VivoChemistryError.invalid("gas-phase LDA result unexpectedly contains solvent") }
                let residual=try VivoHartreeFock.error(fock,r.totalDensity,ao.overlap,r.coefficients).reduce(0.0) { hypot($0,$1) }
                let tolerance=max(1e-8,10*settings.energyToleranceHartree)
                guard r.energyHartree.isFinite, abs(r.energyHartree-energy)<tolerance,
                      abs(r.exchangeCorrelationEnergyHartree-xc.energy)<tolerance,
                      abs(r.integratedElectrons-xc.electrons)<1e-8,
                      abs(xc.electrons-Double(ao.sourceSystem.alphaElectrons+ao.sourceSystem.betaElectrons))<=settings.maximumIntegratedElectronError,
                      r.gridPointCount==grid.count, residual<=1.01*settings.commutatorTolerance,
                      r.finalCommutatorNorm.isFinite, r.finalCommutatorNorm<=settings.commutatorTolerance else {
                    throw VivoChemistryError.invalid("LDA output energy/grid/SCF binding")
                }
            })
    }
    public static func cpcmRHF(implementationFingerprint id:VivoFingerprint) -> VivoChemistryOperation {
        .init(identifier:"vivo.native.cpcm-rhf-from-integrals",version:"1",implementationFingerprint:id,
            outputs:[.init(name:"scf",kind:"vivo.cpcm-rhf")],execute:{ cfg,inputs,budget in
                try slots(inputs,["integrals"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), settings=try decode(VivoCPCMRHFConfiguration.self,cfg)
                return ["scf":try VivoCanonicalJSON.encode(VivoCPCMRHF.solve(system:ao.sourceSystem,integrals:ao,configuration:settings,budget:budget))]
            },validateOutputs:{ cfg,inputs,outputs,budget in
                try slots(inputs,["integrals"]); try slots(outputs,["scf"])
                let ao=try input(VivoAOIntegrals.self,"integrals",inputs), r=try input(VivoCPCMRHFResult.self,"scf",outputs)
                let settings=try decode(VivoCPCMRHFConfiguration.self,cfg)
                try ao.validate(budget:budget); try settings.validate(system:ao.sourceSystem)
                guard ao.sourceSystem.alphaElectrons<=ao.count else { throw VivoChemistryError.invalid("RHF output occupation") }
                let p=try r.alphaDensity.adding(r.betaDensity)
                try checkDensity(p,r.coefficients,ao:ao)
                guard try r.alphaDensity.adding(r.betaDensity,scale:-1).frobeniusNorm<1e-9 else { throw VivoChemistryError.invalid("RHF output spin densities differ") }
                let pcm=try VivoCPCM.evaluate(system:ao.sourceSystem,integrals:ao,totalDensity:p,configuration:settings.cpcm,budget:budget)
                try checkPCM(r.cpcm,pcm)
                let gas=VivoHartreeFock.focks(ao,r.alphaDensity,r.betaDensity)
                let e=VivoHartreeFock.energy(ao,r.alphaDensity,r.betaDensity,gas.0,gas.1)
                let fa=try gas.0.adding(pcm.reactionPotentialMatrix), fb=try gas.1.adding(pcm.reactionPotentialMatrix)
                let errors=try VivoHartreeFock.error(fa,r.alphaDensity,ao.overlap,r.coefficients)+VivoHartreeFock.error(fb,r.betaDensity,ao.overlap,r.coefficients)
                guard abs(r.totalFreeEnergyHartree-e-pcm.polarizationEnergyHartree)<max(1e-8,10*settings.scf.energyToleranceHartree),
                      errors.reduce(0.0,{hypot($0,$1)})<=1.01*settings.scf.commutatorTolerance,
                      r.finalCommutatorNorm.isFinite, r.finalCommutatorNorm<=settings.scf.commutatorTolerance else {
                    throw VivoChemistryError.invalid("C-PCM RHF output energy/residual binding")
                }
            })
    }
}
