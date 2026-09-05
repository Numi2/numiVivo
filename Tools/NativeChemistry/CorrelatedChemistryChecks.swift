import Foundation

/// Portable numerical checks compiled with the production numerical sources.
/// No artifact-store, Metal or provider substitute is introduced by this runner.
@main struct CorrelatedChemistryChecks {
    private struct Observation: Codable {
        let label: String
        let value: Double
        let unit: String
    }
    private struct Report: Codable {
        let schema: String
        let checksPassed: Int
        let backend: String
        let observations: [Observation]
        let limitations: [String]
    }
    private struct Checks {
        var count = 0
        var observations: [Observation] = []
        mutating func require(_ condition: Bool, _ label: String) throws {
            guard condition else { throw VivoChemistryError.invalid("regression failed: \(label)") }
            count += 1
            print("PASS \(label)")
        }
        mutating func observe(_ label: String, _ value: Double, unit: String = "dimensionless") throws {
            guard value.isFinite else { throw VivoChemistryError.invalid("nonfinite observation \(label)") }
            observations.append(.init(label:label,value:value,unit:unit))
            print("VALUE \(label) = \(value) \(unit)")
        }
        mutating func rejects(_ label: String, _ operation: () throws -> Void) throws {
            var rejected = false
            do { try operation() } catch is VivoChemistryError { rejected = true }
            try require(rejected,label)
        }
    }
    private static func write<T: Encodable>(_ value: T, to path: URL) throws {
        let encoder=JSONEncoder()
        encoder.outputFormatting=[.sortedKeys,.prettyPrinted,.withoutEscapingSlashes]
        encoder.nonConformingFloatEncodingStrategy = .throw
        try encoder.encode(value).write(to:path,options:.atomic)
    }
    private static func syntheticHamiltonian() -> VivoEmbeddedHamiltonian {
        let n=4, diagonal=[-2.0,-0.9,0.2,0.8]
        var one=VivoQMMatrix(n,n), two=[Double](repeating:0,count:n*n*n*n)
        for p in 0..<n { for q in 0..<n {
            one[p,q] = p == q ? diagonal[p] : 0.04*cos(Double(p+q+1))
            for r in 0..<n { for s in 0..<n {
                var value=0.0
                for k in 0..<3 {
                    let a=0.15*cos(Double((p+1)*(q+1)*(k+1)))
                    let b=0.15*cos(Double((r+1)*(s+1)*(k+1)))
                    value += a*b
                }
                two[((p*n+q)*n+r)*n+s]=value
            } }
        } }
        return .init(orbitalIdentifiers:(0..<n).map { "synthetic-\($0)" },alphaElectrons:2,betaElectrons:2,
            oneElectron:one,twoElectron:two,constantEnergyHartree:0.42,
            energyReference:"synthetic four-orbital regression; not a molecular reaction",
            provenance:["fixture":"three cosine-factor ERI outer products"])
    }
    static func main() throws {
        guard CommandLine.arguments.count==2 else { throw VivoChemistryError.invalid("expected output directory") }
        let output=URL(fileURLWithPath:CommandLine.arguments[1],isDirectory:true)
        try FileManager.default.createDirectory(at:output,withIntermediateDirectories:true)
        var checks=Checks()
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),
            .init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let ao=try VivoGaussianIntegralEngine.compute(system:system,basis:.hydrogenSTO3G(nucleusIndices:[0,1]))
        let hf=try VivoHartreeFock.solve(system:system,integrals:ao)
        let h2=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,alphaElectrons:1,betaElectrons:1,
            orbitalIdentifiers:["h2-0","h2-1"],energyReference:"isolated H2; distance 1.4 Bohr; STO-3G")
        let cc2=try VivoCoupledCluster.solve(h2), fci2=try VivoConfigurationInteraction.solve(h2)
        let cc2RDM=try VivoCoupledCluster.responseDensityMatrices(cc2)
        try checks.require(cc2.converged,"H2 CCSD converges")
        try checks.require(abs(cc2.energyHartree-fci2.energyHartree)<1e-9,"two-electron CCSD equals FCI")
        try checks.require(abs(try cc2RDM.energy(of:h2)-cc2.energyHartree)<1e-9,"H2 Lambda RDM energy")
        try checks.require(abs(cc2.lambda!.biorthogonalOverlap-1)<1e-10,"CC left/right overlap")
        try checks.observe("h2CCSDEnergy",cc2.energyHartree,unit:"Hartree")
        try checks.observe("h2FCIEnergy",fci2.energyHartree,unit:"Hartree")
        try checks.observe("h2CCSDResidual",cc2.projectedResidualNorm,unit:"Hartree")
        let dispatch=try VivoManyBodySolver.solve(h2,request:.ccsd(configuration:.init()))
        try VivoManyBodySolver.validate(dispatch,hamiltonian:h2,request:.ccsd(configuration:.init()))
        try checks.require(dispatch.informationState==nil,"CC response is not exposed as a QIO probability state")
        let dispatchedFCI=try VivoManyBodySolver.solve(h2,request:.configurationInteraction(method:.fci))
        try VivoManyBodySolver.validate(dispatchedFCI,hamiltonian:h2,request:.configurationInteraction(method:.fci))
        try checks.require(dispatchedFCI.informationState != nil,"FCI exposes a normalized information state")
        try checks.rejects("method mismatch is rejected") {
            try VivoManyBodySolver.validate(dispatch,hamiltonian:h2,request:.configurationInteraction(method:.fci))
        }
        try checks.rejects("CCSD work budget is enforced") {
            _ = try VivoCoupledCluster.solve(h2,budget:.init(maximumOperatorApplications:1))
        }
        try checks.rejects("CCSD rejects a reference in a different particle sector") {
            _ = try VivoCoupledCluster.solve(h2,referenceDeterminant:0)
        }
        let short=try VivoCoupledCluster.solve(h2,configuration:.init(maximumIterations:1))
        try checks.require(!short.converged,"iteration-limited CCSD does not claim convergence")
        try checks.rejects("unsolved CCSD is rejected by common result validation") {
            try VivoManyBodySolver.validate(.ccsd(result:short),hamiltonian:h2,
                request:.ccsd(configuration:.init(maximumIterations:1)))
        }
        let h=syntheticHamiltonian(), partition=VivoActiveSpace(doublyOccupiedCore:[0],active:[1,2])
        let active=try h.frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore)
        let ci=try VivoConfigurationInteraction.solve(active)
        let pairs:[VivoOrbitalRotationPair] = [.init(first:0,second:1),.init(first:0,second:2),.init(first:0,second:3),
            .init(first:1,second:3),.init(first:2,second:3)]
        let gradient=try VivoCASOrbitalGradient.compute(hamiltonian:h,partition:partition,state:ci.state,pairs:pairs)
        var maximumGradientError=0.0
        for (i,pair) in pairs.enumerated() {
            func energy(_ step:Double) throws -> Double {
                var k=VivoQMMatrix(4,4); k[pair.first,pair.second]=step; k[pair.second,pair.first] = -step
                let u=try VivoQMDenseAlgebra.orbitalRotation(generator:k)
                let rotated=try h.rotated(by:u).frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore)
                return try VivoConfigurationInteraction.solve(rotated).energyHartree
            }
            let finiteDifference=try (energy(1e-5)-energy(-1e-5))/2e-5
            let error=abs(finiteDifference-gradient[i]); maximumGradientError=max(maximumGradientError,error)
            try checks.require(error<1e-7,"CAS orbital gradient pair \(pair.first),\(pair.second)")
        }
        let cas=try VivoCASSCF.solve(h,partition:partition)
        try checks.require(cas.converged,"nontrivial core/active/external CASSCF converges")
        try checks.require(cas.activeCI.energyHartree<ci.energyHartree-1e-6,"CASSCF improves fixed-orbital CAS energy")
        try VivoManyBodySolver.validate(.casscf(result:cas),hamiltonian:h,request:.casscf(partition:partition,configuration:.init()))
        try checks.require(cas.orbitalGradient.reduce(0.0,{hypot($0,$1)})<1e-6,"CASSCF residual meets requested tolerance")
        try checks.observe("casMaximumGradientError",maximumGradientError,unit:"Hartree/radian")
        try checks.observe("initialCASCI",ci.energyHartree,unit:"Hartree")
        try checks.observe("optimizedCASSCF",cas.activeCI.energyHartree,unit:"Hartree")
        try checks.observe("casscfGradientNorm",cas.orbitalGradient.reduce(0.0,{hypot($0,$1)}),unit:"Hartree/radian")
        let cc=try VivoCoupledCluster.solve(h), fci=try VivoConfigurationInteraction.solve(h)
        let response=try VivoCoupledCluster.responseDensityMatrices(cc)
        try checks.require(cc.converged,"four-electron CCSD converges")
        try checks.require(abs(try response.energy(of:h)-cc.energyHartree)<1e-8,"four-electron Lambda RDM energy")
        try checks.require(abs(cc.energyHartree-fci.energyHartree)>1e-10,"four-electron CCSD is not replaced by FCI")
        func perturbedEnergy(_ value:Double) throws -> Double {
            var one=h.oneElectron; one[0,1]+=value; one[1,0]+=value
            let perturbed=VivoEmbeddedHamiltonian(orbitalIdentifiers:h.orbitalIdentifiers,alphaElectrons:2,betaElectrons:2,
                oneElectron:one,twoElectron:h.twoElectron,constantEnergyHartree:h.constantEnergyHartree,energyReference:h.energyReference)
            let result=try VivoCoupledCluster.solve(perturbed)
            guard result.converged else { throw VivoChemistryError.convergence("perturbed response check") }
            return result.energyHartree
        }
        let responseDerivative=response.one[0,2]+response.one[2,0]+response.one[1,3]+response.one[3,1]
        let finiteDifference=try (perturbedEnergy(1e-5)-perturbedEnergy(-1e-5))/2e-5
        try checks.require(abs(responseDerivative-finiteDifference)<1e-7,"CC Lambda response matches an energy derivative")
        try checks.observe("fourElectronCCSD",cc.energyHartree,unit:"Hartree")
        try checks.observe("fourElectronFCI",fci.energyHartree,unit:"Hartree")
        try checks.observe("fourElectronCCSDResidual",cc.projectedResidualNorm,unit:"Hartree")
        try checks.observe("ccResponseDerivativeError",abs(responseDerivative-finiteDifference),unit:"dimensionless")
        let fragmentH=VivoEmbeddedHamiltonian(orbitalIdentifiers:["fragment","bath"],alphaElectrons:1,betaElectrons:1,
            oneElectron:try .init(rows:2,columns:2,values:[-1,-0.2,-0.2,0.3]),twoElectron:[Double](repeating:0,count:16),
            constantEnergyHartree:0,energyReference:"synthetic two-level population check")
        let projector=try VivoQMMatrix(rows:2,columns:2,values:[1,0,0,0])
        let matched=try VivoFragmentNumberMatcher.solve(fragments:[.init(identifier:"fragment",hamiltonian:fragmentH,fragmentProjector:projector)],targetPopulation:1)
        try checks.require(matched.converged && abs(matched.chemicalPotentialHartree+1.3)<1e-7,"fragment-only chemical potential matches population")
        let invariant=try VivoFragmentNumberMatcher.solve(fragments:[.init(identifier:"whole",hamiltonian:fragmentH,fragmentProjector:.identity(2))],targetPopulation:1)
        try checks.require(!invariant.converged,"uniform shift does not alter fixed total particle number")
        try checks.observe("matchedChemicalPotential",matched.chemicalPotentialHartree,unit:"Hartree")
        try checks.observe("matchedPopulationResidual",matched.populationResidual,unit:"electrons")
        let gas=try VivoRestrictedLDA.solve(system:system,integrals:ao)
        let solvated=try VivoRestrictedLDA.solve(system:system,integrals:ao,configuration:.init(cpcm:.init()))
        let rhf=try VivoCPCMRHF.solve(system:system,integrals:ao)
        let encoded=[try JSONEncoder().encode(gas),try JSONEncoder().encode(solvated),try JSONEncoder().encode(rhf)]
        try checks.require(encoded.allSatisfy({!$0.isEmpty}),"SCF traces encode without nonfinite JSON sentinels")
        try checks.require(abs(solvated.energyHartree-solvated.gasEnergyFunctionalHartree!-solvated.cpcm!.polarizationEnergyHartree)<1e-12,"C-PCM polarization counted once")
        try checks.require(solvated.integratedElectronError<5e-4 && solvated.finalCommutatorNorm<1e-7,"solvated LDA grid and SCF acceptance")
        try checks.require(solvated.cpcm!.linearRelativeResidual<8e-10,"C-PCM true linear residual")
        try checks.observe("gasLDA",gas.energyHartree,unit:"Hartree")
        try checks.observe("integratedLDAElectrons",gas.integratedElectrons,unit:"electrons")
        try checks.observe("solvatedLDA",solvated.energyHartree,unit:"Hartree")
        try checks.observe("solvatedLDAPolarization",solvated.cpcm!.polarizationEnergyHartree,unit:"Hartree")
        try checks.observe("solvatedLDASCFResidual",solvated.finalCommutatorNorm,unit:"Hartree")
        try checks.observe("cpcmTrueRelativeResidual",solvated.cpcm!.linearRelativeResidual)
        try checks.observe("solvatedRHF",rhf.totalFreeEnergyHartree,unit:"Hartree")
        try write(h,to:output.appendingPathComponent("synthetic-hamiltonian.json"))
        try write(VivoManyBodySolverRequest.ccsd(configuration:.init()),to:output.appendingPathComponent("solver-ccsd.json"))
        try write(VivoManyBodySolverRequest.casscf(partition:partition,configuration:.init()),to:output.appendingPathComponent("solver-casscf.json"))
        try write(cc,to:output.appendingPathComponent("ccsd-result.json"))
        try write(cas,to:output.appendingPathComponent("casscf-result.json"))
        try write(solvated,to:output.appendingPathComponent("lda-cpcm-result.json"))
        #if canImport(Accelerate)
        let backend="CPU FP64; Accelerate matrix multiplication"
        #else
        let backend="CPU FP64; portable matrix multiplication"
        #endif
        let report=Report(schema:"numivivo.org/correlated-regression/v1",checksPassed:checks.count,backend:backend,
            observations:checks.observations,limitations:[
                "H2/STO-3G and explicitly synthetic Hamiltonians only; no paper reaction or BTK calculation.",
                "Internal numerical/derivative consistency checks, not an independent PySCF/ORCA conformance suite.",
                "This executable excludes the artifact store, CLI integration and Metal; it does not qualify the full package.",
                "CCSD uses a bounded determinant evaluation sector; CASSCF uses small active-space FCI.",
                "Continuum-solvent output is electrostatic polarization, not a complete Gibbs free energy."
            ])
        try write(report,to:output.appendingPathComponent("results.json"))
        print("PASS \(checks.count) checks; backend: \(backend)")
    }
}
