import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
@testable import NumiVivoNumerics
#endif

@main struct ECCIntegrationChecks {
    static func synthetic() -> VivoEmbeddedHamiltonian {
        let n=4,diagonal=[-2.0,-0.9,0.2,0.8]
        var h=VivoQMMatrix(n,n),g=[Double](repeating:0,count:256)
        for p in 0..<n{for q in 0..<n{
            h[p,q]=p==q ? diagonal[p]:0.04*cos(Double(p+q+1))
            for r in 0..<n{for s in 0..<n{for k in 0..<3{
                g[((p*n+q)*n+r)*n+s]+=0.0225*cos(Double((p+1)*(q+1)*(k+1)))*cos(Double((r+1)*(s+1)*(k+1)))
            }}}
        }}
        return .init(orbitalIdentifiers:(0..<n).map{"synthetic-\($0)"},alphaElectrons:2,betaElectrons:2,
            oneElectron:h,twoElectron:g,constantEnergyHartree:0.42,energyReference:"synthetic regression; not a molecule")
    }
    static func main() throws {
        let directory=URL(fileURLWithPath:CommandLine.arguments.count>1 ? CommandLine.arguments[1]:"/tmp/ecc-checks")
        try FileManager.default.createDirectory(at:directory,withIntermediateDirectories:true)
        var checks=0,observations:[String:Double]=[:]
        func require(_ condition:Bool,_ label:String) throws {
            guard condition else{throw VivoChemistryError.invalid("regression: \(label)")};checks+=1;print("PASS \(label)")
        }
        func reject(_ label:String,_ operation:()throws->Void) throws {
            var rejected=false;do{try operation()}catch is VivoChemistryError{rejected=true};try require(rejected,label)
        }
        func write<T:Encodable>(_ value:T,_ name:String) throws {
            let encoder=JSONEncoder();encoder.outputFormatting=[.sortedKeys,.prettyPrinted];try encoder.encode(value).write(to:directory.appendingPathComponent(name),options:.atomic)
        }
        let budget=VivoChemistryBudget(maximumBytes:512*1024*1024,maximumBasisFunctions:64,maximumDeterminants:10000,maximumOperatorApplications:1_000_000_000)
        let h=synthetic(),reference=try VivoDirectCI.solve(h,budget:budget).roots[0]
        var k=VivoQMMatrix(4,4);k[0,2]=0.27;k[2,0] = -0.27;k[1,3] = -0.19;k[3,1]=0.19;k[0,3]=0.12;k[3,0] = -0.12
        let u=try VivoQMDenseAlgebra.orbitalRotation(generator:k)
        let rotated=try VivoCIOrbitalFrame.rotated(reference.state,by:u,budget:budget)
        let re=try VivoCIDensityMatrices.compute(rotated,budget:budget).energy(of:h.rotated(by:u,budget:budget),budget:budget)
        observations["frameEnergyError"]=abs(re-reference.energyHartree)
        try require(abs(re-reference.energyHartree)<1e-9,"exact CI frame preserves physical energy")
        let restored=try VivoCIOrbitalFrame.rotated(rotated,by:u.transposed,budget:budget)
        let overlap=zip(reference.state.coefficients,restored.coefficients).reduce(0.0){$0+$1.0*$1.1}
        try require(abs(overlap-1)<1e-10,"CI frame inverse retains fermionic phases")
        var permutation=VivoQMMatrix(4,4);for (new,old) in [3,1,0,2].enumerated(){permutation[old,new]=1}
        let permuted=try VivoCIOrbitalFrame.rotated(reference.state,by:permutation,budget:budget)
        let pe=try VivoCIDensityMatrices.compute(permuted,budget:budget).energy(of:h.rotated(by:permutation,budget:budget),budget:budget)
        try require(abs(pe-reference.energyHartree)<1e-9,"odd occupied-orbital permutations retain spin-order signs")
        let single=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[.init(identifier:"center",orbitals:[0,1],maximumBathOrbitals:2,clusterAlphaElectrons:2,clusterBetaElectrons:2)],bathSelection:.init(minimumBathOrbitals:2))
        let result=try VivoECCDMET.solve(h,configuration:single,budget:budget)
        try require(result.converged,"single-fragment complete bath converges without property matching")
        try require(result.matching==nil && result.frame.numberMatching==nil,"single fragment does not invent a chemical-potential cycle")
        try require(abs(result.energyHartree-reference.energyHartree)<1e-9,"full cluster recovers physical FCI with scalar once")
        try VivoECCDMET.validate(result,hamiltonian:h,configuration:single,budget:budget)
        try require(result.frame.energyContributions.count==1,"single-fragment physical energy decomposition exists")
        let fragments=[VivoECCFragment(identifier:"a",orbitals:[0,1],maximumBathOrbitals:2,clusterAlphaElectrons:2,clusterBetaElectrons:2),
                       VivoECCFragment(identifier:"b",orbitals:[2,3],maximumBathOrbitals:2,clusterAlphaElectrons:2,clusterBetaElectrons:2)]
        var one=VivoQMMatrix(4,4),two=[Double](repeating:0,count:256)
        one[0,0]=1;one[1,1] = -1
        // A genuine local two-body observable, distinct from the fitted 1-RDM.
        two[0]=1;two[((1*4+1)*4+1)*4+1] = -1
        let ops=[VivoECCCorrelationOperator(identifier:"density-contrast",one:one,two:[Double](repeating:0,count:256)),
                 VivoECCCorrelationOperator(identifier:"double-occupancy-contrast",one:VivoQMMatrix(4,4),two:two)]
        let config=VivoECCDMETConfiguration(mode:.selfConsistentPartition,fragments:fragments,correlationOperators:ops,
            matching:.init(referenceMethod:.fci,maximumIterations:60,maximumEvaluations:1000,momentTolerance:1e-8),bathSelection:.init(minimumBathOrbitals:2))
        let cycle=try VivoECCDMET.solve(h,configuration:config,initialPotentialHartree:[0.08,-0.04],budget:budget)
        try require(cycle.converged,"one/two-body correlation feedback converges from nonzero potential")
        try require(cycle.matching!.iterations.count>1,"self-consistency performs actual feedback iterations")
        try require(cycle.matching!.momentResiduals.reduce(0.0,{hypot($0,$1)})<1e-8,"both selected moments meet residual criterion")
        try require(abs(cycle.frame.electronCountDefect)<1e-7,"complete partition particle number closes")
        try require(abs(cycle.energyHartree-reference.energyHartree)<1e-8,"democratic energy recovers exact full-bath limit without impurity double counting")
        try VivoECCDMET.validate(cycle,hamiltonian:h,configuration:config,budget:budget)
        observations["partitionMomentNorm"]=cycle.matching!.momentResiduals.reduce(0.0,{hypot($0,$1)})
        observations["partitionEnergyError"]=abs(cycle.energyHartree-reference.energyHartree)
        observations["matchingIterations"]=Double(cycle.matching!.iterations.count)
        try write(cycle,"partition-result.json");try write(h,"hamiltonian.json");try write(config,"partition-configuration.json")
        // A reduced quantum cluster with a genuinely occupied inactive orbital.
        // The reference is truncated CISD; the impurity FCI correction must
        // restore omitted active-space correlation without recounting the core.
        var embeddedOne=VivoQMMatrix(5,5),embeddedTwo=[Double](repeating:0,count:625)
        embeddedOne[0,0] = -4
        for p in 0..<4{for q in 0..<4{
            embeddedOne[p+1,q+1]=h.oneElectron[p,q]
            for r in 0..<4{for t in 0..<4{embeddedTwo[(((p+1)*5+q+1)*5+r+1)*5+t+1]=h.eri(p,q,r,t)}}
        }}
        for a in 1..<5{let value=0.07/Double(a+1);embeddedTwo[a*5+a]=value;embeddedTwo[(a*5+a)*25]=value}
        let envH=VivoEmbeddedHamiltonian(orbitalIdentifiers:(0..<5).map{"env-\($0)"},alphaElectrons:3,betaElectrons:3,
            oneElectron:embeddedOne,twoElectron:embeddedTwo,constantEnergyHartree:0.37,energyReference:"synthetic four-orbital correlated center plus occupied inactive core")
        let reducedCfg=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[.init(identifier:"reactive",orbitals:[1,2],maximumBathOrbitals:2,clusterAlphaElectrons:2,clusterBetaElectrons:2)],
            matching:.init(referenceMethod:.cisd),bathSelection:.init(minimumBathOrbitals:2))
        let reduced=try VivoECCDMET.solve(envH,configuration:reducedCfg,budget:budget)
        let exactEnv=try VivoDirectCI.solve(envH,budget:budget).roots[0]
        try require(reduced.converged && reduced.frame.clusters[0].coefficients.columns==4,"correlated impurity executes with a nonempty inactive environment")
        try require(abs(reduced.energyHartree-exactEnv.energyHartree)<1e-8,"inactive reference restoration recovers separable-core exact limit")
        try require(abs(reduced.frame.energyContributions[0].correctionHartree)>1e-10,"FCI impurity adds correlation missing from the CISD reference")
        observations["reducedClusterEnergyError"]=abs(reduced.energyHartree-exactEnv.energyHartree)
        observations["reducedClusterCorrection"]=reduced.frame.energyContributions[0].correctionHartree
        try write(envH,"environment-hamiltonian.json");try write(reducedCfg,"environment-configuration.json");try write(reduced,"environment-result.json")
        let system=VivoElectronicSystem(nuclei:[.init(atomicNumber:1,positionBohr:.init(0,0,-0.7)),.init(atomicNumber:1,positionBohr:.init(0,0,0.7))],alphaElectrons:1,betaElectrons:1)
        let ao=try VivoGaussianIntegralEngine.compute(system:system,basis:.hydrogenSTO3G(nucleusIndices:[0,1]))
        let hf=try VivoHartreeFock.solve(system:system,integrals:ao)
        let molecule=try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:hf.alphaCoefficients,alphaElectrons:1,betaElectrons:1,orbitalIdentifiers:["0","1"],energyReference:"H2 STO3G 1.4 Bohr")
        var gen=VivoQMMatrix(2,2);gen[0,1]=0.21;gen[1,0] = -0.21
        let start=try VivoQMDenseAlgebra.orbitalRotation(generator:gen)
        let orbitals=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[.init(identifier:"f",orbitals:[0],maximumBathOrbitals:1,clusterAlphaElectrons:1,clusterBetaElectrons:1)],
            bathSelection:.init(minimumBathOrbitals:1),localityGroups:[[0,1]],maximumOrbitalIterations:100,maximumFrameEvaluations:1000,
            orbitalGradientTolerance:1e-5,maximumOrbitalStep:0.15)
        let optimized=try VivoECCDMET.solve(molecule,configuration:orbitals,initialRotation:start,budget:budget)
        try require(optimized.converged,"locality-constrained QIO converges with all stages recomputed")
        try require(optimized.iterations.count>1 && optimized.iterations.last!.leakageObjective<optimized.iterations.first!.leakageObjective-1e-4,"QIO performs nontrivial objective-reducing orbital updates")
        try VivoECCDMET.validate(optimized,hamiltonian:molecule,configuration:orbitals,budget:budget)
        observations["qioGradientNorm"]=optimized.orbitalGradient.reduce(0.0,{hypot($0,$1)})
        observations["qioInitial"]=optimized.iterations.first!.leakageObjective;observations["qioFinal"]=optimized.iterations.last!.leakageObjective
        try write(optimized,"qio-result.json")
        try reject("configured bath cannot silently discard significant density coupling") {
            _=try VivoECCDMET.solve(h,configuration:.init(mode:.singleFragment,fragments:[.init(identifier:"small",orbitals:[0],maximumBathOrbitals:0,clusterAlphaElectrons:1,clusterBetaElectrons:1)]),budget:budget)
        }
        // Truncating a correlated reference can leave a noninteger inactive
        // population; a fixed-electron impurity is not silently rounded to fit.
        let clipped=VivoECCDMETConfiguration(mode:.singleFragment,fragments:[.init(identifier:"small",orbitals:[0],maximumBathOrbitals:0,clusterAlphaElectrons:1,clusterBetaElectrons:1)],
            matching:.init(referenceMethod:.fci,bathDiscardedWeight:100),electronTolerance:1e-12)
        let mismatch=try VivoECCDMET.solve(h,configuration:clipped,budget:budget)
        try require(!mismatch.converged && mismatch.termination=="inactive-particle-closure-failed","fractional inactive populations do not receive a success certificate")
        try reject("nonconverged ECC result cannot be reused") {try VivoECCDMET.validate(mismatch,hamiltonian:h,configuration:clipped,budget:budget)}
        try reject("self-consistency requires a complete disjoint fragment partition") {
            var bad=config;bad.fragments=Array(fragments.prefix(1));_ = try VivoECCDMET.solve(h,configuration:bad,budget:budget)
        }
        try reject("fixed locality classes reject a rotating restart") {
            _=try VivoECCDMET.solve(h,configuration:single,initialRotation:u,budget:budget)
        }
        var earlyCfg=orbitals;earlyCfg.maximumOrbitalIterations=1
        let early=try VivoECCDMET.solve(molecule,configuration:earlyCfg,initialRotation:start,budget:budget)
        try require(!early.converged,"short QIO execution preserves nonconvergence")
        let encoder=JSONEncoder()
        var forged=try JSONSerialization.jsonObject(with:encoder.encode(early)) as! [String:Any]
        forged["converged"]=true;forged["orbitalGradient"]=[0.0]
        let falseStationarity=try JSONDecoder().decode(VivoECCDMETResult.self,from:JSONSerialization.data(withJSONObject:forged))
        try reject("cache validation recomputes QIO gradient instead of trusting a claimed zero") {
            try VivoECCDMET.validate(falseStationarity,hamiltonian:molecule,configuration:earlyCfg,budget:budget)
        }
        var forgedMoment=try JSONSerialization.jsonObject(with:encoder.encode(cycle)) as! [String:Any]
        forgedMoment["correlationPotentialHartree"]=[0.08,-0.04]
        let falsePotential=try JSONDecoder().decode(VivoECCDMETResult.self,from:JSONSerialization.data(withJSONObject:forgedMoment))
        try reject("cache validation does not silently repair an unconverged correlation potential") {
            try VivoECCDMET.validate(falsePotential,hamiltonian:h,configuration:config,budget:budget)
        }
        struct Report:Codable {let checksPassed:Int;let observations:[String:Double];let scope:String}
        try write(Report(checksPassed:checks,observations:observations,scope:"native finite-basis ECC-DMET algebraic and feedback checks; full-bath exact limits plus locality QIO; not biological validation"),"results.json")
        print("PASS \(checks) ECC integration checks")
    }
}
