import Foundation
#if NUMIVIVO_SCOPED_NUMERICS
@testable import NumiVivoNumerics
#else
@testable import NumiVivoKit
#endif

@main struct ReactionQualificationChecks {
    struct Water:Decodable {
        let positionsBohr:[SIMD3<Double>]
        let massesDa:[Double]
        let hessian:VivoQMMatrix
        let electronicEnergy:Double
    }
    struct Oracle:Decodable {struct Records:Decodable {let water:Water};let records:Records}
    static func main() throws {
        guard CommandLine.arguments.count==3 else {throw VivoChemistryError.invalid("expected oracle file and empty output directory")}
        let out=URL(fileURLWithPath:CommandLine.arguments[2],isDirectory:true)
        try FileManager.default.createDirectory(at:out,withIntermediateDirectories:true)
        let oracle=try JSONDecoder().decode(Oracle.self,from:Data(contentsOf:URL(fileURLWithPath:CommandLine.arguments[1])))
        var labels:[String]=[],observations:[String:Double]=[:]
        func write<T:Encodable>(_ value:T,_ name:String) throws {
            let e=JSONEncoder();e.outputFormatting=[.prettyPrinted,.sortedKeys]
            try e.encode(value).write(to:out.appendingPathComponent(name),options:.atomic)
        }
        func require(_ condition:Bool,_ label:String) throws {
            guard condition else {throw VivoChemistryError.invalid("reaction qualification check: \(label)")}
            labels.append(label);print("PASS \(label)")
        }
        func reject(_ label:String,_ work:() throws -> Void) throws {
            var rejected=false
            do {try work()} catch is VivoChemistryError {rejected=true} catch is VivoECCPathError {rejected=true}
            try require(rejected,label)
        }
        func mutate<T:Codable>(_ value:T,_ transform:(inout [String:Any])->Void) throws -> T {
            var object=try JSONSerialization.jsonObject(with:JSONEncoder().encode(value)) as! [String:Any]
            transform(&object)
            return try JSONDecoder().decode(T.self,from:JSONSerialization.data(withJSONObject:object))
        }
        func hydrogen(_ n:Int,_ distance:Double,_ operation:VivoNuclearOperation,_ kind:VivoStationaryKind)->VivoNuclearQualificationRequest {
            let p:[SIMD3<Double>]=n==1 ? [.zero] : n==2 ? [.init(0,0,-distance/2),.init(0,0,distance/2)] : [.init(0,0,-distance),.zero,.init(0,0,distance)]
            let system=VivoElectronicSystem(nuclei:p.enumerated().map{.init(atomicNumber:1,positionBohr:$0.element,structureAtomIndex:UInt32($0.offset))},alphaElectrons:(n+1)/2,betaElectrons:n/2)
            return .init(model:.init(system:system,basis:.hydrogenSTO3G(nucleusIndices:Array(0..<n)),solver:.fullCI),
                massesDa:[Double](repeating:1.008,count:n),operation:operation,kind:kind,
                thermochemistry:.init(rotationalSymmetryNumber:n==1 ? 1:2,electronicDegeneracy:n%2==0 ? 1:2))
        }
        let h=try VivoNuclearQualification.run(hydrogen(1,0,.characterize,.minimum))
        let h2=try VivoNuclearQualification.run(hydrogen(2,1.5,.minimize,.minimum))
        let h3=try VivoNuclearQualification.run(hydrogen(3,1.85,.refineSaddle,.firstOrderSaddle))
        try write(h,"h.json");try write(h2,"h2.json");try write(h3,"h3.json")
        try require(h.thermochemistry.modes.rigidModeCount==3 && h.thermochemistry.modes.signedFrequenciesCM.isEmpty,"atomic RRHO has no invented vibrational or rotational modes")
        try require(h2.optimizationIterations>0 && abs(h2.finalPositionsBohr[1].z-h2.finalPositionsBohr[0].z-1.388694)<1e-4,"native FCI nuclear minimization reaches the H2 bond minimum")
        try require(h2.thermochemistry.modes.rigidModeCount==5 && h2.thermochemistry.modes.signedFrequenciesCM.count==1,"linear molecule retains all 3N-5 vibrational modes")
        let modes=h3.thermochemistry.modes
        try require(h3.optimizationIterations>0 && abs(h3.finalPositionsBohr[2].z-h3.finalPositionsBohr[1].z-1.7702873)<1e-4,"native FCI saddle refinement changes the H3 seed and reaches stationarity")
        try require(modes.rigidModeCount==5 && modes.signedFrequenciesCM.count==4 && modes.signedFrequenciesCM.filter({$0<0}).count==1,"H3 saddle retains one unstable and three stable vibrational modes")
        try require(modes.signedFrequenciesCM.filter({$0>0}).count==3 && h3.thermochemistry.imaginaryFrequencyCM!<0,"unstable mode is excluded from RRHO rather than made positive")
        try require(h3.gradientStepAgreement<h3.request.differences.maximumGradientDifference && h3.hessianStepAgreement<h3.request.differences.maximumHessianDifference,"nuclear gradients and Hessians satisfy explicit step-halving agreement")
        for point in [h,h2,h3] {try VivoNuclearQualification.validate(point,request:point.request)}
        try require(true,"stationary-point outputs survive numerical reconstruction")
        let water=oracle.records.water
        let wm=try VivoNormalModes.analyze(positions:water.positionsBohr,masses:water.massesDa,
            evaluation:.init(energyHartree:water.electronicEnergy,gradientHartreePerBohr:[SIMD3<Double>](repeating:.zero,count:3)),hessian:water.hessian)
        let wt=try VivoRRHOThermochemistry.evaluate(energyHartree:water.electronicEnergy,modes:wm,kind:.minimum,configuration:.init(rotationalSymmetryNumber:2))
        try require(wm.rigidModeCount==6 && wm.signedFrequenciesCM.count==3 && wt.rotor=="nonlinear","nonlinear water normal analysis removes exactly six rigid modes")
        try write(wt,"water.json")
        // Rigid translation must not change frequencies or ideal rotor entropy.
        let translated=try VivoNormalModes.analyze(positions:water.positionsBohr.map{$0+SIMD3<Double>(2,-1,3)},masses:water.massesDa,
            evaluation:.init(energyHartree:water.electronicEnergy,gradientHartreePerBohr:[SIMD3<Double>](repeating:.zero,count:3)),hessian:water.hessian)
        try require(zip(translated.signedFrequenciesCM,wm.signedFrequenciesCM).allSatisfy({abs($0-$1)<1e-7}),"vibrational spectrum is invariant to overall translation")
        let isotope=try VivoNormalModes.analyze(positions:water.positionsBohr,masses:water.massesDa.map{2*$0},
            evaluation:.init(energyHartree:water.electronicEnergy,gradientHartreePerBohr:[SIMD3<Double>](repeating:.zero,count:3)),hessian:water.hessian)
        try require(zip(isotope.signedFrequenciesCM,wm.signedFrequenciesCM).allSatisfy({abs($0*sqrt(2)-$1)<1e-7}),"mass weighting gives the expected inverse-square-root isotope scaling")
        let barrier=try VivoHarmonicBarrier.estimate(saddle:h3,reactants:[h2,h])
        try require(barrier.reactantMolecularity==2 && barrier.electronicBarrierHartree>0 && barrier.activationGibbsHartree>0,"harmonic activation estimate uses two independently characterized reactants")
        try write(barrier,"barrier.json")
        func concentration(_ point:VivoNuclearQualifiedPoint) throws -> VivoNuclearQualifiedPoint {
            var request=point.request;request.thermochemistry.standardState = .concentration(molPerLitre:1,gasReferencePressurePa:101325)
            let t=try VivoRRHOThermochemistry.evaluate(energyHartree:point.thermochemistry.electronicEnergyHartree,modes:point.thermochemistry.modes,kind:point.request.kind,configuration:request.thermochemistry)
            return .init(schema:point.schema,request:request,finalPositionsBohr:point.finalPositionsBohr,thermochemistry:t,
                gradientStepAgreement:point.gradientStepAgreement,hessianStepAgreement:point.hessianStepAgreement,
                optimizationIterations:point.optimizationIterations,electronicEvaluations:point.electronicEvaluations,qualification:point.qualification)
        }
        let ch=try concentration(h),ch2=try concentration(h2),ch3=try concentration(h3)
        let shifted=try VivoHarmonicBarrier.estimate(saddle:ch3,reactants:[ch2,ch])
        let expectedShift=VivoNuclearUnits.kHartree*298.15*log(1000*VivoAtomicUnits.gasConstantJPerMolK*298.15/101325)
        try require(abs(ch.thermochemistry.gibbsEnergyHartree-h.thermochemistry.gibbsEnergyHartree-expectedShift)<1e-12,"one-molar standard conversion is an explicit per-species Gibbs term")
        try require(abs(shifted.activationGibbsHartree-barrier.activationGibbsHartree+expectedShift)<1e-12,"bimolecular activation standard-state term has the correct stoichiometric sign")
        try write(shifted,"barrier-one-molar.json")
        var cfg=h2.request.thermochemistry;cfg.rotationalSymmetryNumber=1
        let asymmetric=try VivoRRHOThermochemistry.evaluate(energyHartree:h2.thermochemistry.electronicEnergyHartree,modes:h2.thermochemistry.modes,kind:.minimum,configuration:cfg)
        try require(abs(asymmetric.gasEntropyHartreePerK-h2.thermochemistry.gasEntropyHartreePerK-VivoNuclearUnits.kHartree*log(2))<1e-14,"rotational symmetry changes entropy by the declared k ln sigma")
        cfg=h.request.thermochemistry;cfg.electronicDegeneracy=1
        let nondegenerate=try VivoRRHOThermochemistry.evaluate(energyHartree:h.thermochemistry.electronicEnergyHartree,modes:h.thermochemistry.modes,kind:.minimum,configuration:cfg)
        try require(abs(h.thermochemistry.gasEntropyHartreePerK-nondegenerate.gasEntropyHartreePerK-VivoNuclearUnits.kHartree*log(2))<1e-14,"electronic degeneracy remains an explicit entropy input")
        try reject("minimum cannot be relabeled as a first-order saddle") {_=try VivoRRHOThermochemistry.evaluate(energyHartree:wt.electronicEnergyHartree,modes:wm,kind:.firstOrderSaddle)}
        try reject("first-order saddle cannot be labeled a minimum") {_=try VivoRRHOThermochemistry.evaluate(energyHartree:h3.thermochemistry.electronicEnergyHartree,modes:modes,kind:.minimum)}
        try reject("nonstationary geometry cannot acquire RRHO qualification") {_=try VivoNuclearQualification.run(hydrogen(2,2.0,.characterize,.minimum))}
        try reject("unresolved low-frequency modes are not silently deleted") {
            let soft=try mutate(wm){$0["signedFrequenciesCM"]=[0.0,100.0,200.0]}
            _=try VivoRRHOThermochemistry.evaluate(energyHartree:water.electronicEnergy,modes:soft,kind:.minimum)
        }
        try reject("extra unstable modes are rejected") {
            let unstable=try mutate(wm){$0["signedFrequenciesCM"]=[-50.0,-100.0,200.0]}
            _=try VivoRRHOThermochemistry.evaluate(energyHartree:water.electronicEnergy,modes:unstable,kind:.firstOrderSaddle)
        }
        try reject("nonpositive isotope masses are rejected") {
            _=try VivoNormalModes.analyze(positions:water.positionsBohr,masses:[0,1,1],evaluation:.init(energyHartree:0,gradientHartreePerBohr:[.zero,.zero,.zero]),hessian:water.hessian)
        }
        try reject("zero thermochemical temperature is rejected") {var c=cfg;c.temperatureK=0;try c.validate()}
        try reject("energy budget exhaustion cannot publish a stationary point") {
            var r=h2.request;r.differences.maximumEnergyEvaluations=1;_=try VivoNuclearQualification.run(r)
        }
        try reject("forged activation Gibbs output is rejected by stationary reconstruction") {
            let forged=try mutate(h3) {obj in var t=obj["thermochemistry"] as! [String:Any];t["gibbsEnergyHartree"]=99.0;obj["thermochemistry"]=t}
            try VivoNuclearQualification.validate(forged,request:h3.request)
        }
        try reject("wrong reaction isotope inventory is rejected") {_=try VivoHarmonicBarrier.estimate(saddle:h3,reactants:[h2,h2])}
        try reject("basis identity changes cannot pass a balanced reaction comparison") {
            var request=h.request
            request.model.basis.identifier = "unrelated-basis"
            let other=try VivoNuclearQualification.run(request)
            _=try VivoHarmonicBarrier.estimate(saddle:h3,reactants:[h2,other])
        }
        let solventRequest=VivoSolvatedECCPathRequest(path:try VivoMolecularECCPath.hydrogenStretchTemplate(),solvent:.init(dielectricConstant:4,angularPoints:50))
        let solvated=try VivoSolvatedECCPath.solve(solventRequest)
        try VivoSolvatedECCPath.validate(solvated,request:solventRequest)
        try require(solvated.path.converged && solvated.referenceSolventFields.count==5,"smooth reference-polarized field and shared ECC solve at all five geometries")
        try write(solvated,"solvated-path.json")
        let s=solventRequest.path.snapshots[1].system,basis=solventRequest.path.basis
        let ao=try VivoGaussianIntegralEngine.compute(system:s,basis:basis)
        let prepared=try VivoReferencePolarization.prepare(system:s,ao:ao,scf:solventRequest.path.reference,solvent:solventRequest.solvent)
        let state=VivoCIState(orbitalCount:2,alphaElectrons:1,betaElectrons:1,determinants:[3],coefficients:[1])
        let determinantEnergy=try VivoCIDensityMatrices.compute(state).energy(of:prepared.hamiltonian)
        let expectedEnergy=prepared.reference.energyHartree
        try require(abs(determinantEnergy-expectedEnergy)<1e-9,"frozen-field scalar restores the reference energy without double counting polarization")
        let dm=try prepared.reference.alphaDensity.adding(prepared.reference.betaDensity)
        let tr=zip(dm.values,prepared.solvent!.reactionPotentialMatrix.values).reduce(0.0){$0+$1.0*$1.1}
        try require(abs(prepared.frozenFieldConstantHartree-prepared.solvent!.polarizationEnergyHartree+tr)<1e-12,"surface self and nuclear constant equals Gpol minus reference electron-field contraction")
        let anchor=try VivoSolvatedECCPath.nuclearModel(solvated,point:1)
        let surface=try VivoNuclearElectronicSurface(model:anchor)
        let positions=s.nuclei.map(\.positionBohr)
        let anchoredEnergy=try surface.energy(positions)
        try require(abs(anchoredEnergy-solvated.path.pointResults[1].energyHartree)<1e-8,"nuclear anchor retains the actual verified shared orbital rotation")
        let anchoredDerivative=try surface.gradient(positions)
        let fciModel=VivoNuclearElectronicModel(system:s,basis:basis,solver:.fullCI,solvent:solventRequest.solvent)
        let fciSurface=try VivoNuclearElectronicSurface(model:fciModel)
        let independentDerivative=try fciSurface.gradient(positions)
        let derivativeError=zip(anchoredDerivative.gradientHartreePerBohr,independentDerivative.gradientHartreePerBohr).map{vivoQMNorm($0-$1)}.max()!
        try require(derivativeError<1e-7,"full-bath anchored ECC nuclear derivatives agree with a separate frozen-field FCI surface")
        try reject("forged frozen-field scalar is rejected") {
            let forged=try mutate(solvated){$0["frozenFieldConstantsHartree"]=[Double](repeating:1,count:5)}
            try VivoSolvatedECCPath.validate(forged,request:solventRequest)
        }
        let descent=try VivoNuclearDescent.trace(h3,configuration:.init(maximumStepsPerDirection:30))
        try require(descent.forward.count>5 && descent.reverse.count>5,"actual unstable-mode displacements produce two integrated downhill branches")
        func downhill(_ points:[VivoNuclearDescentPoint])->Bool {
            points.allSatisfy{$0.energyHartree<h3.thermochemistry.electronicEnergyHartree} && zip(points,points.dropFirst()).allSatisfy{$1.energyHartree<$0.energyHartree}
        }
        try require(downhill(descent.forward) && downhill(descent.reverse),"both descent branches remain downhill without claiming endpoint connectivity")
        func bondDifference(_ p:VivoNuclearDescentPoint)->Double {vivoQMNorm(p.positionsBohr[0]-p.positionsBohr[1])-vivoQMNorm(p.positionsBohr[1]-p.positionsBohr[2])}
        try require(bondDifference(descent.forward.last!)*bondDifference(descent.reverse.last!)<0,"H3 branches shorten opposite mapped H-H bonds")
        try write(descent,"descent.json")
        observations["H2ElectronicEnergyHartree"]=h2.thermochemistry.electronicEnergyHartree
        observations["H3SaddleEnergyHartree"]=h3.thermochemistry.electronicEnergyHartree
        observations["electronicBarrierHartree"]=barrier.electronicBarrierHartree
        observations["activationGibbsHartree"]=barrier.activationGibbsHartree
        observations["standardStateContributionHartree"]=shifted.standardStateContributionHartree
        observations["H3GradientStepAgreement"]=h3.gradientStepAgreement
        observations["H3HessianStepAgreement"]=h3.hessianStepAgreement
        observations["anchoredFullBathGradientError"]=derivativeError
        struct Report:Codable {let checks:[String];let observations:[String:Double];let passed:Bool;let scope:String}
        try write(Report(checks:labels,observations:observations,passed:true,
            scope:"native H2 FCI minimum and H3 FCI local saddle, water analytic-HF Hessian thermochemistry, reference-polarized five-point full-bath H2 ECC; not correlated self-consistent PCM, converged IRC, paper chemistry or kinetic validation"),"checks.json")
        print("PASS \(labels.count) reaction qualification checks")
    }
}
