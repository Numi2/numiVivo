import Foundation

public enum VivoNuclearOperation: String, Codable, Sendable {case characterize, minimize, refineSaddle}
public struct VivoNuclearQualificationRequest: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/nuclear-qualification/v1"
    public var schema: String
    public var model: VivoNuclearElectronicModel
    public var massesDa: [Double]
    public var operation: VivoNuclearOperation
    public var kind: VivoStationaryKind
    public var differences: VivoNuclearDifferenceConfiguration
    public var thermochemistry: VivoThermochemistryConfiguration
    public var maximumOptimizationIterations: Int
    public var maximumStepBohr: Double
    public var maximumRigidResidual: Double
    public init(model: VivoNuclearElectronicModel,massesDa:[Double],operation:VivoNuclearOperation,kind:VivoStationaryKind,
                differences:VivoNuclearDifferenceConfiguration = .init(),thermochemistry:VivoThermochemistryConfiguration = .init(),
                maximumOptimizationIterations:Int=100,maximumStepBohr:Double=0.15,maximumRigidResidual:Double=1e-4) {
        schema=Self.schema;self.model=model;self.massesDa=massesDa;self.operation=operation;self.kind=kind
        self.differences=differences;self.thermochemistry=thermochemistry;self.maximumOptimizationIterations=maximumOptimizationIterations
        self.maximumStepBohr=maximumStepBohr;self.maximumRigidResidual=maximumRigidResidual
    }
    public func validate() throws {
        try model.validate();try differences.validate();try thermochemistry.validate()
        guard schema==Self.schema,massesDa.count==model.system.nuclei.count,massesDa.allSatisfy({$0.isFinite && $0>0}),
              (1...1000).contains(maximumOptimizationIterations),maximumStepBohr.isFinite,maximumStepBohr>0,maximumStepBohr<=1,
              maximumRigidResidual.isFinite,maximumRigidResidual>0,
              operation != .minimize || kind == .minimum, operation != .refineSaddle || kind == .firstOrderSaddle else {
            throw VivoChemistryError.invalid("nuclear qualification request, masses or optimization kind")
        }
        guard model.system.pointCharges.isEmpty,model.solvent?.extraCavitySpheres.isEmpty ?? true else {
            throw VivoChemistryError.unsupported("free-molecule RRHO cannot remove rigid modes in a frozen external field; constrained/environment thermochemistry requires another model")
        }
    }
}
public struct VivoNuclearQualifiedPoint: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/nuclear-qualified-point/v1"
    public let schema:String
    public let request:VivoNuclearQualificationRequest
    public let finalPositionsBohr:[SIMD3<Double>]
    public let thermochemistry:VivoThermochemistryResult
    public let gradientStepAgreement:Double
    public let hessianStepAgreement:Double
    public let optimizationIterations:Int
    public let electronicEvaluations:Int
    public let qualification:String
}
public enum VivoNuclearQualification {
    private static func flat(_ g:[SIMD3<Double>])->[Double] {g.flatMap{[$0.x,$0.y,$0.z]}}
    private static func norm(_ g:[SIMD3<Double>])->Double {flat(g).reduce(0){hypot($0,$1)}}
    public static func run(_ request:VivoNuclearQualificationRequest) throws -> VivoNuclearQualifiedPoint {
        try request.validate()
        let surface=try VivoNuclearElectronicSurface(model:request.model,differences:request.differences)
        var positions=request.model.system.nuclei.map(\.positionBohr),iterations=0
        let settings=request.thermochemistry
        if request.operation == .minimize {
            let result=try VivoCartesianGeometry.minimize(positionsBohr:positions,
                configuration:.init(maximumIterations:request.maximumOptimizationIterations,maximumEvaluations:request.differences.maximumEnergyEvaluations,
                    gradientRMSTolerance:settings.gradientRMSTolerance,maximumGradientTolerance:settings.maximumGradientTolerance,maximumStepBohr:request.maximumStepBohr),
                budget:request.model.budget,evaluate:{try surface.gradient($0)})
            guard result.converged else {throw VivoChemistryError.convergence("nuclear minimization: \(result.termination)")}
            positions=result.positionsBohr;iterations=result.iterations
        } else if request.operation == .refineSaddle {
            var converged=false
            for iteration in 0...request.maximumOptimizationIterations {
                let value=try surface.gradient(positions),h=try surface.hessian(positions)
                let modes=try VivoNormalModes.analyze(positions:positions,masses:request.massesDa,evaluation:value,hessian:h,budget:request.model.budget)
                guard modes.signedFrequenciesCM.filter({$0 < -settings.minimumResolvedFrequencyCM}).count==1 else {
                    throw VivoChemistryError.convergence("local saddle refiner requires a seed with one resolved unstable vibrational mode")
                }
                iterations=iteration
                if modes.gradientRMS<=settings.gradientRMSTolerance && modes.maximumGradient<=settings.maximumGradientTolerance {converged=true;break}
                if iteration==request.maximumOptimizationIterations {break}
                let n=3*positions.count,mass=request.massesDa,g=flat(value.gradientHartreePerBohr)
                let conversion=sqrt(VivoNuclearUnits.hartreeJ/(VivoNuclearUnits.atomicMassKG*pow(VivoNuclearUnits.bohrM,2)))/(2*Double.pi*VivoNuclearUnits.lightCMPerS)
                var step=[Double](repeating:0,count:n)
                for j in modes.signedFrequenciesCM.indices {
                    let frequency=modes.signedFrequenciesCM[j]
                    guard abs(frequency)>=settings.minimumResolvedFrequencyCM else {throw VivoChemistryError.convergence("unresolved saddle curvature")}
                    let lambda=(frequency<0 ? -1.0:1.0)*pow(frequency/conversion,2)
                    var component=0.0
                    for i in 0..<n {component+=modes.massWeightedModes[i,j]*g[i]/sqrt(mass[i/3])}
                    for i in 0..<n {step[i]-=component/lambda*modes.massWeightedModes[i,j]/sqrt(mass[i/3])}
                }
                let length=step.reduce(0){hypot($0,$1)},bound=min(1,request.maximumStepBohr/max(1e-30,length))
                var accepted=false,scale=bound
                for _ in 0..<20 {
                    var candidate=positions
                    for i in 0..<n {candidate[i/3][i%3]+=scale*step[i]}
                    let trial=try surface.gradient(candidate)
                    if norm(trial.gradientHartreePerBohr)<norm(value.gradientHartreePerBohr)*(1-1e-4*scale) {
                        positions=candidate;accepted=true;break
                    }
                    scale*=0.5
                }
                guard accepted else {throw VivoChemistryError.convergence("saddle residual line search stalled")}
            }
            guard converged else {throw VivoChemistryError.convergence("saddle optimization iteration limit")}
        }
        return try characterize(request,positions:positions,surface:surface,iterations:iterations)
    }
    private static func characterize(_ request:VivoNuclearQualificationRequest,positions:[SIMD3<Double>],surface:VivoNuclearElectronicSurface,
                                     iterations:Int) throws -> VivoNuclearQualifiedPoint {
        let cfg=request.differences,coarse=try surface.gradient(positions),fine=try surface.gradient(positions,step:cfg.gradientStepBohr/2)
        let coarseH=try surface.hessian(positions),fineH=try surface.hessian(positions,step:cfg.hessianStepBohr/2)
        let gradientError=zip(flat(coarse.gradientHartreePerBohr),flat(fine.gradientHartreePerBohr)).map{abs($0-$1)}.max()!
        let hessianError=zip(coarseH.values,fineH.values).map{abs($0-$1)}.max()!
        guard gradientError<=cfg.maximumGradientDifference,hessianError<=cfg.maximumHessianDifference else {
            throw VivoChemistryError.convergence("nuclear derivative step-halving disagreement: gradient \(gradientError), Hessian \(hessianError)")
        }
        let corrected=zip(coarse.gradientHartreePerBohr,fine.gradientHartreePerBohr).map{(4*$1-$0)/3}
        let h=try fineH.scaled(4.0/3).adding(coarseH,scale:-1.0/3)
        let evaluation=VivoGeometryEvaluation(energyHartree:fine.energyHartree,gradientHartreePerBohr:corrected)
        let modes=try VivoNormalModes.analyze(positions:positions,masses:request.massesDa,evaluation:evaluation,hessian:h,budget:request.model.budget)
        guard modes.rigidHessianResidual<=request.maximumRigidResidual else {throw VivoChemistryError.convergence("Hessian does not have the declared free-molecule rigid invariance")}
        let thermo=try VivoRRHOThermochemistry.evaluate(energyHartree:evaluation.energyHartree,modes:modes,kind:request.kind,configuration:request.thermochemistry)
        return .init(schema:VivoNuclearQualifiedPoint.schema,request:request,finalPositionsBohr:positions,thermochemistry:thermo,
            gradientStepAgreement:gradientError,hessianStepAgreement:hessianError,optimizationIterations:iterations,electronicEvaluations:surface.energyEvaluations,
            qualification:"local stationary point and resolved Hessian index; ideal RRHO thermochemistry; reaction connectivity and dynamical rate not certified")
    }
    public static func validate(_ result:VivoNuclearQualifiedPoint,request:VivoNuclearQualificationRequest) throws {
        try request.validate()
        guard result.schema==VivoNuclearQualifiedPoint.schema,result.request==request,
              result.finalPositionsBohr.count==request.model.system.nuclei.count,result.finalPositionsBohr.allSatisfy(vivoQMFinite),
              result.electronicEvaluations>0,result.electronicEvaluations<=request.differences.maximumEnergyEvaluations,
              result.optimizationIterations>=0,result.optimizationIterations<=request.maximumOptimizationIterations else {
            throw VivoChemistryError.invalid("qualified nuclear result request or dimensions")
        }
        let surface=try VivoNuclearElectronicSurface(model:request.model,differences:request.differences)
        let rebuilt=try characterize(request,positions:result.finalPositionsBohr,surface:surface,iterations:result.optimizationIterations)
        guard abs(rebuilt.thermochemistry.gibbsEnergyHartree-result.thermochemistry.gibbsEnergyHartree)<1e-8,
              abs(rebuilt.thermochemistry.electronicEnergyHartree-result.thermochemistry.electronicEnergyHartree)<1e-9,
              abs(rebuilt.thermochemistry.zeroPointEnergyHartree-result.thermochemistry.zeroPointEnergyHartree)<1e-9,
              abs(rebuilt.thermochemistry.enthalpyCorrectionHartree-result.thermochemistry.enthalpyCorrectionHartree)<1e-9,
              abs(rebuilt.thermochemistry.gasEntropyHartreePerK-result.thermochemistry.gasEntropyHartreePerK)<1e-11,
              abs(rebuilt.thermochemistry.standardStateShiftHartree-result.thermochemistry.standardStateShiftHartree)<1e-12,
              result.thermochemistry.modes.positionsBohr==result.finalPositionsBohr,result.thermochemistry.modes.massesDa==request.massesDa,
              try rebuilt.thermochemistry.modes.cartesianHessian.adding(result.thermochemistry.modes.cartesianHessian,scale:-1).frobeniusNorm<1e-6,
              abs(rebuilt.gradientStepAgreement-result.gradientStepAgreement)<1e-10,
              abs(rebuilt.hessianStepAgreement-result.hessianStepAgreement)<1e-8,
              rebuilt.thermochemistry.kind==result.thermochemistry.kind,
              rebuilt.thermochemistry.configuration==result.thermochemistry.configuration,
              rebuilt.thermochemistry.modes.signedFrequenciesCM.count==result.thermochemistry.modes.signedFrequenciesCM.count,
              zip(rebuilt.thermochemistry.modes.signedFrequenciesCM,result.thermochemistry.modes.signedFrequenciesCM).allSatisfy({abs($0-$1)<1e-3}),
              result.qualification==rebuilt.qualification else {throw VivoChemistryError.invalid("nuclear result energy, mode or qualification claim differs on reconstruction")}
    }
}
public struct VivoHarmonicBarrierEstimate: Codable,Sendable,Equatable {
    public let electronicBarrierHartree:Double
    public let zeroPointBarrierHartree:Double
    public let activationGibbsHartree:Double
    public let standardStateContributionHartree:Double
    public let reactantMolecularity:Int
    public let model:String
}
public enum VivoHarmonicBarrier {
    /// Each reactant entry is one molecule, not a precomplex silently standing in
    /// for separated reactants. No rate is published: connectivity is separate.
    public static func estimate(saddle:VivoNuclearQualifiedPoint,reactants:[VivoNuclearQualifiedPoint]) throws -> VivoHarmonicBarrierEstimate {
        guard (1...16).contains(reactants.count),saddle.thermochemistry.kind == .firstOrderSaddle,
              reactants.allSatisfy({$0.thermochemistry.kind == .minimum}) else {throw VivoChemistryError.invalid("barrier requires one saddle and reactant minima")}
        func basisByElement(_ model:VivoNuclearElectronicModel) throws -> [Int:[VivoGaussianShell]] {
            var result:[Int:[VivoGaussianShell]]=[:]
            for (i,atom) in model.system.nuclei.enumerated() {
                let shells=model.basis.shells.filter{$0.nucleusIndex==i}.map { VivoGaussianShell(nucleusIndex:0,angularMomentum:$0.angularMomentum,primitives:$0.primitives) }
                guard !shells.isEmpty else {throw VivoChemistryError.invalid("missing atom basis in barrier")}
                if let previous=result[atom.atomicNumber],previous != shells {throw VivoChemistryError.unsupported("atom-specific unequal basis recipes require mapped fragment comparison")}
                result[atom.atomicNumber]=shells
            };return result
        }
        let saddleBasis=try basisByElement(saddle.request.model)
        let all=[saddle]+reactants
        for point in all {try VivoNuclearQualification.validate(point,request:point.request)}
        let model=saddle.request.model,cfg=saddle.thermochemistry.configuration
        func composition(_ point:VivoNuclearQualifiedPoint)->[String:Int] {
            var result:[String:Int]=[:]
            for (a,m) in zip(point.request.model.system.nuclei,point.request.massesDa) {result["\(a.atomicNumber):\(m.bitPattern)",default:0]+=1}
            return result
        }
        var sum:[String:Int]=[:]
        for point in reactants {
            let other=point.request.model,tc=point.thermochemistry.configuration
            guard other.solver==model.solver,other.basis.identifier==model.basis.identifier,other.basis.representation==model.basis.representation,other.scf==model.scf,other.solvent==model.solvent,
                  other.eccFrame==nil,model.eccFrame==nil,
                  tc.temperatureK==cfg.temperatureK,tc.standardState==cfg.standardState else {
                throw VivoChemistryError.invalid("barrier electronic model/basis/solvent/standard-state mismatch or unqualified independent ECC frames")
            }
            for (element,shells) in try basisByElement(other) {
                guard saddleBasis[element]==shells else {throw VivoChemistryError.invalid("barrier basis primitives differ despite matching basis names")}
            }
            for (key,count) in composition(point) {sum[key,default:0]+=count}
        }
        guard sum==composition(saddle),reactants.reduce(0,{$0+$1.request.model.system.alphaElectrons})==model.system.alphaElectrons,reactants.reduce(0,{$0+$1.request.model.system.betaElectrons})==model.system.betaElectrons else {
            throw VivoChemistryError.invalid("barrier changes element/isotope inventory or total electrons")
        }
        let t=saddle.thermochemistry
        let electronic=t.electronicEnergyHartree-reactants.reduce(0){$0+$1.thermochemistry.electronicEnergyHartree}
        return .init(electronicBarrierHartree:electronic,
            zeroPointBarrierHartree:electronic+t.zeroPointEnergyHartree-reactants.reduce(0){$0+$1.thermochemistry.zeroPointEnergyHartree},
            activationGibbsHartree:t.gibbsEnergyHartree-reactants.reduce(0){$0+$1.thermochemistry.gibbsEnergyHartree},
            standardStateContributionHartree:t.standardStateShiftHartree-reactants.reduce(0){$0+$1.thermochemistry.standardStateShiftHartree},
            reactantMolecularity:reactants.count,model:"local-saddle harmonic TST free-energy estimate; independent reactant molecules; no endpoint-connectivity, tunneling or dynamical-transmission qualification")
    }
}
