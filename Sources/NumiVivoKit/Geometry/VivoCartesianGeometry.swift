import Foundation

public struct VivoGeometryEvaluation: Sendable {
    public let energyHartree: Double
    public let gradientHartreePerBohr: [SIMD3<Double>]
    public init(energyHartree: Double, gradientHartreePerBohr: [SIMD3<Double>]) {
        self.energyHartree=energyHartree;self.gradientHartreePerBohr=gradientHartreePerBohr
    }
}
public struct VivoGeometryConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var maximumEvaluations: Int
    public var gradientRMSTolerance: Double
    public var maximumGradientTolerance: Double
    public var maximumStepBohr: Double
    public init(maximumIterations: Int = 100, maximumEvaluations: Int = 1000,
                gradientRMSTolerance: Double = 1e-5, maximumGradientTolerance: Double = 3e-5,
                maximumStepBohr: Double = 0.2) {
        self.maximumIterations=maximumIterations;self.maximumEvaluations=maximumEvaluations
        self.gradientRMSTolerance=gradientRMSTolerance;self.maximumGradientTolerance=maximumGradientTolerance
        self.maximumStepBohr=maximumStepBohr
    }
}
public struct VivoGeometryResult: Codable, Sendable, Equatable {
    public let positionsBohr: [SIMD3<Double>]
    public let energyHartree: Double
    public let gradientRMSTolerance: Double
    public let finalGradientRMS: Double
    public let maximumGradient: Double
    public let iterations: Int
    public let evaluations: Int
    public let converged: Bool
    public let termination: String
}
public enum VivoCartesianGeometry {
    public static func minimize(positionsBohr: [SIMD3<Double>], frozenCoordinates: Set<Int> = [],
                                configuration cfg: VivoGeometryConfiguration = .init(), budget: VivoChemistryBudget = .init(),
                                evaluate: ([SIMD3<Double>]) throws -> VivoGeometryEvaluation) throws -> VivoGeometryResult {
        let count=positionsBohr.count*3
        guard !positionsBohr.isEmpty,positionsBohr.allSatisfy(vivoQMFinite),count<=3*budget.maximumBasisFunctions,
              frozenCoordinates.allSatisfy({$0>=0 && $0<count}),frozenCoordinates.count<count,
              cfg.maximumIterations>0,cfg.maximumEvaluations>0,
              [cfg.gradientRMSTolerance,cfg.maximumGradientTolerance,cfg.maximumStepBohr].allSatisfy({$0.isFinite && $0>0}) else {
            throw VivoChemistryError.invalid("geometry coordinates/constraints/configuration")
        }
        _ = try budget.elements([count,count],simultaneousArrays:6)
        let free=(0..<count).filter{!frozenCoordinates.contains($0)},n=free.count
        var positions=positionsBohr,h=try VivoQMMatrix.identity(n),calls=0,iterations=0
        func sample(_ p: [SIMD3<Double>]) throws -> (Double,[Double]) {
            guard calls<cfg.maximumEvaluations else { throw VivoChemistryError.resourceLimit("geometry evaluation budget") }
            calls += 1;let result=try evaluate(p)
            guard result.energyHartree.isFinite,result.gradientHartreePerBohr.count==positionsBohr.count,
                  result.gradientHartreePerBohr.allSatisfy(vivoQMFinite) else { throw VivoChemistryError.invalid("geometry evaluator result") }
            return (result.energyHartree,free.map{result.gradientHartreePerBohr[$0/3][$0%3]})
        }
        var (energy,gradient)=try sample(positions),converged=false,termination="iteration-limit"
        for iteration in 0..<cfg.maximumIterations {
            iterations=iteration
            let rms=sqrt(gradient.reduce(0){$0+$1*$1}/Double(n)),maximum=gradient.map(abs).max() ?? 0
            if rms<=cfg.gradientRMSTolerance && maximum<=cfg.maximumGradientTolerance {
                converged=true;termination="gradient-tolerances";break
            }
            var step=(0..<n).map { i in -(0..<n).reduce(0.0){$0+h[i,$1]*gradient[$1]} }
            var derivative=zip(gradient,step).reduce(0){$0+$1.0*$1.1}
            if !derivative.isFinite || derivative>=0 { h=try .identity(n);step=gradient.map{-$0} }
            let norm=sqrt(step.reduce(0){$0+$1*$1})
            if norm>cfg.maximumStepBohr { step=step.map{$0*cfg.maximumStepBohr/norm} }
            derivative=zip(gradient,step).reduce(0){$0+$1.0*$1.1}
            var scale=1.0,accepted=false
            for _ in 0..<20 {
                if calls>=cfg.maximumEvaluations { termination="evaluation-budget";break }
                var trial=positions
                for (i,coordinate) in free.enumerated() { trial[coordinate/3][coordinate%3] += scale*step[i] }
                let (nextEnergy,nextGradient)=try sample(trial)
                if nextEnergy<=energy+1e-4*scale*derivative {
                    let s=step.map{$0*scale},y=zip(nextGradient,gradient).map{$0.0-$0.1}
                    let sy=zip(s,y).reduce(0){$0+$1.0*$1.1},yNorm=sqrt(y.reduce(0){$0+$1*$1})
                    if sy>1e-12*max(1e-20,scale*norm*yNorm) {
                        let hy=(0..<n).map{i in (0..<n).reduce(0.0){$0+h[i,$1]*y[$1]}}
                        let yhy=zip(y,hy).reduce(0){$0+$1.0*$1.1}
                        for i in 0..<n { for j in 0..<n { h[i,j] += (1+yhy/sy)*s[i]*s[j]/sy-(hy[i]*s[j]+s[i]*hy[j])/sy } }
                    } else { h=try .identity(n) }
                    positions=trial;energy=nextEnergy;gradient=nextGradient;accepted=true;iterations=iteration+1;break
                }
                scale *= 0.5
            }
            if !accepted { if termination != "evaluation-budget" { termination="line-search-stalled" };break }
        }
        let rms=sqrt(gradient.reduce(0){$0+$1*$1}/Double(n)),maximum=gradient.map(abs).max() ?? 0
        if rms<=cfg.gradientRMSTolerance && maximum<=cfg.maximumGradientTolerance { converged=true;termination="gradient-tolerances" }
        return .init(positionsBohr:positions,energyHartree:energy,gradientRMSTolerance:cfg.gradientRMSTolerance,
                     finalGradientRMS:rms,maximumGradient:maximum,iterations:iterations,evaluations:calls,converged:converged,termination:termination)
    }
}
public enum VivoHartreeFockGeometry {
    public static func energyAndGradient(system: VivoElectronicSystem, basis: VivoGaussianBasis,
                                         stepBohr: Double = 1e-4, configuration: VivoSCFConfiguration = .init(),
                                         budget: VivoChemistryBudget = .init()) throws -> VivoGeometryEvaluation {
        guard stepBohr.isFinite,stepBohr>=1e-6,stepBohr<=1e-2 else { throw VivoChemistryError.invalid("HF finite-difference step") }
        func energy(_ model: VivoElectronicSystem) throws -> Double {
            let ao=try VivoGaussianIntegralEngine.compute(system:model,basis:basis,budget:budget)
            return try VivoHartreeFock.solve(system:model,integrals:ao,configuration:configuration,budget:budget).energyHartree
        }
        let center=try energy(system)
        var gradient=[SIMD3<Double>](repeating:.zero,count:system.nuclei.count)
        for atom in system.nuclei.indices { for axis in 0..<3 {
            var plus=system,minus=system
            plus.nuclei[atom].positionBohr[axis] += stepBohr;minus.nuclei[atom].positionBohr[axis] -= stepBohr
            gradient[atom][axis]=(try energy(plus)-energy(minus))/(2*stepBohr)
        } }
        return .init(energyHartree:center,gradientHartreePerBohr:gradient)
    }
}
