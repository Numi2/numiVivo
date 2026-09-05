import Foundation

public struct VivoNuclearDescentConfiguration: Codable,Sendable,Equatable {
    public var initialDisplacementMassWeighted:Double
    public var stepMassWeighted:Double
    public var maximumStepsPerDirection:Int
    public var endpointMaximumGradient:Double
    public init(initialDisplacementMassWeighted:Double=0.03,stepMassWeighted:Double=0.08,
                maximumStepsPerDirection:Int=200,endpointMaximumGradient:Double=1e-4) {
        self.initialDisplacementMassWeighted=initialDisplacementMassWeighted;self.stepMassWeighted=stepMassWeighted
        self.maximumStepsPerDirection=maximumStepsPerDirection;self.endpointMaximumGradient=endpointMaximumGradient
    }
}
public struct VivoNuclearDescentPoint:Codable,Sendable,Equatable {
    public let positionsBohr:[SIMD3<Double>]
    public let energyHartree:Double
    public let maximumGradient:Double
    public let arcMassWeighted:Double
}
public struct VivoNuclearDescentResult:Codable,Sendable,Equatable {
    public let reverse:[VivoNuclearDescentPoint]
    public let forward:[VivoNuclearDescentPoint]
    public let reverseStationary:Bool
    public let forwardStationary:Bool
    public let energyEvaluations:Int
    public let interpretation:String
}
public enum VivoNuclearDescent {
    /// Mass-weighted steepest-descent branches with midpoint directions and
    /// downhill acceptance. Endpoints remain unassigned chemical states; low
    /// gradient at separated fragments is not a single-molecule minimum.
    public static func trace(_ saddle:VivoNuclearQualifiedPoint,configuration cfg:VivoNuclearDescentConfiguration = .init()) throws -> VivoNuclearDescentResult {
        try VivoNuclearQualification.validate(saddle,request:saddle.request)
        guard saddle.request.kind == .firstOrderSaddle,
              [cfg.initialDisplacementMassWeighted,cfg.stepMassWeighted,cfg.endpointMaximumGradient].allSatisfy({$0.isFinite && $0>0}),
              cfg.stepMassWeighted<=1,(1...10000).contains(cfg.maximumStepsPerDirection) else {throw VivoChemistryError.invalid("nuclear descent contract")}
        let mass=saddle.request.massesDa,n=3*mass.count
        let surface=try VivoNuclearElectronicSurface(model:saddle.request.model,differences:saddle.request.differences)
        // Reconstruct eigenvectors; a decoded vector/sign is not numerical authority.
        let center=try surface.gradient(saddle.finalPositionsBohr)
        let modes=try VivoNormalModes.analyze(positions:saddle.finalPositionsBohr,masses:mass,evaluation:center,
            hessian:saddle.thermochemistry.modes.cartesianHessian,budget:saddle.request.model.budget)
        let negative=modes.signedFrequenciesCM.firstIndex(where:{$0<0})!
        func direction(_ value:VivoGeometryEvaluation) throws -> [Double] {
            var g=value.gradientHartreePerBohr.flatMap{[$0.x,$0.y,$0.z]}
            for i in 0..<n {g[i]/=sqrt(mass[i/3])}
            let length=g.reduce(0.0){hypot($0,$1)}
            guard length.isFinite,length>1e-16 else {throw VivoChemistryError.convergence("unresolved descent tangent")}
            return g.map{-$0/length}
        }
        func shifted(_ positions:[SIMD3<Double>],_ vector:[Double],_ scale:Double)->[SIMD3<Double>] {
            var x=positions
            for i in 0..<n {x[i/3][i%3]+=scale*vector[i]/sqrt(mass[i/3])};return x
        }
        func branch(_ sign:Double) throws -> ([VivoNuclearDescentPoint],Bool) {
            let mode=(0..<n).map{modes.massWeightedModes[$0,negative]}
            var positions=shifted(saddle.finalPositionsBohr,mode,sign*cfg.initialDisplacementMassWeighted)
            var value=try surface.gradient(positions),arc=cfg.initialDisplacementMassWeighted,history:[VivoNuclearDescentPoint]=[]
            guard value.energyHartree<saddle.thermochemistry.electronicEnergyHartree else {throw VivoChemistryError.convergence("initial unstable-mode displacement is not downhill")}
            for _ in 0..<cfg.maximumStepsPerDirection {
                let maximum=value.gradientHartreePerBohr.flatMap{[$0.x,$0.y,$0.z]}.map(abs).max()!
                history.append(.init(positionsBohr:positions,energyHartree:value.energyHartree,maximumGradient:maximum,arcMassWeighted:sign*arc))
                if maximum<=cfg.endpointMaximumGradient {return (history,true)}
                let tangent=try direction(value)
                var step=cfg.stepMassWeighted,accepted=false
                for _ in 0..<20 {
                    let midpoint=try surface.gradient(shifted(positions,tangent,step/2)),middle=try direction(midpoint)
                    let next=shifted(positions,middle,step),trial=try surface.gradient(next)
                    if trial.energyHartree<value.energyHartree {
                        positions=next;value=trial;arc+=step;accepted=true;break
                    }
                    step*=0.5
                }
                if !accepted {return (history,false)}
            }
            return (history,false)
        }
        let a=try branch(-1),b=try branch(1)
        return .init(reverse:a.0,forward:b.0,reverseStationary:a.1,forwardStationary:b.1,energyEvaluations:surface.energyEvaluations,
            interpretation:"mass-weighted downhill branch evidence; endpoints require chemical assignment and separate characterization; no step-size-converged IRC or rate certification")
    }
}
