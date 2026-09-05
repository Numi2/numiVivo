import Foundation

public struct VivoCPCMCavityConfiguration: Codable, Sendable, Equatable {
    public var pointsPerAtom: Int
    public var radiusScale: Double
    public var radiiAngstromByAtomicNumber: [Int:Double]
    public var maximumTesserae: Int
    public var burialToleranceBohr: Double
    public init(pointsPerAtom: Int = 194, radiusScale: Double = 1.2,
                radiiAngstromByAtomicNumber: [Int:Double] = [1:1.20,6:1.70,7:1.55,8:1.52,9:1.47,15:1.80,16:1.80,17:1.75,35:1.85,53:1.98],
                maximumTesserae: Int = 8192, burialToleranceBohr: Double = 1e-6) {
        self.pointsPerAtom=pointsPerAtom;self.radiusScale=radiusScale;self.radiiAngstromByAtomicNumber=radiiAngstromByAtomicNumber
        self.maximumTesserae=maximumTesserae;self.burialToleranceBohr=burialToleranceBohr
    }
    public func validate(system: VivoElectronicSystem) throws {
        guard pointsPerAtom>=26,pointsPerAtom<=4096,radiusScale.isFinite,radiusScale>0,maximumTesserae>0,
              burialToleranceBohr.isFinite,burialToleranceBohr>=0,
              radiiAngstromByAtomicNumber.allSatisfy({$0.key>0&&$0.value.isFinite&&$0.value>0}),
              system.nuclei.allSatisfy({radiiAngstromByAtomicNumber[$0.atomicNumber] != nil}) else {
            throw VivoChemistryError.invalid("C-PCM cavity radii/tessellation")
        }
        let (count,overflow)=pointsPerAtom.multipliedReportingOverflow(by:system.nuclei.count)
        guard !overflow,count<=maximumTesserae else{throw VivoChemistryError.resourceLimit("C-PCM tessera budget")}
    }
}
public struct VivoCPCMConfiguration: Codable, Sendable, Equatable {
    public var dielectricConstant: Double
    public var cavity: VivoCPCMCavityConfiguration
    public var solverTolerance: Double
    public var maximumSolverIterations: Int
    public var maximumGaussLawRelativeError: Double
    public init(dielectricConstant: Double = 78.3, cavity: VivoCPCMCavityConfiguration = .init(),
                solverTolerance: Double = 1e-10, maximumSolverIterations: Int = 2048,
                maximumGaussLawRelativeError: Double = 0.2) {
        self.dielectricConstant=dielectricConstant;self.cavity=cavity;self.solverTolerance=solverTolerance
        self.maximumSolverIterations=maximumSolverIterations;self.maximumGaussLawRelativeError=maximumGaussLawRelativeError
    }
    public func validate(system: VivoElectronicSystem) throws {
        try cavity.validate(system:system)
        guard dielectricConstant.isFinite,dielectricConstant>1,solverTolerance.isFinite,solverTolerance>0,
              maximumSolverIterations>0,maximumGaussLawRelativeError.isFinite,maximumGaussLawRelativeError>=0 else {
            throw VivoChemistryError.invalid("C-PCM dielectric/solver configuration")
        }
        guard system.pointCharges.isEmpty else {
            throw VivoChemistryError.unsupported("C-PCM with QM/MM embedding charges requires an explicit combined cavity/outlying-charge policy")
        }
    }
}
public struct VivoCPCMTessera: Codable, Sendable, Equatable {
    public let atomIndex:Int
    public let positionBohr:SIMD3<Double>
    public let areaBohr2:Double
    public let selfInteractionPerBohr:Double
}
public struct VivoCPCMResult: Codable, Sendable, Equatable {
    public let dielectricConstant:Double
    public let dielectricScale:Double
    public let tesserae:[VivoCPCMTessera]
    public let solutePotentialHartreePerE:[Double]
    public let apparentSurfaceChargesE:[Double]
    public let polarizationEnergyHartree:Double
    public let reactionPotentialMatrix:VivoQMMatrix
    public let linearResidualNorm:Double
    public let linearRelativeResidual:Double
    public let solverIterations:Int
    public let totalSurfaceChargeE:Double
    public let expectedGaussLawSurfaceChargeE:Double
    public let gaussLawRelativeError:Double
    public let status:String
}
public enum VivoCPCM {
    private static let angstromPerBohr=0.529177210903
    private static func sphere(_ n:Int)->[SIMD3<Double>] {
        let golden=Double.pi*(3-sqrt(5.0))
        return (0..<n).map{i in
            let z=1-2*(Double(i)+0.5)/Double(n),r=sqrt(max(0,1-z*z)),phi=golden*Double(i)
            return .init(r*cos(phi),r*sin(phi),z)
        }
    }
    public static func cavity(system:VivoElectronicSystem,
                              configuration cfg:VivoCPCMCavityConfiguration = .init())throws->[VivoCPCMTessera] {
        try system.validate();try cfg.validate(system:system);let directions=sphere(cfg.pointsPerAtom)
        let radii=try system.nuclei.map{nucleus->Double in
            guard let a=cfg.radiiAngstromByAtomicNumber[nucleus.atomicNumber] else{throw VivoChemistryError.invalid("missing C-PCM radius")}
            return cfg.radiusScale*a/angstromPerBohr
        }
        var result:[VivoCPCMTessera]=[]
        for (atom,nucleus) in system.nuclei.enumerated() {
            let radius=radii[atom],area=4*Double.pi*radius*radius/Double(cfg.pointsPerAtom)
            let selfTerm=1.07*sqrt(4*Double.pi/area)
            for direction in directions {
                let point=nucleus.positionBohr+radius*direction;var buried=false
                for other in system.nuclei.indices where other != atom {
                    if vivoQMNorm(point-system.nuclei[other].positionBohr)<radii[other]-cfg.burialToleranceBohr { buried=true;break }
                }
                if !buried { result.append(.init(atomIndex:atom,positionBohr:point,areaBohr2:area,selfInteractionPerBohr:selfTerm)) }
            }
        }
        guard !result.isEmpty,result.count<=cfg.maximumTesserae else{throw VivoChemistryError.resourceLimit("empty/oversized C-PCM exposed cavity")}
        return result
    }
    private static func applySurfaceMatrix(_ x:[Double],_ t:[VivoCPCMTessera])throws->[Double] {
        guard x.count==t.count,x.allSatisfy(\.isFinite) else{throw VivoChemistryError.invalid("C-PCM surface vector")}
        var y=[Double](repeating:0,count:x.count)
        for i in t.indices {
            y[i]+=t[i].selfInteractionPerBohr*x[i]
            for j in 0..<i {
                let r=vivoQMNorm(t[i].positionBohr-t[j].positionBohr)
                guard r>1e-10 else{throw VivoChemistryError.invalid("coincident C-PCM tesserae")}
                let value=1/r;y[i]+=value*x[j];y[j]+=value*x[i]
            }
        }
        return y
    }
    private static func dot(_ a:[Double],_ b:[Double])->Double{zip(a,b).reduce(0){$0+$1.0*$1.1}}
    private static func solve(rhs:[Double],tesserae:[VivoCPCMTessera],tolerance:Double,
                              maximumIterations:Int)throws->([Double],Double,Int) {
        let n=rhs.count;var x=[Double](repeating:0,count:n),r=rhs
        var z=zip(r,tesserae).map{$0.0/$0.1.selfInteractionPerBohr},p=z,rz=dot(r,z)
        let rhsNorm=sqrt(dot(rhs,rhs));if rhsNorm==0{return(x,0,0)}
        guard rz.isFinite,rz>0 else{throw VivoChemistryError.convergence("C-PCM initial PCG residual")}
        for iteration in 1...maximumIterations {
            let ap=try applySurfaceMatrix(p,tesserae);let denominator=dot(p,ap)
            guard denominator.isFinite,denominator>0 else{
                throw VivoChemistryError.convergence("C-PCM surface matrix is not positive definite under current tessellation")
            }
            let alpha=rz/denominator
            for i in 0..<n{x[i]+=alpha*p[i];r[i]-=alpha*ap[i]}
            let residual=sqrt(dot(r,r));if residual<=tolerance*rhsNorm{return(x,residual,iteration)}
            z=zip(r,tesserae).map{$0.0/$0.1.selfInteractionPerBohr};let next=dot(r,z)
            guard next.isFinite,next>0 else{throw VivoChemistryError.convergence("C-PCM PCG breakdown")}
            let beta=next/rz;for i in 0..<n{p[i]=z[i]+beta*p[i]};rz=next
        }
        throw VivoChemistryError.convergence("C-PCM surface solve iteration limit")
    }
    public static func evaluate(system:VivoElectronicSystem,integrals ao:VivoAOIntegrals,totalDensity:VivoQMMatrix,
                                configuration cfg:VivoCPCMConfiguration = .init(),
                                budget:VivoChemistryBudget = .init())throws->VivoCPCMResult {
        try cfg.validate(system:system);try ao.validate(budget:budget)
        guard system==ao.sourceSystem,totalDensity.rows==ao.count,totalDensity.columns==ao.count,
              totalDensity.values.allSatisfy(\.isFinite) else{throw VivoChemistryError.invalid("C-PCM AO/density binding")}
        let tesserae=try cavity(system:system,configuration:cfg.cavity)
        let (work,overflow)=tesserae.count.multipliedReportingOverflow(by:ao.count*ao.count)
        guard !overflow,work<=budget.maximumOperatorApplications else{throw VivoChemistryError.resourceLimit("C-PCM AO potential work budget")}
        var potential=[Double](repeating:0,count:tesserae.count),matrices:[VivoQMMatrix]=[];matrices.reserveCapacity(tesserae.count)
        for i in tesserae.indices {
            let matrix=try VivoGaussianElectrostaticPotential.matrix(integrals:ao,at:tesserae[i].positionBohr,budget:budget)
            matrices.append(matrix);var value=0.0
            for nucleus in system.nuclei{value += Double(nucleus.atomicNumber)/vivoQMNorm(tesserae[i].positionBohr-nucleus.positionBohr)}
            value -= zip(totalDensity.values,matrix.values).reduce(0){$0+$1.0*$1.1};potential[i]=value
        }
        let scale=(cfg.dielectricConstant-1)/cfg.dielectricConstant,rhs=potential.map{-scale*$0}
        let solution=try solve(rhs:rhs,tesserae:tesserae,tolerance:cfg.solverTolerance,
                               maximumIterations:cfg.maximumSolverIterations),charges=solution.0
        let residual=solution.1,rhsNorm=sqrt(dot(rhs,rhs)),relative=rhsNorm>0 ? residual/rhsNorm:0
        let energy=0.5*dot(charges,potential);var reaction=VivoQMMatrix(ao.count,ao.count)
        for i in charges.indices where charges[i] != 0 { reaction=try reaction.adding(matrices[i],scale:-charges[i]) }
        let soluteCharge=Double(system.nuclei.reduce(0){$0+$1.atomicNumber}-(system.alphaElectrons+system.betaElectrons))
        let expected = -scale * soluteCharge,total=charges.reduce(0,+)
        let gauss=abs(total-expected)/max(1,abs(expected))
        guard gauss<=cfg.maximumGaussLawRelativeError else{
            throw VivoChemistryError.convergence("C-PCM surface charge violates configured Gauss-law diagnostic: \(gauss)")
        }
        guard energy.isFinite,reaction.values.allSatisfy(\.isFinite) else{throw VivoChemistryError.convergence("nonfinite C-PCM reaction field")}
        return .init(dielectricConstant:cfg.dielectricConstant,dielectricScale:scale,tesserae:tesserae,
                     solutePotentialHartreePerE:potential,apparentSurfaceChargesE:charges,
                     polarizationEnergyHartree:energy,reactionPotentialMatrix:reaction,
                     linearResidualNorm:residual,linearRelativeResidual:relative,solverIterations:solution.2,
                     totalSurfaceChargeE:total,expectedGaussLawSurfaceChargeE:expected,gaussLawRelativeError:gauss,
                     status:"C-PCM point-charge tessellation; Aii=1.07*sqrt(4pi/area); x=0 dielectric scaling; serialized cavity radii/scale; no cavitation/dispersion terms")
    }
}
