import Foundation

public struct VivoSolventSphere:Codable,Sendable,Equatable {
    public let centerBohr:SIMD3<Double>
    public let radiusBohr:Double
    public init(centerBohr:SIMD3<Double>,radiusBohr:Double) { self.centerBohr=centerBohr;self.radiusBohr=radiusBohr }
}
public struct VivoSmoothCPCMConfiguration:Codable,Sendable,Equatable {
    public var dielectricConstant:Double
    public var angularPoints:Int
    public var radiusScale:Double
    public var radiiAngstrom:[Int:Double]
    public var extraCavitySpheres:[VivoSolventSphere]
    public var maximumTesserae:Int
    public var tolerance:Double
    public var maximumIterations:Int
    public init(dielectricConstant:Double=78.3,angularPoints:Int=302,radiusScale:Double=1.2,
                radiiAngstrom:[Int:Double]=[1:1.2,6:1.7,7:1.55,8:1.52,9:1.47,15:1.8,16:1.8,17:1.75],
                extraCavitySpheres:[VivoSolventSphere]=[],maximumTesserae:Int=8192,tolerance:Double=1e-10,maximumIterations:Int=2000) {
        self.dielectricConstant=dielectricConstant;self.angularPoints=angularPoints;self.radiusScale=radiusScale
        self.radiiAngstrom=radiiAngstrom;self.extraCavitySpheres=extraCavitySpheres;self.maximumTesserae=maximumTesserae
        self.tolerance=tolerance;self.maximumIterations=maximumIterations
    }
}
public struct VivoSmoothCPCMResult:Codable,Sendable,Equatable {
    public let configuration:VivoSmoothCPCMConfiguration
    public let polarizationEnergyHartree:Double
    public let reactionPotentialMatrix:VivoQMMatrix
    public let apparentSurfaceCharges:[Double]
    public let trueRelativeResidual:Double
    public let surfaceChargeDefect:Double
    public let iterations:Int
    public let tesseraCount:Int
    public let method:String
}
/// Lebedev-Laikov cubic-symmetry orbits. Numerical quadrature coefficients are
/// tabulated independently of any provider runtime; checked against PySCF 2.8.
/// Supported 50/110/302-point rules and Gaussian exponent parameters follow
/// the smooth discretization used in the public PySCF PCM reference.
enum VivoSolventLebedev {
    static func grid(_ count:Int) throws -> [(SIMD3<Double>,Double)] {
        var out:[(SIMD3<Double>,Double)]=[]
        func orbit(_ components:[Double],_ weight:Double) {
            let permutations=[[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]]
            var points=Set<SIMD3<Double>>()
            for p in permutations { for sx in [-1.0,1] { for sy in [-1.0,1] { for sz in [-1.0,1] {
                points.insert(.init(sx*components[p[0]],sy*components[p[1]],sz*components[p[2]]))
            } } } }
            let sorted=points.sorted { a,b in a.x==b.x ? (a.y==b.y ? a.z<b.z:a.y<b.y):a.x<b.x }
            out+=sorted.map { ($0,weight) }
        }
        func axis(_ weight:Double) { orbit([1,0,0],weight) }
        func corner(_ weight:Double) { let a=sqrt(1.0/3);orbit([a,a,a],weight) }
        func aab(_ a:Double,_ weight:Double) { orbit([a,a,sqrt(1-2*a*a)],weight) }
        func ab0(_ a:Double,_ weight:Double) { orbit([a,sqrt(1-a*a),0],weight) }
        func abc(_ a:Double,_ b:Double,_ weight:Double) { orbit([a,b,sqrt(1-a*a-b*b)],weight) }
        switch count {
        case 50:
            axis(0.01269841269841270);orbit([sqrt(0.5),sqrt(0.5),0],0.02257495590828924)
            corner(0.02109375);aab(0.3015113445777636,0.02017333553791887)
        case 110:
            axis(0.003828270494937162);corner(0.009793737512487512)
            aab(0.1851156353447362,0.008211737283191111);aab(0.6904210483822922,0.009942814891178103)
            aab(0.3956894730559419,0.009595471336070963);ab0(0.4783690288121502,0.009694996361663028)
        case 302:
            axis(0.0008545911725128148);corner(0.003599119285025571)
            aab(0.3515640345570105,0.003449788424305883);aab(0.6566329410219612,0.003604822601419882)
            aab(0.4729054132581005,0.003576729661743367);aab(0.09618308522614784,0.002352101413689164)
            aab(0.2219645236294178,0.003108953122413675);aab(0.7011766416089545,0.003650045807677255)
            ab0(0.2644152887060663,0.002982344963171804);ab0(0.5718955891878961,0.003600820932216460)
            abc(0.2510034751770465,0.8000727494073952,0.003571540554273387)
            abc(0.1233548532583327,0.4127724083168531,0.003392312205006170)
        default: throw VivoChemistryError.unsupported("smooth solvent supports Lebedev 50,110,302")
        }
        guard out.count==count,abs(out.reduce(0.0) { $0+$1.1 }-1)<1e-12 else { throw VivoChemistryError.invalid("Lebedev rule normalization") };return out
    }
}
/// Smooth atom-sphere switching and finite Gaussian apparent surface charges.
/// Unlike point tesserae, nuclear/electronic potentials and the surface matrix
/// consistently use the same finite charge distributions. No non-electrostatic
/// free-energy terms or analytic nuclear derivatives are implied.
public struct VivoSmoothCPCMOperator:Sendable {
    public let configuration:VivoSmoothCPCMConfiguration
    private let one:VivoOneElectronIntegrals
    private let surface:VivoQMMatrix
    private let potentials:VivoQMMatrix
    private let nuclear:[Double]
    private let budget:VivoChemistryBudget
    public init(system:VivoElectronicSystem,basis:VivoGaussianBasis,configuration cfg:VivoSmoothCPCMConfiguration = .init(),budget:VivoChemistryBudget = .init()) throws {
        try system.validate();try budget.validate()
        guard cfg.dielectricConstant.isFinite,cfg.dielectricConstant>1,cfg.radiusScale.isFinite,cfg.radiusScale>0,
              cfg.tolerance.isFinite,cfg.tolerance>0,cfg.maximumIterations>0,cfg.maximumTesserae>0,
              cfg.radiiAngstrom.allSatisfy({$0.value.isFinite && $0.value>0}),
              cfg.extraCavitySpheres.allSatisfy({vivoQMFinite($0.centerBohr) && $0.radiusBohr.isFinite && $0.radiusBohr>0}) else { throw VivoChemistryError.invalid("smooth cavity configuration") }
        let grid=try VivoSolventLebedev.grid(cfg.angularPoints)
        var spheres=try system.nuclei.map { atom -> VivoSolventSphere in
            guard let radius=cfg.radiiAngstrom[atom.atomicNumber] else { throw VivoChemistryError.invalid("missing cavity radius") }
            return .init(centerBohr:atom.positionBohr,radiusBohr:radius*cfg.radiusScale/(10*VivoAtomicUnits.bohrInNM))
        };spheres+=cfg.extraCavitySpheres
        guard spheres.allSatisfy({vivoQMFinite($0.centerBohr) && $0.radiusBohr.isFinite && $0.radiusBohr>0}) else { throw VivoChemistryError.invalid("cavity radius overflow") }
        for i in spheres.indices { for j in 0..<i where vivoQMNorm(spheres[i].centerBohr-spheres[j].centerBohr)<1e-10 && abs(spheres[i].radiusBohr-spheres[j].radiusBohr)<1e-10 {
            throw VivoChemistryError.invalid("duplicate solvent cavity spheres")
        } }
        guard spheres.count<=cfg.maximumTesserae/cfg.angularPoints else { throw VivoChemistryError.resourceLimit("smooth cavity tessera count") }
        // External charges must be explicitly enclosed. An outlying charge is
        // never silently treated as if it were inside the solute cavity.
        for charge in system.pointCharges where charge.chargeE != 0 {
            guard spheres.contains(where:{vivoQMNorm(charge.positionBohr-$0.centerBohr)<$0.radiusBohr}) else { throw VivoChemistryError.unsupported("external charge lies outside declared combined cavity") }
        }
        let xiParameter:Double
        switch cfg.angularPoints { case 50:xiParameter=4.89250673295;case 110:xiParameter=4.90101060987;default:xiParameter=4.90498088169 }
        let widths=spheres.map { $0.radiusBohr*sqrt(14.0/Double(cfg.angularPoints)) }
        let inner=spheres.indices.map { i -> Double in
            let ratio=spheres[i].radiusBohr/widths[i],alpha=0.5+ratio-sqrt(ratio*ratio-1.0/28)
            return spheres[i].radiusBohr-alpha*widths[i]
        }
        var points:[SIMD3<Double>]=[],switches:[Double]=[],xis:[Double]=[]
        for (atom,sphere) in spheres.enumerated() { for (direction,weight) in grid {
            let point=sphere.centerBohr+sphere.radiusBohr*direction;var factor=1.0
            for j in spheres.indices where j != atom {
                let x=min(1,max(0,(vivoQMNorm(point-spheres[j].centerBohr)-inner[j])/widths[j]))
                factor*=x*x*x*(10-15*x+6*x*x)
            }
            if factor*weight>1e-16 { points.append(point);switches.append(factor);xis.append(xiParameter/(sphere.radiusBohr*sqrt(weight*4*Double.pi))) }
        } }
        let context=try VivoGaussianIntegralContext(system:system,basis:basis,budget:budget),n=context.orbitals.count,m=points.count
        guard m>0 else { throw VivoChemistryError.invalid("empty smooth cavity") }
        _ = try budget.elements([m,m],simultaneousArrays:3);_ = try budget.elements([m,n,n],simultaneousArrays:3)
        var s=VivoQMMatrix(m,m),p=VivoQMMatrix(m,n*n),nuclear=[Double](repeating:0,count:m),work=0
        for i in 0..<m {
            s[i,i]=xis[i]*sqrt(2/Double.pi)/switches[i]
            for j in 0..<i {
                let distance=vivoQMNorm(points[i]-points[j]),xi=xis[i]*xis[j]/sqrt(xis[i]*xis[i]+xis[j]*xis[j])
                let value=distance<1e-12 ? 2*xi/sqrt(Double.pi):erf(xi*distance)/distance
                s[i,j]=value;s[j,i]=value
            }
        }
        for a in 0..<m {
            let xi=xis[a],exponent=xi*xi,normalization=pow(exponent/Double.pi,1.5)
            func chargePotential(_ charge:Double,_ position:SIMD3<Double>)->Double {
                let d=vivoQMNorm(points[a]-position);return charge*(d<1e-12 ? 2*xi/sqrt(Double.pi):erf(xi*d)/d)
            }
            for nucleus in system.nuclei { nuclear[a]+=chargePotential(Double(nucleus.atomicNumber),nucleus.positionBohr) }
            for charge in system.pointCharges { nuclear[a]+=chargePotential(charge.chargeE,charge.positionBohr) }
            for mu in 0..<n { for nu in 0...mu {
                let left=context.orbitals[mu],right=context.orbitals[nu],ra=system.nuclei[left.nucleusIndex].positionBohr,rb=system.nuclei[right.nucleusIndex].positionBohr
                var value=0.0
                for i in left.weights.indices { for j in right.weights.indices {
                    guard work<budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("smooth solvent Gaussian integral work") };work+=1
                    value+=left.weights[i]*right.weights[j]*normalization*VivoGaussianIntegralEngine.primitiveERI(left.primitiveExponents[i],left.angular,ra,right.primitiveExponents[j],right.angular,rb,exponent,[0,0,0],points[a],0,[0,0,0],points[a])
                } };p[a,mu*n+nu]=value;p[a,nu*n+mu]=value
            } }
        }
        guard s.values.allSatisfy(\.isFinite),p.values.allSatisfy(\.isFinite),nuclear.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("smooth solvent finite integrals") }
        self.configuration=cfg;self.one=try context.oneElectron();self.surface=s;self.potentials=p;self.nuclear=nuclear;self.budget=budget
    }
    public func evaluate(totalDensity d:VivoQMMatrix) throws -> VivoSmoothCPCMResult {
        let n=one.count,m=surface.rows,cfg=configuration
        guard d.rows==n,d.columns==n,d.values.allSatisfy(\.isFinite),try d.adding(d.transposed,scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("smooth PCM density") }
        guard try VivoQMDenseAlgebra.symmetricEigen(d).values.allSatisfy({$0>=(-1e-9)}) else { throw VivoChemistryError.invalid("nonpositive solvent density") }
        let electrons=zip(d.values,one.overlap.values).reduce(0.0) { $0+$1.0*$1.1 }
        guard abs(electrons-Double(one.system.alphaElectrons+one.system.betaElectrons))<1e-6 else { throw VivoChemistryError.invalid("smooth PCM density charge") }
        let electronic=try potentials.multiplied(by:VivoQMMatrix(rows:n*n,columns:1,values:d.values)).values
        let potential=zip(nuclear,electronic).map { $0-$1 },scale=(cfg.dielectricConstant-1)/cfg.dielectricConstant,rhs=potential.map { -scale*$0 }
        func dot(_ a:[Double],_ b:[Double])->Double { zip(a,b).reduce(0.0) { $0+$1.0*$1.1 } }
        func apply(_ q:[Double]) throws -> [Double] { try surface.multiplied(by:VivoQMMatrix(rows:m,columns:1,values:q)).values }
        var q=[Double](repeating:0,count:m),r=rhs,z=(0..<m).map { rhs[$0]/surface[$0,$0] },direction=z,rz=dot(r,z)
        let norm=rhs.reduce(0.0) { hypot($0,$1) }
        var iterations=0
        let maxSteps=min(cfg.maximumIterations,budget.maximumOperatorApplications/max(1,m*m)-1)
        guard maxSteps>0 else { throw VivoChemistryError.resourceLimit("smooth solvent linear solve work") }
        if norm>0 { for step in 1...maxSteps {
            let ap=try apply(direction),denominator=dot(direction,ap)
            guard denominator.isFinite,denominator>0,rz.isFinite,rz>0 else { throw VivoChemistryError.convergence("smooth solvent PCG definiteness") }
            let alpha=rz/denominator
            for i in 0..<m { q[i]+=alpha*direction[i];r[i]-=alpha*ap[i] };iterations=step
            if r.reduce(0.0,{hypot($0,$1)})<=cfg.tolerance*norm { break }
            z=(0..<m).map { r[$0]/surface[$0,$0] };let next=dot(r,z),beta=next/rz
            for i in 0..<m { direction[i]=z[i]+beta*direction[i] };rz=next
        } }
        let trueResidual=zip(try apply(q),rhs).reduce(0.0) { hypot($0,$1.0-$1.1) }/max(norm,1e-300)
        guard trueResidual.isFinite,trueResidual<=5*cfg.tolerance else { throw VivoChemistryError.convergence("smooth solvent true surface residual") }
        let row=try VivoQMMatrix(rows:1,columns:m,values:q.map { -$0 }),response=try row.multiplied(by:potentials)
        let reaction=try VivoQMMatrix(rows:n,columns:n,values:response.values),energy=0.5*dot(q,potential)
        guard energy.isFinite,q.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite smooth solvent energy/charges") }
        let solute=Double(one.system.nuclei.reduce(0) { $0+$1.atomicNumber })+one.system.pointCharges.reduce(0) { $0+$1.chargeE }-electrons
        return .init(configuration:cfg,polarizationEnergyHartree:energy,reactionPotentialMatrix:reaction,apparentSurfaceCharges:q,
            trueRelativeResidual:trueResidual,surfaceChargeDefect:q.reduce(0,+)+scale*solute,iterations:iterations,tesseraCount:m,
            method:"C-PCM; smooth switching/Gaussian surface charges; Lebedev quadrature; electrostatic polarization only")
    }
}
public struct VivoSmoothSolventSCFResult:Codable,Sendable,Equatable {
    public let scf:VivoHartreeFockResult
    public let solvent:VivoSmoothCPCMResult
}
public enum VivoSmoothSolventSCF {
    public static func solve(system:VivoElectronicSystem,integrals ao:VivoAOIntegrals,
                             solvent:VivoSmoothCPCMConfiguration = .init(),scf:VivoSCFConfiguration = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoSmoothSolventSCFResult {
        try ao.validate(budget:budget)
        guard system==ao.sourceSystem else { throw VivoChemistryError.invalid("solvated SCF source binding") }
        let op=try VivoSmoothCPCMOperator(system:system,basis:ao.sourceBasis,configuration:solvent,budget:budget)
        var latest:VivoSmoothCPCMResult?
        let result=try VivoHartreeFock.solveWithFockBuilder(system:system,overlap:ao.overlap,core:ao.coreHamiltonian,
            constant:ao.constantEnergyHartree,configuration:scf,budget:budget,energyFunctional:{ da,db,_,_ in
                guard let latest else { throw VivoChemistryError.invalid("solvent response missing") }
                let gas=VivoHartreeFock.focks(ao,da,db)
                return VivoHartreeFock.energy(ao,da,db,gas.0,gas.1)+latest.polarizationEnergyHartree
            }) { da,db in
                let pcm=try op.evaluate(totalDensity:da.adding(db));latest=pcm
                let gas=VivoHartreeFock.focks(ao,da,db)
                return (try gas.0.adding(pcm.reactionPotentialMatrix),try gas.1.adding(pcm.reactionPotentialMatrix))
            }
        return .init(scf:result,solvent:try op.evaluate(totalDensity:result.alphaDensity.adding(result.betaDensity)))
    }
}
