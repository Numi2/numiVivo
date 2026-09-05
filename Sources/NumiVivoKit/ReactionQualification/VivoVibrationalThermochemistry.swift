import Foundation

public enum VivoStationaryKind: String, Codable, Sendable { case minimum, firstOrderSaddle }
public enum VivoThermoStandardState: Codable, Sendable, Equatable {
    case idealGas(pressurePa: Double)
    /// Formal ideal-gas -> concentration standard-state conversion. This does
    /// not supply non-electrostatic solvent or solution configurational entropy.
    case concentration(molPerLitre: Double, gasReferencePressurePa: Double)
}
public struct VivoThermochemistryConfiguration: Codable, Sendable, Equatable {
    public var temperatureK: Double
    public var standardState: VivoThermoStandardState
    public var rotationalSymmetryNumber: Int
    public var electronicDegeneracy: Int
    public var minimumResolvedFrequencyCM: Double
    public var gradientRMSTolerance: Double
    public var maximumGradientTolerance: Double
    public init(temperatureK: Double = 298.15, standardState: VivoThermoStandardState = .idealGas(pressurePa:101325),
                rotationalSymmetryNumber: Int = 1, electronicDegeneracy: Int = 1,
                minimumResolvedFrequencyCM: Double = 5, gradientRMSTolerance: Double = 1e-5,
                maximumGradientTolerance: Double = 3e-5) {
        self.temperatureK=temperatureK; self.standardState=standardState
        self.rotationalSymmetryNumber=rotationalSymmetryNumber; self.electronicDegeneracy=electronicDegeneracy
        self.minimumResolvedFrequencyCM=minimumResolvedFrequencyCM; self.gradientRMSTolerance=gradientRMSTolerance
        self.maximumGradientTolerance=maximumGradientTolerance
    }
    public func validate() throws {
        guard temperatureK.isFinite,temperatureK>0,temperatureK<=10000,
              (1...1000000).contains(rotationalSymmetryNumber),(1...10000).contains(electronicDegeneracy),
              [minimumResolvedFrequencyCM,gradientRMSTolerance,maximumGradientTolerance].allSatisfy({$0.isFinite && $0>0}) else {
            throw VivoChemistryError.invalid("thermochemical temperature, symmetry, degeneracy or tolerance")
        }
        switch standardState {
        case .idealGas(let p): guard p.isFinite,p>0 else {throw VivoChemistryError.invalid("gas standard pressure")}
        case .concentration(let c,let p): guard c.isFinite,c>0,p.isFinite,p>0 else {throw VivoChemistryError.invalid("concentration standard state")}
        }
    }
}
public struct VivoVibrationalAnalysis: Codable, Sendable, Equatable {
    public let positionsBohr: [SIMD3<Double>]
    public let massesDa: [Double]
    public let cartesianHessian: VivoQMMatrix
    public let rigidModeCount: Int
    /// All 3N-rigid frequencies are retained. Negative means imaginary. Small
    /// modes are NOT silently removed as if they were translations/rotations.
    public let signedFrequenciesCM: [Double]
    public let inertiaEigenvaluesDaBohr2: [Double]
    /// Mass-weighted Cartesian eigenvectors, columns in the vibrational space.
    public let massWeightedModes: VivoQMMatrix
    public let maximumGradient: Double
    public let gradientRMS: Double
    public let rigidHessianResidual: Double
}
public enum VivoNuclearUnits {
    public static let avogadro = 6.02214076e23
    public static let atomicMassKG = 1.66053906892e-27
    public static let bohrM = VivoAtomicUnits.bohrInNM*1e-9
    public static let hartreeJ = VivoAtomicUnits.hartreeInKJPerMol*1000/avogadro
    public static let lightCMPerS = 2.99792458e10
    public static let kHartree = VivoAtomicUnits.boltzmannJPerK/hartreeJ
}
public enum VivoNormalModes {
    private static func dot(_ a: [Double],_ b: [Double]) -> Double {zip(a,b).reduce(0.0){$0+$1.0*$1.1}}
    static func orthogonalized(_ input: [Double], against basis: [[Double]]) -> [Double]? {
        var v=input
        for _ in 0..<2 {for q in basis {let scale=dot(v,q);for i in v.indices {v[i]-=scale*q[i]}}}
        let norm=v.reduce(0.0){hypot($0,$1)}
        return norm>1e-10 ? v.map{$0/norm} : nil
    }
    public static func analyze(positions: [SIMD3<Double>], masses: [Double], evaluation: VivoGeometryEvaluation,
                               hessian: VivoQMMatrix, budget: VivoChemistryBudget = .init()) throws -> VivoVibrationalAnalysis {
        try budget.validate()
        let count=positions.count,n=3*count
        guard count>0,count<=budget.maximumBasisFunctions,masses.count==count,
              positions.allSatisfy(vivoQMFinite),masses.allSatisfy({$0.isFinite && $0>0}),
              evaluation.energyHartree.isFinite,evaluation.gradientHartreePerBohr.count==count,
              evaluation.gradientHartreePerBohr.allSatisfy(vivoQMFinite),hessian.rows==n,hessian.columns==n,
              hessian.values.allSatisfy(\.isFinite),try hessian.adding(hessian.transposed,scale:-1).frobeniusNorm<1e-7 else {
            throw VivoChemistryError.invalid("normal-mode coordinates, isotopic masses, gradient or Hessian")
        }
        _=try budget.elements([n,n],simultaneousArrays:12)
        let total=masses.reduce(0,+)
        var com=SIMD3<Double>.zero
        for i in positions.indices {com+=masses[i]*positions[i]};com/=total
        let centered=positions.map{$0-com}
        var inertia=VivoQMMatrix(3,3)
        for (i,r) in centered.enumerated() {for a in 0..<3 {for b in 0..<3 {inertia[a,b]+=masses[i]*((a==b ? vivoQMDot(r,r):0)-r[a]*r[b])}}}
        let eig=try VivoQMDenseAlgebra.symmetricEigen(inertia,tolerance:1e-13)
        let scale=max(1e-30,eig.values.last!),linear=count>1 && eig.values[0]<=scale*1e-10
        if count>1 {guard eig.values[1]>1e-12 else {throw VivoChemistryError.invalid("collapsed nuclear geometry")}}
        let rigidCount=count==1 ? 3 : (linear ? 5:6)
        var rigid:[[Double]]=[]
        for axis in 0..<3 {
            var v=[Double](repeating:0,count:n)
            for i in 0..<count {v[3*i+axis]=sqrt(masses[i])}
            if let q=orthogonalized(v,against:rigid) {rigid.append(q)}
        }
        if count>1 {for axis in 0..<3 where eig.values[axis]>scale*1e-10 {
            let a=SIMD3<Double>(eig.vectors[0,axis],eig.vectors[1,axis],eig.vectors[2,axis])
            var v=[Double](repeating:0,count:n)
            for (i,r) in centered.enumerated() {
                let cross=SIMD3<Double>(a.y*r.z-a.z*r.y,a.z*r.x-a.x*r.z,a.x*r.y-a.y*r.x)*sqrt(masses[i])
                for j in 0..<3 {v[3*i+j]=cross[j]}
            }
            if let q=orthogonalized(v,against:rigid) {rigid.append(q)}
        }}
        guard rigid.count==rigidCount else {throw VivoChemistryError.convergence("rigid-mode rank ambiguity")}
        var basis=rigid, vibrational:[[Double]]=[]
        for j in 0..<n {
            var e=[Double](repeating:0,count:n);e[j]=1
            if let q=orthogonalized(e,against:basis) {basis.append(q);vibrational.append(q)}
        }
        guard vibrational.count==n-rigidCount else {throw VivoChemistryError.convergence("vibrational complement rank")}
        var mw=hessian
        for i in 0..<n {for j in 0..<n {mw[i,j]/=sqrt(masses[i/3]*masses[j/3])}}
        var rigidResidual=0.0
        for q in rigid {for i in 0..<n {let a=(0..<n).reduce(0.0){$0+mw[i,$1]*q[$1]};rigidResidual=hypot(rigidResidual,a)}}
        var q=VivoQMMatrix(n,vibrational.count)
        for j in vibrational.indices {for i in 0..<n {q[i,j]=vibrational[j][i]}}
        var frequencies:[Double]=[],modes=q
        if !vibrational.isEmpty {
            let reduced=try mw.congruence(q),spectrum=try VivoQMDenseAlgebra.symmetricEigen(reduced,tolerance:1e-13,maximumSweeps:256)
            let conversion=sqrt(VivoNuclearUnits.hartreeJ/(VivoNuclearUnits.atomicMassKG*pow(VivoNuclearUnits.bohrM,2)))/(2*Double.pi*VivoNuclearUnits.lightCMPerS)
            frequencies=spectrum.values.map{($0<0 ? -1:1)*sqrt(abs($0))*conversion}
            modes=try q.multiplied(by:spectrum.vectors)
        }
        let g=evaluation.gradientHartreePerBohr.flatMap{[$0.x,$0.y,$0.z]}
        return .init(positionsBohr:positions,massesDa:masses,cartesianHessian:hessian,rigidModeCount:rigidCount,
            signedFrequenciesCM:frequencies,inertiaEigenvaluesDaBohr2:eig.values,massWeightedModes:modes,
            maximumGradient:g.map(abs).max()!,gradientRMS:sqrt(g.reduce(0){$0+$1*$1}/Double(n)),rigidHessianResidual:rigidResidual)
    }
}
public struct VivoThermochemistryResult: Codable, Sendable, Equatable {
    public let configuration: VivoThermochemistryConfiguration
    public let kind: VivoStationaryKind
    public let electronicEnergyHartree: Double
    public let zeroPointEnergyHartree: Double
    /// Includes ZPE, finite-temperature vibrations, ideal translation/rotation.
    public let enthalpyCorrectionHartree: Double
    public let gasEntropyHartreePerK: Double
    public let standardStateShiftHartree: Double
    public let gibbsEnergyHartree: Double
    public let imaginaryFrequencyCM: Double?
    public let rotor: String
    public let modes: VivoVibrationalAnalysis
    public let model: String
}
public enum VivoRRHOThermochemistry {
    public static func evaluate(energyHartree: Double, modes: VivoVibrationalAnalysis, kind: VivoStationaryKind,
                                configuration cfg: VivoThermochemistryConfiguration = .init()) throws -> VivoThermochemistryResult {
        try cfg.validate()
        guard energyHartree.isFinite,!modes.positionsBohr.isEmpty,modes.positionsBohr.allSatisfy(vivoQMFinite),
              modes.massesDa.count==modes.positionsBohr.count,modes.massesDa.allSatisfy({$0.isFinite && $0>0}),
              [3,5,6].contains(modes.rigidModeCount),modes.inertiaEigenvaluesDaBohr2.count==3,
              modes.inertiaEigenvaluesDaBohr2.allSatisfy(\.isFinite),
              (modes.positionsBohr.count==1 ? modes.rigidModeCount==3 : (modes.rigidModeCount>=5 && modes.inertiaEigenvaluesDaBohr2[1]>0)),
              modes.maximumGradient.isFinite,modes.maximumGradient>=0,modes.gradientRMS.isFinite,modes.gradientRMS>=0,
              modes.maximumGradient<=cfg.maximumGradientTolerance,modes.gradientRMS<=cfg.gradientRMSTolerance,
              modes.signedFrequenciesCM.count==3*modes.positionsBohr.count-modes.rigidModeCount,
              modes.signedFrequenciesCM.allSatisfy({$0.isFinite && abs($0)>=cfg.minimumResolvedFrequencyCM}) else {
            throw VivoChemistryError.convergence("RRHO requires a stationary geometry and every vibrational mode resolved; no soft-mode deletion")
        }
        let negative=modes.signedFrequenciesCM.filter{$0<0}
        guard negative.count==(kind == .minimum ? 0:1) else {throw VivoChemistryError.invalid("wrong Hessian index for requested stationary kind")}
        let k=VivoNuclearUnits.kHartree,T=cfg.temperatureK,h=VivoAtomicUnits.planckJS,kSI=VivoAtomicUnits.boltzmannJPerK
        let pressure:Double,shift:Double
        switch cfg.standardState {
        case .idealGas(let p): pressure=p;shift=0
        case .concentration(let c,let p):
            pressure=p;shift=k*T*log(c*1000*VivoAtomicUnits.gasConstantJPerMolK*T/p)
        }
        let mass=modes.massesDa.reduce(0,+)*VivoNuclearUnits.atomicMassKG
        let logTrans=1.5*log(2*Double.pi*mass*kSI*T/(h*h))+log(kSI*T/pressure)
        var entropy=k*(2.5+logTrans+log(Double(cfg.electronicDegeneracy))),enthalpy=2.5*k*T
        let rotor:String,rotationEnergy:Double,logRotation:Double
        let moments=modes.inertiaEigenvaluesDaBohr2.map{$0*VivoNuclearUnits.atomicMassKG*pow(VivoNuclearUnits.bohrM,2)}
        if modes.positionsBohr.count==1 {rotor="atom";rotationEnergy=0;logRotation=0}
        else if modes.rigidModeCount==5 {
            rotor="linear";rotationEnergy=1
            logRotation=log(8*Double.pi*Double.pi*kSI*T*sqrt(moments[1]*moments[2])/(Double(cfg.rotationalSymmetryNumber)*h*h))
        } else {
            rotor="nonlinear";rotationEnergy=1.5
            logRotation=0.5*log(Double.pi)-log(Double(cfg.rotationalSymmetryNumber))+1.5*log(8*Double.pi*Double.pi*kSI*T/(h*h))+0.5*moments.reduce(0){$0+log($1)}
        }
        entropy+=k*(rotationEnergy+logRotation);enthalpy+=rotationEnergy*k*T
        var zpe=0.0
        for frequency in modes.signedFrequenciesCM where frequency>0 {
            let quantum=h*VivoNuclearUnits.lightCMPerS*frequency/VivoNuclearUnits.hartreeJ,x=quantum/(k*T),decay=exp(-x)
            let denominator = -expm1(-x)
            zpe+=0.5*quantum
            enthalpy+=0.5*quantum+quantum*decay/denominator
            entropy+=k*(x*decay/denominator-log(denominator))
        }
        let gibbs=energyHartree+enthalpy-T*entropy+shift
        guard [zpe,enthalpy,entropy,shift,gibbs].allSatisfy(\.isFinite) else {throw VivoChemistryError.convergence("RRHO numerical overflow")}
        return .init(configuration:cfg,kind:kind,electronicEnergyHartree:energyHartree,zeroPointEnergyHartree:zpe,
            enthalpyCorrectionHartree:enthalpy,gasEntropyHartreePerK:entropy,standardStateShiftHartree:shift,
            gibbsEnergyHartree:gibbs,imaginaryFrequencyCM:negative.first,rotor:rotor,modes:modes,
            model:"ideal RRHO; TS unstable mode excluded, not made positive; no tunneling, hindered rotors, configurational entropy or non-electrostatic solvation")
    }
}
