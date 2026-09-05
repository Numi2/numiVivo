import Foundation

/// One-electron data without allocating the four-index ERI tensor.
public struct VivoOneElectronIntegrals: Codable, Sendable, Equatable {
    public let system: VivoElectronicSystem
    public let basis: VivoGaussianBasis
    public let orbitals: [VivoCartesianOrbital]
    public let overlap: VivoQMMatrix
    public let coreHamiltonian: VivoQMMatrix
    public let constantEnergyHartree: Double
    public var count: Int { overlap.rows }
    public func validate(budget: VivoChemistryBudget = .init()) throws {
        let n=count
        try system.validate();try basis.validate(nucleusCount:system.nuclei.count)
        guard n>0,n<=budget.maximumBasisFunctions,overlap.columns==n,coreHamiltonian.rows==n,coreHamiltonian.columns==n,
              overlap.values.allSatisfy(\.isFinite),coreHamiltonian.values.allSatisfy(\.isFinite),constantEnergyHartree.isFinite,
              try VivoGaussianIntegralEngine.expanded(system:system,basis:basis,budget:budget)==orbitals,
              try overlap.adding(overlap.transposed,scale:-1).frobeniusNorm<1e-10,
              try coreHamiltonian.adding(coreHamiltonian.transposed,scale:-1).frobeniusNorm<1e-10 else {
            throw VivoChemistryError.invalid("one-electron integral source, dimensions or symmetry")
        }
    }
}
public struct VivoGaussianIntegralContext: Sendable {
    public let system: VivoElectronicSystem
    public let basis: VivoGaussianBasis
    public let orbitals: [VivoCartesianOrbital]
    public let budget: VivoChemistryBudget
    public init(system:VivoElectronicSystem,basis:VivoGaussianBasis,budget:VivoChemistryBudget = .init()) throws {
        orbitals=try VivoGaussianIntegralEngine.expanded(system:system,basis:basis,budget:budget)
        self.system=system;self.basis=basis;self.budget=budget
    }
    public func oneElectron() throws -> VivoOneElectronIntegrals {
        let n=orbitals.count;_ = try budget.elements([n,n],simultaneousArrays:8)
        var overlap=VivoQMMatrix(n,n),core=overlap,work=0
        for p in 0..<n { for q in 0...p {
            let a=orbitals[p],b=orbitals[q],ra=system.nuclei[a.nucleusIndex].positionBohr,rb=system.nuclei[b.nucleusIndex].positionBohr
            var s=0.0,h=0.0
            for i in a.weights.indices { for j in b.weights.indices {
                let tasks=system.nuclei.count+system.pointCharges.count+2
                guard tasks<=budget.maximumOperatorApplications-work else { throw VivoChemistryError.resourceLimit("one-electron primitive work") };work+=tasks
                let ea=a.primitiveExponents[i],eb=b.primitiveExponents[j],w=a.weights[i]*b.weights[j]
                s+=w*VivoGaussianIntegralEngine.primitiveOverlap(ea,a.angular,ra,eb,b.angular,rb)
                h+=w*VivoGaussianIntegralEngine.primitiveKinetic(ea,a.angular,ra,eb,b.angular,rb)
                for atom in system.nuclei { h-=w*Double(atom.atomicNumber)*VivoGaussianIntegralEngine.primitivePotential(ea,a.angular,ra,eb,b.angular,rb,atom.positionBohr) }
                for charge in system.pointCharges where charge.chargeE != 0 { h-=w*charge.chargeE*VivoGaussianIntegralEngine.primitivePotential(ea,a.angular,ra,eb,b.angular,rb,charge.positionBohr) }
            } }
            overlap[p,q]=s;overlap[q,p]=s;core[p,q]=h;core[q,p]=h
        } }
        var scalar=0.0
        for (i,a) in system.nuclei.enumerated() {
            for j in 0..<i { scalar+=Double(a.atomicNumber*system.nuclei[j].atomicNumber)/vivoQMNorm(a.positionBohr-system.nuclei[j].positionBohr) }
            for q in system.pointCharges where q.chargeE != 0 { scalar+=Double(a.atomicNumber)*q.chargeE/vivoQMNorm(a.positionBohr-q.positionBohr) }
        }
        let result=VivoOneElectronIntegrals(system:system,basis:basis,orbitals:orbitals,overlap:overlap,coreHamiltonian:core,constantEnergyHartree:scalar)
        try result.validate(budget:budget);return result
    }
    /// A missing leg denotes the constant function, not a zero-exponent basis
    /// function. This is the exact reduction to two/three-centre Coulomb integrals.
    func coulomb(_ a:VivoCartesianOrbital,_ b:VivoCartesianOrbital?,_ c:VivoCartesianOrbital,_ d:VivoCartesianOrbital?,work:inout Int) throws -> Double {
        let ra=system.nuclei[a.nucleusIndex].positionBohr,rb=b.map { system.nuclei[$0.nucleusIndex].positionBohr } ?? ra
        let rc=system.nuclei[c.nucleusIndex].positionBohr,rd=d.map { system.nuclei[$0.nucleusIndex].positionBohr } ?? rc
        let bw=b?.weights ?? [1],dw=d?.weights ?? [1];var value=0.0
        for i in a.weights.indices { for j in bw.indices { for k in c.weights.indices { for l in dw.indices {
            guard work<budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("Gaussian Coulomb work") };work+=1
            value+=a.weights[i]*bw[j]*c.weights[k]*dw[l]*VivoGaussianIntegralEngine.primitiveERI(
                a.primitiveExponents[i],a.angular,ra,b?.primitiveExponents[j] ?? 0,b?.angular ?? [0,0,0],rb,
                c.primitiveExponents[k],c.angular,rc,d?.primitiveExponents[l] ?? 0,d?.angular ?? [0,0,0],rd)
        } } } }
        guard value.isFinite else { throw VivoChemistryError.convergence("Gaussian Coulomb overflow") };return value
    }
}

/// Unweighted packed triangular pairs p>=q. B[pq,L] B[rs,L] represents (pq|rs).
public struct VivoCoulombFactors: Codable, Sendable, Equatable {
    public let orbitalCount:Int
    public let values:VivoQMMatrix
    public let method:String
    public var rank:Int { values.columns }
    public init(orbitalCount:Int,values:VivoQMMatrix,method:String) { self.orbitalCount=orbitalCount;self.values=values;self.method=method }
    public static func pair(_ p:Int,_ q:Int)->Int { let a=max(p,q);return a*(a+1)/2+min(p,q) }
    public func validate(budget:VivoChemistryBudget = .init()) throws {
        try budget.validate();let n=orbitalCount
        guard n>0,n<=budget.maximumBasisFunctions,rank>0,values.rows==n*(n+1)/2,!method.isEmpty,values.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("Coulomb factor shape/identity") }
        _ = try budget.elements([values.rows,rank],simultaneousArrays:2)
    }
    public func eri(_ p:Int,_ q:Int,_ r:Int,_ s:Int)->Double {
        let a=Self.pair(p,q),b=Self.pair(r,s);return (0..<rank).reduce(0.0) { $0+values[a,$1]*values[b,$1] }
    }
    func matrix(_ k:Int)->VivoQMMatrix {
        var result=VivoQMMatrix(orbitalCount,orbitalCount)
        for p in 0..<orbitalCount { for q in 0...p { result[p,q]=values[Self.pair(p,q),k];result[q,p]=result[p,q] } };return result
    }
    public func transformed(by c:VivoQMMatrix,budget:VivoChemistryBudget = .init()) throws -> Self {
        try validate(budget:budget);let m=c.columns
        guard c.rows==orbitalCount,m>0,m<=budget.maximumBasisFunctions,c.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("factor transform") }
        _ = try budget.elements([m*(m+1)/2,rank],simultaneousArrays:3);var output=VivoQMMatrix(m*(m+1)/2,rank)
        for k in 0..<rank { let b=try matrix(k).congruence(c)
            for p in 0..<m { for q in 0...p { output[Self.pair(p,q),k]=b[p,q] } }
        };return .init(orbitalCount:m,values:output,method:method+"; orbital transformed")
    }
    public func coulombExchange(density d:VivoQMMatrix,budget:VivoChemistryBudget = .init()) throws -> (coulomb:VivoQMMatrix,exchange:VivoQMMatrix) {
        try validate(budget:budget);let n=orbitalCount
        guard d.rows==n,d.columns==n,d.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("RI density dimensions") }
        _ = try budget.elements([n,n],simultaneousArrays:12)
        var j=VivoQMMatrix(n,n),k=j
        for l in 0..<rank {
            let b=matrix(l),charge=zip(b.values,d.values).reduce(0.0) { $0+$1.0*$1.1 }
            j=try j.adding(b,scale:charge);k=try k.adding(b.multiplied(by:d).multiplied(by:b))
        };return (j,k)
    }
}
public struct VivoRIConfiguration: Codable, Sendable, Equatable {
    public var metricEigenvalueThreshold:Double
    public var maximumAuxiliaryFunctions:Int
    public var pairBlockSize:Int
    public init(metricEigenvalueThreshold:Double=1e-10,maximumAuxiliaryFunctions:Int=4096,pairBlockSize:Int=256) {
        self.metricEigenvalueThreshold=metricEigenvalueThreshold;self.maximumAuxiliaryFunctions=maximumAuxiliaryFunctions;self.pairBlockSize=pairBlockSize
    }
}
public struct VivoRIIntegrals: Codable, Sendable, Equatable {
    public let oneElectron:VivoOneElectronIntegrals
    public let auxiliaryBasis:VivoGaussianBasis
    public let factors:VivoCoulombFactors
    public let configuration:VivoRIConfiguration
    public let auxiliaryCount:Int
    public let whiteningDefect:Double
    public let primitiveApplications:Int
}
public enum VivoDensityFitting {
    /// Native RI: two/three-centre integrals -> rank-revealing metric whitening
    /// -> packed factor panels. The dense n^4 ERI tensor is never constructed.
    public static func compute(system:VivoElectronicSystem,basis:VivoGaussianBasis,auxiliaryBasis:VivoGaussianBasis,
                               configuration cfg:VivoRIConfiguration = .init(),budget:VivoChemistryBudget = .init()) throws -> VivoRIIntegrals {
        guard cfg.metricEigenvalueThreshold.isFinite,cfg.metricEigenvalueThreshold>0,
              (1...8192).contains(cfg.maximumAuxiliaryFunctions),cfg.pairBlockSize>0 else { throw VivoChemistryError.invalid("RI configuration") }
        let context=try VivoGaussianIntegralContext(system:system,basis:basis,budget:budget)
        var auxiliaryBudget=budget;auxiliaryBudget.maximumBasisFunctions=cfg.maximumAuxiliaryFunctions
        let aux=try VivoGaussianIntegralEngine.expanded(system:system,basis:auxiliaryBasis,budget:auxiliaryBudget)
        let n=context.orbitals.count,m=aux.count,pairs=n*(n+1)/2
        _ = try budget.elements([m,m],simultaneousArrays:6);var metric=VivoQMMatrix(m,m),work=0
        for p in 0..<m { for q in 0...p { metric[p,q]=try context.coulomb(aux[p],nil,aux[q],nil,work:&work);metric[q,p]=metric[p,q] } }
        let eig=try VivoQMDenseAlgebra.symmetricEigen(metric,tolerance:1e-14,maximumSweeps:150)
        guard eig.values.allSatisfy({$0>=(-1e-10)}) else { throw VivoChemistryError.invalid("indefinite RI metric") }
        let keep=eig.values.indices.filter { eig.values[$0]>cfg.metricEigenvalueThreshold }
        guard !keep.isEmpty else { throw VivoChemistryError.invalid("empty RI metric rank") }
        _ = try budget.elements([pairs,keep.count],simultaneousArrays:3)
        var w=VivoQMMatrix(m,keep.count)
        for (j,k) in keep.enumerated() { for i in 0..<m { w[i,j]=eig.vectors[i,k]/sqrt(eig.values[k]) } }
        let defect=try metric.congruence(w).adding(.identity(keep.count),scale:-1).values.map(abs).max() ?? 0
        guard defect<1e-5 else { throw VivoChemistryError.convergence("RI whitening defect") }
        var result=VivoQMMatrix(pairs,keep.count)
        let indices=(0..<n).flatMap { p in (0...p).map { (p,$0) } }
        for start in stride(from:0,to:pairs,by:cfg.pairBlockSize) {
            let length=min(cfg.pairBlockSize,pairs-start);_ = try budget.elements([length,m],simultaneousArrays:3)
            var block=VivoQMMatrix(length,m)
            for i in 0..<length { let (p,q)=indices[start+i]
                for a in 0..<m { block[i,a]=try context.coulomb(context.orbitals[p],context.orbitals[q],aux[a],nil,work:&work) }
            }
            let whitened=try block.multiplied(by:w)
            for i in 0..<length { for k in keep.indices { result[start+i,k]=whitened[i,k] } }
        }
        let factors=VivoCoulombFactors(orbitalCount:n,values:result,method:"RI Coulomb metric; explicit normalized Cartesian auxiliary basis")
        try factors.validate(budget:budget)
        return .init(oneElectron:try context.oneElectron(),auxiliaryBasis:auxiliaryBasis,factors:factors,configuration:cfg,
                     auxiliaryCount:m,whiteningDefect:defect,primitiveApplications:work)
    }
}
public struct VivoRIHartreeFockResult: Codable, Sendable, Equatable {
    public let source:VivoRIIntegrals
    public let configuration:VivoSCFConfiguration
    public let scf:VivoHartreeFockResult
}
public enum VivoFactorizedHartreeFock {
    public static func solve(_ ri:VivoRIIntegrals,configuration:VivoSCFConfiguration = .init(),budget:VivoChemistryBudget = .init()) throws -> VivoRIHartreeFockResult {
        let one=ri.oneElectron;try one.validate(budget:budget);try ri.factors.validate(budget:budget)
        guard ri.factors.orbitalCount==one.count else { throw VivoChemistryError.invalid("RI source dimension") }
        let result=try VivoHartreeFock.solveWithFockBuilder(system:one.system,overlap:one.overlap,core:one.coreHamiltonian,
            constant:one.constantEnergyHartree,configuration:configuration,budget:budget) { da,db in
            let a=try ri.factors.coulombExchange(density:da,budget:budget),b=try ri.factors.coulombExchange(density:db,budget:budget)
            let common=try one.coreHamiltonian.adding(a.coulomb).adding(b.coulomb)
            return (try common.adding(a.exchange,scale:-1),try common.adding(b.exchange,scale:-1))
        };return .init(source:ri,configuration:configuration,scf:result)
    }
}
public enum VivoFactorizedMP2 {
    /// Occupied-virtual factor panels and a single virtual-pair contraction;
    /// neither an MO n^4 ERI nor a complete t2 tensor is allocated.
    public static func solve(_ ref:VivoRIHartreeFockResult,minimumGapHartree:Double=1e-8,budget:VivoChemistryBudget = .init()) throws -> VivoMP2Result {
        let hf=ref.scf,one=ref.source.oneElectron,factors=ref.source.factors
        try one.validate(budget:budget);try factors.validate(budget:budget)
        let n=one.count,o=hf.alphaElectrons,v=n-o,r=factors.rank,c=hf.alphaCoefficients
        guard hf.reference == .restricted,o==hf.betaElectrons,o>=0,o<=n,o==one.system.alphaElectrons,
              hf.betaElectrons==one.system.betaElectrons,c.rows==n,c.columns==n,hf.alphaOrbitalEnergies.count==n,
              hf.alphaOrbitalEnergies.allSatisfy(\.isFinite),minimumGapHartree.isFinite,minimumGapHartree>0,
              try one.overlap.congruence(c).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("DF-MP2 reference") }
        let d=VivoHartreeFock.density(c,occupied:o),jk=try factors.coulombExchange(density:d,budget:budget)
        let f=try one.coreHamiltonian.adding(jk.coulomb,scale:2).adding(jk.exchange,scale:-1),fmo=try f.congruence(c)
        for i in 0..<n { for j in 0..<n {
            guard abs(fmo[i,j]-(i==j ? hf.alphaOrbitalEnergies[i]:0))<1e-7 else { throw VivoChemistryError.invalid("DF-MP2 canonicality") }
        } }
        let energy=one.constantEnergyHartree+zip(d.values,try one.coreHamiltonian.adding(f).values).reduce(0.0) { $0+$1.0*$1.1 }
        guard hf.energyHartree.isFinite,abs(energy-hf.energyHartree)<1e-7 else { throw VivoChemistryError.invalid("DF-MP2 reference energy") }
        _ = try budget.elements([max(1,o*v),r],simultaneousArrays:3)
        _ = try budget.elements([max(1,v),max(1,v)],simultaneousArrays:4)
        let work=try budget.elements([o,o,v,v,max(1,r)])
        guard work<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("DF-MP2 contraction work") }
        var ia=VivoQMMatrix(o*v,r)
        for k in 0..<r { let b=try factors.matrix(k).congruence(c)
            for i in 0..<o { for a in 0..<v { ia[i*v+a,k]=b[i,o+a] } }
        }
        var correlation=0.0,minimum:Double?
        for i in 0..<o { for j in 0..<o {
            let bi=try VivoQMMatrix(rows:v,columns:r,values:Array(ia.values[(i*v*r)..<((i+1)*v*r)]))
            let bj=try VivoQMMatrix(rows:v,columns:r,values:Array(ia.values[(j*v*r)..<((j+1)*v*r)]))
            let pair=try bi.multiplied(by:bj.transposed)
            for a in 0..<v { for b in 0..<v {
                let denominator=hf.alphaOrbitalEnergies[i]+hf.alphaOrbitalEnergies[j]-hf.alphaOrbitalEnergies[o+a]-hf.alphaOrbitalEnergies[o+b]
                guard denominator < -minimumGapHartree else { throw VivoChemistryError.convergence("DF-MP2 small/inverted gap") }
                minimum=min(minimum ?? .infinity,abs(denominator));correlation+=pair[a,b]*(2*pair[a,b]-pair[b,a])/denominator
            } }
        } }
        guard correlation.isFinite else { throw VivoChemistryError.convergence("DF-MP2 overflow") }
        return .init(referenceEnergyHartree:energy,correlationEnergyHartree:correlation,minimumDenominatorMagnitude:minimum)
    }
}
