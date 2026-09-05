import Foundation

public enum VivoCIMethod: String, Codable, Sendable { case fci, cisd }
public struct VivoCIState: Codable, Sendable, Equatable {
    public let orbitalCount: Int
    public let alphaElectrons: Int
    public let betaElectrons: Int
    public let determinants: [UInt64]
    public let coefficients: [Double]
    public init(orbitalCount: Int, alphaElectrons: Int, betaElectrons: Int,
                determinants: [UInt64], coefficients: [Double]) {
        self.orbitalCount=orbitalCount;self.alphaElectrons=alphaElectrons;self.betaElectrons=betaElectrons
        self.determinants=determinants;self.coefficients=coefficients
    }
    public func validate(budget: VivoChemistryBudget = .init()) throws {
        try budget.validate()
        guard (1...31).contains(orbitalCount),orbitalCount<=budget.maximumBasisFunctions,
              (0...orbitalCount).contains(alphaElectrons),(0...orbitalCount).contains(betaElectrons),
              !determinants.isEmpty,determinants.count<=budget.maximumDeterminants,
              Set(determinants).count==determinants.count,coefficients.count==determinants.count,
              coefficients.allSatisfy(\.isFinite),abs(coefficients.reduce(0){$0+$1*$1}-1)<1e-9 else {
            throw VivoChemistryError.invalid("CI wavefunction shape, norm or sector")
        }
        let mask=(UInt64(1) << (2*orbitalCount))-1
        for d in determinants {
            guard d & ~mask == 0 else { throw VivoChemistryError.invalid("CI bit outside orbital range") }
            var na=0,nb=0
            for p in 0..<orbitalCount { na += Int((d >> (2*p)) & 1);nb += Int((d >> (2*p+1)) & 1) }
            guard na==alphaElectrons,nb==betaElectrons else { throw VivoChemistryError.invalid("CI determinant has wrong spin populations") }
        }
    }
}
public struct VivoCIResult: Codable, Sendable, Equatable {
    public let method: VivoCIMethod
    public let energyHartree: Double
    public let state: VivoCIState
    public let eigenResidualNorm: Double
    public let nextStateGapHartree: Double?
}
struct VivoFermionAction: Sendable { let mode: Int; let creation: Bool }
@inline(__always) func vivoApplyFermions(_ input: UInt64, _ actions: [VivoFermionAction]) -> (UInt64,Double)? {
    var determinant=input,sign=1.0
    for action in actions {
        let bit=UInt64(1)<<action.mode,occupied=determinant & bit != 0
        if occupied==action.creation { return nil }
        if (determinant & (bit-1)).nonzeroBitCount % 2 != 0 { sign = -sign }
        if action.creation { determinant |= bit } else { determinant &= ~bit }
    }
    return (determinant,sign)
}
private struct VivoFermionTerm { let coefficient: Double; let actions: [VivoFermionAction] }
public enum VivoConfigurationInteraction {
    private static func combinationCount(_ n: Int, _ k: Int, limit: Int) throws -> Int {
        let k=min(k,n-k); if k==0 { return 1 }; var value=1
        for i in 1...k {
            let (v,overflow)=value.multipliedReportingOverflow(by:n-k+i)
            guard !overflow else { throw VivoChemistryError.resourceLimit("determinant count overflow") }
            value=v/i
            guard value<=limit else { throw VivoChemistryError.resourceLimit("determinant-space budget; use a scalable solver") }
        }
        return value
    }
    private static func combinations(_ n: Int, _ count: Int, spin: Int) -> [UInt64] {
        var out:[UInt64]=[]
        func visit(_ first: Int,_ remaining: Int,_ bits: UInt64) {
            if remaining==0 { out.append(bits);return }
            if first>n-remaining { return }
            for p in first...(n-remaining) { visit(p+1,remaining-1,bits | (UInt64(1) << (2*p+spin))) }
        }
        visit(0,count,0);return out
    }
    public static func solve(_ h: VivoEmbeddedHamiltonian, method: VivoCIMethod = .fci,
                             budget: VivoChemistryBudget = .init()) throws -> VivoCIResult {
        try h.validate(budget:budget)
        let n=h.orbitalCount,na=h.alphaElectrons,nb=h.betaElectrons
        guard n<=31 else { throw VivoChemistryError.resourceLimit("CI determinant representation allows at most 31 spatial orbitals") }
        let ca=try combinationCount(n,na,limit:budget.maximumDeterminants)
        let cb=try combinationCount(n,nb,limit:budget.maximumDeterminants)
        let (fullCount,overflow)=ca.multipliedReportingOverflow(by:cb)
        guard !overflow,fullCount<=budget.maximumDeterminants else { throw VivoChemistryError.resourceLimit("fixed-spin sector exceeds dense CI determinant budget") }
        let alpha=combinations(n,na,spin:0),beta=combinations(n,nb,spin:1)
        var reference:UInt64=0
        for p in 0..<na { reference |= UInt64(1) << (2*p) }
        for p in 0..<nb { reference |= UInt64(1) << (2*p+1) }
        let determinants=alpha.flatMap { a in beta.map { a | $0 } }.filter {
            method == .fci || (($0 ^ reference).nonzeroBitCount / 2 <= 2)
        }.sorted()
        let d=determinants.count
        _ = try budget.elements([d,d],simultaneousArrays:5)
        let maximumTerms=2*n*n+4*n*n*n*n
        _ = try budget.elements([maximumTerms,16])
        let (upperWork,upperOverflow)=maximumTerms.multipliedReportingOverflow(by:d)
        guard !upperOverflow,upperWork<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("CI worst-case work budget") }
        var matrix=VivoQMMatrix(d,d),terms:[VivoFermionTerm]=[]
        for p in 0..<n { for q in 0..<n where h.oneElectron[p,q] != 0 {
            for spin in 0..<2 { terms.append(.init(coefficient:h.oneElectron[p,q],actions:[.init(mode:2*q+spin,creation:false),.init(mode:2*p+spin,creation:true)])) }
        } }
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            let value=0.5*h.eri(p,q,r,s); if value==0 { continue }
            for sigma in 0..<2 { for tau in 0..<2 {
                terms.append(.init(coefficient:value,actions:[.init(mode:2*q+sigma,creation:false),.init(mode:2*s+tau,creation:false),
                                                              .init(mode:2*r+tau,creation:true),.init(mode:2*p+sigma,creation:true)]))
            } }
        } } } }
        let (work,workOverflow)=terms.count.multipliedReportingOverflow(by:d)
        guard !workOverflow,work<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("CI Hamiltonian operator-application budget") }
        let index=Dictionary(uniqueKeysWithValues:determinants.enumerated().map{($0.element,$0.offset)})
        for (j,determinant) in determinants.enumerated() {
            matrix[j,j]=h.constantEnergyHartree
            for term in terms {
                if let (destination,sign)=vivoApplyFermions(determinant,term.actions),let i=index[destination] { matrix[i,j] += sign*term.coefficient }
            }
        }
        for i in 0..<d { for j in 0..<i {
            guard abs(matrix[i,j]-matrix[j,i])<1e-9 else { throw VivoChemistryError.invalid("CI Hamiltonian lost Hermiticity") }
            let mean=0.5*(matrix[i,j]+matrix[j,i]);matrix[i,j]=mean;matrix[j,i]=mean
        } }
        let eig=try VivoQMDenseAlgebra.symmetricEigen(matrix,tolerance:1e-13,maximumSweeps:128)
        var coefficients=(0..<d).map{eig.vectors[$0,0]}
        if let largest=coefficients.indices.max(by:{abs(coefficients[$0])<abs(coefficients[$1])}),coefficients[largest]<0 { coefficients=coefficients.map{-$0} }
        var residualSquared=0.0
        for i in 0..<d {
            let value=(0..<d).reduce(0.0){$0+matrix[i,$1]*coefficients[$1]}-eig.values[0]*coefficients[i]
            residualSquared += value*value
        }
        let residual=sqrt(residualSquared)
        guard residual.isFinite,residual<1e-8 else { throw VivoChemistryError.convergence("CI eigenpair residual") }
        let state=VivoCIState(orbitalCount:n,alphaElectrons:na,betaElectrons:nb,determinants:determinants,coefficients:coefficients)
        try state.validate(budget:budget)
        return .init(method:method,energyHartree:eig.values[0],state:state,eigenResidualNorm:residual,
                     nextStateGapHartree:d>1 ? eig.values[1]-eig.values[0] : nil)
    }
}

public struct VivoSpinRDMs: Codable, Sendable, Equatable {
    public let spinOrbitalCount: Int
    public let electronCount: Int
    public let one: VivoQMMatrix
    public let two: [Double]
    public func gamma2(_ p: Int,_ q: Int,_ r: Int,_ s: Int) -> Double { let n=spinOrbitalCount;return two[((p*n+q)*n+r)*n+s] }
    public func cumulant(_ p: Int,_ q: Int,_ r: Int,_ s: Int) -> Double { gamma2(p,q,r,s)-one[p,r]*one[q,s]+one[p,s]*one[q,r] }
    public var spatialOne: VivoQMMatrix {
        let n=spinOrbitalCount/2;var out=VivoQMMatrix(n,n)
        for p in 0..<n { for q in 0..<n { out[p,q]=one[2*p,2*q]+one[2*p+1,2*q+1] } };return out
    }
    public func energy(of h: VivoEmbeddedHamiltonian) throws -> Double {
        try h.validate()
        guard spinOrbitalCount==2*h.orbitalCount,electronCount==h.alphaElectrons+h.betaElectrons else { throw VivoChemistryError.invalid("RDM/Hamiltonian mismatch") }
        let n=h.orbitalCount;var value=h.constantEnergyHartree
        for p in 0..<n { for q in 0..<n {
            for spin in 0..<2 { value += h.oneElectron[p,q]*one[2*p+spin,2*q+spin] }
            for r in 0..<n { for s in 0..<n { for sigma in 0..<2 { for tau in 0..<2 {
                value += 0.5*h.eri(p,q,r,s)*gamma2(2*p+sigma,2*r+tau,2*q+sigma,2*s+tau)
            } } } }
        } }
        return value
    }
}
public enum VivoCIDensityMatrices {
    static func expectation(_ state: VivoCIState, _ index: [UInt64:Int], _ actions: [VivoFermionAction]) -> Double {
        var value=0.0
        for (j,d) in state.determinants.enumerated() {
            if let (destination,sign)=vivoApplyFermions(d,actions),let i=index[destination] { value += state.coefficients[i]*state.coefficients[j]*sign }
        }
        return value
    }
    public static func compute(_ state: VivoCIState, budget: VivoChemistryBudget = .init()) throws -> VivoSpinRDMs {
        try state.validate(budget:budget)
        let m=2*state.orbitalCount,count=try budget.elements([m,m,m,m],simultaneousArrays:2)
        let (work,overflow)=count.multipliedReportingOverflow(by:state.determinants.count)
        guard !overflow,work<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("RDM operator-application budget") }
        var one=VivoQMMatrix(m,m),two=[Double](repeating:0,count:count)
        let index=Dictionary(uniqueKeysWithValues:state.determinants.enumerated().map{($0.element,$0.offset)})
        for p in 0..<m { for q in 0..<m {
            one[p,q]=expectation(state,index,[.init(mode:q,creation:false),.init(mode:p,creation:true)])
            for r in 0..<m { for s in 0..<m where p != q && r != s {
                two[((p*m+q)*m+r)*m+s]=expectation(state,index,[.init(mode:r,creation:false),.init(mode:s,creation:false),
                                                                              .init(mode:q,creation:true),.init(mode:p,creation:true)])
            } }
        } }
        let electrons=state.alphaElectrons+state.betaElectrons
        guard abs((0..<m).reduce(0.0){$0+one[$1,$1]}-Double(electrons))<1e-8 else { throw VivoChemistryError.invalid("1-RDM trace") }
        let result=VivoSpinRDMs(spinOrbitalCount:m,electronCount:electrons,one:one,two:two)
        for p in 0..<m { for q in 0..<m {
            let contracted=(0..<m).reduce(0.0){$0+result.gamma2(p,$1,q,$1)}
            guard abs(contracted-Double(electrons-1)*one[p,q])<1e-8 else { throw VivoChemistryError.invalid("2-RDM contraction") }
        } }
        return result
    }
    public static func orbitalMarginal(_ state: VivoCIState, orbitals: [Int], budget: VivoChemistryBudget = .init()) throws -> VivoQMMatrix {
        try state.validate(budget:budget)
        guard (1...2).contains(orbitals.count),Set(orbitals).count==orbitals.count,
              orbitals.allSatisfy({$0>=0 && $0<state.orbitalCount}) else { throw VivoChemistryError.invalid("one/two-orbital marginal selection") }
        let modes=orbitals.sorted().flatMap{[2*$0,2*$0+1]},dimension=1<<modes.count
        let selectedMask=modes.reduce(UInt64(0)){$0 | (UInt64(1)<<$1)}
        var groups:[UInt64:[Double]]=[:]
        for (i,determinant) in state.determinants.enumerated() {
            let environment=determinant & ~selectedMask
            var local=0,parity=0
            for (slot,mode) in modes.enumerated() where determinant & (UInt64(1)<<mode) != 0 {
                local |= 1<<slot;parity += (environment & ((UInt64(1)<<mode)-1)).nonzeroBitCount
            }
            var vector=groups[environment] ?? [Double](repeating:0,count:dimension)
            vector[local] += (parity%2==0 ? 1 : -1)*state.coefficients[i];groups[environment]=vector
        }
        var marginal=VivoQMMatrix(dimension,dimension)
        for key in groups.keys.sorted() {
            let vector=groups[key]!
            for i in 0..<dimension { for j in 0..<dimension { marginal[i,j] += vector[i]*vector[j] } }
        }
        return marginal
    }
    public static func entropy(_ density: VivoQMMatrix) throws -> Double {
        let eigen=try VivoQMDenseAlgebra.symmetricEigen(density)
        guard abs(eigen.values.reduce(0,+)-1)<1e-9,eigen.values.allSatisfy({$0 >= -1e-10 && $0<=1+1e-10}) else {
            throw VivoChemistryError.invalid("orbital marginal positivity/normalization")
        }
        return -eigen.values.filter{$0>0}.reduce(0){$0+$1*log($1)}
    }
    public static func orbitalInformation(_ state: VivoCIState, budget: VivoChemistryBudget = .init()) throws -> (entropy: [Double], mutualInformation: VivoQMMatrix) {
        try state.validate(budget:budget)
        let n=state.orbitalCount
        let entropies=try (0..<n).map { try entropy(orbitalMarginal(state,orbitals:[$0],budget:budget)) }
        var mutual=VivoQMMatrix(n,n)
        for i in 0..<n { for j in 0..<i {
            let sij=try entropy(orbitalMarginal(state,orbitals:[i,j],budget:budget)),value=entropies[i]+entropies[j]-sij
            guard value >= -1e-9 else { throw VivoChemistryError.invalid("negative mutual information") }
            mutual[i,j]=max(0,value);mutual[j,i]=mutual[i,j]
        } }
        return (entropies,mutual)
    }
}
