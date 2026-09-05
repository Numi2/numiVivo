import Foundation

/// Portable classical/quantum boundary. Real, orthonormal spatial orbitals;
/// chemist ERIs, interleaved spin order (p-alpha=2p, p-beta=2p+1).
/// H = c + sum h[pq] a†pσ aqσ + 1/2 sum (pq|rs) a†pσ a†rτ asτ aqσ.
public struct VivoEmbeddedHamiltonian: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/embedded-hamiltonian/v1"
    public let schema: String
    public let orbitalIdentifiers: [String]
    public let alphaElectrons: Int
    public let betaElectrons: Int
    public let oneElectron: VivoQMMatrix
    public let twoElectron: [Double]
    public let constantEnergyHartree: Double
    public let energyReference: String
    public let provenance: [String:String]
    public var orbitalCount: Int { orbitalIdentifiers.count }
    public init(orbitalIdentifiers: [String], alphaElectrons: Int, betaElectrons: Int,
                oneElectron: VivoQMMatrix, twoElectron: [Double], constantEnergyHartree: Double,
                energyReference: String, provenance: [String:String] = [:]) {
        schema=Self.schema;self.orbitalIdentifiers=orbitalIdentifiers;self.alphaElectrons=alphaElectrons
        self.betaElectrons=betaElectrons;self.oneElectron=oneElectron;self.twoElectron=twoElectron
        self.constantEnergyHartree=constantEnergyHartree;self.energyReference=energyReference;self.provenance=provenance
    }
    public func eri(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        let n=orbitalCount;return twoElectron[((p*n+q)*n+r)*n+s]
    }
    public func validate(budget: VivoChemistryBudget = .init()) throws {
        let n=orbitalCount
        guard schema==Self.schema,n>0,n<=budget.maximumBasisFunctions,
              Set(orbitalIdentifiers).count==n,orbitalIdentifiers.allSatisfy({!$0.isEmpty}),
              alphaElectrons>=0,betaElectrons>=0,alphaElectrons<=n,betaElectrons<=n,
              oneElectron.rows==n,oneElectron.columns==n,oneElectron.values.allSatisfy(\.isFinite),
              constantEnergyHartree.isFinite,!energyReference.isEmpty else { throw VivoChemistryError.invalid("embedded Hamiltonian identity, spin or dimensions") }
        guard twoElectron.count == (try budget.elements([n,n,n,n],simultaneousArrays:2)),
              twoElectron.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("embedded ERI shape/values") }
        for p in 0..<n { for q in 0..<n {
            guard abs(oneElectron[p,q]-oneElectron[q,p]) < 1e-10 else { throw VivoChemistryError.invalid("one-electron Hermiticity") }
            for r in 0..<n { for s in 0..<n {
                let v=eri(p,q,r,s),tolerance=1e-10*max(1,abs(v))
                guard abs(v-eri(q,p,r,s))<=tolerance,abs(v-eri(p,q,s,r))<=tolerance,
                      abs(v-eri(r,s,p,q))<=tolerance else { throw VivoChemistryError.invalid("embedded ERI convention/symmetry") }
            } }
        } }
    }
    public static func fromAO(_ ao: VivoAOIntegrals, coefficients c: VivoQMMatrix,
                              alphaElectrons: Int, betaElectrons: Int, orbitalIdentifiers: [String],
                              energyReference: String, budget: VivoChemistryBudget = .init()) throws -> Self {
        try ao.validate(budget:budget)
        guard c.rows==ao.count,c.columns==orbitalIdentifiers.count,c.columns>0 else { throw VivoChemistryError.invalid("AO/MO coefficient dimensions") }
        let metric=try ao.overlap.congruence(c),identity=try VivoQMMatrix.identity(c.columns)
        guard try metric.adding(identity,scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("MO columns are not AO-metric orthonormal") }
        let result=Self(orbitalIdentifiers:orbitalIdentifiers,alphaElectrons:alphaElectrons,betaElectrons:betaElectrons,
                        oneElectron:try ao.coreHamiltonian.congruence(c),
                        twoElectron:try VivoOrbitalTransform.transformERI(ao.electronRepulsion,inputCount:ao.count,coefficients:c,budget:budget),
                        constantEnergyHartree:ao.constantEnergyHartree,energyReference:energyReference,
                        provenance:["integralConvention":"real-spatial-chemist","precision":"fp64"])
        try result.validate(budget:budget); return result
    }
    public func frozenCore(active: [Int], doublyOccupiedCore: [Int],
                           budget: VivoChemistryBudget = .init()) throws -> Self {
        try validate(budget:budget)
        let n=orbitalCount,m=active.count,core=doublyOccupiedCore
        guard m>0,Set(active+core).count==m+core.count,(active+core).allSatisfy({$0>=0 && $0<n}),
              core.count<=alphaElectrons,core.count<=betaElectrons,
              alphaElectrons-core.count<=m,betaElectrons-core.count<=m else { throw VivoChemistryError.invalid("frozen core/active partition") }
        var h=VivoQMMatrix(m,m),g=[Double](repeating:0,count:try budget.elements([m,m,m,m],simultaneousArrays:2))
        var constant=constantEnergyHartree
        for i in core { constant += 2*oneElectron[i,i]
            for j in core { constant += 2*eri(i,i,j,j)-eri(i,j,j,i) }
        }
        for p in 0..<m { for q in 0..<m {
            let a=active[p],b=active[q];h[p,q]=oneElectron[a,b]
            for i in core { h[p,q] += 2*eri(a,b,i,i)-eri(a,i,i,b) }
            for r in 0..<m { for s in 0..<m { g[((p*m+q)*m+r)*m+s]=eri(a,b,active[r],active[s]) } }
        } }
        var sources=provenance
        sources["frozenDoublyOccupiedOrbitals"]=core.map{orbitalIdentifiers[$0]}.joined(separator:",")
        let result=Self(orbitalIdentifiers:active.map{orbitalIdentifiers[$0]},alphaElectrons:alphaElectrons-core.count,
                        betaElectrons:betaElectrons-core.count,oneElectron:h,twoElectron:g,constantEnergyHartree:constant,
                        energyReference:energyReference,provenance:sources)
        try result.validate(budget:budget);return result
    }
    public func rotated(by u: VivoQMMatrix, budget: VivoChemistryBudget = .init()) throws -> Self {
        try validate(budget:budget)
        guard u.rows==orbitalCount,u.columns==orbitalCount else { throw VivoChemistryError.invalid("orbital rotation dimension") }
        let metric=try u.transposed.multiplied(by:u)
        guard try metric.adding(.identity(orbitalCount),scale:-1).frobeniusNorm<1e-8 else { throw VivoChemistryError.invalid("nonorthogonal orbital rotation") }
        return .init(orbitalIdentifiers:orbitalIdentifiers,alphaElectrons:alphaElectrons,betaElectrons:betaElectrons,
                     oneElectron:try oneElectron.congruence(u),twoElectron:try VivoOrbitalTransform.transformERI(twoElectron,inputCount:orbitalCount,coefficients:u,budget:budget),
                     constantEnergyHartree:constantEnergyHartree,energyReference:energyReference,provenance:provenance)
    }
}
public enum VivoOrbitalTransform {
    public static func transformERI(_ input: [Double], inputCount n: Int, coefficients c: VivoQMMatrix,
                                    budget: VivoChemistryBudget = .init()) throws -> [Double] {
        let m=c.columns
        guard n>0,m>0,n<=budget.maximumBasisFunctions,m<=budget.maximumBasisFunctions,c.rows==n,
              c.values.allSatisfy(\.isFinite),input.allSatisfy(\.isFinite),
              input.count==(try budget.elements([n,n,n,n],simultaneousArrays:6)) else { throw VivoChemistryError.invalid("ERI transformation shape") }
        var current=input,dims=[n,n,n,n]
        for axis in 0..<4 {
            var nextDims=dims;nextDims[axis]=m
            var output=[Double](repeating:0,count:try budget.elements(nextDims,simultaneousArrays:6))
            let inner=dims[(axis+1)...].reduce(1,*),outer=dims[..<axis].reduce(1,*)
            _ = try budget.elements([outer,inner,n],simultaneousArrays:6)
            var packed=VivoQMMatrix(outer*inner,n)
            for before in 0..<outer { for after in 0..<inner { for k in 0..<n {
                packed[before*inner+after,k]=current[(before*n+k)*inner+after]
            } } }
            let transformed=try packed.multiplied(by:c)
            for before in 0..<outer { for after in 0..<inner { for k in 0..<m {
                output[(before*m+k)*inner+after]=transformed[before*inner+after,k]
            } } }
            current=output;dims=nextDims
        }
        return current
    }
}
public struct VivoMP2Result: Codable, Sendable, Equatable {
    public let referenceEnergyHartree: Double
    public let correlationEnergyHartree: Double
    public var totalEnergyHartree: Double { referenceEnergyHartree+correlationEnergyHartree }
    public let minimumDenominatorMagnitude: Double?
}
public enum VivoRestrictedMP2 {
    public static func solve(_ h: VivoEmbeddedHamiltonian, minimumGapHartree: Double = 1e-8,
                             budget: VivoChemistryBudget = .init()) throws -> VivoMP2Result {
        try h.validate(budget:budget)
        guard h.alphaElectrons==h.betaElectrons,minimumGapHartree.isFinite,minimumGapHartree>0 else {
            throw VivoChemistryError.invalid("restricted MP2 requires a closed-shell reference and positive gap guard")
        }
        let n=h.orbitalCount,occupied=h.alphaElectrons
        var f=h.oneElectron,reference=h.constantEnergyHartree
        for p in 0..<n { for q in 0..<n { for i in 0..<occupied { f[p,q] += 2*h.eri(p,q,i,i)-h.eri(p,i,i,q) } } }
        for p in 0..<n { for q in 0..<p where abs(f[p,q])>1e-7 {
            throw VivoChemistryError.invalid("MP2 input orbitals must be canonical; rotate occupied/virtual blocks first")
        } }
        for i in 0..<occupied { reference += 2*h.oneElectron[i,i]
            for j in 0..<occupied { reference += 2*h.eri(i,i,j,j)-h.eri(i,j,j,i) }
        }
        var correlation=0.0,minimum:Double?
        for i in 0..<occupied { for j in 0..<occupied { for a in occupied..<n { for b in occupied..<n {
            let denominator=f[i,i]+f[j,j]-f[a,a]-f[b,b]
            guard denominator.isFinite,denominator < -minimumGapHartree else { throw VivoChemistryError.convergence("MP2 non-Aufbau or near-zero denominator; no silent regularization") }
            minimum=min(minimum ?? .infinity,abs(denominator))
            let g=h.eri(i,a,j,b)
            correlation += g*(2*g-h.eri(i,b,j,a))/denominator
        } } } }
        guard reference.isFinite,correlation.isFinite else { throw VivoChemistryError.convergence("MP2 nonfinite energy") }
        return .init(referenceEnergyHartree:reference,correlationEnergyHartree:correlation,minimumDenominatorMagnitude:minimum)
    }
}
