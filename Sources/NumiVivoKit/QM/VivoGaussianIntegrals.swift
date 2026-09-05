import Foundation

public struct VivoCartesianOrbital: Codable, Sendable, Equatable {
    public let nucleusIndex: Int
    public let angular: [Int]
    public let shellIndex: Int
    public let primitiveExponents: [Double]
    /// Includes primitive and contracted-function normalization.
    public let weights: [Double]
}

/// Chemist ordering: eri[((p*n+q)*n+r)*n+s] = (pq|rs), in Hartree.
/// Constant energy contains nucleus/nucleus and nucleus/MM-charge interactions,
/// but never the frozen MM self-energy or a classical QM/QM electrostatic term.
public struct VivoAOIntegrals: Codable, Sendable, Equatable {
    public let sourceSystem: VivoElectronicSystem
    public let sourceBasis: VivoGaussianBasis
    public let orbitals: [VivoCartesianOrbital]
    public let overlap: VivoQMMatrix
    public let kinetic: VivoQMMatrix
    public let attraction: VivoQMMatrix
    public let electronRepulsion: [Double]
    public let constantEnergyHartree: Double
    public var count: Int { overlap.rows }
    public var coreHamiltonian: VivoQMMatrix {
        var m = kinetic
        for i in m.values.indices { m.values[i] += attraction.values[i] }; return m
    }
    public func eri(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        electronRepulsion[((p*count+q)*count+r)*count+s]
    }
    public func validate(budget: VivoChemistryBudget = .init()) throws {
        let n = count
        guard try VivoGaussianIntegralEngine.expanded(system:sourceSystem,basis:sourceBasis,budget:budget)==orbitals else { throw VivoChemistryError.invalid("AO orbital/source binding") }
        guard n > 0, n <= budget.maximumBasisFunctions, orbitals.count == n,
              overlap.columns == n, kinetic.rows == n, kinetic.columns == n,
              attraction.rows == n, attraction.columns == n,
              constantEnergyHartree.isFinite else { throw VivoChemistryError.invalid("AO integral dimensions") }
        let size = try budget.elements([n,n,n,n], simultaneousArrays: 2)
        guard electronRepulsion.count == size, electronRepulsion.allSatisfy(\.isFinite),
              overlap.values.allSatisfy(\.isFinite), kinetic.values.allSatisfy(\.isFinite),
              attraction.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("AO integral values") }
        for p in 0..<n { for q in 0..<n {
            guard abs(overlap[p,q]-overlap[q,p]) < 1e-10,
                  abs(kinetic[p,q]-kinetic[q,p]) < 1e-10,
                  abs(attraction[p,q]-attraction[q,p]) < 1e-10 else {
                throw VivoChemistryError.invalid("one-electron integral symmetry")
            }
            for r in 0..<n { for s in 0..<n {
                let v = eri(p,q,r,s), tolerance = 1e-10 * max(1,abs(v))
                guard abs(v-eri(q,p,r,s)) <= tolerance, abs(v-eri(p,q,s,r)) <= tolerance,
                      abs(v-eri(r,s,p,q)) <= tolerance else { throw VivoChemistryError.invalid("ERI permutation symmetry") }
            } }
        } }
    }
}

private struct HermiteIndex: Hashable { let i: Int; let j: Int; let t: Int }
private struct CoulombIndex: Hashable { let t: Int; let u: Int; let v: Int; let n: Int }

/// McMurchie-Davidson Cartesian Gaussian integrals. This is a bounded, dense,
/// all-electron FP64 implementation; it is not RI, ECP or a spherical-basis engine.
public enum VivoGaussianIntegralEngine {
    /// Stable positive-term expansion at small/moderate argument, upward
    /// recurrence only where its cancellation is well-conditioned.
    public static func boys(order n: Int, argument t: Double) throws -> Double {
        guard (0...32).contains(n), t.isFinite, t >= 0 else { throw VivoChemistryError.invalid("Boys function arguments") }
        return boysUnchecked(n,t)
    }
    private static func boysUnchecked(_ n: Int, _ t: Double) -> Double {
        if t < 60 {
            var term = 1 / Double(2*n+1), sum = term
            for k in 1...512 {
                term *= 2*t / Double(2*n+2*k+1); sum += term
                if abs(term) < abs(sum)*2e-16 { break }
            }
            return exp(-t)*sum
        }
        let exponential = exp(-t)
        var f = 0.5 * sqrt(Double.pi/t) * erf(sqrt(t))
        if n > 0 { for m in 0..<n { f = (Double(2*m+1)*f-exponential)/(2*t) } }
        return f
    }
    private static func hermite(_ i: Int, _ j: Int, _ t: Int, _ q: Double,
                                _ a: Double, _ b: Double, _ cache: inout [HermiteIndex:Double]) -> Double {
        if i < 0 || j < 0 || t < 0 || t > i+j { return 0 }
        let key = HermiteIndex(i:i,j:j,t:t)
        if let value = cache[key] { return value }
        let p = a+b
        let value: Double
        if i == 0 && j == 0 { value = exp(-a*b/p*q*q) }
        else if j == 0 {
            value = hermite(i-1,j,t-1,q,a,b,&cache)/(2*p)
                - b*q/p*hermite(i-1,j,t,q,a,b,&cache)
                + Double(t+1)*hermite(i-1,j,t+1,q,a,b,&cache)
        } else {
            value = hermite(i,j-1,t-1,q,a,b,&cache)/(2*p)
                + a*q/p*hermite(i,j-1,t,q,a,b,&cache)
                + Double(t+1)*hermite(i,j-1,t+1,q,a,b,&cache)
        }
        cache[key] = value; return value
    }
    private static func coefficients(_ i: Int, _ j: Int, _ q: Double, _ a: Double, _ b: Double) -> [Double] {
        var cache: [HermiteIndex:Double] = [:]
        return (0...(i+j)).map { hermite(i,j,$0,q,a,b,&cache) }
    }
    private static func coulomb(_ t: Int, _ u: Int, _ v: Int, _ n: Int,
                                _ p: Double, _ pc: SIMD3<Double>, _ cache: inout [CoulombIndex:Double]) -> Double {
        if t < 0 || u < 0 || v < 0 { return 0 }
        let key = CoulombIndex(t:t,u:u,v:v,n:n)
        if let value = cache[key] { return value }
        let value: Double
        if t == 0 && u == 0 && v == 0 {
            value = pow(-2*p,Double(n)) * boysUnchecked(n,p*vivoQMDot(pc,pc))
        } else if t > 0 {
            value = (t > 1 ? Double(t-1)*coulomb(t-2,u,v,n+1,p,pc,&cache) : 0)
                + pc.x*coulomb(t-1,u,v,n+1,p,pc,&cache)
        } else if u > 0 {
            value = (u > 1 ? Double(u-1)*coulomb(t,u-2,v,n+1,p,pc,&cache) : 0)
                + pc.y*coulomb(t,u-1,v,n+1,p,pc,&cache)
        } else {
            value = (v > 1 ? Double(v-1)*coulomb(t,u,v-2,n+1,p,pc,&cache) : 0)
                + pc.z*coulomb(t,u,v-1,n+1,p,pc,&cache)
        }
        cache[key] = value; return value
    }
    private static func primitiveOverlap(_ a: Double, _ la: [Int], _ ra: SIMD3<Double>,
                                         _ b: Double, _ lb: [Int], _ rb: SIMD3<Double>) -> Double {
        if la.contains(where: {$0 < 0}) || lb.contains(where: {$0 < 0}) { return 0 }
        var result = pow(Double.pi/(a+b),1.5)
        for axis in 0..<3 { result *= coefficients(la[axis],lb[axis],ra[axis]-rb[axis],a,b)[0] }
        return result
    }
    private static func primitiveKinetic(_ a: Double, _ la: [Int], _ ra: SIMD3<Double>,
                                         _ b: Double, _ lb: [Int], _ rb: SIMD3<Double>) -> Double {
        var result = b*Double(2*lb.reduce(0,+)+3)*primitiveOverlap(a,la,ra,b,lb,rb)
        for axis in 0..<3 {
            var raised = lb; raised[axis] += 2
            result -= 2*b*b*primitiveOverlap(a,la,ra,b,raised,rb)
            if lb[axis] >= 2 {
                var lowered = lb; lowered[axis] -= 2
                result -= 0.5*Double(lb[axis]*(lb[axis]-1))*primitiveOverlap(a,la,ra,b,lowered,rb)
            }
        }
        return result
    }
    private static func primitivePotential(_ a: Double, _ la: [Int], _ ra: SIMD3<Double>,
                                           _ b: Double, _ lb: [Int], _ rb: SIMD3<Double>,
                                           _ rc: SIMD3<Double>) -> Double {
        let p = a+b, pc = (a*ra+b*rb)/p-rc
        let e = (0..<3).map { coefficients(la[$0],lb[$0],ra[$0]-rb[$0],a,b) }
        var cache: [CoulombIndex:Double] = [:], value = 0.0
        for t in e[0].indices { for u in e[1].indices { for v in e[2].indices {
            value += e[0][t]*e[1][u]*e[2][v]*coulomb(t,u,v,0,p,pc,&cache)
        } } }
        return 2*Double.pi/p*value
    }
    private static func primitiveERI(_ a: Double, _ la: [Int], _ ra: SIMD3<Double>,
                                     _ b: Double, _ lb: [Int], _ rb: SIMD3<Double>,
                                     _ c: Double, _ lc: [Int], _ rc: SIMD3<Double>,
                                     _ d: Double, _ ld: [Int], _ rd: SIMD3<Double>) -> Double {
        let p = a+b, q = c+d, alpha = p*q/(p+q), pq = (a*ra+b*rb)/p-(c*rc+d*rd)/q
        let e = (0..<3).map { coefficients(la[$0],lb[$0],ra[$0]-rb[$0],a,b) }
        let f = (0..<3).map { coefficients(lc[$0],ld[$0],rc[$0]-rd[$0],c,d) }
        var cache: [CoulombIndex:Double] = [:], value = 0.0
        for t in e[0].indices { for u in e[1].indices { for v in e[2].indices {
            let ab = e[0][t]*e[1][u]*e[2][v]
            if ab == 0 { continue }
            for x in f[0].indices { for y in f[1].indices { for z in f[2].indices {
                let cd = f[0][x]*f[1][y]*f[2][z]
                if cd == 0 { continue }
                value += ab*cd*((x+y+z)%2 == 0 ? 1 : -1)*coulomb(t+x,u+y,v+z,0,alpha,pq,&cache)
            } } }
        } } }
        return 2*pow(Double.pi,2.5)/(p*q*sqrt(p+q))*value
    }
    private static func oddDoubleFactorial(_ n: Int) -> Double {
        if n <= 0 { return 1 }; return stride(from:n,through:1,by:-2).reduce(1) { $0*Double($1) }
    }
    private static func normalization(_ exponent: Double, _ angular: [Int]) -> Double {
        let l = angular.reduce(0,+), denominator = angular.reduce(1.0) { $0*oddDoubleFactorial(2*$1-1) }
        return pow(2*exponent/Double.pi,0.75)*sqrt(pow(4*exponent,Double(l))/denominator)
    }
    static func expanded(system: VivoElectronicSystem, basis: VivoGaussianBasis,
                         budget: VivoChemistryBudget) throws -> [VivoCartesianOrbital] {
        try budget.validate()
        guard system.nuclei.count<=budget.maximumBasisFunctions else { throw VivoChemistryError.resourceLimit("QM nucleus budget") }
        try system.validate(); try basis.validate(nucleusCount: system.nuclei.count)
        var result: [VivoCartesianOrbital] = []
        for (shellIndex,shell) in basis.shells.enumerated() {
            let l = shell.angularMomentum
            for x in stride(from:l,through:0,by:-1) { for y in stride(from:l-x,through:0,by:-1) {
                guard result.count < budget.maximumBasisFunctions else { throw VivoChemistryError.resourceLimit("basis-function budget") }
                let angular = [x,y,l-x-y], r = system.nuclei[shell.nucleusIndex].positionBohr
                let exponents = shell.primitives.map(\.exponent)
                var weights = shell.primitives.map { $0.coefficient*normalization($0.exponent,angular) }
                var norm = 0.0
                for i in exponents.indices { for j in exponents.indices {
                    norm += weights[i]*weights[j]*primitiveOverlap(exponents[i],angular,r,exponents[j],angular,r)
                } }
                guard norm.isFinite, norm > 1e-24 else { throw VivoChemistryError.invalid("singular contracted Gaussian") }
                for i in weights.indices { weights[i] /= sqrt(norm) }
                result.append(.init(nucleusIndex:shell.nucleusIndex,angular:angular,shellIndex:shellIndex,
                                    primitiveExponents:exponents,weights:weights))
            } }
        }
        return result
    }
    public static func compute(system: VivoElectronicSystem, basis: VivoGaussianBasis,
                               budget: VivoChemistryBudget = .init()) throws -> VivoAOIntegrals {
        let orbitals = try expanded(system:system,basis:basis,budget:budget), n = orbitals.count
        let size = try budget.elements([n,n,n,n],simultaneousArrays:2)
        _ = try budget.elements([n,n],simultaneousArrays:20)
        var overlap = VivoQMMatrix(n,n), kinetic = overlap, attraction = overlap
        var eri = [Double](repeating:0,count:size)
        var operatorWork = 0
        let potentialCentres = system.nuclei.count + system.pointCharges.filter{$0.chargeE != 0}.count
        func center(_ orbital: VivoCartesianOrbital) -> SIMD3<Double> { system.nuclei[orbital.nucleusIndex].positionBohr }
        for p in 0..<n { for q in 0...p {
            let a = orbitals[p], b = orbitals[q], ra = center(a), rb = center(b)
            var s = 0.0, t = 0.0, v = 0.0
            let (pairWork, overflow) = (a.weights.count*b.weights.count).multipliedReportingOverflow(by:potentialCentres+2)
            guard !overflow,pairWork<=budget.maximumOperatorApplications-operatorWork else { throw VivoChemistryError.resourceLimit("one-electron/point-charge work budget") }
            operatorWork += pairWork
            for i in a.primitiveExponents.indices { for j in b.primitiveExponents.indices {
                let ea = a.primitiveExponents[i], eb = b.primitiveExponents[j], w = a.weights[i]*b.weights[j]
                s += w*primitiveOverlap(ea,a.angular,ra,eb,b.angular,rb)
                t += w*primitiveKinetic(ea,a.angular,ra,eb,b.angular,rb)
                for nucleus in system.nuclei {
                    v -= w*Double(nucleus.atomicNumber)*primitivePotential(ea,a.angular,ra,eb,b.angular,rb,nucleus.positionBohr)
                }
                for charge in system.pointCharges where charge.chargeE != 0 {
                    v -= w*charge.chargeE*primitivePotential(ea,a.angular,ra,eb,b.angular,rb,charge.positionBohr)
                }
            } }
            overlap[p,q] = s; overlap[q,p] = s; kinetic[p,q] = t; kinetic[q,p] = t
            attraction[p,q] = v; attraction[q,p] = v
        } }
        for p in 0..<n { for q in 0...p { for r in 0...p { for s in 0...r {
            if r*(r+1)/2+s > p*(p+1)/2+q { continue }
            let a=orbitals[p],b=orbitals[q],c=orbitals[r],d=orbitals[s]
            let ra=center(a),rb=center(b),rc=center(c),rd=center(d)
            var value = 0.0
            for i in a.weights.indices { for j in b.weights.indices { for k in c.weights.indices { for l in d.weights.indices {
                guard operatorWork < budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("primitive-quartet work budget") }
                operatorWork += 1
                let w = a.weights[i]*b.weights[j]*c.weights[k]*d.weights[l]
                if w != 0 {
                    value += w*primitiveERI(a.primitiveExponents[i],a.angular,ra,b.primitiveExponents[j],b.angular,rb,
                                            c.primitiveExponents[k],c.angular,rc,d.primitiveExponents[l],d.angular,rd)
                }
            } } } }
            for (i,j,k,l) in [(p,q,r,s),(q,p,r,s),(p,q,s,r),(q,p,s,r),(r,s,p,q),(s,r,p,q),(r,s,q,p),(s,r,q,p)] {
                eri[((i*n+j)*n+k)*n+l] = value
            }
        } } } }
        var constant = 0.0
        for i in system.nuclei.indices {
            let a = system.nuclei[i]
            for j in 0..<i {
                let b = system.nuclei[j]; constant += Double(a.atomicNumber*b.atomicNumber)/vivoQMNorm(a.positionBohr-b.positionBohr)
            }
            for q in system.pointCharges where q.chargeE != 0 {
                constant += Double(a.atomicNumber)*q.chargeE/vivoQMNorm(a.positionBohr-q.positionBohr)
            }
        }
        let result = VivoAOIntegrals(sourceSystem:system,sourceBasis:basis,orbitals:orbitals,overlap:overlap,kinetic:kinetic,attraction:attraction,
                                    electronRepulsion:eri,constantEnergyHartree:constant)
        try result.validate(budget:budget); return result
    }
    /// Cross-geometry overlap for path tracking. Orbital indices must originate
    /// from the corresponding basis expansion, never be inferred from energies.
    public static func crossOverlap(leftSystem: VivoElectronicSystem, leftBasis: VivoGaussianBasis,
                                    rightSystem: VivoElectronicSystem, rightBasis: VivoGaussianBasis,
                                    budget: VivoChemistryBudget = .init()) throws -> VivoQMMatrix {
        let a = try expanded(system:leftSystem,basis:leftBasis,budget:budget)
        let b = try expanded(system:rightSystem,basis:rightBasis,budget:budget)
        _ = try budget.elements([a.count,b.count]); var result = VivoQMMatrix(a.count,b.count)
        for p in a.indices { for q in b.indices {
            let ra=leftSystem.nuclei[a[p].nucleusIndex].positionBohr,rb=rightSystem.nuclei[b[q].nucleusIndex].positionBohr
            for i in a[p].weights.indices { for j in b[q].weights.indices {
                result[p,q] += a[p].weights[i]*b[q].weights[j]*primitiveOverlap(a[p].primitiveExponents[i],a[p].angular,ra,
                                                                            b[q].primitiveExponents[j],b[q].angular,rb)
            } }
        } }
        guard result.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("cross-overlap overflow") }; return result
    }
}
