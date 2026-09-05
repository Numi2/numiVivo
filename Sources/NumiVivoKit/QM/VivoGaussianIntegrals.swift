import Foundation

/// An electronic-structure coordinate, explicitly in bohr; not a new molecular
/// identity. sourceAtomID links back to VivoStructure through the QM/MM adapter.
public struct VivoQMNucleus: Codable, Sendable, Equatable {
    public var sourceAtomID: String
    public var nuclearCharge: Int
    public var positionBohr: [Double]
    public init(sourceAtomID: String, nuclearCharge: Int, positionBohr: [Double]) {
        self.sourceAtomID = sourceAtomID; self.nuclearCharge = nuclearCharge
        self.positionBohr = positionBohr
    }
}

public struct VivoQMPointCharge: Codable, Sendable, Equatable {
    public var sourceParticleID: String
    public var chargeE: Double
    public var positionBohr: [Double]
    public init(sourceParticleID: String, chargeE: Double, positionBohr: [Double]) {
        self.sourceParticleID = sourceParticleID; self.chargeE = chargeE
        self.positionBohr = positionBohr
    }
}

public struct VivoGaussianPrimitive: Codable, Sendable, Equatable {
    /// Coefficient multiplies a normalized Cartesian primitive.
    public var exponent: Double
    public var coefficient: Double
    public init(exponent: Double, coefficient: Double) {
        self.exponent = exponent; self.coefficient = coefficient
    }
}

public struct VivoCartesianGaussian: Codable, Sendable, Equatable {
    public var nucleusIndex: Int
    public var angularMomentum: [Int]
    public var primitives: [VivoGaussianPrimitive]
    public init(nucleusIndex: Int, angularMomentum: [Int], primitives: [VivoGaussianPrimitive]) {
        self.nucleusIndex = nucleusIndex; self.angularMomentum = angularMomentum
        self.primitives = primitives
    }
}

public struct VivoGaussianProblem: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/cartesian-gaussian-problem/v1"
    public var schema: String
    public var basisIdentity: String
    public var nuclei: [VivoQMNucleus]
    public var pointCharges: [VivoQMPointCharge]
    public var basis: [VivoCartesianGaussian]
    public var molecularCharge: Int
    /// N_alpha - N_beta; a fixed M_s sector, not a guarantee of spin purity.
    public var spinProjectionTwice: Int

    public init(basisIdentity: String, nuclei: [VivoQMNucleus],
                pointCharges: [VivoQMPointCharge] = [], basis: [VivoCartesianGaussian],
                molecularCharge: Int = 0, spinProjectionTwice: Int = 0) {
        schema = Self.schema; self.basisIdentity = basisIdentity
        self.nuclei = nuclei; self.pointCharges = pointCharges; self.basis = basis
        self.molecularCharge = molecularCharge; self.spinProjectionTwice = spinProjectionTwice
    }
    public var electronCount: Int { nuclei.reduce(0) { $0 + $1.nuclearCharge } - molecularCharge }
    public var alphaElectrons: Int { (electronCount + spinProjectionTwice) / 2 }
    public var betaElectrons: Int { (electronCount - spinProjectionTwice) / 2 }
    public func validate() throws {
        guard schema == Self.schema, !basisIdentity.isEmpty, !nuclei.isEmpty,
              nuclei.count <= 256, (1...48).contains(basis.count),
              (-512...512).contains(molecularCharge), (-96...96).contains(spinProjectionTwice),
              Set(nuclei.map(\.sourceAtomID)).count == nuclei.count else {
            throw VivoQMError.invalid("Gaussian problem schema, identity, or size")
        }
        func position(_ r: [Double]) -> Bool { r.count == 3 && r.allSatisfy(\.isFinite) }
        for atom in nuclei {
            guard !atom.sourceAtomID.isEmpty, (1...118).contains(atom.nuclearCharge),
                  position(atom.positionBohr) else { throw VivoQMError.invalid("QM nucleus") }
        }
        for charge in pointCharges {
            guard !charge.sourceParticleID.isEmpty, charge.chargeE.isFinite,
                  position(charge.positionBohr) else { throw VivoQMError.invalid("MM point charge") }
        }
        guard pointCharges.count <= 100_000,
              Set(pointCharges.map(\.sourceParticleID)).count == pointCharges.count,
              electronCount >= 0, (electronCount + spinProjectionTwice) % 2 == 0,
              (0...basis.count).contains(alphaElectrons), (0...basis.count).contains(betaElectrons) else {
            throw VivoQMError.invalid("electron number / M_s sector / point-charge identities")
        }
        for function in basis {
            guard nuclei.indices.contains(function.nucleusIndex), function.angularMomentum.count == 3,
                  function.angularMomentum.allSatisfy({ (0...3).contains($0) }),
                  function.angularMomentum.reduce(0, +) <= 3,
                  (1...32).contains(function.primitives.count),
                  function.primitives.allSatisfy({ $0.exponent.isFinite && $0.exponent > 0 && $0.coefficient.isFinite }),
                  function.primitives.contains(where: { $0.coefficient != 0 }) else {
                throw VivoQMError.invalid("Cartesian s/p/d/f Gaussian function")
            }
        }
    }
}

public struct VivoGaussianIntegralSet: Codable, Sendable {
    public var problem: VivoGaussianProblem
    public var overlap: VivoQMMatrix
    public var coreHamiltonian: VivoQMMatrix
    public var electronRepulsion: [Double]
    public var constantEnergyHartree: Double
    /// Cartesian AO position matrices <mu|x/y/z|nu>, in bohr.
    public var positionIntegrals: [VivoQMMatrix]
    public var schwarzThreshold: Double

    public func validate() throws {
        try problem.validate(); try overlap.requireSymmetric(); try coreHamiltonian.requireSymmetric()
        let n = problem.basis.count
        guard overlap.rows == n, coreHamiltonian.rows == n,
              electronRepulsion.count == n * n * n * n,
              electronRepulsion.allSatisfy(\.isFinite), constantEnergyHartree.isFinite,
              positionIntegrals.count == 3, schwarzThreshold.isFinite, schwarzThreshold >= 0 else {
            throw VivoQMError.invalid("Gaussian integral set dimensions or values")
        }
        _ = try VivoEmbeddedHamiltonian(orbitalIDs:(0..<n).map { "ao:\($0)" },basisIdentity:problem.basisIdentity,
                                        alphaElectrons:problem.alphaElectrons,betaElectrons:problem.betaElectrons,
                                        constantEnergyHartree:constantEnergyHartree,oneElectron:coreHamiltonian,
                                        electronRepulsion:electronRepulsion)
        for matrix in positionIntegrals {
            try matrix.requireSymmetric()
            guard matrix.rows == n else { throw VivoQMError.invalid("position-integral shape") }
        }
    }
}

/// Native FP64 McMurchie-Davidson/Hermite integrals. Cartesian functions through f;
/// no ECPs, spherical/cartesian guessing, density fitting, or hidden DFT substitutes.
public enum VivoGaussianIntegralEngine {
    private struct Primitive {
        let exponent: Double
        let weight: Double
        let center: [Double]
        let angular: [Int]
    }
    private struct HermiteKey: Hashable { let i: Int; let j: Int; let t: Int }
    private struct CoulombKey: Hashable { let t: Int; let u: Int; let v: Int; let n: Int }

    private static func doubleFactorial(_ n: Int) -> Double {
        if n <= 0 { return 1 }
        return stride(from: n, through: 1, by: -2).reduce(1.0) { $0 * Double($1) }
    }
    private static func normalization(_ exponent: Double, _ angular: [Int]) -> Double {
        let l = angular.reduce(0, +)
        let denominator = angular.reduce(1.0) { $0 * doubleFactorial(2 * $1 - 1) }
        return pow(2 * exponent / .pi, 0.75) * sqrt(pow(4 * exponent, Double(l)) / denominator)
    }
    private static func hermite(_ i: Int, _ j: Int, _ displacement: Double,
                                _ alpha: Double, _ beta: Double) -> [Double] {
        guard i >= 0, j >= 0 else { return [] }
        let p = alpha + beta, reduced = alpha * beta / p
        var cache: [HermiteKey: Double] = [:]
        func evaluate(_ a: Int, _ b: Int, _ t: Int) -> Double {
            if t < 0 || t > a + b || a < 0 || b < 0 { return 0 }
            if a == 0 && b == 0 { return t == 0 ? exp(-reduced * displacement * displacement) : 0 }
            let key = HermiteKey(i: a, j: b, t: t)
            if let value = cache[key] { return value }
            let value: Double
            if b == 0 {
                value = evaluate(a - 1, b, t - 1) / (2 * p)
                    - reduced * displacement / alpha * evaluate(a - 1, b, t)
                    + Double(t + 1) * evaluate(a - 1, b, t + 1)
            } else {
                value = evaluate(a, b - 1, t - 1) / (2 * p)
                    + reduced * displacement / beta * evaluate(a, b - 1, t)
                    + Double(t + 1) * evaluate(a, b - 1, t + 1)
            }
            cache[key] = value
            return value
        }
        return (0...(i + j)).map { evaluate(i, j, $0) }
    }
    /// Positive-series evaluation avoids cancellation at small T; upward
    /// recurrence is used only beyond T=30 (maximum required order is 12).
    private static func boys(_ order: Int, _ t: Double) -> Double {
        if t <= 30 {
            var term = 1 / Double(2 * order + 1), sum = term
            for k in 1...256 {
                term *= 2 * t / Double(2 * order + 2 * k + 1)
                sum += term
                if abs(term) < abs(sum) * 2e-16 { break }
            }
            return exp(-t) * sum
        }
        var value = 0.5 * sqrt(.pi / t) * erf(sqrt(t))
        if order > 0 {
            for n in 0..<order { value = (Double(2 * n + 1) * value - exp(-t)) / (2 * t) }
        }
        return value
    }
    private static func coulomb(_ alpha: Double, _ displacement: [Double]) -> (Int, Int, Int) -> Double {
        let argument = alpha * displacement.reduce(0) { $0 + $1 * $1 }
        var cache: [CoulombKey: Double] = [:]
        func r(_ t: Int, _ u: Int, _ v: Int, _ n: Int) -> Double {
            let key = CoulombKey(t: t, u: u, v: v, n: n)
            if let value = cache[key] { return value }
            let value: Double
            if t == 0 && u == 0 && v == 0 {
                value = pow(-2 * alpha, Double(n)) * boys(n, argument)
            } else if t > 0 {
                value = (t > 1 ? Double(t - 1) * r(t - 2, u, v, n + 1) : 0)
                    + displacement[0] * r(t - 1, u, v, n + 1)
            } else if u > 0 {
                value = (u > 1 ? Double(u - 1) * r(t, u - 2, v, n + 1) : 0)
                    + displacement[1] * r(t, u - 1, v, n + 1)
            } else {
                value = (v > 1 ? Double(v - 1) * r(t, u, v - 2, n + 1) : 0)
                    + displacement[2] * r(t, u, v - 1, n + 1)
            }
            cache[key] = value
            return value
        }
        return { t, u, v in r(t, u, v, 0) }
    }
    private static func overlap(_ a: Primitive, _ b: Primitive, angularB: [Int]? = nil) -> Double {
        let l = angularB ?? b.angular
        guard l.allSatisfy({ $0 >= 0 }) else { return 0 }
        var value = pow(.pi / (a.exponent + b.exponent), 1.5)
        for axis in 0..<3 {
            value *= hermite(a.angular[axis], l[axis], a.center[axis] - b.center[axis], a.exponent, b.exponent)[0]
        }
        return value
    }
    private static func kinetic(_ a: Primitive, _ b: Primitive) -> Double {
        let beta = b.exponent
        var value = beta * Double(2 * b.angular.reduce(0, +) + 3) * overlap(a, b)
        for axis in 0..<3 {
            var raised = b.angular, lowered = b.angular
            raised[axis] += 2; lowered[axis] -= 2
            value -= 2 * beta * beta * overlap(a, b, angularB: raised)
            value -= 0.5 * Double(b.angular[axis] * (b.angular[axis] - 1)) * overlap(a, b, angularB: lowered)
        }
        return value
    }
    private static func attraction(_ a: Primitive, _ b: Primitive, at center: [Double]) -> Double {
        let p = a.exponent + b.exponent
        let product = (0..<3).map { (a.exponent * a.center[$0] + b.exponent * b.center[$0]) / p }
        let coefficients = (0..<3).map {
            hermite(a.angular[$0], b.angular[$0], a.center[$0] - b.center[$0], a.exponent, b.exponent)
        }
        let r = coulomb(p, zip(product, center).map(-))
        var value = 0.0
        for t in coefficients[0].indices { for u in coefficients[1].indices { for v in coefficients[2].indices {
            value += coefficients[0][t] * coefficients[1][u] * coefficients[2][v] * r(t, u, v)
        } } }
        return 2 * .pi / p * value
    }
    private static func repulsion(_ a: Primitive, _ b: Primitive, _ c: Primitive, _ d: Primitive) -> Double {
        let p = a.exponent + b.exponent, q = c.exponent + d.exponent
        let productAB = (0..<3).map { (a.exponent * a.center[$0] + b.exponent * b.center[$0]) / p }
        let productCD = (0..<3).map { (c.exponent * c.center[$0] + d.exponent * d.center[$0]) / q }
        let ab = (0..<3).map { hermite(a.angular[$0], b.angular[$0], a.center[$0] - b.center[$0], a.exponent, b.exponent) }
        let cd = (0..<3).map { hermite(c.angular[$0], d.angular[$0], c.center[$0] - d.center[$0], c.exponent, d.exponent) }
        let r = coulomb(p * q / (p + q), zip(productAB, productCD).map(-))
        var value = 0.0
        for t in ab[0].indices { for u in ab[1].indices { for v in ab[2].indices {
            for x in cd[0].indices { for y in cd[1].indices { for z in cd[2].indices {
                let sign = (x + y + z) % 2 == 0 ? 1.0 : -1.0
                value += ab[0][t] * ab[1][u] * ab[2][v] * cd[0][x] * cd[1][y] * cd[2][z]
                    * sign * r(t + x, u + y, v + z)
            } } }
        } } }
        return 2 * pow(.pi, 2.5) / (p * q * sqrt(p + q)) * value
    }

    private static func normalizedFunctions(_ problem: VivoGaussianProblem) throws -> [[Primitive]] {
        var functions: [[Primitive]] = []
        for function in problem.basis {
            let center = problem.nuclei[function.nucleusIndex].positionBohr
            var primitives = function.primitives.map {
                Primitive(exponent: $0.exponent,
                          weight: $0.coefficient * normalization($0.exponent, function.angularMomentum),
                          center: center, angular: function.angularMomentum)
            }
            var norm = 0.0
            for a in primitives { for b in primitives { norm += a.weight * b.weight * overlap(a, b) } }
            guard norm.isFinite, norm > 1e-20 else { throw VivoQMError.invalid("singular Gaussian contraction") }
            primitives = primitives.map { .init(exponent: $0.exponent, weight: $0.weight / sqrt(norm), center: $0.center, angular: $0.angular) }
            functions.append(primitives)
        }
        return functions
    }

    public static func crossOverlap(reference: VivoGaussianProblem, current: VivoGaussianProblem) throws -> VivoQMMatrix {
        try reference.validate(); try current.validate()
        let left = try normalizedFunctions(reference), right = try normalizedFunctions(current)
        var result = VivoQMMatrix(left.count,right.count)
        for i in left.indices { for j in right.indices {
            for a in left[i] { for b in right[j] { result[i,j] += a.weight*b.weight*overlap(a,b) } }
        } }
        try result.validate()
        return result
    }

    public static func build(_ problem: VivoGaussianProblem, schwarzThreshold: Double = 1e-12,
                             maximumPrimitiveQuartets: Int = 20_000_000,
                             maximumPrimitivePairCenters: Int = 50_000_000) throws -> VivoGaussianIntegralSet {
        try problem.validate()
        guard schwarzThreshold.isFinite, (0...1e-6).contains(schwarzThreshold),
              (1...1_000_000_000).contains(maximumPrimitiveQuartets),
              (1...1_000_000_000).contains(maximumPrimitivePairCenters) else {
            throw VivoQMError.invalid("integral screening or work budget")
        }
        let functions = try normalizedFunctions(problem)
        let n = functions.count
        var s = VivoQMMatrix(n, n), h = VivoQMMatrix(n, n)
        var dipoles = (0..<3).map { _ in VivoQMMatrix(n, n) }
        let charges = problem.nuclei.map { (Double($0.nuclearCharge), $0.positionBohr) }
            + problem.pointCharges.filter { $0.chargeE != 0 }.map { ($0.chargeE, $0.positionBohr) }
        let primitivePairs = functions.indices.reduce(0) { sum,i in
            sum+(0...i).reduce(0) { $0+functions[i].count*functions[$1].count }
        }
        guard primitivePairs <= maximumPrimitivePairCenters/max(1,charges.count) else {
            throw VivoQMError.capacity("one-electron primitive-pair/charge-center budget exhausted")
        }
        var constant = 0.0
        for i in problem.nuclei.indices {
            let atom = problem.nuclei[i]
            for j in 0..<i {
                let other = problem.nuclei[j]
                let distance = sqrt(zip(atom.positionBohr, other.positionBohr).reduce(0) { $0 + pow($1.0 - $1.1, 2) })
                guard distance > 1e-10 else { throw VivoQMError.invalid("coincident nuclei") }
                constant += Double(atom.nuclearCharge * other.nuclearCharge) / distance
            }
            for charge in problem.pointCharges where charge.chargeE != 0 {
                let distance = sqrt(zip(atom.positionBohr, charge.positionBohr).reduce(0) { $0 + pow($1.0 - $1.1, 2) })
                guard distance > 1e-8 else { throw VivoQMError.invalid("MM point charge coincides with a QM nucleus") }
                constant += Double(atom.nuclearCharge) * charge.chargeE / distance
            }
        }
        for i in 0..<n { for j in 0...i {
            for a in functions[i] { for b in functions[j] {
                let weight = a.weight * b.weight
                let sab = overlap(a, b)
                s[i, j] += weight * sab
                var core = kinetic(a, b)
                for (charge, center) in charges { core -= charge * attraction(a, b, at: center) }
                h[i, j] += weight * core
                for axis in 0..<3 {
                    var raised = b.angular; raised[axis] += 1
                    dipoles[axis][i, j] += weight * (overlap(a, b, angularB: raised) + b.center[axis] * sab)
                }
            } }
            s[j, i] = s[i, j]; h[j, i] = h[i, j]
            for axis in 0..<3 { dipoles[axis][j, i] = dipoles[axis][i, j] }
        } }
        var g = Array(repeating: 0.0, count: n * n * n * n)
        let pairs = (0..<n).flatMap { i in (0...i).map { (i, $0) } }
        var work = 0
        func contracted(_ i: Int, _ j: Int, _ k: Int, _ l: Int) throws -> Double {
            let requested = functions[i].count * functions[j].count * functions[k].count * functions[l].count
            guard requested <= maximumPrimitiveQuartets - work else { throw VivoQMError.capacity("primitive ERI quartet budget exhausted") }
            work += requested
            var value = 0.0
            for a in functions[i] { for b in functions[j] { for c in functions[k] { for d in functions[l] {
                value += a.weight * b.weight * c.weight * d.weight * repulsion(a, b, c, d)
            } } } }
            return value
        }
        var bounds: [Double] = []
        for (i, j) in pairs { bounds.append(sqrt(max(0, try contracted(i, j, i, j)))) }
        for (a, pairAB) in pairs.enumerated() {
            let (i, j) = pairAB
            for b in 0...a {
                if bounds[a] * bounds[b] < schwarzThreshold { continue }
                let (k, l) = pairs[b]
                let value = try contracted(i, j, k, l)
                for (p, q, r, t) in [(i,j,k,l),(j,i,k,l),(i,j,l,k),(j,i,l,k),
                                     (k,l,i,j),(l,k,i,j),(k,l,j,i),(l,k,j,i)] {
                    g[((p * n + q) * n + r) * n + t] = value
                }
            }
        }
        let result = VivoGaussianIntegralSet(problem: problem, overlap: s, coreHamiltonian: h,
                                              electronRepulsion: g, constantEnergyHartree: constant,
                                              positionIntegrals: dipoles, schwarzThreshold: schwarzThreshold)
        try result.validate()
        return result
    }
}
