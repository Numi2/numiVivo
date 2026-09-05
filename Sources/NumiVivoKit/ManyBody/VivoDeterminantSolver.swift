import Foundation

public struct VivoFermionLadder: Codable, Sendable, Equatable {
    public let mode: Int
    public let creation: Bool
    public init(mode: Int, creation: Bool) { self.mode = mode; self.creation = creation }
}

public struct VivoFermionTerm: Sendable {
    public let coefficient: Double
    /// Written operator order, left to right. Applied to a ket right to left.
    public let operators: [VivoFermionLadder]
}

public enum VivoFermionAlgebra {
    public static func apply(_ operators: [VivoFermionLadder], to determinant: UInt64) throws -> (UInt64, Double)? {
        guard operators.allSatisfy({ (0..<63).contains($0.mode) }) else { throw VivoQMError.invalid("fermion mode outside UInt64 representation") }
        var state = determinant, phase = 1.0
        for op in operators.reversed() {
            guard let step = ladder(op.mode, creation: op.creation, state: state) else { return nil }
            state = step.0; phase *= step.1
        }
        return (state, phase)
    }
    internal static func ladder(_ mode: Int, creation: Bool, state: UInt64) -> (UInt64, Double)? {
        let bit = UInt64(1) << mode
        let occupied = state & bit != 0
        if creation == occupied { return nil }
        let phase = (state & (bit - 1)).nonzeroBitCount % 2 == 0 ? 1.0 : -1.0
        return (state ^ bit, phase)
    }
    public static func terms(for h: VivoEmbeddedHamiltonian) throws -> [VivoFermionTerm] {
        try h.validate()
        guard h.orbitalCount <= 12 else { throw VivoQMError.capacity("explicit fermion expansion is limited to 12 spatial orbitals") }
        let n = h.orbitalCount
        var terms: [VivoFermionTerm] = []
        for p in 0..<n { for q in 0..<n {
            if h.oneElectron[p, q] != 0 {
                for spin in 0..<2 {
                    terms.append(.init(coefficient: h.oneElectron[p, q], operators: [
                        .init(mode: 2 * p + spin, creation: true), .init(mode: 2 * q + spin, creation: false)]))
                }
            }
            for r in 0..<n { for s in 0..<n {
                let value = 0.5 * h.eri(p, q, r, s)
                if value == 0 { continue }
                for sigma in 0..<2 { for tau in 0..<2 {
                    if 2 * p + sigma == 2 * r + tau || 2 * s + tau == 2 * q + sigma { continue }
                    terms.append(.init(coefficient: value, operators: [
                        .init(mode: 2 * p + sigma, creation: true), .init(mode: 2 * r + tau, creation: true),
                        .init(mode: 2 * s + tau, creation: false), .init(mode: 2 * q + sigma, creation: false)]))
                } }
            } }
        } }
        return terms
    }
}

public enum VivoDeterminantMethod: String, Codable, Sendable { case fci, cisd }
public struct VivoDeterminantConfiguration: Codable, Sendable {
    public var method: VivoDeterminantMethod
    public var maximumDeterminants: Int
    public var maximumOperatorApplications: Int
    public var residualTolerance: Double
    public var referenceDeterminant: UInt64?
    public init(method: VivoDeterminantMethod = .fci, maximumDeterminants: Int = 512,
                maximumOperatorApplications: Int = 50_000_000, residualTolerance: Double = 1e-8,
                referenceDeterminant: UInt64? = nil) {
        self.method = method; self.maximumDeterminants = maximumDeterminants
        self.maximumOperatorApplications = maximumOperatorApplications; self.residualTolerance = residualTolerance
        self.referenceDeterminant = referenceDeterminant
    }
    public func validate() throws {
        guard (1...512).contains(maximumDeterminants), (1...500_000_000).contains(maximumOperatorApplications),
              residualTolerance.isFinite, (1e-12...1e-3).contains(residualTolerance),
              method == .cisd || referenceDeterminant == nil else {
            throw VivoQMError.invalid("determinant solver resource or residual contract")
        }
    }
}

public struct VivoDeterminantState: Codable, Sendable {
    public var spatialOrbitalCount: Int
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public var determinants: [UInt64]
    public var amplitudes: [Double]
    public func validate() throws {
        let n = spatialOrbitalCount
        guard (1...12).contains(n), (0...n).contains(alphaElectrons), (0...n).contains(betaElectrons),
              !determinants.isEmpty, determinants.count <= 512,
              amplitudes.count == determinants.count, amplitudes.allSatisfy(\.isFinite),
              Set(determinants).count == determinants.count else { throw VivoQMError.invalid("determinant-state dimensions") }
        let alphaMask = (0..<n).reduce(UInt64(0)) { $0 | (UInt64(1) << (2 * $1)) }
        let betaMask = alphaMask << 1
        for state in determinants {
            guard state >> (2 * n) == 0, (state & alphaMask).nonzeroBitCount == alphaElectrons,
                  (state & betaMask).nonzeroBitCount == betaElectrons else { throw VivoQMError.invalid("determinant outside the declared electron sector") }
        }
        guard abs(amplitudes.reduce(0) { $0 + $1 * $1 } - 1) < 1e-8 else { throw VivoQMError.invalid("wavefunction is not normalized") }
    }
    /// Exact fermionic partial trace, including the sign from moving selected
    /// modes before their complement. Pair entropy is NOT reconstructed from
    /// the 1-/2-particle RDM alone; that is insufficient for a general CI state.
    public func orbitalDensity(spatialOrbitals: [Int]) throws -> VivoQMMatrix {
        try validate()
        guard (1...2).contains(spatialOrbitals.count), Set(spatialOrbitals).count == spatialOrbitals.count,
              spatialOrbitals.allSatisfy({ (0..<spatialOrbitalCount).contains($0) }) else { throw VivoQMError.invalid("one-/two-orbital partial-trace selection") }
        let modes = spatialOrbitals.sorted().flatMap { [2 * $0, 2 * $0 + 1] }
        let mask = modes.reduce(UInt64(0)) { $0 | (UInt64(1) << $1) }
        let dimension = 1 << modes.count
        var groups: [UInt64: [Double]] = [:]
        for (i, state) in determinants.enumerated() {
            let environment = state & ~mask
            var local = 0, parity = 0
            for (j, mode) in modes.enumerated() where state & (UInt64(1) << mode) != 0 {
                local |= 1 << j
                parity += (environment & ((UInt64(1) << mode) - 1)).nonzeroBitCount
            }
            if groups[environment] == nil { groups[environment] = Array(repeating: 0, count: dimension) }
            groups[environment]![local] += (parity % 2 == 0 ? 1 : -1) * amplitudes[i]
        }
        var rho = VivoQMMatrix(dimension, dimension)
        for key in groups.keys.sorted() {
            let vector = groups[key]!
            for p in 0..<dimension { for q in 0..<dimension { rho[p, q] += vector[p] * vector[q] } }
        }
        return rho
    }
}

public struct VivoDeterminantResult: Codable, Sendable {
    public let method: VivoDeterminantMethod
    public let energyHartree: Double
    public let residualNormHartree: Double
    public let excitationGapInSectorHartree: Double?
    public let referenceDeterminant: UInt64?
    public let state: VivoDeterminantState
}

public enum VivoDeterminantSolver {
    private static func combinations(_ n: Int, _ k: Int, spin: Int) -> [UInt64] {
        var result: [UInt64] = []
        func visit(_ start: Int, _ remaining: Int, _ bits: UInt64) {
            if remaining == 0 { result.append(bits); return }
            if n - start < remaining { return }
            for orbital in start...(n - remaining) { visit(orbital + 1, remaining - 1, bits | (UInt64(1) << (2 * orbital + spin))) }
        }
        visit(0, k, 0)
        return result
    }
    public static func solve(_ h: VivoEmbeddedHamiltonian,
                             configuration config: VivoDeterminantConfiguration = .init()) throws -> VivoDeterminantResult {
        try h.validate(); try config.validate()
        let n = h.orbitalCount
        guard n <= 12 else { throw VivoQMError.capacity("small CI solver supports at most 12 spatial orbitals") }
        let aa = combinations(n, h.alphaElectrons, spin: 0), bb = combinations(n, h.betaElectrons, spin: 1)
        if config.method == .fci && aa.count * bb.count > config.maximumDeterminants {
            throw VivoQMError.capacity("FCI sector has \(aa.count * bb.count) determinants; requested limit is \(config.maximumDeterminants)")
        }
        let referenceA = (0..<h.alphaElectrons).reduce(UInt64(0)) { $0 | (UInt64(1) << (2 * $1)) }
        let referenceB = (0..<h.betaElectrons).reduce(UInt64(0)) { $0 | (UInt64(1) << (2 * $1 + 1)) }
        let reference = config.referenceDeterminant ?? (referenceA | referenceB)
        let alphaMask = (0..<n).reduce(UInt64(0)) { $0 | (UInt64(1) << (2*$1)) }
        guard reference >> (2*n) == 0, (reference&alphaMask).nonzeroBitCount == h.alphaElectrons,
              (reference&(alphaMask<<1)).nonzeroBitCount == h.betaElectrons else {
            throw VivoQMError.invalid("CISD reference determinant is outside the requested electron sector")
        }
        var states: [UInt64] = []
        for a in aa { for b in bb {
            let state = a | b
            if config.method == .fci || (state ^ reference).nonzeroBitCount <= 4 { states.append(state) }
        } }
        states.sort()
        guard !states.isEmpty, states.count <= config.maximumDeterminants else { throw VivoQMError.capacity("CI determinant subspace exceeds the configured dense-solver budget") }
        let index = Dictionary(uniqueKeysWithValues: states.enumerated().map { ($1, $0) })
        let terms = try VivoFermionAlgebra.terms(for: h)
        guard terms.count <= config.maximumOperatorApplications / states.count else { throw VivoQMError.capacity("CI operator-application budget") }
        let dimension = states.count
        var matrix = VivoQMMatrix(dimension, dimension)
        for (column, determinant) in states.enumerated() {
            matrix[column, column] = h.constantEnergyHartree
            for term in terms {
                if let (target, phase) = try VivoFermionAlgebra.apply(term.operators, to: determinant), let row = index[target] {
                    matrix[row, column] += phase * term.coefficient
                }
            }
        }
        try matrix.requireSymmetric(tolerance: 1e-10)
        let eigen = try matrix.symmetricEigen(tolerance: min(1e-12, config.residualTolerance / Double(10 * dimension)))
        let vector = (0..<dimension).map { eigen.vectors[$0, 0] }
        var residual = 0.0
        for i in 0..<dimension {
            var value = -eigen.values[0] * vector[i]
            for j in 0..<dimension { value += matrix[i, j] * vector[j] }
            residual += value * value
        }
        residual = sqrt(residual)
        guard residual <= config.residualTolerance else { throw VivoQMError.convergence("CI eigenpair failed its physical residual check") }
        let state = VivoDeterminantState(spatialOrbitalCount: n, alphaElectrons: h.alphaElectrons,
                                         betaElectrons: h.betaElectrons, determinants: states, amplitudes: vector)
        try state.validate()
        return .init(method: config.method, energyHartree: eigen.values[0], residualNormHartree: residual,
                     excitationGapInSectorHartree: dimension > 1 ? eigen.values[1] - eigen.values[0] : nil,
                     referenceDeterminant:config.method == .cisd ? reference : nil, state: state)
    }
}

/// Spin-mode tensors: gamma[p,r]=<a†p ar>, Gamma[p,q,r,s]=<a†p a†q as ar>.
/// This differs from a spatial chemists'-ordered 2-RDM. Conversion is explicit.
public struct VivoFermionRDM: Codable, Sendable {
    public let electronCount: Int
    public let oneParticle: VivoQMMatrix
    public let twoParticle: [Double]
    public var modeCount: Int { oneParticle.rows }
    public func gamma2(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        let n = modeCount; return twoParticle[((p * n + q) * n + r) * n + s]
    }
    public func cumulant(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        gamma2(p, q, r, s) - oneParticle[p, r] * oneParticle[q, s] + oneParticle[p, s] * oneParticle[q, r]
    }
    public func validate(tolerance: Double = 1e-7) throws {
        try oneParticle.requireSymmetric()
        let n = modeCount
        guard (2...24).contains(n), n % 2 == 0, (0...n).contains(electronCount),
              twoParticle.count == n * n * n * n, twoParticle.allSatisfy(\.isFinite) else { throw VivoQMError.invalid("RDM dimensions") }
        let trace = (0..<n).reduce(0.0) { $0 + oneParticle[$1, $1] }
        guard abs(trace - Double(electronCount)) <= tolerance else { throw VivoQMError.invalid("1-RDM trace violates particle number") }
        for p in 0..<n { for r in 0..<n {
            let contracted = (0..<n).reduce(0.0) { $0 + gamma2(p, $1, r, $1) }
            guard abs(contracted - Double(electronCount - 1) * oneParticle[p, r]) < tolerance else {
                throw VivoQMError.invalid("2-RDM contraction violates (N-1) gamma")
            }
        } }
    }
    public var spatialDensity: VivoQMMatrix {
        let n = modeCount / 2
        var result = VivoQMMatrix(n, n)
        for p in 0..<n { for q in 0..<n { result[p, q] = oneParticle[2*p, 2*q] + oneParticle[2*p+1, 2*q+1] } }
        return result
    }
    public func energy(of h: VivoEmbeddedHamiltonian) throws -> Double {
        try validate(); try h.validate()
        guard 2 * h.orbitalCount == modeCount, h.electronCount == electronCount else { throw VivoQMError.invalid("RDM/Hamiltonian sector mismatch") }
        let n = h.orbitalCount
        var value = h.constantEnergyHartree
        for p in 0..<n { for q in 0..<n {
            for sigma in 0..<2 { value += h.oneElectron[p, q] * oneParticle[2*p+sigma, 2*q+sigma] }
            for r in 0..<n { for s in 0..<n { for sigma in 0..<2 { for tau in 0..<2 {
                value += 0.5 * h.eri(p,q,r,s) * gamma2(2*p+sigma,2*r+tau,2*q+sigma,2*s+tau)
            } } } }
        } }
        return value
    }
    public static func from(_ wavefunction: VivoDeterminantState,
                            maximumApplications: Int = 100_000_000) throws -> Self {
        try wavefunction.validate()
        guard (1...500_000_000).contains(maximumApplications) else { throw VivoQMError.invalid("RDM work budget") }
        let n = 2 * wavefunction.spatialOrbitalCount
        let electrons = wavefunction.alphaElectrons + wavefunction.betaElectrons
        let estimated = wavefunction.determinants.count * max(1, electrons) * max(1, electrons - 1)
            * (n - electrons + 2) * (n - electrons + 1)
        guard estimated <= maximumApplications else { throw VivoQMError.capacity("2-RDM construction work budget") }
        let index = Dictionary(uniqueKeysWithValues: wavefunction.determinants.enumerated().map { ($1, $0) })
        var one = VivoQMMatrix(n, n), two = Array(repeating: 0.0, count: n*n*n*n)
        for (column, determinant) in wavefunction.determinants.enumerated() {
            let coefficient = wavefunction.amplitudes[column]
            for r in 0..<n {
                guard let (d1, phase1) = VivoFermionAlgebra.ladder(r, creation: false, state: determinant) else { continue }
                for p in 0..<n {
                    if let (target, phase2) = VivoFermionAlgebra.ladder(p, creation: true, state: d1), let row = index[target] {
                        one[p,r] += coefficient * wavefunction.amplitudes[row] * phase1 * phase2
                    }
                }
                for s in 0..<n {
                    guard let (d2, phase2) = VivoFermionAlgebra.ladder(s, creation: false, state: d1) else { continue }
                    for q in 0..<n {
                        guard let (d3, phase3) = VivoFermionAlgebra.ladder(q, creation: true, state: d2) else { continue }
                        for p in 0..<n {
                            if let (target, phase4) = VivoFermionAlgebra.ladder(p, creation: true, state: d3), let row = index[target] {
                                two[((p*n+q)*n+r)*n+s] += coefficient * wavefunction.amplitudes[row] * phase1 * phase2 * phase3 * phase4
                            }
                        }
                    }
                }
            }
        }
        let result = Self(electronCount: electrons, oneParticle: one, twoParticle: two)
        try result.validate()
        return result
    }
}

public struct VivoOrbitalInformation: Codable, Sendable {
    public let singleOrbitalEntropy: [Double]
    /// Natural logarithms, I_ij = S_i + S_j - S_ij, with zero diagonal.
    /// No implicit 1/2 normalization; objective weights use this convention.
    public let mutualInformation: VivoQMMatrix
    private static func entropy(_ rho: VivoQMMatrix) throws -> Double {
        let eigen = try rho.symmetricEigen()
        guard abs(eigen.values.reduce(0,+) - 1) < 1e-8,
              eigen.values.allSatisfy({ $0 >= -1e-10 && $0 <= 1 + 1e-10 }) else { throw VivoQMError.invalid("orbital RDM is not a normalized positive density operator") }
        return eigen.values.reduce(0) { $0 - ($1 > 1e-15 ? $1 * log($1) : 0) }
    }
    public static func from(_ state: VivoDeterminantState) throws -> Self {
        try state.validate()
        let n = state.spatialOrbitalCount
        let singles = try (0..<n).map { try entropy(state.orbitalDensity(spatialOrbitals: [$0])) }
        var mi = VivoQMMatrix(n,n)
        for i in 0..<n { for j in 0..<i {
            let pair = try entropy(state.orbitalDensity(spatialOrbitals: [i,j]))
            let value = singles[i] + singles[j] - pair
            guard value >= -1e-8 else { throw VivoQMError.invalid("negative mutual information beyond roundoff") }
            mi[i,j] = max(0,value); mi[j,i] = mi[i,j]
        } }
        return .init(singleOrbitalEntropy: singles, mutualInformation: mi)
    }
}
