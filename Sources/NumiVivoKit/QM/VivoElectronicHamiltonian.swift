import Foundation

/// Real spatial orbitals, interleaved spin modes (2p=alpha, 2p+1=beta).
/// g[p,q,r,s] = (pq|rs), in chemists' notation, row-major.
/// H = E0 + sum h[p,q] a†[p,sigma] a[q,sigma]
///       + 1/2 sum g[p,q,r,s] a†[p,sigma] a†[r,tau] a[s,tau] a[q,sigma].
/// E0 includes every declared constant exactly once. Solvers return E0 + <H_e>.
public struct VivoEmbeddedHamiltonian: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/embedded-hamiltonian/v1"
    public var schema: String
    public var orbitalIDs: [String]
    public var basisIdentity: String
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public var constantEnergyHartree: Double
    public var oneElectron: VivoQMMatrix
    public var electronRepulsion: [Double]
    public var provenance: [String: String]

    public init(orbitalIDs: [String], basisIdentity: String,
                alphaElectrons: Int, betaElectrons: Int,
                constantEnergyHartree: Double, oneElectron: VivoQMMatrix,
                electronRepulsion: [Double], provenance: [String: String] = [:]) throws {
        schema = Self.schema; self.orbitalIDs = orbitalIDs; self.basisIdentity = basisIdentity
        self.alphaElectrons = alphaElectrons; self.betaElectrons = betaElectrons
        self.constantEnergyHartree = constantEnergyHartree; self.oneElectron = oneElectron
        self.electronRepulsion = electronRepulsion; self.provenance = provenance
        try validate()
    }
    public var orbitalCount: Int { orbitalIDs.count }
    public var electronCount: Int { alphaElectrons + betaElectrons }
    public func eri(_ p: Int, _ q: Int, _ r: Int, _ s: Int) -> Double {
        let n = orbitalCount
        return electronRepulsion[((p * n + q) * n + r) * n + s]
    }
    public func validate() throws {
        let n = orbitalCount
        guard schema == Self.schema, (1...48).contains(n),
              Set(orbitalIDs).count == n, orbitalIDs.allSatisfy({ !$0.isEmpty }),
              !basisIdentity.isEmpty, (0...n).contains(alphaElectrons),
              (0...n).contains(betaElectrons), constantEnergyHartree.isFinite,
              oneElectron.rows == n, oneElectron.columns == n,
              electronRepulsion.count == n * n * n * n,
              electronRepulsion.allSatisfy(\.isFinite) else {
            throw VivoQMError.invalid("Hamiltonian schema, identities, electron sector, or tensor shape")
        }
        try oneElectron.requireSymmetric()
        let scale = max(1, electronRepulsion.map(abs).max() ?? 0)
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            let value = eri(p, q, r, s)
            guard abs(value - eri(q, p, r, s)) <= 1e-9 * scale,
                  abs(value - eri(p, q, s, r)) <= 1e-9 * scale,
                  abs(value - eri(r, s, p, q)) <= 1e-9 * scale else {
                throw VivoQMError.invalid("ERI tensor is not real chemists' eightfold-symmetric")
            }
        } } } }
    }

    /// Sequential O(n^5) FP64 AO/MO transform; no O(n^8) contraction.
    public func rotated(by c: VivoQMMatrix, orbitalIDs newIDs: [String]? = nil) throws -> Self {
        try validate(); try c.validate()
        let n = orbitalCount
        guard c.rows == n, c.columns == n else { throw VivoQMError.invalid("orbital rotation shape") }
        let gram = try c.transposed.multiplied(by: c)
        guard try gram.adding(.identity(n), scale: -1).frobeniusNorm < 1e-8 else {
            throw VivoQMError.invalid("embedded orbital rotation is not orthogonal")
        }
        return try .init(orbitalIDs: newIDs ?? orbitalIDs, basisIdentity: basisIdentity,
                         alphaElectrons: alphaElectrons, betaElectrons: betaElectrons,
                         constantEnergyHartree: constantEnergyHartree,
                         oneElectron: oneElectron.transformed(by: c),
                         electronRepulsion: Self.transformERI(electronRepulsion, dimension: n, coefficients: c),
                         provenance: provenance)
    }

    internal static func transformERI(_ source: [Double], dimension n: Int,
                                      coefficients c: VivoQMMatrix) throws -> [Double] {
        try c.validate()
        let m = c.columns
        guard (1...48).contains(n), (1...48).contains(m), c.rows == n,
              source.count == n * n * n * n, source.allSatisfy(\.isFinite) else {
            throw VivoQMError.invalid("ERI transformation dimensions")
        }
        var data = source
        var dims = [n, n, n, n]
        for axis in 0..<4 {
            let before = dims.prefix(axis).reduce(1, *)
            let after = dims.suffix(3 - axis).reduce(1, *)
            var next = Array(repeating: 0.0, count: before * m * after)
            for b in 0..<before { for j in 0..<m { for i in 0..<n {
                let coefficient = c[i, j]
                for a in 0..<after {
                    next[(b * m + j) * after + a] += coefficient * data[(b * n + i) * after + a]
                }
            } } }
            data = next; dims[axis] = m
        }
        guard data.allSatisfy(\.isFinite) else { throw VivoQMError.convergence("non-finite ERI transform") }
        return data
    }

    /// Closed-shell frozen-core contraction. Unlisted orbitals are empty virtuals,
    /// not an implicitly occupied environment. This is CASCI preprocessing, not DMET.
    public func activeSpace(active: [Int], doublyOccupiedCore core: [Int]) throws -> Self {
        try validate()
        let n = orbitalCount
        guard !active.isEmpty, Set(active).count == active.count,
              Set(core).count == core.count, Set(active).isDisjoint(with: core),
              (active + core).allSatisfy({ (0..<n).contains($0) }),
              alphaElectrons >= core.count, betaElectrons >= core.count,
              alphaElectrons - core.count <= active.count,
              betaElectrons - core.count <= active.count else {
            throw VivoQMError.invalid("active/frozen-core partition or electron count")
        }
        var constant = constantEnergyHartree
        for i in core {
            constant += 2 * oneElectron[i, i]
            for j in core { constant += 2 * eri(i, i, j, j) - eri(i, j, j, i) }
        }
        let m = active.count
        var h = VivoQMMatrix(m, m)
        var g = Array(repeating: 0.0, count: m * m * m * m)
        for (p, a) in active.enumerated() { for (q, b) in active.enumerated() {
            h[p, q] = oneElectron[a, b]
            for i in core { h[p, q] += 2 * eri(a, b, i, i) - eri(a, i, i, b) }
            for (r, c) in active.enumerated() { for (s, d) in active.enumerated() {
                g[((p * m + q) * m + r) * m + s] = eri(a, b, c, d)
            } }
        } }
        var origin = provenance
        origin["partition"] = "closed-shell-frozen-core"
        origin["frozenCoreOrbitals"] = core.map { orbitalIDs[$0] }.joined(separator: ",")
        return try .init(orbitalIDs: active.map { orbitalIDs[$0] }, basisIdentity: basisIdentity,
                         alphaElectrons: alphaElectrons - core.count,
                         betaElectrons: betaElectrons - core.count,
                         constantEnergyHartree: constant, oneElectron: h,
                         electronRepulsion: g, provenance: origin)
    }
}
