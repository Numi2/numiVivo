import Foundation

public struct VivoOrbitalLocalizationResult: Codable, Sendable {
    public let rotation: VivoQMMatrix
    public let objective: Double
    public let sweeps: Int
    public let converged: Bool
}

public enum VivoOrbitalSpace {
    /// Boys: supply three MO position-integral matrices. Pipek-Mezey: supply
    /// one symmetric MO atomic-population matrix per atom. Both maximize the
    /// squared diagonals; no orbital occupation is inferred from an index.
    public static func localize(metrics: [VivoQMMatrix], groups: [[Int]],
                                maximumSweeps: Int = 80, tolerance: Double = 1e-10) throws -> VivoOrbitalLocalizationResult {
        guard let first = metrics.first, metrics.count <= 256,
              (1...200).contains(maximumSweeps), tolerance.isFinite,
              (1e-14...1e-3).contains(tolerance) else { throw VivoQMError.invalid("localization settings") }
        let n = first.rows
        guard n <= 48 else { throw VivoQMError.capacity("localization orbital count") }
        for metric in metrics {
            try metric.requireSymmetric()
            guard metric.rows == n else { throw VivoQMError.invalid("localization metric shape") }
        }
        let indices = groups.flatMap { $0 }
        guard !groups.isEmpty, groups.allSatisfy({ !$0.isEmpty }), Set(indices).count == indices.count,
              indices.allSatisfy({ (0..<n).contains($0) }) else { throw VivoQMError.invalid("localization groups overlap or contain absent orbitals") }
        var matrices = metrics, rotation = try VivoQMMatrix.identity(n)
        func objective(_ values: [VivoQMMatrix]) -> Double {
            values.reduce(0) { sum, matrix in sum + (0..<n).reduce(0) { $0 + matrix[$1,$1] * matrix[$1,$1] } }
        }
        var previous = objective(matrices)
        for sweep in 1...maximumSweeps {
            for group in groups { for (offset, p) in group.enumerated() { for q in group.dropFirst(offset + 1) {
                var aa = 0.0, bb = 0.0, ab = 0.0
                for m in matrices {
                    let x = 0.5 * (m[p,p] - m[q,q]), y = m[p,q]
                    aa += x*x; bb += y*y; ab += x*y
                }
                if hypot(aa-bb, 2*ab) < tolerance * 0.01 { continue }
                let theta = 0.25 * atan2(2*ab, aa-bb), c = cos(theta), s = sin(theta)
                for k in 0..<n {
                    let a = rotation[k,p], b = rotation[k,q]
                    rotation[k,p] = c*a+s*b; rotation[k,q] = -s*a+c*b
                }
                for index in matrices.indices {
                    var m = matrices[index]
                    let a = m[p,p], b = m[p,q], d = m[q,q]
                    for k in 0..<n where k != p && k != q {
                        let x = m[k,p], y = m[k,q]
                        m[k,p] = c*x+s*y; m[p,k] = m[k,p]
                        m[k,q] = -s*x+c*y; m[q,k] = m[k,q]
                    }
                    m[p,p] = c*c*a+2*c*s*b+s*s*d
                    m[q,q] = s*s*a-2*c*s*b+c*c*d
                    m[p,q] = (c*c-s*s)*b+c*s*(d-a); m[q,p] = m[p,q]
                    matrices[index] = m
                }
            } } }
            let value = objective(matrices)
            guard value + 1e-9 * max(1,abs(previous)) >= previous else {
                throw VivoQMError.convergence("localization objective decreased")
            }
            if abs(value-previous) <= tolerance * max(1,abs(value)) {
                return .init(rotation: rotation, objective: value, sweeps: sweep, converged: true)
            }
            previous = value
        }
        return .init(rotation: rotation, objective: previous, sweeps: maximumSweeps, converged: false)
    }

    public static func mullikenPopulationMetrics(coefficients c: VivoQMMatrix,
                                                 overlap s: VivoQMMatrix,
                                                 basisAtomIndices: [Int]) throws -> [VivoQMMatrix] {
        try c.validate(); try s.requireSymmetric()
        guard c.rows == s.rows, basisAtomIndices.count == c.rows,
              basisAtomIndices.allSatisfy({ (0..<256).contains($0) }), c.columns <= 48 else {
            throw VivoQMError.invalid("Mulliken population dimensions")
        }
        let gram = try s.transformed(by: c)
        guard try gram.adding(.identity(c.columns),scale:-1).frobeniusNorm < 1e-8 else {
            throw VivoQMError.invalid("population orbitals are not S-orthonormal")
        }
        let sc = try s.multiplied(by:c)
        return Set(basisAtomIndices).sorted().map { atom in
            var m = VivoQMMatrix(c.columns,c.columns)
            for mu in basisAtomIndices.indices where basisAtomIndices[mu] == atom {
                for p in 0..<c.columns { for q in 0..<c.columns {
                    m[p,q] += 0.5 * (c[mu,p]*sc[mu,q] + sc[mu,p]*c[mu,q])
                } }
            }
            return m
        }
    }

    /// U=exp(-K), K real skew-symmetric. Scaling/squaring in FP64; the
    /// orthogonality check is mandatory before any Hamiltonian transform.
    public static func rotation(dimension n: Int, pairs: [[Int]], parameters: [Double]) throws -> VivoQMMatrix {
        guard (1...48).contains(n), pairs.count == parameters.count,
              pairs.count <= n*(n-1)/2, parameters.allSatisfy({ $0.isFinite && abs($0) <= 100 }),
              pairs.allSatisfy({ $0.count == 2 && $0[0] >= 0 && $0[0] < $0[1] && $0[1] < n }),
              Set(pairs).count == pairs.count else { throw VivoQMError.invalid("skew orbital generator") }
        var generator = VivoQMMatrix(n,n)
        for (pair,value) in zip(pairs,parameters) { generator[pair[0],pair[1]] = -value; generator[pair[1],pair[0]] = value }
        let norm = generator.frobeniusNorm
        let squarings = norm <= 0.5 ? 0 : Int(ceil(log2(norm/0.5)))
        let scaled = try generator.scaled(pow(2,-Double(squarings)))
        var result = try VivoQMMatrix.identity(n), term = result
        var converged = false
        for order in 1...64 {
            term = try term.multiplied(by:scaled).scaled(1/Double(order))
            result = try result.adding(term)
            if term.frobeniusNorm < 1e-16 { converged = true; break }
        }
        guard converged else { throw VivoQMError.convergence("orbital exponential series") }
        for _ in 0..<squarings { result = try result.multiplied(by:result) }
        let gram = try result.transposed.multiplied(by:result)
        guard try gram.adding(.identity(n),scale:-1).frobeniusNorm < 1e-9 else {
            throw VivoQMError.convergence("orbital exponential lost orthogonality")
        }
        return result
    }

    /// Real orthogonal Procrustes transport. overlap[i,j]=<reference_i|current_j>.
    /// A low singular value is an active-space continuity failure, not a reason
    /// to silently relabel unrelated orbitals. Supply cross-geometry AO overlap,
    /// not C_ref^T C_current when the AO basis is nonorthogonal or has moved.
    public static func alignment(overlap: VivoQMMatrix,
                                 minimumSingularValue: Double = 0.5) throws -> VivoOrbitalAlignment {
        try overlap.validate()
        let n = overlap.rows
        guard n == overlap.columns, n <= 48, minimumSingularValue.isFinite,
              (1e-6...1).contains(minimumSingularValue) else { throw VivoQMError.invalid("path overlap transport contract") }
        let gram = try overlap.multiplied(by:overlap.transposed)
        let eigen = try gram.symmetricEigen()
        let singular = eigen.values.map { sqrt(max(0,$0)) }
        guard singular.allSatisfy({ $0 >= minimumSingularValue && $0 <= 1 + 1e-7 }) else {
            throw VivoQMError.invalid("reaction-path orbital subspaces differ beyond the overlap threshold")
        }
        var inverseRoot = VivoQMMatrix(n,n)
        for k in 0..<n { for i in 0..<n { for j in 0..<n {
            inverseRoot[i,j] += eigen.vectors[i,k]*eigen.vectors[j,k]/singular[k]
        } } }
        let rotation = try overlap.transposed.multiplied(by:inverseRoot)
        let check = try rotation.transposed.multiplied(by:rotation)
        guard try check.adding(.identity(n),scale:-1).frobeniusNorm < 1e-8 else {
            throw VivoQMError.convergence("path orbital transport is not orthogonal")
        }
        return .init(rotation:rotation, singularValues:singular)
    }
}

public struct VivoOrbitalAlignment: Codable, Sendable {
    public let rotation: VivoQMMatrix
    public let singularValues: [Double]
}

public struct VivoSchmidtBath: Codable, Sendable {
    public let fragment: [Int]
    public let environment: [Int]
    /// n x retained rank. Zero rank is represented by nil, not a fake orbital.
    public let coefficients: VivoQMMatrix?
    public let singularValues: [Double]
    public let discardedSquaredWeight: Double
}

public enum VivoDMETBathBuilder {
    /// SVD of gamma[E,F]. Its rank cannot exceed |F|, even when gamma is
    /// correlated. Extra correlated bath directions need a different explicit
    /// construction and are never fabricated by changing this rank contract.
    public static func schmidt(density: VivoQMMatrix, fragment: [Int],
                               singularTolerance: Double = 1e-8) throws -> VivoSchmidtBath {
        try density.requireSymmetric()
        let n = density.rows
        guard n <= 48, !fragment.isEmpty, fragment.count < n,
              Set(fragment).count == fragment.count,
              fragment.allSatisfy({ (0..<n).contains($0) }), singularTolerance.isFinite,
              (1e-10...0.1).contains(singularTolerance) else { throw VivoQMError.invalid("Schmidt bath selection or tolerance") }
        let environment = (0..<n).filter { !fragment.contains($0) }
        var cross = VivoQMMatrix(environment.count,fragment.count)
        for (i,a) in environment.enumerated() { for (j,b) in fragment.enumerated() { cross[i,j] = density[a,b] } }
        let eigen = try cross.transposed.multiplied(by:cross).symmetricEigen(tolerance:1e-15)
        let order = eigen.values.indices.reversed()
        let retained = order.filter { eigen.values[$0] > singularTolerance*singularTolerance }
        guard retained.count <= min(environment.count,fragment.count) else { throw VivoQMError.convergence("Schmidt SVD rank exceeds matrix dimensions") }
        let singular = retained.map { sqrt(eigen.values[$0]) }
        var coefficients: VivoQMMatrix? = retained.isEmpty ? nil : VivoQMMatrix(n,retained.count)
        for (column,k) in retained.enumerated() { for (i,a) in environment.enumerated() {
            var value = 0.0
            for j in fragment.indices { value += cross[i,j]*eigen.vectors[j,k] }
            coefficients![a,column] = value/singular[column]
        } }
        if let c = coefficients {
            let gram = try c.transposed.multiplied(by:c)
            guard try gram.adding(.identity(c.columns),scale:-1).frobeniusNorm < 1e-6 else {
                throw VivoQMError.convergence("ill-conditioned Schmidt bath; increase the singular tolerance")
            }
        }
        return .init(fragment:fragment,environment:environment,coefficients:coefficients,
                     singularValues:singular,
                     discardedSquaredWeight:eigen.values.indices.filter { !retained.contains($0) }.reduce(0) { $0+max(0,eigen.values[$1]) })
    }

    /// Single-shot closed-shell determinant DMET construction. The unentangled
    /// environment must have integer occupations. Not self-consistent ECC-DMET.
    public static func closedShellCluster(hamiltonian h: VivoEmbeddedHamiltonian,
                                          spinDensity d: VivoQMMatrix,
                                          fragment: [Int]) throws -> VivoDMETCluster {
        try h.validate(); try d.requireSymmetric()
        let n = h.orbitalCount
        guard h.alphaElectrons == h.betaElectrons, d.rows == n,
              abs((0..<n).reduce(0.0) { $0+d[$1,$1] } - Double(h.alphaElectrons)) < 1e-7,
              try d.multiplied(by:d).adding(d,scale:-1).frobeniusNorm < 1e-7 else {
            throw VivoQMError.invalid("determinant DMET requires a closed-shell idempotent spin density")
        }
        let bath = try schmidt(density:d,fragment:fragment)
        var columns: [[Double]] = fragment.map { index in (0..<n).map { $0 == index ? 1 : 0 } }
        if let b = bath.coefficients { for j in 0..<b.columns { columns.append((0..<n).map { b[$0,j] }) } }
        let clusterCount = columns.count
        for index in bath.environment {
            var vector = (0..<n).map { $0 == index ? 1.0 : 0.0 }
            // Twice-reorthogonalized modified Gram-Schmidt.
            for _ in 0..<2 { for column in columns {
                let dot = zip(vector,column).reduce(0) { $0+$1.0*$1.1 }
                for i in 0..<n { vector[i] -= dot*column[i] }
            } }
            let norm = sqrt(vector.reduce(0) { $0+$1*$1 })
            if norm > 1e-8 { columns.append(vector.map { $0/norm }) }
        }
        guard columns.count == n else { throw VivoQMError.convergence("DMET environment complement rank") }
        var rotation = VivoQMMatrix(n,n)
        for j in 0..<n { for i in 0..<n { rotation[i,j] = columns[j][i] } }
        var occupied: [Int] = []
        if clusterCount < n {
            let remainder = n-clusterCount
            var complement = VivoQMMatrix(n,remainder)
            for j in 0..<remainder { for i in 0..<n { complement[i,j] = rotation[i,j+clusterCount] } }
            let eigen = try d.transformed(by:complement).symmetricEigen()
            guard eigen.values.allSatisfy({ abs($0) < 1e-7 || abs($0-1) < 1e-7 }) else {
                throw VivoQMError.invalid("DMET residual environment is correlated or the bath is incomplete")
            }
            let natural = try complement.multiplied(by:eigen.vectors)
            for j in 0..<remainder {
                if eigen.values[j] > 0.5 { occupied.append(j+clusterCount) }
                for i in 0..<n { rotation[i,j+clusterCount] = natural[i,j] }
            }
        }
        let ids = fragment.map { h.orbitalIDs[$0] } + (fragment.count..<n).map { "dmet-\($0)" }
        let rotated = try h.rotated(by:rotation,orbitalIDs:ids)
        let cluster = try rotated.activeSpace(active:Array(0..<clusterCount),doublyOccupiedCore:occupied)
        return .init(hamiltonian:cluster,rotation:rotation,bath:bath,frozenCoreColumns:occupied)
    }
}

public struct VivoDMETCluster: Codable, Sendable {
    public let hamiltonian: VivoEmbeddedHamiltonian
    public let rotation: VivoQMMatrix
    public let bath: VivoSchmidtBath
    public let frozenCoreColumns: [Int]
}
