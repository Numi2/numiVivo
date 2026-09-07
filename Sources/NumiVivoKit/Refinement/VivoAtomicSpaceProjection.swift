import Foundation

public struct VivoAtomicSpaceProjectionConfiguration: Codable, Sendable, Equatable {
    public let minimumActiveWeight: Double
    public let minimumTargetOverlapEigenvalue: Double
    public let degeneracyTolerance: Double
    public init(minimumActiveWeight: Double = 0.15, minimumTargetOverlapEigenvalue: Double = 1e-9,
                degeneracyTolerance: Double = 1e-8) {
        self.minimumActiveWeight = minimumActiveWeight
        self.minimumTargetOverlapEigenvalue = minimumTargetOverlapEigenvalue
        self.degeneracyTolerance = degeneracyTolerance
    }
    public func validate() throws {
        guard minimumActiveWeight.isFinite, minimumActiveWeight > 0, minimumActiveWeight <= 1,
              minimumTargetOverlapEigenvalue.isFinite, minimumTargetOverlapEigenvalue > 0,
              minimumTargetOverlapEigenvalue < 1, degeneracyTolerance.isFinite,
              degeneracyTolerance > 0, degeneracyTolerance < minimumActiveWeight else {
            throw VivoChemistryError.invalid("atomic-space projection thresholds")
        }
    }
}
public struct VivoAtomicSpaceProjectionResult: Codable, Sendable, Equatable {
    public let rotation: VivoQMMatrix
    public let weights: [Double]
    public let occupiedOrbitals: Int
    public let targetAOIndices: [Int]
    public let retainedTargetRank: Int
    public let discardedTargetEigenvalues: [Double]
    public let groups: [[Int]]
    public let initialSpace: VivoActiveSpace
    public let candidateBlocks: [VivoOrbitalRefinementBlock]
    public let interpretation: String
}
/// Occupied/virtual atomic-subspace projection, using the physical AO metric.
/// No principal-shell or atomic-valence character is guessed from exponents.
/// All selected atomic AOs are the conservative default; explicit shell targets
/// can define a smaller valence space. Projection scores are not correlation.
public enum VivoAtomicSpaceProjection {
    public static func select(overlap s: VivoQMMatrix, coefficients c: VivoQMMatrix,
                              occupied: Int, targetAOIndices targets: [Int],
                              configuration cfg: VivoAtomicSpaceProjectionConfiguration = .init(),
                              budget: VivoChemistryBudget = .init()) throws -> VivoAtomicSpaceProjectionResult {
        try cfg.validate(); try budget.validate()
        let n = s.rows
        guard (1...31).contains(n), s.columns == n, c.rows == n, c.columns == n,
              (0...n).contains(occupied), !targets.isEmpty, targets == targets.sorted(),
              Set(targets).count == targets.count, targets.allSatisfy({ (0..<n).contains($0) }),
              s.values.allSatisfy(\.isFinite), c.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("atomic-space AO metric, coefficients, occupation or targets")
        }
        _ = try budget.elements([n,n],simultaneousArrays: 20)
        for i in 0..<n { for j in 0..<i where abs(s[i,j]-s[j,i]) > 1e-10 {
            throw VivoChemistryError.invalid("asymmetric AO overlap")
        } }
        let metric = try VivoQMDenseAlgebra.symmetricEigen(s)
        guard metric.values.allSatisfy({ $0 > 0 }),
              try s.congruence(c).adding(.identity(n),scale: -1).frobeniusNorm < 1e-8 else {
            throw VivoChemistryError.invalid("nonorthonormal molecular orbitals or nonpositive AO metric")
        }
        let m = targets.count
        var st = VivoQMMatrix(m,m), targetCross = VivoQMMatrix(m,n)
        let sc = try s.multiplied(by: c)
        for (i,p) in targets.enumerated() {
            for (j,q) in targets.enumerated() { st[i,j] = s[p,q] }
            for j in 0..<n { targetCross[i,j] = sc[p,j] }
        }
        let eigen = try VivoQMDenseAlgebra.symmetricEigen(st,tolerance: 1e-14)
        let cutoff = cfg.minimumTargetOverlapEigenvalue * max(1,eigen.values.last ?? 1)
        let retained = eigen.values.indices.filter { eigen.values[$0] >= cutoff }
        guard !retained.isEmpty else { throw VivoChemistryError.convergence("atomic target subspace has no retained metric rank") }
        var inverse = VivoQMMatrix(m,m)
        for i in 0..<m { for j in 0..<m { for k in retained {
            inverse[i,j] += eigen.vectors[i,k]*eigen.vectors[j,k]/eigen.values[k]
        } } }
        let projected = try targetCross.transposed.multiplied(by: inverse).multiplied(by: targetCross)
        var rotation = try VivoQMMatrix.identity(n), weights = [Double](repeating: 0,count: n), groups: [[Int]] = []
        for range in [Array(0..<occupied),Array(occupied..<n)] where !range.isEmpty {
            var block = VivoQMMatrix(range.count,range.count)
            for (i,p) in range.enumerated() { for (j,q) in range.enumerated() { block[i,j] = projected[p,q] } }
            let decomposition = try VivoQMDenseAlgebra.symmetricEigen(block,tolerance: 1e-14)
            guard decomposition.values.allSatisfy({ $0 >= -1e-8 && $0 <= 1+1e-8 }) else {
                throw VivoChemistryError.convergence("atomic projection weights outside [0,1]")
            }
            let order = Array(decomposition.values.indices.reversed())
            for (j,k) in order.enumerated() {
                weights[range[j]] = max(0,min(1,decomposition.values[k]))
                // Fix a deterministic sign without making a degenerate eigenvector
                // identity claim. Complete near-degenerate groups are indivisible.
                let pivot = range.indices.max { abs(decomposition.vectors[$0,k]) < abs(decomposition.vectors[$1,k]) }!
                let sign = decomposition.vectors[pivot,k] < 0 ? -1.0 : 1.0
                for (i,p) in range.enumerated() { rotation[p,range[j]] = sign*decomposition.vectors[i,k] }
            }
            var group = [range[0]]
            for p in range.dropFirst() {
                if abs(weights[p]-weights[group[0]]) <= cfg.degeneracyTolerance { group.append(p) }
                else { groups.append(group); group = [p] }
            }
            groups.append(group)
        }
        var selected = Set(groups.filter { group in group.contains { weights[$0] >= cfg.minimumActiveWeight } }.flatMap { $0 })
        // Preserve an occupied and an empty response subspace where both exist.
        // This is a conservative seed, not evidence that the seed is sufficient.
        for range in [Array(0..<occupied),Array(occupied..<n)] where !range.isEmpty && selected.isDisjoint(with: range) {
            if let group = groups.first(where: { $0.contains(range[0]) }) { selected.formUnion(group) }
        }
        let active = selected.sorted()
        guard !active.isEmpty else { throw VivoChemistryError.invalid("empty atomic-space seed") }
        let core = (0..<occupied).filter { !selected.contains($0) }
        let candidates = groups.filter { Set($0).isDisjoint(with: selected) }.map { group in
            VivoOrbitalRefinementBlock(identifier: "atomic-response-"+group.map(String.init).joined(separator: "-"),orbitals: group)
        }
        return .init(rotation: rotation,weights: weights,occupiedOrbitals: occupied,targetAOIndices: targets,
            retainedTargetRank: retained.count,discardedTargetEigenvalues: eigen.values.filter { $0 < cutoff },
            groups: groups,initialSpace: .init(doublyOccupiedCore: core,active: active),candidateBlocks: candidates,
            interpretation: "AO-metric occupied/virtual atomic projection; all nonseed orbitals remain in the declared candidate universe; projection is not entanglement or a chemical-accuracy certificate")
    }
}
