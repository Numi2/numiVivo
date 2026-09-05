import Foundation

public struct VivoActiveSpace: Codable, Sendable, Equatable {
    public let doublyOccupiedCore: [Int]
    public let active: [Int]
    /// Orbital columns excluded from every optimization rotation.
    public let frozenOrbitals: [Int]
    public init(doublyOccupiedCore: [Int] = [], active: [Int], frozenOrbitals: [Int] = []) {
        self.doublyOccupiedCore = doublyOccupiedCore; self.active = active; self.frozenOrbitals = frozenOrbitals
    }
    public func validate(for h: VivoEmbeddedHamiltonian, budget: VivoChemistryBudget = .init()) throws {
        try h.validate(budget: budget)
        let selected = doublyOccupiedCore+active
        guard !active.isEmpty, Set(selected).count == selected.count,
              selected.allSatisfy({ $0 >= 0 && $0 < h.orbitalCount }),
              Set(frozenOrbitals).count == frozenOrbitals.count,
              frozenOrbitals.allSatisfy({ $0 >= 0 && $0 < h.orbitalCount }),
              doublyOccupiedCore.count <= min(h.alphaElectrons,h.betaElectrons),
              h.alphaElectrons-doublyOccupiedCore.count <= active.count,
              h.betaElectrons-doublyOccupiedCore.count <= active.count else {
            throw VivoChemistryError.invalid("CASSCF core/active/frozen partition or spin population")
        }
    }
}
public struct VivoOrbitalRotationPair: Codable, Sendable, Equatable {
    public let first: Int
    public let second: Int
}
public struct VivoCASSCFConfiguration: Codable, Sendable, Equatable {
    public var maximumMacroIterations: Int
    public var maximumEnergyEvaluations: Int
    public var gradientTolerance: Double
    public var energyToleranceHartree: Double
    public var maximumRotationNorm: Double
    public var maximumLineSearchSteps: Int
    public var historySize: Int
    public var minimumStateGapHartree: Double
    public init(maximumMacroIterations: Int = 100, maximumEnergyEvaluations: Int = 1000,
                gradientTolerance: Double = 1e-6, energyToleranceHartree: Double = 1e-10,
                maximumRotationNorm: Double = 0.2, maximumLineSearchSteps: Int = 16,
                historySize: Int = 7, minimumStateGapHartree: Double = 1e-9) {
        self.maximumMacroIterations = maximumMacroIterations; self.maximumEnergyEvaluations = maximumEnergyEvaluations
        self.gradientTolerance = gradientTolerance; self.energyToleranceHartree = energyToleranceHartree
        self.maximumRotationNorm = maximumRotationNorm; self.maximumLineSearchSteps = maximumLineSearchSteps
        self.historySize = historySize; self.minimumStateGapHartree = minimumStateGapHartree
    }
    public func validate() throws {
        guard (1...1000).contains(maximumMacroIterations), (1...100_000).contains(maximumEnergyEvaluations),
              (1...64).contains(maximumLineSearchSteps), (0...32).contains(historySize),
              gradientTolerance.isFinite, gradientTolerance > 0, energyToleranceHartree.isFinite,
              energyToleranceHartree > 0, maximumRotationNorm.isFinite, maximumRotationNorm > 0,
              maximumRotationNorm <= 1, minimumStateGapHartree.isFinite, minimumStateGapHartree >= 0 else {
            throw VivoChemistryError.invalid("CASSCF optimization settings")
        }
    }
}
public struct VivoCASSCFIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let orbitalGradientNorm: Double
    public let acceptedRotationNorm: Double
    public let activeCIGapHartree: Double?
}
public enum VivoCASSCFTermination: String, Codable, Sendable { case converged, iterationLimit, lineSearchFailed, evaluationLimit }
public struct VivoCASSCFResult: Codable, Sendable, Equable {
    public let partition: VivoActiveSpace
    /// Columns expressed in the original orthonormal spatial-orbital basis.
    public let orbitalRotation: VivoQMMatrix
    public let activeHamiltonian: VivoEmbeddedHamiltonian
    public let activeCI: VivoCIResult
    public let rotationPairs: [VivoOrbitalRotationPair]
    public let orbitalGradient: [Double]
    public let termination: VivoCASSCFTermination
    public var converged: Bool { termination == .converged }
    public let energyEvaluations: Int
    public let iterations: [VivoCASSCFIteration]
}

/// Analytic orbital derivatives of a real CAS state. Full-space one/two-particle
/// densities are reconstructed from the active RDM cumulant plus the doubly
/// occupied core. Empty external orbitals carry no density. No large full-space
/// CI expansion is constructed. This is an electronic, not nuclear, gradient.
public enum VivoCASOrbitalGradient {
    public static func compute(hamiltonian h: VivoEmbeddedHamiltonian, partition: VivoActiveSpace,
                               state: VivoCIState, pairs: [VivoOrbitalRotationPair],
                               budget: VivoChemistryBudget = .init()) throws -> [Double] {
        try partition.validate(for: h, budget: budget); try state.validate(budget: budget)
        let n = h.orbitalCount, a = partition.active.count, core = partition.doublyOccupiedCore.count
        guard state.orbitalCount == a, state.alphaElectrons == h.alphaElectrons-core,
              state.betaElectrons == h.betaElectrons-core,
              pairs.allSatisfy({ $0.first >= 0 && $0.second > $0.first && $0.second < n }) else {
            throw VivoChemistryError.invalid("CAS orbital derivative state/partition")
        }
        let count = try budget.elements([n,n,n,n], simultaneousArrays: 8)
        let (work, overflow) = count.multipliedReportingOverflow(by: n)
        guard !overflow, work <= budget.maximumOperatorApplications/4 else {
            throw VivoChemistryError.resourceLimit("CAS orbital-gradient contraction work")
        }
        let rdm = try VivoCIDensityMatrices.compute(state, budget: budget)
        var alpha = VivoQMMatrix(n,n), beta = VivoQMMatrix(n,n)
        for i in partition.doublyOccupiedCore { alpha[i,i] = 1; beta[i,i] = 1 }
        for p in 0..<a { for q in 0..<a {
            alpha[partition.active[p],partition.active[q]] = rdm.one[2*p,2*q]
            beta[partition.active[p],partition.active[q]] = rdm.one[2*p+1,2*q+1]
        } }
        let one = try alpha.adding(beta)
        var local = [Int](repeating: -1, count: n)
        for (i,p) in partition.active.enumerated() { local[p] = i }
        var two = [Double](repeating: 0, count: count)
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            var value = one[p,q]*one[r,s]-alpha[p,s]*alpha[r,q]-beta[p,s]*beta[r,q]
            if local[p] >= 0, local[q] >= 0, local[r] >= 0, local[s] >= 0 {
                for sigma in 0..<2 { for tau in 0..<2 {
                    value += rdm.cumulant(2*local[p]+sigma,2*local[r]+tau,2*local[q]+sigma,2*local[s]+tau)
                } }
            }
            two[((p*n+q)*n+r)*n+s] = value
        } } } }
        // Symmetrization in each ERI pair is exact in energy and its orbital
        // derivative because real chemist ERIs retain these symmetries under U.
        func gamma(_ p: Int,_ q: Int,_ r: Int,_ s: Int) -> Double { two[((p*n+q)*n+r)*n+s] }
        var symmetric = VivoQMMatrix(n,n*n*n)
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            symmetric[p,(q*n+r)*n+s] = 0.25*(gamma(p,q,r,s)+gamma(q,p,r,s)+gamma(p,q,s,r)+gamma(q,p,s,r))
        } } } }
        let integrals = try VivoQMMatrix(rows: n, columns: n*n*n, values: h.twoElectron)
        // Both contractions use the existing FP64 GEMM path (Accelerate on Apple).
        let g1 = try h.oneElectron.multiplied(by: one.transposed).scaled(2)
        let g2 = try integrals.multiplied(by: symmetric.transposed).scaled(2)
        let derivative = try g1.adding(g2)
        let gradient = pairs.map { derivative[$0.second,$0.first]-derivative[$0.first,$0.second] }
        guard gradient.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite CAS orbital gradient") }
        return gradient
    }
}

public enum VivoCASSCF {
    private struct Correction { var s: [Double]; var y: [Double] }
    private static func dot(_ a: [Double], _ b: [Double]) -> Double { zip(a,b).reduce(0) { $0+$1.0*$1.1 } }
    private static func norm(_ a: [Double]) -> Double { a.reduce(0) { hypot($0,$1) } }
    private static func generator(_ vector: [Double], pairs: [VivoOrbitalRotationPair], n: Int) -> VivoQMMatrix {
        var k = VivoQMMatrix(n,n)
        for (i,pair) in pairs.enumerated() { k[pair.first,pair.second] = vector[i]; k[pair.second,pair.first] = -vector[i] }
        return k
    }
    private static func transport(_ vector: [Double], by u: VivoQMMatrix,
                                  pairs: [VivoOrbitalRotationPair]) throws -> [Double] {
        let k = generator(vector, pairs: pairs, n: u.rows)
        let rotated = try k.congruence(u)
        return pairs.map { rotated[$0.first,$0.second] }
    }
    /// Ground-state, fixed-(Nalpha,Nbeta) CAS orbital optimization. Does not
    /// assert total-spin purity or state-averaged/root-followed CASSCF. Redundant
    /// within-core, within-active and within-external rotations are excluded.
    public static func solve(_ input: VivoEmbeddedHamiltonian, partition: VivoActiveSpace,
                             configuration: VivoCASSCFConfiguration = .init(),
                             initialRotation: VivoQMMatrix? = nil,
                             budget: VivoChemistryBudget = .init()) throws -> VivoCASSCFResult {
        try configuration.validate(); try partition.validate(for: input, budget: budget)
        let n = input.orbitalCount, frozen = Set(partition.frozenOrbitals)
        var group = [Int](repeating: 2, count: n)
        for p in partition.doublyOccupiedCore { group[p] = 0 }
        for p in partition.active { group[p] = 1 }
        var pairs: [VivoOrbitalRotationPair] = []
        for p in 0..<n { for q in (p+1)..<n where group[p] != group[q] && !frozen.contains(p) && !frozen.contains(q) {
            pairs.append(.init(first: p, second: q))
        } }
        var rotation = try initialRotation ?? VivoQMMatrix.identity(n)
        if let initialRotation {
            guard initialRotation.rows == n, initialRotation.columns == n else { throw VivoChemistryError.invalid("CASSCF initial rotation") }
            for p in frozen { for q in 0..<n where abs(initialRotation[q,p]-(p == q ? 1 : 0)) > 1e-12 {
                throw VivoChemistryError.invalid("CASSCF restart rotates an explicitly frozen orbital")
            } }
        }
        var full = try input.rotated(by: rotation, budget: budget), evaluations = 0
        func evaluate(_ h: VivoEmbeddedHamiltonian) throws -> (VivoEmbeddedHamiltonian,VivoCIResult) {
            guard evaluations < configuration.maximumEnergyEvaluations else { throw VivoChemistryError.resourceLimit("CASSCF energy evaluations") }
            evaluations += 1
            let active = try h.frozenCore(active: partition.active, doublyOccupiedCore: partition.doublyOccupiedCore, budget: budget)
            let ci = try VivoConfigurationInteraction.solve(active, budget: budget)
            return (active,ci)
        }
        var (active,ci) = try evaluate(full)
        var gradient: [Double] = [], corrections: [Correction] = [], history: [VivoCASSCFIteration] = []
        var previousEnergy: Double?, stepNorm = 0.0, termination = VivoCASSCFTermination.iterationLimit
        for iteration in 0...configuration.maximumMacroIterations {
            if let gap = ci.nextStateGapHartree, gap < configuration.minimumStateGapHartree, !pairs.isEmpty {
                throw VivoChemistryError.convergence("CAS ground-root degeneracy requires state averaging or an explicit root-following policy")
            }
            gradient = pairs.isEmpty ? [] : try VivoCASOrbitalGradient.compute(hamiltonian: full, partition: partition,
                state: ci.state, pairs: pairs, budget: budget)
            let error = norm(gradient)
            history.append(.init(iteration: iteration, energyHartree: ci.energyHartree, orbitalGradientNorm: error,
                                 acceptedRotationNorm: stepNorm, activeCIGapHartree: ci.nextStateGapHartree))
            let settled = previousEnergy.map { abs(ci.energyHartree-$0) <= configuration.energyToleranceHartree } ?? true
            if error <= configuration.gradientTolerance {
                if settled { termination = .converged; break }
                previousEnergy = ci.energyHartree; continue
            }
            if iteration == configuration.maximumMacroIterations { break }
            if evaluations == configuration.maximumEnergyEvaluations { termination = .evaluationLimit; break }
            var direction = gradient, weights = [Double](repeating: 0, count: corrections.count)
            for i in corrections.indices.reversed() {
                weights[i] = dot(corrections[i].s,direction)/dot(corrections[i].s,corrections[i].y)
                for j in direction.indices { direction[j] -= weights[i]*corrections[i].y[j] }
            }
            if let last = corrections.last {
                let scale = dot(last.s,last.y)/dot(last.y,last.y)
                direction = direction.map { $0*scale }
            }
            for i in corrections.indices {
                let beta = dot(corrections[i].y,direction)/dot(corrections[i].s,corrections[i].y)
                for j in direction.indices { direction[j] += corrections[i].s[j]*(weights[i]-beta) }
            }
            direction = direction.map { -$0 }
            if !direction.allSatisfy(\.isFinite) || dot(direction,gradient) >= 0 {
                corrections.removeAll(); direction = gradient.map { -$0 }
            }
            let bound = min(1, configuration.maximumRotationNorm/max(norm(direction),1e-300))
            direction = direction.map { $0*bound }
            let slope = dot(gradient,direction)
            var scale = 1.0, accepted = false
            for _ in 0..<configuration.maximumLineSearchSteps {
                if evaluations == configuration.maximumEnergyEvaluations { break }
                let step = direction.map { scale*$0 }
                let u = try VivoQMDenseAlgebra.orbitalRotation(generator: generator(step,pairs:pairs,n:n))
                let candidateRotation = try rotation.multiplied(by:u)
                // Always transform original integrals to avoid accumulating
                // four-index roundoff across accepted and rejected trials.
                let candidateFull = try input.rotated(by:candidateRotation,budget:budget)
                let (candidateActive,candidateCI) = try evaluate(candidateFull)
                if candidateCI.energyHartree <= ci.energyHartree+1e-4*scale*slope {
                    let nextGradient = try VivoCASOrbitalGradient.compute(hamiltonian:candidateFull,partition:partition,
                        state:candidateCI.state,pairs:pairs,budget:budget)
                    let oldGradient = try transport(gradient,by:u,pairs:pairs)
                    var transported: [Correction] = []
                    for correction in corrections {
                        let s = try transport(correction.s,by:u,pairs:pairs), y = try transport(correction.y,by:u,pairs:pairs)
                        if dot(s,y) > 1e-12*norm(s)*norm(y) { transported.append(.init(s:s,y:y)) }
                    }
                    let s = try transport(step,by:u,pairs:pairs), y = zip(nextGradient,oldGradient).map { $0-$1 }
                    if dot(s,y) > 1e-12*norm(s)*norm(y), configuration.historySize > 0 { transported.append(.init(s:s,y:y)) }
                    corrections = Array(transported.suffix(configuration.historySize))
                    previousEnergy = ci.energyHartree; stepNorm = norm(step)
                    rotation = candidateRotation; full = candidateFull; active = candidateActive; ci = candidateCI
                    accepted = true; break
                }
                scale *= 0.5
            }
            if !accepted {
                termination = evaluations == configuration.maximumEnergyEvaluations ? .evaluationLimit : .lineSearchFailed
                break
            }
        }
        return .init(partition: partition, orbitalRotation: rotation, activeHamiltonian: active, activeCI: ci,
                     rotationPairs: pairs, orbitalGradient: gradient, termination: termination,
                     energyEvaluations: evaluations, iterations: history)
    }
}
