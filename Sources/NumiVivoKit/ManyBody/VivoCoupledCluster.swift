import Foundation

/// Small-system projective CCSD, not a variational CI substitute. All arithmetic
/// is FP64. The complete fixed-(Nalpha,Nbeta) sector is an evaluation space for
/// exp(-T) H exp(T); it is never diagonalized by this solver. This deliberately
/// bounded implementation is a conformance implementation, not scalable CCSD.
public struct VivoCCSDConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var residualTolerance: Double
    public var energyToleranceHartree: Double
    public var maximumStepNorm: Double
    public var maximumLineSearchSteps: Int
    public var solveLambda: Bool
    public init(maximumIterations: Int = 50, residualTolerance: Double = 1e-9,
                energyToleranceHartree: Double = 1e-10, maximumStepNorm: Double = 1,
                maximumLineSearchSteps: Int = 16, solveLambda: Bool = true) {
        self.maximumIterations = maximumIterations; self.residualTolerance = residualTolerance
        self.energyToleranceHartree = energyToleranceHartree; self.maximumStepNorm = maximumStepNorm
        self.maximumLineSearchSteps = maximumLineSearchSteps; self.solveLambda = solveLambda
    }
    public func validate() throws {
        guard (1...1000).contains(maximumIterations), (1...64).contains(maximumLineSearchSteps),
              residualTolerance.isFinite, residualTolerance > 0, energyToleranceHartree.isFinite,
              energyToleranceHartree > 0, maximumStepNorm.isFinite, maximumStepNorm > 0 else {
            throw VivoChemistryError.invalid("CCSD iteration, residual, energy or step contract")
        }
    }
}
public struct VivoCCSDExcitation: Codable, Sendable, Equatable {
    /// Interleaved spin-orbital indices. The application order is annihilations
    /// ascending, then creations descending. referencePhase normalizes tau|0>
    /// to the positive canonical determinant. These are not raw tensor t2 slots.
    public let annihilatedModes: [Int]
    public let createdModes: [Int]
    public let referencePhase: Double
    public let determinant: UInt64
}
public enum VivoCCSDTermination: String, Codable, Sendable {
    case converged, iterationLimit, lineSearchFailed, singularJacobian
}
public struct VivoCCSDIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let projectedResidualNorm: Double
    public let acceptedStepScale: Double
}
public struct VivoCCSDLambdaResponse: Codable, Sendable, Equatable {
    public let multipliers: [Double]
    public let residualNorm: Double
    /// Biorthogonal, NOT probability amplitudes. <left|right> = 1.
    public let leftCoefficients: [Double]
    public let rightCoefficients: [Double]
    public let biorthogonalOverlap: Double
}
public struct VivoCCSDResult: Codable, Sendable, Equatable {
    public let configuration: VivoCCSDConfiguration
    public let orbitalIdentifiers: [String]
    public let referenceDeterminant: UInt64
    public let excitations: [VivoCCSDExcitation]
    public let amplitudes: [Double]
    public let referenceEnergyHartree: Double
    public let energyHartree: Double
    public var correlationEnergyHartree: Double { energyHartree - referenceEnergyHartree }
    public let projectedResidualNorm: Double
    /// Normalized exp(T)|0>, for inspection only. Its Hermitian expectation and
    /// RDMs are not CC response properties. Use lambda for CC response RDMs.
    public let normalizedRightState: VivoCIState
    public let rightExpectationEnergyHartree: Double
    public let lambda: VivoCCSDLambdaResponse?
    public let termination: VivoCCSDTermination
    public var converged: Bool { termination == .converged }
    public let iterations: [VivoCCSDIteration]
    public let operatorApplications: Int
    public let energyReference: String
}

private struct VivoCCWork {
    var remaining: Int
    let limit: Int
    init(_ limit: Int) { self.limit = limit; remaining = limit }
    mutating func charge(_ count: Int) throws {
        guard count >= 0, count <= remaining else {
            throw VivoChemistryError.resourceLimit("CCSD aggregate operator-application budget")
        }
        remaining -= count
    }
    var used: Int { limit - remaining }
}
private struct VivoCCLink { let source: Int; let destination: Int; let sign: Double }
private struct VivoCCSector {
    let determinants: [UInt64]
    let referenceIndex: Int
    let excitations: [VivoCCSDExcitation]
    let excitedIndices: [Int]
    let links: [[VivoCCLink]]
    let h: VivoQMMatrix // Scalar energy offset intentionally omitted.
    let nilpotenceOrder: Int

    static func make(_ input: VivoEmbeddedHamiltonian, reference requested: UInt64?,
                     budget: VivoChemistryBudget, work: inout VivoCCWork) throws -> Self {
        try budget.validate(); try input.validate(budget: budget)
        let n = input.orbitalCount
        guard n <= 31 else { throw VivoChemistryError.resourceLimit("CCSD determinant width") }
        func count(_ k: Int) throws -> Int {
            let k = min(k, n-k)
            if k == 0 { return 1 }
            var value = 1
            for i in 1...k {
                let (product, overflow) = value.multipliedReportingOverflow(by: n-k+i)
                guard !overflow else { throw VivoChemistryError.resourceLimit("CCSD sector count overflow") }
                value = product/i
                guard value <= budget.maximumDeterminants else {
                    throw VivoChemistryError.resourceLimit("CCSD full fixed-spin sector exceeds determinant limit")
                }
            }
            return value
        }
        let ca = try count(input.alphaElectrons), cb = try count(input.betaElectrons)
        let (dimension, overflow) = ca.multipliedReportingOverflow(by: cb)
        guard !overflow, dimension <= budget.maximumDeterminants else {
            throw VivoChemistryError.resourceLimit("CCSD full fixed-spin sector exceeds determinant limit")
        }
        // Conservative live storage: H, Jacobian, elimination, sparse links and
        // vectors. Covers both dense workspaces and link allocation before use.
        _ = try budget.elements([dimension, dimension], simultaneousArrays: 16)
        func combinations(_ population: Int, spin: Int) -> [UInt64] {
            var result: [UInt64] = []
            func visit(_ start: Int, _ left: Int, _ bits: UInt64) {
                if left == 0 { result.append(bits); return }
                guard start <= n-left else { return }
                for p in start...(n-left) { visit(p+1, left-1, bits | (UInt64(1) << (2*p+spin))) }
            }
            visit(0, population, 0); return result
        }
        let alpha = combinations(input.alphaElectrons, spin: 0)
        let beta = combinations(input.betaElectrons, spin: 1)
        let determinants = alpha.flatMap { a in beta.map { a | $0 } }.sorted()
        let index = Dictionary(uniqueKeysWithValues: determinants.enumerated().map { ($0.element, $0.offset) })
        var defaultReference: UInt64 = 0
        for p in 0..<input.alphaElectrons { defaultReference |= UInt64(1) << (2*p) }
        for p in 0..<input.betaElectrons { defaultReference |= UInt64(1) << (2*p+1) }
        let reference = requested ?? defaultReference
        guard let referenceIndex = index[reference] else {
            throw VivoChemistryError.invalid("CCSD reference lies outside requested spin sector")
        }
        var h = VivoQMMatrix(dimension, dimension)
        func add(_ coefficient: Double, _ actions: [VivoFermionAction],
                 work: inout VivoCCWork) throws {
            if coefficient == 0 { return }
            try work.charge(dimension)
            for (j, determinant) in determinants.enumerated() {
                if let (target, phase) = vivoApplyFermions(determinant, actions), let i = index[target] {
                    h[i,j] += coefficient*phase
                }
            }
        }
        for p in 0..<n { for q in 0..<n {
            for spin in 0..<2 {
                try add(input.oneElectron[p,q], [.init(mode: 2*q+spin, creation: false),
                                                .init(mode: 2*p+spin, creation: true)], work: &work)
            }
            for r in 0..<n { for s in 0..<n {
                for sigma in 0..<2 { for tau in 0..<2 {
                    try add(0.5*input.eri(p,q,r,s), [.init(mode: 2*q+sigma, creation: false),
                        .init(mode: 2*s+tau, creation: false), .init(mode: 2*r+tau, creation: true),
                        .init(mode: 2*p+sigma, creation: true)], work: &work)
                } }
            } }
        } }
        guard h.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("nonfinite CCSD sector Hamiltonian") }
        for i in 0..<dimension { for j in 0..<i {
            guard abs(h[i,j]-h[j,i]) <= 1e-10*max(1, abs(h[i,j]), abs(h[j,i])) else {
                throw VivoChemistryError.invalid("CCSD sector Hamiltonian is not Hermitian")
            }
            let average = 0.5*(h[i,j]+h[j,i]); h[i,j] = average; h[j,i] = average
        } }
        var excitations: [VivoCCSDExcitation] = [], excited: [Int] = [], links: [[VivoCCLink]] = []
        for (i, determinant) in determinants.enumerated() {
            let rank = (reference ^ determinant).nonzeroBitCount/2
            guard rank == 1 || rank == 2 else { continue }
            let removed = (0..<(2*n)).filter { reference & ~determinant & (UInt64(1) << $0) != 0 }
            let added = (0..<(2*n)).filter { determinant & ~reference & (UInt64(1) << $0) != 0 }
            let actions = removed.map { VivoFermionAction(mode: $0, creation: false) }
                + added.reversed().map { VivoFermionAction(mode: $0, creation: true) }
            guard let (target, phase) = vivoApplyFermions(reference, actions), target == determinant else {
                throw VivoChemistryError.invalid("CCSD excitation convention")
            }
            excitations.append(.init(annihilatedModes: removed, createdModes: added,
                                     referencePhase: phase, determinant: determinant))
            excited.append(i)
            var transitions: [VivoCCLink] = []
            try work.charge(dimension)
            for (j, source) in determinants.enumerated() {
                if let (target, sign) = vivoApplyFermions(source, actions), let k = index[target] {
                    transitions.append(.init(source: j, destination: k, sign: sign*phase))
                }
            }
            links.append(transitions)
        }
        return .init(determinants: determinants, referenceIndex: referenceIndex,
                     excitations: excitations, excitedIndices: excited, links: links, h: h,
                     nilpotenceOrder: input.alphaElectrons+input.betaElectrons)
    }
    func cluster(_ vector: [Double], amplitudes: [Double], transpose: Bool = false,
                 work: inout VivoCCWork) throws -> [Double] {
        var result = [Double](repeating: 0, count: determinants.count)
        for mu in links.indices where amplitudes[mu] != 0 {
            try work.charge(links[mu].count)
            for link in links[mu] {
                let source = transpose ? link.destination : link.source
                let destination = transpose ? link.source : link.destination
                result[destination] += amplitudes[mu]*link.sign*vector[source]
            }
        }
        return result
    }
    func exponential(_ vector: [Double], amplitudes: [Double], sign: Double,
                     transpose: Bool = false, work: inout VivoCCWork) throws -> [Double] {
        var result = vector, term = vector
        if nilpotenceOrder > 0 {
            for order in 1...nilpotenceOrder {
                term = try cluster(term, amplitudes: amplitudes, transpose: transpose, work: &work)
                for i in term.indices { term[i] *= sign/Double(order); result[i] += term[i] }
            }
        }
        guard result.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("CCSD exponential overflow") }
        return result
    }
    func hTimes(_ vector: [Double], work: inout VivoCCWork) throws -> [Double] {
        let d = determinants.count
        try work.charge(d*d)
        let column = try VivoQMMatrix(rows: d, columns: 1, values: vector)
        return try h.multiplied(by: column).values
    }
    func similarity(_ vector: [Double], amplitudes: [Double], work: inout VivoCCWork) throws -> [Double] {
        let right = try exponential(vector, amplitudes: amplitudes, sign: 1, work: &work)
        let applied = try hTimes(right, work: &work)
        return try exponential(applied, amplitudes: amplitudes, sign: -1, work: &work)
    }
    /// Excitation operators commute. d(Hbar)/dt_mu = [Hbar,tau_mu].
    /// This Jacobian is analytic; no energy finite difference is used for CCSD.
    func derivatives(_ amplitudes: [Double], referenceColumn: [Double],
                     work: inout VivoCCWork) throws -> (jacobian: VivoQMMatrix, energy: [Double]) {
        let m = links.count, d = determinants.count
        var jacobian = VivoQMMatrix(m,m), energy = [Double](repeating: 0, count: m)
        for mu in 0..<m {
            var ket = [Double](repeating: 0, count: d); ket[excitedIndices[mu]] = 1
            var column = try similarity(ket, amplitudes: amplitudes, work: &work)
            try work.charge(links[mu].count)
            for link in links[mu] { column[link.destination] -= link.sign*referenceColumn[link.source] }
            energy[mu] = column[referenceIndex]
            for nu in 0..<m { jacobian[nu,mu] = column[excitedIndices[nu]] }
        }
        guard jacobian.values.allSatisfy(\.isFinite), energy.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.convergence("nonfinite CCSD Jacobian")
        }
        return (jacobian,energy)
    }
}

public enum VivoCoupledCluster {
    public static func solve(_ h: VivoEmbeddedHamiltonian,
                             configuration: VivoCCSDConfiguration = .init(),
                             referenceDeterminant: UInt64? = nil, initialAmplitudes: [Double]? = nil,
                             budget: VivoChemistryBudget = .init()) throws -> VivoCCSDResult {
        try configuration.validate(); try budget.validate()
        var work = VivoCCWork(budget.maximumOperatorApplications)
        let sector = try VivoCCSector.make(h, reference: referenceDeterminant, budget: budget, work: &work)
        let m = sector.excitations.count, d = sector.determinants.count, ref = sector.referenceIndex
        var t = initialAmplitudes ?? [Double](repeating: 0, count: m)
        guard t.count == m, t.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("CCSD restart amplitudes") }
        var reference = [Double](repeating: 0, count: d); reference[ref] = 1
        var column = try sector.similarity(reference, amplitudes: t, work: &work)
        func norm(_ values: [Double]) -> Double { values.reduce(0.0) { hypot($0,$1) } }
        var residual = sector.excitedIndices.map { column[$0] }
        var previousEnergy: Double?, history: [VivoCCSDIteration] = []
        var termination = VivoCCSDTermination.iterationLimit, stepScale = 0.0
        for iteration in 0...configuration.maximumIterations {
            let energy = column[ref]+h.constantEnergyHartree, error = norm(residual)
            history.append(.init(iteration: iteration, energyHartree: energy,
                                 projectedResidualNorm: error, acceptedStepScale: stepScale))
            let energySettled = previousEnergy.map { abs(energy-$0) <= configuration.energyToleranceHartree } ?? true
            if error <= configuration.residualTolerance {
                if energySettled { termination = .converged; break }
                // A stationary confirmation uses unchanged amplitudes; it does
                // not force a spurious Newton step on an already solved root.
                previousEnergy = energy; continue
            }
            if iteration == configuration.maximumIterations { break }
            let derivatives = try sector.derivatives(t, referenceColumn: column, work: &work)
            let step: [Double]
            do { step = try VivoQMDenseAlgebra.solve(derivatives.jacobian, rhs: residual.map { -$0 }) }
            catch VivoChemistryError.convergence(_) { termination = .singularJacobian; break }
            let length = norm(step), bound = min(1, configuration.maximumStepNorm/max(length, 1e-300))
            var scale = bound, accepted = false
            for _ in 0..<configuration.maximumLineSearchSteps {
                let candidate = zip(t,step).map { $0+scale*$1 }
                let trial = try sector.similarity(reference, amplitudes: candidate, work: &work)
                let trialResidual = sector.excitedIndices.map { trial[$0] }, trialNorm = norm(trialResidual)
                if trialNorm < error*(1-1e-4*scale) || trialNorm <= configuration.residualTolerance {
                    previousEnergy = energy; t = candidate; column = trial; residual = trialResidual
                    stepScale = scale; accepted = true; break
                }
                scale *= 0.5
            }
            if !accepted { termination = .lineSearchFailed; break }
        }
        let right = try sector.exponential(reference, amplitudes: t, sign: 1, work: &work)
        let rightNorm = norm(right)
        guard rightNorm.isFinite, rightNorm > 0 else { throw VivoChemistryError.convergence("CCSD right-state norm") }
        let normalized = right.map { $0/rightNorm }
        let state = VivoCIState(orbitalCount: h.orbitalCount, alphaElectrons: h.alphaElectrons,
                               betaElectrons: h.betaElectrons, determinants: sector.determinants, coefficients: normalized)
        try state.validate(budget: budget)
        let hRight = try sector.hTimes(normalized, work: &work)
        let expectation = zip(normalized,hRight).reduce(h.constantEnergyHartree) { $0+$1.0*$1.1 }
        var response: VivoCCSDLambdaResponse?
        if termination == .converged, configuration.solveLambda {
            let derivatives = try sector.derivatives(t, referenceColumn: column, work: &work)
            let multipliers = m == 0 ? [] : try VivoQMDenseAlgebra.solve(derivatives.jacobian.transposed,
                                                                       rhs: derivatives.energy.map { -$0 })
            var lambdaResidual = derivatives.energy
            for j in 0..<m { for i in 0..<m { lambdaResidual[j] += derivatives.jacobian[i,j]*multipliers[i] } }
            let responseNorm = norm(lambdaResidual)
            guard responseNorm.isFinite, responseNorm <= max(1e-9, 10*configuration.residualTolerance) else {
                throw VivoChemistryError.convergence("CCSD lambda residual")
            }
            var bra = reference
            for mu in 0..<m { bra[sector.excitedIndices[mu]] = multipliers[mu] }
            let left = try sector.exponential(bra, amplitudes: t, sign: -1, transpose: true, work: &work)
            let overlap = zip(left,right).reduce(0.0) { $0+$1.0*$1.1 }
            guard overlap.isFinite, abs(overlap-1) <= 1e-9 else {
                throw VivoChemistryError.convergence("CCSD biorthogonal normalization")
            }
            response = .init(multipliers: multipliers, residualNorm: responseNorm, leftCoefficients: left,
                             rightCoefficients: right, biorthogonalOverlap: overlap)
        }
        return .init(configuration: configuration, orbitalIdentifiers: h.orbitalIdentifiers,
                     referenceDeterminant: sector.determinants[ref], excitations: sector.excitations, amplitudes: t,
                     referenceEnergyHartree: sector.h[ref,ref]+h.constantEnergyHartree,
                     energyHartree: column[ref]+h.constantEnergyHartree, projectedResidualNorm: norm(residual),
                     normalizedRightState: state, rightExpectationEnergyHartree: expectation,
                     lambda: response, termination: termination, iterations: history,
                     operatorApplications: work.used, energyReference: h.energyReference)
    }

    /// Unrelaxed CCSD response RDMs from <0|(1+Lambda)e^-T O e^T|0>.
    /// These can be non-Hermitian/nonpositive. They MUST NOT be used as physical
    /// density operators for orbital entropy or QIO mutual-information scoring.
    public static func responseDensityMatrices(_ result: VivoCCSDResult,
                                               budget: VivoChemistryBudget = .init()) throws -> VivoSpinRDMs {
        try budget.validate(); try result.configuration.validate()
        try result.normalizedRightState.validate(budget: budget)
        guard result.converged, result.projectedResidualNorm.isFinite,
              result.projectedResidualNorm <= result.configuration.residualTolerance,
              let response = result.lambda else {
            throw VivoChemistryError.invalid("CC response requires converged CCSD and lambda")
        }
        let state = result.normalizedRightState, d = state.determinants.count, m = 2*state.orbitalCount
        guard response.leftCoefficients.count == d, response.rightCoefficients.count == d,
              response.leftCoefficients.allSatisfy(\.isFinite), response.rightCoefficients.allSatisfy(\.isFinite),
              response.residualNorm.isFinite, response.residualNorm <= max(1e-9, 10*result.configuration.residualTolerance) else {
            throw VivoChemistryError.invalid("CC response vector dimensions or residual")
        }
        let overlap = zip(response.leftCoefficients,response.rightCoefficients).reduce(0.0) { $0+$1.0*$1.1 }
        guard overlap.isFinite, abs(overlap-1) < 1e-9 else { throw VivoChemistryError.invalid("CC response overlap") }
        let count = try budget.elements([m,m,m,m], simultaneousArrays: 3)
        var work = VivoCCWork(budget.maximumOperatorApplications)
        let index = Dictionary(uniqueKeysWithValues: state.determinants.enumerated().map { ($0.element,$0.offset) })
        func expectation(_ actions: [VivoFermionAction], work: inout VivoCCWork) throws -> Double {
            try work.charge(d)
            var value = 0.0
            for (j, determinant) in state.determinants.enumerated() {
                if let (target, sign) = vivoApplyFermions(determinant,actions), let i = index[target] {
                    value += response.leftCoefficients[i]*response.rightCoefficients[j]*sign
                }
            }
            return value
        }
        var one = VivoQMMatrix(m,m), two = [Double](repeating: 0, count: count)
        for p in 0..<m { for q in 0..<m {
            one[p,q] = try expectation([.init(mode: q, creation: false), .init(mode: p, creation: true)], work: &work)
            for r in 0..<m { for s in 0..<m where p != q && r != s {
                two[((p*m+q)*m+r)*m+s] = try expectation([.init(mode: r, creation: false), .init(mode: s, creation: false),
                    .init(mode: q, creation: true), .init(mode: p, creation: true)], work: &work)
            } }
        } }
        let electrons = state.alphaElectrons+state.betaElectrons
        guard one.values.allSatisfy(\.isFinite), two.allSatisfy(\.isFinite),
              abs((0..<m).reduce(0.0) { $0+one[$1,$1] }-Double(electrons)) < 1e-8 else {
            throw VivoChemistryError.convergence("CC response RDM trace or finite values")
        }
        let rdm = VivoSpinRDMs(spinOrbitalCount: m, electronCount: electrons, one: one, two: two)
        for p in 0..<m { for q in 0..<m {
            let trace = (0..<m).reduce(0.0) { $0+rdm.gamma2(p,$1,q,$1) }
            guard abs(trace-Double(electrons-1)*one[p,q]) < 1e-8 else {
                throw VivoChemistryError.convergence("CC response 2-RDM contraction")
            }
        } }
        return rdm
    }
}
