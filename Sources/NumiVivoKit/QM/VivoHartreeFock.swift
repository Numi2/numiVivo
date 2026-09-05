import Foundation

public enum VivoHartreeFockMethod: String, Codable, Sendable { case rhf, uhf }

public struct VivoSCFConfiguration: Codable, Sendable, Equatable {
    public var method: VivoHartreeFockMethod
    public var maximumIterations: Int
    public var energyToleranceHartree: Double
    public var densityTolerance: Double
    public var residualTolerance: Double
    public var overlapEigenvalueFloor: Double
    public var diisHistory: Int
    public var initialDamping: Double
    public init(method: VivoHartreeFockMethod = .rhf, maximumIterations: Int = 128,
                energyToleranceHartree: Double = 1e-10, densityTolerance: Double = 1e-8,
                residualTolerance: Double = 1e-8, overlapEigenvalueFloor: Double = 1e-9,
                diisHistory: Int = 8, initialDamping: Double = 0.2) {
        self.method = method; self.maximumIterations = maximumIterations
        self.energyToleranceHartree = energyToleranceHartree; self.densityTolerance = densityTolerance
        self.residualTolerance = residualTolerance; self.overlapEigenvalueFloor = overlapEigenvalueFloor
        self.diisHistory = diisHistory; self.initialDamping = initialDamping
    }
    public func validate() throws {
        guard (2...10_000).contains(maximumIterations), (2...16).contains(diisHistory),
              [energyToleranceHartree, densityTolerance, residualTolerance, overlapEigenvalueFloor]
                .allSatisfy({ $0.isFinite && $0 > 0 && $0 < 1 }),
              initialDamping.isFinite, (0...0.8).contains(initialDamping) else {
            throw VivoQMError.invalid("SCF convergence settings")
        }
    }
}

public struct VivoSCFIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let energyHartree: Double
    public let densityChangeRMS: Double
    public let commutatorRMS: Double
}

public struct VivoHartreeFockResult: Codable, Sendable {
    public var method: VivoHartreeFockMethod
    public var basisIdentity: String
    public var alphaElectrons: Int
    public var betaElectrons: Int
    public var energyHartree: Double
    public var coefficientsAlpha: VivoQMMatrix
    public var coefficientsBeta: VivoQMMatrix
    public var densityAlpha: VivoQMMatrix
    public var densityBeta: VivoQMMatrix
    public var orbitalEnergiesAlpha: [Double]
    public var orbitalEnergiesBeta: [Double]
    public var spinSquared: Double
    public var iterations: [VivoSCFIteration]

    public func embeddedHamiltonian(from integrals: VivoGaussianIntegralSet) throws -> VivoEmbeddedHamiltonian {
        try integrals.validate()
        guard integrals.problem.basisIdentity == basisIdentity,
              alphaElectrons == integrals.problem.alphaElectrons,
              betaElectrons == integrals.problem.betaElectrons else {
            throw VivoQMError.invalid("SCF result and Gaussian problem identity differ")
        }
        let c = coefficientsAlpha
        let n = integrals.problem.basis.count
        try c.validate()
        guard c.rows == n, c.columns == n else { throw VivoQMError.invalid("SCF coefficient shape") }
        let metric = try integrals.overlap.transformed(by: c)
        guard try metric.adding(.identity(n), scale: -1).frobeniusNorm < 1e-7 else {
            throw VivoQMError.invalid("SCF orbitals are not orthonormal in the supplied AO metric")
        }
        return try .init(orbitalIDs: (0..<n).map { "mo:\($0)" },
                         basisIdentity: basisIdentity + ":canonical-alpha",
                         alphaElectrons: alphaElectrons, betaElectrons: betaElectrons,
                         constantEnergyHartree: integrals.constantEnergyHartree,
                         oneElectron: integrals.coreHamiltonian.transformed(by: c),
                         electronRepulsion: VivoEmbeddedHamiltonian.transformERI(integrals.electronRepulsion,
                                                                                 dimension: n, coefficients: c),
                         provenance: ["referenceMethod": method.rawValue, "precision": "FP64",
                                      "orbitalMetric": "orthonormal", "spinModeOrder": "interleaved-alpha-beta"])
    }
}

public enum VivoHartreeFock {
    public static func solve(_ integrals: VivoGaussianIntegralSet,
                             configuration: VivoSCFConfiguration = .init()) throws -> VivoHartreeFockResult {
        try integrals.validate()
        return try solveCore(h: integrals.coreHamiltonian, g: integrals.electronRepulsion,
                             overlap: integrals.overlap, constant: integrals.constantEnergyHartree,
                             alpha: integrals.problem.alphaElectrons, beta: integrals.problem.betaElectrons,
                             basisIdentity: integrals.problem.basisIdentity, configuration: configuration)
    }
    public static func solve(_ hamiltonian: VivoEmbeddedHamiltonian,
                             configuration: VivoSCFConfiguration = .init()) throws -> VivoHartreeFockResult {
        try hamiltonian.validate()
        return try solveCore(h: hamiltonian.oneElectron, g: hamiltonian.electronRepulsion,
                             overlap: .identity(hamiltonian.orbitalCount),
                             constant: hamiltonian.constantEnergyHartree,
                             alpha: hamiltonian.alphaElectrons, beta: hamiltonian.betaElectrons,
                             basisIdentity: hamiltonian.basisIdentity, configuration: configuration)
    }
    private static func density(_ coefficients: VivoQMMatrix, occupied: Int) -> VivoQMMatrix {
        let n = coefficients.rows
        var result = VivoQMMatrix(n, n)
        for p in 0..<n { for q in 0..<n { for i in 0..<occupied {
            result[p, q] += coefficients[p, i] * coefficients[q, i]
        } } }
        return result
    }
    private static func solveCore(h: VivoQMMatrix, g: [Double], overlap s: VivoQMMatrix,
                                  constant: Double, alpha: Int, beta: Int,
                                  basisIdentity: String, configuration config: VivoSCFConfiguration) throws -> VivoHartreeFockResult {
        try config.validate()
        let n = h.rows
        if config.method == .rhf && alpha != beta { throw VivoQMError.unsupported("RHF requires equal alpha/beta electron counts; select UHF") }
        let metric = try s.symmetricEigen()
        guard metric.values.allSatisfy({ $0 > config.overlapEigenvalueFloor }) else {
            throw VivoQMError.convergence("linearly dependent AO basis; no orbitals were silently deleted")
        }
        var x = metric.vectors
        for i in 0..<n { for j in 0..<n { x[i, j] /= sqrt(metric.values[j]) } }
        func diagonalize(_ f: VivoQMMatrix) throws -> VivoQMMatrix {
            try x.multiplied(by: f.transformed(by: x).symmetricEigen().vectors)
        }
        func fock(_ pa: VivoQMMatrix, _ pb: VivoQMMatrix) -> (VivoQMMatrix, VivoQMMatrix) {
            var fa = h, fb = h
            for p in 0..<n { for q in 0..<n { for r in 0..<n { for t in 0..<n {
                let coulomb = g[((p * n + q) * n + r) * n + t]
                let exchange = g[((p * n + r) * n + q) * n + t]
                let j = (pa[r, t] + pb[r, t]) * coulomb
                fa[p, q] += j - pa[r, t] * exchange
                fb[p, q] += j - pb[r, t] * exchange
            } } } }
            return (fa, fb)
        }
        func energy(_ pa: VivoQMMatrix, _ pb: VivoQMMatrix,
                    _ fa: VivoQMMatrix, _ fb: VivoQMMatrix) -> Double {
            var value = constant
            for p in 0..<n { for q in 0..<n {
                value += 0.5 * ((pa[p, q] + pb[p, q]) * h[p, q] + pa[p, q] * fa[p, q] + pb[p, q] * fb[p, q])
            } }
            return value
        }
        func residual(_ f: VivoQMMatrix, _ p: VivoQMMatrix) throws -> [Double] {
            let fps = try f.multiplied(by: p).multiplied(by: s)
            let spf = try s.multiplied(by: p).multiplied(by: f)
            return try fps.adding(spf, scale: -1).transformed(by: x).values
        }
        let initial = try diagonalize(h)
        var pa = density(initial, occupied: alpha), pb = density(initial, occupied: beta)
        var previousEnergy = Double.infinity
        var histories: [(VivoQMMatrix, VivoQMMatrix, [Double])] = []
        var trace: [VivoSCFIteration] = []
        for iteration in 1...config.maximumIterations {
            let (rawA, rawB) = fock(pa, pb)
            let error = try residual(rawA, pa) + residual(rawB, pb)
            histories.append((rawA, rawB, error))
            if histories.count > config.diisHistory { histories.removeFirst() }
            var fa = rawA, fb = rawB
            if histories.count >= 2 {
                let count = histories.count
                var equations = VivoQMMatrix(count + 1, count + 1)
                var rhs = Array(repeating: 0.0, count: count + 1); rhs[count] = -1
                var maximum = 0.0
                for i in 0..<count { for j in 0..<count {
                    let dot = zip(histories[i].2, histories[j].2).reduce(0) { $0 + $1.0 * $1.1 }
                    equations[i, j] = dot; maximum = max(maximum, abs(dot))
                } }
                // Scaling the error block leaves the DIIS coefficients invariant.
                if maximum > 1e-28 {
                    for i in 0..<count { for j in 0..<count { equations[i, j] /= maximum } }
                    for i in 0..<count { equations[i, count] = -1; equations[count, i] = -1 }
                    if let weights = try? equations.solve(rhs), weights.prefix(count).allSatisfy({ abs($0) < 1e5 }) {
                        fa = VivoQMMatrix(n, n); fb = VivoQMMatrix(n, n)
                        for i in 0..<count {
                            fa = try fa.adding(histories[i].0, scale: weights[i])
                            fb = try fb.adding(histories[i].1, scale: weights[i])
                        }
                    }
                }
            }
            let ca = try diagonalize(fa)
            let cb = config.method == .rhf ? ca : try diagonalize(fb)
            var nextA = density(ca, occupied: alpha), nextB = density(cb, occupied: beta)
            if iteration <= 3 && config.initialDamping > 0 {
                nextA = try nextA.scaled(1 - config.initialDamping).adding(pa, scale: config.initialDamping)
                nextB = try nextB.scaled(1 - config.initialDamping).adding(pb, scale: config.initialDamping)
            }
            let changeA = try nextA.adding(pa, scale: -1).frobeniusNorm
            let changeB = try nextB.adding(pb, scale: -1).frobeniusNorm
            let change = sqrt((changeA * changeA + changeB * changeB) / Double(2 * n * n))
            let (physicalA, physicalB) = fock(nextA, nextB)
            let nextEnergy = energy(nextA, nextB, physicalA, physicalB)
            let errors = try residual(physicalA, nextA) + residual(physicalB, nextB)
            let rms = sqrt(errors.reduce(0) { $0 + $1 * $1 } / Double(errors.count))
            guard nextEnergy.isFinite, rms.isFinite, change.isFinite else { throw VivoQMError.convergence("SCF produced non-finite values") }
            trace.append(.init(iteration: iteration, energyHartree: nextEnergy, densityChangeRMS: change, commutatorRMS: rms))
            if iteration > 3 && abs(nextEnergy - previousEnergy) < config.energyToleranceHartree,
               change < config.densityTolerance, rms < config.residualTolerance {
                // Canonicalize with the *physical* Fock, never the DIIS extrapolate.
                let canonicalA = try diagonalize(physicalA)
                let canonicalB = config.method == .rhf ? canonicalA : try diagonalize(physicalB)
                let finalA = density(canonicalA, occupied: alpha), finalB = density(canonicalB, occupied: beta)
                let (finalFA, finalFB) = fock(finalA, finalB)
                let finalEnergy = energy(finalA, finalB, finalFA, finalFB)
                let finalErrors = try residual(finalFA, finalA) + residual(finalFB, finalB)
                let finalRMS = sqrt(finalErrors.reduce(0) { $0 + $1 * $1 } / Double(finalErrors.count))
                if abs(finalEnergy - nextEnergy) < config.energyToleranceHartree && finalRMS < config.residualTolerance {
                    let moa = try finalFA.transformed(by: canonicalA), mob = try finalFB.transformed(by: canonicalB)
                    let spinOverlap = try canonicalA.transposed.multiplied(by: s).multiplied(by: canonicalB)
                    let sz = Double(alpha - beta) / 2
                    var s2 = sz * sz + Double(alpha + beta) / 2
                    for i in 0..<alpha { for j in 0..<beta { s2 -= spinOverlap[i, j] * spinOverlap[i, j] } }
                    return .init(method: config.method, basisIdentity: basisIdentity,
                                 alphaElectrons: alpha, betaElectrons: beta, energyHartree: finalEnergy,
                                 coefficientsAlpha: canonicalA, coefficientsBeta: canonicalB,
                                 densityAlpha: finalA, densityBeta: finalB,
                                 orbitalEnergiesAlpha: (0..<n).map { moa[$0, $0] },
                                 orbitalEnergiesBeta: (0..<n).map { mob[$0, $0] },
                                 spinSquared: max(0, s2), iterations: trace)
                }
            }
            pa = nextA; pb = nextB; previousEnergy = nextEnergy
        }
        throw VivoQMError.convergence("SCF exhausted \(config.maximumIterations) iterations; last record: \(String(describing: trace.last))")
    }
}

public struct VivoMP2Result: Codable, Sendable {
    public let referenceEnergyHartree: Double
    public let correlationEnergyHartree: Double
    public let totalEnergyHartree: Double
    public let minimumDenominatorMagnitudeHartree: Double?
}

public enum VivoRestrictedMP2 {
    public static func solve(integrals: VivoGaussianIntegralSet, reference: VivoHartreeFockResult,
                             denominatorFloorHartree: Double = 1e-8) throws -> VivoMP2Result {
        guard reference.method == .rhf, reference.alphaElectrons == reference.betaElectrons,
              denominatorFloorHartree.isFinite, denominatorFloorHartree > 0 else {
            throw VivoQMError.unsupported("canonical restricted MP2 requires a closed-shell RHF reference")
        }
        let h = try reference.embeddedHamiltonian(from: integrals)
        let n = h.orbitalCount, occupied = h.alphaElectrons
        guard reference.orbitalEnergiesAlpha.count == n,
              reference.orbitalEnergiesAlpha.allSatisfy(\.isFinite) else { throw VivoQMError.invalid("MP2 orbital energies") }
        var hfEnergy = h.constantEnergyHartree
        for i in 0..<occupied {
            hfEnergy += 2 * h.oneElectron[i, i]
            for j in 0..<occupied { hfEnergy += 2 * h.eri(i, i, j, j) - h.eri(i, j, j, i) }
        }
        guard reference.energyHartree.isFinite, abs(hfEnergy - reference.energyHartree) < 1e-7 else {
            throw VivoQMError.invalid("MP2 reference energy does not match the supplied integrals")
        }
        var correlation = 0.0, minimum: Double?
        let epsilon = reference.orbitalEnergiesAlpha
        for i in 0..<occupied { for j in 0..<occupied { for a in occupied..<n { for b in occupied..<n {
            let denominator = epsilon[i] + epsilon[j] - epsilon[a] - epsilon[b]
            guard denominator < -denominatorFloorHartree else {
                throw VivoQMError.convergence("MP2 near-zero or inverted energy denominator; select a multireference solver")
            }
            minimum = min(minimum ?? .infinity, abs(denominator))
            let direct = h.eri(i, a, j, b), exchange = h.eri(i, b, j, a)
            correlation += direct * (2 * direct - exchange) / denominator
        } } } }
        guard correlation.isFinite else { throw VivoQMError.convergence("non-finite MP2 correlation energy") }
        return .init(referenceEnergyHartree: hfEnergy, correlationEnergyHartree: correlation,
                     totalEnergyHartree: hfEnergy + correlation,
                     minimumDenominatorMagnitudeHartree: minimum)
    }
}
