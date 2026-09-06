import Foundation

/// Exactly one approximation family is used in a ladder. CAS and ECC results
/// are not presented as successive refinements of the same approximation.
public enum VivoBarrierLevelMethod: Codable, Sendable, Equatable {
    case casci(partition: VivoActiveSpace)
    case eccPath(configuration: VivoECCPathConfiguration)
}
public struct VivoBarrierLevel: Codable, Sendable, Equatable {
    public let identifier: String
    public let method: VivoBarrierLevelMethod
    public init(identifier: String, method: VivoBarrierLevelMethod) {
        self.identifier = identifier; self.method = method
    }
}
public struct VivoBarrierAcceptance: Codable, Sendable, Equatable {
    public var maximumBarrierErrorHartree: Double
    public var maximumRelativeProfileErrorHartree: Double
    public var maximumSuccessiveChangeHartree: Double
    /// At least two genuinely reduced levels must agree with the reference and
    /// one another. A full-space endpoint is never counted in this window.
    public var minimumStableReducedLevels: Int
    public init(maximumBarrierErrorHartree: Double = 1e-3,
                maximumRelativeProfileErrorHartree: Double = 1e-3,
                maximumSuccessiveChangeHartree: Double = 1e-3,
                minimumStableReducedLevels: Int = 2) {
        self.maximumBarrierErrorHartree = maximumBarrierErrorHartree
        self.maximumRelativeProfileErrorHartree = maximumRelativeProfileErrorHartree
        self.maximumSuccessiveChangeHartree = maximumSuccessiveChangeHartree
        self.minimumStableReducedLevels = minimumStableReducedLevels
    }
    public func validate() throws {
        guard [maximumBarrierErrorHartree, maximumRelativeProfileErrorHartree,
               maximumSuccessiveChangeHartree].allSatisfy({ $0.isFinite && $0 > 0 }),
              (2...16).contains(minimumStableReducedLevels) else {
            throw VivoChemistryError.invalid("barrier accuracy tolerances or reduced-level stability window")
        }
    }
}
public struct VivoBarrierConvergenceRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/barrier-convergence/v1"
    public let schema: String
    public let identifier: String
    public let atomIdentifiers: [String]
    public let coordinateUnit: String
    public let snapshots: [VivoMolecularPathSnapshot]
    public let basis: VivoGaussianBasis
    public let anchorPointIdentifier: String
    /// This is a declared evaluation point, NOT an assertion of a saddle.
    public let barrierPointIdentifier: String
    public let transportGroups: [[Int]]
    public let anchorCoefficients: VivoQMMatrix?
    public let minimumTransportSingularValue: Double
    public let minimumReferenceOverlapSquared: Double
    public let referenceResidualTolerance: Double
    public let levels: [VivoBarrierLevel]
    public let acceptance: VivoBarrierAcceptance
    /// Bounds point solves plus reserved outer ECC point evaluations. Nested
    /// ECC iterations keep their own explicit budgets; this is not a FLOP count.
    public let maximumPointEvaluations: Int
    public let budget: VivoChemistryBudget
    public init(identifier: String, atomIdentifiers: [String], coordinateUnit: String,
                snapshots: [VivoMolecularPathSnapshot], basis: VivoGaussianBasis,
                anchorPointIdentifier: String, barrierPointIdentifier: String,
                transportGroups: [[Int]], anchorCoefficients: VivoQMMatrix? = nil,
                minimumTransportSingularValue: Double = 0.5,
                minimumReferenceOverlapSquared: Double = 0.1,
                referenceResidualTolerance: Double = 1e-11,
                levels: [VivoBarrierLevel], acceptance: VivoBarrierAcceptance = .init(),
                maximumPointEvaluations: Int = 10000, budget: VivoChemistryBudget = .init()) {
        schema = Self.schema; self.identifier = identifier; self.atomIdentifiers = atomIdentifiers
        self.coordinateUnit = coordinateUnit; self.snapshots = snapshots; self.basis = basis
        self.anchorPointIdentifier = anchorPointIdentifier; self.barrierPointIdentifier = barrierPointIdentifier
        self.transportGroups = transportGroups; self.anchorCoefficients = anchorCoefficients
        self.minimumTransportSingularValue = minimumTransportSingularValue
        self.minimumReferenceOverlapSquared = minimumReferenceOverlapSquared
        self.referenceResidualTolerance = referenceResidualTolerance; self.levels = levels
        self.acceptance = acceptance; self.maximumPointEvaluations = maximumPointEvaluations; self.budget = budget
    }
    public func validate() throws {
        try budget.validate(); try acceptance.validate()
        guard schema == Self.schema, !identifier.isEmpty, identifier.utf8.count <= 1024,
              !coordinateUnit.isEmpty, coordinateUnit.utf8.count <= 128,
              (3...64).contains(snapshots.count), (2...32).contains(levels.count),
              Set(snapshots.map(\.identifier)).count == snapshots.count,
              Set(levels.map(\.identifier)).count == levels.count,
              levels.allSatisfy({ !$0.identifier.isEmpty && $0.identifier.utf8.count <= 1024 }),
              snapshots.contains(where: { $0.identifier == anchorPointIdentifier }),
              snapshots.dropFirst().dropLast().contains(where: { $0.identifier == barrierPointIdentifier }),
              (1...1000000).contains(maximumPointEvaluations),
              [minimumTransportSingularValue, minimumReferenceOverlapSquared].allSatisfy({ $0.isFinite && $0 > 0 && $0 <= 1 }),
              referenceResidualTolerance.isFinite, referenceResidualTolerance > 0,
              referenceResidualTolerance <= min(1e-8, acceptance.maximumBarrierErrorHartree / 100) else {
            throw VivoChemistryError.invalid("barrier campaign identities, points, reference or work contract")
        }
        let first = snapshots[0].system
        guard atomIdentifiers.count == first.nuclei.count, Set(atomIdentifiers).count == atomIdentifiers.count,
              atomIdentifiers.allSatisfy({ !$0.isEmpty && $0.utf8.count <= 1024 }) else {
            throw VivoChemistryError.invalid("barrier campaign requires explicit unique mapped atoms")
        }
        try basis.validate(nucleusCount: first.nuclei.count)
        var coordinate = -Double.infinity
        for point in snapshots {
            try point.system.validate()
            guard !point.identifier.isEmpty, point.identifier.utf8.count <= 1024,
                  point.coordinate.isFinite, point.coordinate > coordinate,
                  point.system.nuclei.map(\.atomicNumber) == first.nuclei.map(\.atomicNumber),
                  point.system.nuclei.map(\.structureAtomIndex) == first.nuclei.map(\.structureAtomIndex),
                  point.system.alphaElectrons == first.alphaElectrons, point.system.betaElectrons == first.betaElectrons,
                  point.system.pointCharges == first.pointCharges else {
                throw VivoChemistryError.invalid("campaign changes atom mapping, charge/spin, frozen environment or coordinate order")
            }
            coordinate = point.coordinate
        }
        let n = try VivoGaussianIntegralEngine.expanded(system: first, basis: basis, budget: budget).count
        let groups = transportGroups.flatMap { $0 }
        guard n <= 31, n <= budget.maximumBasisFunctions, groups.count == n, Set(groups) == Set(0..<n),
              transportGroups.allSatisfy({ !$0.isEmpty }), first.alphaElectrons <= n, first.betaElectrons <= n else {
            throw VivoChemistryError.invalid("campaign basis capacity, orbital transport groups or populations")
        }
        // Check complete reference-sector capacity before expensive integral work.
        _ = try VivoDirectCI.determinants(n: n, na: first.alphaElectrons, nb: first.betaElectrons, budget: budget)
        _ = try budget.elements([snapshots.count, n, n, n, n], simultaneousArrays: 8)
        if let c = anchorCoefficients {
            guard c.rows == n, c.columns == n, c.values.allSatisfy(\.isFinite) else {
                throw VivoChemistryError.invalid("campaign anchor coefficient shape")
            }
        }
        var previous: VivoBarrierLevelMethod?, reservation = snapshots.count
        for level in levels {
            switch level.method {
            case .casci(let p):
                let all = p.doublyOccupiedCore + p.active
                guard !p.active.isEmpty, p.frozenOrbitals.isEmpty, Set(all).count == all.count,
                      all.allSatisfy({ $0 >= 0 && $0 < n }), p.doublyOccupiedCore.count <= min(first.alphaElectrons, first.betaElectrons),
                      first.alphaElectrons - p.doublyOccupiedCore.count <= p.active.count,
                      first.betaElectrons - p.doublyOccupiedCore.count <= p.active.count else {
                    throw VivoChemistryError.invalid("campaign CAS core, active, external or electron partition")
                }
                if let previous {
                    guard case .casci(let old) = previous,
                          Set(old.active).isStrictSubset(of: Set(p.active)),
                          Set(p.doublyOccupiedCore).isSubset(of: Set(old.doublyOccupiedCore)),
                          Set(old.doublyOccupiedCore).subtracting(p.doublyOccupiedCore).isSubset(of: Set(p.active)) else {
                        throw VivoChemistryError.invalid("CAS ladder must enlarge a nested active space without changing families or demoting occupied core to empty external orbitals")
                    }
                }
                reservation += snapshots.count
            case .eccPath(let cfg):
                guard first.alphaElectrons == first.betaElectrons,
                      cfg.expectedBathOrbitals.count == cfg.embedding.fragments.count,
                      cfg.maximumPointEvaluations >= snapshots.count, !cfg.embedding.fragments.isEmpty else { throw VivoChemistryError.invalid("ECC campaign spin or bath contract") }
                // Structural preflight only: no model energies are inferred
                // from this zero-valued shape carrier.
                let shape = VivoEmbeddedHamiltonian(orbitalIdentifiers: (0..<n).map { "shape-\($0)" },
                    alphaElectrons: first.alphaElectrons, betaElectrons: first.betaElectrons,
                    oneElectron: VivoQMMatrix(n,n), twoElectron: [Double](repeating: 0, count: n*n*n*n),
                    constantEnergyHartree: 0, energyReference: "shape-validation-only")
                let shapePoints = try snapshots.indices.map { i in try VivoECCPathPoint(identifier: snapshots[i].identifier,
                    hamiltonian: shape, overlapWithPrevious: i == 0 ? nil : VivoQMMatrix.identity(n)) }
                try cfg.validate(points: shapePoints, budget: budget)
                if let previous {
                    guard case .eccPath(let old) = previous, old.embedding.mode == cfg.embedding.mode,
                          old.embedding.matching == cfg.embedding.matching,
                          old.embedding.qioWeights == cfg.embedding.qioWeights,
                          old.embedding.localityGroups == cfg.embedding.localityGroups,
                          old.pointWeights == cfg.pointWeights, old.transportGroups == cfg.transportGroups,
                          old.embedding.fragments.map(\.identifier) == cfg.embedding.fragments.map(\.identifier) else {
                        throw VivoChemistryError.invalid("ECC ladder changes solver family, reference, matching, objective or fragment identities")
                    }
                    var expanded = false
                    for i in cfg.embedding.fragments.indices {
                        let a = old.embedding.fragments[i], b = cfg.embedding.fragments[i]
                        let oldSize = a.orbitals.count + old.expectedBathOrbitals[i]
                        let newSize = b.orbitals.count + cfg.expectedBathOrbitals[i]
                        guard Set(a.orbitals).isSubset(of: Set(b.orbitals)), newSize >= oldSize else {
                            throw VivoChemistryError.invalid("ECC ladder shrinks a fragment or its declared cluster")
                        }
                        expanded = expanded || newSize > oldSize
                    }
                    guard expanded else { throw VivoChemistryError.invalid("duplicate-sized ECC level is not an enlargement") }
                }
                guard cfg.maximumPointEvaluations <= maximumPointEvaluations - reservation else {
                    throw VivoChemistryError.resourceLimit("campaign ECC point-evaluation reservation")
                }
                reservation += cfg.maximumPointEvaluations
            }
            guard reservation <= maximumPointEvaluations else { throw VivoChemistryError.resourceLimit("campaign point-evaluation reservation") }
            previous = level.method
        }
    }
}

public struct VivoBarrierDifferences: Codable, Sendable, Equatable {
    public let forwardHartree: Double
    public let reverseHartree: Double
    public let reactionHartree: Double
    public let relativeProfileHartree: [Double]
}
public struct VivoBarrierLevelEvaluation: Codable, Sendable, Equatable {
    public let energiesHartree: [Double]
    public let differences: VivoBarrierDifferences
    public let maximumBarrierErrorHartree: Double
    public let maximumRelativeProfileErrorHartree: Double
    public let maximumAbsoluteEnergyErrorHartree: Double
    /// Nil only for the first level or after a failed immediate predecessor.
    public let maximumSuccessiveChangeHartree: Double?
    public let maximumImpurityOrbitals: Int
    public let isGenuinelyReduced: Bool
    public let maximumSolverResidual: Double
    public let chargedPointEvaluations: Int
    /// Explicitly diagnoses incorrect error cancellation rather than concealing it.
    public let meetsReferenceAccuracy: Bool
}
public enum VivoBarrierLevelOutcome: Codable, Sendable, Equatable {
    case evaluated(result: VivoBarrierLevelEvaluation)
    case failed(category: String, reason: String)
    public var evaluation: VivoBarrierLevelEvaluation? {
        if case .evaluated(let result) = self { return result }; return nil
    }
}
public struct VivoBarrierLevelReport: Codable, Sendable, Equatable {
    public let identifier: String
    public let outcome: VivoBarrierLevelOutcome
}
public enum VivoBarrierAssessment: String, Codable, Sendable {
    case reducedAccuracyEstablished, reducedAccuracyNotEstablished, incompleteLevelExecution
}
public struct VivoBarrierConvergenceResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/barrier-convergence-result/v1"
    public let schema: String
    public let request: VivoBarrierConvergenceRequest
    /// The actual complete AO frames used for every level and reference.
    public let orbitalCoefficients: [VivoQMMatrix]
    public let transportMinimumSingularValues: [Double]
    public let adjacentReferenceOverlapsSquared: [Double]
    public let referenceEnergiesHartree: [Double]
    public let referenceResiduals: [Double]
    public let referenceDifferences: VivoBarrierDifferences
    public let levels: [VivoBarrierLevelReport]
    public let assessment: VivoBarrierAssessment
    public let acceptedReducedLevelIdentifier: String?
    public let chargedPointEvaluations: Int
    public let meaning: String
}

/// Accuracy assessment at identical, explicitly declared geometries. No nuclear
/// optimization, stationary-point, thermochemical or rate claim is made here.
/// A failed approximation is retained, not filtered out of a favorable report.
public enum VivoBarrierConvergence {
    public static let meaning = "fixed-geometry electronic barrier/profile convergence against full-space FCI; first and last supplied structures define references; not a saddle, connectivity, Gibbs energy, rate or paper-reproduction certificate"
    private struct Prepared {
        let hamiltonians: [VivoEmbeddedHamiltonian]
        let coefficients: [VivoQMMatrix]
        let minima: [Double]
        let overlaps: [VivoQMMatrix]
    }
    private static func prepare(_ r: VivoBarrierConvergenceRequest) throws -> Prepared {
        try r.validate()
        let integrals = try r.snapshots.map { try VivoGaussianIntegralEngine.compute(system: $0.system, basis: r.basis, budget: r.budget) }
        let anchor = r.snapshots.firstIndex { $0.identifier == r.anchorPointIdentifier }!, n = integrals[0].count
        var candidates: [VivoQMMatrix] = []
        for ao in integrals {
            let spectrum = try VivoQMDenseAlgebra.symmetricEigen(ao.overlap)
            guard spectrum.values[0] >= 1e-8 else { throw VivoChemistryError.invalid("campaign AO space lost rank") }
            var x = spectrum.vectors
            for i in 0..<n { for j in 0..<n { x[i,j] /= sqrt(spectrum.values[j]) } }
            let orbitals = try VivoQMDenseAlgebra.symmetricEigen(ao.coreHamiltonian.congruence(x))
            candidates.append(try x.multiplied(by: orbitals.vectors))
        }
        var coefficients = candidates, minima = [Double](repeating: 1, count: candidates.count)
        coefficients[anchor] = r.anchorCoefficients ?? candidates[anchor]
        guard try integrals[anchor].overlap.congruence(coefficients[anchor]).adding(.identity(n), scale: -1).frobeniusNorm < 1e-8 else {
            throw VivoChemistryError.invalid("campaign anchor is not AO-metric orthonormal")
        }
        let cross = try (0..<(r.snapshots.count-1)).map { i in
            try VivoGaussianIntegralEngine.crossOverlap(leftSystem: r.snapshots[i].system, leftBasis: r.basis,
                rightSystem: r.snapshots[i+1].system, rightBasis: r.basis, budget: r.budget)
        }
        func align(_ from: Int, _ to: Int, _ metric: VivoQMMatrix) throws {
            let transport = try VivoOrbitalPathTransport.align(referenceCoefficients: coefficients[from],
                currentCoefficients: candidates[to], crossAOOverlap: metric, groups: r.transportGroups,
                minimumSingularValue: r.minimumTransportSingularValue)
            coefficients[to] = try candidates[to].multiplied(by: transport.rotation)
            minima[to] = transport.minimumSubspaceSingularValue
        }
        if anchor > 0 { for i in stride(from: anchor-1, through: 0, by: -1) { try align(i+1, i, cross[i].transposed) } }
        if anchor+1 < r.snapshots.count { for i in (anchor+1)..<r.snapshots.count { try align(i-1, i, cross[i-1]) } }
        let h = try integrals.indices.map { i in
            try VivoEmbeddedHamiltonian.fromAO(integrals[i], coefficients: coefficients[i],
                alphaElectrons: r.snapshots[i].system.alphaElectrons, betaElectrons: r.snapshots[i].system.betaElectrons,
                orbitalIdentifiers: (0..<n).map { "campaign-orbital-\($0)" },
                energyReference: "physical electronic Hamiltonian; geometry-specific AO scalar once; no thermal or standard-state correction", budget: r.budget)
        }
        let overlaps = try cross.indices.map { i in try coefficients[i].transposed.multiplied(by: cross[i]).multiplied(by: coefficients[i+1]) }
        return .init(hamiltonians: h, coefficients: coefficients, minima: minima, overlaps: overlaps)
    }
    private static func differences(_ energies: [Double], barrier: Int) -> VivoBarrierDifferences {
        .init(forwardHartree: energies[barrier]-energies[0], reverseHartree: energies[barrier]-energies.last!,
              reactionHartree: energies.last!-energies[0], relativeProfileHartree: energies.map { $0-energies[0] })
    }
    private static func barrierChange(_ a: VivoBarrierDifferences, _ b: VivoBarrierDifferences) -> Double {
        max(abs(a.forwardHartree-b.forwardHartree), abs(a.reverseHartree-b.reverseHartree), abs(a.reactionHartree-b.reactionHartree))
    }
    private static func profileChange(_ a: VivoBarrierDifferences, _ b: VivoBarrierDifferences) -> Double {
        zip(a.relativeProfileHartree, b.relativeProfileHartree).map { abs($0-$1) }.max()!
    }
    public static func run(_ request: VivoBarrierConvergenceRequest) throws -> VivoBarrierConvergenceResult {
        let p = try prepare(request), n = p.hamiltonians[0].orbitalCount
        let barrier = request.snapshots.firstIndex { $0.identifier == request.barrierPointIdentifier }!
        let reference = try p.hamiltonians.map { try VivoDirectCI.solve($0,
            configuration: .init(residualTolerance: request.referenceResidualTolerance), budget: request.budget).roots[0] }
        // Physical CI residuals are measured again with the shared Hamiltonian
        // action, including all fixed-sector components, not a serialized flag.
        let residuals = try reference.indices.map { i in try VivoDirectCI.residualNorm(hamiltonian: p.hamiltonians[i],
            state: reference[i].state, energyHartree: reference[i].energyHartree, budget: request.budget) }
        guard residuals.allSatisfy({ $0 <= 1.01*request.referenceResidualTolerance }) else {
            throw VivoChemistryError.convergence("campaign full-space reference physical residual")
        }
        let overlaps = try p.overlaps.indices.map { i -> Double in
            let value = try VivoStateFollowing.overlaps(previous: [reference[i].state], candidates: [reference[i+1].state],
                orbitalOverlap: p.overlaps[i], budget: request.budget)[0,0]
            guard value.isFinite, value*value <= 1+1e-6, value*value >= request.minimumReferenceOverlapSquared else {
                throw VivoChemistryError.convergence("campaign ground reference loses physical continuity; refine the path or specify a multistate method")
            }
            return value*value
        }
        let energies = reference.map(\.energyHartree), truth = differences(energies, barrier: barrier)
        var reports: [VivoBarrierLevelReport] = [], calls = reference.count
        let eccPoints = p.hamiltonians.indices.map { i in VivoECCPathPoint(identifier: request.snapshots[i].identifier,
            hamiltonian: p.hamiltonians[i], overlapWithPrevious: i == 0 ? nil : p.overlaps[i-1]) }
        for level in request.levels {
            do {
                let values: [Double], size: Int, solverResidual: Double, chargeForLevel: Int
                switch level.method {
                case .casci(let partition):
                    var states: [VivoCIResult] = []
                    for h in p.hamiltonians {
                        calls += 1
                        guard calls <= request.maximumPointEvaluations else { throw VivoChemistryError.resourceLimit("campaign point solve limit") }
                        let active = try h.frozenCore(active: partition.active, doublyOccupiedCore: partition.doublyOccupiedCore, budget: request.budget)
                        let ci = try VivoDirectCI.solve(active, configuration: .init(residualTolerance: request.referenceResidualTolerance), budget: request.budget).roots[0]
                        let residual = try VivoDirectCI.residualNorm(hamiltonian: active, state: ci.state, energyHartree: ci.energyHartree, budget: request.budget)
                        guard residual <= 1.01*request.referenceResidualTolerance else { throw VivoChemistryError.convergence("campaign CAS physical residual") }
                        states.append(ci)
                    }
                    chargeForLevel = states.count
                    values = states.map(\.energyHartree); size = partition.active.count
                    solverResidual = states.map(\.eigenResidualNorm).max()!
                case .eccPath(let cfg):
                    // Reserve all work before entry, including a failing nested
                    // solve whose partial count is intentionally not inferred.
                    chargeForLevel = cfg.maximumPointEvaluations
                    calls += cfg.maximumPointEvaluations
                    guard calls <= request.maximumPointEvaluations else { throw VivoChemistryError.resourceLimit("campaign ECC reservation") }
                    let result = try VivoECCReactionPath.solve(points: eccPoints, configuration: cfg, budget: request.budget)
                    guard result.converged else { throw VivoChemistryError.convergence("campaign ECC path: \(result.termination)") }
                    try VivoECCReactionPath.validate(result, points: eccPoints, configuration: cfg, budget: request.budget)
                    values = result.pointResults.map(\.energyHartree)
                    size = result.pointResults.flatMap { $0.frame.clusters.map { $0.coefficients.columns } }.max()!
                    solverResidual = result.pointResults.flatMap { $0.frame.states.map { $0.biasedCI.eigenResidualNorm } }.max()!
                }
                let diff = differences(values, barrier: barrier), bError = barrierChange(diff, truth), pError = profileChange(diff, truth)
                let previous = reports.last?.outcome.evaluation?.differences
                let successive = previous.map { max(barrierChange(diff, $0), profileChange(diff, $0)) }
                let evaluation = VivoBarrierLevelEvaluation(energiesHartree: values, differences: diff,
                    maximumBarrierErrorHartree: bError, maximumRelativeProfileErrorHartree: pError,
                    maximumAbsoluteEnergyErrorHartree: zip(values, energies).map { abs($0-$1) }.max()!,
                    maximumSuccessiveChangeHartree: successive, maximumImpurityOrbitals: size,
                    isGenuinelyReduced: size < n, maximumSolverResidual: solverResidual,
                    chargedPointEvaluations: chargeForLevel,
                    meetsReferenceAccuracy: bError <= request.acceptance.maximumBarrierErrorHartree && pError <= request.acceptance.maximumRelativeProfileErrorHartree)
                reports.append(.init(identifier: level.identifier, outcome: .evaluated(result: evaluation)))
            } catch let error as VivoChemistryError {
                switch error {
                case .convergence(let message): reports.append(.init(identifier: level.identifier, outcome: .failed(category: "numerical-convergence", reason: message)))
                case .resourceLimit(let message): reports.append(.init(identifier: level.identifier, outcome: .failed(category: "resource-limit", reason: message)))
                default: throw error // malformed methods are not scientific observations
                }
            } catch let error as VivoECCPathError {
                reports.append(.init(identifier: level.identifier, outcome: .failed(category: "subspace-continuity", reason: String(describing: error))))
            }
        }
        let decision = assessment(reports, acceptance: request.acceptance)
        return .init(schema: VivoBarrierConvergenceResult.schema, request: request, orbitalCoefficients: p.coefficients,
            transportMinimumSingularValues: p.minima, adjacentReferenceOverlapsSquared: overlaps,
            referenceEnergiesHartree: energies, referenceResiduals: residuals, referenceDifferences: truth, levels: reports,
            assessment: decision.0, acceptedReducedLevelIdentifier: decision.1, chargedPointEvaluations: calls, meaning: meaning)
    }
    private static func assessment(_ reports: [VivoBarrierLevelReport], acceptance a: VivoBarrierAcceptance) -> (VivoBarrierAssessment, String?) {
        guard reports.allSatisfy({ $0.outcome.evaluation != nil }) else { return (.incompleteLevelExecution, nil) }
        let values = reports.map { $0.outcome.evaluation! }
        let reduced = values.indices.filter { values[$0].isGenuinelyReduced }
        guard let last = reduced.last, reduced.count >= a.minimumStableReducedLevels else { return (.reducedAccuracyNotEstablished, nil) }
        let window = Array(reduced.suffix(a.minimumStableReducedLevels))
        guard window == Array((last-a.minimumStableReducedLevels+1)...last),
              window.allSatisfy({ values[$0].meetsReferenceAccuracy }),
              window.dropFirst().allSatisfy({ (values[$0].maximumSuccessiveChangeHartree ?? .infinity) <= a.maximumSuccessiveChangeHartree }),
              values[(last+1)...].allSatisfy({ $0.meetsReferenceAccuracy && ($0.maximumSuccessiveChangeHartree ?? .infinity) <= a.maximumSuccessiveChangeHartree }) else {
            return (.reducedAccuracyNotEstablished, nil)
        }
        return (.reducedAccuracyEstablished, reports[last].identifier)
    }
    /// Re-executes the bounded numerical campaign. Cache hits cannot promote a
    /// claimed acceptance, hide a failed level or replace the reference energies.
    public static func validate(_ result: VivoBarrierConvergenceResult, request: VivoBarrierConvergenceRequest) throws {
        guard result.schema == VivoBarrierConvergenceResult.schema, result.request == request, result.meaning == meaning else {
            throw VivoChemistryError.invalid("campaign result identity, configuration or scientific meaning")
        }
        let rebuilt = try run(request)
        func close(_ a: [Double], _ b: [Double], tolerance: Double = 1e-9) -> Bool {
            a.count == b.count && zip(a,b).allSatisfy { $0.isFinite && $1.isFinite && abs($0-$1) <= tolerance }
        }
        func differencesClose(_ a: VivoBarrierDifferences, _ b: VivoBarrierDifferences) -> Bool {
            close([a.forwardHartree,a.reverseHartree,a.reactionHartree], [b.forwardHartree,b.reverseHartree,b.reactionHartree]) && close(a.relativeProfileHartree,b.relativeProfileHartree)
        }
        guard result.assessment == rebuilt.assessment, result.acceptedReducedLevelIdentifier == rebuilt.acceptedReducedLevelIdentifier,
              result.chargedPointEvaluations == rebuilt.chargedPointEvaluations, result.levels.count == rebuilt.levels.count,
              result.orbitalCoefficients.count == rebuilt.orbitalCoefficients.count,
              close(result.transportMinimumSingularValues, rebuilt.transportMinimumSingularValues),
              close(result.adjacentReferenceOverlapsSquared, rebuilt.adjacentReferenceOverlapsSquared),
              close(result.referenceEnergiesHartree, rebuilt.referenceEnergiesHartree),
              close(result.referenceResiduals, rebuilt.referenceResiduals, tolerance: request.referenceResidualTolerance),
              differencesClose(result.referenceDifferences, rebuilt.referenceDifferences) else {
            throw VivoChemistryError.invalid("campaign reference, work or assessment differs on reconstruction")
        }
        for (a,b) in zip(result.orbitalCoefficients, rebuilt.orbitalCoefficients) {
            guard a.rows == b.rows, a.columns == b.columns, close(a.values,b.values) else { throw VivoChemistryError.invalid("campaign orbital-frame reconstruction") }
        }
        for (a,b) in zip(result.levels, rebuilt.levels) {
            guard a.identifier == b.identifier else { throw VivoChemistryError.invalid("campaign level identity") }
            switch (a.outcome,b.outcome) {
            case (.failed(let c,let why),.failed(let d,let other)):
                guard c == d, why == other else { throw VivoChemistryError.invalid("campaign failure observation differs") }
            case (.evaluated(let x),.evaluated(let y)):
                guard close(x.energiesHartree,y.energiesHartree), differencesClose(x.differences,y.differences),
                      close([x.maximumBarrierErrorHartree,x.maximumRelativeProfileErrorHartree,x.maximumAbsoluteEnergyErrorHartree,x.maximumSolverResidual],
                            [y.maximumBarrierErrorHartree,y.maximumRelativeProfileErrorHartree,y.maximumAbsoluteEnergyErrorHartree,y.maximumSolverResidual]),
                      x.maximumImpurityOrbitals == y.maximumImpurityOrbitals, x.isGenuinelyReduced == y.isGenuinelyReduced,
                      x.chargedPointEvaluations == y.chargedPointEvaluations, x.meetsReferenceAccuracy == y.meetsReferenceAccuracy,
                      (x.maximumSuccessiveChangeHartree == nil) == (y.maximumSuccessiveChangeHartree == nil),
                      close([x.maximumSuccessiveChangeHartree ?? 0], [y.maximumSuccessiveChangeHartree ?? 0]) else {
                    throw VivoChemistryError.invalid("campaign approximation energy, residual, size or accuracy claim")
                }
            default: throw VivoChemistryError.invalid("a failed campaign level was removed or presented as computed")
            }
        }
    }
}
