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
/// Reference-assisted natural orbitals of the weighted path density in the
/// physically transported common frame. Weights are inputs, never fitted to
/// barrier energies. Full reference CI is still required; no speedup is implied.
public struct VivoBarrierEnsembleOrbitals: Codable, Sendable, Equatable {
    public let pointWeights: [Double]
    public let minimumOccupationBoundaryGap: Double
    public init(pointWeights: [Double], minimumOccupationBoundaryGap: Double = 1e-7) {
        self.pointWeights = pointWeights; self.minimumOccupationBoundaryGap = minimumOccupationBoundaryGap
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
    public let ensembleOrbitals: VivoBarrierEnsembleOrbitals?
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
                maximumPointEvaluations: Int = 10000, budget: VivoChemistryBudget = .init(),
                ensembleOrbitals: VivoBarrierEnsembleOrbitals? = nil) {
        schema = Self.schema; self.identifier = identifier; self.atomIdentifiers = atomIdentifiers
        self.coordinateUnit = coordinateUnit; self.snapshots = snapshots; self.basis = basis
        self.anchorPointIdentifier = anchorPointIdentifier; self.barrierPointIdentifier = barrierPointIdentifier
        self.transportGroups = transportGroups; self.anchorCoefficients = anchorCoefficients
        self.ensembleOrbitals = ensembleOrbitals
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
              referenceResidualTolerance <= min(1e-8, min(acceptance.maximumBarrierErrorHartree, acceptance.maximumRelativeProfileErrorHartree) / 100) else {
            throw VivoChemistryError.invalid("barrier campaign identities, points, reference or work contract")
        }
        if let ensemble = ensembleOrbitals {
            guard ensemble.pointWeights.count == snapshots.count,
                  ensemble.pointWeights.allSatisfy({ $0.isFinite && $0 > 0 }),
                  abs(ensemble.pointWeights.reduce(0,+)-1) < 1e-12,
                  ensemble.minimumOccupationBoundaryGap.isFinite, ensemble.minimumOccupationBoundaryGap > 0 else {
                throw VivoChemistryError.invalid("ensemble density weights or occupation-gap contract")
            }
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
                guard 2*cfg.maximumPointEvaluations <= maximumPointEvaluations - reservation else {
                    throw VivoChemistryError.resourceLimit("campaign ECC point-evaluation reservation")
                }
                reservation += 2*cfg.maximumPointEvaluations
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
    public let ensembleOccupations: [Double]?
    public let ensembleRotation: VivoQMMatrix?
    public let referenceEnergiesHartree: [Double]
    public let referenceResiduals: [Double]
    public let referenceDifferences: VivoBarrierDifferences
    public let levels: [VivoBarrierLevelReport]
    public let assessment: VivoBarrierAssessment
    public let acceptedReducedLevelIdentifier: String?
    public let chargedPointEvaluations: Int
    public let meaning: String
}

