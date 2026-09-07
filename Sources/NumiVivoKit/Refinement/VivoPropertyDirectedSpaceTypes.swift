import Foundation

public struct VivoOrbitalRefinementBlock: Codable, Sendable, Equatable {
    public let identifier: String
    public let orbitals: [Int]
    public init(identifier: String, orbitals: [Int]) { self.identifier = identifier; self.orbitals = orbitals }
}
public enum VivoSpaceRefinementSolver: Codable, Sendable, Equatable {
    case directCI(configuration: VivoDavidsonConfiguration)
    case selectedCI(configuration: VivoSelectedCIConfiguration)
    case multistateCI(configuration: VivoDavidsonConfiguration, states: VivoRefinementStatePolicy)
    case stateAveragedCASSCF(configuration: VivoMultiStateCASSCFConfiguration, states: VivoRefinementStatePolicy)
    public var statePolicy: VivoRefinementStatePolicy {
        switch self {
        case .directCI, .selectedCI: return .groundState
        case .multistateCI(_, let states), .stateAveragedCASSCF(_, let states): return states
        }
    }
    public var optimizesOrbitals: Bool {
        if case .stateAveragedCASSCF = self { return true }; return false
    }
}
public struct VivoElectronicProfileTarget: Codable, Sendable, Equatable {
    public let barrierPointIdentifier: String
    public let maximumBarrierShiftHartree: Double
    public let maximumRelativeProfileShiftHartree: Double
    /// V2 only. Bounds every included pairwise electronic excitation-gap change.
    public let maximumStateGapShiftHartree: Double?
    public init(barrierPointIdentifier: String, maximumBarrierShiftHartree: Double = 1e-4,
                maximumRelativeProfileShiftHartree: Double = 1e-4,
                maximumStateGapShiftHartree: Double? = nil) {
        self.barrierPointIdentifier = barrierPointIdentifier
        self.maximumBarrierShiftHartree = maximumBarrierShiftHartree
        self.maximumRelativeProfileShiftHartree = maximumRelativeProfileShiftHartree
        self.maximumStateGapShiftHartree = maximumStateGapShiftHartree
    }
}
/// Fixed-geometry exploration, not an activation Gibbs free energy or a rate.
/// The first and last points define the two endpoints. The declared interior
/// barrier point is an evaluation point, not a saddle characterization.
public struct VivoPropertyDirectedSpaceRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/property-directed-space/v2"
    public static let legacySchema = "numivivo.org/property-directed-space/v1"
    public let schema: String
    public let identifier: String
    public let points: [VivoECCPathPoint]
    public let transportGroups: [[Int]]
    public let initialSpace: VivoActiveSpace
    public let mandatoryActiveOrbitals: [Int]
    public let candidateBlocks: [VivoOrbitalRefinementBlock]
    public let discoveryPointIdentifiers: [String]
    public let confirmationPointIdentifiers: [String]
    public let target: VivoElectronicProfileTarget
    public let solver: VivoSpaceRefinementSolver
    public let maximumRounds: Int
    public let maximumActiveOrbitals: Int
    public let maximumPointEvaluations: Int
    public let minimumTransportSingularValue: Double
    public let minimumStateOverlapSquared: Double
    public let collectOrbitalInformation: Bool
    public let budget: VivoChemistryBudget
    public init(identifier: String, points: [VivoECCPathPoint], transportGroups: [[Int]],
                initialSpace: VivoActiveSpace, mandatoryActiveOrbitals: [Int],
                candidateBlocks: [VivoOrbitalRefinementBlock], discoveryPointIdentifiers: [String],
                confirmationPointIdentifiers: [String], target: VivoElectronicProfileTarget,
                solver: VivoSpaceRefinementSolver = .directCI(configuration: .init()),
                maximumRounds: Int = 8, maximumActiveOrbitals: Int = 31,
                maximumPointEvaluations: Int = 512, minimumTransportSingularValue: Double = 0.5,
                minimumStateOverlapSquared: Double = 0.1, collectOrbitalInformation: Bool = false,
                budget: VivoChemistryBudget = .init()) {
        schema = Self.schema; self.identifier = identifier; self.points = points; self.transportGroups = transportGroups
        self.initialSpace = initialSpace; self.mandatoryActiveOrbitals = mandatoryActiveOrbitals
        self.candidateBlocks = candidateBlocks; self.discoveryPointIdentifiers = discoveryPointIdentifiers
        self.confirmationPointIdentifiers = confirmationPointIdentifiers; self.target = target; self.solver = solver
        self.maximumRounds = maximumRounds; self.maximumActiveOrbitals = maximumActiveOrbitals
        self.maximumPointEvaluations = maximumPointEvaluations; self.minimumTransportSingularValue = minimumTransportSingularValue
        self.minimumStateOverlapSquared = minimumStateOverlapSquared; self.collectOrbitalInformation = collectOrbitalInformation
        self.budget = budget
    }
    public func validate() throws {
        try budget.validate()
        guard [Self.schema, Self.legacySchema].contains(schema), !identifier.isEmpty, identifier.utf8.count <= 1024,
              (4...32).contains(points.count), (0...32).contains(candidateBlocks.count),
              (1...32).contains(maximumRounds), (1...31).contains(maximumActiveOrbitals),
              (1...100_000).contains(maximumPointEvaluations),
              [minimumTransportSingularValue,minimumStateOverlapSquared].allSatisfy({ $0.isFinite && $0 > 0 && $0 <= 1 }),
              [target.maximumBarrierShiftHartree,target.maximumRelativeProfileShiftHartree].allSatisfy({ $0.isFinite && $0 > 0 }) else {
            throw VivoChemistryError.invalid("property-directed space identity, dimensions, tolerances or capacities")
        }
        let statePolicy = solver.statePolicy
        try statePolicy.validate()
        if let tolerance = target.maximumStateGapShiftHartree {
            guard tolerance.isFinite, tolerance > 0 else { throw VivoChemistryError.invalid("state-gap tolerance") }
        }
        if statePolicy.labels.count > 1 {
            guard schema == Self.schema, target.maximumStateGapShiftHartree != nil else {
                throw VivoChemistryError.invalid("multistate refinement requires a v2 request and an explicit gap tolerance")
            }
        }
        if schema == Self.legacySchema {
            guard !solver.optimizesOrbitals, statePolicy.labels.count == 1,
                  target.maximumStateGapShiftHartree == nil, !candidateBlocks.isEmpty else {
                throw VivoChemistryError.invalid("v2 refinement features in a v1 request")
            }
        }
        let first = points[0].hamiltonian, n = first.orbitalCount
        try first.validate(budget: budget); try initialSpace.validate(for: first,budget: budget)
        guard !candidateBlocks.isEmpty || initialSpace.active.count == n else {
            throw VivoChemistryError.invalid("an empty candidate universe is allowed only when the entire orbital space is active")
        }
        let pointIDs = points.map(\.identifier), discovery = Set(discoveryPointIdentifiers), confirmation = Set(confirmationPointIdentifiers)
        guard n <= 31, Set(pointIDs).count == pointIDs.count,
              points.allSatisfy({ !$0.identifier.isEmpty && $0.identifier.utf8.count <= 1024 }),
              points[0].overlapWithPrevious == nil,
              discovery.count == discoveryPointIdentifiers.count, confirmation.count == confirmationPointIdentifiers.count,
              discovery.count >= 3, !confirmation.isEmpty, discovery.isDisjoint(with: confirmation),
              discovery.union(confirmation) == Set(pointIDs), discovery.contains(pointIDs[0]), discovery.contains(pointIDs.last!),
              discovery.contains(target.barrierPointIdentifier), pointIDs.dropFirst().dropLast().contains(target.barrierPointIdentifier),
              initialSpace.frozenOrbitals.isEmpty, initialSpace.active == initialSpace.active.sorted(),
              initialSpace.doublyOccupiedCore == initialSpace.doublyOccupiedCore.sorted(),
              initialSpace.active.count <= maximumActiveOrbitals,
              Set(mandatoryActiveOrbitals).count == mandatoryActiveOrbitals.count,
              Set(mandatoryActiveOrbitals).isSubset(of: Set(initialSpace.active)) else {
            throw VivoChemistryError.invalid("property-directed discovery/confirmation split, active space or mandatory seeds")
        }
        let flattened = transportGroups.flatMap { $0 }
        guard transportGroups.allSatisfy({ !$0.isEmpty }), flattened.count == n, Set(flattened) == Set(0..<n) else {
            throw VivoChemistryError.invalid("complete disjoint orbital-transport groups required")
        }
        let active = Set(initialSpace.active), core = Set(initialSpace.doublyOccupiedCore)
        for group in transportGroups {
            let g = Set(group)
            guard g.isSubset(of: active) || g.isSubset(of: core) || g.isDisjoint(with: active.union(core)) else {
                throw VivoChemistryError.invalid("initial space splits a transported orbital subspace")
            }
        }
        guard Set(candidateBlocks.map(\.identifier)).count == candidateBlocks.count,
              Set(candidateBlocks.map(\.orbitals)).count == candidateBlocks.count else {
            throw VivoChemistryError.invalid("duplicate orbital-refinement blocks")
        }
        for block in candidateBlocks {
            let b = Set(block.orbitals)
            guard !block.identifier.isEmpty, block.identifier.utf8.count <= 1024, !block.identifier.hasPrefix("__"),
                  !block.orbitals.isEmpty, block.orbitals.count <= n, block.orbitals == block.orbitals.sorted(),
                  b.count == block.orbitals.count, b.allSatisfy({ (0..<n).contains($0) }),
                  !b.isSubset(of: active) else { throw VivoChemistryError.invalid("orbital-refinement block") }
            for group in transportGroups {
                let g = Set(group)
                guard g.isDisjoint(with: b) || g.isSubset(of: b) else {
                    throw VivoChemistryError.invalid("refinement must promote whole transported subspaces")
                }
            }
        }
        let accuracy = min(target.maximumBarrierShiftHartree,target.maximumRelativeProfileShiftHartree,
                           target.maximumStateGapShiftHartree ?? Double.greatestFiniteMagnitude)/100
        switch solver {
        case .directCI(let cfg):
            try VivoDirectCI.validate(cfg)
            guard cfg.roots == 1, cfg.residualTolerance <= accuracy else {
                throw VivoChemistryError.invalid("space exploration currently requires one tightly solved CI root")
            }
        case .selectedCI(let cfg):
            try cfg.validate(budget: budget)
            guard cfg.effectiveFullResidualTolerance <= accuracy, cfg.pt2ToleranceHartree <= accuracy else {
                throw VivoChemistryError.invalid("selected-CI probe tolerances are too loose for the requested property shifts")
            }
        case .multistateCI(let cfg, let states):
            try VivoDirectCI.validate(cfg)
            guard cfg.roots == states.labels.count, cfg.residualTolerance <= accuracy else {
                throw VivoChemistryError.invalid("multistate CI root count or accuracy")
            }
        case .stateAveragedCASSCF(let cfg, let states):
            try cfg.optimization.validate(); try VivoDirectCI.validate(cfg.davidson)
            guard schema == Self.schema, cfg.rootLabels == states.labels,
                  cfg.weights.count == states.labels.count, cfg.davidson.roots == states.labels.count,
                  cfg.weights.allSatisfy({ $0.isFinite && $0 >= 0 }), abs(cfg.weights.reduce(0,+)-1) < 1e-12,
                  !cfg.followRoots, cfg.optimization.energyToleranceHartree <= accuracy,
                  cfg.davidson.residualTolerance <= accuracy else {
                throw VivoChemistryError.invalid("refinement SA-CASSCF requires tightly converged energy-ordered roots and explicit weights")
            }
            for group in states.groups {
                guard group.allSatisfy({ abs(cfg.weights[$0] - cfg.weights[group[0]]) < 1e-12 }) else {
                    throw VivoChemistryError.invalid("weights must be equal inside an interchangeable state subspace")
                }
            }
        }
        _ = try budget.elements([points.count,n,n,n,n],simultaneousArrays: 8)
        _ = try budget.elements([points.count,budget.maximumDeterminants,16,statePolicy.labels.count],simultaneousArrays: 6)
        for (i,point) in points.enumerated() {
            let h = point.hamiltonian; try h.validate(budget: budget)
            guard h.orbitalIdentifiers == first.orbitalIdentifiers, h.alphaElectrons == first.alphaElectrons,
                  h.betaElectrons == first.betaElectrons, h.energyReference == first.energyReference else {
                throw VivoChemistryError.invalid("refinement path changes orbital identity, spin sector or energy convention")
            }
            if i > 0 {
                guard let s = point.overlapWithPrevious, s.rows == n, s.columns == n,
                      s.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("physical cross-geometry orbital overlap required") }
                let spectrum = try VivoQMDenseAlgebra.symmetricEigen(s.transposed.multiplied(by: s))
                guard spectrum.values.allSatisfy({ $0 >= -1e-10 && $0 <= 1+1e-6 }) else {
                    throw VivoChemistryError.invalid("cross-geometry overlap is not a subspace contraction")
                }
            }
        }
    }
}
public struct VivoSpacePointObservation: Codable, Sendable, Equatable {
    public let pointIdentifier: String
    public let energyHartree: Double
    public let determinantCount: Int
    public let ciResidualNorm: Double
    public let pt2CorrectionHartree: Double?
    /// Local active-column indices; the enclosing evaluation's partition maps
    /// them back to the common full orbital frame. Not environmental coverage.
    public let orbitalInformation: VivoSelectiveOrbitalInformationResult?
    /// Nil for legacy single-root evidence. Energies are in adiabatic energy order.
    public var roots: [VivoRefinementRootObservation]? = nil
    /// Relative to the transported input frame; exported anchors must retain it.
    public var optimizedOrbitalRotation: VivoQMMatrix? = nil
}
public struct VivoSpaceEvaluation: Codable, Sendable, Equatable {
    public let partition: VivoActiveSpace
    public let points: [VivoSpacePointObservation]
    public let minimumAdjacentStateOverlapSquared: Double?
    public let assessedAdjacentEdges: Int
    public let hamiltonianOperatorApplications: Int
    public let reservedAuxiliaryWork: Int
}
public struct VivoSpacePropertyShift: Codable, Sendable, Equatable {
    public let forwardHartree: Double
    public let reverseHartree: Double
    public let reactionHartree: Double
    public let maximumRelativeProfileHartree: Double
    public let maximumAbsoluteEnergyHartree: Double
    public let minimumExpansionStateOverlapSquared: Double
    public let normalizedPropertyImpact: Double
    public var rootShifts: [VivoRefinementRootShift]? = nil
    public var maximumStateGapShiftHartree: Double? = nil
}
public struct VivoSpaceRefinementTrial: Codable, Sendable, Equatable {
    public let blockIdentifier: String
    public let evaluation: VivoSpaceEvaluation
    public let shift: VivoSpacePropertyShift
    /// Observed property sensitivity / charged work, not an uncertainty bound.
    public let priority: Double
}
public struct VivoSpaceRefinementRound: Codable, Sendable, Equatable {
    public let round: Int
    public let baseline: VivoSpaceEvaluation
    public let trials: [VivoSpaceRefinementTrial]
    public let selectedBlockIdentifier: String?
}
public enum VivoSpaceRefinementTermination: String, Codable, Sendable {
    case sensitivityStable, candidatePoolIncluded, confirmationFailed
    case roundLimit, activeSpaceLimit, incompleteExecution
}
public struct VivoSpaceConfirmation: Codable, Sendable, Equatable {
    public let chosen: VivoSpaceEvaluation
    public let candidateUnion: VivoSpaceEvaluation
    public let shift: VivoSpacePropertyShift
    public let passed: Bool
}
public struct VivoPropertyDirectedSpaceResult: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/property-directed-space-result/v2"
    public let schema: String
    public let request: VivoPropertyDirectedSpaceRequest
    public let transportRotations: [VivoQMMatrix]
    public let transportMinimumSingularValues: [Double]
    public let rounds: [VivoSpaceRefinementRound]
    public let finalSpace: VivoActiveSpace
    public let candidatePoolOrbitals: [Int]
    public let unexaminedOrbitals: [Int]
    public let confirmation: VivoSpaceConfirmation?
    public let termination: VivoSpaceRefinementTermination
    public let failure: String?
    public let pointEvaluations: Int
    public let hamiltonianOperatorApplications: Int
    public let reservedAuxiliaryWork: Int
    public let reservedFailedSolverWork: Int
    public let meaning: String
    public var sensitivityEstablishedWithinDeclaredPool: Bool {
        (termination == .sensitivityStable || termination == .candidatePoolIncluded) && confirmation?.passed == true
    }
}
