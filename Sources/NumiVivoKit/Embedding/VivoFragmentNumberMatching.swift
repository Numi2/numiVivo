import Foundation

/// A fragment number operator in the cluster's orthonormal spatial basis.
/// Bath occupations are not counted. Different clusters may overlap: this type
/// does not by itself define a non-overlapping global density/energy partition.
public struct VivoNumberMatchedFragment: Codable, Sendable, Equatable {
    public let identifier: String
    public let hamiltonian: VivoEmbeddedHamiltonian
    public let fragmentProjector: VivoQMMatrix
    public init(identifier: String, hamiltonian: VivoEmbeddedHamiltonian, fragmentProjector: VivoQMMatrix) {
        self.identifier = identifier; self.hamiltonian = hamiltonian; self.fragmentProjector = fragmentProjector
    }
    public func validate(budget: VivoChemistryBudget = .init()) throws {
        try hamiltonian.validate(budget: budget)
        let p = fragmentProjector, n = hamiltonian.orbitalCount
        guard !identifier.isEmpty, identifier.utf8.count <= 1024, p.rows == n, p.columns == n,
              p.values.allSatisfy(\.isFinite), try p.adding(p.transposed,scale:-1).frobeniusNorm < 1e-9,
              try p.multiplied(by:p).adding(p,scale:-1).frobeniusNorm < 1e-9 else {
            throw VivoChemistryError.invalid("fragment number operator must be a finite symmetric orthogonal projector")
        }
        let rank = (0..<n).reduce(0.0) { $0+p[$1,$1] }
        guard rank >= 1-1e-9, rank <= Double(n)+1e-9 else {
            throw VivoChemistryError.invalid("empty or invalid fragment projector")
        }
    }
}
public struct VivoNumberMatchingConfiguration: Codable, Sendable, Equatable {
    public var populationTolerance: Double
    public var initialChemicalPotentialHartree: Double
    public var initialBracketStepHartree: Double
    public var maximumAbsoluteChemicalPotentialHartree: Double
    public var maximumEvaluations: Int
    public init(populationTolerance: Double = 1e-8, initialChemicalPotentialHartree: Double = 0,
                initialBracketStepHartree: Double = 0.1, maximumAbsoluteChemicalPotentialHartree: Double = 32,
                maximumEvaluations: Int = 128) {
        self.populationTolerance = populationTolerance; self.initialChemicalPotentialHartree = initialChemicalPotentialHartree
        self.initialBracketStepHartree = initialBracketStepHartree
        self.maximumAbsoluteChemicalPotentialHartree = maximumAbsoluteChemicalPotentialHartree
        self.maximumEvaluations = maximumEvaluations
    }
    public func validate() throws {
        guard populationTolerance.isFinite, populationTolerance > 0, populationTolerance < 0.1,
              initialChemicalPotentialHartree.isFinite, initialBracketStepHartree.isFinite,
              initialBracketStepHartree > 0, maximumAbsoluteChemicalPotentialHartree.isFinite,
              maximumAbsoluteChemicalPotentialHartree > 0,
              abs(initialChemicalPotentialHartree) <= maximumAbsoluteChemicalPotentialHartree,
              (1...4096).contains(maximumEvaluations) else {
            throw VivoChemistryError.invalid("fragment number-matching settings")
        }
    }
}
public struct VivoNumberMatchedState: Codable, Sendable, Equatable {
    public let identifier: String
    public let fragmentPopulation: Double
    public let biasedCI: VivoCIResult
    /// <Psi(mu)|H_cluster|Psi(mu)>, with the artificial -mu*N_F removed.
    /// Summing these overlapping cluster energies is NOT a DMET total energy.
    public let physicalClusterEnergyHartree: Double
}
public struct VivoNumberMatchingIteration: Codable, Sendable, Equatable {
    public let chemicalPotentialHartree: Double
    public let fragmentPopulationSum: Double
    public let populationResidual: Double
}
public enum VivoNumberMatchingTermination: String, Codable, Sendable {
    case converged, unbracketed, evaluationLimit, discontinuousOrUnresponsive
}
public struct VivoNumberMatchingResult: Codable, Sendable, Equatable {
    public let targetFragmentPopulationSum: Double
    public let chemicalPotentialHartree: Double
    public let states: [VivoNumberMatchedState]
    public let populationResidual: Double
    public let termination: VivoNumberMatchingTermination
    public var converged: Bool { termination == .converged }
    public let history: [VivoNumberMatchingIteration]
}

public enum VivoFragmentNumberMatcher {
    private struct Evaluation { let mu: Double; let population: Double; let states: [VivoNumberMatchedState] }
    /// Solve sum_A <N_F,A> = target by shifting fragment projectors only:
    /// H_A(mu) = H_A - mu*N_F,A. Each cluster retains its fixed electron sector.
    /// Exact CI gives a monotone population response; a jump/plateau is not
    /// accepted merely because the chemical-potential bracket becomes narrow.
    /// This is a multi-fragment building block, not a full correlation-potential
    /// DMET cycle. The paper's single-fragment workflow need not invoke it.
    public static func solve(fragments: [VivoNumberMatchedFragment], targetPopulation: Double,
                             configuration cfg: VivoNumberMatchingConfiguration = .init(),
                             budget: VivoChemistryBudget = .init()) throws -> VivoNumberMatchingResult {
        try cfg.validate(); try budget.validate()
        guard (1...64).contains(fragments.count), Set(fragments.map(\.identifier)).count == fragments.count,
              targetPopulation.isFinite, targetPopulation >= 0 else {
            throw VivoChemistryError.invalid("number-matching fragment identity, count or target")
        }
        var maximumPopulation = 0.0
        for fragment in fragments {
            try fragment.validate(budget: budget)
            let h = fragment.hamiltonian
            let rank = (0..<h.orbitalCount).reduce(0.0) { $0+fragment.fragmentProjector[$1,$1] }.rounded()
            maximumPopulation += min(rank,Double(h.alphaElectrons))+min(rank,Double(h.betaElectrons))
        }
        guard targetPopulation <= maximumPopulation+cfg.populationTolerance else {
            throw VivoChemistryError.invalid("target exceeds fragment capacity in the selected cluster spin sectors")
        }
        // A global evaluation bound complements each local solver's work budget.
        // CI/RDM work limits are per fragment solve, not advertised as total work.
        var history: [VivoNumberMatchingIteration] = []
        func evaluate(_ mu: Double) throws -> Evaluation {
            guard history.count < cfg.maximumEvaluations else { throw VivoChemistryError.resourceLimit("fragment number evaluations") }
            var states: [VivoNumberMatchedState] = [], population = 0.0
            for fragment in fragments {
                let h = fragment.hamiltonian, n = h.orbitalCount
                let shifted = VivoEmbeddedHamiltonian(orbitalIdentifiers:h.orbitalIdentifiers,
                    alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
                    oneElectron:try h.oneElectron.adding(fragment.fragmentProjector,scale:-mu),
                    twoElectron:h.twoElectron,constantEnergyHartree:h.constantEnergyHartree,
                    energyReference:h.energyReference+"; auxiliary fragment chemical potential",provenance:h.provenance)
                let ci = try VivoDirectCI.solve(shifted,budget:budget).roots[0]
                let state = ci.state
                let index = Dictionary(uniqueKeysWithValues:state.determinants.enumerated().map { ($0.element,$0.offset) })
                let (work,overflow) = (2*n*n).multipliedReportingOverflow(by:state.determinants.count)
                guard !overflow, work <= budget.maximumOperatorApplications else {
                    throw VivoChemistryError.resourceLimit("fragment one-RDM contraction work")
                }
                var electrons = 0.0
                for p in 0..<n { for q in 0..<n where fragment.fragmentProjector[p,q] != 0 {
                    for spin in 0..<2 {
                        electrons += fragment.fragmentProjector[p,q]*VivoCIDensityMatrices.expectation(state,index,
                            [.init(mode:2*q+spin,creation:false),.init(mode:2*p+spin,creation:true)])
                    }
                } }
                guard electrons.isFinite, electrons >= -1e-8 else { throw VivoChemistryError.convergence("fragment population") }
                let physical = ci.energyHartree+mu*electrons
                guard physical.isFinite else { throw VivoChemistryError.convergence("physical cluster energy overflow") }
                population += electrons
                states.append(.init(identifier:fragment.identifier,fragmentPopulation:electrons,
                                    biasedCI:ci,physicalClusterEnergyHartree:physical))
            }
            guard population.isFinite else { throw VivoChemistryError.convergence("fragment population sum") }
            for previous in history {
                let monotonicityTolerance = max(1e-7,10*cfg.populationTolerance)
                if mu > previous.chemicalPotentialHartree, population < previous.fragmentPopulationSum-monotonicityTolerance {
                    throw VivoChemistryError.convergence("exact ground-state fragment response lost monotonicity")
                }
                if mu < previous.chemicalPotentialHartree, population > previous.fragmentPopulationSum+monotonicityTolerance {
                    throw VivoChemistryError.convergence("exact ground-state fragment response lost monotonicity")
                }
            }
            history.append(.init(chemicalPotentialHartree:mu,fragmentPopulationSum:population,
                                 populationResidual:population-targetPopulation))
            return .init(mu:mu,population:population,states:states)
        }
        var best = try evaluate(cfg.initialChemicalPotentialHartree)
        func result(_ termination: VivoNumberMatchingTermination) -> VivoNumberMatchingResult {
            .init(targetFragmentPopulationSum:targetPopulation,chemicalPotentialHartree:best.mu,states:best.states,
                  populationResidual:best.population-targetPopulation,termination:termination,history:history)
        }
        func keep(_ evaluation: Evaluation) {
            if abs(evaluation.population-targetPopulation) < abs(best.population-targetPopulation) { best = evaluation }
        }
        if abs(best.population-targetPopulation) <= cfg.populationTolerance { return result(.converged) }
        var lower: Evaluation? = best.population < targetPopulation ? best : nil
        var upper: Evaluation? = best.population > targetPopulation ? best : nil
        let direction = best.population < targetPopulation ? 1.0 : -1.0
        var step = cfg.initialBracketStepHartree, previous = best
        while lower == nil || upper == nil {
            if history.count == cfg.maximumEvaluations { return result(.evaluationLimit) }
            let limit = cfg.maximumAbsoluteChemicalPotentialHartree
            let candidate = min(limit,max(-limit,previous.mu+direction*min(step,2*limit)))
            if candidate == previous.mu { return result(.unbracketed) }
            let current = try evaluate(candidate); keep(current)
            if abs(best.population-targetPopulation) <= cfg.populationTolerance { return result(.converged) }
            if current.population < targetPopulation { lower = current } else { upper = current }
            previous = current; step = min(limit,step*2)
        }
        var low = lower!, high = upper!
        while history.count < cfg.maximumEvaluations {
            let midpoint = low.mu/2+high.mu/2
            if midpoint == low.mu || midpoint == high.mu { return result(.discontinuousOrUnresponsive) }
            let current = try evaluate(midpoint); keep(current)
            if abs(best.population-targetPopulation) <= cfg.populationTolerance { return result(.converged) }
            if current.population < targetPopulation { low = current } else { high = current }
        }
        return result(.evaluationLimit)
    }
}
