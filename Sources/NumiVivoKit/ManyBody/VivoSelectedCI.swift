import Foundation

/// Existing selected-CI request fields are retained. The two optional guards
/// decode older requests, but new results use a stricter, versioned contract.
public struct VivoSelectedCIConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var maximumDeterminants: Int
    public var selectionBatchSize: Int
    public var minimumSelectionContributionHartree: Double
    public var pt2ToleranceHartree: Double
    public var eigenResidualTolerance: Double
    public var minimumDenominatorHartree: Double
    public var maximumDavidsonSubspace: Int
    public var maximumExternalDeterminants: Int?
    public var fullResidualTolerance: Double?
    public init(maximumIterations: Int = 32, maximumDeterminants: Int = 512,
                selectionBatchSize: Int = 128, minimumSelectionContributionHartree: Double = 1e-9,
                pt2ToleranceHartree: Double = 1e-6, eigenResidualTolerance: Double = 1e-9,
                minimumDenominatorHartree: Double = 1e-5, maximumDavidsonSubspace: Int = 48,
                maximumExternalDeterminants: Int? = nil, fullResidualTolerance: Double? = nil) {
        self.maximumIterations = maximumIterations; self.maximumDeterminants = maximumDeterminants
        self.selectionBatchSize = selectionBatchSize; self.minimumSelectionContributionHartree = minimumSelectionContributionHartree
        self.pt2ToleranceHartree = pt2ToleranceHartree; self.eigenResidualTolerance = eigenResidualTolerance
        self.minimumDenominatorHartree = minimumDenominatorHartree; self.maximumDavidsonSubspace = maximumDavidsonSubspace
        self.maximumExternalDeterminants = maximumExternalDeterminants; self.fullResidualTolerance = fullResidualTolerance
    }
    public var effectiveFullResidualTolerance: Double { fullResidualTolerance ?? max(1e-7,10*eigenResidualTolerance) }
    public func externalCapacity(budget: VivoChemistryBudget) -> Int {
        maximumExternalDeterminants ?? min(50_000,max(1,budget.maximumBytes/4096))
    }
    public func validate(budget: VivoChemistryBudget) throws {
        try budget.validate()
        guard (1...1000).contains(maximumIterations), (1...budget.maximumDeterminants).contains(maximumDeterminants),
              (1...maximumDeterminants).contains(selectionBatchSize),
              minimumSelectionContributionHartree.isFinite, minimumSelectionContributionHartree >= 0,
              pt2ToleranceHartree.isFinite, pt2ToleranceHartree > 0,
              eigenResidualTolerance.isFinite, eigenResidualTolerance > 0,
              minimumDenominatorHartree.isFinite, minimumDenominatorHartree > 0,
              (4...512).contains(maximumDavidsonSubspace), externalCapacity(budget: budget) > 0,
              effectiveFullResidualTolerance.isFinite, effectiveFullResidualTolerance > 0,
              eigenResidualTolerance <= effectiveFullResidualTolerance/10 else {
            throw VivoChemistryError.invalid("selected-CI capacities, thresholds or inner/full residual contract")
        }
    }
    var davidson: VivoDavidsonConfiguration {
        .init(maximumIterations: max(150,4*maximumIterations),maximumSubspace: maximumDavidsonSubspace,
              residualTolerance: eigenResidualTolerance)
    }
}
public enum VivoSelectedCITermination: String, Codable, Sendable {
    case residualConverged, determinantLimit, iterationLimit, selectionThreshold
}
public struct VivoSelectedCIExternalContribution: Codable, Sendable, Equatable {
    public let determinant: UInt64
    /// <D|H|Psi> after summing all connected contributions, not max |H_Di c_i|.
    public let residualAmplitudeHartree: Double
    public let denominatorHartree: Double
    public let intruder: Bool
    public let selectionAmplitude: Double
    public let selectionContributionHartree: Double
}
public struct VivoSelectedCIDiagnostics: Codable, Sendable, Equatable {
    public let projectedResidualNorm: Double
    public let externalResidualNorm: Double
    public var fullResidualNorm: Double { hypot(projectedResidualNorm,externalResidualNorm) }
    /// Complete EN2 sum, or nil if any nonzero external contribution has a
    /// nonnegative/near-zero denominator. This estimate is never an error bound.
    public let epsteinNesbetCorrectionHartree: Double?
    public let intruderCount: Int
    public let externalDeterminantCount: Int
    public let operatorApplications: Int
    public let externalContributions: [VivoSelectedCIExternalContribution]
}
public struct VivoSelectedCIIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let determinantCount: Int
    public let variationalEnergyHartree: Double
    public let eigenResidual: Double
    public let externalCandidateCount: Int
    public let selectedCount: Int
    public let pt2CorrectionHartree: Double?
    public let largestExternalContributionHartree: Double
    public let intruderCandidateCount: Int
    public let fullResidualNorm: Double
    public let addedDeterminants: [UInt64]
}
public struct VivoSelectedCIResult: Codable, Sendable, Equatable {
    public let schema: String
    public let converged: Bool
    public let variationalEnergyHartree: Double
    /// Nil when a valid unregularized EN2 correction is unavailable.
    public let pt2CorrectedEnergyHartree: Double?
    public let pt2CorrectionHartree: Double?
    public let state: VivoCIState
    public let fullSectorDimension: Int
    public let selectedDeterminantCount: Int
    public let eigenResidual: Double
    public let iterations: [VivoSelectedCIIteration]
    public let operatorApplications: Int
    public let method: String
    public let configuration: VivoSelectedCIConfiguration
    public let seedDeterminants: [UInt64]
    public let diagnostics: VivoSelectedCIDiagnostics
    public let termination: VivoSelectedCITermination
    public let matrixVectorProducts: Int
}
/// Consolidated deterministic CIPSI-style selected CI. The existing Davidson
/// and Slater-Condon implementations serve both full and selected spaces.
/// PT2 is a diagnostic, never an error bound or global-ground-root certificate.
public enum VivoSelectedCI {
    public static let method = "numivivo.selected-ci.v2: deterministic residual/EN2 selection with shared Davidson; PT2 is a diagnostic remainder, not a certified error bound; full connected residual required"
    private static func sectorDimension(_ n: Int, _ a: Int, _ b: Int) throws -> Int {
        func choose(_ population: Int) throws -> Int {
            let k = min(population,n-population); var c = 1
            if k > 0 { for i in 1...k {
                let x = c.multipliedReportingOverflow(by: n-k+i)
                guard !x.overflow else { throw VivoChemistryError.resourceLimit("selected-CI sector dimension overflow") }
                c = x.partialValue/i
            } }
            return c
        }
        let ca = try choose(a), cb = try choose(b), x = ca.multipliedReportingOverflow(by: cb)
        guard !x.overflow else { throw VivoChemistryError.resourceLimit("selected-CI sector dimension overflow") }
        return x.partialValue
    }
    private static func defaultSeed(_ h: VivoEmbeddedHamiltonian) -> UInt64 {
        var d: UInt64 = 0
        for p in 0..<h.alphaElectrons { d |= UInt64(1) << (2*p) }
        for p in 0..<h.betaElectrons { d |= UInt64(1) << (2*p+1) }
        return d
    }
    /// Conservative accounting for live sparse maps, sort buffers, occupied/
    /// virtual lists and result records; not an OS resident-memory measurement.
    private static func reservedBytes(external: Int, selected: Int, budget: VivoChemistryBudget) throws -> Int {
        let a = external.multipliedReportingOverflow(by: 1024)
        let b = selected.multipliedReportingOverflow(by: 2048)
        let sum = a.partialValue.addingReportingOverflow(b.partialValue)
        guard !a.overflow, !b.overflow, !sum.overflow, sum.partialValue < budget.maximumBytes else {
            throw VivoChemistryError.resourceLimit("selected CI aggregate sparse workspace")
        }
        return sum.partialValue
    }
    /// Recomputes the complete connected residual without enumerating the sector.
    /// Capacity exhaustion throws rather than publishing a screened 'full' norm.
    public static func diagnose(_ h: VivoEmbeddedHamiltonian, state: VivoCIState, energyHartree: Double,
                                configuration cfg: VivoSelectedCIConfiguration = .init(),
                                budget: VivoChemistryBudget = .init()) throws -> VivoSelectedCIDiagnostics {
        try cfg.validate(budget: budget); try h.validate(budget: budget); try state.validate(budget: budget)
        guard state.orbitalCount == h.orbitalCount, state.alphaElectrons == h.alphaElectrons,
              state.betaElectrons == h.betaElectrons, energyHartree.isFinite else {
            throw VivoChemistryError.invalid("selected CI diagnostic Hamiltonian/state binding")
        }
        _ = try reservedBytes(external: cfg.externalCapacity(budget: budget), selected: state.determinants.count, budget: budget)
        let action = try VivoDirectHamiltonian(h,determinants: state.determinants,budget: budget)
        var projected = [Double](repeating: 0,count: state.determinants.count)
        var external: [UInt64: Double] = [:], work = 0
        // Ascending determinant order fixes accumulation order. No small
        // coefficient is discarded from a purported complete residual.
        for k in state.determinants.indices.sorted(by: { state.determinants[$0] < state.determinants[$1] }) where state.coefficients[k] != 0 {
            try action.connections(from: k,work: &work) { destination,value in
                let contribution = value*state.coefficients[k]
                guard contribution.isFinite else { throw VivoChemistryError.convergence("selected CI residual overflow") }
                if let i = action.index[destination] { projected[i] += contribution }
                else if contribution != 0 {
                    if external[destination] == nil, external.count >= cfg.externalCapacity(budget: budget) {
                        throw VivoChemistryError.resourceLimit("selected CI external frontier capacity; residual is not complete")
                    }
                    external[destination,default: 0] += contribution
                }
            }
        }
        let electronic = energyHartree-h.constantEnergyHartree
        var inside = 0.0, outside = 0.0, correction = 0.0, compensation = 0.0, intruders = 0
        for i in projected.indices { inside = hypot(inside,projected[i]-electronic*state.coefficients[i]) }
        var contributions: [VivoSelectedCIExternalContribution] = []
        for d in external.keys.sorted() {
            let r = external[d]!
            if r == 0 { continue }
            guard work < budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("selected CI diagnostic work") }
            work += 1
            let denominator = electronic-action.diagonalElement(d)
            let intruder = denominator >= -cfg.minimumDenominatorHartree
            let amplitude = abs(r)/max(abs(denominator),cfg.minimumDenominatorHartree)
            let selection = amplitude*abs(r)
            guard r.isFinite, denominator.isFinite, amplitude.isFinite, selection.isFinite else {
                throw VivoChemistryError.convergence("selected CI diagnostic nonfinite")
            }
            outside = hypot(outside,r)
            if intruder { intruders += 1 }
            else {
                let term = (r/denominator)*r
                guard term.isFinite else { throw VivoChemistryError.convergence("selected CI EN2 overflow") }
                let y = term-compensation, t = correction+y
                compensation = (t-correction)-y; correction = t
            }
            contributions.append(.init(determinant: d,residualAmplitudeHartree: r,denominatorHartree: denominator,
                intruder: intruder,selectionAmplitude: amplitude,selectionContributionHartree: selection))
        }
        guard inside.isFinite, outside.isFinite, correction.isFinite else { throw VivoChemistryError.convergence("selected CI norm overflow") }
        return .init(projectedResidualNorm: inside,externalResidualNorm: outside,
            epsteinNesbetCorrectionHartree: intruders == 0 ? correction : nil,
            intruderCount: intruders,externalDeterminantCount: contributions.count,
            operatorApplications: work,externalContributions: contributions)
    }
    public static func solve(_ h: VivoEmbeddedHamiltonian, configuration cfg: VivoSelectedCIConfiguration = .init(),
                             initialDeterminants suppliedSeeds: [UInt64]? = nil,
                             budget: VivoChemistryBudget = .init()) throws -> VivoSelectedCIResult {
        try cfg.validate(budget: budget); try h.validate(budget: budget)
        guard h.orbitalCount <= 31 else { throw VivoChemistryError.resourceLimit("selected CI supports at most 31 spatial orbitals") }
        let seeds = (suppliedSeeds ?? [defaultSeed(h)]).sorted()
        guard !seeds.isEmpty, seeds.count <= cfg.maximumDeterminants, Set(seeds).count == seeds.count else {
            throw VivoChemistryError.invalid("selected CI seed set")
        }
        var dets = seeds, history: [VivoSelectedCIIteration] = [], work = 0, products = 0
        var previousEnergy: Double?
        for iteration in 1...cfg.maximumIterations {
            let reserved = try reservedBytes(external: cfg.externalCapacity(budget: budget), selected: dets.count, budget: budget)
            var innerBudget = budget
            innerBudget.maximumBytes -= reserved; innerBudget.maximumOperatorApplications -= work
            guard innerBudget.maximumOperatorApplications > 0 else { throw VivoChemistryError.resourceLimit("selected CI aggregate work") }
            let eig = try VivoDirectCI.diagonalize(h,determinants: dets,configuration: cfg.davidson,budget: innerBudget)
            work += eig.operatorApplications; products += eig.matrixVectorProducts
            let state = eig.states[0], energy = eig.energiesHartree[0]
            if let previousEnergy, energy > previousEnergy+max(1e-10,10*cfg.davidson.residualTolerance) {
                throw VivoChemistryError.convergence("selected CI variational energy increased on nested expansion")
            }
            var diagnosticBudget = budget; diagnosticBudget.maximumOperatorApplications -= work
            guard diagnosticBudget.maximumOperatorApplications > 0 else { throw VivoChemistryError.resourceLimit("selected CI aggregate diagnostic work") }
            let diagnostics = try diagnose(h,state: state,energyHartree: energy,configuration: cfg,budget: diagnosticBudget)
            work += diagnostics.operatorApplications
            let remaining = cfg.maximumDeterminants-dets.count
            let ranked = diagnostics.externalContributions.filter {
                $0.intruder || $0.selectionContributionHartree >= cfg.minimumSelectionContributionHartree
            }.sorted {
                if $0.intruder != $1.intruder { return $0.intruder }
                if $0.selectionContributionHartree != $1.selectionContributionHartree {
                    return $0.selectionContributionHartree > $1.selectionContributionHartree
                }
                return $0.determinant < $1.determinant
            }
            var termination: VivoSelectedCITermination?
            if diagnostics.fullResidualNorm <= cfg.effectiveFullResidualTolerance,
               let pt2 = diagnostics.epsteinNesbetCorrectionHartree, abs(pt2) <= cfg.pt2ToleranceHartree { termination = .residualConverged }
            else if remaining == 0 { termination = .determinantLimit }
            else if iteration == cfg.maximumIterations { termination = .iterationLimit }
            else if ranked.isEmpty { termination = .selectionThreshold }
            let added = termination == nil ? Array(ranked.prefix(min(remaining,cfg.selectionBatchSize))).map(\.determinant) : []
            history.append(.init(iteration: iteration,determinantCount: dets.count,variationalEnergyHartree: energy,
                eigenResidual: diagnostics.projectedResidualNorm,externalCandidateCount: diagnostics.externalDeterminantCount,
                selectedCount: added.count,pt2CorrectionHartree: diagnostics.epsteinNesbetCorrectionHartree,
                largestExternalContributionHartree: diagnostics.externalContributions.map(\.selectionContributionHartree).max() ?? 0,
                intruderCandidateCount: diagnostics.intruderCount,fullResidualNorm: diagnostics.fullResidualNorm,
                addedDeterminants: added))
            if let termination {
                return .init(schema: "numivivo.selected-ci.v2",converged: termination == .residualConverged,
                    variationalEnergyHartree: energy,
                    pt2CorrectedEnergyHartree: diagnostics.epsteinNesbetCorrectionHartree.map { energy+$0 },
                    pt2CorrectionHartree: diagnostics.epsteinNesbetCorrectionHartree,state: state,
                    fullSectorDimension: try sectorDimension(h.orbitalCount,h.alphaElectrons,h.betaElectrons),
                    selectedDeterminantCount: dets.count,eigenResidual: diagnostics.projectedResidualNorm,
                    iterations: history,operatorApplications: work,method: method,configuration: cfg,
                    seedDeterminants: seeds,diagnostics: diagnostics,termination: termination,matrixVectorProducts: products)
            }
            previousEnergy = energy; dets = (dets+added).sorted()
        }
        throw VivoChemistryError.convergence("unreachable selected CI termination")
    }
    /// Recompute eigenpair evidence, not global-ground-root certification.
    /// Workflow receipts separately bind the exact request and implementation.
    public static func validate(_ result: VivoSelectedCIResult, hamiltonian h: VivoEmbeddedHamiltonian,
                                budget: VivoChemistryBudget = .init()) throws {
        guard result.schema == "numivivo.selected-ci.v2", result.method == method, !result.iterations.isEmpty,
              result.iterations.count <= result.configuration.maximumIterations,
              result.operatorApplications > 0, result.matrixVectorProducts > 0,
              result.seedDeterminants == result.seedDeterminants.sorted(), !result.seedDeterminants.isEmpty,
              Set(result.seedDeterminants).count == result.seedDeterminants.count else {
            throw VivoChemistryError.invalid("selected CI result schema/provenance")
        }
        var reconstructed = result.seedDeterminants
        for (index,step) in result.iterations.enumerated() {
            guard step.iteration == index+1, step.determinantCount == reconstructed.count,
                  step.variationalEnergyHartree.isFinite, step.fullResidualNorm.isFinite, step.fullResidualNorm >= 0,
                  step.selectedCount == step.addedDeterminants.count,
                  step.addedDeterminants.count <= result.configuration.selectionBatchSize,
                  Set(step.addedDeterminants).count == step.addedDeterminants.count,
                  Set(step.addedDeterminants).isDisjoint(with: Set(reconstructed)),
                  (index+1 == result.iterations.count ? step.addedDeterminants.isEmpty : !step.addedDeterminants.isEmpty) else {
                throw VivoChemistryError.invalid("selected CI nested iteration record")
            }
            reconstructed = (reconstructed+step.addedDeterminants).sorted()
        }
        guard reconstructed == result.state.determinants,
              result.iterations.last?.variationalEnergyHartree == result.variationalEnergyHartree,
              result.iterations.last?.fullResidualNorm == result.diagnostics.fullResidualNorm else {
            throw VivoChemistryError.invalid("selected CI final iteration binding")
        }
        let checked = try diagnose(h,state: result.state,energyHartree: result.variationalEnergyHartree,
                                   configuration: result.configuration,budget: budget)
        func close(_ a: Double, _ b: Double) -> Bool { a.isFinite && b.isFinite && abs(a-b) <= 1e-11*max(1,abs(a),abs(b)) }
        let residualConverged = checked.fullResidualNorm <= result.configuration.effectiveFullResidualTolerance &&
            (checked.epsteinNesbetCorrectionHartree.map { abs($0) <= result.configuration.pt2ToleranceHartree } ?? false)
        guard close(checked.fullResidualNorm,result.diagnostics.fullResidualNorm),
              close(checked.projectedResidualNorm,result.diagnostics.projectedResidualNorm),
              checked.intruderCount == result.diagnostics.intruderCount,
              checked.externalDeterminantCount == result.diagnostics.externalDeterminantCount,
              checked.operatorApplications == result.diagnostics.operatorApplications,
              checked.externalContributions == result.diagnostics.externalContributions,
              checked.epsteinNesbetCorrectionHartree == result.diagnostics.epsteinNesbetCorrectionHartree,
              result.converged == (result.termination == .residualConverged), result.converged == residualConverged,
              result.selectedDeterminantCount == result.state.determinants.count,
              result.fullSectorDimension == (try sectorDimension(h.orbitalCount,h.alphaElectrons,h.betaElectrons)),
              result.pt2CorrectionHartree == checked.epsteinNesbetCorrectionHartree,
              result.pt2CorrectedEnergyHartree == checked.epsteinNesbetCorrectionHartree.map({ result.variationalEnergyHartree+$0 }),
              close(result.eigenResidual,checked.projectedResidualNorm),
              checked.projectedResidualNorm <= 1.1*result.configuration.eigenResidualTolerance else {
            throw VivoChemistryError.invalid("selected CI result fails independent residual reconstruction")
        }
    }
}
