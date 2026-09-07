import Foundation

/// Property-directed exploration around the existing native CI, orbital-frame
/// transport and state-overlap implementations. All candidate comparisons use
/// exactly matched geometries and a fixed common orbital frame. No production
/// trajectory or Hamiltonian is mutated by this calculation.
public enum VivoPropertyDirectedSpace {
    public static let meaning = "fixed-geometry electronic sensitivity within the declared candidate pool, confirmed on predeclared held-out geometries and against the joint candidate union; not a reference-error bound, global-ground-root certificate, activation Gibbs energy, kinetic rate or production-Hamiltonian update"

    public static func run(_ request: VivoPropertyDirectedSpaceRequest) throws -> VivoPropertyDirectedSpaceResult {
        try request.validate()
        return try Runner(request).execute()
    }
    public static func validate(_ result: VivoPropertyDirectedSpaceResult, request: VivoPropertyDirectedSpaceRequest) throws {
        guard result.schema == VivoPropertyDirectedSpaceResult.schema, result.request == request,
              result.meaning == meaning, result == (try run(request)) else {
            throw VivoChemistryError.invalid("property-directed refinement fails exact bounded reconstruction")
        }
    }
    /// Promote complete groups. A promoted occupied orbital leaves the inactive
    /// core; no occupied core orbital is silently turned into an empty external.
    public static func expanding(_ space: VivoActiveSpace, with orbitals: [Int]) -> VivoActiveSpace {
        let active = Set(space.active).union(orbitals)
        return .init(doublyOccupiedCore: space.doublyOccupiedCore.filter { !active.contains($0) }.sorted(),
                     active: active.sorted(),frozenOrbitals: space.frozenOrbitals)
    }
    private enum Limit: Error { case activeSpace }
    private struct Evaluation {
        let summary: VivoSpaceEvaluation
        let indices: [Int]
        let states: [VivoCIState] // lifted into the full common spatial-orbital frame
    }
    private final class Runner {
        let request: VivoPropertyDirectedSpaceRequest
        let pool: [Int]
        let discovery: [Int]
        var hamiltonians: [VivoEmbeddedHamiltonian] = []
        var rotations: [VivoQMMatrix] = []
        var minima: [Double] = []
        var edges: [VivoQMMatrix] = []
        var pointEvaluations = 0, operatorWork = 0, auxiliaryWork = 0, failedWorkReservation = 0
        var rounds: [VivoSpaceRefinementRound] = []
        var current: VivoActiveSpace
        init(_ request: VivoPropertyDirectedSpaceRequest) {
            self.request = request; current = request.initialSpace
            pool = Set(request.initialSpace.active + request.candidateBlocks.flatMap(\.orbitals)).sorted()
            let ids = Set(request.discoveryPointIdentifiers)
            discovery = request.points.indices.filter { ids.contains(request.points[$0].identifier) }
        }
        var chargedWork: Int { operatorWork+auxiliaryWork+failedWorkReservation }
        func remainingBudget() throws -> VivoChemistryBudget {
            var budget = request.budget
            budget.maximumOperatorApplications -= chargedWork
            guard budget.maximumOperatorApplications > 0 else { throw VivoChemistryError.resourceLimit("refinement aggregate work budget") }
            return budget
        }
        func reserve(_ amount: Int) throws {
            guard amount >= 0, amount <= request.budget.maximumOperatorApplications-chargedWork else {
                throw VivoChemistryError.resourceLimit("refinement aggregate auxiliary-work reservation")
            }
            auxiliaryWork += amount
        }
        func prepare() throws {
            let n = request.points[0].hamiltonian.orbitalCount, identity = try VivoQMMatrix.identity(n)
            rotations = [identity]; minima = [1]; hamiltonians = [request.points[0].hamiltonian]
            for i in 1..<request.points.count {
                // Charged estimates for transport/eigensolves and four-index
                // transformation. Hamiltonian action counts remain separate.
                try reserve(256*n*n*n+8*n*n*n*n*n)
                let s = request.points[i].overlapWithPrevious!
                let alignment = try VivoOrbitalPathTransport.align(referenceCoefficients: rotations[i-1],
                    currentCoefficients: identity,crossAOOverlap: s,groups: request.transportGroups,
                    minimumSingularValue: request.minimumTransportSingularValue)
                rotations.append(alignment.rotation); minima.append(alignment.minimumSubspaceSingularValue)
                edges.append(try rotations[i-1].transposed.multiplied(by: s).multiplied(by: alignment.rotation))
                hamiltonians.append(try request.points[i].hamiltonian.rotated(by: alignment.rotation,budget: remainingBudget()))
            }
        }
        func lift(_ state: VivoCIState, partition: VivoActiveSpace, full h: VivoEmbeddedHamiltonian) throws -> VivoCIState {
            // Active columns are sorted. Inserting a complete doubly occupied
            // spatial core crosses an even number of fermion modes, so no
            // additional sign is introduced by this ordered embedding.
            var core: UInt64 = 0
            for p in partition.doublyOccupiedCore { core |= UInt64(3) << (2*p) }
            let determinants = state.determinants.map { det -> UInt64 in
                var full = core
                for (local,p) in partition.active.enumerated() { full |= ((det >> (2*local)) & 3) << (2*p) }
                return full
            }
            let result = VivoCIState(orbitalCount: h.orbitalCount,alphaElectrons: h.alphaElectrons,betaElectrons: h.betaElectrons,
                                    determinants: determinants,coefficients: state.coefficients)
            try result.validate(budget: request.budget)
            return result
        }
        func evaluate(_ partition: VivoActiveSpace, indices: [Int]) throws -> Evaluation {
            guard partition.active.count <= request.maximumActiveOrbitals else { throw Limit.activeSpace }
            let startOperators = operatorWork, startAuxiliary = auxiliaryWork
            var observations: [VivoSpacePointObservation] = [], states: [VivoCIState] = []
            for index in indices {
                guard pointEvaluations < request.maximumPointEvaluations else {
                    throw VivoChemistryError.resourceLimit("refinement point-evaluation budget")
                }
                pointEvaluations += 1
                let h = hamiltonians[index], n = h.orbitalCount
                try reserve(16*n*n*n*n)
                let reduced = try h.frozenCore(active: partition.active,doublyOccupiedCore: partition.doublyOccupiedCore,
                                               budget: remainingBudget())
                let budget = try remainingBudget()
                let state: VivoCIState, energy: Double, residual: Double, pt2: Double?
                var returned = false
                do {
                    switch request.solver {
                    case .directCI(let cfg):
                        let result = try VivoDirectCI.solve(reduced,configuration: cfg,budget: budget)
                        operatorWork += result.operatorApplications; returned = true
                        state = result.roots[0].state; energy = result.roots[0].energyHartree
                        residual = result.roots[0].eigenResidualNorm; pt2 = nil
                    case .selectedCI(let cfg):
                        let result = try VivoSelectedCI.solve(reduced,configuration: cfg,budget: budget)
                        operatorWork += result.operatorApplications; returned = true
                        guard result.converged else {
                            throw VivoChemistryError.convergence("selected CI at \(request.points[index].identifier): \(result.termination.rawValue); no property sensitivity inferred from an unconverged probe")
                        }
                        state = result.state; energy = result.variationalEnergyHartree
                        residual = result.diagnostics.fullResidualNorm; pt2 = result.pt2CorrectionHartree
                    }
                } catch {
                    // An exception does not report its partial inner work. Charge
                    // the entire remaining allowance instead of inventing a count.
                    if !returned { failedWorkReservation += budget.maximumOperatorApplications }
                    throw error
                }
                var information: VivoSelectiveOrbitalInformationResult?
                if request.collectOrbitalInformation {
                    let selection = try VivoOrbitalInformationSelection.all(orbitalCount: state.orbitalCount)
                    let report = try VivoSelectiveOrbitalInformation.analyze(state,selection: selection,budget: remainingBudget())
                    try reserve(report.reservedPrimitiveWork); information = report
                }
                states.append(try lift(state,partition: partition,full: h))
                observations.append(.init(pointIdentifier: request.points[index].identifier,energyHartree: energy,
                    determinantCount: state.determinants.count,ciResidualNorm: residual,
                    pt2CorrectionHartree: pt2,orbitalInformation: information))
            }
            var minimum: Double?, assessed = 0
            if indices.count > 1 {
                for j in 1..<indices.count where indices[j] == indices[j-1]+1 {
                    let a = states[j-1], b = states[j], electrons = max(1,a.alphaElectrons+a.betaElectrons)
                    let cost = try request.budget.elements([a.determinants.count,b.determinants.count,electrons,electrons,electrons])
                    let budget = try remainingBudget(); try reserve(cost)
                    let overlap = try VivoStateFollowing.overlaps(previous: [a],candidates: [b],orbitalOverlap: edges[indices[j]-1],budget: budget)[0,0]
                    let square = overlap*overlap
                    guard square.isFinite, square <= 1+1e-6, square >= request.minimumStateOverlapSquared else {
                        throw VivoChemistryError.convergence("refinement path state discontinuity at \(request.points[indices[j]].identifier)")
                    }
                    minimum = min(minimum ?? 1,square); assessed += 1
                }
            }
            // Nonadjacent physical overlaps are never approximated by multiplying
            // adjacent overlap matrices. Unchecked discovery edges stay explicit;
            // the final confirmation evaluates every adjacent physical edge.
            return .init(summary: .init(partition: partition,points: observations,
                minimumAdjacentStateOverlapSquared: minimum,assessedAdjacentEdges: assessed,
                hamiltonianOperatorApplications: operatorWork-startOperators,
                reservedAuxiliaryWork: auxiliaryWork-startAuxiliary),indices: indices,states: states)
        }
        func contrast(_ base: Evaluation, _ trial: Evaluation) throws -> VivoSpacePropertyShift {
            guard base.indices == trial.indices, Set(base.summary.partition.active).isSubset(of: Set(trial.summary.partition.active)) else {
                throw VivoChemistryError.invalid("paired refinement must use identical geometries and nested active spaces")
            }
            let a = base.summary.points.map(\.energyHartree), b = trial.summary.points.map(\.energyHartree)
            let r = base.indices.firstIndex(of: 0)!, p = base.indices.firstIndex(of: request.points.count-1)!
            let t = base.summary.points.firstIndex { $0.pointIdentifier == request.target.barrierPointIdentifier }!
            var minOverlap = 1.0
            for i in base.indices.indices {
                let left = base.states[i], right = trial.states[i]
                try reserve(left.determinants.count+right.determinants.count)
                // A coefficient overlap is exact HERE: both states are lifted
                // into the same complete orbital frame at the same geometry.
                let map = Dictionary(uniqueKeysWithValues: zip(left.determinants,left.coefficients))
                let overlap = zip(right.determinants,right.coefficients).reduce(0.0) { $0+(map[$1.0] ?? 0)*$1.1 }
                let square = overlap*overlap
                guard square.isFinite, square <= 1+1e-6, square >= request.minimumStateOverlapSquared else {
                    throw VivoChemistryError.convergence("active-space expansion lost the tracked electronic state")
                }
                guard b[i] <= a[i]+1e-8*max(1,abs(a[i])) else {
                    throw VivoChemistryError.convergence("nested electronic space increased variational energy; state/solver review required")
                }
                minOverlap = min(minOverlap,square)
            }
            let forward = (b[t]-b[r])-(a[t]-a[r]), reverse = (b[t]-b[p])-(a[t]-a[p])
            let reaction = (b[p]-b[r])-(a[p]-a[r])
            let relative = a.indices.map { abs((b[$0]-b[r])-(a[$0]-a[r])) }.max()!
            let absolute = zip(a,b).map { abs($0-$1) }.max()!
            let impact = max(max(abs(forward),abs(reverse),abs(reaction))/request.target.maximumBarrierShiftHartree,
                             relative/request.target.maximumRelativeProfileShiftHartree)
            guard [forward,reverse,reaction,relative,absolute,impact].allSatisfy(\.isFinite) else {
                throw VivoChemistryError.convergence("property sensitivity overflow")
            }
            return .init(forwardHartree: forward,reverseHartree: reverse,reactionHartree: reaction,
                maximumRelativeProfileHartree: relative,maximumAbsoluteEnergyHartree: absolute,
                minimumExpansionStateOverlapSquared: minOverlap,normalizedPropertyImpact: impact)
        }
        func confirmation(_ partition: VivoActiveSpace) throws -> VivoSpaceConfirmation {
            // Confirmation geometries have not supplied energies to selection.
            // If this fails, stop; do not tune on them and still call them held out.
            let indices = Array(request.points.indices)
            let chosen = try evaluate(partition,indices: indices)
            let unionSpace = VivoPropertyDirectedSpace.expanding(partition,with: pool)
            let union = unionSpace == partition ? chosen : try evaluate(unionSpace,indices: indices)
            let shift = try contrast(chosen,union)
            return .init(chosen: chosen.summary,candidateUnion: union.summary,shift: shift,
                         passed: shift.normalizedPropertyImpact <= 1 && chosen.summary.assessedAdjacentEdges == request.points.count-1 &&
                            union.summary.assessedAdjacentEdges == request.points.count-1)
        }
        func result(_ termination: VivoSpaceRefinementTermination, confirmation: VivoSpaceConfirmation? = nil,
                    failure: String? = nil) -> VivoPropertyDirectedSpaceResult {
            .init(schema: VivoPropertyDirectedSpaceResult.schema,request: request,transportRotations: rotations,
                transportMinimumSingularValues: minima,rounds: rounds,finalSpace: current,candidatePoolOrbitals: pool,
                unexaminedOrbitals: (0..<request.points[0].hamiltonian.orbitalCount).filter { !pool.contains($0) },
                confirmation: confirmation,termination: termination,failure: failure,pointEvaluations: pointEvaluations,
                hamiltonianOperatorApplications: operatorWork,reservedAuxiliaryWork: auxiliaryWork,
                reservedFailedSolverWork: failedWorkReservation,meaning: VivoPropertyDirectedSpace.meaning)
        }
        func execute() throws -> VivoPropertyDirectedSpaceResult {
            do { try prepare() } catch { return result(.incompleteExecution,failure: "orbital preparation: \(error)") }
            var baseline: Evaluation
            do { baseline = try evaluate(current,indices: discovery) }
            catch { return result(.incompleteExecution,failure: "initial electronic space: \(error)") }
            for round in 1...request.maximumRounds {
                var trials: [VivoSpaceRefinementTrial] = []
                var best: (String,Double,Evaluation)?, unionEvaluation: Evaluation?
                var activeBlock = "__initial__"
                do {
                    if current.active == pool {
                        rounds.append(.init(round: round,baseline: baseline.summary,trials: [],selectedBlockIdentifier: nil))
                        let checked = try confirmation(current)
                        return result(checked.passed ? .candidatePoolIncluded : .confirmationFailed,confirmation: checked)
                    }
                    var visited: Set<[Int]> = []
                    for block in request.candidateBlocks.sorted(by: { $0.identifier < $1.identifier }) {
                        let partition = VivoPropertyDirectedSpace.expanding(current,with: block.orbitals)
                        if partition == current || !visited.insert(partition.active).inserted { continue }
                        activeBlock = block.identifier
                        let before = chargedWork, candidate = try evaluate(partition,indices: discovery)
                        let shift = try contrast(baseline,candidate)
                        let priority = shift.normalizedPropertyImpact/Double(max(1,chargedWork-before))
                        trials.append(.init(blockIdentifier: block.identifier,evaluation: candidate.summary,shift: shift,priority: priority))
                        if partition.active == pool { unionEvaluation = candidate }
                        if shift.normalizedPropertyImpact > 1 && (best == nil || priority > best!.1) {
                            best = (block.identifier,priority,candidate)
                        }
                    }
                    if let best {
                        rounds.append(.init(round: round,baseline: baseline.summary,trials: trials,selectedBlockIdentifier: best.0))
                        current = best.2.summary.partition; baseline = best.2
                        continue
                    }
                    // Individually small effects cannot qualify the joint space.
                    // The complete PREDECLARED candidate union is an obligation.
                    activeBlock = "__candidate_union__"
                    let before = chargedWork
                    let union = try unionEvaluation ?? evaluate(VivoPropertyDirectedSpace.expanding(current,with: pool),indices: discovery)
                    let shift = try contrast(baseline,union)
                    if unionEvaluation == nil {
                        trials.append(.init(blockIdentifier: activeBlock,evaluation: union.summary,shift: shift,
                                            priority: shift.normalizedPropertyImpact/Double(max(1,chargedWork-before))))
                    }
                    if shift.normalizedPropertyImpact > 1 {
                        rounds.append(.init(round: round,baseline: baseline.summary,trials: trials,selectedBlockIdentifier: activeBlock))
                        current = union.summary.partition; baseline = union
                        continue
                    }
                    rounds.append(.init(round: round,baseline: baseline.summary,trials: trials,selectedBlockIdentifier: nil))
                    activeBlock = "__held_out_confirmation__"
                    let checked = try confirmation(current)
                    return result(checked.passed ? .sensitivityStable : .confirmationFailed,confirmation: checked)
                } catch {
                    // Keep completed trials even when a later mandatory trial
                    // fails. A budget limit is not evidence of low sensitivity.
                    if rounds.last?.round != round {
                        rounds.append(.init(round: round,baseline: baseline.summary,trials: trials,selectedBlockIdentifier: nil))
                    }
                    return result(error is Limit ? .activeSpaceLimit : .incompleteExecution,
                                  failure: "\(activeBlock): \(error)")
                }
            }
            return result(.roundLimit)
        }
    }
}
