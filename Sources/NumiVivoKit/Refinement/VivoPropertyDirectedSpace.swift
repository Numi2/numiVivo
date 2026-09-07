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
        let states: [[VivoCIState]] // [point][energy-ordered root], lifted into each full orbital frame
        let orbitalFrames: [VivoQMMatrix] // relative to the transported input frame
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
            let policy = request.solver.statePolicy
            var observations: [VivoSpacePointObservation] = [], states: [[VivoCIState]] = []
            var frames: [VivoQMMatrix] = []
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
                let rootStates: [VivoCIState], energies: [Double], residuals: [Double], pt2: Double?
                let frame: VivoQMMatrix
                var returned = false
                do {
                    switch request.solver {
                    case .directCI(let cfg), .multistateCI(let cfg, _):
                        let result = try VivoDirectCI.solve(reduced,configuration: cfg,budget: budget)
                        operatorWork += result.operatorApplications; returned = true
                        rootStates = result.roots.map(\.state); energies = result.roots.map(\.energyHartree)
                        residuals = result.roots.map(\.eigenResidualNorm); pt2 = nil
                        frame = try .identity(n)
                    case .selectedCI(let cfg):
                        let result = try VivoSelectedCI.solve(reduced,configuration: cfg,budget: budget)
                        operatorWork += result.operatorApplications; returned = true
                        guard result.converged else {
                            throw VivoChemistryError.convergence("selected CI at \(request.points[index].identifier): \(result.termination.rawValue); no sensitivity from an unconverged probe")
                        }
                        rootStates = [result.state]; energies = [result.variationalEnergyHartree]
                        residuals = [result.diagnostics.fullResidualNorm]; pt2 = result.pt2CorrectionHartree
                        frame = try .identity(n)
                    case .stateAveragedCASSCF(let cfg, _):
                        let result = try VivoMultiStateCASSCF.solve(h,partition: partition,configuration: cfg,budget: budget)
                        guard let work = result.numericalWork else {
                            throw VivoChemistryError.invalid("optimized refinement requires aggregate solver work accounting")
                        }
                        operatorWork += work.hamiltonianOperatorApplications
                        try reserve(work.reservedAuxiliaryWork); returned = true
                        guard result.converged else {
                            throw VivoChemistryError.convergence("SA-CASSCF refinement did not converge: \(result.termination.rawValue)")
                        }
                        rootStates = result.states.map(\.state); energies = result.states.map(\.energyHartree)
                        residuals = result.states.map(\.eigenResidualNorm); pt2 = nil
                        frame = result.orbitalRotation
                    }
                } catch {
                    // A failed inner call does not report all its partial work.
                    // Reserve the remaining allowance rather than inventing it.
                    if !returned { failedWorkReservation += max(0,request.budget.maximumOperatorApplications-chargedWork) }
                    throw error
                }
                try policy.validateSpectrum(energies)
                var rootReports: [VivoRefinementRootObservation] = []
                for root in rootStates.indices {
                    var information: VivoSelectiveOrbitalInformationResult?
                    if request.collectOrbitalInformation {
                        let selection = try VivoOrbitalInformationSelection.all(orbitalCount: rootStates[root].orbitalCount)
                        let report = try VivoSelectiveOrbitalInformation.analyze(rootStates[root],selection: selection,budget: remainingBudget())
                        try reserve(report.reservedPrimitiveWork); information = report
                    }
                    rootReports.append(.init(label: policy.labels[root],energyHartree: energies[root],
                        determinantCount: rootStates[root].determinants.count,residualNorm: residuals[root],orbitalInformation: information))
                }
                states.append(try rootStates.map { try lift($0,partition: partition,full: h) })
                frames.append(frame)
                observations.append(.init(pointIdentifier: request.points[index].identifier,energyHartree: energies[0],
                    determinantCount: rootStates[0].determinants.count,ciResidualNorm: residuals[0],
                    pt2CorrectionHartree: pt2,orbitalInformation: rootReports[0].orbitalInformation,
                    roots: policy.labels.count > 1 || request.solver.optimizesOrbitals ? rootReports : nil,
                    optimizedOrbitalRotation: request.solver.optimizesOrbitals ? frame : nil))
            }
            var minimum: Double?, assessed = 0
            if indices.count > 1 {
                for j in 1..<indices.count where indices[j] == indices[j-1]+1 {
                    let a = states[j-1], b = states[j], electrons = max(1,a[0].alphaElectrons+a[0].betaElectrons)
                    let cost = try request.budget.elements([a[0].determinants.count,b[0].determinants.count,electrons,electrons,electrons])
                    let budget = try remainingBudget(); try reserve(cost)
                    let physicalOverlap = try frames[j-1].transposed.multiplied(by: edges[indices[j]-1]).multiplied(by: frames[j])
                    let overlap = try VivoStateFollowing.overlaps(previous: a,candidates: b,orbitalOverlap: physicalOverlap,budget: budget)
                    let retained = try policy.minimumRetainedOverlapSquared(overlap,threshold: request.minimumStateOverlapSquared)
                    minimum = min(minimum ?? 1,retained); assessed += 1
                }
            }
            // Never invent nonadjacent overlaps by multiplying adjacent ones.
            // Confirmation evaluates every physical adjacent edge.
            return .init(summary: .init(partition: partition,points: observations,
                minimumAdjacentStateOverlapSquared: minimum,assessedAdjacentEdges: assessed,
                hamiltonianOperatorApplications: operatorWork-startOperators,
                reservedAuxiliaryWork: auxiliaryWork-startAuxiliary),indices: indices,states: states,orbitalFrames: frames)
        }
        func contrast(_ base: Evaluation, _ trial: Evaluation) throws -> VivoSpacePropertyShift {
            guard base.indices == trial.indices, Set(base.summary.partition.active).isSubset(of: Set(trial.summary.partition.active)) else {
                throw VivoChemistryError.invalid("paired refinement must use identical geometries and nested active spaces")
            }
            let policy = request.solver.statePolicy, roots = policy.labels.count
            let r = base.indices.firstIndex(of: 0)!, p = base.indices.firstIndex(of: request.points.count-1)!
            let t = base.summary.points.firstIndex { $0.pointIdentifier == request.target.barrierPointIdentifier }!
            func energies(_ evaluation: Evaluation, _ point: Int) -> [Double] {
                evaluation.summary.points[point].roots?.map(\.energyHartree) ?? [evaluation.summary.points[point].energyHartree]
            }
            var minOverlap = 1.0, maxGap = 0.0
            for i in base.indices.indices {
                let left = base.states[i], right = trial.states[i]
                let overlap: VivoQMMatrix
                if !request.solver.optimizesOrbitals {
                    // Exact HERE: both states are lifted into the same complete
                    // fixed orbital frame at the same geometry.
                    try reserve((left[0].determinants.count+right[0].determinants.count)*roots*roots)
                    var values = VivoQMMatrix(roots,roots)
                    for a in 0..<roots {
                        let map = Dictionary(uniqueKeysWithValues: zip(left[a].determinants,left[a].coefficients))
                        for b in 0..<roots {
                            values[a,b] = zip(right[b].determinants,right[b].coefficients).reduce(0.0) { $0+(map[$1.0] ?? 0)*$1.1 }
                        }
                    }
                    overlap = values
                } else {
                    let electrons = max(1,left[0].alphaElectrons+left[0].betaElectrons)
                    let cost = try request.budget.elements([left[0].determinants.count,right[0].determinants.count,electrons,electrons,electrons])
                    let budget = try remainingBudget(); try reserve(cost)
                    overlap = try VivoStateFollowing.overlaps(previous: left,candidates: right,
                        orbitalOverlap: base.orbitalFrames[i].transposed.multiplied(by: trial.orbitalFrames[i]),budget: budget)
                }
                minOverlap = min(minOverlap,try policy.minimumRetainedOverlapSquared(overlap,threshold: request.minimumStateOverlapSquared))
                let a = energies(base,i), b = energies(trial,i)
                for root in 0..<roots {
                    // Min-max ordering applies to fixed-orbital nested CI spaces,
                    // not to separately stationary state-averaged orbital optima.
                    if !request.solver.optimizesOrbitals && b[root] > a[root]+1e-8*max(1,abs(a[root])) {
                        throw VivoChemistryError.convergence("nested CI increased an energy-ordered variational root")
                    }
                    for other in 0..<root { maxGap = max(maxGap,abs((b[root]-b[other])-(a[root]-a[other]))) }
                }
            }
            var shifts: [VivoRefinementRootShift] = [], impact = 0.0
            for root in 0..<roots {
                let a = base.indices.indices.map { energies(base,$0)[root] }
                let b = base.indices.indices.map { energies(trial,$0)[root] }
                let forward = (b[t]-b[r])-(a[t]-a[r]), reverse = (b[t]-b[p])-(a[t]-a[p])
                let reaction = (b[p]-b[r])-(a[p]-a[r])
                let relative = a.indices.map { abs((b[$0]-b[r])-(a[$0]-a[r])) }.max()!
                let absolute = zip(a,b).map { abs($0-$1) }.max()!
                impact = max(impact,max(abs(forward),abs(reverse),abs(reaction))/request.target.maximumBarrierShiftHartree,
                             relative/request.target.maximumRelativeProfileShiftHartree)
                shifts.append(.init(label: policy.labels[root],forwardHartree: forward,reverseHartree: reverse,reactionHartree: reaction,
                    maximumRelativeProfileHartree: relative,maximumAbsoluteEnergyHartree: absolute))
            }
            if let tolerance = request.target.maximumStateGapShiftHartree { impact = max(impact,maxGap/tolerance) }
            guard impact.isFinite, maxGap.isFinite else { throw VivoChemistryError.convergence("property sensitivity overflow") }
            let first = shifts[0]
            return .init(forwardHartree: first.forwardHartree,reverseHartree: first.reverseHartree,reactionHartree: first.reactionHartree,
                maximumRelativeProfileHartree: first.maximumRelativeProfileHartree,maximumAbsoluteEnergyHartree: first.maximumAbsoluteEnergyHartree,
                minimumExpansionStateOverlapSquared: minOverlap,normalizedPropertyImpact: impact,
                rootShifts: roots > 1 ? shifts : nil,maximumStateGapShiftHartree: roots > 1 ? maxGap : nil)
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
