import Foundation

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
        var p = try prepare(request)
        let n = p.hamiltonians[0].orbitalCount
        let barrier = request.snapshots.firstIndex { $0.identifier == request.barrierPointIdentifier }!
        var reference = try p.hamiltonians.map { try VivoDirectCI.solve($0,
            configuration: .init(residualTolerance: request.referenceResidualTolerance), budget: request.budget).roots[0] }
        var occupations: [Double]?, ensembleRotation: VivoQMMatrix?
        if let policy = request.ensembleOrbitals {
            var density = VivoQMMatrix(n,n), work = 0
            for (point,ci) in reference.enumerated() {
                let index = Dictionary(uniqueKeysWithValues: ci.state.determinants.enumerated().map { ($0.element,$0.offset) })
                let needed = 2*n*n*ci.state.determinants.count
                guard needed <= request.budget.maximumOperatorApplications-work else {
                    throw VivoChemistryError.resourceLimit("ensemble one-particle-density work")
                }
                work += needed
                for a in 0..<n { for b in 0..<n { for spin in 0..<2 {
                    let value = VivoCIDensityMatrices.expectation(ci.state, index,
                        [.init(mode:2*b+spin,creation:false),.init(mode:2*a+spin,creation:true)])
                    density[a,b] += policy.pointWeights[point]*value
                } } }
            }
            let spectrum = try VivoQMDenseAlgebra.symmetricEigen(density, tolerance: 1e-13)
            guard spectrum.values.allSatisfy({ $0 >= -1e-8 && $0 <= 2+1e-8 }),
                  abs(spectrum.values.reduce(0,+)-Double(reference[0].state.alphaElectrons+reference[0].state.betaElectrons)) < 1e-8 else {
                throw VivoChemistryError.convergence("ensemble density particle trace or representability")
            }
            var rotation = VivoQMMatrix(n,n)
            for i in 0..<n { for j in 0..<n { rotation[i,j] = spectrum.vectors[i,n-1-j] } }
            occupations = Array(spectrum.values.reversed()); ensembleRotation = rotation
            p = .init(hamiltonians: try p.hamiltonians.map { try $0.rotated(by: rotation, budget: request.budget) },
                coefficients: try p.coefficients.map { try $0.multiplied(by: rotation) }, minima: p.minima,
                overlaps: try p.overlaps.map { try $0.congruence(rotation) })
            reference = try reference.map { ci in
                .init(method: ci.method, energyHartree: ci.energyHartree,
                    state: try VivoCIOrbitalFrame.rotated(ci.state, by: rotation, budget: request.budget),
                    eigenResidualNorm: ci.eigenResidualNorm, nextStateGapHartree: ci.nextStateGapHartree)
            }
        }
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
                    if let occupations, let policy = request.ensembleOrbitals {
                        var owner = [Int](repeating: 2, count: n)
                        for p in partition.active { owner[p] = 1 }
                        for p in partition.doublyOccupiedCore { owner[p] = 0 }
                        for i in 0..<n { for j in (i+1)..<n where owner[i] != owner[j] {
                            guard abs(occupations[i]-occupations[j]) >= policy.minimumOccupationBoundaryGap else {
                                throw VivoChemistryError.convergence("CAS boundary splits a near-degenerate ensemble occupation subspace")
                            }
                        } }
                    }
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
                    // Reserve both execution and path reconstruction, including
                    // a failing solve whose partial count cannot be inferred.
                    chargeForLevel = 2*cfg.maximumPointEvaluations
                    calls += chargeForLevel
                    guard calls <= request.maximumPointEvaluations else { throw VivoChemistryError.resourceLimit("campaign ECC reservation") }
                    let result = try VivoECCReactionPath.solve(points: eccPoints, configuration: cfg, budget: request.budget)
                    guard result.converged else { throw VivoChemistryError.convergence("campaign ECC path: \(result.termination)") }
                    try VivoECCReactionPath.validate(result, points: eccPoints, configuration: cfg, budget: request.budget)
                    values = result.pointResults.map(\.energyHartree)
                    size = result.pointResults.flatMap { $0.frame.clusters.map { $0.coefficients.columns } }.max()!
                    solverResidual = result.pointResults.flatMap { $0.frame.states.map { $0.biasedCI.eigenResidualNorm } }.max()!
                }
                guard solverResidual.isFinite, solverResidual >= 0,
                      solverResidual <= min(1e-8, min(request.acceptance.maximumBarrierErrorHartree,
                          request.acceptance.maximumRelativeProfileErrorHartree)/100) else {
                    throw VivoChemistryError.convergence("level solver accuracy is insufficient for the declared barrier target")
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
            ensembleOccupations: occupations, ensembleRotation: ensembleRotation,
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
        guard (result.ensembleOccupations == nil) == (rebuilt.ensembleOccupations == nil),
              (result.ensembleRotation == nil) == (rebuilt.ensembleRotation == nil),
              close(result.ensembleOccupations ?? [], rebuilt.ensembleOccupations ?? []) else {
            throw VivoChemistryError.invalid("campaign ensemble occupation or policy binding")
        }
        if let a = result.ensembleRotation, let b = rebuilt.ensembleRotation {
            guard a.rows == b.rows, a.columns == b.columns, close(a.values,b.values) else {
                throw VivoChemistryError.invalid("campaign ensemble orbital-frame reconstruction")
            }
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
