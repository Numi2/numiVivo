import Foundation

public struct VivoResidualBarrierRequest: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/shared-residual-barrier/v1"
    public let schema: String
    /// The original campaign remains executable and is retained in the output.
    /// Its full CI reference is used for assessment, not residual directions.
    public let baseline: VivoBarrierConvergenceRequest
    public let seed: VivoActiveSpace
    public let maximumRefinementRounds: Int
    public let space: VivoVariationalSpaceConfiguration
    public init(baseline: VivoBarrierConvergenceRequest, seed: VivoActiveSpace,
                maximumRefinementRounds: Int = 6, space: VivoVariationalSpaceConfiguration = .init()) {
        schema=Self.schema; self.baseline=baseline; self.seed=seed
        self.maximumRefinementRounds=maximumRefinementRounds; self.space=space
    }
    public func validate() throws {
        try baseline.validate(); try space.validate()
        guard schema==Self.schema, (1...32).contains(maximumRefinementRounds), seed.frozenOrbitals.isEmpty,
              space.projectedResidualToleranceHartree <= min(baseline.acceptance.maximumBarrierErrorHartree,
                  baseline.acceptance.maximumRelativeProfileErrorHartree)/100 else {
            throw VivoChemistryError.invalid("residual campaign rounds, seed or projected eigenvalue accuracy")
        }
    }
}
public struct VivoResidualBarrierLevel: Codable, Sendable, Equatable {
    public let round: Int
    public let variationalDimension: Int
    public let fullSectorDimension: Int
    public let newDirections: Int
    public let stateSpaceIsReduced: Bool
    public let points: [VivoVariationalCIResult]
    public let maximumBarrierErrorHartree: Double
    public let maximumProfileErrorHartree: Double
    public let maximumAbsoluteEnergyErrorHartree: Double
    public let maximumSuccessiveChangeHartree: Double?
    public let minimumPhysicalStateOverlapSquared: Double
    public let meetsReferenceAccuracy: Bool
}
public struct VivoResidualBarrierResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/shared-residual-barrier-result/v1"
    public let schema: String
    public let request: VivoResidualBarrierRequest
    public let baseline: VivoBarrierConvergenceResult
    public let levels: [VivoResidualBarrierLevel]
    public let finalSubspace: VivoQMMatrix
    public let hamiltonianOperatorApplications: Int
    public let reducedAccuracyEstablished: Bool
    public let acceptedRound: Int?
    public let meaning: String
}

/// Systematically enlarges ONE nested global variational space at all path
/// points. Shared physical residual vectors provide external correlation;
/// no PT2 denominator, fitted energy shift or full-reference eigenvector is
/// substituted for the requested approximation. Its cost remains finite-sector.
public enum VivoResidualBarrierCampaign {
    public static let meaning="reference-assessed fixed-geometry electronic profile; shared residual-enriched global CI, reduced eigenproblem dimension not fewer occupied orbital modes; complete determinant-sector vectors retained; no nuclear, Gibbs, kinetic or paper-reproduction certificate"
    public static func run(_ request: VivoResidualBarrierRequest) throws -> VivoResidualBarrierResult {
        try request.validate()
        let r=request.baseline
        let baseline=try VivoBarrierConvergence.run(r)
        let hs=try molecularHamiltonians(request:r,baseline:baseline), h=hs[0]
        let overlaps=try (0..<(r.snapshots.count-1)).map { i in
            let s=try VivoGaussianIntegralEngine.crossOverlap(leftSystem:r.snapshots[i].system,leftBasis:r.basis,
                rightSystem:r.snapshots[i+1].system,rightBasis:r.basis,budget:r.budget)
            return try baseline.orbitalCoefficients[i].transposed.multiplied(by:s).multiplied(by:baseline.orbitalCoefficients[i+1])
        }
        try request.seed.validate(for:h,budget:r.budget)
        var space=try VivoFockSpace(hamiltonian:h,configuration:request.space,budget:r.budget)
        try space.add(.init(identifier:"shared-CAS-seed",partition:request.seed),hamiltonian:h)
        let barrier=r.snapshots.firstIndex { $0.identifier==r.barrierPointIdentifier }!
        let truth=baseline.referenceEnergiesHartree, tForward=truth[barrier]-truth[0], tReverse=truth[barrier]-truth.last!
        var work=0, levels:[VivoResidualBarrierLevel]=[], newDirections=space.dimension
        for round in 0...request.maximumRefinementRounds {
            let points=try hs.map { try space.solve($0,work:&work) }
            let energies=points.map(\.energyHartree)
            guard zip(energies,truth).allSatisfy({ $0 >= $1-1e-8 }) else {
                throw VivoChemistryError.convergence("variational energy below the independent full-sector reference")
            }
            if let previous=levels.last {
                guard zip(energies,previous.points).allSatisfy({ $0 <= $1.energyHartree+1e-8 }) else {
                    throw VivoChemistryError.convergence("nested variational space increased an energy")
                }
            }
            let bError=max(abs(energies[barrier]-energies[0]-tForward),
                abs(energies[barrier]-energies.last!-tReverse),abs(energies.last!-energies[0]-truth.last!+truth[0]))
            let profile=energies.map { $0-energies[0] }, tProfile=truth.map { $0-truth[0] }
            let pError=zip(profile,tProfile).map { abs($0-$1) }.max()!
            let successive=levels.last.map { previous -> Double in
                let e=previous.points.map(\.energyHartree), p=e.map { $0-e[0] }
                return max(zip(profile,p).map { abs($0-$1) }.max()!,
                    abs(energies[barrier]-energies.last!-e[barrier]+e.last!))
            }
            var minimumOverlap=1.0
            for i in 0..<(points.count-1) {
                let overlap=try VivoStateFollowing.overlaps(previous:[points[i].state],candidates:[points[i+1].state],
                    orbitalOverlap:overlaps[i],budget:r.budget)[0,0]
                guard overlap.isFinite, overlap*overlap <= 1+1e-6 else { throw VivoChemistryError.convergence("residual path state overlap") }
                minimumOverlap=min(minimumOverlap,overlap*overlap)
            }
            guard minimumOverlap >= r.minimumReferenceOverlapSquared else {
                throw VivoChemistryError.convergence("residual model state discontinuity; refinement cannot silently switch roots")
            }
            levels.append(.init(round:round,variationalDimension:space.dimension,fullSectorDimension:space.determinants.count,
                newDirections:newDirections,stateSpaceIsReduced:space.dimension<space.determinants.count,points:points,
                maximumBarrierErrorHartree:bError,maximumProfileErrorHartree:pError,
                maximumAbsoluteEnergyErrorHartree:zip(energies,truth).map { abs($0-$1) }.max()!,
                maximumSuccessiveChangeHartree:successive,minimumPhysicalStateOverlapSquared:minimumOverlap,
                meetsReferenceAccuracy:bError<=r.acceptance.maximumBarrierErrorHartree && pError<=r.acceptance.maximumRelativeProfileErrorHartree))
            if round==request.maximumRefinementRounds || space.dimension==space.determinants.count { break }
            newDirections=try space.enrich(hamiltonians:hs,states:points,work:&work)
            if newDirections==0 { break } // no duplicated level to fabricate stability
        }
        let window=r.acceptance.minimumStableReducedLevels
        let reduced=levels.indices.filter { levels[$0].stateSpaceIsReduced }
        var accepted:Int?
        if reduced.count>=window, let last=reduced.last {
            let indices=Array(reduced.suffix(window))
            let stable=indices.allSatisfy { levels[$0].meetsReferenceAccuracy } && indices.dropFirst().allSatisfy {
                levels[$0].variationalDimension>levels[$0-1].variationalDimension &&
                (levels[$0].maximumSuccessiveChangeHartree ?? .infinity)<=r.acceptance.maximumSuccessiveChangeHartree
            }
            let later=levels[(last+1)...].allSatisfy { $0.meetsReferenceAccuracy && ($0.maximumSuccessiveChangeHartree ?? .infinity)<=r.acceptance.maximumSuccessiveChangeHartree }
            if stable && later { accepted=levels[last].round }
        }
        return .init(schema:VivoResidualBarrierResult.schema,request:request,baseline:baseline,levels:levels,
            finalSubspace:space.matrix,hamiltonianOperatorApplications:work,reducedAccuracyEstablished:accepted != nil,
            acceptedRound:accepted,meaning:meaning)
    }
    /// Reuse the already verified orbital frames; no energy-derived fitting.
    /// AO integral generation remains the existing authoritative implementation.
    static func molecularHamiltonians(request r:VivoBarrierConvergenceRequest,baseline:VivoBarrierConvergenceResult) throws -> [VivoEmbeddedHamiltonian] {
        guard baseline.request==r,baseline.orbitalCoefficients.count==r.snapshots.count else {
            throw VivoChemistryError.invalid("residual molecular preparation identity")
        }
        return try r.snapshots.indices.map { i in
            let ao=try VivoGaussianIntegralEngine.compute(system:r.snapshots[i].system,basis:r.basis,budget:r.budget)
            return try VivoEmbeddedHamiltonian.fromAO(ao,coefficients:baseline.orbitalCoefficients[i],
                alphaElectrons:r.snapshots[i].system.alphaElectrons,betaElectrons:r.snapshots[i].system.betaElectrons,
                orbitalIdentifiers:(0..<ao.count).map { "campaign-orbital-\($0)" },
                energyReference:"physical electronic Hamiltonian; geometry-specific AO scalar once; no thermal or standard-state correction",budget:r.budget)
        }
    }
    public static func validate(_ result:VivoResidualBarrierResult,request:VivoResidualBarrierRequest) throws {
        guard result.schema==VivoResidualBarrierResult.schema,result.request==request,result.meaning==meaning else {
            throw VivoChemistryError.invalid("residual campaign request or scientific meaning")
        }
        // Same bounded and deterministic operations reconstruct projector rank,
        // all physical CI vectors/residuals, energy differences and acceptance.
        // An error outside the subspace is never replaced by the projected one.
        let rebuilt=try run(request)
        guard result==rebuilt else { throw VivoChemistryError.invalid("residual campaign projector, energies or assessment differs on reconstruction") }
    }
    public static func template() -> VivoResidualBarrierRequest {
        .init(baseline:VivoBarrierBenchmarks.hydrogenExchange631G(ensemble:true),seed:.init(active:[0,1,2]),
              maximumRefinementRounds:6,space:.init(maximumDimension:90))
    }
}
