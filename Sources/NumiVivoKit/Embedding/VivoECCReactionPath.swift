import Foundation

/// An ordered electronic path. overlapWithPrevious is <previous orbital|current
/// orbital>, not an identity inferred from equal orbital counts or energies.
public struct VivoECCPathPoint: Codable, Sendable, Equatable {
    public let identifier: String
    public let hamiltonian: VivoEmbeddedHamiltonian
    public let overlapWithPrevious: VivoQMMatrix?
    public init(identifier: String, hamiltonian: VivoEmbeddedHamiltonian,
                overlapWithPrevious: VivoQMMatrix? = nil) {
        self.identifier=identifier; self.hamiltonian=hamiltonian; self.overlapWithPrevious=overlapWithPrevious
    }
}
public struct VivoECCPathConfiguration: Codable, Sendable, Equatable {
    public var embedding: VivoECCDMETConfiguration
    public var transportGroups: [[Int]]
    public var expectedBathOrbitals: [Int]
    public var pointWeights: [Double]
    public var minimumTransportSingularValue: Double
    public var minimumPathSubspaceSingularValue: Double
    public var minimumReferenceOverlapSquared: Double
    public var minimumTrialSubspaceSingularValue: Double
    /// Counts complete ECC point evaluations, not just outer path evaluations.
    /// Each point evaluation has the nested limits in embedding.matching.
    public var maximumPointEvaluations: Int
    public init(embedding: VivoECCDMETConfiguration, transportGroups: [[Int]],
                expectedBathOrbitals: [Int], pointWeights: [Double],
                minimumTransportSingularValue: Double = 0.5,
                minimumPathSubspaceSingularValue: Double = 0.5,
                minimumTrialSubspaceSingularValue: Double = 0.95,
                minimumReferenceOverlapSquared: Double = 0.1,
                maximumPointEvaluations: Int = 20000) {
        self.embedding=embedding; self.transportGroups=transportGroups
        self.expectedBathOrbitals=expectedBathOrbitals; self.pointWeights=pointWeights
        self.minimumTransportSingularValue=minimumTransportSingularValue
        self.minimumPathSubspaceSingularValue=minimumPathSubspaceSingularValue
        self.minimumTrialSubspaceSingularValue=minimumTrialSubspaceSingularValue
        self.minimumReferenceOverlapSquared=minimumReferenceOverlapSquared
        self.maximumPointEvaluations=maximumPointEvaluations
    }
    public func validate(points: [VivoECCPathPoint], budget: VivoChemistryBudget = .init()) throws {
        try budget.validate()
        guard (2...128).contains(points.count), let first=points.first,
              Set(points.map(\.identifier)).count==points.count,
              points.allSatisfy({ !$0.identifier.isEmpty && $0.identifier.utf8.count<=1024 }),
              points[0].overlapWithPrevious==nil, pointWeights.count==points.count,
              pointWeights.allSatisfy({ $0.isFinite && $0>0 }),
              abs(pointWeights.reduce(0,+)-1)<1e-12,
              (points.count...1_000_000).contains(maximumPointEvaluations),
              [minimumTransportSingularValue,minimumPathSubspaceSingularValue,minimumTrialSubspaceSingularValue,minimumReferenceOverlapSquared]
                .allSatisfy({ $0.isFinite && $0>0 && $0<=1 }) else {
            throw VivoChemistryError.invalid("ECC path identity, weights, transport or evaluation limit")
        }
        let n=first.hamiltonian.orbitalCount, groups=transportGroups.flatMap { $0 }
        guard groups.count==n, Set(groups)==Set(0..<n), transportGroups.allSatisfy({ !$0.isEmpty }),
              expectedBathOrbitals.count==embedding.fragments.count,
              zip(expectedBathOrbitals,embedding.fragments).allSatisfy({ $0>=0 && $0<=$1.maximumBathOrbitals }),
              embedding.bathSelection.maximumCandidates==nil else {
            throw VivoChemistryError.invalid("path requires fixed bath dimensions, complete transport groups and untruncated environment candidates")
        }
        _=try budget.elements([points.count,n,n,n,n],simultaneousArrays:64*embedding.fragments.count+16)
        for (index,point) in points.enumerated() {
            let h=point.hamiltonian
            try embedding.validate(hamiltonian:h,budget:budget)
            guard h.orbitalIdentifiers==first.hamiltonian.orbitalIdentifiers,
                  h.alphaElectrons==first.hamiltonian.alphaElectrons,
                  h.betaElectrons==first.hamiltonian.betaElectrons,
                  h.energyReference==first.hamiltonian.energyReference else {
                throw VivoChemistryError.invalid("ECC path changes orbital identity, spin sector or energy convention")
            }
            if index>0 {
                guard let s=point.overlapWithPrevious, s.rows==n, s.columns==n,
                      s.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("missing or malformed physical cross-geometry overlap") }
                // The overlap between two orthonormal finite subspaces is a
                // contraction. Equal dimensions alone do not make it unitary.
                let eigen=try VivoQMDenseAlgebra.symmetricEigen(s.transposed.multiplied(by:s))
                guard eigen.values.allSatisfy({ $0>=(-1e-10) && $0<=1+1e-6 }) else {
                    throw VivoChemistryError.invalid("path overlap is not a contraction between orthonormal orbital spaces")
                }
            }
        }
    }
}
public struct VivoECCPathContinuity: Codable, Sendable, Equatable {
    public let leftIdentifier: String
    public let rightIdentifier: String
    public let fragmentIdentifier: String
    public let referenceStateOverlapSquared: Double
    public let fragmentMinimumSingularValue: Double
    /// Nil when the declared bath is empty; no fictitious singular value.
    public let bathMinimumSingularValue: Double?
}
public struct VivoECCPathIteration: Codable, Sendable, Equatable {
    public let iteration: Int
    public let objective: Double
    public let gradientNorm: Double
    public let energyHartree: [Double]
    public let maximumEnergyChangeHartree: Double?
    public let acceptedRotationNorm: Double
}
public struct VivoECCPathResult: Codable, Sendable, Equatable {
    public static let schema="numivivo.org/ecc-path-result/v1"
    public let schema: String
    public let configuration: VivoECCPathConfiguration
    public let pointIdentifiers: [String]
    /// Raw point orbitals -> common path gauge, before the shared rotation.
    public let transportRotations: [VivoQMMatrix]
    public let transportMinimumSingularValues: [Double]
    public let sharedRotation: VivoQMMatrix
    /// Each ECC result is in its already transported/shared orbital frame and
    /// has no independent pointwise orbital optimization enabled.
    public let pointResults: [VivoECCDMETResult]
    public let continuity: [VivoECCPathContinuity]
    public let objective: Double
    public let sharedGradient: [Double]
    public let relativeElectronicEnergiesHartree: [Double]
    public let iterations: [VivoECCPathIteration]
    public let pointEvaluations: Int
    public let converged: Bool
    public let termination: String
}
public enum VivoECCPathError: Error, Sendable {
    case discontinuousSubspace(String)
}

/// Shared-orbital ECC-DMET. Each path trial runs the full single-fragment or
/// selected-moment partition cycle at every geometry. Only the orbital rotation
/// is shared; each geometry retains its own self-consistent correlation potential.
/// An electronic profile is not a TS characterization or a Gibbs barrier.
public enum VivoECCReactionPath {
    private struct Prepared {
        let h: [VivoEmbeddedHamiltonian]
        let transport: [VivoQMMatrix]
        let transportMinimum: [Double]
        let adjacent: [VivoQMMatrix] // edge i joins point i to i+1
    }
    private struct Evaluation {
        let points: [VivoECCDMETResult]
        let continuity: [VivoECCPathContinuity]
        let objective: Double
        var energies: [Double] { points.map(\.energyHartree) }
    }
    private static func fixed(_ cfg: VivoECCPathConfiguration) -> VivoECCDMETConfiguration {
        var result=cfg.embedding
        result.localityGroups=[]
        return result
    }
    private static func pairs(_ cfg: VivoECCPathConfiguration) -> [(Int,Int)] {
        cfg.embedding.localityGroups.flatMap { g in g.indices.flatMap { i in
            ((i+1)..<g.count).map { (g[i],g[$0]) }
        } }
    }
    private static func prepare(_ points: [VivoECCPathPoint], _ cfg: VivoECCPathConfiguration,
                                _ budget: VivoChemistryBudget) throws -> Prepared {
        try cfg.validate(points:points,budget:budget)
        let n=points[0].hamiltonian.orbitalCount, identity=try VivoQMMatrix.identity(n)
        var rotations=[identity], minima=[1.0], h=[points[0].hamiltonian], edges:[VivoQMMatrix]=[]
        for i in 1..<points.count {
            let s=points[i].overlapWithPrevious!
            let alignment=try VivoOrbitalPathTransport.align(referenceCoefficients:rotations[i-1],
                currentCoefficients:identity,crossAOOverlap:s,groups:cfg.transportGroups,
                minimumSingularValue:cfg.minimumTransportSingularValue)
            rotations.append(alignment.rotation); minima.append(alignment.minimumSubspaceSingularValue)
            edges.append(try rotations[i-1].transposed.multiplied(by:s).multiplied(by:alignment.rotation))
            h.append(try points[i].hamiltonian.rotated(by:alignment.rotation,budget:budget))
        }
        return .init(h:h,transport:rotations,transportMinimum:minima,adjacent:edges)
    }
    private static func rotation(_ values: [Double], _ pairs: [(Int,Int)], _ n: Int) throws -> VivoQMMatrix {
        var k=VivoQMMatrix(n,n)
        for (i,pair) in pairs.enumerated() { k[pair.0,pair.1]=values[i]; k[pair.1,pair.0] = -values[i] }
        return try VivoQMDenseAlgebra.orbitalRotation(generator:k)
    }
    private static func checkRotation(_ u: VivoQMMatrix, _ cfg: VivoECCPathConfiguration, _ n: Int) throws {
        guard u.rows==n,u.columns==n,u.values.allSatisfy(\.isFinite),
              try u.transposed.multiplied(by:u).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {
            throw VivoChemistryError.invalid("shared path rotation is not orthogonal")
        }
        var owner=Array(0..<n)
        if !cfg.embedding.localityGroups.isEmpty {
            for (i,group) in cfg.embedding.localityGroups.enumerated() { for p in group { owner[p]=i } }
        }
        for p in 0..<n { for q in 0..<n where owner[p] != owner[q] && abs(u[p,q])>1e-10 {
            throw VivoChemistryError.invalid("shared path rotation violates declared locality")
        } }
        // Singleton fixed orbitals may not secretly undergo sign/reflection restarts.
        if cfg.embedding.localityGroups.isEmpty,
           try u.adding(.identity(n),scale:-1).frobeniusNorm>1e-10 {
            throw VivoChemistryError.invalid("fixed-orbital path requires identity shared rotation")
        }
    }
    private static func columns(_ c: VivoQMMatrix, _ range: Range<Int>) -> VivoQMMatrix {
        var out=VivoQMMatrix(c.rows,range.count)
        for p in 0..<c.rows { for (q,k) in range.enumerated() { out[p,q]=c[p,k] } }
        return out
    }
    private static func minimum(_ a: VivoQMMatrix, _ b: VivoQMMatrix, _ metric: VivoQMMatrix) throws -> Double {
        let overlap=try a.transposed.multiplied(by:metric).multiplied(by:b)
        let eigen=try VivoQMDenseAlgebra.symmetricEigen(overlap.transposed.multiplied(by:overlap),tolerance:1e-13)
        guard eigen.values.allSatisfy({ $0>=(-1e-10) && $0<=1+1e-6 }) else {
            throw VivoChemistryError.invalid("invalid path projector overlap spectrum")
        }
        return sqrt(max(0,eigen.values[0]))
    }
    private static func compare(_ a: VivoECCDMETResult, _ b: VivoECCDMETResult,
                                metric: VivoQMMatrix, threshold: Double,
                                left: String, right: String, stateThreshold: Double,
                                budget: VivoChemistryBudget) throws -> [VivoECCPathContinuity] {
        guard a.frame.clusters.count==b.frame.clusters.count else { throw VivoECCPathError.discontinuousSubspace("fragment count") }
        let stateOverlap=try VivoStateFollowing.overlaps(previous:[a.frame.reference.state],candidates:[b.frame.reference.state],
            orbitalOverlap:metric,budget:budget)[0,0]
        let stateSquare=stateOverlap*stateOverlap
        guard stateSquare.isFinite,stateSquare<=1+1e-6,stateSquare>=stateThreshold else {
            throw VivoECCPathError.discontinuousSubspace("reference state changes discontinuously; refine the path or use an explicit multistate model")
        }
        var checks:[VivoECCPathContinuity]=[]
        for (ac,bc) in zip(a.frame.clusters,b.frame.clusters) {
            let f=ac.fragment.orbitals.count, k=ac.coefficients.columns
            guard ac.fragment==bc.fragment,k==bc.coefficients.columns else { throw VivoECCPathError.discontinuousSubspace("cluster identity/dimension") }
            let fragment=try minimum(columns(ac.coefficients,0..<f),columns(bc.coefficients,0..<f),metric)
            let bath=try k>f ? minimum(columns(ac.coefficients,f..<k),columns(bc.coefficients,f..<k),metric) : nil
            guard fragment>=threshold, bath.map({$0>=threshold}) ?? true else {
                throw VivoECCPathError.discontinuousSubspace("\(left) -> \(right), fragment \(ac.fragment.identifier): fragment or bath subspace lost")
            }
            checks.append(.init(leftIdentifier:left,rightIdentifier:right,fragmentIdentifier:ac.fragment.identifier,
                referenceStateOverlapSquared:stateSquare,fragmentMinimumSingularValue:fragment,bathMinimumSingularValue:bath))
        }
        return checks
    }
    private static func evaluate(_ prepared: Prepared, _ points: [VivoECCPathPoint], _ cfg: VivoECCPathConfiguration,
                                 _ u: VivoQMMatrix, _ budget: VivoChemistryBudget, _ calls: inout Int) throws -> Evaluation {
        guard calls<=cfg.maximumPointEvaluations-points.count else {
            throw VivoChemistryError.resourceLimit("shared ECC path point-evaluation budget")
        }
        var results:[VivoECCDMETResult]=[]
        for h in prepared.h {
            calls+=1
            // Deterministic zero potential initialization at each geometry/trial:
            // no path-history-dependent choice of an inner matching root.
            let point=try VivoECCDMET.solve(h.rotated(by:u,budget:budget),configuration:fixed(cfg),budget:budget)
            guard point.converged else { throw VivoChemistryError.convergence("path inner ECC: \(point.termination)") }
            let dimensions=point.frame.clusters.map { $0.coefficients.columns-$0.fragment.orbitals.count }
            guard dimensions==cfg.expectedBathOrbitals else { throw VivoECCPathError.discontinuousSubspace("declared bath dimension differs from constructed cluster") }
            results.append(point)
        }
        var continuity:[VivoECCPathContinuity]=[]
        for i in prepared.adjacent.indices {
            continuity += try compare(results[i],results[i+1],metric:prepared.adjacent[i].congruence(u),
                threshold:cfg.minimumPathSubspaceSingularValue,left:points[i].identifier,right:points[i+1].identifier,
                stateThreshold:cfg.minimumReferenceOverlapSquared,budget:budget)
        }
        let objective=zip(results,cfg.pointWeights).reduce(0.0) { sum,pair in
            sum+pair.1*pair.0.frame.clusters.reduce(0.0) { $0+$1.qio.objective }/Double(pair.0.frame.clusters.count)
        }
        guard objective.isFinite else { throw VivoChemistryError.convergence("nonfinite shared path objective") }
        return .init(points:results,continuity:continuity,objective:objective)
    }
    private static func trialContinuity(_ a: Evaluation, _ b: Evaluation, _ old: VivoQMMatrix,
                                        _ new: VivoQMMatrix, _ cfg: VivoECCPathConfiguration, _ budget: VivoChemistryBudget) throws {
        let metric=try old.transposed.multiplied(by:new)
        for i in a.points.indices {
            _=try compare(a.points[i],b.points[i],metric:metric,threshold:cfg.minimumTrialSubspaceSingularValue,
                left:"accepted-frame",right:"trial-frame",stateThreshold:cfg.minimumTrialSubspaceSingularValue*cfg.minimumTrialSubspaceSingularValue,budget:budget)
        }
    }
    private static func gradient(_ prepared: Prepared, _ points: [VivoECCPathPoint], _ cfg: VivoECCPathConfiguration,
                                 _ u: VivoQMMatrix, _ base: Evaluation, _ budget: VivoChemistryBudget,
                                 _ calls: inout Int) throws -> [Double] {
        let indices=pairs(cfg), n=u.rows, step=cfg.embedding.orbitalDifferenceStep
        var result:[Double]=[]
        for i in indices.indices {
            var delta=[Double](repeating:0,count:indices.count); delta[i]=step
            let plus=try u.multiplied(by:rotation(delta,indices,n)); delta[i] = -step
            let minus=try u.multiplied(by:rotation(delta,indices,n))
            let a=try evaluate(prepared,points,cfg,plus,budget,&calls)
            let b=try evaluate(prepared,points,cfg,minus,budget,&calls)
            try trialContinuity(base,a,u,plus,cfg,budget); try trialContinuity(base,b,u,minus,cfg,budget)
            result.append((a.objective-b.objective)/(2*step))
        }
        guard result.allSatisfy(\.isFinite) else { throw VivoChemistryError.convergence("nonfinite shared ECC derivative") }
        return result
    }
    public static func solve(points: [VivoECCPathPoint], configuration cfg: VivoECCPathConfiguration,
                             initialSharedRotation: VivoQMMatrix? = nil,
                             budget: VivoChemistryBudget = .init()) throws -> VivoECCPathResult {
        let prepared=try prepare(points,cfg,budget), n=prepared.h[0].orbitalCount, indices=pairs(cfg)
        var u=try initialSharedRotation ?? VivoQMMatrix.identity(n)
        try checkRotation(u,cfg,n)
        var calls=0, current=try evaluate(prepared,points,cfg,u,budget,&calls)
        var history:[VivoECCPathIteration]=[], g:[Double]=[], previous:[Double]?, stepNorm=0.0
        var converged=false, termination="iteration-limit"
        for iteration in 0...cfg.embedding.maximumOrbitalIterations {
            let required=2*indices.count*points.count
            guard required<=cfg.maximumPointEvaluations-calls else { termination="point-evaluation-limit"; break }
            g=try gradient(prepared,points,cfg,u,current,budget,&calls)
            let norm=g.reduce(0.0) { hypot($0,$1) }
            let energyChange=previous.map { zip($0,current.energies).map { abs($0-$1) }.max()! }
            history.append(.init(iteration:iteration,objective:current.objective,gradientNorm:norm,
                energyHartree:current.energies,maximumEnergyChangeHartree:energyChange,acceptedRotationNorm:stepNorm))
            if norm<=cfg.embedding.orbitalGradientTolerance {
                if energyChange == nil || energyChange!<=cfg.embedding.energyToleranceHartree {
                    converged=true; termination=indices.isEmpty ? "fixed-shared-frame" : "shared-orbital-stationarity"; break
                }
                previous=current.energies; continue
            }
            if iteration==cfg.embedding.maximumOrbitalIterations { break }
            let bound=min(1,cfg.embedding.maximumOrbitalStep/max(norm,1e-300)), direction=g.map { -bound*$0 }
            var scale=1.0, accepted=false
            for _ in 0..<20 {
                if points.count>cfg.maximumPointEvaluations-calls { termination="point-evaluation-limit"; break }
                let next=try u.multiplied(by:rotation(direction.map { scale*$0 },indices,n))
                do {
                    let candidate=try evaluate(prepared,points,cfg,next,budget,&calls)
                    try trialContinuity(current,candidate,u,next,cfg,budget)
                    if candidate.objective<=current.objective-1e-4*scale*bound*norm*norm {
                        previous=current.energies; current=candidate; u=next; stepNorm=scale*bound*norm; accepted=true; break
                    }
                } catch is VivoECCPathError { /* reject this discontinuous step, never accept changed spaces */ }
                scale*=0.5
            }
            if !accepted { if termination != "point-evaluation-limit" { termination="shared-line-search-stalled" }; break }
        }
        return .init(schema:VivoECCPathResult.schema,configuration:cfg,pointIdentifiers:points.map(\.identifier),
            transportRotations:prepared.transport,transportMinimumSingularValues:prepared.transportMinimum,sharedRotation:u,
            pointResults:current.points,continuity:current.continuity,objective:current.objective,sharedGradient:g,
            relativeElectronicEnergiesHartree:current.energies.map { $0-current.energies[0] },iterations:history,
            pointEvaluations:calls,converged:converged,termination:termination)
    }
    /// Reconstruct physical transports and every point's final ECC state, then
    /// independently recompute shared stationarity. Pointwise stationary flags
    /// do not prove the shared orbital frame is stationary.
    public static func validate(_ result: VivoECCPathResult, points: [VivoECCPathPoint],
                                configuration cfg: VivoECCPathConfiguration, budget: VivoChemistryBudget = .init()) throws {
        let prepared=try prepare(points,cfg,budget), n=prepared.h[0].orbitalCount
        try checkRotation(result.sharedRotation,cfg,n)
        guard result.schema==VivoECCPathResult.schema,result.configuration==cfg,result.converged,
              result.pointIdentifiers==points.map(\.identifier),result.pointResults.count==points.count,
              result.transportRotations.count==points.count,result.transportMinimumSingularValues.count==points.count,
              result.relativeElectronicEnergiesHartree.count==points.count,
              result.sharedGradient.count==pairs(cfg).count,result.sharedGradient.allSatisfy(\.isFinite),
              result.sharedGradient.reduce(0.0,{hypot($0,$1)})<=cfg.embedding.orbitalGradientTolerance,
              result.pointEvaluations>0,result.pointEvaluations<=cfg.maximumPointEvaluations,
              !result.iterations.isEmpty,result.iterations.count<=cfg.embedding.maximumOrbitalIterations+1 else {
            throw VivoChemistryError.invalid("ECC path result method, dimensions or convergence contract")
        }
        for i in points.indices {
            guard try result.transportRotations[i].adding(prepared.transport[i],scale:-1).frobeniusNorm<1e-8,
                  abs(result.transportMinimumSingularValues[i]-prepared.transportMinimum[i])<1e-8 else {
                throw VivoChemistryError.invalid("path transport differs from physical cross overlaps")
            }
            try VivoECCDMET.validate(result.pointResults[i],hamiltonian:prepared.h[i].rotated(by:result.sharedRotation,budget:budget),
                configuration:fixed(cfg),budget:budget)
        }
        var calls=0
        let rebuilt=try evaluate(prepared,points,cfg,result.sharedRotation,budget,&calls)
        let actual=try gradient(prepared,points,cfg,result.sharedRotation,rebuilt,budget,&calls)
        let tolerance=cfg.embedding.orbitalGradientTolerance
        guard actual.reduce(0.0,{hypot($0,$1)})<=1.01*tolerance,
              zip(actual,result.sharedGradient).reduce(0.0,{hypot($0,$1.0-$1.1)})<=max(1e-8,0.1*tolerance),
              result.objective.isFinite,abs(rebuilt.objective-result.objective)<1e-7,
              result.continuity.count==rebuilt.continuity.count else { throw VivoChemistryError.invalid("path shared stationarity or continuity claim") }
        for (a,b) in zip(result.continuity,rebuilt.continuity) {
            guard a.leftIdentifier==b.leftIdentifier,a.rightIdentifier==b.rightIdentifier,a.fragmentIdentifier==b.fragmentIdentifier,
                  abs(a.fragmentMinimumSingularValue-b.fragmentMinimumSingularValue)<1e-8,
                  abs(a.referenceStateOverlapSquared-b.referenceStateOverlapSquared)<1e-8,
                  (a.bathMinimumSingularValue==nil)==(b.bathMinimumSingularValue==nil),
                  abs((a.bathMinimumSingularValue ?? 0)-(b.bathMinimumSingularValue ?? 0))<1e-8 else {
                throw VivoChemistryError.invalid("path fragment/bath continuity diagnostic differs")
            }
        }
        for i in points.indices {
            let expected=rebuilt.energies[i]-rebuilt.energies[0]
            guard result.relativeElectronicEnergiesHartree[i].isFinite,
                  abs(expected-result.relativeElectronicEnergiesHartree[i])<=max(1e-8,10*cfg.embedding.energyToleranceHartree) else {
                throw VivoChemistryError.invalid("path relative energy/reference mismatch")
            }
        }
    }
}
