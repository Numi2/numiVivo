import Foundation

public struct VivoRootFollowingResult:Codable,Sendable,Equatable {
    public let states:[VivoCIResult]
    public let candidateIndices:[Int]
    public let overlapsSquared:[Double]
}
public enum VivoStateFollowing {
    /// Exact Slater overlaps including inactive doubly occupied core. Coefficient
    /// dot products alone are invalid after changing the orbital frame. The
    /// supplied matrix is <previous spatial orbitals | candidate orbitals>.
    public static func overlaps(previous:[VivoCIState],candidates:[VivoCIState],orbitalOverlap:VivoQMMatrix,
                                coreCount:Int=0,budget:VivoChemistryBudget = .init()) throws -> VivoQMMatrix {
        guard let left=previous.first,let right=candidates.first,coreCount>=0 else { throw VivoChemistryError.invalid("empty root tracking states") }
        for state in previous+candidates { try state.validate(budget:budget) }
        let n=left.orbitalCount+coreCount
        guard right.orbitalCount==left.orbitalCount,orbitalOverlap.rows==n,orbitalOverlap.columns==n,
              orbitalOverlap.values.allSatisfy(\.isFinite),n<=31,
              previous.allSatisfy({$0.determinants==left.determinants}),candidates.allSatisfy({$0.determinants==right.determinants}),
              left.alphaElectrons==right.alphaElectrons,left.betaElectrons==right.betaElectrons else { throw VivoChemistryError.invalid("root overlap sectors or frames") }
        let electrons=left.alphaElectrons+left.betaElectrons+2*coreCount
        _ = try budget.elements([left.determinants.count,right.determinants.count],simultaneousArrays:3)
        let work=try budget.elements([left.determinants.count,right.determinants.count,max(1,electrons),max(1,electrons),max(1,electrons)])
        guard work<=budget.maximumOperatorApplications else { throw VivoChemistryError.resourceLimit("determinant-overlap work") }
        func modes(_ d:UInt64)->[Int] {
            let core=coreCount==0 ? UInt64(0):(UInt64(1)<<(2*coreCount))-1
            let bits=(d<<(2*coreCount))|core
            return (0..<(2*n)).filter { bits & (UInt64(1)<<$0) != 0 }
        }
        func determinant(_ a:VivoQMMatrix)->Double {
            if a.rows==0 { return 1 };var a=a,value=1.0
            for k in 0..<a.rows {
                var pivot=k
                for i in (k+1)..<a.rows where abs(a[i,k])>abs(a[pivot,k]) { pivot=i }
                if abs(a[pivot,k])<1e-15 { return 0 }
                if pivot != k { for j in k..<a.rows { let temp=a[k,j];a[k,j]=a[pivot,j];a[pivot,j]=temp };value = -value }
                value*=a[k,k]
                for i in (k+1)..<a.rows { let factor=a[i,k]/a[k,k];for j in (k+1)..<a.rows { a[i,j]-=factor*a[k,j] } }
            };return value
        }
        let lm=left.determinants.map(modes),rm=right.determinants.map(modes)
        var determinants=VivoQMMatrix(lm.count,rm.count)
        for i in lm.indices { for j in rm.indices {
            var a=VivoQMMatrix(electrons,electrons)
            for p in 0..<electrons { for q in 0..<electrons {
                let x=lm[i][p],y=rm[j][q];a[p,q]=x%2==y%2 ? orbitalOverlap[x/2,y/2]:0
            } };determinants[i,j]=determinant(a)
        } }
        var l=VivoQMMatrix(previous.count,lm.count),r=VivoQMMatrix(rm.count,candidates.count)
        for i in previous.indices { for j in lm.indices { l[i,j]=previous[i].coefficients[j] } }
        for i in candidates.indices { for j in rm.indices { r[j,i]=candidates[i].coefficients[j] } }
        return try l.multiplied(by:determinants).multiplied(by:r)
    }
    /// One-to-one maximum-overlap assignment by bounded subset dynamic program.
    /// Root loss rejects instead of silently replacing an excited state.
    public static func match(previous:[VivoCIResult],candidates:[VivoCIResult],orbitalOverlap:VivoQMMatrix,
                             coreCount:Int=0,minimumOverlapSquared:Double=0.25,budget:VivoChemistryBudget = .init()) throws -> VivoRootFollowingResult {
        guard previous.count==candidates.count,(1...12).contains(previous.count),minimumOverlapSquared.isFinite,
              minimumOverlapSquared>0,minimumOverlapSquared<=1 else { throw VivoChemistryError.invalid("root assignment contract") }
        let overlap=try overlaps(previous:previous.map(\.state),candidates:candidates.map(\.state),orbitalOverlap:orbitalOverlap,coreCount:coreCount,budget:budget)
        let n=previous.count,total=1<<n
        var best=[Double](repeating:-Double.infinity,count:total),parent=[Int](repeating:-1,count:total),chosen=parent;best[0]=0
        for mask in 0..<total {
            let i=mask.nonzeroBitCount;if i>=n { continue }
            for j in 0..<n where mask & (1<<j)==0 {
                let next=mask|(1<<j),score=best[mask]+overlap[i,j]*overlap[i,j]
                if score>best[next] { best[next]=score;parent[next]=mask;chosen[next]=j }
            }
        }
        var assignment=[Int](repeating:0,count:n),mask=total-1
        for i in (0..<n).reversed() { assignment[i]=chosen[mask];mask=parent[mask] }
        var states:[VivoCIResult]=[],squares:[Double]=[]
        for i in 0..<n {
            let j=assignment[i],value=overlap[i,j],square=value*value
            guard square>=minimumOverlapSquared else { throw VivoChemistryError.convergence("tracked root left the retained state/orbital subspace") }
            let c=candidates[j],s=c.state,sign=value<0 ? -1.0:1.0
            let state=VivoCIState(orbitalCount:s.orbitalCount,alphaElectrons:s.alphaElectrons,betaElectrons:s.betaElectrons,
                determinants:s.determinants,coefficients:s.coefficients.map { sign*$0 })
            states.append(.init(method:c.method,energyHartree:c.energyHartree,state:state,eigenResidualNorm:c.eigenResidualNorm,nextStateGapHartree:c.nextStateGapHartree));squares.append(square)
        };return .init(states:states,candidateIndices:assignment,overlapsSquared:squares)
    }
}
public struct VivoMultiStateCASSCFConfiguration:Codable,Sendable,Equatable {
    public var weights:[Double]
    public var rootLabels:[String]
    public var optimization:VivoCASSCFConfiguration
    public var davidson:VivoDavidsonConfiguration
    public var followRoots:Bool
    public var minimumOverlapSquared:Double
    public init(weights:[Double],rootLabels:[String],optimization:VivoCASSCFConfiguration = .init(),followRoots:Bool=true,minimumOverlapSquared:Double=0.25) {
        self.weights=weights;self.rootLabels=rootLabels;self.optimization=optimization;self.followRoots=followRoots;self.minimumOverlapSquared=minimumOverlapSquared
        davidson = .init(roots:weights.count,maximumSubspace:max(40,4*weights.count+2))
    }
}
public struct VivoMultiStateCASSCFResult:Codable,Sendable,Equatable {
    public let configuration:VivoMultiStateCASSCFConfiguration
    public let orbitalRotation:VivoQMMatrix
    public let states:[VivoCIResult]
    public let weightedEnergyHartree:Double
    public let orbitalGradient:[Double]
    public let termination:VivoCASSCFTermination
    public var converged:Bool { termination == .converged }
    public let iterations:[VivoCASSCFIteration]
    public let finalRootAssignment:[Int]
    public let finalOverlapSquared:[Double]
}
public enum VivoMultiStateCASSCF {
    /// State-averaged energy with analytic weighted CI-RDM orbital gradients.
    /// A weight vector with one nonzero component implements state-specific
    /// root-followed optimization. Individual total-spin purity is not imposed.
    public static func solve(_ input:VivoEmbeddedHamiltonian,partition:VivoActiveSpace,
                             configuration cfg:VivoMultiStateCASSCFConfiguration,initialRotation:VivoQMMatrix?=nil,
                             budget:VivoChemistryBudget = .init()) throws -> VivoMultiStateCASSCFResult {
        try partition.validate(for:input,budget:budget);try cfg.optimization.validate()
        guard (1...12).contains(cfg.weights.count),cfg.weights.count==cfg.rootLabels.count,
              Set(cfg.rootLabels).count==cfg.rootLabels.count,cfg.rootLabels.allSatisfy({!$0.isEmpty}),
              cfg.weights.allSatisfy({$0.isFinite && $0>=0}),abs(cfg.weights.reduce(0,+)-1)<1e-12,
              cfg.davidson.roots==cfg.weights.count,cfg.minimumOverlapSquared.isFinite,cfg.minimumOverlapSquared>0,cfg.minimumOverlapSquared<=1 else {
            throw VivoChemistryError.invalid("state-average weights, labels or root policy")
        }
        let n=input.orbitalCount,frozen=Set(partition.frozenOrbitals),selected=partition.doublyOccupiedCore+partition.active,opt=cfg.optimization
        var group=[Int](repeating:2,count:n)
        for p in partition.doublyOccupiedCore { group[p]=0 };for p in partition.active { group[p]=1 }
        var pairs:[VivoOrbitalRotationPair]=[]
        for p in 0..<n { for q in (p+1)..<n where group[p] != group[q] && !frozen.contains(p) && !frozen.contains(q) { pairs.append(.init(first:p,second:q)) } }
        var u=try initialRotation ?? .identity(n),evaluations=0
        guard u.rows==n,u.columns==n,u.values.allSatisfy(\.isFinite) else { throw VivoChemistryError.invalid("multistate initial rotation dimensions") }
        for p in frozen { for q in 0..<n where abs(u[q,p]-(p==q ? 1:0))>1e-12 { throw VivoChemistryError.invalid("initial rotation changes frozen orbital") } }
        func frame(_ rotation:VivoQMMatrix)->VivoQMMatrix {
            var c=VivoQMMatrix(n,selected.count);for i in 0..<n { for j in selected.indices { c[i,j]=rotation[i,selected[j]] } };return c
        }
        func spectrum(_ rotation:VivoQMMatrix) throws -> (VivoEmbeddedHamiltonian,[VivoCIResult]) {
            guard evaluations<opt.maximumEnergyEvaluations else { throw VivoChemistryError.resourceLimit("multistate energy evaluations") };evaluations+=1
            let full=try input.rotated(by:rotation,budget:budget),active=try full.frozenCore(active:partition.active,doublyOccupiedCore:partition.doublyOccupiedCore,budget:budget)
            return (full,try VivoDirectCI.solve(active,configuration:cfg.davidson,budget:budget).roots)
        }
        func gradient(_ full:VivoEmbeddedHamiltonian,_ states:[VivoCIResult]) throws -> [Double] {
            var g=[Double](repeating:0,count:pairs.count)
            for i in states.indices where cfg.weights[i]>0 {
                let part=try VivoCASOrbitalGradient.compute(hamiltonian:full,partition:partition,state:states[i].state,pairs:pairs,budget:budget)
                for j in g.indices { g[j]+=cfg.weights[i]*part[j] }
            };return g
        }
        func energy(_ states:[VivoCIResult])->Double { zip(cfg.weights,states).reduce(0.0) { $0+$1.0*$1.1.energyHartree } }
        func norm(_ x:[Double])->Double { x.reduce(0.0) { hypot($0,$1) } }
        var (full,states)=try spectrum(u),value=energy(states),history:[VivoCASSCFIteration]=[],g:[Double]=[]
        var previous:Double?,stepNorm=0.0,termination=VivoCASSCFTermination.iterationLimit
        var assignment=Array(states.indices),overlaps=[Double](repeating:1,count:states.count),stepScale=1.0
        for iteration in 0...opt.maximumMacroIterations {
            g=try gradient(full,states);let gn=norm(g)
            history.append(.init(iteration:iteration,energyHartree:value,orbitalGradientNorm:gn,acceptedRotationNorm:stepNorm,activeCIGapHartree:nil))
            if gn<=opt.gradientTolerance {
                if previous == nil || abs(value-previous!)<=opt.energyToleranceHartree { termination = .converged;break }
                previous=value;continue
            }
            if iteration==opt.maximumMacroIterations { break }
            if evaluations>=opt.maximumEnergyEvaluations { termination = .evaluationLimit;break }
            let scale=min(stepScale,opt.maximumRotationNorm/max(gn,1e-300)),direction=g.map { -scale*$0 },slope = -scale*gn*gn
            var alpha=1.0,accepted=false
            for _ in 0..<opt.maximumLineSearchSteps {
                if evaluations>=opt.maximumEnergyEvaluations { break }
                var k=VivoQMMatrix(n,n)
                for (i,pair) in pairs.enumerated() { k[pair.first,pair.second]=alpha*direction[i];k[pair.second,pair.first] = -alpha*direction[i] }
                let delta=try VivoQMDenseAlgebra.orbitalRotation(generator:k),candidateU=try u.multiplied(by:delta)
                let (candidateFull,raw)=try spectrum(candidateU)
                let match:VivoRootFollowingResult
                if cfg.followRoots {
                    do { match=try VivoStateFollowing.match(previous:states,candidates:raw,orbitalOverlap:frame(u).transposed.multiplied(by:frame(candidateU)),coreCount:partition.doublyOccupiedCore.count,minimumOverlapSquared:cfg.minimumOverlapSquared,budget:budget) }
                    catch VivoChemistryError.convergence(_) { alpha*=0.5;continue }
                } else { match = .init(states:raw,candidateIndices:Array(raw.indices),overlapsSquared:[]) }
                let candidateValue=energy(match.states)
                if candidateValue<=value+1e-4*alpha*slope {
                    let nextGradient=try gradient(candidateFull,match.states)
                    var oldGenerator=VivoQMMatrix(n,n)
                    for (i,pair) in pairs.enumerated() { oldGenerator[pair.first,pair.second]=g[i];oldGenerator[pair.second,pair.first] = -g[i] }
                    let transported=try oldGenerator.congruence(delta),oldGradient=pairs.map { transported[$0.first,$0.second] }
                    let step=direction.map { alpha*$0 },difference=zip(nextGradient,oldGradient).map { $0-$1 }
                    let sy=zip(step,difference).reduce(0.0) { $0+$1.0*$1.1 },ss=step.reduce(0.0) { $0+$1*$1 }
                    stepScale=sy>1e-14 ? min(10,max(1e-4,ss/sy)):1
                    previous=value;value=candidateValue;u=candidateU;full=candidateFull;states=match.states
                    assignment=match.candidateIndices;overlaps=match.overlapsSquared;stepNorm=norm(step);accepted=true;break
                };alpha*=0.5
            }
            if !accepted { termination=evaluations>=opt.maximumEnergyEvaluations ? .evaluationLimit:.lineSearchFailed;break }
        }
        return .init(configuration:cfg,orbitalRotation:u,states:states,weightedEnergyHartree:value,orbitalGradient:g,
            termination:termination,iterations:history,finalRootAssignment:assignment,finalOverlapSquared:overlaps)
    }
}
