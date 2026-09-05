import Foundation

public struct VivoECCBathSelection: Codable, Sendable, Equatable {
    public var maximumCandidates: Int?
    public var minimumBathOrbitals: Int
    public var mutualInformationWeight: Double
    public var entropyWeight: Double
    public var cumulantWeight: Double
    public var maximumOmittedCouplingSquared: Double
    public init(maximumCandidates: Int? = nil, minimumBathOrbitals: Int = 0,
                mutualInformationWeight: Double = 1, entropyWeight: Double = 0, cumulantWeight: Double = 0,
                maximumOmittedCouplingSquared: Double = 1e-10) {
        self.maximumCandidates=maximumCandidates; self.minimumBathOrbitals=minimumBathOrbitals
        self.mutualInformationWeight=mutualInformationWeight; self.entropyWeight=entropyWeight
        self.cumulantWeight=cumulantWeight; self.maximumOmittedCouplingSquared=maximumOmittedCouplingSquared
    }
    public func validate() throws {
        guard maximumCandidates.map({$0>0}) ?? true, minimumBathOrbitals>=0,
              [mutualInformationWeight,entropyWeight,cumulantWeight,maximumOmittedCouplingSquared].allSatisfy({$0.isFinite && $0>=0}) else {
            throw VivoChemistryError.invalid("ECC bath-selection settings")
        }
    }
}
public struct VivoECCBathCandidate: Codable, Sendable, Equatable {
    public let orbital: Int
    public let score: Double
}

/// Exact change of orbital frame for an explicit real CI wavefunction. Alpha and
/// beta exterior-power matrices are contracted separately, retaining the signs
/// required by the runtime's interleaved spin order. No eigensolver is called.
public enum VivoCIOrbitalFrame {
    public static func rotated(_ state: VivoCIState, by u: VivoQMMatrix,
                               budget: VivoChemistryBudget = .init()) throws -> VivoCIState {
        try state.validate(budget:budget)
        let n=state.orbitalCount
        guard u.rows==n,u.columns==n,u.values.allSatisfy(\.isFinite),
              try u.transposed.multiplied(by:u).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {
            throw VivoChemistryError.invalid("CI orbital-frame rotation")
        }
        let dets=try VivoDirectCI.determinants(n:n,na:state.alphaElectrons,nb:state.betaElectrons,budget:budget)
        func spinMask(_ d:UInt64,_ spin:Int)->UInt64 {
            (0..<n).reduce(UInt64(0)) { $0 | (((d>>(2*$1+spin))&1)<<$1) }
        }
        func phase(_ d:UInt64)->Double {
            var parity=0
            for p in 0..<n where d & (UInt64(1)<<(2*p+1)) != 0 {
                for q in (p+1)..<n where d & (UInt64(1)<<(2*q)) != 0 { parity ^= 1 }
            }
            return parity==0 ? 1 : -1
        }
        let alpha=Set(dets.map{spinMask($0,0)}).sorted(),beta=Set(dets.map{spinMask($0,1)}).sorted()
        _ = try budget.elements([max(alpha.count,beta.count),max(alpha.count,beta.count)],simultaneousArrays:6)
        _ = try budget.elements([alpha.count,beta.count],simultaneousArrays:6)
        let ai=Dictionary(uniqueKeysWithValues:alpha.enumerated().map{($0.element,$0.offset)})
        let bi=Dictionary(uniqueKeysWithValues:beta.enumerated().map{($0.element,$0.offset)})
        var work=0
        func charge(_ count:Int) throws {
            guard count>=0,count<=budget.maximumOperatorApplications-work else { throw VivoChemistryError.resourceLimit("CI exterior-power rotation work") }
            work+=count
        }
        func exterior(_ masks:[UInt64],population:Int) throws -> VivoQMMatrix {
            let occupied=masks.map{d in (0..<n).filter{d & (UInt64(1)<<$0) != 0}}
            var result=VivoQMMatrix(masks.count,masks.count)
            for i in masks.indices { for j in masks.indices {
                try charge(max(1,population*population*population))
                if population==0 { result[i,j]=1;continue }
                var a=VivoQMMatrix(population,population)
                for p in 0..<population { for q in 0..<population { a[p,q]=u[occupied[i][p],occupied[j][q]] } }
                var value=1.0
                for k in 0..<population {
                    var pivot=k
                    for p in (k+1)..<population where abs(a[p,k])>abs(a[pivot,k]) { pivot=p }
                    if abs(a[pivot,k])<1e-15 { value=0;break }
                    if pivot != k { for q in k..<population { let v=a[k,q];a[k,q]=a[pivot,q];a[pivot,q]=v };value = -value }
                    value*=a[k,k]
                    for p in (k+1)..<population { let f=a[p,k]/a[k,k];for q in (k+1)..<population { a[p,q]-=f*a[k,q] } }
                }
                result[i,j]=value
            } }
            return result
        }
        var amplitudes=VivoQMMatrix(alpha.count,beta.count)
        for (i,d) in state.determinants.enumerated() { amplitudes[ai[spinMask(d,0)]!,bi[spinMask(d,1)]!]=phase(d)*state.coefficients[i] }
        let a=try exterior(alpha,population:state.alphaElectrons),b=try exterior(beta,population:state.betaElectrons)
        try charge(alpha.count*beta.count*(alpha.count+beta.count))
        let transformed=try a.transposed.multiplied(by:amplitudes).multiplied(by:b)
        let coefficients=dets.map{phase($0)*transformed[ai[spinMask($0,0)]!,bi[spinMask($0,1)]!]}
        let result=VivoCIState(orbitalCount:n,alphaElectrons:state.alphaElectrons,betaElectrons:state.betaElectrons,
            determinants:dets,coefficients:coefficients)
        try result.validate(budget:budget);return result
    }
}

public struct VivoECCBuiltCluster: Codable, Sendable, Equatable {
    public let fragment: VivoECCFragment
    public let coefficients: VivoQMMatrix
    public let completeFrame: VivoQMMatrix
    public let bathCandidates: [VivoECCBathCandidate]
    public let spectralBathRank: Int
    public let supplementalBathRank: Int
    public let discardedCouplingSquared: Double
    public let omittedCandidateCouplingSquared: Double
    public let physicalHamiltonian: VivoEmbeddedHamiltonian
    public let bareOneElectron: VivoQMMatrix
    public let environmentPotential: VivoQMMatrix
    public let referenceStateInFrame: VivoCIState
    public let referenceDensityInFrame: VivoSpinRDMs
    public let qio: VivoQIOEvaluation
}
public struct VivoECCEnergyContribution: Codable, Sendable, Equatable {
    public let fragmentIdentifier: String
    public let referenceContributionHartree: Double
    public let impurityContributionHartree: Double
    public var correctionHartree: Double { impurityContributionHartree-referenceContributionHartree }
}
public struct VivoECCFrameResult: Codable, Sendable, Equatable {
    public let reference: VivoCIResult
    public let referencePhysicalEnergyHartree: Double
    public let clusters: [VivoECCBuiltCluster]
    public let states: [VivoNumberMatchedState]
    public let numberMatching: VivoNumberMatchingResult?
    public let energyContributions: [VivoECCEnergyContribution]
    public let energyHartree: Double
    public let energyConvention: String
    /// For a complete partition this is the sum over non-overlapping fragments.
    /// In single-fragment mode it is the reference N plus the cluster N change;
    /// fractional inactive occupations can produce a reported nonzero defect.
    public let reconstructedElectronCount: Double
    public let electronCountDefect: Double
}

public enum VivoECCClusterConstruction {
    /// Trace of a spatial operator against a spin RDM. Supports a projected
    /// reference RDM with fractional cluster population; no fictitious integer
    /// electron count is assigned to that projection.
    static func expectation(one:VivoQMMatrix,two:[Double],rdm:VivoSpinRDMs)->Double {
        let n=one.rows
        var value=0.0
        for p in 0..<n { for q in 0..<n {
            for spin in 0..<2 { value+=one[p,q]*rdm.one[2*p+spin,2*q+spin] }
            for r in 0..<n { for s in 0..<n {
                let g=0.5*two[((p*n+q)*n+r)*n+s]
                for sigma in 0..<2 { for tau in 0..<2 { value+=g*rdm.gamma2(2*p+sigma,2*r+tau,2*q+sigma,2*s+tau) } }
            } }
        } }
        return value
    }
    public static func build(_ physical:VivoEmbeddedHamiltonian,reference:VivoCIState,fragment f:VivoECCFragment,
                             selection:VivoECCBathSelection = .init(),discardedWeight:Double=1e-10,
                             qioWeights:VivoQIOWeights = .init(),budget:VivoChemistryBudget = .init()) throws -> VivoECCBuiltCluster {
        try physical.validate(budget:budget);try reference.validate(budget:budget);try selection.validate();try qioWeights.validate()
        let n=physical.orbitalCount
        guard reference.orbitalCount==n,reference.alphaElectrons==physical.alphaElectrons,reference.betaElectrons==physical.betaElectrons,
              !f.identifier.isEmpty,!f.orbitals.isEmpty,Set(f.orbitals).count==f.orbitals.count,
              f.orbitals.allSatisfy({$0>=0 && $0<n}),f.maximumBathOrbitals>=0,f.maximumBathOrbitals<=n-f.orbitals.count,
              selection.minimumBathOrbitals<=f.maximumBathOrbitals,discardedWeight.isFinite,discardedWeight>=0 else {
            throw VivoChemistryError.invalid("ECC correlated reference/fragment/bath contract")
        }
        let rdm=try VivoCIDensityMatrices.compute(reference,budget:budget),density=rdm.spatialOne
        for p in 0..<n { for q in 0..<n where abs(rdm.one[2*p,2*q]-rdm.one[2*p+1,2*q+1])>1e-7 {
            throw VivoChemistryError.unsupported("ECC spin-polarized environment needs a spin-dependent Hamiltonian boundary")
        } }
        let environment=(0..<n).filter{!f.orbitals.contains($0)}
        let info=try VivoCIDensityMatrices.orbitalInformation(reference,budget:budget)
        // Explicit loop types avoid an expensive map/sort inference chain in
        // the Apple Swift compiler. Accumulation and ordering stay unchanged.
        var ranked: [VivoECCBathCandidate] = []
        ranked.reserveCapacity(environment.count)
        for j in environment {
            var information = 0.0
            for i in f.orbitals { information += info.mutualInformation[i,j] }
            var score = selection.mutualInformationWeight * information
            score += selection.entropyWeight * info.entropy[j]
            if selection.cumulantWeight > 0 {
                for i in f.orbitals { for k in f.orbitals {
                    for spin in 0..<2 { for tau in 0..<2 {
                        let cumulant = rdm.cumulant(2*i+spin,2*j+tau,2*k+spin,2*j+tau)
                        score += selection.cumulantWeight * abs(cumulant)
                    } }
                } }
            }
            guard score.isFinite else { throw VivoChemistryError.invalid("ECC bath score overflow") }
            ranked.append(VivoECCBathCandidate(orbital: j, score: score))
        }
        ranked.sort { left, right in
            if left.score == right.score { return left.orbital < right.orbital }
            return left.score > right.score
        }
        let maximumCandidates = selection.maximumCandidates ?? ranked.count
        let candidates: [Int] = ranked.prefix(maximumCandidates).map { $0.orbital }
        var omitted=0.0
        for e in environment where !candidates.contains(e) { for p in f.orbitals { omitted+=density[e,p]*density[e,p] } }
        guard omitted<=selection.maximumOmittedCouplingSquared else { throw VivoChemistryError.convergence("ECC candidate truncation exceeds omitted 1-RDM coupling bound") }
        let bath:VivoBathResult?
        if candidates.isEmpty { bath=nil }
        else { bath=try VivoBathBuilder.fromOneRDM(density,fragment:f.orbitals,candidateEnvironment:candidates,
            maximumRank:f.maximumBathOrbitals,discardedSquaredWeight:discardedWeight,budget:budget) }
        var columns:[[Double]]=[]
        for p in f.orbitals { var v=[Double](repeating:0,count:n);v[p]=1;columns.append(v) }
        func orthogonalize(_ vector:[Double])->[Double]? {
            var v=vector
            for _ in 0..<2 { for c in columns { let d=zip(c,v).reduce(0.0){$0+$1.0*$1.1};for i in 0..<n { v[i]-=d*c[i] } } }
            let norm=v.reduce(0.0){hypot($0,$1)}
            guard norm>1e-10 else{return nil}
            v=v.map{$0/norm}
            if let pivot=v.indices.max(by:{abs(v[$0])<abs(v[$1])}),v[pivot]<0 { v=v.map{-$0} }
            return v
        }
        if let bath { for j in 0..<bath.retainedRank {
            var v=[Double](repeating:0,count:n)
            for (i,p) in candidates.enumerated(){v[p]=bath.coefficients[i,j]}
            guard let c=orthogonalize(v) else {throw VivoChemistryError.convergence("ECC bath lost rank")};columns.append(c)
        } }
        let spectralRank=bath?.retainedRank ?? 0
        for candidate in ranked where columns.count<f.orbitals.count+selection.minimumBathOrbitals {
            var v=[Double](repeating:0,count:n);v[candidate.orbital]=1
            if let c=orthogonalize(v){columns.append(c)}
        }
        let k=columns.count,supplemental=k-f.orbitals.count-spectralRank
        guard k>=f.orbitals.count+selection.minimumBathOrbitals,f.clusterAlphaElectrons>=0,f.clusterBetaElectrons>=0,
              f.clusterAlphaElectrons<=k,f.clusterBetaElectrons<=k,f.clusterAlphaElectrons==f.clusterBetaElectrons else {
            throw VivoChemistryError.invalid("ECC declared impurity sector or supplemental bath capacity")
        }
        for p in 0..<n where columns.count<n {
            var v=[Double](repeating:0,count:n);v[p]=1
            if let c=orthogonalize(v){columns.append(c)}
        }
        guard columns.count==n else {throw VivoChemistryError.convergence("ECC complete orbital frame")}
        var complete=VivoQMMatrix(n,n),c=VivoQMMatrix(n,k)
        for p in 0..<n { for q in 0..<n {complete[p,q]=columns[q][p];if q<k{c[p,q]=columns[q][p]}} }
        let complement=try VivoQMMatrix.identity(n).adding(c.multiplied(by:c.transposed),scale:-1)
        let env=try density.congruence(complement)
        var potential=VivoQMMatrix(n,n)
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            potential[p,q]+=env[r,s]*(physical.eri(p,q,r,s)-0.5*physical.eri(p,r,q,s))
        } } } }
        let bare=try physical.oneElectron.congruence(c),localPotential=try potential.congruence(c)
        let h=VivoEmbeddedHamiltonian(orbitalIdentifiers:(0..<k).map{f.identifier+"-cluster-\($0)"},
            alphaElectrons:f.clusterAlphaElectrons,betaElectrons:f.clusterBetaElectrons,
            oneElectron:try bare.adding(localPotential),twoElectron:try VivoOrbitalTransform.transformERI(physical.twoElectron,inputCount:n,coefficients:c,budget:budget),
            constantEnergyHartree:0,energyReference:"physical impurity plus spin-unpolarized correlated-density environment field; scalar excluded")
        try h.validate(budget:budget)
        let localState=try VivoCIOrbitalFrame.rotated(reference,by:complete,budget:budget)
        let localRDM=try VivoCIDensityMatrices.compute(localState,budget:budget)
        let qio=try VivoQIOObjective.evaluate(state:localState,partition:.init(fragment:Array(0..<f.orbitals.count),
            bath:Array(f.orbitals.count..<k),environment:Array(k..<n)),weights:qioWeights,budget:budget)
        return .init(fragment:f,coefficients:c,completeFrame:complete,bathCandidates:ranked,spectralBathRank:spectralRank,
            supplementalBathRank:supplemental,discardedCouplingSquared:bath?.discardedSquaredWeight ?? 0,
            omittedCandidateCouplingSquared:omitted,physicalHamiltonian:h,bareOneElectron:bare,environmentPotential:localPotential,
            referenceStateInFrame:localState,referenceDensityInFrame:localRDM,qio:qio)
    }
    public static func evaluate(_ physical:VivoEmbeddedHamiltonian,reference:VivoCIResult,fragments:[VivoECCFragment],
                                partitionMatching:VivoNumberMatchingConfiguration?,selection:VivoECCBathSelection = .init(),
                                discardedWeight:Double=1e-10,qioWeights:VivoQIOWeights = .init(),
                                budget:VivoChemistryBudget = .init()) throws -> VivoECCFrameResult {
        try physical.validate(budget:budget);try reference.state.validate(budget:budget)
        guard !fragments.isEmpty,Set(fragments.map(\.identifier)).count==fragments.count else {throw VivoChemistryError.invalid("ECC fragment identifiers")}
        if partitionMatching != nil {
            let indices=fragments.flatMap(\.orbitals)
            guard indices.count==physical.orbitalCount,Set(indices)==Set(0..<physical.orbitalCount) else {throw VivoChemistryError.invalid("ECC democratic partition must cover every orbital once")}
        } else if fragments.count != 1 {throw VivoChemistryError.invalid("single-fragment ECC requires exactly one fragment")}
        let clusters=try fragments.map{try build(physical,reference:reference.state,fragment:$0,selection:selection,
            discardedWeight:discardedWeight,qioWeights:qioWeights,budget:budget)}
        let matching:VivoNumberMatchingResult?,states:[VivoNumberMatchedState]
        if let cfg=partitionMatching {
            let inputs=clusters.map{cluster -> VivoNumberMatchedFragment in
                let k=cluster.coefficients.columns;var p=VivoQMMatrix(k,k)
                for i in cluster.fragment.orbitals.indices{p[i,i]=1}
                return .init(identifier:cluster.fragment.identifier,hamiltonian:cluster.physicalHamiltonian,fragmentProjector:p)
            }
            let result=try VivoFragmentNumberMatcher.solve(fragments:inputs,targetPopulation:Double(physical.alphaElectrons+physical.betaElectrons),configuration:cfg,budget:budget)
            guard result.converged else {throw VivoChemistryError.convergence("ECC partition particle-number feedback: \(result.termination)")}
            matching=result;states=result.states
        } else {
            let cluster=clusters[0]
            let ci=try VivoDirectCI.solve(cluster.physicalHamiltonian,budget:budget).roots[0]
            let rdm=try VivoCIDensityMatrices.compute(ci.state,budget:budget)
            let population=cluster.fragment.orbitals.indices.reduce(0.0){$0+rdm.one[2*$1,2*$1]+rdm.one[2*$1+1,2*$1+1]}
            states=[.init(identifier:cluster.fragment.identifier,fragmentPopulation:population,biasedCI:ci,physicalClusterEnergyHartree:ci.energyHartree)]
            matching=nil
        }
        let referenceRDM=try VivoCIDensityMatrices.compute(reference.state,budget:budget)
        let base=try referenceRDM.energy(of:physical,budget:budget)
        var contributions:[VivoECCEnergyContribution]=[]
        for (i,cluster) in clusters.enumerated() {
            let k=cluster.coefficients.columns,nf=cluster.fragment.orbitals.count
            let density=try VivoCIDensityMatrices.compute(states[i].biasedCI.state,budget:budget)
            var one=cluster.physicalHamiltonian.oneElectron,two=cluster.physicalHamiltonian.twoElectron
            if partitionMatching != nil {
                let effective=try cluster.bareOneElectron.adding(cluster.environmentPotential,scale:0.5)
                for p in 0..<k { for q in 0..<k {
                    one[p,q]=0.5*Double((p<nf ? 1:0)+(q<nf ? 1:0))*effective[p,q]
                    for r in 0..<k { for s in 0..<k {
                        two[((p*k+q)*k+r)*k+s]*=0.25*Double((p<nf ? 1:0)+(q<nf ? 1:0)+(r<nf ? 1:0)+(s<nf ? 1:0))
                    } }
                } }
            }
            let low=expectation(one:one,two:two,rdm:cluster.referenceDensityInFrame)
            let high=expectation(one:one,two:two,rdm:density)
            guard low.isFinite,high.isFinite else {throw VivoChemistryError.convergence("ECC partition-energy contraction")}
            contributions.append(.init(fragmentIdentifier:cluster.fragment.identifier,referenceContributionHartree:low,impurityContributionHartree:high))
        }
        let energy=base+contributions.reduce(0.0){$0+$1.correctionHartree},ne=Double(physical.alphaElectrons+physical.betaElectrons)
        let reconstructed:Double
        if let matching {reconstructed=ne+matching.populationResidual}
        else {
            let cluster=clusters[0],k=cluster.coefficients.columns
            let projected=(0..<(2*k)).reduce(0.0){$0+cluster.referenceDensityInFrame.one[$1,$1]}
            reconstructed=ne+Double(cluster.physicalHamiltonian.alphaElectrons+cluster.physicalHamiltonian.betaElectrons)-projected
        }
        guard energy.isFinite,reconstructed.isFinite else {throw VivoChemistryError.convergence("ECC reconstructed total")}
        return .init(reference:reference,referencePhysicalEnergyHartree:base,clusters:clusters,states:states,numberMatching:matching,
            energyContributions:contributions,energyHartree:energy,
            energyConvention:partitionMatching == nil ? "reference physical energy + full impurity expectation correction; inactive reference restored" :
                "reference-restored democratic cluster correction; physical scalar once; half environment field in fragment one-body energy; no auxiliary potentials",
            reconstructedElectronCount:reconstructed,electronCountDefect:reconstructed-ne)
    }
}
