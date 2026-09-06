import Foundation

public struct VivoSelectedCIConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var maximumDeterminants: Int
    public var selectionBatchSize: Int
    public var minimumSelectionContributionHartree: Double
    public var pt2ToleranceHartree: Double
    public var eigenResidualTolerance: Double
    public var minimumDenominatorHartree: Double
    public var maximumDavidsonSubspace: Int
    public init(maximumIterations: Int = 32, maximumDeterminants: Int = 4096,
                selectionBatchSize: Int = 128, minimumSelectionContributionHartree: Double = 1e-9,
                pt2ToleranceHartree: Double = 1e-6, eigenResidualTolerance: Double = 1e-9,
                minimumDenominatorHartree: Double = 1e-5, maximumDavidsonSubspace: Int = 48) {
        self.maximumIterations=maximumIterations;self.maximumDeterminants=maximumDeterminants
        self.selectionBatchSize=selectionBatchSize;self.minimumSelectionContributionHartree=minimumSelectionContributionHartree
        self.pt2ToleranceHartree=pt2ToleranceHartree;self.eigenResidualTolerance=eigenResidualTolerance
        self.minimumDenominatorHartree=minimumDenominatorHartree;self.maximumDavidsonSubspace=maximumDavidsonSubspace
    }
    public func validate(budget: VivoChemistryBudget) throws {
        guard (1...1000).contains(maximumIterations),(1...budget.maximumDeterminants).contains(maximumDeterminants),
              (1...maximumDeterminants).contains(selectionBatchSize),minimumSelectionContributionHartree.isFinite,
              minimumSelectionContributionHartree>=0,pt2ToleranceHartree.isFinite,pt2ToleranceHartree>0,
              eigenResidualTolerance.isFinite,eigenResidualTolerance>0,minimumDenominatorHartree.isFinite,
              minimumDenominatorHartree>0,(4...512).contains(maximumDavidsonSubspace) else {
            throw VivoChemistryError.invalid("selected-CI iteration, determinant, threshold, denominator or Davidson contract")
        }
    }
}

public struct VivoSelectedCIIteration: Codable, Sendable, Equatable {
    public let iteration:Int
    public let determinantCount:Int
    public let variationalEnergyHartree:Double
    public let eigenResidual:Double
    public let externalCandidateCount:Int
    public let selectedCount:Int
    public let pt2CorrectionHartree:Double
    public let largestExternalContributionHartree:Double
    public let intruderCandidateCount:Int
}

public struct VivoSelectedCIResult: Codable, Sendable, Equatable {
    public let converged:Bool
    public let variationalEnergyHartree:Double
    public let pt2CorrectedEnergyHartree:Double
    public let pt2CorrectionHartree:Double
    public let state:VivoCIState
    public let fullSectorDimension:Int
    public let selectedDeterminantCount:Int
    public let eigenResidual:Double
    public let iterations:[VivoSelectedCIIteration]
    public let operatorApplications:Int
    public let method:String
}

/// Deterministic CIPSI-style selected CI. Only selected determinants are stored
/// in the variational Hamiltonian. Connected external determinants are generated
/// from the wavefunction, ranked by Epstein-Nesbet second-order contributions,
/// and admitted in bounded batches. PT2 is a remainder diagnostic, not a rigorous
/// error bar. Near-zero denominators are selected preferentially and reported.
public enum VivoSelectedCI {
    public static let method="adaptive determinant-selected variational CI with matrix-free Davidson and Epstein-Nesbet PT2 selection; PT2 is a diagnostic remainder, not a certified error bound or DMRG replacement"

    private static func sectorDimension(_ n:Int,_ a:Int,_ b:Int) throws -> Int {
        func choose(_ n:Int,_ input:Int) throws -> Int {
            guard input>=0,input<=n else{return 0};let k=min(input,n-input);if k==0{return 1};var c=1
            for i in 1...k {let x=c.multipliedReportingOverflow(by:n-k+i);guard !x.overflow else{throw VivoChemistryError.resourceLimit("selected-CI sector dimension overflow")};c=x.partialValue/i}
            return c
        }
        let ca=try choose(n,a),cb=try choose(n,b),x=ca.multipliedReportingOverflow(by:cb)
        guard !x.overflow else{throw VivoChemistryError.resourceLimit("selected-CI sector dimension overflow")};return x.partialValue
    }
    private static func hartreeFockDeterminant(_ h:VivoEmbeddedHamiltonian)->UInt64 {
        var d:UInt64=0
        for p in 0..<h.alphaElectrons {d |= UInt64(1) << (2*p)}
        for p in 0..<h.betaElectrons {d |= UInt64(1) << (2*p+1)}
        return d
    }
    private static func diagonal(_ det:UInt64,_ h:VivoEmbeddedHamiltonian)->Double {
        let occ=(0..<(2*h.orbitalCount)).filter{det & (UInt64(1)<<$0) != 0};var e=0.0
        func integral(_ p:Int,_ q:Int,_ r:Int,_ s:Int)->Double {
            (p%2==r%2 && q%2==s%2 ? h.eri(p/2,r/2,q/2,s/2):0) -
            (p%2==s%2 && q%2==r%2 ? h.eri(p/2,s/2,q/2,r/2):0)
        }
        for p in occ {e += h.oneElectron[p/2,p/2];for q in occ {e += 0.5*integral(p,q,p,q)}}
        return e
    }
    private static func connected(_ det:UInt64,_ h:VivoEmbeddedHamiltonian,
                                  visit:(UInt64,Double)throws->Void) throws {
        let modes=2*h.orbitalCount,occ=(0..<modes).filter{det & (UInt64(1)<<$0) != 0},vir=(0..<modes).filter{det & (UInt64(1)<<$0) == 0}
        func integral(_ p:Int,_ q:Int,_ r:Int,_ s:Int)->Double {
            (p%2==r%2 && q%2==s%2 ? h.eri(p/2,r/2,q/2,s/2):0) -
            (p%2==s%2 && q%2==r%2 ? h.eri(p/2,s/2,q/2,r/2):0)
        }
        for i in occ {for a in vir where i%2==a%2 {
            guard let action=vivoApplyFermions(det,[.init(mode:i,creation:false),.init(mode:a,creation:true)]) else{continue}
            var value=h.oneElectron[a/2,i/2]
            for j in occ where j != i {value += integral(a,j,i,j)}
            if value != 0 {try visit(action.0,Double(action.1)*value)}
        }}
        if occ.count>=2 && vir.count>=2 {for ii in 0..<(occ.count-1) {for jj in (ii+1)..<occ.count {
            let i=occ[ii],j=occ[jj]
            for aa in 0..<(vir.count-1) {for bb in (aa+1)..<vir.count {
                let a=vir[aa],b=vir[bb]
                guard i%2+j%2==a%2+b%2,
                      let action=vivoApplyFermions(det,[.init(mode:i,creation:false),.init(mode:j,creation:false),.init(mode:b,creation:true),.init(mode:a,creation:true)]) else{continue}
                let value=Double(action.1)*integral(a,b,i,j);if value != 0 {try visit(action.0,value)}
            }}
        }}}
    }

    private static func solveProjected(_ h:VivoEmbeddedHamiltonian,determinants:[UInt64],initial:[Double]?,
                                       cfg:VivoSelectedCIConfiguration,budget:VivoChemistryBudget,
                                       work:inout Int) throws -> (energy:Double,coefficients:[Double],residual:Double) {
        let action=try VivoDirectHamiltonian(h,determinants:determinants,budget:budget),d=determinants.count
        let capacity=min(d,cfg.maximumDavidsonSubspace)
        func dot(_ a:[Double],_ b:[Double])->Double{zip(a,b).reduce(0){$0+$1.0*$1.1}}
        func orth(_ source:[Double],_ basis:[[Double]])->[Double]?{var q=source;for _ in 0..<2{for b in basis{let p=dot(q,b);for i in q.indices{q[i]-=p*b[i]}}};let n=q.reduce(0){hypot($0,$1)};return n>1e-12 ? q.map{$0/n}:nil}
        var seed=initial ?? [Double](repeating:0,count:d)
        if initial==nil {seed[action.diagonal.indices.min(by:{action.diagonal[$0]<action.diagonal[$1]})!]=1}
        guard seed.count==d,let q=orth(seed,[]) else{throw VivoChemistryError.invalid("selected-CI Davidson seed")}
        var basis=[q],images=[try action.apply(q,work:&work)]
        for _ in 0..<max(20,cfg.maximumIterations*4) {
            let m=basis.count;var p=VivoQMMatrix(m,m)
            for i in 0..<m {for j in 0...i {let x=0.5*(dot(basis[i],images[j])+dot(basis[j],images[i]));p[i,j]=x;p[j,i]=x}}
            let eig=try VivoQMDenseAlgebra.symmetricEigen(p,tolerance:1e-14,maximumSweeps:128)
            var v=[Double](repeating:0,count:d),av=v
            for j in 0..<m {for i in 0..<d {v[i]+=eig.vectors[j,0]*basis[j][i];av[i]+=eig.vectors[j,0]*images[j][i]}}
            let r=(0..<d).map{av[$0]-eig.values[0]*v[$0]},rn=r.reduce(0){hypot($0,$1)}
            if rn<=cfg.eigenResidualTolerance {if let k=v.indices.max(by:{abs(v[$0])<abs(v[$1])}),v[k]<0{v=v.map { -$0 }};return(eig.values[0]+h.constantEnergyHartree,v,rn)}
            let correction=(0..<d).map{i->Double in let gap=eig.values[0]-action.diagonal[i];return r[i]/((gap<0 ? -1.0:1.0)*max(abs(gap),cfg.minimumDenominatorHartree))}
            guard let next=orth(correction,basis) else{throw VivoChemistryError.convergence("selected-CI Davidson stagnation")}
            if basis.count>=capacity {basis=[v];images=[av]}
            basis.append(next);images.append(try action.apply(next,work:&work))
        }
        throw VivoChemistryError.convergence("selected-CI Davidson iteration bound")
    }

    public static func solve(_ h:VivoEmbeddedHamiltonian,configuration cfg:VivoSelectedCIConfiguration = .init(),
                             initialDeterminants:[UInt64]? = nil,budget:VivoChemistryBudget = .init()) throws -> VivoSelectedCIResult {
        try h.validate(budget:budget);try cfg.validate(budget:budget)
        guard h.orbitalCount<=31 else{throw VivoChemistryError.resourceLimit("selected CI uses UInt64 spin determinants")}
        let full=try sectorDimension(h.orbitalCount,h.alphaElectrons,h.betaElectrons)
        var determinants=initialDeterminants ?? [hartreeFockDeterminant(h)]
        determinants=Array(Set(determinants)).sorted()
        guard !determinants.isEmpty,determinants.count<=cfg.maximumDeterminants else{throw VivoChemistryError.invalid("selected-CI initial determinant population")}
        let probe=VivoCIState(orbitalCount:h.orbitalCount,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
            determinants:determinants,coefficients:[Double](repeating:1/sqrt(Double(determinants.count)),count:determinants.count))
        try probe.validate(budget:budget)
        var work=0,previous:[UInt64:Double]=[:],history:[VivoSelectedCIIteration]=[],final:(Double,[Double],Double)?
        var lastPT2=Double.infinity
        for iteration in 1...cfg.maximumIterations {
            let initial=determinants.map{previous[$0] ?? 0},seed=initial.contains(where:{$0 != 0}) ? initial:nil
            let solved=try solveProjected(h,determinants:determinants,initial:seed,cfg:cfg,budget:budget,work:&work)
            final=solved;previous=Dictionary(uniqueKeysWithValues:zip(determinants,solved.1))
            let selected=Set(determinants),electronic=solved.0-h.constantEnergyHartree
            var coupling:[UInt64:Double]=[:]
            for (det,c) in zip(determinants,solved.1) where abs(c)>1e-15 {
                try connected(det,h) { external,value in
                    guard work<budget.maximumOperatorApplications else{throw VivoChemistryError.resourceLimit("selected-CI external coupling work")};work+=1
                    if !selected.contains(external){coupling[external,default:0]+=c*value}
                }
            }
            var ranked:[(det:UInt64,contribution:Double,intruder:Bool)]=[];var pt2=0.0,intruders=0,largest=0.0
            ranked.reserveCapacity(coupling.count)
            for (det,v) in coupling where v != 0 {
                let raw=electronic-diagonal(det,h),intruder=abs(raw)<cfg.minimumDenominatorHartree
                let denom=intruder ? (raw<0 ? -cfg.minimumDenominatorHartree:cfg.minimumDenominatorHartree):raw
                let contribution=v*v/denom;pt2+=contribution;largest=max(largest,abs(contribution));if intruder{intruders+=1}
                ranked.append((det,contribution,intruder))
            }
            ranked.sort{a,b in let x=abs(a.contribution),y=abs(b.contribution);return x==y ? a.det<b.det:x>y}
            let room=cfg.maximumDeterminants-determinants.count
            let eligible=ranked.filter{$0.intruder || abs($0.contribution)>=cfg.minimumSelectionContributionHartree}
            let chosen=Array(eligible.prefix(min(room,cfg.selectionBatchSize)))
            history.append(.init(iteration:iteration,determinantCount:determinants.count,variationalEnergyHartree:solved.0,
                eigenResidual:solved.2,externalCandidateCount:coupling.count,selectedCount:chosen.count,
                pt2CorrectionHartree:pt2,largestExternalContributionHartree:largest,intruderCandidateCount:intruders))
            lastPT2=pt2
            if abs(pt2)<=cfg.pt2ToleranceHartree && intruders==0 {
                let state=VivoCIState(orbitalCount:h.orbitalCount,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
                    determinants:determinants,coefficients:solved.1);try state.validate(budget:budget)
                return .init(converged:true,variationalEnergyHartree:solved.0,pt2CorrectedEnergyHartree:solved.0+pt2,
                    pt2CorrectionHartree:pt2,state:state,fullSectorDimension:full,selectedDeterminantCount:determinants.count,
                    eigenResidual:solved.2,iterations:history,operatorApplications:work,method:method)
            }
            guard iteration<cfg.maximumIterations,!chosen.isEmpty,room>0 else{break}
            determinants.append(contentsOf:chosen.map(\.det));determinants.sort()
        }
        guard let solved=final else{throw VivoChemistryError.convergence("selected-CI produced no variational state")}
        let state=VivoCIState(orbitalCount:h.orbitalCount,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
            determinants:determinants,coefficients:determinants.map{previous[$0] ?? 0});try state.validate(budget:budget)
        return .init(converged:false,variationalEnergyHartree:solved.0,pt2CorrectedEnergyHartree:solved.0+lastPT2,
            pt2CorrectionHartree:lastPT2,state:state,fullSectorDimension:full,selectedDeterminantCount:determinants.count,
            eigenResidual:solved.2,iterations:history,operatorApplications:work,method:method)
    }
}
