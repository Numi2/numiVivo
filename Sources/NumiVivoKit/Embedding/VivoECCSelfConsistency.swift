import Foundation

public struct VivoECCFragment:Codable,Sendable,Equatable {
    public let identifier:String
    public let orbitals:[Int]
    public let maximumBathOrbitals:Int
    public let clusterAlphaElectrons:Int
    public let clusterBetaElectrons:Int
    public init(identifier:String,orbitals:[Int],maximumBathOrbitals:Int,clusterAlphaElectrons:Int,clusterBetaElectrons:Int) {
        self.identifier=identifier;self.orbitals=orbitals;self.maximumBathOrbitals=maximumBathOrbitals
        self.clusterAlphaElectrons=clusterAlphaElectrons;self.clusterBetaElectrons=clusterBetaElectrons
    }
}
/// Dimensionless Hermitian one-/two-body operator; its multiplier is in Hartree.
/// Operators must be sums of fragment-local terms in the declared orbital frame.
public struct VivoECCCorrelationOperator:Codable,Sendable,Equatable {
    public let identifier:String
    public let one:VivoQMMatrix
    public let two:[Double]
    public init(identifier:String,one:VivoQMMatrix,two:[Double]) { self.identifier=identifier;self.one=one;self.two=two }
}
public struct VivoECCSelfConsistencyConfiguration:Codable,Sendable,Equatable {
    public var referenceMethod:VivoCIMethod
    public var maximumIterations:Int
    public var maximumEvaluations:Int
    public var momentTolerance:Double
    public var bathDiscardedWeight:Double
    public var differenceStepHartree:Double
    public var maximumPotentialHartree:Double
    public var numberMatching:VivoNumberMatchingConfiguration
    public init(referenceMethod:VivoCIMethod = .cisd,maximumIterations:Int=30,maximumEvaluations:Int=500,
                momentTolerance:Double=1e-7,bathDiscardedWeight:Double=1e-10,differenceStepHartree:Double=1e-4,
                maximumPotentialHartree:Double=10,numberMatching:VivoNumberMatchingConfiguration = .init()) {
        self.referenceMethod=referenceMethod;self.maximumIterations=maximumIterations;self.maximumEvaluations=maximumEvaluations
        self.momentTolerance=momentTolerance;self.bathDiscardedWeight=bathDiscardedWeight
        self.differenceStepHartree=differenceStepHartree;self.maximumPotentialHartree=maximumPotentialHartree;self.numberMatching=numberMatching
    }
}
public struct VivoECCSelfConsistencyIteration:Codable,Sendable,Equatable {
    public let iteration:Int
    public let momentResidualNorm:Double
    public let electronResidual:Double
    public let chemicalPotentialHartree:Double
    public let referenceEnergyHartree:Double
}
public struct VivoECCSelfConsistencyResult:Codable,Sendable,Equatable {
    public let converged:Bool
    public let termination:String
    public let correlationPotentialHartree:[Double]
    public let momentResiduals:[Double]
    public let referenceMoments:[Double]
    public let fragmentMoments:[Double]
    public let numberMatching:VivoNumberMatchingResult
    public let clusterCoefficients:[VivoQMMatrix]
    public let iterations:[VivoECCSelfConsistencyIteration]
    public let evaluations:Int
    public let scope:String
    /// Present for newly executed cycles; absent in older archival results.
    public let frame:VivoECCFrameResult?
}
/// Selected-moment, multi-fragment correlated-reference self-consistency.
/// Each trial recomputes the reference, correlated-density baths, environment
/// mean-field contraction, impurity solutions and global chemical potential.
/// This extends the paper's single-fragment workflow, which omits matching.
/// It is NOT a claim that a finite selected operator set matches every 2-RDM
/// element, nor that converged matching certifies a reaction barrier.
public enum VivoECCSelfConsistency {
    private struct Evaluation {
        let reference:[Double];let fragment:[Double];let residual:[Double]
        let number:VivoNumberMatchingResult;let coefficients:[VivoQMMatrix];let referenceEnergy:Double
        let frame:VivoECCFrameResult
    }
    public static func solve(_ physical:VivoEmbeddedHamiltonian,fragments:[VivoECCFragment],
                             operators:[VivoECCCorrelationOperator],initialPotentialHartree:[Double]?=nil,
                             configuration cfg:VivoECCSelfConsistencyConfiguration = .init(),
                             bathSelection:VivoECCBathSelection = .init(),qioWeights:VivoQIOWeights = .init(),
                             budget:VivoChemistryBudget = .init()) throws -> VivoECCSelfConsistencyResult {
        try physical.validate(budget:budget);try cfg.numberMatching.validate();try bathSelection.validate();try qioWeights.validate()
        let n=physical.orbitalCount,m=operators.count
        guard physical.alphaElectrons==physical.betaElectrons,(1...16).contains(fragments.count),(1...256).contains(m),
              Set(fragments.map(\.identifier)).count==fragments.count,Set(operators.map(\.identifier)).count==m,
              fragments.allSatisfy({!$0.identifier.isEmpty && !$0.orbitals.isEmpty && $0.maximumBathOrbitals>=0 && $0.maximumBathOrbitals<=n && $0.clusterAlphaElectrons>=0 && $0.clusterAlphaElectrons==$0.clusterBetaElectrons}),
              fragments.flatMap(\.orbitals).count==n,Set(fragments.flatMap(\.orbitals))==Set(0..<n),
              (1...500).contains(cfg.maximumIterations),(1...100000).contains(cfg.maximumEvaluations),
              [cfg.momentTolerance,cfg.differenceStepHartree,cfg.maximumPotentialHartree].allSatisfy({$0.isFinite && $0>0}),
              cfg.bathDiscardedWeight.isFinite,cfg.bathDiscardedWeight>=0 else { throw VivoChemistryError.invalid("ECC fragment partition, spin or iteration contract") }
        _ = try budget.elements([m,m],simultaneousArrays:8)
        var owner=[Int](repeating:0,count:n)
        for (k,f) in fragments.enumerated() { for p in f.orbitals { owner[p]=k } }
        func operatorHamiltonian(_ op:VivoECCCorrelationOperator)->VivoEmbeddedHamiltonian {
            .init(orbitalIdentifiers:physical.orbitalIdentifiers,alphaElectrons:physical.alphaElectrons,betaElectrons:physical.betaElectrons,
                oneElectron:op.one,twoElectron:op.two,constantEnergyHartree:0,energyReference:"dimensionless selected correlation observable")
        }
        for op in operators {
            try operatorHamiltonian(op).validate(budget:budget)
            guard !op.identifier.isEmpty,op.one.values.contains(where:{$0 != 0}) || op.two.contains(where:{$0 != 0}) else { throw VivoChemistryError.invalid("empty correlation operator") }
            // Remove the uniform one-body gauge; mu has its own number equation.
            guard abs((0..<n).reduce(0.0) { $0+op.one[$1,$1] })<1e-10 else { throw VivoChemistryError.invalid("one-body correlation potential must be traceless") }
            for p in 0..<n { for q in 0..<n {
                if op.one[p,q] != 0,owner[p] != owner[q] { throw VivoChemistryError.unsupported("cross-fragment one-body fitting operator") }
                for r in 0..<n { for s in 0..<n where op.two[((p*n+q)*n+r)*n+s] != 0 {
                    guard owner[p]==owner[q],owner[p]==owner[r],owner[p]==owner[s] else { throw VivoChemistryError.unsupported("cross-fragment two-body fitting operator") }
                } }
            } }
        }
        var potential=initialPotentialHartree ?? [Double](repeating:0,count:m),evaluations=0
        guard potential.count==m,potential.allSatisfy({$0.isFinite && abs($0)<=cfg.maximumPotentialHartree}) else { throw VivoChemistryError.invalid("initial correlation potential") }
        func evaluate(_ x:[Double]) throws -> Evaluation {
            guard evaluations<cfg.maximumEvaluations else { throw VivoChemistryError.resourceLimit("ECC aggregate evaluation count") };evaluations+=1
            var one=physical.oneElectron,two=physical.twoElectron
            for i in 0..<m {
                one=try one.adding(operators[i].one,scale:x[i])
                for k in two.indices { two[k]+=x[i]*operators[i].two[k] }
            }
            let referenceHamiltonian=VivoEmbeddedHamiltonian(orbitalIdentifiers:physical.orbitalIdentifiers,
                alphaElectrons:physical.alphaElectrons,betaElectrons:physical.betaElectrons,oneElectron:one,twoElectron:two,
                constantEnergyHartree:physical.constantEnergyHartree,energyReference:physical.energyReference+"; auxiliary correlation potential")
            let reference:VivoCIResult
            if cfg.referenceMethod == .fci { reference=try VivoDirectCI.solve(referenceHamiltonian,budget:budget).roots[0] }
            else { reference=try VivoConfigurationInteraction.solve(referenceHamiltonian,method:.cisd,budget:budget) }
            let rdm=try VivoCIDensityMatrices.compute(reference.state,budget:budget)
            let frame=try VivoECCClusterConstruction.evaluate(physical,reference:reference,fragments:fragments,
                partitionMatching:cfg.numberMatching,selection:bathSelection,discardedWeight:cfg.bathDiscardedWeight,
                qioWeights:qioWeights,budget:budget)
            guard let number=frame.numberMatching else { throw VivoChemistryError.invalid("partition cycle omitted particle matching") }
            let referenceMoments=try operators.map { try rdm.energy(of:operatorHamiltonian($0),budget:budget) }
            var fragmentMoments=[Double](repeating:0,count:m)
            for (index,cluster) in frame.clusters.enumerated() {
                let f=cluster.fragment.orbitals,k=cluster.coefficients.columns
                let localRDM=try VivoCIDensityMatrices.compute(frame.states[index].biasedCI.state,budget:budget)
                for (i,op) in operators.enumerated() {
                    var one=VivoQMMatrix(k,k),two=[Double](repeating:0,count:try budget.elements([k,k,k,k],simultaneousArrays:3))
                    for p in f.indices { for q in f.indices {
                        one[p,q]=op.one[f[p],f[q]]
                        for r in f.indices { for s in f.indices { two[((p*k+q)*k+r)*k+s]=op.two[((f[p]*n+f[q])*n+f[r])*n+f[s]] } }
                    } }
                    fragmentMoments[i]+=VivoECCClusterConstruction.expectation(one:one,two:two,rdm:localRDM)
                }
            }
            return .init(reference:referenceMoments,fragment:fragmentMoments,residual:zip(referenceMoments,fragmentMoments).map { $0-$1 },
                number:number,coefficients:frame.clusters.map(\.coefficients),referenceEnergy:reference.energyHartree,frame:frame)
        }
        func norm(_ x:[Double])->Double { x.reduce(0.0) { hypot($0,$1) } }
        var current=try evaluate(potential),history:[VivoECCSelfConsistencyIteration]=[],converged=false,termination="iterationLimit"
        for iteration in 0...cfg.maximumIterations {
            let error=norm(current.residual)
            history.append(.init(iteration:iteration,momentResidualNorm:error,electronResidual:current.number.populationResidual,
                chemicalPotentialHartree:current.number.chemicalPotentialHartree,referenceEnergyHartree:current.referenceEnergy))
            if error<=cfg.momentTolerance { converged=true;termination="convergedSelectedMoments";break }
            if iteration==cfg.maximumIterations { break }
            if evaluations+2*m+1>cfg.maximumEvaluations { termination="evaluationLimit";break }
            var jacobian=VivoQMMatrix(m,m)
            for j in 0..<m {
                var plus=potential,minus=potential;plus[j]=min(cfg.maximumPotentialHartree,plus[j]+cfg.differenceStepHartree);minus[j]=max(-cfg.maximumPotentialHartree,minus[j]-cfg.differenceStepHartree)
                let a=try evaluate(plus),b=try evaluate(minus)
                for i in 0..<m { jacobian[i,j]=(a.residual[i]-b.residual[i])/(plus[j]-minus[j]) }
            }
            var normal=try jacobian.transposed.multiplied(by:jacobian)
            let gradient=try jacobian.transposed.multiplied(by:VivoQMMatrix(rows:m,columns:1,values:current.residual)).values
            // Bounded Levenberg step handles ill-conditioned selected moments.
            // Regularization does not relax the final physical residual test.
            for i in 0..<m { normal[i,i]+=1e-8 }
            var step=try VivoQMDenseAlgebra.solve(normal,rhs:gradient.map { -$0 })
            let bound=min(1,0.5/max(norm(step),1e-300));step=step.map { bound*$0 }
            var alpha=1.0,accepted=false
            for _ in 0..<16 {
                if evaluations>=cfg.maximumEvaluations { break }
                let candidate=zip(potential,step).map { $0+alpha*$1 }
                if candidate.contains(where:{abs($0)>cfg.maximumPotentialHartree}) { alpha*=0.5;continue }
                let trial=try evaluate(candidate)
                if norm(trial.residual)<error*(1-1e-4*alpha) {
                    potential=candidate;current=trial;accepted=true;break
                };alpha*=0.5
            }
            if !accepted { termination=evaluations>=cfg.maximumEvaluations ? "evaluationLimit":"lineSearchFailed";break }
        }
        return .init(converged:converged,termination:termination,correlationPotentialHartree:potential,momentResiduals:current.residual,
            referenceMoments:current.reference,fragmentMoments:current.fragment,numberMatching:current.number,clusterCoefficients:current.coefficients,
            iterations:history,evaluations:evaluations,
            scope:"correlated-reference and impurity feedback; selected fragment-local one/two-body matching; reference-restored democratic energy; no claim of a globally N-representable reconstructed 2-RDM",frame:current.frame)
    }
}
