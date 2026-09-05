import Foundation

public struct VivoOrbitalTransport: Codable, Sendable, Equatable {
    public let rotation: VivoQMMatrix
    public let minimumSubspaceSingularValue: Double
    public let groups: [[Int]]
}
public enum VivoOrbitalPathTransport {
    public static func align(referenceCoefficients a: VivoQMMatrix, currentCoefficients b: VivoQMMatrix,
                             crossAOOverlap s: VivoQMMatrix, groups: [[Int]], minimumSingularValue: Double = 0.5) throws -> VivoOrbitalTransport {
        let n=a.columns,indices=groups.flatMap{$0}
        guard n>0,b.columns==n,s.rows==a.rows,s.columns==b.rows,
              indices.count==n,Set(indices)==Set(0..<n),minimumSingularValue.isFinite,
              minimumSingularValue>0,minimumSingularValue<=1 else { throw VivoChemistryError.invalid("orbital path transport dimensions/groups") }
        let overlap=try a.transposed.multiplied(by:s).multiplied(by:b)
        var rotation=try VivoQMMatrix.identity(n),minimum=Double.infinity
        for group in groups where !group.isEmpty {
            let m=group.count;var block=VivoQMMatrix(m,m)
            for i in 0..<m { for j in 0..<m { block[i,j]=overlap[group[i],group[j]] } }
            let gram=try block.transposed.multiplied(by:block),eig=try VivoQMDenseAlgebra.symmetricEigen(gram,tolerance:1e-13)
            guard eig.values.allSatisfy({$0>=minimumSingularValue*minimumSingularValue && $0<=1+1e-6}) else {
                throw VivoChemistryError.convergence("reaction path lost an active/fragment/environment subspace; expand a shared path partition")
            }
            minimum=min(minimum,sqrt(eig.values[0]))
            var inverseRoot=VivoQMMatrix(m,m)
            for i in 0..<m { for j in 0..<m { for k in 0..<m {
                inverseRoot[i,j] += eig.vectors[i,k]*eig.vectors[j,k]/sqrt(eig.values[k])
            } } }
            let local=try inverseRoot.multiplied(by:block.transposed)
            for i in 0..<m { for j in 0..<m { rotation[group[i],group[j]]=local[i,j] } }
        }
        guard try rotation.transposed.multiplied(by:rotation).adding(.identity(n),scale:-1).frobeniusNorm<1e-8 else {
            throw VivoChemistryError.convergence("path transport orthogonality")
        }
        return .init(rotation:rotation,minimumSubspaceSingularValue:minimum,groups:groups)
    }
}
public struct VivoPathQIOConfiguration: Codable, Sendable, Equatable {
    public var maximumIterations: Int
    public var maximumObjectiveEvaluations: Int
    public var gradientTolerance: Double
    public var finiteDifferenceStep: Double
    public var initialStep: Double
    public var minimumCIStateGapHartree: Double
    public init(maximumIterations: Int = 24, maximumObjectiveEvaluations: Int = 512,
                gradientTolerance: Double = 1e-5, finiteDifferenceStep: Double = 1e-4, initialStep: Double = 0.2, minimumCIStateGapHartree: Double = 1e-8) {
        self.maximumIterations=maximumIterations;self.maximumObjectiveEvaluations=maximumObjectiveEvaluations
        self.gradientTolerance=gradientTolerance;self.finiteDifferenceStep=finiteDifferenceStep;self.initialStep=initialStep
        self.minimumCIStateGapHartree=minimumCIStateGapHartree
    }
}
public struct VivoPathQIOResult: Codable, Sendable, Equatable {
    public let sharedRotation: VivoQMMatrix
    public let partition: VivoOrbitalPartition
    public let evaluations: [VivoQIOEvaluation]
    public let acceptedObjectives: [Double]
    public let objectiveEvaluations: Int
    public let converged: Bool
    public let termination: String
}
public enum VivoPathConsistentQIO {
    public static func optimize(alignedPoints points: [VivoEmbeddedHamiltonian], partition: VivoOrbitalPartition,
                                weights: VivoQIOWeights = .init(), configuration cfg: VivoPathQIOConfiguration = .init(),
                                budget: VivoChemistryBudget = .init()) throws -> VivoPathQIOResult {
        guard let first=points.first,points.count<=128,
              cfg.maximumIterations>0,cfg.maximumObjectiveEvaluations>0,
              [cfg.gradientTolerance,cfg.finiteDifferenceStep,cfg.initialStep,cfg.minimumCIStateGapHartree].allSatisfy({$0.isFinite && $0>0}) else {
            throw VivoChemistryError.invalid("joint path QIO configuration")
        }
        let n=first.orbitalCount
        try partition.validate(orbitalCount:n);try weights.validate()
        for point in points {
            try point.validate(budget:budget)
            guard point.orbitalIdentifiers==first.orbitalIdentifiers,point.alphaElectrons==first.alphaElectrons,
                  point.betaElectrons==first.betaElectrons,point.energyReference==first.energyReference else {
                throw VivoChemistryError.invalid("path QIO requires aligned orbital identities, equal spin sector and common energy convention")
            }
        }
        let pairs=(0..<n).flatMap { i in ((i+1)..<n).map{(i,$0)} }
        var parameters=[Double](repeating:0,count:pairs.count),calls=0
        func rotation(_ p: [Double]) throws -> VivoQMMatrix {
            var k=VivoQMMatrix(n,n)
            for (index,pair) in pairs.enumerated() { k[pair.0,pair.1]=p[index];k[pair.1,pair.0] = -p[index] }
            return try VivoQMDenseAlgebra.orbitalRotation(generator:k)
        }
        func evaluate(_ p: [Double]) throws -> (Double,[VivoQIOEvaluation]) {
            guard calls<cfg.maximumObjectiveEvaluations else { throw VivoChemistryError.resourceLimit("joint QIO objective-evaluation budget") }
            calls += 1;let u=try rotation(p)
            let evaluations=try points.map { point in
                let h=try point.rotated(by:u,budget:budget),ci=try VivoConfigurationInteraction.solve(h,budget:budget)
                if let gap=ci.nextStateGapHartree,gap<cfg.minimumCIStateGapHartree {
                    throw VivoChemistryError.unsupported("near-degenerate path state requires a root-tracked/state-averaged reference")
                }
                return try VivoQIOObjective.evaluate(state:ci.state,partition:partition,weights:weights,budget:budget)
            }
            return (evaluations.reduce(0){$0+$1.objective}/Double(evaluations.count),evaluations)
        }
        var (value,evaluations)=try evaluate(parameters),accepted=[value],converged=false,termination="iteration-limit"
        for _ in 0..<cfg.maximumIterations {
            if calls+2*pairs.count+1>cfg.maximumObjectiveEvaluations { termination="evaluation-budget";break }
            var gradient=[Double](repeating:0,count:pairs.count)
            for i in parameters.indices {
                var plus=parameters,minus=parameters;plus[i] += cfg.finiteDifferenceStep;minus[i] -= cfg.finiteDifferenceStep
                gradient[i]=(try evaluate(plus).0-evaluate(minus).0)/(2*cfg.finiteDifferenceStep)
            }
            let normSquared=gradient.reduce(0){$0+$1*$1}
            if sqrt(normSquared)<cfg.gradientTolerance { converged=true;termination="gradient-tolerance";break }
            var step=cfg.initialStep,acceptedStep=false
            for _ in 0..<16 {
                if calls>=cfg.maximumObjectiveEvaluations { termination="evaluation-budget";break }
                let candidate=zip(parameters,gradient).map{$0.0-step*$0.1}
                let (next,nextEvaluations)=try evaluate(candidate)
                if next<=value-1e-4*step*normSquared {
                    parameters=candidate;value=next;evaluations=nextEvaluations;accepted.append(next);acceptedStep=true;break
                }
                step *= 0.5
            }
            if !acceptedStep { if termination != "evaluation-budget" { termination="line-search-stalled" };break }
        }
        return .init(sharedRotation:try rotation(parameters),partition:partition,evaluations:evaluations,acceptedObjectives:accepted,
                     objectiveEvaluations:calls,converged:converged,termination:termination)
    }
}
