import Foundation

public struct VivoOrbitalPartition: Codable, Sendable, Equatable {
    public let fragment: [Int]
    public let bath: [Int]
    public let environment: [Int]
    public var cluster: [Int] { fragment+bath }
    public init(fragment: [Int], bath: [Int], environment: [Int]) { self.fragment=fragment;self.bath=bath;self.environment=environment }
    public func validate(orbitalCount n: Int) throws {
        let all=fragment+bath+environment
        guard !fragment.isEmpty,all.count==n,Set(all)==Set(0..<n) else { throw VivoChemistryError.invalid("F/B/E must be a disjoint complete orbital partition") }
    }
}
public struct VivoBathResult: Codable, Sendable, Equatable {
    public let candidateEnvironment: [Int]
    public let coefficients: VivoQMMatrix
    public let singularValues: [Double]
    public let retainedRank: Int
    public let discardedSquaredWeight: Double
    public let requestedDiscardedSquaredWeight: Double
}
public enum VivoBathBuilder {
    public static func fromCouplings(_ x: VivoQMMatrix, candidateEnvironment: [Int], maximumRank: Int,
                                     discardedSquaredWeight: Double = 1e-10,
                                     budget: VivoChemistryBudget = .init()) throws -> VivoBathResult {
        try budget.validate()
        guard x.rows>0,x.columns>0,x.rows==candidateEnvironment.count,
              Set(candidateEnvironment).count==x.rows,candidateEnvironment.allSatisfy({$0>=0}),
              x.rows<=budget.maximumBasisFunctions,x.columns<=budget.maximumBasisFunctions,
              x.values.allSatisfy(\.isFinite),maximumRank>=0,
              discardedSquaredWeight.isFinite,discardedSquaredWeight>=0 else { throw VivoChemistryError.invalid("bath coupling matrix/configuration") }
        _ = try budget.elements([max(x.rows,x.columns),max(x.rows,x.columns)],simultaneousArrays:8)
        let gram=try x.transposed.multiplied(by:x)
        let eig=try VivoQMDenseAlgebra.symmetricEigen(gram,tolerance:1e-14)
        let order=Array(eig.values.indices.reversed())
        let weights=order.map { max(0,eig.values[$0]) }
        let total=x.values.reduce(0){$0+$1*$1}
        let numericalFloor=1e-13*max(1,total)
        var retained=0,remaining=total
        while remaining>discardedSquaredWeight && retained<min(maximumRank,min(x.rows,x.columns)) {
            if weights[retained]<=numericalFloor { break }
            remaining=max(0,remaining-weights[retained]);retained += 1
        }
        guard remaining<=discardedSquaredWeight+1e-14*max(1,total) else {
            throw VivoChemistryError.convergence("bath discarded-weight bound cannot be met under rank/FP64 Gram precision limits")
        }
        var q=VivoQMMatrix(x.rows,retained)
        for column in 0..<retained {
            let index=order[column],sigma=sqrt(weights[column])
            for row in 0..<x.rows { for k in 0..<x.columns { q[row,column] += x[row,k]*eig.vectors[k,index]/sigma } }
        }
        if retained>0 {
            let orthogonality=try q.transposed.multiplied(by:q).adding(.identity(retained),scale:-1).frobeniusNorm
            guard orthogonality<1e-8 else { throw VivoChemistryError.convergence("bath orthogonality") }
        }
        return .init(candidateEnvironment:candidateEnvironment,coefficients:q,singularValues:weights.map(sqrt),
                     retainedRank:retained,discardedSquaredWeight:remaining,requestedDiscardedSquaredWeight:discardedSquaredWeight)
    }
    public static func fromOneRDM(_ gamma: VivoQMMatrix, fragment: [Int], candidateEnvironment: [Int],
                                  maximumRank: Int, discardedSquaredWeight: Double = 1e-10,
                                  budget: VivoChemistryBudget = .init()) throws -> VivoBathResult {
        let n=gamma.rows
        guard n>0,gamma.columns==n,!fragment.isEmpty,!candidateEnvironment.isEmpty,
              Set(fragment+candidateEnvironment).count==fragment.count+candidateEnvironment.count,
              (fragment+candidateEnvironment).allSatisfy({$0>=0 && $0<n}),gamma.values.allSatisfy(\.isFinite) else {
            throw VivoChemistryError.invalid("1-RDM bath partition")
        }
        var x=VivoQMMatrix(candidateEnvironment.count,fragment.count)
        for i in candidateEnvironment.indices { for j in fragment.indices { x[i,j]=gamma[candidateEnvironment[i],fragment[j]] } }
        return try fromCouplings(x,candidateEnvironment:candidateEnvironment,maximumRank:maximumRank,
                                 discardedSquaredWeight:discardedSquaredWeight,budget:budget)
    }
}
public struct VivoQIOWeights: Codable, Sendable, Equatable {
    public var fragmentEntropy: Double
    public var fragmentBathInformation: Double
    public var clusterEnvironmentInformation: Double
    public var oneRDMLeakage: Double
    public var cumulantPairTransfer: Double
    public init(fragmentEntropy: Double = 1, fragmentBathInformation: Double = 1,
                clusterEnvironmentInformation: Double = 1, oneRDMLeakage: Double = 1,
                cumulantPairTransfer: Double = 1) {
        self.fragmentEntropy=fragmentEntropy;self.fragmentBathInformation=fragmentBathInformation
        self.clusterEnvironmentInformation=clusterEnvironmentInformation;self.oneRDMLeakage=oneRDMLeakage
        self.cumulantPairTransfer=cumulantPairTransfer
    }
    public func validate() throws {
        guard [fragmentEntropy,fragmentBathInformation,clusterEnvironmentInformation,oneRDMLeakage,cumulantPairTransfer].allSatisfy({$0.isFinite && $0>=0}) else { throw VivoChemistryError.invalid("QIO weights") }
    }
}
public struct VivoQIOEvaluation: Codable, Sendable, Equatable {
    public let objective: Double
    public let fragmentEntropy: Double
    public let fragmentBathInformation: Double
    public let clusterEnvironmentInformation: Double
    public let oneRDMLeakageSquared: Double
    public let cumulantPairTransferSquared: Double
    public let entropyNats: [Double]
    public let mutualInformation: VivoQMMatrix
}
public enum VivoQIOObjective {
    public static func evaluate(state: VivoCIState, partition: VivoOrbitalPartition,
                                weights: VivoQIOWeights = .init(), budget: VivoChemistryBudget = .init()) throws -> VivoQIOEvaluation {
        try partition.validate(orbitalCount:state.orbitalCount);try weights.validate()
        let rdm=try VivoCIDensityMatrices.compute(state,budget:budget)
        let information=try VivoCIDensityMatrices.orbitalInformation(state,budget:budget)
        let entropy=partition.fragment.reduce(0){$0+information.entropy[$1]}
        var fb=0.0,ce=0.0,gamma=0.0,cumulant=0.0
        for f in partition.fragment { for b in partition.bath { fb += information.mutualInformation[f,b] } }
        for c in partition.cluster { for e in partition.environment { ce += information.mutualInformation[c,e] } }
        let cluster=partition.cluster.flatMap{[2*$0,2*$0+1]},environment=partition.environment.flatMap{[2*$0,2*$0+1]}
        for p in cluster { for r in environment { gamma += rdm.one[p,r]*rdm.one[p,r] } }
        for p in cluster { for q in cluster { for r in environment { for s in environment {
            let value=rdm.cumulant(p,q,r,s);cumulant += value*value
        } } } }
        let objective = -weights.fragmentEntropy*entropy-weights.fragmentBathInformation*fb
            + weights.clusterEnvironmentInformation*ce+weights.oneRDMLeakage*gamma+weights.cumulantPairTransfer*cumulant
        guard objective.isFinite else { throw VivoChemistryError.invalid("QIO nonfinite objective") }
        return .init(objective:objective,fragmentEntropy:entropy,fragmentBathInformation:fb,clusterEnvironmentInformation:ce,
                     oneRDMLeakageSquared:gamma,cumulantPairTransferSquared:cumulant,
                     entropyNats:information.entropy,mutualInformation:information.mutualInformation)
    }
}
