import Foundation

public struct VivoDoubleFactor: Codable, Sendable, Equatable {
    public let weight: Double
    public let matrix: VivoQMMatrix
    public let orbitalRotation: VivoQMMatrix
    public let orbitalEigenvalues: [Double]
    public let numberShift: Double
}

public struct VivoDoubleFactorizedHamiltonian: Codable, Sendable, Equatable {
    public static let schema = "numivivo.org/double-factorized-hamiltonian/v1"
    public let schema: String
    public let orbitalIdentifiers: [String]
    public let alphaElectrons: Int
    public let betaElectrons: Int
    public let effectiveOneElectron: VivoQMMatrix
    public let constantEnergyHartree: Double
    public let factors: [VivoDoubleFactor]
    public let discardedERIFrobeniusNorm: Double
    public let sourceEnergyReference: String
    public let symmetryOptimized: Bool
    public var orbitalCount: Int { orbitalIdentifiers.count }
    public var electronCount: Int { alphaElectrons + betaElectrons }

    public func validate(budget: VivoChemistryBudget = .init()) throws {
        let n = orbitalCount
        guard schema == Self.schema, n > 0, n <= budget.maximumBasisFunctions,
              Set(orbitalIdentifiers).count == n, alphaElectrons >= 0, betaElectrons >= 0,
              alphaElectrons <= n, betaElectrons <= n, effectiveOneElectron.rows == n,
              effectiveOneElectron.columns == n, effectiveOneElectron.values.allSatisfy(\.isFinite),
              constantEnergyHartree.isFinite, discardedERIFrobeniusNorm.isFinite,
              discardedERIFrobeniusNorm >= 0, !sourceEnergyReference.isEmpty else {
            throw VivoChemistryError.invalid("double-factorized Hamiltonian metadata")
        }
        for factor in factors {
            guard factor.weight.isFinite, factor.matrix.rows == n, factor.matrix.columns == n,
                  factor.matrix.values.allSatisfy(\.isFinite), factor.orbitalRotation.rows == n,
                  factor.orbitalRotation.columns == n, factor.orbitalEigenvalues.count == n,
                  factor.orbitalEigenvalues.allSatisfy(\.isFinite), factor.numberShift.isFinite else {
                throw VivoChemistryError.invalid("double-factor dimensions/values")
            }
            for i in 0..<n { for j in 0..<i where abs(factor.matrix[i,j]-factor.matrix[j,i]) > 1e-9 {
                throw VivoChemistryError.invalid("double-factor matrix symmetry")
            } }
            let identity = try VivoQMMatrix.identity(n)
            guard try factor.orbitalRotation.transposed.multiplied(by:factor.orbitalRotation)
                .adding(identity,scale:-1).frobeniusNorm < 1e-8 else {
                throw VivoChemistryError.invalid("double-factor orbital rotation")
            }
        }
    }
}

public struct VivoDoubleFactorizationConfiguration: Codable, Sendable, Equatable {
    public var maximumOuterFactors: Int
    public var maximumDiscardedERIFrobeniusNorm: Double
    public var eigenvalueTolerance: Double
    public init(maximumOuterFactors: Int = 4096, maximumDiscardedERIFrobeniusNorm: Double = 1e-10,
                eigenvalueTolerance: Double = 1e-12) {
        self.maximumOuterFactors=maximumOuterFactors
        self.maximumDiscardedERIFrobeniusNorm=maximumDiscardedERIFrobeniusNorm
        self.eigenvalueTolerance=eigenvalueTolerance
    }
    public func validate() throws {
        guard maximumOuterFactors > 0, maximumDiscardedERIFrobeniusNorm.isFinite,
              maximumDiscardedERIFrobeniusNorm >= 0, eigenvalueTolerance.isFinite,
              eigenvalueTolerance > 0 else { throw VivoChemistryError.invalid("double-factorization configuration") }
    }
}

public enum VivoDoubleFactorization {
    public static func factorize(_ h: VivoEmbeddedHamiltonian,
                                 configuration cfg: VivoDoubleFactorizationConfiguration = .init(),
                                 budget: VivoChemistryBudget = .init()) throws -> VivoDoubleFactorizedHamiltonian {
        try h.validate(budget:budget); try cfg.validate()
        let n=h.orbitalCount,pairs=n*n
        _ = try budget.elements([pairs,pairs],simultaneousArrays:4)
        var supermatrix=VivoQMMatrix(pairs,pairs)
        for p in 0..<n { for q in 0..<n { for r in 0..<n { for s in 0..<n {
            supermatrix[p*n+q,r*n+s]=h.eri(p,q,r,s)
        } } } }
        let eig=try VivoQMDenseAlgebra.symmetricEigen(supermatrix,tolerance:cfg.eigenvalueTolerance,maximumSweeps:256)
        let order=eig.values.indices.sorted { abs(eig.values[$0]) > abs(eig.values[$1]) }
        var discardedSquared=eig.values.reduce(0.0){$0+$1*$1},retained:[Int]=[]
        let targetSquared=cfg.maximumDiscardedERIFrobeniusNorm*cfg.maximumDiscardedERIFrobeniusNorm
        for index in order where retained.count < cfg.maximumOuterFactors {
            if discardedSquared <= targetSquared { break }
            if abs(eig.values[index]) <= cfg.eigenvalueTolerance { break }
            retained.append(index); discardedSquared=max(0,discardedSquared-eig.values[index]*eig.values[index])
        }
        let discarded=sqrt(discardedSquared)
        guard discarded <= cfg.maximumDiscardedERIFrobeniusNorm + 10*cfg.eigenvalueTolerance else {
            throw VivoChemistryError.resourceLimit("outer-factor limit cannot meet ERI Frobenius error target")
        }
        var effective=h.oneElectron
        for p in 0..<n { for s in 0..<n { for q in 0..<n {
            effective[p,s] -= 0.5*h.eri(p,q,q,s)
        } } }
        var factors:[VivoDoubleFactor]=[]
        for index in retained {
            var matrix=VivoQMMatrix(n,n)
            for p in 0..<n { for q in 0..<n { matrix[p,q]=eig.vectors[p*n+q,index] } }
            let antisymmetry=try matrix.adding(matrix.transposed,scale:-1).frobeniusNorm
            guard antisymmetry < 1e-7 else { throw VivoChemistryError.invalid("ERI supermatrix eigenvector is not a symmetric density factor") }
            matrix=try matrix.adding(matrix.transposed).scaled(0.5)
            let inner=try VivoQMDenseAlgebra.symmetricEigen(matrix,tolerance:cfg.eigenvalueTolerance,maximumSweeps:128)
            factors.append(.init(weight:eig.values[index],matrix:matrix,orbitalRotation:inner.vectors,
                                 orbitalEigenvalues:inner.values,numberShift:0))
        }
        let result=VivoDoubleFactorizedHamiltonian(schema:VivoDoubleFactorizedHamiltonian.schema,
            orbitalIdentifiers:h.orbitalIdentifiers,alphaElectrons:h.alphaElectrons,betaElectrons:h.betaElectrons,
            effectiveOneElectron:effective,constantEnergyHartree:h.constantEnergyHartree,factors:factors,
            discardedERIFrobeniusNorm:discarded,sourceEnergyReference:h.energyReference,symmetryOptimized:false)
        try result.validate(budget:budget); return result
    }

    public static func optimizeNumberSymmetry(_ input: VivoDoubleFactorizedHamiltonian,
                                              budget: VivoChemistryBudget = .init()) throws -> VivoDoubleFactorizedHamiltonian {
        try input.validate(budget:budget)
        let n=input.orbitalCount,N=Double(input.electronCount)
        guard N>0 else { throw VivoChemistryError.invalid("number-symmetry optimization requires nonzero fixed particle number") }
        var one=input.effectiveOneElectron,constant=input.constantEnergyHartree,newFactors=input.factors

        func objective(oneBody: VivoQMMatrix, factors: [VivoDoubleFactor]) -> Double {
            let oneNorm=oneBody.values.reduce(0){$0+abs($1)}
            let twoNorm=0.5*factors.reduce(0.0) { total,factor in
                let l1=factor.orbitalEigenvalues.reduce(0){$0+abs($1)}
                return total+abs(factor.weight)*l1*l1
            }
            return oneNorm+twoNorm
        }

        var currentObjective=objective(oneBody:one,factors:newFactors)
        for factorIndex in newFactors.indices {
            let original=newFactors[factorIndex]
            let sorted=original.orbitalEigenvalues.sorted()
            var candidates=Array(Set([0.0,sorted[sorted.count/2]]+sorted)).sorted()
            if sorted.count>1 {
                for i in 0..<(sorted.count-1) { candidates.append(0.5*(sorted[i]+sorted[i+1])) }
            }
            var bestMu=0.0,bestObjective=currentObjective,bestOne=one,bestConstant=constant,bestFactor=original
            for mu in candidates where mu.isFinite {
                if abs(mu)<1e-18 { continue }
                var shifted=original.matrix
                for p in 0..<n { shifted[p,p] -= mu }
                var candidateOne=one
                for p in 0..<n { for q in 0..<n {
                    candidateOne[p,q] += original.weight*mu*N*shifted[p,q]
                } }
                let candidateConstant=constant+0.5*original.weight*mu*mu*N*N
                let candidateFactor=VivoDoubleFactor(weight:original.weight,matrix:shifted,
                    orbitalRotation:original.orbitalRotation,
                    orbitalEigenvalues:original.orbitalEigenvalues.map{$0-mu},numberShift:mu)
                var factorSet=newFactors;factorSet[factorIndex]=candidateFactor
                let candidateObjective=objective(oneBody:candidateOne,factors:factorSet)
                if candidateObjective < bestObjective-1e-12*max(1,bestObjective) {
                    bestMu=mu;bestObjective=candidateObjective;bestOne=candidateOne
                    bestConstant=candidateConstant;bestFactor=candidateFactor
                }
            }
            if bestMu != 0 {
                one=bestOne;constant=bestConstant;newFactors[factorIndex]=bestFactor;currentObjective=bestObjective
            }
        }
        let result=VivoDoubleFactorizedHamiltonian(schema:VivoDoubleFactorizedHamiltonian.schema,
            orbitalIdentifiers:input.orbitalIdentifiers,alphaElectrons:input.alphaElectrons,betaElectrons:input.betaElectrons,
            effectiveOneElectron:one,constantEnergyHartree:constant,factors:newFactors,
            discardedERIFrobeniusNorm:input.discardedERIFrobeniusNorm,sourceEnergyReference:input.sourceEnergyReference,
            symmetryOptimized:true)
        try result.validate(budget:budget)
        let before=try VivoDoubleFactorizationNorm.estimate(input).normalization
        let after=try VivoDoubleFactorizationNorm.estimate(result).normalization
        guard after<=before+1e-10*max(1,before) else { throw VivoChemistryError.convergence("number-symmetry optimization increased global normalization") }
        return result
    }
}

public struct VivoBlockEncodingNorm: Codable, Sendable, Equatable {
    public let oneBodyOneNorm: Double
    public let factorizedTwoBodyOneNorm: Double
    public let normalization: Double
}
public enum VivoDoubleFactorizationNorm {
    public static func estimate(_ h: VivoDoubleFactorizedHamiltonian) throws -> VivoBlockEncodingNorm {
        try h.validate()
        let one=h.effectiveOneElectron.values.reduce(0){$0+abs($1)}
        let two=0.5*h.factors.reduce(0.0) { total,factor in
            let l1=factor.orbitalEigenvalues.reduce(0){$0+abs($1)}
            return total+abs(factor.weight)*l1*l1
        }
        let normalization=one+two
        guard normalization.isFinite,normalization>0 else { throw VivoChemistryError.invalid("block-encoding normalization") }
        return .init(oneBodyOneNorm:one,factorizedTwoBodyOneNorm:two,normalization:normalization)
    }
}

public struct VivoQubitizationCostAssumptions: Codable, Sendable, Equatable {
    public let tGatesPerQuery: UInt64
    public let blockEncodingAncillaQubits: Int
    public let phaseEstimationAncillaQubits: Int
    public let modelIdentifier: String
    public init(tGatesPerQuery: UInt64, blockEncodingAncillaQubits: Int,
                phaseEstimationAncillaQubits: Int, modelIdentifier: String) {
        self.tGatesPerQuery=tGatesPerQuery;self.blockEncodingAncillaQubits=blockEncodingAncillaQubits
        self.phaseEstimationAncillaQubits=phaseEstimationAncillaQubits;self.modelIdentifier=modelIdentifier
    }
}
public struct VivoQPEResourceEstimate: Codable, Sendable, Equatable {
    public let targetEnergyErrorHartree: Double
    public let blockEncodingNormalization: Double
    public let walkQueries: UInt64
    public let logicalTCount: UInt64
    public let logicalQubits: Int
    public let assumptions: VivoQubitizationCostAssumptions
    public let status: String
}
public enum VivoQubitizedQPE {
    public static func estimate(_ h: VivoDoubleFactorizedHamiltonian, targetEnergyErrorHartree epsilon: Double,
                                assumptions: VivoQubitizationCostAssumptions) throws -> VivoQPEResourceEstimate {
        try h.validate()
        guard epsilon.isFinite,epsilon>0,assumptions.tGatesPerQuery>0,
              assumptions.blockEncodingAncillaQubits>=0,assumptions.phaseEstimationAncillaQubits>=0,
              !assumptions.modelIdentifier.isEmpty else { throw VivoChemistryError.invalid("QPE target/cost assumptions") }
        let norm=try VivoDoubleFactorizationNorm.estimate(h)
        let raw=ceil(Double.pi*norm.normalization/(2*epsilon))
        guard raw.isFinite,raw>0,raw<=Double(UInt64.max) else { throw VivoChemistryError.resourceLimit("QPE query count overflow") }
        let queries=UInt64(raw),product=queries.multipliedReportingOverflow(by:assumptions.tGatesPerQuery)
        guard !product.overflow else { throw VivoChemistryError.resourceLimit("QPE T-count overflow") }
        let qubits=2*h.orbitalCount+assumptions.blockEncodingAncillaQubits+assumptions.phaseEstimationAncillaQubits
        return .init(targetEnergyErrorHartree:epsilon,blockEncodingNormalization:norm.normalization,walkQueries:queries,
                     logicalTCount:product.partialValue,logicalQubits:qubits,assumptions:assumptions,
                     status:"algorithmic logical-resource model; excludes surface-code overhead, routing, magic-state factory and hardware runtime")
    }
}
