import Foundation

public struct VivoQIOPartition: Codable, Sendable, Equatable {
    public var fragment: [Int]
    public var bath: [Int]
    public var environment: [Int]
    public init(fragment: [Int], bath: [Int], environment: [Int]) {
        self.fragment = fragment; self.bath = bath; self.environment = environment
    }
    public func validate(orbitalCount n: Int) throws {
        let all = fragment+bath+environment
        guard !fragment.isEmpty, all.count == n, Set(all) == Set(0..<n) else {
            throw VivoQMError.invalid("QIO requires one disjoint, exhaustive F/B/E partition shared by all path points")
        }
    }
}

public struct VivoQIOWeights: Codable, Sendable {
    public var fragmentEntropy: Double
    public var fragmentBathInformation: Double
    public var clusterEnvironmentInformation: Double
    public var densityLeakage: Double
    public var cumulantLeakage: Double
    public init(fragmentEntropy: Double = 1, fragmentBathInformation: Double = 1,
                clusterEnvironmentInformation: Double = 1, densityLeakage: Double = 1,
                cumulantLeakage: Double = 1) {
        self.fragmentEntropy = fragmentEntropy; self.fragmentBathInformation = fragmentBathInformation
        self.clusterEnvironmentInformation = clusterEnvironmentInformation
        self.densityLeakage = densityLeakage; self.cumulantLeakage = cumulantLeakage
    }
    public func validate() throws {
        let all = [fragmentEntropy,fragmentBathInformation,clusterEnvironmentInformation,densityLeakage,cumulantLeakage]
        guard all.allSatisfy({ $0.isFinite && (0...1e6).contains($0) }), all.contains(where: { $0 > 0 }) else {
            throw VivoQMError.invalid("nonnegative QIO weights with at least one positive term required")
        }
    }
}

public struct VivoQIOLeakage: Codable, Sendable {
    public let fragmentEntropy: Double
    public let fragmentBathInformation: Double
    public let clusterEnvironmentInformation: Double
    public let densitySquaredNorm: Double
    /// Explicit convention: sum |Lambda[p,q,r,s]|^2 for spin modes p,r in C
    /// and q,s in E. This is a defined cross-pair block, not an unspecified
    /// spatial-RDM permutation or every mixed-index cumulant block.
    public let crossPairCumulantSquaredNorm: Double
    public let objective: Double
    public static func evaluate(rdm: VivoFermionRDM, information: VivoOrbitalInformation,
                                partition p: VivoQIOPartition, weights w: VivoQIOWeights) throws -> Self {
        try rdm.validate(); try w.validate()
        let n = rdm.modeCount/2
        try p.validate(orbitalCount:n)
        try information.mutualInformation.requireSymmetric()
        guard information.singleOrbitalEntropy.count == n,
              information.singleOrbitalEntropy.allSatisfy({ $0.isFinite && $0 >= 0 }),
              information.mutualInformation.rows == n,
              information.mutualInformation.values.allSatisfy({ $0 >= 0 }) else { throw VivoQMError.invalid("QIO information tensor dimensions") }
        let cluster = p.fragment+p.bath
        let f = p.fragment.reduce(0.0) { $0+information.singleOrbitalEntropy[$1] }
        var fb = 0.0, ce = 0.0
        for i in p.fragment { for j in p.bath { fb += information.mutualInformation[i,j] } }
        for i in cluster { for j in p.environment { ce += information.mutualInformation[i,j] } }
        let cModes = cluster.flatMap { [2*$0,2*$0+1] }, eModes = p.environment.flatMap { [2*$0,2*$0+1] }
        var gamma = 0.0, lambda = 0.0
        for i in cModes { for j in eModes {
            gamma += pow(rdm.oneParticle[i,j],2)
            for k in cModes { for l in eModes { lambda += pow(rdm.cumulant(i,j,k,l),2) } }
        } }
        let objective = -w.fragmentEntropy*f-w.fragmentBathInformation*fb
            + w.clusterEnvironmentInformation*ce+w.densityLeakage*gamma+w.cumulantLeakage*lambda
        guard objective.isFinite else { throw VivoQMError.convergence("non-finite leakage objective") }
        return .init(fragmentEntropy:f,fragmentBathInformation:fb,clusterEnvironmentInformation:ce,
                     densitySquaredNorm:gamma,crossPairCumulantSquaredNorm:lambda,objective:objective)
    }
}

public struct VivoQIOPathPoint: Codable, Sendable {
    public var identifier: String
    public var hamiltonian: VivoEmbeddedHamiltonian
    public var weight: Double
    /// Required except on the first/reference point. Matrix entries must come
    /// from actual cross-geometry orbital overlaps, not matching array lengths.
    public var overlapToReference: VivoQMMatrix?
    public init(identifier: String, hamiltonian: VivoEmbeddedHamiltonian,
                weight: Double = 1, overlapToReference: VivoQMMatrix? = nil) {
        self.identifier = identifier; self.hamiltonian = hamiltonian
        self.weight = weight; self.overlapToReference = overlapToReference
    }
}

public struct VivoQIOConfiguration: Codable, Sendable {
    public var maximumIterations: Int
    public var maximumEvaluations: Int
    public var finiteDifferenceStep: Double
    public var gradientTolerance: Double
    public var maximumParameterStep: Double
    public var minimumOverlapSingularValue: Double
    public var minimumGroundStateGapHartree: Double
    public var solver: VivoDeterminantConfiguration
    public init(maximumIterations: Int = 20, maximumEvaluations: Int = 256,
                finiteDifferenceStep: Double = 1e-4, gradientTolerance: Double = 1e-6,
                maximumParameterStep: Double = 0.2, minimumOverlapSingularValue: Double = 0.5,
                minimumGroundStateGapHartree: Double = 1e-7,
                solver: VivoDeterminantConfiguration = .init()) {
        self.maximumIterations = maximumIterations; self.maximumEvaluations = maximumEvaluations
        self.finiteDifferenceStep = finiteDifferenceStep; self.gradientTolerance = gradientTolerance
        self.maximumParameterStep = maximumParameterStep; self.minimumOverlapSingularValue = minimumOverlapSingularValue
        self.minimumGroundStateGapHartree = minimumGroundStateGapHartree; self.solver = solver
    }
    public func validate() throws {
        try solver.validate()
        guard solver.method == .fci, (0...200).contains(maximumIterations), (1...4096).contains(maximumEvaluations),
              finiteDifferenceStep.isFinite, (1e-6...1e-2).contains(finiteDifferenceStep),
              gradientTolerance.isFinite, (1e-10...1e-2).contains(gradientTolerance),
              maximumParameterStep.isFinite, (1e-4...1).contains(maximumParameterStep),
              minimumOverlapSingularValue.isFinite, (1e-6...1).contains(minimumOverlapSingularValue),
              minimumGroundStateGapHartree.isFinite, (1e-10...0.1).contains(minimumGroundStateGapHartree) else {
            throw VivoQMError.invalid("path-QIO budget, tolerances, or solver; this implementation requires small parent-space FCI")
        }
    }
}

public struct VivoQIOPathRequest: Codable, Sendable {
    public var points: [VivoQIOPathPoint]
    public var partition: VivoQIOPartition
    public var weights: VivoQIOWeights
    public var rotationPairs: [[Int]]
    public var configuration: VivoQIOConfiguration
    public init(points: [VivoQIOPathPoint], partition: VivoQIOPartition,
                weights: VivoQIOWeights = .init(), rotationPairs: [[Int]],
                configuration: VivoQIOConfiguration = .init()) {
        self.points = points; self.partition = partition; self.weights = weights
        self.rotationPairs = rotationPairs; self.configuration = configuration
    }
}

public enum VivoQIOStopReason: String, Codable, Sendable {
    case gradientConverged, evaluationBudget, iterationBudget, lineSearchStalled, noRotationDegreesOfFreedom
}
public struct VivoQIOPointResult: Codable, Sendable {
    public let identifier: String
    public let hamiltonian: VivoEmbeddedHamiltonian
    public let fullParentEnergyHartree: Double
    public let leakage: VivoQIOLeakage
    public let transportSingularValues: [Double]
}
public struct VivoQIOPathResult: Codable, Sendable {
    public let partition: VivoQIOPartition
    public let sharedRotation: VivoQMMatrix
    public let parameters: [Double]
    public let initialObjective: Double
    public let acceptedObjectives: [Double]
    public let evaluations: Int
    public let stopReason: VivoQIOStopReason
    public let points: [VivoQIOPointResult]
    public var converged: Bool { stopReason == .gradientConverged }
}

public enum VivoPathConsistentQIO {
    /// A bounded, native small-space implementation, with central-difference
    /// gradients and Armijo backtracking, NOT CMA-ES. All points use one shared
    /// U after overlap transport. Full-parent FCI is used to obtain exact orbital
    /// entropies. It is not a scalable ECC-DMET implementation or a barrier result.
    public static func optimize(_ request: VivoQIOPathRequest) throws -> VivoQIOPathResult {
        let config = request.configuration
        try config.validate(); try request.weights.validate()
        guard let first = request.points.first, (1...16).contains(request.points.count),
              Set(request.points.map(\.identifier)).count == request.points.count,
              request.points.allSatisfy({ !$0.identifier.isEmpty && $0.weight.isFinite && (1e-12...1e12).contains($0.weight) }),
              first.overlapToReference == nil else { throw VivoQMError.invalid("QIO path identity, weights, or reference point") }
        try first.hamiltonian.validate()
        let n = first.hamiltonian.orbitalCount
        guard n <= 12 else { throw VivoQMError.capacity("path-QIO parent is limited to the small CI solver") }
        try request.partition.validate(orbitalCount:n)
        var aligned: [VivoEmbeddedHamiltonian] = [], singularValues: [[Double]] = []
        for (index,point) in request.points.enumerated() {
            try point.hamiltonian.validate()
            guard point.hamiltonian.orbitalCount == n,
                  point.hamiltonian.alphaElectrons == first.hamiltonian.alphaElectrons,
                  point.hamiltonian.betaElectrons == first.hamiltonian.betaElectrons else {
                throw VivoQMError.invalid("path points use different orbital dimensions or electron sectors")
            }
            if index == 0 { aligned.append(point.hamiltonian); singularValues.append(Array(repeating:1,count:n)); continue }
            guard let overlap = point.overlapToReference, overlap.rows == n, overlap.columns == n else {
                throw VivoQMError.invalid("every non-reference point requires cross-geometry orbital overlaps")
            }
            let transport = try VivoOrbitalSpace.alignment(overlap:overlap,minimumSingularValue:config.minimumOverlapSingularValue)
            aligned.append(try point.hamiltonian.rotated(by:transport.rotation,orbitalIDs:first.hamiltonian.orbitalIDs))
            singularValues.append(transport.singularValues)
        }
        let pairs = request.rotationPairs, count = pairs.count
        var parameters = Array(repeating:0.0,count:count), evaluations = 0
        let totalWeight = request.points.reduce(0.0) { $0+$1.weight }
        func evaluate(_ parameters: [Double]) throws -> (Double,VivoQMMatrix,[VivoQIOPointResult]) {
            guard evaluations < config.maximumEvaluations else { throw VivoQMError.capacity("internal QIO evaluation budget") }
            evaluations += 1
            let rotation = try VivoOrbitalSpace.rotation(dimension:n,pairs:pairs,parameters:parameters)
            var results: [VivoQIOPointResult] = [], objective = 0.0
            for (index,h) in aligned.enumerated() {
                let rotated = try h.rotated(by:rotation)
                let ci = try VivoDeterminantSolver.solve(rotated,configuration:config.solver)
                if let gap = ci.excitationGapInSectorHartree, gap < config.minimumGroundStateGapHartree {
                    throw VivoQMError.convergence("near-degenerate parent state at \(request.points[index].identifier); state tracking is required for QIO")
                }
                let rdm = try VivoFermionRDM.from(ci.state)
                let information = try VivoOrbitalInformation.from(ci.state)
                let leakage = try VivoQIOLeakage.evaluate(rdm:rdm,information:information,partition:request.partition,weights:request.weights)
                objective += request.points[index].weight/totalWeight*leakage.objective
                results.append(.init(identifier:request.points[index].identifier,hamiltonian:rotated,
                                     fullParentEnergyHartree:ci.energyHartree,leakage:leakage,transportSingularValues:singularValues[index]))
            }
            return (objective,rotation,results)
        }
        var accepted = try evaluate(parameters)
        let initial = accepted.0
        var history = [initial]
        var reason: VivoQIOStopReason = count == 0 ? .noRotationDegreesOfFreedom : .iterationBudget
        if count > 0 {
            for _ in 0..<config.maximumIterations {
                if evaluations+2*count > config.maximumEvaluations { reason = .evaluationBudget; break }
                var gradient = Array(repeating:0.0,count:count)
                for i in 0..<count {
                    var plus = parameters, minus = parameters
                    plus[i] += config.finiteDifferenceStep; minus[i] -= config.finiteDifferenceStep
                    gradient[i] = try (evaluate(plus).0-evaluate(minus).0)/(2*config.finiteDifferenceStep)
                }
                let norm = sqrt(gradient.reduce(0) { $0+$1*$1 })
                guard norm.isFinite else { throw VivoQMError.convergence("non-finite path-QIO gradient") }
                if norm < config.gradientTolerance { reason = .gradientConverged; break }
                var step = min(1,config.maximumParameterStep/norm), improved = false
                for _ in 0..<16 {
                    if evaluations >= config.maximumEvaluations { reason = .evaluationBudget; break }
                    let trial = zip(parameters,gradient).map { $0-step*$1 }
                    if trial.contains(where:{ abs($0)>99 }) { step *= 0.5; continue }
                    let evaluated = try evaluate(trial)
                    if evaluated.0 <= accepted.0-1e-4*step*norm*norm {
                        parameters = trial; accepted = evaluated; history.append(accepted.0); improved = true; break
                    }
                    step *= 0.5
                }
                if !improved { if reason != .evaluationBudget { reason = .lineSearchStalled }; break }
            }
        }
        return .init(partition:request.partition,sharedRotation:accepted.1,parameters:parameters,
                     initialObjective:initial,acceptedObjectives:history,evaluations:evaluations,stopReason:reason,points:accepted.2)
    }
}
